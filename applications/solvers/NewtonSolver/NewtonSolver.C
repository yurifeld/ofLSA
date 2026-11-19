#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"
#include "Switch.H"
#include "fvMatrices.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"
#include "fvcLaplacian.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "Pstream.H"
#include "SolverPerformance.H"
#include "boolList.H"

#include "MatrixFreeFlow.H"
#include "NewtonKrylov.H"
#include "fvcSurfaceIntegrate.H"
#include "fixedValueFvPatchFields.H"
#include <cmath>

using namespace Foam;
using namespace Foam::matrixFreeNewton;

namespace
{
inline bool finiteScalarLocal(const scalar x)
{
    return std::isfinite(static_cast<double>(x));
}

inline scalar robustHypotLocal(const scalar a, const scalar b)
{
    return scalar(std::hypot(static_cast<double>(a), static_cast<double>(b)));
}
}


// Output-only preservation of prescribed velocity boundary data.
// This does NOT participate in the Newton iteration or projection.
static void capturePrescribedVelocityBoundary
(
    const volVectorField& U,
    boolList& prescribedPatch,
    List<vectorField>& prescribedValue
)
{
    const auto& Ub = U.boundaryField();

    prescribedPatch.setSize(Ub.size());
    prescribedValue.setSize(Ub.size());

    forAll(Ub, patchi)
    {
        prescribedPatch[patchi] = Ub[patchi].fixesValue();

        if (prescribedPatch[patchi])
        {
            prescribedValue[patchi] = vectorField(Ub[patchi]);
        }
        else
        {
            prescribedValue[patchi].clear();
        }
    }
}


static void restorePrescribedVelocityBoundary
(
    volVectorField& U,
    const boolList& prescribedPatch,
    const List<vectorField>& prescribedValue
)
{
    auto& Ub = U.boundaryFieldRef();

    forAll(Ub, patchi)
    {
        if (!prescribedPatch[patchi]) continue;

        vectorField& dst = Ub[patchi];
        const vectorField& src = prescribedValue[patchi];

        if (dst.size() != src.size())
        {
            FatalErrorInFunction
                << "Velocity boundary size changed on patch "
                << Ub[patchi].patch().name()
                << exit(FatalError);
        }

        forAll(dst, facei)
        {
            dst[facei] = src[facei];
        }
    }
}


// Zero a velocity field without evaluating boundary conditions.  The boundary
// values are still cleared so repeated accumulation cannot drift, but the
// (communicating) coupled-patch evaluation is left to the single
// setHomogeneousVelocityBoundary() that every caller ends with.
static void zeroVolVectorNoEvaluate(volVectorField& f)
{
    vectorField& fi = f.primitiveFieldRef();
    forAll(fi, celli) fi[celli] = vector::zero;

    forAll(f.boundaryFieldRef(), patchi)
    {
        auto& bf = f.boundaryFieldRef()[patchi];
        forAll(bf, i) bf[i] = vector::zero;
    }
}

// Impose the homogeneous flux boundary condition without communication.
static void zeroPhysicalSurfaceBoundary(surfaceScalarField& f)
{
    forAll(f.boundaryFieldRef(), patchi)
    {
        if (f.boundaryField()[patchi].coupled()) continue;
        auto& bf = f.boundaryFieldRef()[patchi];
        forAll(bf, i) bf[i] = 0.0;
    }
}

static void multiplySurfaceInto
(
    surfaceScalarField& dst,
    const surfaceScalarField& a,
    const surfaceScalarField& b
)
{
    scalarField& d = dst.primitiveFieldRef();
    const scalarField& ai = a.primitiveField();
    const scalarField& bi = b.primitiveField();
    forAll(d, facei) d[facei] = ai[facei]*bi[facei];

    forAll(dst.boundaryFieldRef(), patchi)
    {
        scalarField& db = dst.boundaryFieldRef()[patchi];
        const scalarField& ab = a.boundaryField()[patchi];
        const scalarField& bb = b.boundaryField()[patchi];
        forAll(db, facei) db[facei] = ab[facei]*bb[facei];
    }
}

static void addScaledProductToSurface
(
    surfaceScalarField& dst,
    const scalar coeff,
    const surfaceScalarField& a,
    const surfaceScalarField& b
)
{
    scalarField& d = dst.primitiveFieldRef();
    const scalarField& ai = a.primitiveField();
    const scalarField& bi = b.primitiveField();
    forAll(d, facei) d[facei] += coeff*ai[facei]*bi[facei];

    forAll(dst.boundaryFieldRef(), patchi)
    {
        scalarField& db = dst.boundaryFieldRef()[patchi];
        const scalarField& ab = a.boundaryField()[patchi];
        const scalarField& bb = b.boundaryField()[patchi];
        forAll(db, facei) db[facei] += coeff*ab[facei]*bb[facei];
    }
}

// Patch types of a scalar field representing one Cartesian component of the
// homogeneous velocity perturbation.  A velocity patch type with no scalar
// counterpart could not be reproduced by the single-component operator.
static wordList velocityComponentPatchTypes(const volVectorField& U)
{
    const fvMesh& mesh = U.mesh();
    wordList patchTypes(mesh.boundary().size());

    forAll(mesh.boundary(), patchi)
    {
        const word& type = U.boundaryField()[patchi].type();

        if (!fvPatchField<scalar>::patchConstructorTablePtr_->found(type))
        {
            FatalErrorInFunction
                << "Velocity patch type " << type << " on patch "
                << mesh.boundary()[patchi].name()
                << " has no scalar counterpart, so the mAL single-component "
                << "velocity operator cannot represent it." << nl
                << "Use component-wise velocity boundary conditions "
                << "(e.g. fixedValue/zeroGradient)."
                << exit(FatalError);
        }

        patchTypes[patchi] = type;
    }

    return patchTypes;
}

// Zero every boundary value once.  Physical patch values of the perturbation
// component are identically zero for the whole run, so this never needs
// repeating; coupled patches are refreshed by correctBoundaryConditions().
static void initialiseComponentField(volScalarField& f)
{
    scalarField& fi = f.primitiveFieldRef();
    forAll(fi, celli) fi[celli] = 0.0;

    forAll(f.boundaryFieldRef(), patchi)
    {
        auto& bf = f.boundaryFieldRef()[patchi];
        forAll(bf, i) bf[i] = 0.0;
    }

    f.correctBoundaryConditions();
}

// Load one Cartesian velocity component into its scalar carrier.  Only the
// coupled patches need refreshing: the physical patch values of the
// homogeneous perturbation are identically zero for the whole run.
static void loadComponentField(volScalarField& f, const scalarField& src)
{
    scalarField& fi = f.primitiveFieldRef();
    forAll(fi, celli) fi[celli] = src[celli];
    f.correctBoundaryConditions();
}

static void enforceZeroPressureReference
(
    volScalarField& p,
    const label pRefCell
)
{
    scalar localRef = 0.0;
    if (pRefCell >= 0) localRef = p.primitiveField()[pRefCell];

    const scalar refValue = returnReduce(localRef, sumOp<scalar>());
    scalarField& pi = p.primitiveFieldRef();
    forAll(pi, celli) pi[celli] -= refValue;
    p.correctBoundaryConditions();
}

static label localPressureReferenceCell
(
    const fvMesh& mesh,
    const IOdictionary& dict
)
{
    if (!Pstream::parRun())
    {
        return dict.lookupOrDefault<label>("pRefCell", 0);
    }

    if (!dict.found("pRefPoint"))
    {
        FatalErrorInFunction
            << "Parallel run requires pRefPoint in system/newtonProperties"
            << exit(FatalError);
    }

    const point refPoint(dict.lookup("pRefPoint"));
    scalar localMinDist = GREAT;
    label localCell = -1;

    const vectorField& C = mesh.C();
    forAll(C, celli)
    {
        const scalar d = magSqr(C[celli] - refPoint);
        if (d < localMinDist)
        {
            localMinDist = d;
            localCell = celli;
        }
    }

    const scalar globalMinDist = returnReduce(localMinDist, minOp<scalar>());
    const bool ownsRef =
        localCell >= 0
     && mag(localMinDist - globalMinDist)
        <= SMALL*max(scalar(1), globalMinDist);

    const label nOwners = returnReduce(ownsRef ? 1 : 0, sumOp<label>());
    if (nOwners != 1)
    {
        FatalErrorInFunction
            << "pRefPoint " << refPoint << " selected " << nOwners
            << " reference owners; choose a point clearly inside one cell"
            << exit(FatalError);
    }

    return ownsRef ? localCell : -1;
}

static scalar maxMagVolVectorField(const volVectorField& U)
{
    scalar localMax = 0.0;
    const vectorField& Ui = U.primitiveField();
    forAll(Ui, celli) localMax = max(localMax, mag(Ui[celli]));
    return returnReduce(localMax, maxOp<scalar>());
}

static scalar maxAbsVolScalarField(const volScalarField& f)
{
    scalar localMax = 0.0;
    const scalarField& fi = f.primitiveField();
    forAll(fi, celli) localMax = max(localMax, mag(fi[celli]));
    return returnReduce(localMax, maxOp<scalar>());
}

static void enforcePhysicalPhiFromU(Flow& a)
{
    const surfaceVectorField& Sf = a.U.mesh().Sf();

    forAll(a.phi.boundaryField(), patchi)
    {
        if (a.phi.boundaryField()[patchi].coupled()) continue;

        scalarField& phibf = a.phi.boundaryFieldRef()[patchi];
        const vectorField& Ubf = a.U.boundaryField()[patchi];
        const vectorField& Sfbf = Sf.boundaryField()[patchi];

        forAll(phibf, i) phibf[i] = Ubf[i] & Sfbf[i];
    }
}


// Canonicalise only derived/prescribed volField patch storage.
// The independent nonlinear base unknowns are:
//   U, p  : cell values
//   phi   : face values (including processor-interface faces)
static void canonicaliseBaseVolFields(Flow& base)
{
    base.U.correctBoundaryConditions();
    base.p.correctBoundaryConditions();
    enforcePhysicalPhiFromU(base);
}


// Internal-cell-only U update.
// Used for the nonlinear BASE state, not for generic Krylov-vector algebra.
static void axpyBaseInternalU
(
    volVectorField& y,
    const scalar a,
    const volVectorField& x
)
{
    vectorField& yi = y.primitiveFieldRef();
    const vectorField& xi = x.primitiveField();

    forAll(yi, celli)
    {
        yi[celli] += a*xi[celli];
    }
}


// Internal-cell-only p update.
// Processor/zeroGradient patch values are reconstructed afterwards.
static void axpyBaseInternalP
(
    volScalarField& y,
    const scalar a,
    const volScalarField& x
)
{
    scalarField& yi = y.primitiveFieldRef();
    const scalarField& xi = x.primitiveField();

    forAll(yi, celli)
    {
        yi[celli] += a*xi[celli];
    }
}


// Apply a tmp<volVectorField> to BASE U internal cells only.
static void addBaseInternalFromTmp
(
    volVectorField& y,
    const tmp<volVectorField>& tx,
    const scalar a
)
{
    const volVectorField& x = tx();
    axpyBaseInternalU(y, a, x);
}


// Newton update of the nonlinear BASE state.
//
// IMPORTANT:
// Generic axpyFlow() is still used for Krylov/preconditioner vectors.
// It is deliberately NOT used here, because U and p patch values are not
// independent nonlinear base-state unknowns in a decomposed OpenFOAM case.
static void axpyBaseState
(
    Flow& base,
    const scalar a,
    const Flow& correction
)
{
    axpyBaseInternalU(base.U, a, correction.U);
    axpyBaseInternalP(base.p, a, correction.p);

    // phi is a genuine face unknown.  Keep processor-interface face values
    // in the algebra, then reset only physical boundary flux from U.
    axpySurfScalar(base.phi, a, correction.phi);

    canonicaliseBaseVolFields(base);
}

static void copyInternalAndCorrect
(
    volScalarField& dst,
    const volScalarField& src
)
{
    scalarField& d = dst.primitiveFieldRef();
    const scalarField& s = src.primitiveField();
    forAll(d, celli) d[celli] = s[celli];
    dst.correctBoundaryConditions();
}

static void zeroInternalAndCorrect(volScalarField& dst)
{
    scalarField& d = dst.primitiveFieldRef();
    forAll(d, celli) d[celli] = 0.0;
    dst.correctBoundaryConditions();
}

static void copySurfVector
(
    surfaceVectorField& dst,
    const surfaceVectorField& src
)
{
    vectorField& dI = dst.primitiveFieldRef();
    const vectorField& sI = src.primitiveField();
    forAll(dI, facei) dI[facei] = sI[facei];

    forAll(dst.boundaryField(), patchi)
    {
        vectorField& dB = dst.boundaryFieldRef()[patchi];
        const vectorField& sB = src.boundaryField()[patchi];
        forAll(dB, facei) dB[facei] = sB[facei];
    }
}

static tmp<volVectorField> frozenBaseVelocityFluxDivergence
(
    const surfaceScalarField& perturbationPhi,
    const surfaceVectorField& baseUFaceFrozen,
    surfaceVectorField& momentumFluxWork
)
{
    vectorField& mI = momentumFluxWork.primitiveFieldRef();
    const scalarField& phiI = perturbationPhi.primitiveField();
    const vectorField& uFI = baseUFaceFrozen.primitiveField();
    forAll(mI, facei) mI[facei] = phiI[facei]*uFI[facei];

    forAll(momentumFluxWork.boundaryField(), patchi)
    {
        vectorField& mB = momentumFluxWork.boundaryFieldRef()[patchi];
        const scalarField& phiB = perturbationPhi.boundaryField()[patchi];
        const vectorField& uFB = baseUFaceFrozen.boundaryField()[patchi];
        forAll(mB, facei) mB[facei] = phiB[facei]*uFB[facei];
    }

    return fvc::div(momentumFluxWork);
}

static scalar projectBaseDivergenceFree
(
    Flow& base,
    volScalarField& pProjection,
    const label pRefCell,
    const scalar pRefValue,
    const scalar tolerance,
    const label maxIter
)
{
    canonicaliseBaseVolFields(base);

    scalar divNow = maxAbsVolScalarField(fvc::div(base.phi)());

    for (label iter=0; iter<maxIter && divNow>tolerance; ++iter)
    {
        zeroVolScalar(pProjection);

        tmp<volScalarField> tDivPhi = fvc::div(base.phi);
        fvScalarMatrix pEqn
        (
            fvm::laplacian(pProjection) == tDivPhi()
        );
        pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve();

        // The conservative FV pressure flux is the correction required by
        // the discrete continuity equation.
        tmp<surfaceScalarField> tPhiCorr = pEqn.flux();
        axpySurfScalar(base.phi, -1.0, tPhiCorr());

        // Keep cell-centred velocity and pressure consistent with that flux
        // correction.  For the nonlinear BASE state, U and p are updated on
        // internal cells only.  Their patch/coupled values are reconstructed
        // from the BC objects so the in-memory state is identical in form to
        // a freshly read state on the same decomposition.
        tmp<volVectorField> tGradP = fvc::grad(pProjection);
        addBaseInternalFromTmp(base.U, tGradP, -1.0);
        axpyBaseInternalP(base.p, 1.0, pProjection);

        canonicaliseBaseVolFields(base);

        tmp<volScalarField> tDivNow = fvc::div(base.phi);
        divNow = maxAbsVolScalarField(tDivNow());
    }

    return divNow;
}

class NewtonMALSolver
{
    const Flow& base_;
    const dimensionedScalar& nu_;
    volScalarField& continuityWork_;
    surfaceVectorField& baseUFaceFrozen_;
    surfaceVectorField& momentumFluxWork_;
    label pRefCell_;

    scalar gmresVolumeWeightScale_;
    label itMax_;
    scalar epsilon_;
    scalar absoluteTolerance_;

    scalar malGamma_;
    label malVelocityBlockMaxIter_;
    scalar malVelocityBlockTolerance_;
    scalar malVelocityBlockAbsTolerance_;

    // ---------------------------------------------------------------------- //
    // Cached constants.
    //
    // Sf_i and the directional viscosity nu + gamma (Sf_i/|Sf|)^2 depend only
    // on the mesh, nu and gamma, so they are run constants.  The frozen base
    // face velocity component and the per-component Oseen matrices depend on
    // the base state and are rebuilt once per Newton iteration.
    //
    // The original code recomputed Sf_i/|Sf|, the directional viscosity and a
    // complete fvVectorMatrix on EVERY preconditioner application.
    // ---------------------------------------------------------------------- //
    PtrList<surfaceScalarField> sfComponent_;
    PtrList<surfaceScalarField> directionalNu_;
    PtrList<surfaceScalarField> baseUFaceComponent_;

    // Scalar matrices acting on the single active velocity component.  A
    // vector matrix made OpenFOAM run three segregated component solves, two
    // of them on an identically zero system, and still pay their norm-factor
    // reductions and matrix-vector products.
    PtrList<fvScalarMatrix> componentMatrix_;
    bool componentMatricesValid_;

    ComponentWorkspace malComponentWorkspace_;

    Flow malCanonicalRhsWork_;
    Flow malVelocityActionWork_;

    // Scalar work fields of the single-component operator and of the Oseen
    // component solve.  The single-component action never touches the two
    // inactive velocity components, so it costs a third of the vector action
    // it replaces.
    volScalarField malOseenU_;
    volScalarField componentInput_;
    surfaceScalarField componentPhi_;
    surfaceScalarField componentFluxWork_;

    // Compact right-hand side and solution of the active component block.
    // These must not alias the inner workspace.
    scalarField componentRhsWork_;
    scalarField componentSolutionWork_;

public:
    NewtonMALSolver
    (
        const Flow& base,
        const dimensionedScalar& nu,
        volScalarField& continuityWork,
        surfaceVectorField& baseUFaceFrozen,
        surfaceVectorField& momentumFluxWork,
        const label pRefCell,
        const label gmresMaxIter,
        const scalar gmresTolerance,
        const scalar gmresAbsTolerance,
        const scalar malGamma,
        const label malVelocityBlockMaxIter,
        const scalar malVelocityBlockTolerance,
        const scalar malVelocityBlockAbsTolerance,
        const label malVelocityBlockRestart
    )
    :
        base_(base),
        nu_(nu),
        continuityWork_(continuityWork),
        baseUFaceFrozen_(baseUFaceFrozen),
        momentumFluxWork_(momentumFluxWork),
        pRefCell_(pRefCell),
        gmresVolumeWeightScale_(1.0),
        itMax_(gmresMaxIter),
        epsilon_(gmresTolerance),
        absoluteTolerance_(max(gmresAbsTolerance, SMALL)),
        malGamma_(malGamma),
        malVelocityBlockMaxIter_(malVelocityBlockMaxIter),
        malVelocityBlockTolerance_(malVelocityBlockTolerance),
        malVelocityBlockAbsTolerance_(malVelocityBlockAbsTolerance),
        sfComponent_(3),
        directionalNu_(3),
        baseUFaceComponent_(3),
        componentMatrix_(),
        componentMatricesValid_(false),
        malComponentWorkspace_
        (
            base.U.mesh().nCells(),
            min(malVelocityBlockRestart, malVelocityBlockMaxIter)
        ),
        malCanonicalRhsWork_("mAL_canonicalRhs", base.U, base.phi, base.p),
        malVelocityActionWork_("mAL_velocityAction", base.U, base.phi, base.p),
        malOseenU_
        (
            IOobject
            (
                // The fvSolution entry is looked up by field name.
                "mALOseen_U",
                base.U.mesh().time().name(),
                base.U.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            base.U.mesh(),
            dimensionedScalar("mALOseen_U", base.U.dimensions(), 0.0),
            velocityComponentPatchTypes(base.U)
        ),
        componentInput_
        (
            IOobject
            (
                "mALComponentInput",
                base.U.mesh().time().name(),
                base.U.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            base.U.mesh(),
            dimensionedScalar("mALComponentInput", base.U.dimensions(), 0.0),
            velocityComponentPatchTypes(base.U)
        ),
        componentPhi_
        (
            IOobject
            (
                "mALComponentPhi",
                base.U.mesh().time().name(),
                base.U.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            base.phi
        ),
        componentFluxWork_
        (
            IOobject
            (
                "mALComponentFluxWork",
                base.U.mesh().time().name(),
                base.U.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            base.phi
        ),
        componentRhsWork_(base.U.mesh().nCells(), 0.0),
        componentSolutionWork_(base.U.mesh().nCells(), 0.0)
    {
        if (malGamma_ < 0.0)
        {
            FatalErrorInFunction << "mAL gamma must be >= 0" << exit(FatalError);
        }
        if
        (
            malVelocityBlockMaxIter_ < 1
         || malVelocityBlockTolerance_ <= 0.0
         || malVelocityBlockAbsTolerance_ <= 0.0
         || malComponentWorkspace_.restart < 1
        )
        {
            FatalErrorInFunction << "Invalid mAL velocity-block controls"
                << exit(FatalError);
        }

        const fvMesh& mesh = base_.U.mesh();
        const scalarField& V = mesh.V();

        scalar localVolume = 0.0;
        forAll(V, celli) localVolume += V[celli];

        scalar reduced[2] = {localVolume, scalar(mesh.nCells())};
        globalSumInPlace(reduced, 2);

        if
        (
            reduced[1] <= 0
         || !finiteScalarLocal(reduced[0])
         || reduced[0] <= VSMALL
        )
        {
            FatalErrorInFunction
                << "Cannot construct volume-weighted Krylov metric"
                << exit(FatalError);
        }

        gmresVolumeWeightScale_ = reduced[1]/reduced[0];

        zeroFlow(malCanonicalRhsWork_);
        zeroFlow(malVelocityActionWork_);

        initialiseComponentField(malOseenU_);
        initialiseComponentField(componentInput_);
        zeroSurfScalar(componentPhi_);
        zeroSurfScalar(componentFluxWork_);

        // Sf_i is pure geometry; the directional viscosity depends only on the
        // mesh, nu and gamma.  Both are constant for the whole run.
        const dimensionedScalar gammaNu
        (
            "mALDirectionalGamma",
            nu_.dimensions(),
            malGamma_
        );

        for (direction cmpt=0; cmpt<3; ++cmpt)
        {
            sfComponent_.set
            (
                cmpt,
                new surfaceScalarField
                (
                    IOobject
                    (
                        "mALSf" + name(label(cmpt)),
                        mesh.time().name(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh.Sf().component(cmpt)
                )
            );

            directionalNu_.set
            (
                cmpt,
                new surfaceScalarField
                (
                    IOobject
                    (
                        "mALDirectionalNu" + name(label(cmpt)),
                        mesh.time().name(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    nu_ + gammaNu*sqr(sfComponent_[cmpt]/mesh.magSf())
                )
            );
        }
    }

    // ---------------------------------------------------------------------- //
    // Base-state dependent caches
    // ---------------------------------------------------------------------- //

    void updateFrozenBaseFaceVelocity()
    {
        tmp<surfaceVectorField> tBaseUFace = fvc::interpolate(base_.U);
        copySurfVector(baseUFaceFrozen_, tBaseUFace());
        tBaseUFace.clear();

        updateBaseFaceComponents();

        // base_.phi changed with the Newton update, so the cached Oseen
        // matrices are stale.
        componentMatricesValid_ = false;
    }

    void computeNewtonRHS
    (
        Flow& rhs,
        volScalarField& rhsContinuity
    ) const
    {
        zeroFlow(rhs);
        zeroInternalAndCorrect(rhsContinuity);

        tmp<volVectorField> tGradP = fvc::grad(base_.p);
        tmp<volVectorField> tLapU = fvc::laplacian(nu_, base_.U);
        tmp<volVectorField> tConvU = fvc::div(base_.phi, base_.U);

        addFromTmp(rhs.U, tGradP, 1.0);
        addFromTmp(rhs.U, tLapU, -1.0);
        addFromTmp(rhs.U, tConvU, 1.0);

        tmp<volScalarField> tDivBase = fvc::div(base_.phi);
        copyInternalAndCorrect(rhsContinuity, tDivBase());
        copyInternalAndCorrect(rhs.p, rhsContinuity);
        zeroSurfScalar(rhs.phi);
        setHomogeneousVelocityBoundary(rhs);
    }

    // ---------------------------------------------------------------------- //
    // Outer Newton FGMRES
    // ---------------------------------------------------------------------- //

    bool solveGMRES
    (
        Flow& x,
        const Flow& rhs,
        NewtonStateWorkspace& ws,
        label& maxIterUsed,
        const bool useRelativeTolerance,
        const label printInterval
    )
    {
        const label restart = ws.restart;
        const label printEvery = max(label(1), printInterval);

        auto H = [&](const label row, const label col) -> scalar&
        {
            return ws.H[row + (restart + 1)*col];
        };

        // main() zeroes the Newton correction before every solve, and the
        // action is linear, so A*x is exactly zero and r0 is exactly rhs.
        loadState(ws.W, rhs);

        const scalar initialResidual = residualNorm(ws.W);
        const scalar rawTarget = useRelativeTolerance
          ? epsilon_*max(initialResidual, SMALL)
          : epsilon_;
        const scalar target = max(rawTarget, absoluteTolerance_);

        Info<< "FGMRES: initialResidual=" << initialResidual
            << " target=" << target
            << " rawTarget=" << rawTarget
            << " absFloor=" << absoluteTolerance_
            << " restart=" << restart
            << " maxIter=" << itMax_ << nl;

        if (initialResidual <= target)
        {
            maxIterUsed = 0;
            return true;
        }

        label globalIter = 0;
        scalar trueResidual = initialResidual;

        while (globalIter < itMax_ && trueResidual > target)
        {
            forAll(ws.H, i) ws.H[i] = 0.0;
            forAll(ws.cs, i) ws.cs[i] = 0.0;
            forAll(ws.sn, i) ws.sn[i] = 0.0;
            forAll(ws.g, i) ws.g[i] = 0.0;
            forAll(ws.y, i) ws.y[i] = 0.0;

            const scalar beta = residualNorm(ws.W);
            if (!finiteScalarLocal(beta) || beta < VSMALL) break;

            copyState(ws.V[0], ws.W);
            scaleState(ws.V[0], 1.0/beta);
            ws.g[0] = beta;

            label m = 0;

            for (label j=0; j<restart && globalIter<itMax_; ++j)
            {
                ++globalIter;

                // The mAL preconditioner reads the internal velocity and
                // pressure of its argument only, so no boundary evaluation is
                // needed here.
                storeState(ws.precondInput, ws.V[j]);
                applyMAL(ws.precondInput, ws.precondOutput);

                // Keep the compact image for the solution update, then apply
                // the operator to the full field just produced.
                loadState(ws.Z[j], ws.precondOutput);
                rawMatrixAction(ws.precondOutput, ws.AX);
                loadState(ws.W, ws.AX);

                // Gram-Schmidt with one reorthogonalization (CGS2): one
                // global reduction per pass instead of one per basis vector.
                for (label pass=0; pass<2; ++pass)
                {
                    basisInnerProducts(ws.W, ws.V, j+1, ws.dots.begin());

                    for (label i=0; i<=j; ++i)
                    {
                        scalar h = ws.dots[i];
                        if (!finiteScalarLocal(h)) h = 0.0;

                        H(i,j) += h;
                        axpyState(ws.W, -h, ws.V[i]);
                    }
                }

                H(j+1,j) = residualNorm(ws.W);
                if (finiteScalarLocal(H(j+1,j)) && H(j+1,j) > VSMALL)
                {
                    copyState(ws.V[j+1], ws.W);
                    scaleState(ws.V[j+1], 1.0/H(j+1,j));
                }
                else
                {
                    H(j+1,j) = 0.0;
                    zeroState(ws.V[j+1]);
                }

                for (label i=0; i<j; ++i)
                {
                    const scalar hij = H(i,j);
                    const scalar hi1j = H(i+1,j);
                    H(i,j)   =  ws.cs[i]*hij + ws.sn[i]*hi1j;
                    H(i+1,j) = -ws.sn[i]*hij + ws.cs[i]*hi1j;
                }

                const scalar h1 = H(j,j);
                const scalar h2 = H(j+1,j);
                const scalar den = robustHypotLocal(h1, h2);

                if (finiteScalarLocal(den) && den > VSMALL)
                {
                    ws.cs[j] = h1/den;
                    ws.sn[j] = h2/den;
                }
                else
                {
                    ws.cs[j] = 1.0;
                    ws.sn[j] = 0.0;
                }

                H(j,j) = ws.cs[j]*h1 + ws.sn[j]*h2;
                H(j+1,j) = 0.0;

                const scalar gj = ws.g[j];
                ws.g[j]   =  ws.cs[j]*gj;
                ws.g[j+1] = -ws.sn[j]*gj;

                m = j + 1;
                const scalar estimatedResidual = mag(ws.g[j+1]);

                if (!finiteScalarLocal(estimatedResidual))
                {
                    Info<< "FGMRES: non-finite estimated residual, "
                        << "leaving current restart at iter=" << globalIter
                        << nl;
                    break;
                }

                if
                (
                    globalIter == 1
                 || globalIter % printEvery == 0
                 || estimatedResidual < target
                )
                {
                    Info<< "FGMRES: iter=" << globalIter
                        << " estimatedResidual=" << estimatedResidual
                        << " target=" << target << nl;
                }

                if (estimatedResidual < target) break;
            }

            for (label i=m-1; i>=0; --i)
            {
                scalar sum = ws.g[i];
                for (label k=i+1; k<m; ++k) sum -= H(i,k)*ws.y[k];

                if
                (
                    !finiteScalarLocal(sum)
                 || !finiteScalarLocal(H(i,i))
                 || mag(H(i,i)) < VSMALL
                )
                {
                    ws.y[i] = 0.0;
                }
                else
                {
                    ws.y[i] = sum/H(i,i);
                    if (!finiteScalarLocal(ws.y[i])) ws.y[i] = 0.0;
                }
            }

            for (label i=0; i<m; ++i)
            {
                if (finiteScalarLocal(ws.y[i]))
                {
                    axpyStateIntoFlow(x, ws.y[i], ws.Z[i]);
                }
            }

            // x is the only field whose boundary values and conservative flux
            // are consumed downstream, so rebuild them from the state.
            finaliseFlowState(x);

            rawMatrixAction(x, ws.AX);
            stateEqualsFlowMinusFlow(ws.W, rhs, ws.AX);
            trueResidual = residualNorm(ws.W);

            Info<< "FGMRES: restart end iter=" << globalIter
                << " trueResidual=" << trueResidual
                << " target=" << target << nl;
        }

        maxIterUsed = globalIter;

        const bool converged =
            finiteScalarLocal(trueResidual) && trueResidual <= target;

        if (!converged)
        {
            Info<< "FGMRES: no convergence in " << globalIter
                << " iterations. trueResidual=" << trueResidual
                << " target=" << target << nl;
        }

        return converged;
    }

private:

    // ---------------------------------------------------------------------- //
    // Metric: one global reduction per dot product, one per Gram-Schmidt pass
    // ---------------------------------------------------------------------- //

    scalar innerProduct(const StateVector& A, const StateVector& B) const
    {
        const scalarField& V = base_.U.mesh().V();

        scalar local = 0.0;
        forAll(A.U, celli)
        {
            const scalar w = gmresVolumeWeightScale_*V[celli];
            local += w*((A.U[celli] & B.U[celli]) + A.p[celli]*B.p[celli]);
        }

        globalSumInPlace(&local, 1);
        return local;
    }

    scalar residualNorm(const StateVector& A) const
    {
        const scalar v = innerProduct(A, A);
        if (!finiteScalarLocal(v) || v < 0.0) return GREAT;
        return Foam::sqrt(v);
    }

    void basisInnerProducts
    (
        const StateVector& W,
        const PtrList<StateVector>& V,
        const label n,
        scalar* dots
    ) const
    {
        const scalarField& cellV = base_.U.mesh().V();

        for (label i=0; i<n; ++i)
        {
            const StateVector& Vi = V[i];

            scalar local = 0.0;
            forAll(W.U, celli)
            {
                const scalar w = gmresVolumeWeightScale_*cellV[celli];
                local += w*((W.U[celli] & Vi.U[celli]) + W.p[celli]*Vi.p[celli]);
            }

            dots[i] = local;
        }

        globalSumInPlace(dots, n);
    }

    scalar componentInnerProduct
    (
        const scalarField& A,
        const scalarField& B
    ) const
    {
        const scalarField& V = base_.U.mesh().V();

        scalar local = 0.0;
        forAll(A, celli)
        {
            local += gmresVolumeWeightScale_*V[celli]*A[celli]*B[celli];
        }

        globalSumInPlace(&local, 1);
        return local;
    }

    scalar componentNorm(const scalarField& A) const
    {
        const scalar v = componentInnerProduct(A, A);
        if (!finiteScalarLocal(v) || v < 0.0) return GREAT;
        return Foam::sqrt(v);
    }

    void componentBasisInnerProducts
    (
        const scalarField& W,
        const PtrList<scalarField>& V,
        const label n,
        scalar* dots
    ) const
    {
        const scalarField& cellV = base_.U.mesh().V();

        for (label i=0; i<n; ++i)
        {
            const scalarField& Vi = V[i];

            scalar local = 0.0;
            forAll(W, celli)
            {
                local += gmresVolumeWeightScale_*cellV[celli]*W[celli]*Vi[celli];
            }

            dots[i] = local;
        }

        globalSumInPlace(dots, n);
    }

    // Restore the boundary/flux consistency that the operator needs after the
    // internal Krylov state of a full flow field has been updated.  These are
    // the exact invariants the mAL preconditioner output satisfies.
    static void finaliseFlowState(Flow& x)
    {
        setHomogeneousVelocityBoundary(x);
        x.p.correctBoundaryConditions();
        replacePhiWithFluxOfU(x);
    }

    // ---------------------------------------------------------------------- //
    // Operator
    // ---------------------------------------------------------------------- //

    void rawMatrixAction(const Flow& x, Flow& Ax) const
    {
        applyLinearAction(x, Ax);
        tmp<volScalarField> tDiv = fvc::div(Ax.phi);
        copyInternalAndCorrect(Ax.p, tDiv());
    }

    void applyLinearAction(const Flow& x, Flow& Ax) const
    {
        // Ax.phi and Ax.p are overwritten in full by the caller chain, so only
        // the velocity row needs clearing, and its coupled boundary values are
        // never read before the next full evaluation.
        zeroVolVectorNoEvaluate(Ax.U);

        tmp<volVectorField> tGradP = fvc::grad(x.p);
        tmp<volVectorField> tLapU = fvc::laplacian(nu_, x.U);
        tmp<volVectorField> tConvBase = fvc::div(base_.phi, x.U);
        tmp<volVectorField> tConvNewton =
            frozenBaseVelocityFluxDivergence
            (
                x.phi,
                baseUFaceFrozen_,
                momentumFluxWork_
            );

        addFromTmp(Ax.U, tGradP, 1.0);
        addFromTmp(Ax.U, tLapU, -1.0);
        addFromTmp(Ax.U, tConvBase, 1.0);
        addFromTmp(Ax.U, tConvNewton, 1.0);

        copySurfScalar(Ax.phi, x.phi);
        zeroPhysicalSurfaceBoundary(Ax.phi);
    }

    // Velocity-row-only action used inside mAL, where the input pressure is
    // identically zero and the output pressure/flux rows are discarded.
    void applyMALVelocityRowAction
    (
        const Flow& x,
        volVectorField& AxU
    ) const
    {
        zeroVolVectorNoEvaluate(AxU);

        tmp<volVectorField> tLapU = fvc::laplacian(nu_, x.U);
        tmp<volVectorField> tConvBase = fvc::div(base_.phi, x.U);
        tmp<volVectorField> tConvNewton =
            frozenBaseVelocityFluxDivergence
            (
                x.phi,
                baseUFaceFrozen_,
                momentumFluxWork_
            );

        addFromTmp(AxU, tLapU, -1.0);
        addFromTmp(AxU, tConvBase, 1.0);
        addFromTmp(AxU, tConvNewton, 1.0);
    }

    // Full augmented velocity action.  Used only for the off-diagonal coupling
    // of the block Gauss-Seidel sweep, where the input has more than one
    // active velocity component.  x must already carry homogeneous velocity
    // boundaries and phi=fvc::flux(U); only the velocity rows of Ax are
    // written, its pressure/flux rows are never read.
    void applyMALAugmentedVelocityAction(const Flow& x, Flow& Ax)
    {
        applyMALVelocityRowAction(x, Ax.U);

        if (malGamma_ > SMALL)
        {
            tmp<volScalarField> tDiv = fvc::div(x.phi);
            tmp<volVectorField> tGradDiv = fvc::grad(tDiv());
            tDiv.clear();
            addFromTmp(Ax.U, tGradDiv, -malGamma_);
        }
    }

    // Diagonal block of the augmented velocity operator for one Cartesian
    // component.  The two inactive components of the input are identically
    // zero and the two inactive output components are discarded, so the whole
    // action is evaluated on scalar fields.  Every operator involved (Gauss
    // div/laplacian/grad, linear interpolation) acts component-wise with the
    // same face weights, so this reproduces the corresponding component of the
    // full vector action; the two surface integrals are summed before
    // integration, which is equivalent up to rounding.
    void applyMALDiagonalComponentAction
    (
        const direction cmpt,
        const scalarField& x,
        scalarField& Ax
    )
    {
        loadComponentField(componentInput_, x);

        const surfaceScalarField& Sfi = sfComponent_[cmpt];
        const surfaceScalarField& U0i = baseUFaceComponent_[cmpt];

        // phi' = fvc::flux(u e_i) = interpolate(u) Sf_i
        {
            tmp<surfaceScalarField> tI = fvc::interpolate(componentInput_);
            multiplySurfaceInto(componentPhi_, tI(), Sfi);
        }

        // -laplacian(nu,u) + div(phi0,u)
        {
            tmp<volScalarField> tLap = fvc::laplacian(nu_, componentInput_);
            tmp<volScalarField> tConv = fvc::div(base_.phi, componentInput_);

            const scalarField& lap = tLap().primitiveField();
            const scalarField& conv = tConv().primitiveField();
            forAll(Ax, celli)
            {
                Ax[celli] = -lap[celli] + conv[celli];
            }
        }

        // + div(phi' U0_i(face)) - gamma grad(div(phi'))_i
        //
        // Both terms are surface integrals, so they are assembled into one
        // face field and integrated once.
        {
            multiplySurfaceInto(componentFluxWork_, componentPhi_, U0i);

            if (malGamma_ > SMALL)
            {
                tmp<volScalarField> tDiv = fvc::div(componentPhi_);
                tmp<surfaceScalarField> tI = fvc::interpolate(tDiv());
                tDiv.clear();

                addScaledProductToSurface
                (
                    componentFluxWork_,
                    -malGamma_,
                    Sfi,
                    tI()
                );
            }

            tmp<volScalarField::Internal> tIntegral =
                fvc::surfaceIntegrate(componentFluxWork_);

            const scalarField& si = tIntegral();
            forAll(Ax, celli) Ax[celli] += si[celli];
        }
    }

    // ---------------------------------------------------------------------- //
    // Deepest velocity preconditioner
    // ---------------------------------------------------------------------- //

    void updateBaseFaceComponents()
    {
        const fvMesh& mesh = base_.U.mesh();

        for (direction cmpt=0; cmpt<3; ++cmpt)
        {
            if (!baseUFaceComponent_.set(cmpt))
            {
                baseUFaceComponent_.set
                (
                    cmpt,
                    new surfaceScalarField
                    (
                        IOobject
                        (
                            "mALBaseUFace" + name(label(cmpt)),
                            mesh.time().name(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        baseUFaceFrozen_.component(cmpt)
                    )
                );
            }
            else
            {
                baseUFaceComponent_[cmpt] = baseUFaceFrozen_.component(cmpt);
            }
        }
    }

    // Assemble the constant-per-Newton-iteration Oseen operators once instead
    // of on every preconditioner application.  Only the source changes between
    // solves.
    void updateComponentMatrices()
    {
        // The boundary coefficients are taken from the (identically zero)
        // homogeneous perturbation boundary values, so assemble from a zero
        // field.
        initialiseComponentField(malOseenU_);

        componentMatrix_.clear();
        componentMatrix_.setSize(3);

        for (direction cmpt=0; cmpt<3; ++cmpt)
        {
            // Retain the diagonal conservative Newton cross-convection in the
            // deepest component preconditioner: div(Sf_i U0_i,f * u_i).
            tmp<surfaceScalarField> tCross
            (
                sfComponent_[cmpt]*baseUFaceComponent_[cmpt]
            );

            // The implicit div/laplacian operators contribute nothing to
            // source_, so the freshly built source_ is exactly zero; the solve
            // then writes V*rhs into it, reproducing "== momentumSource".
            tmp<fvScalarMatrix> tUEqn
            (
                fvm::div(base_.phi, malOseenU_)
              + fvm::div(tCross(), malOseenU_)
              - fvm::laplacian(directionalNu_[cmpt], malOseenU_)
            );

            componentMatrix_.set(cmpt, tUEqn.ptr());
        }

        componentMatricesValid_ = true;
    }

    void applyMALOseenComponentPreconditioner
    (
        const direction cmpt,
        const scalarField& rhs,
        scalarField& z
    )
    {
        if (!componentMatricesValid_) updateComponentMatrices();

        fvScalarMatrix& UEqn = componentMatrix_[cmpt];

        const scalarField& V = base_.U.mesh().V();
        scalarField& psi = malOseenU_.primitiveFieldRef();
        scalarField& source = UEqn.source();

        // Deterministic zero initial guess.  Only the internal field is used
        // by the solve: the physical boundary values stay zero for the whole
        // run and the coupled ones are exchanged by the solver itself.
        forAll(psi, celli) psi[celli] = 0.0;
        forAll(source, celli) source[celli] = V[celli]*rhs[celli];

        const int oldDebug = SolverPerformance<scalar>::debug;
        SolverPerformance<scalar>::debug = 0;
        UEqn.solve();
        SolverPerformance<scalar>::debug = oldDebug;

        forAll(z, celli) z[celli] = psi[celli];
    }

    // ---------------------------------------------------------------------- //
    // Inner one-component FGMRES, entirely on compact scalar fields
    // ---------------------------------------------------------------------- //

    bool solveMALDiagonalComponent
    (
        const direction cmpt,
        const scalarField& rhs,
        scalarField& x,
        label& iterations
    )
    {
        ComponentWorkspace& ws = malComponentWorkspace_;
        const label restart = ws.restart;

        auto H = [&](const label row, const label col) -> scalar&
        {
            return ws.H[row + (restart + 1)*col];
        };

        // x starts from zero and the action is linear, so A*x is exactly zero
        // and r0 is exactly rhs.
        ws.X = scalar(0);
        forAll(ws.RHS, celli)
        {
            ws.RHS[celli] = rhs[celli];
            ws.W[celli] = rhs[celli];
        }

        const scalar initialResidual = componentNorm(ws.W);
        const scalar target = max
        (
            malVelocityBlockTolerance_*max(initialResidual, SMALL),
            malVelocityBlockAbsTolerance_
        );

        if (initialResidual <= target)
        {
            iterations = 0;
            forAll(x, celli) x[celli] = ws.X[celli];
            return true;
        }

        label globalIter = 0;
        scalar trueResidual = initialResidual;

        while (globalIter < malVelocityBlockMaxIter_ && trueResidual > target)
        {
            forAll(ws.H, i) ws.H[i] = 0.0;
            forAll(ws.cs, i) ws.cs[i] = 0.0;
            forAll(ws.sn, i) ws.sn[i] = 0.0;
            forAll(ws.g, i) ws.g[i] = 0.0;
            forAll(ws.y, i) ws.y[i] = 0.0;

            const scalar beta = componentNorm(ws.W);
            if (!finiteScalarLocal(beta) || beta < VSMALL) break;

            {
                scalarField& V0 = ws.V[0];
                const scalar s0 = 1.0/beta;
                forAll(V0, celli) V0[celli] = s0*ws.W[celli];
            }
            ws.g[0] = beta;

            label m = 0;
            for
            (
                label j=0;
                j<restart && globalIter<malVelocityBlockMaxIter_;
                ++j
            )
            {
                ++globalIter;

                applyMALOseenComponentPreconditioner(cmpt, ws.V[j], ws.Z[j]);
                applyMALDiagonalComponentAction(cmpt, ws.Z[j], ws.W);

                // Always reorthogonalize (CGS2): one global reduction per pass
                // instead of one per basis vector.
                for (label pass=0; pass<2; ++pass)
                {
                    componentBasisInnerProducts
                    (
                        ws.W,
                        ws.V,
                        j+1,
                        ws.dots.begin()
                    );

                    for (label i=0; i<=j; ++i)
                    {
                        scalar h = ws.dots[i];
                        if (!finiteScalarLocal(h)) h = 0.0;

                        H(i,j) += h;

                        const scalarField& Vi = ws.V[i];
                        forAll(ws.W, celli) ws.W[celli] -= h*Vi[celli];
                    }
                }

                H(j+1,j) = componentNorm(ws.W);
                if (finiteScalarLocal(H(j+1,j)) && H(j+1,j) > VSMALL)
                {
                    scalarField& Vn = ws.V[j+1];
                    const scalar sn = 1.0/H(j+1,j);
                    forAll(Vn, celli) Vn[celli] = sn*ws.W[celli];
                }
                else
                {
                    H(j+1,j) = 0.0;
                    ws.V[j+1] = scalar(0);
                }

                for (label i=0; i<j; ++i)
                {
                    const scalar hij = H(i,j);
                    const scalar hi1j = H(i+1,j);
                    H(i,j)   =  ws.cs[i]*hij + ws.sn[i]*hi1j;
                    H(i+1,j) = -ws.sn[i]*hij + ws.cs[i]*hi1j;
                }

                const scalar h1 = H(j,j);
                const scalar h2 = H(j+1,j);
                const scalar den = robustHypotLocal(h1, h2);

                if (finiteScalarLocal(den) && den > VSMALL)
                {
                    ws.cs[j] = h1/den;
                    ws.sn[j] = h2/den;
                }
                else
                {
                    ws.cs[j] = 1.0;
                    ws.sn[j] = 0.0;
                }

                H(j,j) = ws.cs[j]*h1 + ws.sn[j]*h2;
                H(j+1,j) = 0.0;
                const scalar gj = ws.g[j];
                ws.g[j]   =  ws.cs[j]*gj;
                ws.g[j+1] = -ws.sn[j]*gj;

                m = j+1;
                if (!finiteScalarLocal(ws.g[j+1]) || mag(ws.g[j+1]) <= target)
                {
                    break;
                }
            }

            for (label i=m-1; i>=0; --i)
            {
                scalar sum = ws.g[i];
                for (label k=i+1; k<m; ++k) sum -= H(i,k)*ws.y[k];

                if
                (
                    !finiteScalarLocal(sum)
                 || !finiteScalarLocal(H(i,i))
                 || mag(H(i,i)) < VSMALL
                )
                {
                    ws.y[i] = 0.0;
                }
                else
                {
                    ws.y[i] = sum/H(i,i);
                    if (!finiteScalarLocal(ws.y[i])) ws.y[i] = 0.0;
                }
            }

            for (label i=0; i<m; ++i)
            {
                if (finiteScalarLocal(ws.y[i]))
                {
                    const scalar yi = ws.y[i];
                    const scalarField& Zi = ws.Z[i];
                    forAll(ws.X, celli) ws.X[celli] += yi*Zi[celli];
                }
            }

            applyMALDiagonalComponentAction(cmpt, ws.X, ws.AX);
            forAll(ws.W, celli) ws.W[celli] = ws.RHS[celli] - ws.AX[celli];
            trueResidual = componentNorm(ws.W);
        }

        iterations = globalIter;
        forAll(x, celli) x[celli] = ws.X[celli];
        return trueResidual <= target;
    }

    // ---------------------------------------------------------------------- //
    // mAL right preconditioner
    // ---------------------------------------------------------------------- //

    void applyMAL(const Flow& r, Flow& z)
    {
        // Only the (U,p) rows of r carry information.
        copyVolVector(malCanonicalRhsWork_.U, r.U);
        copyVolScalar(malCanonicalRhsWork_.p, r.p);

        // Consistent augmented RHS: f_gamma = f - gamma grad(g).
        copyInternalAndCorrect(continuityWork_, r.p);
        tmp<volVectorField> tGradG = fvc::grad(continuityWork_);
        addFromTmp(malCanonicalRhsWork_.U, tGradG, -malGamma_);
        tGradG.clear();

        // Block Gauss-Seidel over Ux, Uy, Uz.  z accumulates the velocity
        // solution in place, so the coupling action can read it directly: the
        // components not yet solved are still exactly zero.
        zeroVolVectorNoEvaluate(z.U);

        scalarField& rhsC = componentRhsWork_;
        label iterations = 0;

        for (direction cmpt=0; cmpt<3; ++cmpt)
        {
            const vectorField& canon =
                malCanonicalRhsWork_.U.primitiveField();

            if (cmpt == 0)
            {
                forAll(rhsC, celli) rhsC[celli] = canon[celli][cmpt];
            }
            else
            {
                // Off-diagonal coupling of the already solved components.
                setHomogeneousVelocityBoundary(z);
                replacePhiWithFluxOfU(z);
                applyMALAugmentedVelocityAction(z, malVelocityActionWork_);

                const vectorField& act =
                    malVelocityActionWork_.U.primitiveField();
                forAll(rhsC, celli)
                {
                    rhsC[celli] = canon[celli][cmpt] - act[celli][cmpt];
                }
            }

            solveMALDiagonalComponent
            (
                cmpt,
                rhsC,
                componentSolutionWork_,
                iterations
            );

            vectorField& zu = z.U.primitiveFieldRef();
            forAll(zu, celli)
            {
                zu[celli][cmpt] = componentSolutionWork_[celli];
            }
        }

        setHomogeneousVelocityBoundary(z);
        replacePhiWithFluxOfU(z);

        // Schur approximation: Sp^{-1} ~= -(gamma+nu) I applied to
        // h = div(u_tilde) - g, then fix the pressure gauge.
        tmp<volScalarField> tDiv = fvc::div(z.phi);
        {
            const scalar a = -(malGamma_ + nu_.value());
            const scalarField& d = tDiv().primitiveField();
            const scalarField& g = malCanonicalRhsWork_.p.primitiveField();
            scalarField& zp = z.p.primitiveFieldRef();
            forAll(zp, celli) zp[celli] = a*(d[celli] - g[celli]);
        }
        tDiv.clear();

        enforceZeroPressureReference(z.p, pRefCell_);
    }
};

int main(int argc, char *argv[])
{
    argList::addNote("Steady Newton solver with mAL right preconditioning");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    boolList prescribedVelocityPatch;
    List<vectorField> prescribedVelocityValue;
    capturePrescribedVelocityBoundary
    (
        U0,
        prescribedVelocityPatch,
        prescribedVelocityValue
    );

    NewtonMALSolver solver
    (
        base,
        nu,
        continuityWork,
        baseUFaceFrozen,
        momentumFluxWork,
        pRefCell,
        gmresMaxIter,
        gmresTolerance,
        gmresAbsTolerance,
        malGamma,
        malVelocityBlockMaxIter,
        malVelocityBlockTolerance,
        malVelocityBlockAbsTolerance,
        malVelocityBlockRestart
    );

    // Projection is part of the nonlinear algorithm, not a Newton
    // preconditioner.  It keeps the stored FV flux discretely solenoidal.
    projectBaseDivergenceFree
    (
        base,
        pProjection,
        pRefCell,
        pRefValue,
        projectionTolerance,
        projectionMaxIter
    );

    scalar nonlinearMomentum = GREAT;
    scalar nonlinearContinuity = GREAT;
    bool converged = false;

    for (label iterNewton=0; iterNewton<maxNewtonIter; ++iterNewton)
    {
        // The nonlinear residual is always evaluated on the same canonical
        // representation that will be obtained after writing and rereading:
        // U,p patch storage is reconstructed; physical phi follows U.
        canonicaliseBaseVolFields(base);

        solver.updateFrozenBaseFaceVelocity();
        solver.computeNewtonRHS(rhsNewton, rhsContinuity);

        nonlinearMomentum = maxMagVolVectorField(rhsNewton.U);
        nonlinearContinuity = maxAbsVolScalarField(rhsContinuity);

        Info<< "Newton " << iterNewton + 1
            << ": momentum=" << nonlinearMomentum
            << " continuity=" << nonlinearContinuity << nl;

        if
        (
            nonlinearMomentum < newtonTolerance
         && nonlinearContinuity < newtonContinuityTolerance
        )
        {
            converged = true;
            break;
        }

        zeroFlow(newtonCorrection);
        label nLinearIter = 0;
        const bool linearConverged = solver.solveGMRES
        (
            newtonCorrection,
            rhsNewton,
            gmresWorkspace,
            nLinearIter,
            true,
            gmresPrintInterval
        );

        if (!linearConverged)
        {
            FatalErrorInFunction
                << "Newton FGMRES failed after " << nLinearIter
                << " iterations. Correction not applied."
                << exit(FatalError);
        }

        // Apply damping directly to the BASE update.  This avoids using
        // generic Flow algebra on base U/p boundary storage.
        const scalar baseUpdateScale = -newtonDamping;
        axpyBaseState(base, baseUpdateScale, newtonCorrection);

        projectBaseDivergenceFree
        (
            base,
            pProjection,
            pRefCell,
            pRefValue,
            projectionTolerance,
            projectionMaxIter
        );
    }

    if (!converged)
    {
        WarningInFunction
            << "Newton did not reach the requested nonlinear tolerance. "
            << "Final momentum=" << nonlinearMomentum
            << " continuity=" << nonlinearContinuity << nl;
    }

    // At convergence the base is already canonical.  Refresh once more
    // before copying to disk; this should not change the nonlinear residual.
    canonicaliseBaseVolFields(base);

    copyVolVector(U0, base.U);

    // Safety-only restoration of the originally prescribed fixed-value U.
    // With internal-only base algebra this should now be a no-op.
    restorePrescribedVelocityBoundary
    (
        U0,
        prescribedVelocityPatch,
        prescribedVelocityValue
    );

    copyVolScalar(p0, base.p);
    copySurfScalar(phi0, base.phi);

    // Results go to time 1: the initial conditions in 0/ stay untouched, and
    // ParaView (which skips time 0 by default) shows the base flow right away.
    // A continuation run starts from time 1 (startFrom latestTime) and
    // overwrites it.
    if (runTime.value() == 0)
    {
        runTime.setTime(scalar(1), runTime.timeIndex() + 1);
    }
    Info<< "Writing base flow to time " << runTime.name() << nl;

    U0.write();
    p0.write();
    phi0.write();

    Info<< "End" << nl;
    return 0;
}
