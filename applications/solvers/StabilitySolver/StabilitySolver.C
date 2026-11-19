#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"
#include "Switch.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"
#include "Pstream.H"

#include "ShiftedMALSolver.H"
#include "ArpackShiftInvertSolver.H"

using namespace Foam;
using namespace Foam::matrixFreeNewton;
using namespace Foam::matrixFreeStability;

namespace
{
void zeroInternalAndCorrectLocal(volScalarField& f)
{
    scalarField& fi = f.primitiveFieldRef();
    forAll(fi, celli) fi[celli] = 0.0;
    f.correctBoundaryConditions();
}

label localPressureReferenceCell
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
            << "Parallel run requires pRefPoint in system/stabilityProperties"
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
            << " reference owners" << exit(FatalError);
    }

    return ownsRef ? localCell : -1;
}

void writeLeadingMode(const Time& runTime, const ComplexFlow& mode)
{
    const fvMesh& mesh = mode.real.U.mesh();

    volVectorField realU
    (
        IOobject("leadingEigenModeRealU", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mode.real.U
    );
    volVectorField imagU
    (
        IOobject("leadingEigenModeImagU", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mode.imag.U
    );
    volScalarField realP
    (
        IOobject("leadingEigenModeRealP", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mode.real.p
    );
    volScalarField imagP
    (
        IOobject("leadingEigenModeImagP", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mode.imag.p
    );
    surfaceScalarField realPhi
    (
        IOobject("leadingEigenModeRealPhi", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mode.real.phi
    );
    surfaceScalarField imagPhi
    (
        IOobject("leadingEigenModeImagPhi", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mode.imag.phi
    );

    realU.write();
    imagU.write();
    realP.write();
    imagP.write();
    realPhi.write();
    imagPhi.write();
}
}

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "PARPACK shift-invert stability using matrix-free J and mAL/native-GAMG pressure preconditioning"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    const scalar baseDiv = Foam::matrixFreeNewton::maxAbsDivPhi(base.phi);

    Info<< "Stability PARPACK run:" << nl
        << "  Re              = " << Reynolds << nl
        << "  nu              = " << nu.value() << nl
        << "  jacobianSign    = " << jacobianSign << nl
        << "  sigma           = (" << sigmaReal << ',' << sigmaImag << ')' << nl
        << "  base max|divPhi|= " << baseDiv << nl
        << "  mAL gamma       = " << malGamma << nl
        << "  velocity Oseen  = OpenFOAM PBiCGStab/DILU, shift-aware block GS, sweeps="
        << malOseenShiftSweeps << nl
        << "  pressure L^-1   = native OpenFOAM GAMG (one cycle)" << nl
        << "  inner restart   = " << shiftInvertRestart << nl
        << "  inner maxIter   = " << shiftInvertMaxIter << nl
        << "  inner tolerance = " << shiftInvertTolerance << nl;

    if (baseDiv > baseContinuityTolerance)
    {
        FatalErrorInFunction
            << "Base state is not sufficiently divergence-free: max|div(phi)|="
            << baseDiv << " tolerance=" << baseContinuityTolerance
            << exit(FatalError);
    }

    ShiftedMALSolver solver
    (
        base,
        nu,
        continuityWork,
        baseUFaceFrozen,
        momentumFluxWork,
        mALPressurePoisson,
        mALPressurePoissonStore,
        pRefCell,
        pRefValue,
        jacobianSign,
        sigmaReal,
        sigmaImag,
        continuityWeight,
        shiftInvertMaxIter,
        shiftInvertRestart,
        shiftInvertTolerance,
        shiftInvertAbsTolerance,
        malGamma,
        malVelocityBlockMaxIter,
        malVelocityBlockRestart,
        malVelocityBlockTolerance,
        malVelocityBlockAbsTolerance,
        malVelocityBlockPrintInterval,
        malOseenShiftSweeps,
        malRequireVelocityBlockConvergence,
        malReportVelocityBlock
    );
    solver.updateFrozenBaseFaceVelocity();

    ArpackArnoldiControls controls;
    controls.nev = arnoldiNev;
    controls.ncv = arnoldiNcv;
    controls.maxIterations = arnoldiMaxIterations;
    controls.tolerance = arnoldiTolerance;
    controls.which = arnoldiWhich;
    controls.deterministicStart = arnoldiDeterministicStart;
    controls.volumeWeighting = arnoldiVolumeWeighting;
    controls.printInterval = arnoldiPrintInterval;
    controls.innerPrintInterval = shiftInvertPrintInterval;
    controls.leadingSelection = leadingSelection;
    controls.transformedResidualTolerance = arnoldiTransformedResidualTolerance;
    controls.generalizedResidualTolerance = arnoldiGeneralizedResidualTolerance;
    controls.constraintTolerance = arnoldiConstraintTolerance;
    controls.abortOnConstraintViolation = arnoldiAbortOnConstraintViolation;

    ArpackShiftInvertSolver arnoldi
    (
        mesh,
        solver,
        U0,
        phi0,
        p0,
        controls
    );

    ComplexFlow leadingMode("leadingEigenMode", U0, phi0, p0);
    zeroComplexFlow(leadingMode);

    const ArpackArnoldiResult result = arnoldi.solve(leadingMode);

    if (writeEigenModeFields)
    {
        // The base flow is read from the latest time (1/, written by
        // NewtonSolver); the eigenmode is written next to it.  Never into 0/,
        // which holds only the initial conditions.
        if (runTime.value() == 0)
        {
            runTime.setTime(scalar(1), runTime.timeIndex() + 1);
        }
        Info<< "Writing leading eigenmode to time " << runTime.name() << nl;
        writeLeadingMode(runTime, leadingMode);
        Info<< "Wrote leadingEigenMode Real/Imag U,p,phi" << nl;
    }

    Info<< "End" << nl;
    return result.passed ? 0 : 2;
}
