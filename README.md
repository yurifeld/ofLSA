# ofLSA: Preconditioned Shift-and-Invert Global Linear Stability Analysis in OpenFOAM

ofLSA is a parallel, matrix-free computational framework for the steady-state solution and
three-dimensional global linear stability analysis (TriGlobal LSA) of incompressible flows.
Developed by D. Poliak and Y. Feldman at the Department of Mechanical Engineering,
Ben-Gurion University of the Negev, the framework is implemented in OpenFOAM 14 (OpenFOAM Foundation) and uses its
finite-volume discretization, domain decomposition, MPI communication and native iterative
linear solvers.

The repository accompanies the manuscript *Preconditioned Shift-and-Invert Global Linear
Stability Analysis in OpenFOAM* by D. Poliak and Y. Feldman, currently under review. It
provides the solver implementations, installation,
compilation and execution scripts, a reduced-resolution tutorial for initial verification,
and the high-resolution straight and diagonally lid-driven cavity configurations used in the
study.

---

## 1. Numerical methodology

The framework comprises two solvers sharing the same discrete linearization of the
incompressible Navier–Stokes equations and a common modified augmented-Lagrangian (mAL)
preconditioning strategy.

### Steady base-flow computation (`NewtonSolver`)

The steady finite-volume equations are solved directly by Newton iteration, allowing the
computation of steady base states irrespective of their temporal stability, provided the
nonlinear iteration converges. At each Newton step, the coupled velocity–pressure correction
satisfies

$$
\mathbf{J}(\mathbf{Q}^{k})\,\delta\mathbf{q}^{k} = -\mathbf{r}(\mathbf{Q}^{k}),
\qquad
\mathbf{J} = \begin{bmatrix} \mathbf{F} & \mathbf{G} \\ \mathbf{D} & 0 \end{bmatrix},
$$

where $\mathbf{Q}^k$ is the current base-state approximation, $\mathbf{r}$ is the discrete
residual, $\mathbf{F}$ is the full Newton velocity block, and $\mathbf{G}$ and $\mathbf{D}$
are the discrete gradient and divergence operators, respectively. The Jacobian action
includes both contributions from the linearization of convection. The coupled system is
solved using right-preconditioned flexible GMRES (FGMRES). A pressure-Poisson projection is
used to control the discrete continuity error before the Newton iteration and following each
Newton update.

### Global linear stability analysis (`StabilitySolver`)

Linearization about the converged steady base state yields the generalized eigenvalue problem

$$
\mathbf{L}\widehat{\mathbf{q}} = \lambda\mathbf{B}\widehat{\mathbf{q}},
\qquad
\mathbf{L} = -\mathbf{J},
\qquad
\mathbf{B} = \begin{bmatrix} \mathbf{I}_u & 0 \\ 0 & 0 \end{bmatrix},
$$

where $\widehat{\mathbf{q}}$ contains the velocity and pressure eigenmode components and
$\lambda = \lambda_r + i\omega$ denotes the temporal growth rate and angular frequency of the
perturbation. The finite eigenvalues near a prescribed complex shift $\sigma$ are extracted
using distributed Arnoldi iteration provided by PARPACK, applied to the shift-and-invert
operator

$$
\mathbf{T}_{\sigma} = (\mathbf{L} - \sigma\mathbf{B})^{-1}\mathbf{B},
\qquad
\mu = \frac{1}{\lambda - \sigma}.
$$

Each application of this operator requires the iterative solution of a complex coupled
velocity–pressure system. The shifted operator is evaluated using a doubled-real
representation. Only the velocity components are used to construct the Arnoldi basis; the
associated pressure eigenmode is recovered following convergence by an additional
shifted-system solve.

### Modified augmented-Lagrangian preconditioning

Both solvers employ an mAL right preconditioner based on the augmented velocity operator

$$
\mathbf{F}_{\gamma,\sigma} = \mathbf{F} - \gamma\mathbf{G}\mathbf{D} + \sigma\mathbf{I}_u,
$$

with $\sigma = 0$ for the steady Newton system. A lower-triangular approximation of the
augmented velocity block permits sequential solution of the scalar component systems. The
inverse pressure Schur complement is approximated by

$$
\widetilde{\mathbf{S}}_{\gamma,\sigma}^{-1} = -(Re^{-1} + \gamma)\mathbf{I}_p - \sigma\mathbf{L}_p^{-1},
$$

where $\mathbf{L}_p$ denotes the positive pressure-Laplacian operator with the pressure
reference imposed.

The nested iterative solution procedure uses OpenFOAM's PBiCGStab/DILU infrastructure for
scalar velocity preconditioning and GAMG for pressure-Poisson operations. The coupled Newton
Jacobian and shifted stability operators are evaluated directly from the finite-volume
discretization without assembling or factorizing the global coupled matrices. Sparse scalar
auxiliary matrices are assembled within the preconditioning and projection procedures. The
approximations introduced by mAL are confined to the preconditioner; the original discrete
operators are retained in the outer Krylov and Arnoldi iterations.

The mathematical derivation and details of the nested solver hierarchy are provided in the
accompanying manuscript.

---

## 2. Reference configurations and reproducibility

The framework was verified for the straight and diagonally lid-driven cubic cavities using
structured, wall-clustered grids of $256^3$ control volumes. The grid mapping

$$
x = \xi - 0.0975\sin(2\pi\xi), \qquad 0 \le \xi \le 1,
$$

is applied in each of the three spatial directions. The same augmentation parameter,
$\gamma = 0.02$, was used for both cavity configurations and for the steady and shifted
systems.

| Configuration | Grid | Near-critical Reynolds number | Angular frequency |
|---|---|---|---|
| Straight lid-driven cavity | $256^3$ | 1928.75 | 0.5841 |
| Diagonally lid-driven cavity | $256^3$ | 2320.8 | 0.2507 |

The reported Reynolds numbers correspond to computed near-neutral states on the specified
finite-volume grids. The critical eigenmode of the straight cavity preserves midplane
reflection symmetry, whereas that of the diagonally driven cavity is antisymmetric with
respect to its diagonal symmetry plane. The latter implies half-period spatio-temporal
reflection symmetry of the neutral linear disturbance.

The reference computations employed 256 CPU cores on a Lenovo ThinkSystem SR645 V3 server
equipped with two 128-core AMD EPYC 9745 processors (nominal frequency 2.4 GHz) and 1536 GB
of DDR5 memory (rated at 6400 MT/s). On the $256^3$ grids, the steady base-flow computations
required at most approximately two hours, while the complete shift-and-invert Arnoldi
calculations required approximately 14–30 hours, depending on the configuration and
operating point. These timings characterize the reported computational setup and are not
general performance guarantees.

The repository provides two levels of verification:

- `tutorials/lidDrivenCavity3D/`, `tutorials/diagonalCavity3D/`: $50^3$ straight- and
  diagonal-cavity tutorials for verifying
  compilation, parallel execution and the numerical workflow (Section 5). Its eigenvalues are
  resolution-dependent and should not be expected to reproduce the high-resolution reference
  values.
- `tutorials/paper/`: the two $256^3$ configurations corresponding to the manuscript,
  including the prescribed Reynolds numbers, mesh and shift settings. Follow the
  case-specific instructions in [`tutorials/paper/README.md`](tutorials/paper/README.md) to
  reproduce the reported calculations.

---

## 3. Software requirements

- Linux workstation or HPC system.
- **OpenFOAM 14 from the OpenFOAM Foundation** (openfoam.org). Compatibility with the
  separately maintained OpenFOAM distribution from openfoam.com is not established.
- A C++ toolchain compatible with OpenFOAM, a Fortran compiler and an MPI implementation
  (OpenMPI, Intel MPI or another MPI supported by the OpenFOAM installation).
- ARPACK and its parallel extension PARPACK, built with an MPI and compiler configuration
  compatible with OpenFOAM. If not available, the build procedure compiles them locally.
- Python 3 for the supplied mesh-generation utility and ParaView for visualization of the
  computed fields.

Two installation routes are documented: an Ubuntu workstation with administrator privileges
(Section 4.1), and an HPC system without administrator privileges (Section 4.2).
Alternatively, a container is provided (Section 9).

---

## 4. Installation and compilation

### 4.1 Ubuntu workstation with administrator (`sudo`) privileges

Supported: Ubuntu 22.04, 24.04 and 26.04 (64-bit).

```bash
git clone https://github.com/yurifeld/ofLSA.git
cd ofLSA
scripts/install_ubuntu.sh
```

The script installs the missing system components — OpenFOAM 14 from the OpenFOAM Foundation
package repository, OpenMPI, gfortran, ARPACK/PARPACK — and subsequently compiles the
framework. Components already present are left unchanged, so the script may be re-run
safely; `scripts/install_ubuntu.sh --check` only reports the installation status. The build
log is written to `log.Allwmake`.

To initialize the environment in a new terminal:

```bash
scripts/oflsa-shell
```

This opens a shell with OpenFOAM 14 and the project environment loaded (prompt
`(ofLSA) ...$`). Execute the remaining commands in this shell. The utility also avoids two
known problems of OpenFOAM environment initialization on Ubuntu desktops (Section 8).

### 4.2 HPC system without administrator privileges

On shared HPC systems the operating system, compilers and MPI are site-specific, and
OpenFOAM must be compiled for that environment. Its installation should therefore be
performed by the system administrator. A suitable request is:

> Please install **OpenFOAM 14 (OpenFOAM Foundation, openfoam.org)**, compiled with the
> system's MPI, and indicate how to load it (an environment module, or the path to
> `OpenFOAM-14/etc/bashrc`). An MPI Fortran compiler (`mpif90`) and `cmake` are also
> required. If ARPACK-ng with MPI support (PARPACK) is available as a module, please
> indicate its name.

All remaining components are installed in the user's home directory. The commands below
assume the `bash` shell; users of `tcsh`/`csh` should start `bash` first.

```bash
bash                                          # only if the login shell is tcsh/csh
module load openfoam/14                       # as specified by the administrator, or:
source /path/to/OpenFOAM-14/etc/bashrc
git clone https://github.com/yurifeld/ofLSA.git
cd ofLSA
source etc/bashrc
./Allwmake -j
```

`Allwmake` uses an existing PARPACK if found (`ARPACK_DIR` set by an environment module, or
system packages); otherwise it downloads and compiles ARPACK-ng — and OpenBLAS, if no
BLAS/LAPACK is available — into `~/software/arpack-ng-install`. This is required only once.

In each new session (and in each scheduler job script), initialize the environment before
execution:

```bash
module load openfoam/14                       # or: source /path/to/OpenFOAM-14/etc/bashrc
source ~/ofLSA/etc/bashrc
```

Example SLURM and SGE job scripts are given in
[`tutorials/paper/README.md`](tutorials/paper/README.md).

### Successful compilation

Both routes compile the shared matrix-free FGMRES library and the two solver applications.
Successful compilation ends with:

```
Built:
.../bin/NewtonSolver
.../bin/StabilitySolver
```

---

## 5. Execution of the tutorial

With the project environment initialized:

```bash
cd tutorials/lidDrivenCavity3D
./Allrun -np 4 > log.Allrun 2>&1
```

The `-np` argument specifies the number of MPI processes; use at most the number of physical
cores (not hardware threads). For a detached local execution:

```bash
nohup ./Allrun -np 4 > log.Allrun 2>&1 &
tail -f log.Allrun
```

On managed HPC systems, execute the workflow within a resource allocation using the site's
scheduler (Section 4.2).

The `Allrun` script performs the following operations:

1. Removal of earlier results (time directories other than `0/`). The initial and boundary
   conditions in `0/` are never modified: both solvers write their results into time `1/`.
2. Mesh generation using `blockMesh`.
3. Domain decomposition using `decomposePar`.
4. Computation of the steady base flow with `NewtonSolver`, using continuation in Reynolds
   number (200, 500, 1000, then increments of 250, up to the target) to improve nonlinear
   convergence.
5. Global linear stability analysis of the converged base flow at the target Reynolds number
   using `StabilitySolver`.
6. Reconstruction of the distributed solution into `1/` using `reconstructPar -latestTime`, and
   creation of `case.foam` for ParaView. ParaView opens the case directly at time 1; the
   paper figures are produced from it by the scripts in `figures/` (Section 7a).

`./Allclean` removes all generated data and restores the case to its original state.

With the supplied settings (target $Re = 2100$, shift $\sigma = 0.5i$) the complete tutorial
requires approximately 20 minutes on four cores of a current workstation, and yields a
leading eigenvalue of approximately $\lambda \approx -0.0049 + 0.5479\,i$. Small
differences in the last digits are expected between platforms and process counts. The
preliminary tutorial is intended to verify the complete workflow before undertaking the
substantially more expensive high-resolution calculations.

---

## 6. Case configuration

The tutorial follows the standard OpenFOAM case organization:

```
tutorials/lidDrivenCavity3D/
├── 0/
│   ├── U
│   └── p
├── system/
│   ├── blockMeshDict
│   ├── newtonProperties
│   ├── stabilityProperties
│   ├── decomposeParDict
│   ├── controlDict
│   ├── fvSchemes
│   └── fvSolution
├── Allrun
└── Allclean
```

| File | Purpose |
|---|---|
| `0/U`, `0/p` | Initial and boundary conditions (not modified by the solvers) |
| `1/` | Generated: base flow (`U`, `p`, `phi`) and eigenmode fields (`leadingEigenMode*`) |
| `system/blockMeshDict` | Mesh geometry and resolution |
| `system/newtonProperties` | Reynolds number, nonlinear and linear solver controls, and Newton preconditioning parameters |
| `system/stabilityProperties` | Reynolds number, spectral shift, Arnoldi controls, and shifted-system solver parameters |
| `system/decomposeParDict` | MPI domain decomposition (set by `Allrun -np`) |
| `system/fvSchemes`, `system/fvSolution`, `system/controlDict` | Finite-volume discretization and OpenFOAM execution controls |
| `Allrun`, `Allclean` | Execution and removal of generated case data, respectively |

### Reynolds number and mesh resolution

For the nondimensional cavity cases, the characteristic lid speed and cavity side length are
unity, so that the nondimensional kinematic viscosity is $Re^{-1}$.

To modify a case while preserving the supplied tutorial, first create a working copy from the
repository root:

```bash
cp -r tutorials/lidDrivenCavity3D ~/myCavity
cd ~/myCavity
./Allclean
```

Set the target Reynolds number and execute the case:

```bash
setRe 1900
./Allrun -np 8
```

The `setRe` utility updates both `newtonProperties` and `stabilityProperties`. The Reynolds
number specified for stability analysis must coincide with that of the converged base flow.

To generate a $100^3$ grid with the wall-clustering parameter used in the paper:

```bash
makeCavityMesh.py 100 0.0975
./Allrun -np 8
```

The mesh-generation utility overwrites `system/blockMeshDict`. Its arguments specify the
number of cells in each direction and the mesh-stretching coefficient (zero for a uniform
grid). A mesh modification requires recomputation of the base state and a consistent domain
decomposition, both performed by `Allrun`. Select the resolution and process count according
to the available memory and computational resources.

### Spectral targeting and Arnoldi parameters

The principal eigenvalue extraction controls in `system/stabilityProperties` include:

| Parameter | Description |
|---|---|
| `sigmaReal`, `sigmaImag` | Real and imaginary components of the spectral shift $\sigma$ |
| `arnoldiNev` | Number of requested eigenvalues |
| `arnoldiNcv` | Dimension of the Arnoldi search subspace |
| `shiftInvertTolerance` | Convergence tolerance for the shifted linear systems |
| `arnoldiConstraintTolerance` | Tolerance associated with monitoring the divergence constraint |

Choose the shift in the vicinity of the spectral region of interest. For example,
`sigmaImag 0.5` targets oscillatory eigenvalues near an angular frequency of 0.5 when the
real shift component is zero. An Arnoldi subspace dimension of at least
`2*arnoldiNev + 1` is recommended as an initial choice, subject to PARPACK's dimensional
requirements and the available memory.

If only the spectral or stability-solver controls are modified, the converged base flow can
be reused. For an existing eight-process decomposition:

```bash
mpirun -np 8 StabilitySolver -parallel > log.StabilitySolver 2>&1
reconstructPar -latestTime > log.reconstructPar 2>&1
touch case.foam
```

The MPI process count must match the existing decomposition. Both solvers write into time
`1`; a continued Newton run (`startFrom latestTime`) starts from it and overwrites it. Eigenvalues returned by shift-and-invert Arnoldi are targeted by
proximity to the shift; establishing the globally leading eigenvalue may require additional
shifts and a sufficiently broad investigation of the relevant spectrum.

### Continuation and restarting

The tutorial script uses continuation in Reynolds number to obtain suitable initial
approximations for Newton iteration. If a continuation step fails, the script halves the
increment and retries, down to an increment of 10 (option `-minstep`). The last converged
Reynolds number is recorded in `lastConvergedRe`.

Provided the existing mesh, decomposed fields and process count are consistent, continuation
may be resumed using:

```bash
./Allrun -np 8 -continue
```

or from a specified converged state already stored on disk:

```bash
./Allrun -np 8 -continue -from 1500
```

The continuation steps can also be prescribed explicitly, e.g.
`./Allrun -np 8 -ramp "500 1000 1500 1800"`, or omitted with `-ramp none`.

---

## 7. Analysis of the computed results

### Base-flow convergence

The workflow records the overall execution in `log.Allrun` and the Newton solver output at
individual continuation points in files of the form `log.NewtonSolver.Re<Re>`. Inspect both
momentum and continuity residuals and the reported convergence status before using a base
state for stability analysis.

### Eigenvalue convergence and critical parameters

The stability solver records its output in `log.StabilitySolver`. The reported growth rate
and angular frequency correspond to the real and imaginary parts of the recovered eigenvalue,
respectively:

- $\lambda_r < 0$: exponential decay of the selected linear mode.
- $\lambda_r > 0$: exponential growth of the selected linear mode.
- $\lambda_r = 0$: neutral stability of that mode.

For a Hopf instability, the critical Reynolds number is identified by the crossing of the
leading eigenvalue through the imaginary axis at nonzero angular frequency. Compute the same
eigenvalue branch at several Reynolds numbers near the crossing and refine the estimate of
the zero growth rate. A negative growth rate for a single targeted mode does not by itself
establish stability of the complete flow.

Assess the shifted-system convergence, the eigensolver convergence, the generalized
eigenpair residual and the continuity constraint when accepting an eigenpair. Inner
linear-solver tolerances should be sufficiently stringent to support the intended eigenvalue
accuracy. For quantitative comparisons, also establish the sensitivity of the computed
eigenpairs to mesh resolution and the spectral extraction parameters.

### Eigenmode visualization

After reconstruction, open the case in ParaView:

```bash
paraview case.foam
```

The fields `leadingEigenModeRealU` and `leadingEigenModeImagU` contain the real and
imaginary parts of the selected velocity eigenmode (the corresponding pressure fields are
`leadingEigenModeRealP` and `leadingEigenModeImagP`). Signed component isosurfaces and
componentwise amplitudes can be used to characterize its spatial structure and reflection
parity.

The physical linear velocity disturbance is proportional to

$$
\mathbf{u}'(\mathbf{x},t) = \operatorname{Re}\left\{ \widehat{\mathbf{u}}(\mathbf{x})\, e^{\lambda t} \right\}.
$$

The complex eigenmode has arbitrary normalization and global phase; comparisons of signed
fields between independent calculations should therefore employ a consistent normalization
and phase convention.

### Paper figures in one click

The folder [`figures/`](figures/README.md) rebuilds Figs. 1–4 of the manuscript in ParaView:
isosurfaces at ±25 % of the maximum of each velocity component, cavity outline, lid arrow,
the same camera in all panels. After `./Allrun` in a case:

| Figure | 256³ (paper case) | 50³ (tutorial) |
|---|---|---|
| Fig. 1, straight cavity, amplitudes | `figures/fig1.sh` / `fig1.bat` | `figures/fig1_50.sh` / `fig1_50.bat` |
| Fig. 2, straight cavity, Re / Im | `figures/fig2.sh` / `fig2.bat` | `figures/fig2_50.sh` / `fig2_50.bat` |
| Fig. 3, diagonal cavity, amplitudes | `figures/fig3.sh` / `fig3.bat` | `figures/fig3_50.sh` / `fig3_50.bat` |
| Fig. 4, diagonal cavity, Re / Im | `figures/fig4.sh` / `fig4.bat` | `figures/fig4_50.sh` / `fig4_50.bat` |

On Windows a double-click on the `.bat` file opens ParaView with the figure; on Linux and
macOS run the `.sh` file. Add `-png` to render the panels without a window.

---

## 8. Troubleshooting

| Symptom | Suggested action |
|---|---|
| OpenFOAM environment initialization hangs or reports `_foamSource: command not found` | Use `scripts/oflsa-shell`; a ParaView version check during environment initialization may be responsible. |
| `Allwmake` reports `Load the OpenFOAM-14 environment first` | Initialize OpenFOAM in the current shell (Section 4). The environment scripts require `bash`; under `tcsh`/`csh` start `bash` first. |
| `Missing mpif90` while building PARPACK | Load the MPI environment module specified by the system administrator (Section 4.2). |
| `exe-not-found` from `mpirun` | Verify that the project environment is initialized and that both solver executables were built and are accessible on the MPI ranks. |
| Newton FGMRES exceeds its iteration limit | Inspect the residual history, solver tolerances and initial base state. Resume from the last converged solution with smaller continuation increments, e.g. `./Allrun -np 8 -continue -ramp "1600 1700 1800"`. |
| PARPACK reports that an Arnoldi basis vector left the solenoidal space | Inspect the base-flow and shifted-system continuity residuals. Test a tighter inner tolerance, e.g. `shiftInvertTolerance 1e-6`. Adjust `arnoldiConstraintTolerance` only if the resulting eigenpair accuracy and continuity error remain acceptable. |
| Linker errors involving `-liomp5` or `-lifport` | Check compiler consistency and use the repository's supplied `Make/options` rather than settings from a different toolchain. |
| PARPACK is not detected | Run `scripts/install_arpack.sh` or install the appropriate ARPACK/PARPACK development packages, ensuring MPI and compiler compatibility. |
| ParaView shows no `leadingEigenMode*` fields | The results were not reconstructed: run `reconstructPar -latestTime` and reopen `case.foam` (time 1). |
| ParaView aborts at start-up (`could not load the Qt platform plugin`, `createPlatformIntegration`) | No usable display for Qt. On a Wayland desktop use `QT_QPA_PLATFORM=xcb paraview case.foam`; over SSH use X forwarding (`ssh -X`) or visualize the reconstructed fields on a separate workstation. |

---

## 9. Container-based execution

A container image including OpenFOAM 14, PARPACK and the solver executables is provided for
systems without a native OpenFOAM installation. The documented execution environment uses
Apptainer (formerly Singularity).

On a single node, the complete workflow of Section 5 — Reynolds-number continuation,
stability analysis and reconstruction — is executed with the container's MPI:

```bash
apptainer pull "$HOME/oflsa.sif" docker://ghcr.io/yurifeld/oflsa:latest
cd ~/myCavity
apptainer exec "$HOME/oflsa.sif" oflsa-run ./Allrun -np 32
```

For multi-node runs the host MPI launches the container on each node, and the individual
steps are invoked explicitly. Newton continuation in Reynolds number (Section 6) remains
necessary:

```bash
SIF="$HOME/oflsa.sif"
apptainer exec "$SIF" oflsa-run foamListTimes -rm      # old results; 0/ = initial conditions
apptainer exec "$SIF" oflsa-run blockMesh
apptainer exec "$SIF" oflsa-run decomposePar -force    # numberOfSubdomains = 128 in system/decomposeParDict
for re in 200 500 1000 1250 1500 1750 2000; do         # ... up to the target Reynolds number
    apptainer exec "$SIF" oflsa-run setRe $re
    mpirun -np 128 apptainer exec "$SIF" oflsa-run NewtonSolver -parallel > log.NewtonSolver.Re$re 2>&1 || break
done
mpirun -np 128 apptainer exec "$SIF" oflsa-run StabilitySolver -parallel > log.StabilitySolver 2>&1
apptainer exec "$SIF" oflsa-run reconstructPar -latestTime
touch case.foam
```

The host MPI implementation and container MPI libraries must be compatible. Follow the local
HPC system's guidance for containerized MPI applications and scheduler integration. For
reproducible research, record the exact container image version or digest rather than
relying solely on the mutable `latest` tag.

---

## 10. Repository structure

| Path | Contents |
|---|---|
| `src/matrixFreeFGMRESSafe/` | Matrix-free flexible GMRES library |
| `applications/solvers/NewtonSolver/` | Steady Newton solver and associated preconditioner |
| `applications/solvers/StabilitySolver/` | Shift-and-invert global stability solver and associated preconditioner |
| `tutorials/lidDrivenCavity3D/` | Reduced-resolution verification case (straight cavity) |
| `tutorials/diagonalCavity3D/` | Reduced-resolution verification case (diagonal cavity) |
| `figures/` | ParaView scripts and one-click launchers for Figs. 1–4 |
| `tutorials/paper/` | Reference configurations from the manuscript |
| `scripts/install_ubuntu.sh` | Installation of all prerequisites and compilation on Ubuntu (`sudo`) |
| `scripts/install_arpack.sh` | Local ARPACK/PARPACK installation (no `sudo`) |
| `scripts/oflsa-shell` | Project environment initialization |
| `scripts/setRe` | Consistent Reynolds-number configuration |
| `scripts/makeCavityMesh.py` | Structured cavity mesh generation |
| `etc/bashrc` | Project environment and library configuration |
| `container/` | Container build configuration |

---

## 11. Citation and license

If the software or reference configurations are used in research, please cite the
accompanying manuscript:

> D. Poliak and Y. Feldman, *Preconditioned Shift-and-Invert Global Linear Stability
> Analysis in OpenFOAM*, under review (2026).

The citation will be updated with the journal reference upon publication.

For reproducibility, report the repository commit or release, OpenFOAM and MPI versions,
mesh characteristics, solver parameters and computational resources used in the study.

The software is distributed under the GNU General Public License, version 3 or later
(GPL-3.0-or-later), consistent with OpenFOAM.
