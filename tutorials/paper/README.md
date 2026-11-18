# Cases from the paper

Ready-to-run setups that reproduce Tables 1 and 2 of

> D. Polyak and Y. Feldman, *Preconditioned Shift-and-Invert Global Linear Stability
> Analysis in OpenFOAM*, under review (2026).

| Case | Lid | Re | Shift σ | Paper result |
|---|---|---|---|---|
| `straightCavity256` | along x, \|U\| = 1 | 1928.75 | 0.5 i | λ = −2.81·10⁻⁶ + 0.584058 i |
| `diagonalCavity256` | along the horizontal diagonal, \|U\| = 1 | 2320.8 | 0.22 i | ω = 0.2507, near neutral |

Both use the 256³ grid with wall clustering `x = ξ − 0.0975 sin(2πξ)` in all directions
(Gelfgat, TCFD 2019), γ = 0.02 in the mAL preconditioner.

**These are cluster runs.** 16.7 M cells; the paper used 256 cores: Newton up to ~2 h,
shift-invert Arnoldi 14–30 h. Run on a node or through the scheduler, e.g.

```bash
./Allrun -np 256
```

Through a scheduler, a job script loads the environment (README section 2b) and calls
`Allrun`. SLURM:
```bash
#!/bin/bash
#SBATCH --job-name=straight256
#SBATCH --ntasks=256
#SBATCH --time=48:00:00
module load openfoam/14                  # or: source /path/to/OpenFOAM-14/etc/bashrc
source ~/ofLSA/etc/bashrc
cd ~/ofLSA/tutorials/paper/straightCavity256
./Allrun -np $SLURM_NTASKS > log.Allrun 2>&1
```
SGE (the parallel environment name, here `mpi`, differs between clusters — ask your administrator):
```bash
#!/bin/bash
#$ -N straight256
#$ -pe mpi 256
#$ -l h_rt=48:00:00
#$ -S /bin/bash
#$ -cwd
module load openfoam/14                  # or: source /path/to/OpenFOAM-14/etc/bashrc
source ~/ofLSA/etc/bashrc
./Allrun -np $NSLOTS > log.Allrun 2>&1
```

`Allrun` reaches the target Re by continuation (200, 500, 1000, then every 250), halving the
step if Newton does not converge. Each Newton step on this grid is expensive, so if you already
have a converged base flow at a nearby Re, copy it into the case and start from it instead:

```bash
./Allrun -np 256 -continue -from <Re of that flow>
```

To try the same physics on a workstation, coarsen the grid inside the case:
`makeCavityMesh.py 100 0.0975` (1 M cells). Results then differ from the paper by the grid error.

Shift σ, Arnoldi settings and the diagonal lid velocity (1/√2, 1/√2, 0) are taken from the
authors' original run directories (`Stability_LDC_grid_256`, `Stability_LDC_grid_256_45_deg`).
