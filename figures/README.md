# Paper figures in ParaView

Each launcher opens ParaView with one figure of the manuscript, built from the computed
eigenmode: one panel per velocity component, all cameras linked (rotate one panel and the
others follow).

| Figure | Content | 256³ (paper case) | 50³ (tutorial) |
|---|---|---|---|
| 1 | straight cavity, amplitudes \|û_x\|, \|û_y\|, \|û_z\| | `fig1.sh`, `fig1.bat` | `fig1_50.sh`, `fig1_50.bat` |
| 2 | straight cavity, Re / Im of û_x, û_y, û_z | `fig2.sh`, `fig2.bat` | `fig2_50.sh`, `fig2_50.bat` |
| 3 | diagonal cavity, amplitudes \|û_∥\|, \|û_⊥\|, \|û_z\| | `fig3.sh`, `fig3.bat` | `fig3_50.sh`, `fig3_50.bat` |
| 4 | diagonal cavity, Re / Im of û_∥, û_⊥, û_z | `fig4.sh`, `fig4.bat` | `fig4_50.sh`, `fig4_50.bat` |

The cases are `tutorials/paper/straightCavity256`, `tutorials/paper/diagonalCavity256`
(256³) and `tutorials/lidDrivenCavity3D`, `tutorials/diagonalCavity3D` (50³).

## Use

1. Run the case (`./Allrun`). The results are written into time `1/`; `0/` holds only the
   initial conditions. A decomposed case (`processor*/1/`) also works.
2. Start the figure:
   - **Windows:** double-click `fig2.bat` (or `fig2_50.bat`). ParaView is found in `PATH`
     or in `C:\Program Files\ParaView*`. To view results computed on the cluster, copy the
     whole case folder (with `1\`, `constant\polyMesh\`, `system\`) into the same place in
     the repository, e.g. with WinSCP.
   - **Linux / macOS:** `./fig2.sh` (or `./fig2_50.sh`).
   - Any other case: `./fig2.sh /path/to/case` or `fig2.bat C:\path\to\case`.
3. Without a window (cluster, batch): `./fig2.sh -png` renders every panel with `pvbatch`
   into `<case>/figures/`.

Every run also writes into the case directory:

- `figN_levels.txt` — the isosurface levels of each panel, for the figure caption;
- `figN.pvsm` — a regular ParaView state (*File → Load State*) for this machine.

## What is drawn

- Signed parts (Figs. 2, 4): isosurfaces at **+25 %** (red) and **−25 %** (blue) of
  max|component|, real and imaginary parts each with their own level.
- Amplitudes (Figs. 1, 3): one isosurface at 25 % of the maximum amplitude
  √(Re² + Im²) of the component.
- Diagonal cavity: û_∥ = (û_x + û_y)/√2, û_⊥ = (−û_x + û_y)/√2 (Eq. 40 of the paper).
- Cavity outline, an arrow along the lid velocity on the top face, axes grid.

Levels are rounded to three significant digits.

## Tuning

All settings are at the top of `oflsa_figures.py`: `LEVEL_FRACTION` (0.25),
`SHARED_RE_IM_LEVEL` (one level for Re and Im, as in the paper captions), colours,
camera direction, PNG size. A view adjusted by hand in ParaView can be kept with
*File → Save State*.
