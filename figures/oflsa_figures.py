# -*- coding: utf-8 -*-
"""
Builders for the eigenmode figures of the paper (Figs. 1-4), for ParaView.

Not run directly: the small scripts fig1.py ... fig4.py call build_figure(),
and the launchers fig*.sh / fig*.bat start ParaView with them.

Environment (set by the launchers):
  OFLSA_CASE     case directory with case results in time 1/ (required)
  OFLSA_FIGDIR   this folder (to find the module)
  OFLSA_PNG      "1": also save PNGs (used with pvbatch, see README.md)

Everything a user may want to tune is in the SETTINGS block below.
"""

import math
import os

from paraview.simple import *  # noqa: F401,F403

# ----------------------------------------------------------------- SETTINGS --
# Isosurface level as a fraction of max|component|.  0.25 = 25 % of the max.
LEVEL_FRACTION = 0.25
# True:  one level per component, from max over the real AND imaginary parts
#        (as in the captions of Figs. 2 and 4 of the paper).
# False: real and imaginary parts get their own 25 % levels.
SHARED_RE_IM_LEVEL = False
# Round levels to this many significant digits (nicer captions: 2.63e-4).
LEVEL_DIGITS = 3

POSITIVE = [0.80, 0.08, 0.08]     # red
NEGATIVE = [0.25, 0.45, 0.85]     # blue
AMPLITUDE = [0.00, 0.55, 0.60]    # teal, Figs. 1 and 3
OPACITY = 1.0

# Camera direction (from the focal point towards the eye), like the paper:
# x along the front edge, y going back to the right, z up.
CAMERA_DIRECTION = (-0.55, -1.0, 0.70)
PANEL_PIXELS = (900, 900)         # PNG size of one panel
# ------------------------------------------------------------------------------

REAL_U = "leadingEigenModeRealU"
IMAG_U = "leadingEigenModeImagU"


# --- components ------------------------------------------------------------------
def _component(field, comp):
    """Calculator expression of one velocity component of `field`."""
    x, y, z = field + "_X", field + "_Y", field + "_Z"
    return {
        "x": x,
        "y": y,
        "z": z,
        # Eq. (40) of the paper: along / across the lid diagonal
        "par": "(%s+%s)/sqrt(2)" % (x, y),
        "perp": "(-%s+%s)/sqrt(2)" % (x, y),
    }[comp]


LABEL = {"x": "u_x", "y": "u_y", "z": "u_z", "par": "u_par", "perp": "u_perp"}

# Paper figures: kind, lid direction (deg from x), components (one per row)
FIGURES = {
    1: dict(kind="amplitude", lid_deg=0.0, comps=["x", "y", "z"],
            title="Fig. 1 - straight cavity, mode amplitudes"),
    2: dict(kind="signed", lid_deg=0.0, comps=["x", "y", "z"],
            title="Fig. 2 - straight cavity, real / imaginary parts"),
    3: dict(kind="amplitude", lid_deg=45.0, comps=["par", "perp", "z"],
            title="Fig. 3 - diagonal cavity, mode amplitudes"),
    4: dict(kind="signed", lid_deg=45.0, comps=["par", "perp", "z"],
            title="Fig. 4 - diagonal cavity, real / imaginary parts"),
}


# --- helpers ---------------------------------------------------------------------
def _round_sig(v, digits):
    if v == 0:
        return 0.0
    return round(v, digits - 1 - int(math.floor(math.log10(abs(v)))))


def _try_set(obj, name, value):
    try:
        setattr(obj, name, value)
        return True
    except Exception:
        return False


def _find_case():
    case = os.environ.get("OFLSA_CASE", "")
    if not case:
        raise RuntimeError(
            "OFLSA_CASE is not set: start the figure with figures/fig*.sh or "
            "figures/fig*.bat, or set OFLSA_CASE to the case directory.")
    case = os.path.abspath(case)
    reconstructed = os.path.isfile(os.path.join(case, "1", REAL_U))
    decomposed = os.path.isfile(os.path.join(case, "processor0", "1", REAL_U))
    if not (reconstructed or decomposed):
        raise RuntimeError(
            "No eigenmode in %s/1 (or processor*/1). Run the case first "
            "(./Allrun) - the results are written into time 1." % case)
    return case, reconstructed


def _open_case(case, reconstructed):
    foam = os.path.join(case, "case.foam")
    if not os.path.exists(foam):
        open(foam, "a").close()
    reader = OpenFOAMReader(registrationName="case.foam", FileName=foam)
    _try_set(reader, "CaseType",
             "Reconstructed Case" if reconstructed else "Decomposed Case")
    reader.MeshRegions = ["internalMesh"]
    reader.CellArrays = [REAL_U, IMAG_U]
    reader.UpdatePipelineInformation()
    times = reader.TimestepValues
    if not isinstance(times, (list, tuple)):
        times = [times]
    t = max(times) if len(times) else 0.0
    reader.UpdatePipeline(t)
    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()
    scene.AnimationTime = t
    points = CellDatatoPointData(registrationName="pointData", Input=reader)
    return reader, points, t


def _calc(points, name, expr, t):
    c = Calculator(registrationName=name, Input=points)
    c.ResultArrayName = name
    c.Function = expr
    c.UpdatePipeline(t)
    lo, hi = c.PointData[name].GetRange()
    return c, max(abs(lo), abs(hi))


def _iso(calc, name, level, color, view, t):
    if level <= 0:
        return
    s = Contour(registrationName="iso_%s_%+.3g" % (name, level), Input=calc)
    s.ContourBy = ["POINTS", name]
    s.Isosurfaces = [level]
    _try_set(s, "ComputeScalars", 0)
    s.UpdatePipeline(t)
    d = Show(s, view)
    ColorBy(d, None)
    d.AmbientColor = color
    d.DiffuseColor = color
    d.Opacity = OPACITY
    _try_set(d, "Specular", 0.3)


def _new_view():
    v = CreateView("RenderView")
    v.OrientationAxesVisibility = 0
    _try_set(v, "UseColorPaletteForBackground", 0)
    v.Background = [1, 1, 1]
    ag = v.AxesGrid
    ag.Visibility = 1
    ag.XTitle, ag.YTitle, ag.ZTitle = "x", "y", "z"
    for a in ("XTitleColor", "YTitleColor", "ZTitleColor",
              "XLabelColor", "YLabelColor", "ZLabelColor"):
        _try_set(ag, a, [0, 0, 0])
    _try_set(ag, "ShowGrid", 1)
    _try_set(ag, "GridColor", [0.85, 0.85, 0.85])
    return v


def _lid_sources(reader, bounds, lid_deg):
    """Cavity outline and lid arrow, built once and shown in every panel."""
    outline = Outline(registrationName="outline", Input=reader)
    x0, x1, y0, y1, z0, z1 = bounds
    L = 0.4 * (x1 - x0)
    a = math.radians(lid_deg)
    cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
    arrow = Arrow(registrationName="lidArrowShape")
    arrow.TipLength, arrow.TipRadius, arrow.ShaftRadius = 0.25, 0.07, 0.022
    tr = Transform(registrationName="lidArrow", Input=arrow)
    _try_set(tr, "Transform", "Transform")
    tr.Transform.Translate = [cx - 0.5 * L * math.cos(a),
                              cy - 0.5 * L * math.sin(a),
                              z1 + 0.02 * (z1 - z0)]
    tr.Transform.Rotate = [0.0, 0.0, lid_deg]
    tr.Transform.Scale = [L, L, L]
    return outline, tr


def _decorate(view, sources, bounds, label, name):
    outline, arrow = sources
    o = Show(outline, view)
    ColorBy(o, None)
    o.AmbientColor = o.DiffuseColor = [0, 0, 0]
    _try_set(o, "LineWidth", 1.5)
    d = Show(arrow, view)
    ColorBy(d, None)
    d.AmbientColor = d.DiffuseColor = [0, 0, 0]
    x0, x1, y0, y1, z0, z1 = bounds
    cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)

    # panel label (component, part, isolevel)
    txt = Text(registrationName="label_" + name)
    txt.Text = label
    td = Show(txt, view)
    if not _try_set(td, "WindowLocation", "Upper Left Corner"):
        _try_set(td, "WindowLocation", "UpperLeftCorner")
    _try_set(td, "Color", [0, 0, 0])
    _try_set(td, "FontSize", 14)

    # camera
    c = [cx, cy, 0.5 * (z0 + z1)]
    dx, dy, dz = CAMERA_DIRECTION
    view.CameraFocalPoint = c
    view.CameraPosition = [c[0] + 3 * dx, c[1] + 3 * dy, c[2] + 3 * dz]
    view.CameraViewUp = [0, 0, 1]
    ResetCamera(view)


def _layout(name, rows, cols):
    """Layout with rows x cols panels; returns the cell indices row by row."""
    lay = CreateLayout(name)
    if cols == 1:
        cells = [0]
    else:                                  # two columns: left / right
        lay.SplitHorizontal(0, 0.5)
        cells = [1, 2]
    columns = []
    for c in cells:
        if rows == 1:
            columns.append([c])
        else:                              # three rows
            lay.SplitVertical(c, 1.0 / 3.0)
            top, rest = 2 * c + 1, 2 * c + 2
            lay.SplitVertical(rest, 0.5)
            columns.append([top, 2 * rest + 1, 2 * rest + 2])
    return lay, [[columns[j][i] for j in range(cols)] for i in range(rows)]


def _layout_row(name):
    """One row of three panels (Figs. 1 and 3)."""
    lay = CreateLayout(name)
    lay.SplitHorizontal(0, 1.0 / 3.0)
    lay.SplitHorizontal(2, 0.5)
    return lay, [1, 5, 6]


# --- main entry --------------------------------------------------------------------
def build_figure(number):
    spec = FIGURES[number]
    case, reconstructed = _find_case()
    reader, points, t = _open_case(case, reconstructed)
    bounds = reader.GetDataInformation().GetBounds()
    grid = os.path.basename(os.path.normpath(case))
    sources = _lid_sources(reader, bounds, spec["lid_deg"])
    lines = ["%s" % spec["title"], "case: %s, time %g" % (case, t),
             "level = %g x max|component|%s" % (
                 LEVEL_FRACTION,
                 ", shared by Re/Im" if SHARED_RE_IM_LEVEL else "")]

    views = []
    if spec["kind"] == "amplitude":
        lay, cells = _layout_row("Fig%d %s" % (number, grid))
        for comp, cell in zip(spec["comps"], cells):
            re_ = _component(REAL_U, comp)
            im_ = _component(IMAG_U, comp)
            name = "amp_" + comp
            calc, amax = _calc(points, name,
                               "sqrt((%s)^2+(%s)^2)" % (re_, im_), t)
            level = _round_sig(LEVEL_FRACTION * amax, LEVEL_DIGITS)
            v = _new_view()
            AssignViewToLayout(view=v, layout=lay, hint=cell)
            _iso(calc, name, level, AMPLITUDE, v, t)
            label = "|%s| = %.*e" % (LABEL[comp], LEVEL_DIGITS - 1, level)
            _decorate(v, sources, bounds, label, name)
            lines.append(label)
            views.append((v, name))
    else:
        lay, grid_cells = _layout("Fig%d %s" % (number, grid), 3, 2)
        for comp, row in zip(spec["comps"], grid_cells):
            parts = []
            for part, field in (("Re", REAL_U), ("Im", IMAG_U)):
                name = "%s_%s" % (part, comp)
                calc, amax = _calc(points, name, _component(field, comp), t)
                parts.append((part, name, calc, amax))
            if SHARED_RE_IM_LEVEL:
                m = max(p[3] for p in parts)
                parts = [(p[0], p[1], p[2], m) for p in parts]
            for (part, name, calc, amax), cell in zip(parts, row):
                level = _round_sig(LEVEL_FRACTION * amax, LEVEL_DIGITS)
                v = _new_view()
                AssignViewToLayout(view=v, layout=lay, hint=cell)
                _iso(calc, name, +level, POSITIVE, v, t)
                _iso(calc, name, -level, NEGATIVE, v, t)
                label = "%s(%s) = +-%.*e" % (part, LABEL[comp],
                                             LEVEL_DIGITS - 1, level)
                _decorate(v, sources, bounds, label, name)
                lines.append(label)
                views.append((v, name))

    # rotate one panel -> all panels follow
    for i, (v, _) in enumerate(views[1:]):
        AddCameraLink(views[0][0], v, "oflsaCam%d_%d" % (number, i))
    SetActiveView(views[0][0])
    for v, _ in views:
        Render(v)

    # isolevels for the caption + a .pvsm state next to the case
    tag = "fig%d" % number
    with open(os.path.join(case, tag + "_levels.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    try:
        SaveState(os.path.join(case, tag + ".pvsm"))
    except Exception as e:  # not fatal
        print("Could not save the state: %s" % e)

    if os.environ.get("OFLSA_PNG") == "1":
        outdir = os.path.join(case, "figures")
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        for v, name in views:
            SaveScreenshot(os.path.join(outdir, "%s_%s.png" % (tag, name)), v,
                           ImageResolution=list(PANEL_PIXELS))
        print("PNGs written to %s" % outdir)
