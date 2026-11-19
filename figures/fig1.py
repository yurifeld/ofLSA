# Fig. 1 of the paper: straight cavity, mode amplitudes.
# Start it with fig1.sh (Linux/macOS) or fig1.bat (Windows): they find ParaView
# and the case. Manually: OFLSA_CASE=<case> paraview --script=fig1.py
import os
import sys

_here = os.environ.get("OFLSA_FIGDIR") or os.path.dirname(
    os.path.abspath(globals().get("__file__", "fig1.py")))
if _here not in sys.path:
    sys.path.insert(0, _here)

import oflsa_figures  # noqa: E402

oflsa_figures.build_figure(1)
