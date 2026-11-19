# Fig. 4 of the paper: diagonal cavity, real/imaginary parts.
# Start it with fig4.sh (Linux/macOS) or fig4.bat (Windows): they find ParaView
# and the case. Manually: OFLSA_CASE=<case> paraview --script=fig4.py
import os
import sys

_here = os.environ.get("OFLSA_FIGDIR") or os.path.dirname(
    os.path.abspath(globals().get("__file__", "fig4.py")))
if _here not in sys.path:
    sys.path.insert(0, _here)

import oflsa_figures  # noqa: E402

oflsa_figures.build_figure(4)
