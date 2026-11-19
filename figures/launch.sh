#!/bin/sh
# Common launcher for fig*.sh:  launch.sh <figure 1-4> <256 | 50 | case dir> [-png]
#   256 / 50   the paper case or the 50^3 tutorial of that cavity
#   -png       also save every panel as PNG into <case>/figures/
here=$(cd "$(dirname "$0")" && pwd)
root=$(dirname "$here")
fig="$1"; grid="${2:-256}"; png="$3"

case "$fig:$grid" in
    [12]:256) caseDir="$root/tutorials/paper/straightCavity256" ;;
    [12]:50)  caseDir="$root/tutorials/lidDrivenCavity3D" ;;
    [34]:256) caseDir="$root/tutorials/paper/diagonalCavity256" ;;
    [34]:50)  caseDir="$root/tutorials/diagonalCavity3D" ;;
    *)        caseDir=$(cd "$grid" 2>/dev/null && pwd) \
                  || { echo "No such case directory: $grid"; exit 1; } ;;
esac

if [ ! -f "$caseDir/1/leadingEigenModeRealU" ] \
   && [ ! -f "$caseDir/processor0/1/leadingEigenModeRealU" ]; then
    echo "No results in $caseDir/1 - run ./Allrun in that case first."
    exit 1
fi

# ParaView: PATH, then the usual install places (Linux tarballs, macOS app).
# PNGs are also made with the GUI: pvbatch often cannot render without a display.
exe=paraview
pv=$(command -v "$exe" 2>/dev/null)
if [ -z "$pv" ]; then
    for d in /opt/ParaView*/bin "$HOME"/ParaView*/bin "$HOME"/opt/ParaView*/bin \
             /Applications/ParaView*.app/Contents/bin /Applications/ParaView*.app/Contents/MacOS; do
        [ -x "$d/$exe" ] && pv="$d/$exe"
    done
fi
[ -n "$pv" ] || { echo "$exe not found: install ParaView or put its bin/ in PATH"; exit 1; }

export OFLSA_CASE="$caseDir" OFLSA_FIGDIR="$here"
echo "Fig. $fig from $caseDir"
[ "$png" = "-png" ] && export OFLSA_PNG=1 && echo "PNGs -> $caseDir/figures/"
# Detach ParaView from the terminal: ParaView 6 started from a terminal waits
# for it (stdin) and does not show its window until the terminal is closed.
log="$caseDir/log.figure$fig"
if command -v setsid >/dev/null 2>&1; then
    setsid "$pv" --script="$here/fig$fig.py" < /dev/null > "$log" 2>&1 &
else
    nohup "$pv" --script="$here/fig$fig.py" < /dev/null > "$log" 2>&1 &
fi
echo "ParaView started (messages: $log)"
