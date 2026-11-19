#!/bin/sh
# Common launcher for fig*.sh:  launch.sh <figure 1-4> <256 | 50 | case dir> [-png]
#   256 / 50   the paper case or the 50^3 tutorial of that cavity
#   -png       no window: render the panels with pvbatch into <case>/figures/
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

# ParaView: PATH, then the usual install places (Linux tarballs, macOS app)
exe=paraview; [ "$png" = "-png" ] && exe=pvbatch
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
if [ "$png" = "-png" ]; then
    OFLSA_PNG=1 exec "$pv" "$here/fig$fig.py"
else
    exec "$pv" --script="$here/fig$fig.py"
fi
