#!/bin/sh
# Fig. 2 on the 50^3 grid. Usage: ./fig2_50.sh [<case dir>] [-png]
grid=50; png=
for a in "$@"; do case "$a" in -png) png=-png ;; *) grid="$a" ;; esac; done
exec "$(dirname "$0")/launch.sh" 2 "$grid" $png
