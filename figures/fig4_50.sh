#!/bin/sh
# Fig. 4 on the 50^3 grid. Usage: ./fig4_50.sh [<case dir>] [-png]
grid=50; png=
for a in "$@"; do case "$a" in -png) png=-png ;; *) grid="$a" ;; esac; done
exec "$(dirname "$0")/launch.sh" 4 "$grid" $png
