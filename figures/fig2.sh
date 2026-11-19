#!/bin/sh
# Fig. 2 on the 256^3 grid. Usage: ./fig2.sh [<case dir>] [-png]
grid=256; png=
for a in "$@"; do case "$a" in -png) png=-png ;; *) grid="$a" ;; esac; done
exec "$(dirname "$0")/launch.sh" 2 "$grid" $png
