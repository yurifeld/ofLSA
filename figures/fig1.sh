#!/bin/sh
# Fig. 1 on the 256^3 grid. Usage: ./fig1.sh [<case dir>] [-png]
grid=256; png=
for a in "$@"; do case "$a" in -png) png=-png ;; *) grid="$a" ;; esac; done
exec "$(dirname "$0")/launch.sh" 1 "$grid" $png
