#!/bin/bash
# Build 2D and 3D as C++20 and fail on a C++20 deprecation in first-party code.
# usage: cxx20-check.sh <make args...>   (run from the repository root)
#
# It changes the standard in the Makefile itself, so it edits the checkout:
# run it only in CI or a scratch copy.
set -euo pipefail

n=$(grep -c -- '-std=c++17' Makefile || true)
test "$n" -gt 0 || { echo "no -std=c++17 in the Makefile: update this script"; exit 1; }
# -Wall leaves out the deprecated-copy warnings, so they are added here.
perl -pi -e 's/-std=c\+\+17/-std=c++20 -Wdeprecated -Wdeprecated-copy -Wdeprecated-copy-dtor/g' Makefile
if grep -q -- '-std=c++17' Makefile; then echo "-std=c++17 left in the Makefile"; exit 1; fi
echo "Makefile: $n line(s) now ask for -std=c++20"

# Vendored code carries its own license and is not fixed here.
vendored='^(/|(\./)?(3x3-C|tetgen|triangle|nanoflann|knn-bvh|mmg|third_party)/)'
for nd in 2 3; do
    log=build-cxx20-${nd}d.log
    make -j4 ndims=$nd "$@" 2>&1 | tee $log
    # proves the patched flag reached the compile lines
    grep -q -- '-std=c++20' $log || { echo "no -std=c++20 in $log"; exit 1; }
    grep -E '^[^ :]+\.(cxx|hpp|cpp|c|h):[0-9]+:[0-9]+: warning:' $log \
        | grep -v -E "$vendored" \
        | grep -i -E 'deprecat|c\+\+20|reversed' | sort -u > cxx20-${nd}d.txt || true
    if [ -s cxx20-${nd}d.txt ]; then
        echo "C++20 deprecations in first-party code (${nd}D):"
        cat cxx20-${nd}d.txt
        exit 1
    fi
    echo "${nd}D: builds as C++20, no first-party deprecation warnings"
done
