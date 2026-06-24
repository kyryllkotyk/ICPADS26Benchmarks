#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$SCRIPT_DIR/../FinalSource"
OUTPUT="$SCRIPT_DIR/pm2heat3d"

CXX="${CXX:-}"

if [ -z "$CXX" ]; then
    if command -v mpicxx.madmpi >/dev/null 2>&1; then
        CXX="mpicxx.madmpi"
    elif command -v mpic++.madmpi >/dev/null 2>&1; then
        CXX="mpic++.madmpi"
    elif command -v mpicxx >/dev/null 2>&1; then
        CXX="mpicxx"
    else
        CXX="g++"
    fi
fi

"$CXX" -O3 -flto=8 \
    "$SOURCE_DIR/mainPM2.cpp" \
    "$SOURCE_DIR/pm2heat3d.cpp" \
    -o "$OUTPUT"

echo "Built $OUTPUT"
