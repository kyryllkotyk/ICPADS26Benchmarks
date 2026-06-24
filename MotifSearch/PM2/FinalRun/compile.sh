#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$SCRIPT_DIR/../FinalSource"

MOTIF_OUTPUT="$SCRIPT_DIR/pm2motifsearch"
GENERATOR_SOURCE="$SCRIPT_DIR/generateMotifGraph.cpp"
GENERATOR_OUTPUT="$SCRIPT_DIR/generateMotifGraph"

EDGE_FILE="$SCRIPT_DIR/triangle250edgelist.edgelist"
VERTICES=250000
EDGE_PROB=0.0011
SEED=42

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

echo "Building MotifSearch executable"

"$CXX" -O3 -flto=8 \
    "$SOURCE_DIR/main.cpp" \
    "$SOURCE_DIR/pm2motifsearch.cpp" \
    -o "$MOTIF_OUTPUT"

echo "Built $MOTIF_OUTPUT"

if [ ! -f "$GENERATOR_SOURCE" ]; then
    echo "ERROR: missing graph generator: $GENERATOR_SOURCE" >&2
    exit 1
fi

echo "Building graph generator"

g++ -O3 -flto=8 \
    "$GENERATOR_SOURCE" \
    -o "$GENERATOR_OUTPUT"

echo "Built $GENERATOR_OUTPUT"

if [ -f "$EDGE_FILE" ]; then
    echo "Edge file already exists: $EDGE_FILE"
    echo "Skipping regeneration"
else
    echo "Generating edge file: $EDGE_FILE"

    "$GENERATOR_OUTPUT" \
        --vertices "$VERTICES" \
        --edge-prob "$EDGE_PROB" \
        --seed "$SEED" \
        --output "$EDGE_FILE"

    echo "Generated $EDGE_FILE"
fi
