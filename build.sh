#!/usr/bin/env bash
# Build driver for the rast project. Invoked by `pixi run build`.
# Compiles into ./_build/, using the vendored colib headers in src/colib.
set -euo pipefail

mkdir -p _build
PYINC=$(python -c 'import sysconfig; print(sysconfig.get_path("include"))')
cd _build
exec make -f ../Makefile \
    VPATH=../src:../tests:../bindings/python \
    SRCDIR=../src \
    CXX="$CXX -g -Wall -I../src -O3 -fPIC" \
    PYINC="$PYINC" \
    "$@"
