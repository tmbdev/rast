#!/usr/bin/env bash
# Build driver for the rast project. Invoked by `pixi run build`.
# Uses the vendored colib headers in ./colib.
set -euo pipefail

PYINC=$(python -c 'import sysconfig; print(sysconfig.get_path("include"))')
exec make \
    CXX="$CXX -g -Wall -I. -Icolib -O3 -fPIC" \
    PYINC="$PYINC" \
    "$@"
