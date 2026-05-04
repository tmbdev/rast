#!/usr/bin/env bash
# Build driver for the rast project. Invoked by `pixi run build*`.
# Three configurations from the same source tree, per cpp_guidelines.md:
#   debug      -- -O0 -g, ASan/UBSan, _GLIBCXX_DEBUG, -Werror, assertions on
#   production -- -O2 -g, assertions on, no sanitizers
#   ultra      -- -O3 -DNDEBUG, assertions stripped, -flto where available
# Each config builds into _build/<config>/. The required warning set is
# applied in all three; -Werror is debug-only so production / ultra still
# build through warnings.
set -euo pipefail

BUILD_TYPE="${BUILD_TYPE:-debug}"
case "$BUILD_TYPE" in
debug | production | ultra) ;;
*)
    echo "BUILD_TYPE must be debug, production, or ultra (got: $BUILD_TYPE)" >&2
    exit 1
    ;;
esac

WARNINGS="-Wall -Wextra -Wpedantic -Wshadow -Wnon-virtual-dtor"
WARNINGS="$WARNINGS -Wold-style-cast -Wcast-align -Woverloaded-virtual -Wconversion"
# -Wdouble-promotion is intentionally omitted pending a separate cleanup
# sweep: the matchers do most arithmetic in float but use double literals
# (0.0 / 1.0); each std::max(0.0,...) and `q == 0.0` triggers it. M_PI uses
# already migrated to the kPi family in util.h, but the literal-suffix sweep
# is still pending.

COMMON="-std=c++20 $WARNINGS -fPIC -fvisibility=hidden -fvisibility-inlines-hidden"

AR_OVERRIDE="ar"
case "$BUILD_TYPE" in
debug)
    CONFIG_FLAGS="-O0 -g -fsanitize=address,undefined -D_GLIBCXX_DEBUG -Werror"
    ;;
production)
    CONFIG_FLAGS="-O2 -g"
    ;;
ultra)
    # -flto produces LTO bitcode in each .o; librast.a must therefore be
    # built with gcc-ar so the archive index covers the bitcode. The linker
    # then resolves make*() factories correctly.
    CONFIG_FLAGS="-O3 -DNDEBUG -flto"
    AR_OVERRIDE="gcc-ar"
    ;;
esac

mkdir -p "_build/$BUILD_TYPE"
PYINC=$(python -c 'import sysconfig; print(sysconfig.get_path("include"))')
cd "_build/$BUILD_TYPE"
exec make -f ../../Makefile \
    VPATH=../../src:../../tests:../../bindings/python \
    SRCDIR=../../src \
    CXX="$CXX $COMMON $CONFIG_FLAGS -I../../src" \
    AR="$AR_OVERRIDE" \
    PYINC="$PYINC" \
    "$@"
