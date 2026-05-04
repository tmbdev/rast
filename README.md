# rast

A C++ library implementing the **RAST** family of branch-and-bound feature
matching algorithms (Recognition by Adaptive Subdivisions of Transformation
space). Given a set of model features and a set of image features, RAST finds
globally optimal correspondences and the rigid 2D transform between them by
recursively subdividing the parameter space and bounding the match quality
in each subregion.

The repository ships:

- a static C++ library `librast.a`,
- a command-line driver `rast` exposing the algorithms as subcommands,
- a self-contained Canny-style edge detector `cedges` (file: `src/cedges.cc`),
- a randomized regression test binary `rast-test`,
- a doctest-based unit-test binary `tests`,
- Python bindings via pybind11 (`rast.so`).

## Public interfaces

All defined in [src/rast.h](src/rast.h). Each interface comes with a `make*()`
factory that returns a heap-allocated implementation; clients hold the result
in a smart pointer.

| Interface         | Purpose                                                          |
| ----------------- | ---------------------------------------------------------------- |
| `AlignmentP2D`    | Brute-force 2-point alignment of point sets                      |
| `InstanceP2D`     | Synthetic instance generator for the 2D recognition problem      |
| `LinesP2D`        | Line finding from oriented points (RAST + Hough-style)           |
| `LinesS2D`        | Line finding from line segments                                  |
| `RastP2D`         | Point-to-point matching, branch-and-bound over (tx, ty, θ, s)    |
| `RastS2D`         | Segment-to-segment matching                                      |
| `RastSS2D`        | Segment-to-segment matching with sampling                        |
| `RastRS2D`        | Line+blob segment matching with sampling                         |

The standalone edge detector `EdgeDetector` is in [src/cedges.h](src/cedges.h).

## Build

The build is driven by [pixi](https://pixi.sh/), which provisions the entire
toolchain (C++ compiler, pybind11, Python, NumPy, doctest, clang-format) from
conda-forge into a project-local environment.

```sh
# one-time: install pixi (Rust-based, single binary)
curl -fsSL https://pixi.sh/install.sh | bash

# in this repo
pixi install        # downloads and pins the toolchain into .pixi/
pixi run build      # builds everything into ./_build/
pixi run test       # runs the doctest suite
pixi run clean      # removes ./_build/
```

All build artifacts (`.o` files, `librast.a`, `rast`, `rast-test`, `cedges`,
`tests`, `rast.so`) land under `_build/`. The project root stays
clean.

The C++ standard library is the only runtime dependency; the previously-vendored
`iulib`/`colib` code has been removed and all of its container, smart-pointer,
heap, and I/O abstractions replaced with their std equivalents.

### What pixi runs under the hood

- [build.sh](build.sh) creates `_build/`, then runs `make -f ../Makefile` from
  inside it with `VPATH` and include paths shifted by `..`. This keeps the
  project root free of generated files.
- [Makefile](Makefile) drives the actual compilation. It is a plain GNU Make
  file (the project's chosen build system per the coding guidelines).
- The conda-forge compiler's flags are auto-augmented via `$CXX`; build.sh
  appends `-I../src` so headers in `src/` are found.

### Python module

After `pixi run build`, the pybind11 extension module is at `_build/rast.so`.
The pixi environment sets `PYTHONPATH=$PIXI_PROJECT_ROOT/_build` so:

```sh
pixi run python -c "import rast; print(rast.makeRastP2D())"
```

works without any extra setup.

## Layout

```text
rast/
├── src/                  # C++ sources and headers
│   ├── rast.h            # public umbrella header (all factories)
│   ├── cedges.h          # edge-detector public header
│   ├── vec2.h            # local 2D vector type (no std equivalent)
│   ├── heap.h            # priority-queue wrapper used by the matchers
│   ├── util.h            # math helpers + env-var parameter readers
│   ├── trie.h            # spatial bucket index used by AlignmentP2D
│   ├── *.cc              # implementation files (one per impl class)
│   └── *_test.cc         # doctest unit tests (alongside the impls)
├── tests/
│   └── rast-test.cc      # randomized regression-test driver
├── bindings/python/
│   ├── rast_pybind.cc    # pybind11 module
│   ├── rast-test.py      # Python smoke test
│   └── setup.py          # legacy distutils setup (kept for reference)
├── examples/
│   └── RastP2D-Demo.ipynb
├── doc/                  # Texinfo documentation
├── Makefile              # canonical build rules
├── build.sh              # build driver invoked by pixi
├── pixi.toml             # pixi environment + tasks
├── pixi.lock             # pixi lockfile
├── .clang-format         # formatter config
├── cpp_guidelines.md     # project coding conventions
├── LICENSE               # Apache 2.0
└── README.md             # this file
```

## CLI quickstart

The `rast` binary exposes the algorithms as subcommands. Parameters are passed
through the environment, in the lowercase-env-var style described in the coding
guidelines:

```sh
# generate a synthetic instance
nclutter=20 nmodel_total=15 nmodel_unoccluded=10 \
    pixi run -- _build/rast instance model.points image.points

# match it back
maxresults=1 minscale=0.95 maxscale=1.05 \
    pixi run -- _build/rast rast model.points image.points

# find lines through a set of oriented points
pixi run -- _build/rast lines data.points
```

`pixi run -- _build/rast` with no arguments prints the full subcommand list
and parameter inventory.

## License

Apache License 2.0; see [LICENSE](LICENSE). All source files carry a
two-line header pointing at it.

Copyright 1990-2026 Thomas M. Breuel.
