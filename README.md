# TalisMan 2.0

TalisMan Version 2.0
is a verification and analysis tool for AIG-based circuits.
It combines algebraic reasoning and SAT-based techniques using
CaDiCaL, FLINT, PBLib, and cmsgen.
PBLib and cmsgen are included in this repository,
CaDiCaL and FLINT have to be provided (see below).

In the corresponding IJCAR'26 paper 
C. Hofstadler, D.Kaufmann, and C. Chen
Avoiding Big Integers: Parallel Multimodular Algebraic Verification of Arithmetic Circuits

we have used the version of TalisMan tagged ijcar26
and the following versions of the dependencies
CaDiCaL 2.2.0
FLINT 3.4.0
PBLib 1.2.1
cmsgen 1.0



------------------------------------------------------------
Prerequisites
------------------------------------------------------------

You need:

- CMake >= 3.16
- GCC or Clang (C++20)
- GMP (including its C++ interface gmpxx)
- FLINT 3
- CaDiCaL (built from source, see below)


Linux (Ubuntu):

    sudo apt install \
    build-essential \
    cmake \
    libgmp-dev \
    libflint-dev


macOS:

    brew install cmake git gmp flint


CaDiCaL:

    git clone https://github.com/arminbiere/cadical.git
    cd cadical && ./configure && make -j && cd ..


------------------------------------------------------------
Building
------------------------------------------------------------

Build cmsgen (included in includes/cmsgen):

    cmake -S includes/cmsgen -B build/cmsgen \
      -DBUILD_SHARED_LIBS=OFF -DMANPAGE=OFF -DNOZLIB=ON
    cmake --build build/cmsgen --target cmsgen CopyPublicHeaders -j

Configure (adapt the paths to CaDiCaL):

    cmake -S . -B build \
      -DCADICAL_INCLUDE_DIR=/path/to/cadical/src \
      -DCADICAL_LIBRARY_DIR=/path/to/cadical/build \
      -Dcmsgen_DIR=$PWD/build/cmsgen

FLINT is searched in the default locations. If it is installed
elsewhere, add -DFLINT_INCLUDE_DIR=/path/to/include (the directory
containing flint/) and -DFLINT_LIBRARY=/path/to/libflint.so.

Build:

    cmake --build build -j

The executable will be located at:

    build/talisman

Test:

    ./build/talisman examples/abc6.aig -mult-spec
    cd build && ctest && cd ..   # library API test



If something goes wrong, remove the build directory
and repeat the steps above (including cmsgen):

    rm -rf build

------------------------------------------------------------
Usage
------------------------------------------------------------

    talisman <input-file> <spec-mode> [options]

Positional arguments:

    <input-file>      File containing the AIG graph
    <spec_file>       Custom specification file
                      (used when not using built-in specs)

Specification modes (choose one):

    -miter-spec       Use the miter specification
    -mult-spec        Use the multiplier specification
    -assert-spec      Use the assertion specification

If no built-in specification mode is used,
a custom specification file must be provided.

------------------------------------------------------------
General Options
------------------------------------------------------------

    -h, --help                Show help message and exit
    -v <int>                  Verbosity level (0..4)

    --custom-linalg <0|1>     Use custom linear algebra routines
                              Default: 1

    --primes <0|1>            Prime selection strategy
                              0 = structured
                              1 = random

    --reduction <0|1|2>       Reduction strategy
                              0 = hybrid (default)
                              1 = linear
                              2 = nonlinear

    --sc-depth <int>          Subcircuit depth (default: 2)

    --trust-abc               Trust the linear polynomial cache
                              (skip its SAT checks)

    --uniform-sampling <0|1>  Use uniform sampling (default: 1)

------------------------------------------------------------
Examples
------------------------------------------------------------

Run multiplier specification:

    ./build/talisman input.aig -mult-spec

Use a custom specification file:

    ./build/talisman input.aig spec.txt

Increase verbosity:

    ./build/talisman input.aig -mult-spec -v 3

Trust the linear polynomial cache:

    ./build/talisman input.aig -mult-spec --trust-abc

------------------------------------------------------------
Notes
------------------------------------------------------------

- Adders are detected by TalisMan itself, ABC is no longer needed.
- If unexpected behaviour occurs, delete the build directory and rebuild.

------------------------------------------------------------
Library
------------------------------------------------------------

Besides the executable, the build produces the static library
build/libtalisman.a (CMake target libtalisman), which is used, e.g.,
by the TalisMan application engine of the distributed SAT solving
platform Mallob. Its interface is include/talisman/talisman.hpp:
talisman::run(options, hooks) verifies one AIG, the hooks select the
SAT backend (default: CaDiCaL) and a terminate callback. All state is
reset after each run, errors are returned instead of exiting the process.
The library API test is test/api_test.cpp (repeated runs, custom SAT
backend, interruption).


