# talisman

TalisMan Version 2.0
is a verification and analysis tool for AIG-based circuits.
It combines algebraic reasoning and SAT-based techniques using
CaDiCaL, FLINT, PBLib, and cmsgen.
These tools will automatically be fetched and locally installed.

In the corresponding IJCAR'26 paper 
C. Hofstadler, D.Kaufmann, and C. Chen
Avoiding Big Integers: Parallel Multimodular Algebraic Verification of Arithmetic Circuits

we have used the following versions of the dependencies
CaDiCaL 2.2.0
FLINT 3.4.0
PBLib 1.2.1
cmsgen 1.0



------------------------------------------------------------
Prerequisites
------------------------------------------------------------

You need:

- CMake >= 3.16
- GCC or Clang
- Autoconf, Automake, and Libtool (required for FLINT and CaDiCaL)


Linux (Ubuntu):

    sudo apt install \
    build-essential \
    cmake \
    autoconf \
    automake \
    libtool 


macOS:

    brew install cmake git autoconf automake libtool


------------------------------------------------------------
Building
------------------------------------------------------------

Configure:

    cmake -S . -B build

Build:

    cmake --build build -j

The executable will be located at:

    build/talisman



If something goes wrong, remove the build directory and rebuild:

    rm -rf build
    cmake -S . -B build
    cmake --build build -j

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

    -primes <0|1>             Prime selection strategy
                              0 = structured
                              1 = random

    --reduction <0|1|2>       Reduction strategy
                              0 = hybrid (default)
                              1 = linear
                              2 = nonlinear

    --sc-depth <int>          Subcircuit depth (default: 2)

    --no-trust-abc            Do not trust ABC's linear polynomial cache
                              (perform SAT checks instead)

    --uniform-sampling <0|1>  Disable uniform sampling (default: 1)

------------------------------------------------------------
Examples
------------------------------------------------------------

Run multiplier specification:

    ./build/talisman input.aig -mult-spec

Use a custom specification file:

    ./build/talisman input.aig spec.txt

Increase verbosity:

    ./build/talisman input.aig -mult-spec -v 3

Disable trusting ABC:

    ./build/talisman input.aig -mult-spec --no-trust-abc

------------------------------------------------------------
Notes
------------------------------------------------------------

- Dependencies are fetched and built automatically by CMake.
- The first build may take several minutes.
- If unexpected behaviour occurs, delete the build directory and rebuild.



