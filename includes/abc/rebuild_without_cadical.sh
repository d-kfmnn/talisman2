#!/bin/bash
# Script to rebuild ABC without built-in CaDiCaL (to avoid symbol conflicts)
# This allows linking with external CaDiCaL library

set -e

cd "$(dirname "$0")"

echo "Cleaning ABC build..."
make clean 2>/dev/null || true

# Remove the cadical object files from the archive
# We'll rebuild libabc.a excluding cadical

echo "Rebuilding ABC without built-in CaDiCaL..."

# Backup and disable cadical module
if [ ! -f src/sat/cadical/module.make.bak ]; then
    cp src/sat/cadical/module.make src/sat/cadical/module.make.bak
fi
echo "# CaDiCaL disabled - use external library instead" > src/sat/cadical/module.make

# Backup and disable kissat module
if [ ! -f src/sat/kissat/module.make.bak ]; then
    cp src/sat/kissat/module.make src/sat/kissat/module.make.bak
fi
echo "# Kissat disabled - kitten conflicts with external CaDiCaL" > src/sat/kissat/module.make

# Backup and disable eslim module (depends on kissat)
if [ ! -f src/opt/eslim/module.make.bak ]; then
    cp src/opt/eslim/module.make src/opt/eslim/module.make.bak
fi
echo "# eSLIM disabled - depends on kissat" > src/opt/eslim/module.make

# Build ABC
echo "Building ABC (this may take a while)..."
make libabc.a ABC_USE_NO_READLINE=1 -j$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu)


echo ""
echo "Done! libabc.a has been rebuilt without built-in CaDiCaL, Kissat, and eSLIM."
echo "You can now rebuild talisman."
