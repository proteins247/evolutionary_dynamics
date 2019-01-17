#!/bin/bash

# install_folding_evolution.sh

# stand-in for a makefile at the moment

# Example
# DESCRIPTION="this is my description" ./install_folding_evolution.sh 0.0.1


if [ $# -lt 1 ]; then
    echo "Need version number argument: ./install_folding_evolution.sh VERSION_NUMBER"
    echo "exiting"
    exit 1
fi

: ${VERSION:=$1}

: ${FOLDEVODIR:=/n/home00/vzhao/code/evolutionary_dynamics/src/app_folding_evolution/}
: ${GENEDIR:=/n/home00/vzhao/code/evolutionary_dynamics/src/app_amy_LP_structure/}
: ${PKGDIR:=/n/home00/vzhao/pkg/FoldEvo}
: ${MODULEFILEDIR:=/n/home00/vzhao/modulefiles/FoldEvo}

if [ -f /usr/local/bin/centos7-modules.sh ]; then
    source centos7-modules.sh
    module add gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02 hdf5/1.10.1-fasrc01
else
    echo "Run this on CentOS 7."
    echo "exiting"
    exit 1
fi


# ------------------------------------------------
# Make the install path
# ------------------------------------------------
BINDIR="${PKGDIR}/${VERSION}/bin"
SHAREDIR="${PKGDIR}/${VERSION}/share"

if [ -d "${BINDIR}" ]; then
    echo "Path already exists (${PKGDIR}/${VERSION})"
    echo "Install anyway (will overwrite)? ('yes' to proceed)"
    read response
    if [ "$response" != "yes" ]; then
	echo "Exiting"
	exit
    fi
fi
echo "Creating paths: ${BINDIR}, ${SHAREDIR}"
mkdir -p "${BINDIR}"
mkdir -p "${SHAREDIR}"
echo


# ------------------------------------------------
# Install
# ------------------------------------------------
cd $FOLDEVODIR
echo "Compiling folding_evolution..."
mpic++ -Wall -O3 -std=c++11 folding_evolution.cpp -o folding_evolution -lhdf5 \
       ../rng.o ../gencode.o
if [ $? -ne 0 ]; then
    echo "Compile FoldEvo failed"
    exit -1
fi
echo "Compiling generateStableGene..."
cd $GENEDIR
gcc -Wall -g -o generateStableGene generateStableGene.c general.c \
    structurelib.c ../gencode.c ../latticelib.c ../rng.c -lz -lm
if [ $? -ne 0 ]; then
    echo "Compile generateStableGene failed"
    exit -1
fi
echo

echo "Installing..."
cd $FOLDEVODIR
cp -f folding_evolution "${BINDIR}"
cd $GENEDIR
cp -f generateStableGene "${BINDIR}"
/bin/cp -rf ../../share/* "${SHAREDIR}"

echo


# ------------------------------------------------
# Make a lua file
# ------------------------------------------------
LUAFILEPATH="${MODULEFILEDIR}/${VERSION}.lua"
echo ".lua file at $LUAFILEPATH"
cat >"${LUAFILEPATH}" <<EOF
local helpstr = [[
FoldEvo ${VERSION} -- Protein folding evolution code.

$(echo ${DESCRIPTION} | fold -s -w 70)

Installed on $(date).

Changelog:

]]

help(helpstr)

whatis("Name: FoldEvo")
whatis("Version: ${VERSION}")
whatis("Description: Evolutionary simulations of protein folding.")

setenv("FOLDEVO_SHARE", "${SHAREDIR}")
prepend_path("PATH", "${BINDIR}")
EOF

echo
echo "Printing contents of .lua file:"
cat "$LUAFILEPATH"
