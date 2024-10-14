#!/bin/bash

# Written by Andreas Luttens

# Verbose
echo "UniverseGenerator location is: ${UniverseGenerator}"

# RDKit root dir
if [ -z "${RDBASE}" ]; then
    echo "RDBASE is not set!\nLooking in UniverseGenerator folders"
    export RDBASE=${UniverseGenerator}/rdkit
fi

# Boost root dir
if [ -z "${BOOSTBASE}" ]; then 
    echo "BOOSTBASE is not set!\nLooking in conda env"
    export BOOSTBASE=${CONDA_PREFIX}
fi

# Link RDKit and conda
if [[ "$(uname)" == "Darwin" ]]; then
    echo "Working on MacOS"
    export DYLD_LIBRARY_PATH=$RDBASE/lib:$CONDA_PREFIX/lib
elif [[ "$(uname)" == "Linux" ]]; then
    echo "Working on Linux"
    export LD_LIBRARY_PATH=$RDBASE/lib:$CONDA_PREFIX/lib
else
    echo "Unknown OS."
fi

# Create shortcut variable for demo folder
export demo=${UniverseGenerator}/demo