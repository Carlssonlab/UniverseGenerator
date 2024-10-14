#!/bin/bash

# Written by Andreas Luttens

home=$PWD

# Iterate over all steps in the process, compile the codes
for directory in $(cat ${UniverseGenerator}/code/steps.list)
do

    # Enter the step's directory
    cd ${UniverseGenerator}/code/${directory}

    echo "Compiling" ${directory}

    # Compile the code
    make

    # Find all source codes
    for code in $(ls *.cpp)
    do
        # Fine the name of the executable
        executable=$(echo ${code} | sed 's|.cpp||g')

	# Put key executables in UniverseGenerator's bin folder
	cp ${executable} ${UniverseGenerator}/bin/${executable}

    done
	
    # Remove remaining object files
    if ls *.o 1> /dev/null 2>&1; then
        rm *.o
    fi

    # Go back to previous directory
    cd ${home}

    # Newline
    echo ""

done

