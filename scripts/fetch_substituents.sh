#!/bin/bash

# Written by Andreas Luttens

home=$PWD

# Iterate over all substituent files
for sub_file in $(cat ${UniverseGenerator}/auxiliaries/substituents/substituents.list)
do

    # Create a symbolic link
    ln -s ${UniverseGenerator}/auxiliaries/substituents/${sub_file}

done