#!/bin/bash

# Written by Andreas Luttens

home=$PWD

# Iterate over all substituent files
for sub_file in $(cat ${UniverseGenerator}/auxiliaries/substituent.list)
do

    # Create a symbolic link
    ln -s ${sub_file}

done