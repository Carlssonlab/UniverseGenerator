#!/bin/bash

# Written by Andreas Luttens

# Download nauty and Traces
wget http://users.cecs.anu.edu.au/~bdm/nauty/nauty2_8_8.tar.gz

# Install the package
tar -xvzf nauty2_8_8.tar.gz
cd nauty2_8_8
./configure
make

# Put key executables in UniverseGenerator's bin folder
cp ${UniverseGenerator}/nauty2_8_8/geng ${UniverseGenerator}/bin/geng
cp ${UniverseGenerator}/nauty2_8_8/planarg ${UniverseGenerator}/bin/planarg
cp ${UniverseGenerator}/nauty2_8_8/showg ${UniverseGenerator}/bin/showg
