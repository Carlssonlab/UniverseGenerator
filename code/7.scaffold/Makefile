# Author: Andreas Luttens
# Date: June 22, 2023

print-%  : ; @echo $* = $($*)

# define the C++ compiler to use (gcc, g++)
CC = g++

# define any compile-time flags
CFLAGS = -std=c++17 -Wall -Ofast

# include and library flags for compilation
INCLUDES = -I$(RDBASE)/Code/ -I$(BOOSTBASE)/include
LFLAGS = -L$(RDBASE)/lib/ -L$(BOOSTBASE)/lib

# RDKit libraries to include for this compilation
LIBS = -lstdc++ -lRDKitGraphMol -lRDKitAlignment -lRDKitMolAlign -lRDKitSmilesParse -lRDKitRDGeometryLib -lRDKitRDGeneral -lRDKitFileParsers -lRDKitMolTransforms -lRDKitSubstructMatch -lRDKitChemReactions

# name of the main executable
MAIN=activate-scaffold

all:$(MAIN)

# compile the source code
$(MAIN):$(MAIN).cpp
	$(CC) $(CFLAGS) $(INCLUDES) $(MAIN).cpp -o $(MAIN) $(LFLAGS) $(LIBS)

# clean the directory
clean:
	$(RM) $(TARGET)
