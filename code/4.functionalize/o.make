print-%  : ; @echo $* = $($*)

# define the C++ compiler to use (gcc, g++)
CC = g++

# define any compile-time flags
CFLAGS = -std=c++17 -Wall -Ofast

# location of rdkit installation
RDBASE=/Users/andreas/software/rdkit-master

# location of boost installation
BOOSTBASE=/opt/homebrew/Cellar/boost/1.83.0

# include and library flags for compilation
INCLUDES = -I$(RDBASE)/Code/ -I$(BOOSTBASE)/include
LFLAGS = -L$(RDBASE)/lib/ -L$(BOOSTBASE)/lib

# RDKit libraries to include for this compilation
LIBS = -lstdc++ -lRDKitGraphMol -lRDKitOptimizer -lRDKitAlignment -lRDKitForceFieldHelpers -lRDKitMolAlign -lRDKitSmilesParse -lRDKitForceField -lRDKitRDGeometryLib -lRDKitRDGeneral -lRDKitDistGeomHelpers -lRDKitDepictor -lRDKitFileParsers -lRDKitMolTransforms -lRDKitSubstructMatch -lRDKitFilterCatalog

# name of the main executable
MAIN=unsaturate


all:$(MAIN)

# compile the source code
$(MAIN):$(MAIN).cpp
	$(CC) $(CFLAGS) $(INCLUDES) $(MAIN).cpp -o $(MAIN) $(LFLAGS) $(LIBS)

# clean the directory
clean:
	$(RM) $(TARGET)

