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
LIBS = -lstdc++ -lRDKitGraphMol -lRDKitAlignment -lRDKitMolAlign -lRDKitSmilesParse -lRDKitRDGeometryLib -lRDKitRDGeneral -lRDKitFileParsers -lRDKitMolTransforms -lRDKitSubstructMatch -lRDKitChemReactions -lRDKitChemTransforms

# Source files for program 1
SRCS1 = build-instructions.cpp

# Object files for program 1
OBJS1 = $(SRCS1:.cpp=.o)

# Executable name for program 1
EXEC1 = build-instructions

# Source files for program 2
SRCS2 = build-superstructures.cpp

# Object files for program 2
OBJS2 = $(SRCS2:.cpp=.o)

# Executable name for program 2
EXEC2 = build-superstructures

# Default target
all: $(EXEC1) $(EXEC2)

# Compile each source file for program 1 into an object file
%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Compile each source file for program 2 into an object file
%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Link object files for program 1 to create executable
$(EXEC1): $(OBJS1)
	$(CC) $(OBJS1) $(LFLAGS) $(LIBS) -o $@

# Link object files for program 2 to create executable
$(EXEC2): $(OBJS2)
	$(CC) $(OBJS2) $(LFLAGS) $(LIBS) -o $@

# Clean target to remove object files and executables
clean:
	rm -f $(OBJS1) $(EXEC1).o $(OBJS2) $(EXEC2).o
