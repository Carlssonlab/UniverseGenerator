#!/Users/andreasluttens/anaconda2/bin/python

"""
Translate a file into another.
[Andreas Luttens - 2018]
"""

# basic modules
import argparse

# openeye modules
from openeye.oechem import *

def ParseArgs():
    """
    Parse the arguments.
    """
    parser = argparse.ArgumentParser(description='A script for translating one file into another.', epilog= 'GitHub URL')
    parser.add_argument('-i1','--inputfile1', help = 'path to input file', required=True)
    parser.add_argument('-i2','--inputfile2', help = 'path to input file', required=True)
    parser.add_argument('-o','--outputfile', help = 'name of the outputfile', required=True)
    parser.set_defaults(verbose=False)
    args = parser.parse_args()
    return args

def main():
    """
    Run the main script.
    """

    # retrieve the parsed arguments
    args = ParseArgs()

    inputfileStream = oemolistream(args.inputfile1)
    if not inputfileStream.open(args.inputfile1):
        OEThrow.Fatal("Unable to open %s inputfile for reading" % args.inputfile1)

    inputfileStream2 = oemolistream(args.inputfile2)
    if not inputfileStream2.open(args.inputfile2):
        OEThrow.Fatal("Unable to open %s inputfile for reading" % args.inputfile2)

    # create the outputfile stream
    outputfileStream = oemolostream(args.outputfile)
    if not outputfileStream.open(args.outputfile):
        OEThrow.Fatal("Unable to open %s outputfile for writing" % args.outputfile)

    seen = set()

    # iterate over all molecules in the inputfile, write them to the outputfile
    for molecule in inputfileStream.GetOEGraphMols():

        smiles = OECreateIsoSmiString(molecule)

        seen.add(smiles)
    
    for molecule in inputfileStream2.GetOEGraphMols():

        smiles = OECreateIsoSmiString(molecule)

        if smiles not in seen:
            OEWriteMolecule(outputfileStream, molecule)

    
if __name__ == '__main__':
    main()



