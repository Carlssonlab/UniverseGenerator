#!/usr/bin/env python

# find the matches

import argparse
from openeye import oedepict, oechem

def parseArgs():
    """
    Parse the arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputfile', required=True)
    parser.add_argument('-o', '--outputfile', required=True)
    parser.parse_args()
    parser.set_defaults(verbose=False)
    args = parser.parse_args()
    return args

def depict(inputfile, outputfile):

    inputfileStream = oechem.oemolistream(inputfile)
    if not inputfileStream.open(inputfile):
        OEThrow.Fatal("Unable to open %s inputfile for reading" % inputfile)

    multi = oedepict.OEMultiPageImageFile(oedepict.OEPageOrientation_Landscape, oedepict.OEPageSize_US_Letter)
    image = multi.NewPage()

    opts = oedepict.OE2DMolDisplayOptions()
    opts.SetAromaticStyle(oedepict.OEAromaticStyle_Circle)
    #opts.SetBondLineGapScale(5)
    opts.SetTitleFontScale(2.0)
    # SHOW ME ALL HYDROGENS
    #opts.SetHydrogenStyle(oedepict.OEHydrogenStyle_ImplicitAll)

    # add lines to change the title here.
    rows, cols = 4, 3
    grid = oedepict.OEImageGrid(image, rows, cols)
    grid.SetCellGap(0)
    grid.SetMargins(0)
    citer = grid.GetCells()

    for mol in inputfileStream.GetOEGraphMols():
        mol.SetTitle("%s\n%s"%(oechem.OEMolToSmiles(mol), mol.GetTitle()))

        if not citer.IsValid():
            # go to next page
            image = multi.NewPage()
            grid = oedepict.OEImageGrid(image, rows, cols)
            grid.SetCellGap(0)
            grid.SetMargins(0)
            citer = grid.GetCells()

        cell = citer.Target()
        oedepict.OEPrepareDepiction(mol)
        opts.SetDimensions(cell.GetWidth(), cell.GetHeight(), oedepict.OEScale_AutoScale)
        #opts.SetMargin(oedepict.OEMargin_Bottom, 25.0)
        disp = oedepict.OE2DMolDisplay(mol, opts)
        clearbackground = False
        oedepict.OERenderMolecule(cell, disp, clearbackground)
        oedepict.OEDrawBorder(cell, oedepict.OEPen(oedepict.OEBlackPen))
        citer.Next()

    oedepict.OEWriteMultiPageImage(outputfile, multi)

def main():

    args = parseArgs()

    depict(args.inputfile, args.outputfile)

if __name__ == '__main__':
    main()
