#!/usr/bin/env python

## This example loads a mmCIF or PDB file, then changes the temperature factor
## of all Lysine atoms to 50.0, all others are set to 1.0.  The result is
## written to standard output as a PDB file.

import sys
from mmLib.FileLoader import LoadStructure, SaveStructure


def main(path):
    structure = LoadStructure(path)

    for aa_res in structure.aminoAcidResidueIterator():
        temp_factor = 1.0
        if aa_res.getName() == "LYS":
            temp_factor = 50.0

        for atm in aa_res.atomIterator():
            atm.setTemperatureFactor(temp_factor)

    SaveStructure(sys.stdout, structure, ".pdb")


if __name__ == '__main__':
    try:
        path = sys.argv[1]
    except IndexError:
        print "usage %s <file path>" % (sys.argv[0])
        sys.exit(1)

    main(path)
