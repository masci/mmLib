#!/usr/bin/env python

## This example loads a mmCIF or PDB file, then changes the temperature factor
## of all Lysine atoms to 50.0, all others are set to 1.0.  The result is
## written to standard output as a PDB file.

import sys
from mmLib.FileLoader import LoadStructure, SaveStructure


def main(path):
    struct = LoadStructure(path)

    for res in struct.iter_amino_acids():

        temp_factor = 1.0
        if res.res_name == "LYS": temp_factor = 50.0

        for atm in res.iter_atoms():
            atm.temp_factor = temp_factor

    SaveStructure(sys.stdout, structure = struct, format = "PDB")


if __name__ == '__main__':
    try:
        path = sys.argv[1]
    except IndexError:
        print "usage %s <file path>" % (sys.argv[0])
        sys.exit(1)

    main(path)
