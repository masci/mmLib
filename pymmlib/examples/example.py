#!/usr/bin/env python

## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

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
