#!/usr/bin/env python

## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

## This example program calculates the adverage B value of 
## the protein atoms of a structure.

import sys
from mmLib.FileLoader import LoadStructure, SaveStructure
from mmLib.Structure import *

def main(path):

    ## load structure
    struct = LoadStructure(
        fil = path,
        build_properties = ("no_bonds",))

    atom_list = AtomList()

    for res in struct.iter_amino_acids():
        for atm in res.iter_atoms():
            if atm.element == "H":
                continue
            atom_list.append(atm)

    print "ADV B: ",atom_list.calc_adv_temp_factor()

try:
    path = sys.argv[1]
except IndexError:
    print "usage: badv.py <PDB/mmCIF file>"
else:
    main(path)
