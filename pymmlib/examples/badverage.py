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

    num_atoms = 0

    mean_B  = 0.0
    sigma_B = 0.0

    min_B = 1000.0
    max_B = 0.0


    for res in struct.iter_amino_acids():
        for atm in res.iter_atoms():
            if atm.element == "H":
                continue

            num_atoms += 1
            mean_B    += atm.temp_factor

            min_B = min(min_B, atm.temp_factor)
            max_B = max(max_B, atm.temp_factor)

    
            
    mean_B = mean_B / num_atoms

    print "mean B: ",mean_B
    print "max  B: ",max_B
    print "min  B: ",min_B

try:
    path = sys.argv[1]
except IndexError:
    print "usage: badv.py <PDB/mmCIF file>"
else:
    main(path)
