#!/usr/bin/env python
## Copyright 2002-2010 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

## This example program calculates the adverage B value of 
## the protein atoms of a structure.

import sys
from mmLib import FileIO

def main(path):
    ## load structure
    struct = FileIO.LoadStructure(file = path)

    num_atoms = 0

    mean_B  = 0.0
    sigma_B = 0.0

    min_B = 1000.0
    max_B = 0.0

    for atm in struct.iter_all_atoms():
        num_atoms += 1
        mean_B    += atm.temp_factor
        
        min_B = min(min_B, atm.temp_factor)
        max_B = max(max_B, atm.temp_factor)

    if num_atoms > 0:
        mean_B = mean_B / num_atoms
        print "mean B: ",mean_B
        print "max  B: ",max_B
        print "min  B: ",min_B

    else:
        print "No Atoms Found"

try:
    path = sys.argv[1]
except IndexError:
    print "usage: calc_Bmean.py <PDB/mmCIF file>"
else:
    main(path)
