#!/usr/bin/env python

## This example program calculates the distribution of ansotropy
## for a protein structure

import sys
from mmLib.FileLoader import LoadStructure, SaveStructure

try:
    path = sys.argv[1]
except IndexError:
    print "usage: ansotropy.py <PDB/mmCIF file>"

def main(path):
    ## load structure
    print "Loading ",path
    struct = LoadStructure(path, build_properties=("polymers",))

    ## list of anisotropic values
    anisou_list = []

    ## iterate over all atoms in the structure
    for res in struct.iter_amino_acids():
        for atm in res.iter_atoms():
            anisou = atm.calc_anisotropy()	
            anisou_list.append(anisou)

    ## this function counts the number of anisotropic values in
    ## anisou_list in the range amin to amax
    def count(amin, amax):
        c = 0
        for x in anisou_list:
            if x >= amin and x < amax: c += 1 
        return c

    bin_size = 0.05
    amax = bin_size
    while amax <= 1.001:
        amin = amax - bin_size
        print "Ansotropy Range (%f,%f): %d atoms " % (
            amin, amax, count(amin, amax))
        amax += bin_size


try:
    path = sys.argv[1]
except IndexError:
    print "usage: ansotropy.py <PDB/mmCIF file>"
else:
    main(path)

print "END"
