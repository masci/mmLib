#!/usr/bin/env python

import sys
from mmLib.FileLoader import LoadStructure, SaveStructure

## load structure
structure = LoadStructure(sys.argv[1])

## list of anisotropic values
anisou_list = []

## iterate over all atoms in the structure
for atm in structure.atomIterator():
    ## get the atom's parent molecule, if it is part of a
    ## AminoAcidResidue, then calculate the atom's anisotropy
    ## and add it to the list of anisotropic values
    aa_res = atm.getAminoAcidResidue()
    if aa_res:
        anisou = atm.calcAnisotropy()	
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
    print "(%f,%f) " % (amin, amax), count(amin, amax)
    amax += bin_size
