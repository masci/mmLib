#!/usr/bin/env python
## Copyright 2002-2010 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os, sys
from mmLib import PDB, Library

AMINO_ACID13_DICT = {
    'A': 'ALA', 'C': 'CYS', 'E': 'GLU', 'D': 'ASP', 'G': 'GLY',
    'F': 'PHE', 'I': 'ILE', 'H': 'HIS', 'K': 'LYS', 'M': 'MET',
    'L': 'LEU', 'N': 'ASN', 'Q': 'GLN', 'P': 'PRO', 'S': 'SER',
    'R': 'ARG', 'T': 'THR', 'W': 'TRP', 'V': 'VAL', 'Y': 'TYR'}


def usage():
    print "seq2seqres.py: convert a 1-letter protein sequence into"
    print "               PDB file SEQRES records"
    print "usage: seq2seqres.py <chain_id> <1-leter-code sequence>"
    print
    sys.exit(-1)


def main(chain_id, sequence):
    pdb_file = PDB.PDBFile()

    ## parse sequence into 
    sequence3  = []
    i = 0
    while i < len(sequence):
        c = sequence[i]
        if c=="(":
            sequence3.append(sequence[i+1:i+4])
            i += 5
        else:
            sequence3.append(sequence[i])
            i += 1

    for i in range(len(sequence3)):
        res = sequence3[i]        
        if AMINO_ACID13_DICT.has_key(res):
            sequence3[i] = AMINO_ACID13_DICT[res]

    ## generate PDB records
    serial  = 0
    res_num = 0
    newrec  = True

    for res in sequence3:

        if newrec:
            newrec = False

            serial += 1
            res_num = 1
            
            seqres = PDB.SEQRES()
            pdb_file.append(seqres)
            
            seqres["serNum"] = serial
            seqres["chainID"] = chain_id
            seqres["numRes"] = len(sequence3)

        seqres["resName%d"%(res_num)] = res

        res_num += 1
        if res_num>13:
            newrec = True

    pdb_file.save_file(sys.stdout)


if __name__ == "__main__":
    import os

    try:
        chain_id = sys.argv[1]
        sequence = sys.argv[2]
    except IndexError:
        usage()

    if len(chain_id)>1:
        usage()

    main(chain_id, sequence)
