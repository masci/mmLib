#!/usr/bin/env python

## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

## PDB->mmCIF or mmCIF->PDB file translator

import sys
import getopt
import string

from mmLib.FileLoader import LoadStructure, SaveStructure

def usage():
    print "NAME"
    print "  pdb2cif.py - convert a PDB file to a mmCIF file"
    print
    print "SYNOPSIS"
    print "  pdb2cif.py [-h|--help]" 
    print "  pdb2cif.py [input] [output]"
    print
    print "DESCRIPTION"
    print "  input     path of input file, or '-' for stdin"
    print "  output    path of output file, or '-' for stdout"
    print "  -h, --help  display this help and exit" 
    print
    print "EXAMPLE"
    print "  Translate mmCIF file lysozyme.cif to lysozyme.pdb:"
    print "    # python pdb2cif.py lysozyme.pdb lysozyme.cif"
    print
    print "  Translate a mmCIF file from standard input to a PDB file"
    print "  written to standard output (2 equivelent forms)"
    print "    # cat lysozyme.cif | ./cif2pdb.py > lysozyme.pdb"
    print "    # cat lysozyme.cif | ./cif2pdb.py - - > lysozyme.pdb"
    print
    print "AUTHORS"
    print "  Jay Painter <jpaint@u.washington.edu>"


if __name__ == '__main__':
    if "-h" in sys.argv or "--help" in sys.argv:
        usage()
        sys.exit(0)

    cif = "-"
    pdb = "-"

    if   len(sys.argv) == 1:
        pass
    elif len(sys.argv) == 2:
        cif = sys.argv[1]
    elif len(sys.argv) == 3:
        cif = sys.argv[1]
        pdb = sys.argv[2]
    else:
        usage()
        sys.exit(1)

    if cif == "-": cif = sys.stdin
    if pdb == "-": pdb = sys.stdout

    struct = LoadStructure(fil = cif, format = "PDB")
    SaveStructure(fil = pdb, structure = struct, format = "CIF")
