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
    print "  filtrans.py - read in a PDB or mmCIF file, and write it"
    print "                out as a mmCIF or PDB file"
    print
    print "SYNOPSIS"
    print "  filtrans.py [OPTION] <input> <output>"
    print
    print "DESCRIPTION"
    print "  <input>     path of input file, if the path is suffixed"
    print "              by a '-' it is read from standard input, but"
    print "              the path extention is still needed to determine"
    print "              the file type"
    print "  <output>    same as above, but '-' is for standard output"
    print "  -h, --help  display this help and exit" 
    print
    print "EXAMPLE"
    print "  Translate mmCIF file 102l.cif to 102l.pdb:"
    print "    # python filtrans.py 102l.cif 102l.pdb"
    print
    print "  Translate a mmCIF file from standard input to a PDB file"
    print "  written to standard output:"
    print "    # python filtrans.py .pdb- .cif-"
    print
    print "AUTHORS"
    print "  Jay Painter <jpaint@u.washington.edu>"


if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", ['help'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)

    if len(args) != 2:
        usage()
        sys.exit(0)

    iformat = ""
    oformat = ""
    ifil = args[0]
    ofil = args[1]
    
    if ifil[-1] == "-":
        iformat = ifil[:-1]
        ifil = sys.stdin

    if ofil[-1] == "-":
        oformat = ofil[:-1]
        ofil = sys.stdout

    structure = LoadStructure(ifil, iformat)
    SaveStructure(ofil, structure, oformat)
