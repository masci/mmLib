#!/usr/bin/env python
## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
## PDB->mmCIF or mmCIF->PDB file translator
import sys
import os
import getopt
import string
from mmLib.FileLoader import LoadStructure, SaveStructure

if __name__ == '__main__':
    dir = "/home/jpaint/mmCIF"

    go = 0
    for path in os.listdir(dir):
	path = os.path.join(dir, path)

        print "[LOADING] ",path
        struct = LoadStructure(
            fil = path,
            build_properties = ("sequence","bonds"))

        print struct
        SaveStructure(sys.stdout, struct, "PDB")


