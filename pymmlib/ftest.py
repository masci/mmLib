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

    mmCIFDirectory = "/home/jpaint/mmCIF"

    for path in os.listdir(mmCIFDirectory):
        path = os.path.join(mmCIFDirectory, path)


        print "[SCANNING] ",path
        structure = LoadStructure(path)

        print "[SAVING] ",structure.getTitle()
        SaveStructure("current.pdb", structure)

