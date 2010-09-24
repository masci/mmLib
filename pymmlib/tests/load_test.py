#!/usr/bin/env python
## Copyright 2002-2010 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""This example uses the Python code profiler to load a file with mmLib.
We use this to help performance tune mmLib.
"""

## Python
import sys

## pymmlib
import test_util
from mmLib import FileIO


def main(path):
    print "mmLib.LoadStructure(fil=%s)" % (path)
    struct = FileIO.LoadStructure(fil=path)

if __name__ == "__main__":
    import os

    try:
        path = sys.argv[1]
    except IndexError:
        print "usage: load_test.py <PDB/mmCIF file or directory of files>"
        sys.exit(1)

    for pathx in test_util.walk_pdb_cif(path):
        main(pathx)
