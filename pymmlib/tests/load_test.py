#!/usr/bin/env python

## This example uses the Python code profiler to load a file with mmLib
## I use this to help performance tune mmLib

import sys
import test_util
from mmLib.FileLoader import LoadStructure,SaveStructure


def main(path):
    print "mmLib.LoadStructure(fil=%s)" % (path)
    struct = LoadStructure(fil=path)
    SaveStructure(fil=sys.stdout, struct=struct)

if __name__ == "__main__":
    import os

    try:
        path = sys.argv[1]
    except IndexError:
        print "usage: load_test.py <PDB/mmCIF file or directory of files>"
        sys.exit(1)

    for pathx in test_util.walk_pdb_cif(path):
        main(pathx)
