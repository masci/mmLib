#!/usr/bin/env python

## This example uses the Python code profiler to load a file with mmLib
## I use this to help performance tune mmLib

import sys
import profile
from mmLib.FileLoader import LoadStructure


def main(path):
    print "PERFORMANCE PROFILE: mmLib.LoadStructure(fil=%s)" % (path)
    profile.run("LoadStructure(fil=path)")

if __name__ == "__main__":
    import os

    try:
        path = sys.argv[1]
    except IndexError:
        print "usage: profile_load.py <PDB/mmCIF file or directory of files>"
        sys.exit(1)

    if os.path.isfile(path):
        main(path)
    elif os.path.isdir(path):
        for name in os.listdir(path):
            name = os.path.join(path, name)
            if not os.path.isfile(name):
                continue
            try:
                main(name)
            except:
                print "ERROR: ",name
                raise
