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
import profile

## pymmlib
from mmLib import FileIO


def main(path):
    print "PERFORMANCE PROFILE: mmLib.LoadStructure(fil=%s)" % (path)
    profile.run("FileIO.LoadStructure(fil=path)")

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
