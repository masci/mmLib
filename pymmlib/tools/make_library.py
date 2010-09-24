#!/usr/bin/env python
## Copyright 2002-2010 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

## Python
import os
import sys
import copy

## pymmlib
from mmLib.mmCIF import mmCIFFile

def usage():
    print """
    usage: python make_library.py <component.cif file>

    Utility to create the mmLib monomer library located
    in mmLib/Data/Monomers. It uses the RCSB public component dictionary
    available on the web at: http://pdb.rutgers.edu/public-component-erf.cif

    When this utility is run, it creates the Monomers directory in the
    current working directory. To use it, it must be moved to its
    loading location at mmLib/Data/Monomers.
    """

def main(path):
    print "parsing %s" % (path)

    cif_file = mmCIFFile()
    cif_file.load_file(path)

    os.mkdir("Monomers")

    for cif_data in cif_file:

        mkdir_path = os.path.join("Monomers", cif_data.name[0])
        if not os.path.isdir(mkdir_path):
            os.mkdir(mkdir_path)

        save_path = os.path.join(mkdir_path, "%s.cif" % (cif_data.name))

        print "saving %s" % (save_path)

        cf = mmCIFFile()
        cf.append(copy.deepcopy(cif_data))
        cf.save_file(save_path)


if __name__ == "__main__":
    try:
        path = sys.argv[1]
    except:
        usage()
        sys.exit(1)

    main(path)
