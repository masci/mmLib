#!/usr/bin/env python

## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import sys
import getopt
import string

from mmLib.FileLoader import LoadStructure, SaveStructure

def usage():
    print """
    NAME
       sec_struct.py - print out secondary structure

    SYNOPSIS
       sec_struct.py <file>

    """

if __name__ == '__main__':
    if "-h" in sys.argv or "--help" in sys.argv:
        usage()
        sys.exit(1)

    try:
        path = sys.argv[1]
    except KeyError:
        usage()
        sys.exit(1)

    if path == "-":
        path = sys.stdin

    struct = LoadStructure(fil = path)
    print struct

    print "Alpha Helicies:"
    for ahelix in struct.iter_alpha_helicies():
        print "    ",ahelix

    print "Beta Sheets:"
    for bsheet in struct.iter_beta_sheets():
        print "    ",bsheet

    print "Sites:"
    for site in struct.iter_sites():
        print "    ",site
        
