#!/usr/bin/env python
## Copyright 2002-2010 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import sys
import getopt
import string

from mmLib import FileIO

def usage():
    print """
    NAME
       sec_struct.py - print out secondary structure

    SYNOPSIS
       sec_struct.py <file>
       sec_struct.py -v <file> (verbose output

    """
    

def main(path, verbose):
    
    struct = FileIO.LoadStructure(file = path)
    print struct

    print "Alpha Helicies:"
    for ahelix in struct.iter_alpha_helicies():
        print "    ",ahelix

        if verbose == True:
            print "        ",ahelix.segment

    print "Beta Sheets:"
    for bsheet in struct.iter_beta_sheets():

        print "    ",bsheet

        if verbose == True:
            for strand in bsheet.iter_strands():
                print "        ",strand
                print "            ",strand.segment

    print "Sites:"
    for site in struct.iter_sites():
        print "    ",site

        if verbose == True:
            for frag in site.iter_fragments():
                print "        ",frag


if __name__ == '__main__':
   ## parse options
    (opts, args) = getopt.getopt(sys.argv[1:], "h?v")

    verbose = False

    for (opt, item) in opts:
        if opt=="-h" or opt=="-?":
            usage()
            sys.exit(1)

        elif opt=="-v":
            verbose = True

    try:
        path = args[0]
    except KeyError:
        usage()
        sys.exit(1)

    if path == "-":
        path = sys.stdin

    main(path, verbose)
