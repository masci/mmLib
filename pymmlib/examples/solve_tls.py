#!/usr/bin/env python

## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import sys
from mmLib.Structure import *
from mmLib.FileLoader import *
from mmLib.Extensions.TLS import *

def usage():
    print "solve_tls.py <file path>"
    print
    print "description:"
    print "    compute TLS tensor from anisotropic U values"
    print "    in a structure input file using all amino acid atoms"
    print

def main(path):
    print """
    Calculating 1 TLS group per chain by least squares.
    """
    
    struct = LoadStructure(fil = path)
    tls_group = TLSGroup()

    for chain in struct.iter_chains():
        print chain
        
        tls_group = TLSGroup()

        for aa in chain.iter_amino_acids():
            for atm in aa.iter_atoms():
                tls_group.append(atm)

        print "    TLS GROUP ATOMS: ",len(tls_group)

        if len(tls_group)<20:
            print "    NOT ENOUGH ATOMS IN CHAIN"
            continue
            

        tls_group.origin = tls_group.calc_centroid()

        print "    CENTROID: ", tls_group.origin
        print

        tls_group.calc_TLS_least_squares_fit()
        tls_group.write()


if __name__ == "__main__":
    try:
        path = sys.argv[1]
    except IndexError:
        usage()
        sys.exit(1)

    main(path)
