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
    Calculating TLS parameters for a single rigid body group composed of
    all the amino acids
    """
    
    struct = LoadStructure(fil = path)
    tls_group = TLSGroup()

    for res in struct.iter_amino_acids():
        for atm in res.iter_atoms():
            tls_group.append(atm)

    print "Atoms: ",len(tls_group)

    tls_group.origin = tls_group.calc_centroid()

    print "Centroid of Atoms:"
    print tls_group.origin
    print

    tls_group.calc_tls_tensors()
    tls_group.write()


if __name__ == "__main__":
    try:
        path = sys.argv[1]
    except IndexError:
        usage()
        sys.exit(1)

    main(path)
