#!/usr/bin/env python
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
    struct = LoadStructure(fil = path)
    tls_group = TLSGroup()

    for res in struct.iter_amino_acids():
        for atm in res.iter_atoms():
            tls_group.append(atm)

    print len(tls_group),"atoms"

    tls_group.origin = tls_group.calc_centroid()
    print "centroid=",tls_group.origin

    tls_group.calc_tls_tensors()
    print tls_group


if __name__ == "__main__":
    try:
        path = sys.argv[1]
    except IndexError:
        usage()
        sys.exit(1)

    main(path)
