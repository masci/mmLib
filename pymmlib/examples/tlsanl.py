#!/usr/bin/env python

## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import sys
import re
import string

from mmLib.PDB import PDBFile
from mmLib.Structure import *
from mmLib.FileLoader import LoadStructure, SaveStructure
from mmLib.Extensions.TLS import *


def usage():
    print "tlsanl.py <xxx.pdb> [REFMAC TLS File]"
    print
    print "description:"
    print "    Compute anisotropic ADP records from the given TLS"
    print "    description.  The TLS description is taken from the"
    print "    REMARK fields in the PDB file, or from the TLSOUT file"
    print "    written by REFMAC."
    print


def astr(a):
    a0 = "%.3f" % (a[0])
    a1 = "%.3f" % (a[1])
    a2 = "%.3f" % (a[2])

    return a0.rjust(7) + a1.rjust(7) + a2.rjust(7)

def print_TLS(text, Tt, T, Lt, L, St, S):
    print text

    s0 = "%s TENSOR" % (Tt)
    s1 = "%s TENSOR" % (Lt)
    s2 = "%s TENSOR" % (St)

    print "".ljust(5) + s0.ljust(24) + s1.ljust(24) + s2
    print "     (A^2)                   (DEG^2)                 (A DEG)"

    L = L * rad2deg2
    S = S * rad2deg

    for i in range(3):
        print "   %s   %s   %s" % (astr(T[i]), astr(L[i]), astr(S[i]))


def main(pdb_path, tls_out_path = None):

    struct     = LoadStructure(fil = pdb_path)
    tls_groups = TLSGroupFile()

    ## extract TLS tensors from PDB REMARK records
    if tls_out_path == None:
        pdb_file = PDBFile()
        pdb_file.load_file(pdb_path)
        pdb_file.record_processor(tls_groups)

    ## or get TLS tensors from REMAC tlsout file
    else:
        try:
            fil = open(tls_out_path)
        except IOError, e:
            print "[Error] %s: %s" % (str(e), tls_out_path)
        tls_groups.load_refmac_tlsout_file(fil)

    print "<TLSGroups>"
    print tls_groups
    print "</TLSGroups>"
    print
    print

    for tls_info in tls_groups.tls_info_list:
        tls = tls_info.make_tls_group(struct)
        calcs = tls.calc_COR()

        print_TLS(
            "INPUT TENSOR MATRICES WRT ORTHOGONAL AXES USING ORIGIN "\
            "OF CALCULATIONS",
            "T", tls.T,
            "L", tls.L,
            "S", tls.S)

        print
        print "TRACE OF TRANSLATION TENSOR               %.3f" % (
            trace(tls.T))
        print "MEAN TRANSLATION (TRACE/3)                %.3f" % (
            trace(tls.T)/3.0)
        print "MEAN LIBRATION   (TRACE/3)                %.3f" % (
            trace(tls.L * rad2deg2)/3.0)
        print

        print_TLS(
            "TENSOR MATRICES WRT LIBRATION AXES USING ORIGIN OF CALCULATIONS",
            "T^", calcs["T^"],
            "L^", calcs["L^"],
            "S^", calcs["S^"])

        print
        print "ORIGIN SHIFT RHO(O)^ TO CENTRE WRT LIBRATION AXES (A): "+\
              astr(calcs["RHO^"])
        print "ORIGIN SHIFT TO CENTRE WRT ORTHOGONAL AXES        (A): "+\
              astr(calcs["RHO"])
        print "TLS CENTRE OF REACTION WRT ORTHOGONAL AXES        (A): "+\
              astr(calcs["COR"])
        print

        print_TLS(
            "TENSOR MATRICES WRT LIBRATION AXES USING CENTRE OF REACTION",
            "T'^", calcs["T'^"],
            "L'^", calcs["L'^"],
            "S'^", calcs["S'^"])


if __name__ == "__main__":
    try:
        pdb_path = sys.argv[1]
    except IndexError:
        usage()
        sys.exit(1)

    try:
        tls_out_path = sys.argv[2]
    except IndexError:
        tls_out_path = None

    main(pdb_path, tls_out_path)
