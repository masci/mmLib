#!/usr/bin/env python

## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import sys
import getopt

from mmLib.Structure      import *
from mmLib.FileLoader     import *
from mmLib.Extensions.TLS import *

def usage():
    print "lsq_tls_fit.py [-t <tlsin>] [-p <pdbout>] <structure file>"
    print
    print "description:"
    print "    Performs a least squares fit of TLS tensors to the"
    print "    tempature factors of the given structure file.  If"
    print "    no TLS groups are defined by the TLSIN file, then"
    print "    one group is created per chain."
    print


def main(path, opt_dict):
    struct = LoadStructure(fil = path)

    tls_group_list = []

    ## make the TLS groups
    if opt_dict.has_key("-t"):

        try:
            fil = open(opt_dict["-t"], "r")
        except IOError:
            print "ERROR: TLSIN File not found %s" % (opt_dict["-t"])
            sys.exit(-1)
        
        tilist = TLSInfoList()
        tilist.load_refmac_tlsout_file(fil)
        
        for tls_info in tilist:
            tls = tls_info.make_tls_group(struct)
            tls.tls_info = tls_info
            tls_group_list.append(tls)

    else:
        ## create one TLS group per chain by default
        for chain in struct.iter_chains():
            tls = TLSGroup()
            tls.tls_info  = TLSInfo()

            tls.tls_info.range_list = [
                (chain.chain_id, chain[0].fragment_id,
                 chain.chain_id, chain[-1].fragment_id)]

            for aa in chain.iter_amino_acids():
                for atm in chain.iter_atoms():
                    tls.append(atm)


    ## fit TLS groups and write output
    print "REFMAC"
    print

    for tls in tls_group_list:

        if len(tls)<20:
            print "NOT ENOUGH ATOMS IN TLSGROUP"
            continue

        tls.origin = tls.calc_centroid()
        tls.calc_TLS_least_squares_fit()

        T = tls.T
        L = tls.L * rad2deg2
        S = tls.S * rad2deg

        tls.tls_info.origin = tls.origin.copy()
        tls.tls_info.T = (T[0,0], T[1,1], T[2,2], T[0,1], T[0,2], T[1,2])
        tls.tls_info.L = (L[0,0], L[1,1], L[2,2], L[0,1], L[0,2], L[1,2])
        tls.tls_info.S = (
            S[1,1]-S[0,0], S[0,0]-S[2,2], S[0,1], S[0,2], S[1,2], S[1,0],
            S[2,0], S[2,1])

        print
        print tls.tls_info.refmac_description()

    ## write out a PDB file with 0.0 tempature factors for all
    ## atoms in TLS groups
    if opt_dict.has_key("-p"):

        ## dictionary of all atoms in TLS groups
        tls_atoms = {}
        for tls in tls_group_list:
            for atm in tls:
                tls_atoms[atm] = True

        ## set the temp_factor of TLS atoms to 0.0
        for atm in struct.iter_all_atoms():
            if tls_atoms.has_key(atm):

                for aatm in atm.iter_alt_loc():
                    aatm.temp_factor = 0.0


        ## save the struct
        SaveStructure(fil=opt_dict["-p"], struct=struct)

    print


if __name__ == "__main__":
    try:
        (opts, args) = getopt.getopt(sys.argv[1:], "t:p:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    opt_dict = {}
    for (flag, data) in opts:
        opt_dict[flag] = data

    try:
        path = args[0]
    except IndexError:
        usage()
        sys.exit(1)

    main(path, opt_dict)
