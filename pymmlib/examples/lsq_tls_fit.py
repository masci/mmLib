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
        
        tls_file = TLSFile()
        tls_file.set_file_format(TLSFileFormatTLSOUT())
        tls_file.load(fil)
        
        for tls_desc in tls_file.tls_desc_list:
            tls = tls_desc.generate_tls_group(struct)
            tls.tls_desc = tls_desc
            tls_group_list.append(tls)

    else:
        ## create one TLS group per chain by default
        for chain in struct.iter_chains():
            try:
                chain_id1 = chain.chain_id
                frag_id1  = chain[0].fragment_id
                frag_id2  = chain[-1].fragment_id
            except IndexError:
                continue
            
            tls_desc = TLSGroupDesc()
            tls_desc.add_range(chain_id1, frag_id1, chain_id1, frag_id2, "ALL")
            tls = tls_desc.generate_tls_group(struct)
            tls_group_list.append(tls)
            tls.tls_desc = tls_desc


    ## fit TLS groups and write output
    print "REFMAC"
    print

    for tls in tls_group_list:

        if len(tls)<20:
            print "NOT ENOUGH ATOMS IN TLSGROUP"
            continue

        tls.origin = tls.calc_centroid()
        tls.calc_TLS_least_squares_fit()
        tls.shift_COR()

        tls.tls_desc.set_tls_group(tls)

        tls_file = TLSFile()
        tls_file.set_file_format(TLSFileFormatTLSOUT())
        tls_file.tls_desc_list.append(tls.tls_desc)
        tls_file.save(sys.stdout)


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
