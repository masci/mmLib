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
    print "refmac_tls_prep.py [-s ] [-t <tlsin>] [-o <tlsout>]"
    print "                   [-p <pdbout>] <structure file>"
    print
    print "description:"
    print "    Performs a least squares fit of TLS tensors to the"
    print "    tempature factors of the given structure file.  If"
    print "    no TLS groups are defined by the TLSIN file, then"
    print "    one group is created per chain."
    print
    print "    -t <tlsin>  Use the TLSIN file for defining groups and"
    print "                calculating atomic ADPs.  If no tensors exist"
    print "                in the group, then a least squares fit"
    print "                is automatically preformed."
    print
    print
    print "    -o <tlsout> Write a TLSOUT file with the LSQ fit tensors."
    print
    print "    -p <pdbout> Write a PDB file to the given path containing"
    print "                the TLS predicted anisotropic tempature factors"
    print "                unless the -s option is given (see below)."
    print
    print "    -s          Set temperature factors of all TLS atoms to"
    print "                0.0, leaving all ADP magnitude in the TLS model."
    print

def main(path, opt_dict):
    struct = LoadStructure(fil = path)

    tls_group_list = []

    ## make the TLS groups
    if opt_dict.has_key("-t"):

        try:
            fil = open(opt_dict["-t"], "r")
        except IOError:
            print "[ERROR]: TLSIN File not found %s" % (opt_dict["-t"])
            sys.exit(-1)
        
        tls_file = TLSFile()
        tls_file.set_file_format(TLSFileFormatTLSOUT())
        tls_file.load(fil)

        for tls_desc in tls_file.tls_desc_list:
            tls = tls_desc.construct_tls_group_with_atoms(struct)
            tls.tls_desc = tls_desc
            tls_group_list.append(tls)

    else:
        ## create one TLS group per chain by default
        for chain in struct.iter_chains():
            if chain.count_amino_acids()<10:
                continue

            chain_id1 = chain.chain_id

            ## find beginning and ending amino acids
            frag_id1 = None
            frag_id2 = None

            for frag in chain.iter_amino_acids():
                if frag_id1==None:
                    frag_id1 = frag.fragment_id
                frag_id2 = frag.fragment_id

            tls_desc = TLSGroupDesc()
            tls_desc.add_range(chain_id1, frag_id1, chain_id1, frag_id2, "ALL")
            tls = tls_desc.construct_tls_group_with_atoms(struct)
            tls_group_list.append(tls)
            tls.tls_desc = tls_desc

            print "Creating TLS Group: %s" % (tls.name)

    ## fit TLS groups and write output
    tls_file = TLSFile()
    tls_file.set_file_format(TLSFileFormatTLSOUT())

    ## preform a LSQ fit if necessary
    for tls in tls_group_list:

        print "[TLS GROUP  %s]" % (tls.name)

        ## if the TLS group is null, then perform a LSQ-TLS fit
        if tls.is_null():
            print "Null Group: Running TLS-LSQ"
            
            if len(tls)<20:
                print "ERROR: Not Enough Atoms in TLS Group."
                continue
                
            tls.origin = tls.calc_centroid()
            lsq_residual = tls.calc_TLS_least_squares_fit()

        tls.tls_desc.set_tls_group(tls)
        tls_file.tls_desc_list.append(tls.tls_desc)

    if opt_dict.has_key("-s"):
        ## strip out all ADP magnitude from atoms
        for tls_group in tls_group_list:
            for atm, U in tls_group.iter_atm_Utls():
                assert min(eigenvalues(U))>0.0
                atm.temp_factor = 0.0
                atm.U = None

    else:
        ## shift some Uiso displacement from the TLS T tensor to the
        ## individual atoms
        for tls_group in tls_group_list:

            for atm, U in tls_group.iter_atm_Utls():
                assert min(eigenvalues(U))>0.0

            ## leave some B magnitude in the file for refinement
            (tevals, R) = eigenvectors(tls_group.T)
            tmin = min(tevals)
            T = matrixmultiply(R, matrixmultiply(tls_group.T, transpose(R)))
            T = T - (tmin * identity(3, Float))
            tls_group.T = matrixmultiply(transpose(R), matrixmultiply(T, R))

            bmin = U2B * tmin
            for atm, U in tls_group.iter_atm_Utls():
                atm.temp_factor = bmin
                atm.U = None

                assert min(eigenvalues(U + (B2U*bmin*identity(3, Float))))>0.0

    ## write TLSOUT file with new tensor values
    if opt_dict.has_key("-o"):
        for tls_group in tls_group_list:
            tls_group.tls_desc.set_tls_group(tls_group)
        
        print "Saving TLSIN: %s" % (opt_dict["-o"])
        tls_file.save(open(opt_dict["-o"], "w"))

    ## write out a PDB file with 0.0 tempature factors for all
    ## atoms in TLS groups
    if opt_dict.has_key("-p"):
        ## save the struct
        print "Saving XYZIN: %s" % (opt_dict["-p"])
        SaveStructure(fil=opt_dict["-p"], struct=struct)


if __name__ == "__main__":
    try:
        (opts, args) = getopt.getopt(sys.argv[1:], "st:p:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    opt_dict = {}
    for (flag, data) in opts:
        opt_dict[flag] = data

    try:
        pdb_path = args[0]
    except IndexError:
        usage()
        sys.exit(1)


    main(pdb_path, opt_dict)
