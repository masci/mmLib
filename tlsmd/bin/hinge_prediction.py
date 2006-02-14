#!/usr/bin/env python

## TLS Motion Determination (TLSMD)
## Copyright 2002 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import sys
import math
import numpy

from mmLib import Constants, FileLoader, Structure
from tlsmdlib import tlsmd_analysis, fit_engine, lineartls, nonlineartls


class TLSHingePredictionHypothosis(object):
    """Use a fixed width sliding window to detect hinge points in a protein chain.
    """

    def __init__(self, chain, residue_window_with):
        self.chain = chain
        self.residue_window_width = residue_window_with

        for atm in chain.iter_all_atoms():
            atm.include = True
    
    def hinge_analysis(self, data_filename):
        fil = open(data_filename, "w")
        xchain = fit_engine.XChain(tlsmd_analysis.chain_to_xmlrpc_list(self.chain))

        tls_model = lineartls.LinearTLSModel()
        #tls_model = nonlineartls.NLTLSModel()

        tls_model.set_xmlrpc_chain(xchain.xmlrpc_chain)

        win = self.residue_window_width
        ifrag = 0
        ifrag_end = len(self.chain) - 2*win + 1

        for ifrag in range(ifrag, ifrag_end):

            ifrag1a = ifrag
            ifrag2a = ifrag+win-1

            ifrag1b = ifrag+win   
            ifrag2b = ifrag+2*win-1
            
            frag_id1a = self.chain[ifrag1a].fragment_id
            frag_id2a = self.chain[ifrag2a].fragment_id
            frag_id1b = self.chain[ifrag1b].fragment_id
            frag_id2b = self.chain[ifrag2b].fragment_id

            istarta = xchain.get_istart(frag_id1a)
            ienda   = xchain.get_iend(frag_id2a)
            istartb = xchain.get_istart(frag_id1b)
            iendb   = xchain.get_iend(frag_id2b)

            print "HINGE WINDOW: %s %s-%s:%s-%s" % (
                self.chain.chain_id, frag_id1a, frag_id2a, frag_id1b, frag_id2b)
            
            hdict = tls_model.calc_isotropic_hinge_delta(istarta, ienda, istartb, iendb)

            msd_ab  = hdict["hdelta_ab"]
            msd_abo = hdict["hdelta_abo"]
            msd_Lab = hdict["msd_c"]

            print "SEGMENT A(%d-%d): msd=%f rmsd=%f" % (
                istarta, ienda, Constants.U2B**2 * hdict["msd_a"],
                Constants.U2B * math.sqrt(hdict["msd_a"]))
            print "SEGMENT B(%d-%d): msd=%f rmsd=%f" % (
                istartb, iendb, Constants.U2B**2 * hdict["msd_b"],
                Constants.U2B * math.sqrt(hdict["msd_b"]))
            print "HINGE VALS: msd_ab=%f rmsd_ab=%f" % (
                Constants.U2B**2 * msd_ab, Constants.U2B * math.sqrt(msd_ab)) 
            print "HINGE VALS: msd_abo=%f rmsd_abo=%f" % (
                Constants.U2B**2 * msd_abo, Constants.U2B * math.sqrt(msd_abo))
            print "L TENSOR: msd_Lab=%f rmsd_Lab=%f" % (
                Constants.RAD2DEG2**2 *  msd_Lab, Constants.RAD2DEG2 * math.sqrt(msd_Lab))
            print

            fil.write("%s %f %f\n" % (
                frag_id2a, Constants.U2B**2 * msd_abo, Constants.U2B * math.sqrt(msd_abo)))

        fil.close()

def usage():
    print "hinge_prediction.py <structure file> <chain ID> <data output file>"
    sys.exit(-1)

def main():
    try:
        struct_file = sys.argv[1]
        chain_id = sys.argv[2]
        plot_file = sys.argv[3]
    except IndexError:
        usage()

    struct = FileLoader.LoadStructure(fil = struct_file)

    chain = struct.get_chain(chain_id)
    if chain == None:
        print "No such chain ID in structure"
        sys.exit(-1)

    segment = Structure.Segment(chain_id = chain_id)
    for frag in chain.iter_amino_acids():
        segment.add_fragment(frag)
    
    hinge_hyp = TLSHingePredictionHypothosis(segment, 14)
    hinge_hyp.hinge_analysis(plot_file)

if __name__ == "__main__":
    main()
