#!/usr/bin/env python

## TLS Motion Determination (TLSMD)
## Copyright 2002 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import sys
import math
import numpy

from mmLib import Constants, FileIO, Structure, TLS, AtomMath
from tlsmdlib import conf, atom_selection, tlsmdmodule, tlsmd_analysis, table, tls_calcs

def iter_fragment_atoms(fragiter):
    for frag in fragiter:
        for atm in frag.iter_all_atoms():
            if atm.include:
                yield atm

def symmtensor33diff(T1, T2):
    delta = numpy.trace(T1) - numpy.trace(T2)
    return delta**2

    msd_sum = 0.0
    for j in xrange(3):
        for i in xrange(j, 3):
            delta = T1[i, j] - T2[i, j]
            msd_sum += delta**2
    return msd_sum / 6

def ResidualInfo(chain, range, tlsdict):
    IT, IL, IS, IO = tls_calcs.isotlsdict2tensors(tlsdict)

    num_atoms = 0
    weight_sum = 0.0
    msd_sum = 0.0
    atomiter = iter_fragment_atoms(chain.iter_fragments(*range))
    for atm, uiso_tls in TLS.iter_itls_uiso(atomiter, IT, IL, IS, IO):
        num_atoms += 1
        delta = atm.temp_factor - (Constants.U2B * uiso_tls)
        msd_sum += atm.occupancy * delta**2
        weight_sum += atm.occupancy
    msd = msd_sum / weight_sum
    return msd
    
def CrossPrectionResidual(chain, range1, tlsdict1, range2, tlsdict2):
    """
    """
    msd1 = ResidualInfo(chain, range1, tlsdict1)
    msd2 = ResidualInfo(chain, range2, tlsdict2)
    msd1_predict2 = ResidualInfo(chain, range2, tlsdict1)
    msd2_predict1 = ResidualInfo(chain, range1, tlsdict2)
    return msd1, msd2, msd1_predict2, msd2_predict1


class TLSHingePredictionHypothosis(object):
    """Use a fixed width sliding window to detect hinge points in a protein chain.
    """
    def __init__(self, structfile, run_tlsmd = True):
        self.structfile = structfile
        self.run_tlsmd = run_tlsmd
        
        self.residue_window_width = 14
        self.hinge_analysis_column = 1
        self.tlsmd_start_column = 2

        self.tlsmd = tlsmd_analysis.TLSMDAnalysis(
            struct_file_path = self.structfile)

        self.chain = None
        self.tbl_row_index = None
        self.tbl = None

    def run(self, plotfileprefix):
        if self.run_tlsmd:
            self.tlsmd.run_optimization(True)
        
        for chain in self.tlsmd.iter_chains():
            self.set_chain(chain)
            
            if self.run_tlsmd:
                self.extract_tlsmd_partitions(chain)

            self.hinge_analysis(chain)

            plotfilename = "%s_CHAIN%s.txt" % (plotfileprefix, chain.chain_id)
            self.write_data_file(plotfilename)

    def set_chain(self, chain):
        """Set the chain to this instance's current chain.  Construct output
        tables to hold information on hinge calculations.
        """
        first_last_frag_id = (chain[0].fragment_id, chain[-1].fragment_id)

        if self.run_tlsmd:
            cols = 2 + conf.globalconf.nparts - 1
        else:
            cols = 2

        tbl_row_index = {}
        tbl = table.StringTable(len(chain), cols, "0.0")
        for i, frag in enumerate(chain.iter_fragments()):
            tbl[i, 0] = frag.fragment_id
            tbl_row_index[frag.fragment_id] = i
            if cols > 2:
                for j in xrange(2, cols):
                    tbl[i, j] = 0

        self.chain = chain
        self.tbl_row_index = tbl_row_index
        self.tbl = tbl
        self.first_last_frag_id = first_last_frag_id

    def extract_tlsmd_partitions(self, chain):
        """Run the TLSMD partitioning algorithm and add it to the
        output data table.
        """
        presult = []
        partition_collection = chain.partition_collection

        col = self.tlsmd_start_column
        for cpartition in partition_collection.iter_chain_partitions():
            if cpartition.num_tls_segments() == 1:
                continue
            for tls in cpartition.iter_tls_segments():
                for frag_id1, frag_id2 in tls.iter_segment_ranges():
                    for frag_id in (frag_id1, frag_id2):
                        if frag_id not in self.first_last_frag_id:
                            row = self.tbl_row_index[frag_id]
                            self.tbl[row, col] = 1
            col += 1

    def hinge_analysis(self, chain):
        atom_list = atom_selection.chain_to_xmlrpc_list(chain.iter_all_atoms())

        tls_model = tlsmdmodule.TLSModelAnalyzer()
        tls_model.set_xmlrpc_chain(atom_list)

        win = self.residue_window_width
        ifrag = 0
        ifrag_end = len(chain) - 2*win + 1

        residual_dict = {}

        for ibegin in xrange(ifrag, ifrag_end):

            ifrag1a = ibegin
            ifrag2a = ibegin + win - 1

            ifrag1b = ibegin + win
            ifrag2b = ibegin + 2*win - 1
            
            frag_id1a = chain[ifrag1a].fragment_id
            frag_id2a = chain[ifrag2a].fragment_id
            frag_id1b = chain[ifrag1b].fragment_id
            frag_id2b = chain[ifrag2b].fragment_id

            print "HINGE WINDOW: %s %s-%s:%s-%s" % (
                chain.chain_id, frag_id1a, frag_id2a, frag_id1b, frag_id2b)

            range1 = (frag_id1a, frag_id2a)
            range2 = (frag_id1b, frag_id2b)
            
            tlsdict1 = tls_model.constrained_isotropic_fit([range1])
            tlsdict2 = tls_model.constrained_isotropic_fit([range2])

            residual = self.cross_prediction_residual(
                chain, range1, tlsdict1, range2, tlsdict2)
            
            for frag_id in (frag_id2a, ):
                if residual_dict.has_key(frag_id):
                    residual_dict[frag_id] += residual
                else:
                    residual_dict[frag_id] = residual

        for i, frag in enumerate(chain.iter_fragments()):
            if residual_dict.has_key(frag.fragment_id):
                self.tbl[i, 1] = residual_dict[frag.fragment_id]

    def cross_prediction_residual(self, chain, range1, tlsdict1, range2, tlsdict2):
        msd1, msd2, msd1_predict2, msd2_predict1 = CrossPrectionResidual(
            chain, range1, tlsdict1, range2, tlsdict2)
        return msd1_predict2 + msd2_predict1

        IT1, IL1, IS1, IO1 = tls_calcs.isotlsdict2tensors(tlsdict1)
        IT2, IL2, IS2, IO2 = tls_calcs.isotlsdict2tensors(tlsdict2)
        return Constants.RAD2DEG2 * symmtensor33diff(IL1, IL2)

    def write_data_file(self, filename):
        open(filename, "w").write(str(self.tbl))
        

def usage():
    print "hinge_prediction.py <structure file> <chain ID> <data output file>"
    raise SystemExit

def main():
    try:
        structfile = sys.argv[1]
        plotfileprefix = sys.argv[2]
    except IndexError:
        usage()

    hinge_hyp = TLSHingePredictionHypothosis(structfile)
    hinge_hyp.run(plotfileprefix)
    
if __name__ == "__main__":
    main()
