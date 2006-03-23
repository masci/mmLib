## TLS Motion Determination (TLSMD)
## Copyright 2002-2006 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import sys
import math
import numpy
import gc

from mmLib import Constants, FileIO, TLS

import misc
import const
import conf
import atom_selection
import tlsmdmodule
import opt_containers
import independent_segment_opt


class TLSMDAnalysis(object):
    """Central object for a whole-structure TLS analysis.
    """
    def __init__(self,
                 struct_file_path  = None,
                 struct2_file_path = None,
                 struct2_chain_id  = None,
                 sel_chain_ids     = None,
                 tlsdb_file        = None,
                 tlsdb_complete    = False):

        conf.globalconf.prnt()

        self.struct_file_path     = struct_file_path
        self.struct2_file_path    = struct2_file_path
        self.struct2_chain_id     = struct2_chain_id
        self.target_chain         = None

        if sel_chain_ids is not None:
            self.sel_chain_ids = sel_chain_ids.split(",")
        else:
            self.sel_chain_ids = None

        self.tlsdb_file      = tlsdb_file
        self.tlsdb_complete  = tlsdb_complete

        self.struct          = None
        self.struct_id       = None
        self.chains          = None

    def run_optimization(self):
        """Run the TLSMD optimization on the structure.
        """
        self.load_struct()

        ## auto name of tlsdb file then open
        if self.tlsdb_file is None:
            self.tlsdb_file = "%s_%s_%s.db" % (
                self.struct_id, conf.globalconf.tls_model,
                conf.globalconf.weight_model)

        self.select_chains()
        self.prnt_settings()
        self.calc_chain_minimization()
        #self.calc_chain_partition_recombination()
        
        if self.struct2_file_path != None and self.struct2_chain_id != None:
            self.calc_chain_partition_superposition()

    def prnt_settings(self):
        chain_ids = []
        for chain in self.chains:
            chain_ids.append(chain.chain_id)
        cids = ",".join(chain_ids)
        
        print "STRUCTURE FILE.....................: %s" % (self.struct_file_path)
        print "STRUCTURE ID.......................: %s" % (self.struct_id)
        print "CHAIN IDs SELECTED FOR ANALYSIS....: %s" % (cids)
        print "DATABASE FILE PATH.................: %s" % (self.tlsdb_file)
        print
        
    def load_struct(self):
        """Loads Structure, chooses a unique struct_id string.
        """
        print "LOADING STRUCTURE..................: %s" % (self.struct_file_path)

        ## load struct
        self.struct = FileIO.LoadStructure(
            fil = self.struct_file_path,
            build_properties = ("library_bonds","distance_bonds"))

        ## set the structure ID
        if conf.globalconf.struct_id != None:
            struct_id = conf.globalconf.struct_id
        else:
            struct_id = self.struct.structure_id
            conf.globalconf.struct_id = struct_id

        self.struct.structure_id = struct_id
        self.struct_id = struct_id

        print

        ## if there are REFMAC5 TLS groups in the REMARK records of
        ## the PDB file, then add those in
        tls_file = TLS.TLSFile()
        tls_file.set_file_format(TLS.TLSFileFormatPDB())

        fil = open(self.struct_file_path, "r")
        tls_file.load(fil)

        if len(tls_file.tls_desc_list)>0:
            print "MERGING TLS GROUPS"
            print "    NUM TLS GROUPS: %d" % (len(tls_file.tls_desc_list))

            ## assume REFMAC5 groups where Utotal = Utls + Biso(temp_factor)
            for tls_desc in tls_file.tls_desc_list:
                tls_group = tls_desc.construct_tls_group_with_atoms(self.struct)

                print "    TLS GROUP: %s" % (tls_group.name)
                
                for atm, Utls in tls_group.iter_atm_Utls():
                    bresi = atm.temp_factor
                    atm.temp_factor = bresi + (Constants.U2B * numpy.trace(Utls) / 3.0)
                    atm.U = (Constants.B2U * bresi * numpy.identity(3, float)) + Utls
            
            print

    def select_chains(self):
        """Selects chains for analysis.
        """
        ## select viable chains for TLS analysis
        segments = []
        
        for chain in self.struct.iter_chains():
            ## if self.sel_chain_ids is set, then only use those
            ## selected chain ids
            if self.sel_chain_ids is not None:
                if chain.chain_id not in self.sel_chain_ids:
                    continue

            ## count the number of amino acid residues in the chain
            if chain.count_amino_acids() < 10:
                continue

            ## ok, use the chain but use a segment and cut off
            ## any leading and trailing non-amino acid residues
            frag_id1 = None
            for aa in chain.iter_amino_acids():
                if frag_id1 is None:
                    frag_id1 = aa.fragment_id
                frag_id2 = aa.fragment_id
                
            segment = chain[frag_id1:frag_id2]
            segments.append(segment)

            ## Sets that atm.include attribute for each atom in the chains
            ## being analyzed by tlsmd
            for atm in segment.iter_all_atoms():
                atm.include = atom_selection.calc_include_atom(atm)
            
            segment.struct = self.struct
            segment.tls_analyzer = tlsmdmodule.TLSModelAnalyzer()
            xlist = atom_selection.chain_to_xmlrpc_list(segment)
            segment.tls_analyzer.set_xmlrpc_chain(xlist)
        
        self.chains = segments

    def calc_chain_minimization(self):
        """Performs the TLS graph minimization on all TLSGraphs.
        """
        for chain in self.chains:
            isopt = independent_segment_opt.ISOptimization(
                chain,
                conf.globalconf.min_subsegment_size,
                conf.globalconf.nparts)
            
            isopt.run_minimization()
            if not isopt.minimized:
                continue

            print
            print "="*79
            print "MINIMIZING CHAIN %s" % (chain)
            isopt.prnt_detailed_paths()

            chain.partition_collection = isopt.construct_partition_collection(conf.globalconf.nparts)
            chain.partition_collection.struct = self.struct
            chain.partition_collection.struct_file_path = self.struct_file_path
            
    def calc_chain_partition_recombination(self):
        print
        print "TLS SEGMENT RECOMBINATION: N INTO N-1"

        for chain in self.chains:

            improvements = []
            prev_cpartition = None
            for ntls, cpartition in chain.partition_collection.iter_ntls_chain_partitions():
                if prev_cpartition is None:
                    prev_cpartition = cpartition
                    continue

                print
                print "%d->%d Segment Recombination beat(%8.6f)" % (ntls, ntls-1, prev_cpartition.residual())
                combined_cp = ChainPartitionRecombination(cpartition)
                print "%s %8.6f" % (combined_cp, combined_cp.residual())
                
                if combined_cp.residual() < prev_cpartition.residual():
                    improvements.append(combined_cp)

                prev_cpartition = cpartition

            ## insert replacement ChainPartitions
            for cpartition in improvements:
                chain.partition_collection.insert_chain_partition(cpartition)

                
    def calc_chain_partition_superposition(self):
        import structcmp

        target_struct = FileIO.LoadStructure(fil = self.struct2_file_path)
        target_chain = target_struct.get_chain(self.struct2_chain_id)

        if target_chain == None:
            print "UNABLE TO LOAD TARGET STRUCTURE/CHAIN: %s:%s" % (
                target_struct, target_chain)
            return

        self.target_chain = target_chain
        
        for chain in self.iter_chains():
            print
            print "Superimposing Chain.............%s" % (chain.chain_id)
            hyp = structcmp.TLSConformationPredctionHypothosis(chain, target_chain)
            for ntls, cpartition in chain.partition_collection.iter_ntls_chain_partitions():
                print
                print "Number of TLS Segments........: %d" % (ntls)
                hyp.add_conformation_prediction_to_chain_partition(cpartition)

    def iter_chains(self):
        return iter(self.chains)

    def num_chains(self):
        return len(self.chains)
