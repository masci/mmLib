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
import adp_smoothing
import independent_segment_opt
import cpartition_recombination


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
        self.calc_chain_partition_recombination()
        self.calc_visualization_tls_models()

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
            fil = self.struct_file_path, distance_bonds = True)

        print "HEADER..............................: %s" % (self.struct.header)
        print "TITLE...............................: %s" % (self.struct.title)
        print "EXPERIMENTAL METHOD.................: %s" % (self.struct.experimental_method)

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
            print "ADDING TLS GROUP Bequiv TO ATOM TEMPERATURE FACTORS"
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
            
            segment = ConstructSegmentForAnalysis(chain)
            segments.append(segment)
            segment.struct = self.struct
        
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
        if not conf.globalconf.recombination:
            return

        print
        print "TLS SEGMENT RECOMBINATION"
        for chain in self.chains:
            cpartition_recombination.ChainPartitionRecombinationOptimization(chain)

    def calc_visualization_tls_models(self):
        print
        print "CALCULATING CONSTRAINED TLS MODEL FOR VISUALIZATION"
        for chain in self.chains:
            print "CHAIN %s" % (chain.chain_id)

            for cpartition in chain.partition_collection.iter_chain_partitions():
                print "TLS GROUPS: %d" % (cpartition.num_tls_segments())
                
                for tls in cpartition.iter_tls_segments():
                    tls.fit_to_chain(cpartition.chain)
        print
            
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


def ConstructSegmentForAnalysis(raw_chain):
    """Returns a list of Segment instance from the
    Chain instance which is properly modified for use in
    the this application.
    """
    ## ok, use the chain but use a segment and cut off
    ## any leading and trailing non-amino acid residues
    frag_id1 = None
    for aa in raw_chain.iter_amino_acids():
        if frag_id1 is None:
            frag_id1 = aa.fragment_id
        frag_id2 = aa.fragment_id

    segment = raw_chain[frag_id1:frag_id2]
    print "SELECTING RANGE %s-%s; GOT RANGE %s-%s" % (frag_id1, frag_id2, segment[0].fragment_id, segment[-1].fragment_id)
    
    ## this is useful: for each fragment in the minimization
    ## set a attribute for its index position
    for i, frag in enumerate(segment.iter_fragments()):
        frag.ifrag = i

    ## Sets that atm.include attribute for each atom in the chains
    ## being analyzed by tlsmd
    for atm in segment.iter_all_atoms():
        atm.include = atom_selection.calc_include_atom(atm)

    ## apply data smooth if desired
    if conf.globalconf.adp_smoothing > 0:
        adp_smoothing.IsotropicADPDataSmoother(segment, conf.globalconf.adp_smoothing)

    ## create a TLSModelAnalyzer instance for the chain, and
    ## attach the instance to the chain for use by the rest of the
    ## program
    segment.tls_analyzer = tlsmdmodule.TLSModelAnalyzer()
    xlist = atom_selection.chain_to_xmlrpc_list(segment.iter_all_atoms())
    segment.tls_analyzer.set_xmlrpc_chain(xlist)

    return segment

