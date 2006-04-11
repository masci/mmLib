## TLS Motion Determination (TLSMD)
## Copyright 2002-2006 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import sys
import math
import numpy
import gc

from mmLib import Constants, FileIO, TLS, Structure

import misc
import const
import conf
import console
import atom_selection
import tlsmdmodule
import opt_containers
import adp_smoothing
import independent_segment_opt
import cpartition_recombination
import html


def TLSMD_Main(struct_file_path  = None,
               sel_chain_ids     = None,
               html_report_dir   = None):

    ## create the analysis processor and load the structure, select chains
    analysis = TLSMDAnalysis(
        struct_file_path    = struct_file_path,
        sel_chain_ids       = sel_chain_ids,
        struct2_file_path   = conf.globalconf.target_struct_path,
        struct2_chain_id    = conf.globalconf.target_struct_chain_id)
    
    IndependentTLSSegmentOptimization(analysis)
    RecombineIndependentTLSSegments(analysis)

    if analysis.struct2_file_path is not None and analysis.struct2_chain_id is not None:
        SumperimposeHomologousStructure(analysis)

    if html_report_dir is not None and analysis.num_chains() > 0:
        FitConstrainedTLSModel(analysis)
        report = html.HTMLReport(analysis)
        report.write(html_report_dir)


class TLSMDAnalysis(object):
    """Central object for a whole-structure TLS analysis.
    """
    def __init__(self,
                 struct_file_path  = None,
                 struct2_file_path = None,
                 struct2_chain_id  = None,
                 sel_chain_ids     = None):

        conf.globalconf.prnt()

        self.struct_file_path = struct_file_path
        self.struct2_file_path = struct2_file_path
        self.struct2_chain_id = struct2_chain_id
        self.target_chain = None

        if isinstance(sel_chain_ids, str):
            self.sel_chain_ids = sel_chain_ids.split(",")
        elif isinstance(sel_chain_ids, list):
            self.sel_chain_ids = sel_chain_ids
        else:
            self.sel_chain_ids = None

        self.struct = LoadStructure(self.struct_file_path)
        self.struct_id = self.struct.structure_id
        self.chains = None

        self.select_chains()
        self.prnt_settings()

    def prnt_settings(self):
        chain_ids = []
        for chain in self.chains:
            chain_ids.append(chain.chain_id)
        cids = ",".join(chain_ids)

        console.kvformat("STRUCTURE FILE", self.struct_file_path)
        console.kvformat("STRUCTURE ID", self.struct_id)
        console.kvformat("CHAIN IDs SELECTED FOR ANALYSIS", cids)
        console.endln()
        
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
            naa = chain.count_amino_acids()
            nna = chain.count_nucleic_acids()
            num_frags = max(naa, nna)
            if num_frags < 10:
                continue
            
            segment = ConstructSegmentForAnalysis(chain)
            segments.append(segment)
            segment.struct = self.struct
        
        self.chains = segments

    def iter_chains(self):
        return iter(self.chains)

    def num_chains(self):
        return len(self.chains)

    def get_chain(self, chain_id):
        for chain in self.iter_chains():
            if chain.chain_id == chain_id:
                return chain
        return None


def LoadStructure(struct_file_path):
    """Loads Structure, chooses a unique struct_id string.
    Also, search the REMARK records for TLS group records.  If they
    are found, then add the TLS group ADP magnitude to the B facors of
    the ATOM records.
    """
    console.kvformat("LOADING STRUCTURE", struct_file_path)

    ## load struct
    struct = FileIO.LoadStructure(
        fil = struct_file_path, distance_bonds = True)

    console.kvformat("HEADER", struct.header)
    console.kvformat("TITLE", struct.title)
    console.kvformat("EXPERIMENTAL METHOD", struct.experimental_method)
    
    ## set the structure ID
    if conf.globalconf.struct_id is not None:
        struct_id = conf.globalconf.struct_id
    else:
        struct_id = struct.structure_id
        conf.globalconf.struct_id = struct_id
    struct.structure_id = struct_id

    console.endln()

    ## if there are REFMAC5 TLS groups in the REMARK records of
    ## the PDB file, then add those in
    tls_file = TLS.TLSFile()
    tls_file.set_file_format(TLS.TLSFileFormatPDB())

    fil = open(struct_file_path, "r")
    tls_file.load(fil)

    if len(tls_file.tls_desc_list) > 0:
        console.stdoutln("ADDING TLS GROUP Bequiv TO ATOM TEMPERATURE FACTORS")
        console.stdoutln("    NUM TLS GROUPS: %d" % (len(tls_file.tls_desc_list)))

        ## assume REFMAC5 groups where Utotal = Utls + Biso(temp_factor)
        for tls_desc in tls_file.tls_desc_list:
            tls_group = tls_desc.construct_tls_group_with_atoms(struct)
            console.stdoutln("    TLS GROUP: %s" % (tls_group.name))
            for atm, Utls in tls_group.iter_atm_Utls():
                bresi = atm.temp_factor
                atm.temp_factor = bresi + (Constants.U2B * numpy.trace(Utls) / 3.0)
                atm.U = (Constants.B2U * bresi * numpy.identity(3, float)) + Utls

        console.endln()

    return struct

        
def ConstructSegmentForAnalysis(raw_chain):
    """Returns a list of Segment instance from the
    Chain instance which is properly modified for use in
    the this application.
    """
    ## Sets that atm.include attribute for each atom in the chains
    ## being analyzed by tlsmd
    for atm in raw_chain.iter_all_atoms():
        atm.include = atom_selection.calc_include_atom(atm)

    ## ok, use the chain but use a segment and cut off
    ## any leading and trailing non-amino acid residues
    ## do not include a fragment with no included atoms
    naa = raw_chain.count_amino_acids()
    nna = raw_chain.count_nucleic_acids()

    if naa > nna:
        iter_residues = raw_chain.iter_amino_acids()
    elif nna > naa:
        iter_residues = raw_chain.iter_nucleic_acids()
        
    segment = Structure.Segment(chain_id = raw_chain.chain_id)
    for frag in iter_residues:
        for atm in frag.iter_all_atoms():
            if atm.include:
                segment.add_fragment(frag)
                break

    ## apply data smooth if desired
    if conf.globalconf.adp_smoothing > 0:
        adp_smoothing.IsotropicADPDataSmoother(segment, conf.globalconf.adp_smoothing)

    ## this is useful: for each fragment in the minimization
    ## set a attribute for its index position
    for i, frag in enumerate(segment.iter_fragments()):
        frag.ifrag = i

    ## create a TLSModelAnalyzer instance for the chain, and
    ## attach the instance to the chain for use by the rest of the
    ## program
    segment.tls_analyzer = tlsmdmodule.TLSModelAnalyzer()
    xlist = atom_selection.chain_to_xmlrpc_list(segment.iter_all_atoms())
    segment.tls_analyzer.set_xmlrpc_chain(xlist)

    return segment


def IndependentTLSSegmentOptimization(analysis):
    """Performs the TLS graph minimization on all TLSGraphs.
    """
    for chain in analysis.chains:
        isopt = independent_segment_opt.ISOptimization(
            chain,
            conf.globalconf.min_subsegment_size,
            conf.globalconf.nparts)

        isopt.run_minimization()
        if not isopt.minimized:
            continue

        console.endln()
        console.stdoutln("="*79)
        console.stdoutln("MINIMIZING CHAIN %s" % (chain))
        isopt.prnt_detailed_paths()

        chain.partition_collection = isopt.construct_partition_collection(conf.globalconf.nparts)
        chain.partition_collection.struct = analysis.struct
        chain.partition_collection.struct_file_path = analysis.struct_file_path


def RecombineIndependentTLSSegments(analysis):
    console.endln()
    console.stdoutln("TLS SEGMENT RECOMBINATION")
    for chain in analysis.chains:
        cpartition_recombination.ChainPartitionRecombinationOptimization(chain)


def FitConstrainedTLSModel(analysis):
    """
    """
    console.endln()
    console.stdoutln("CALCULATING CONSTRAINED TLS MODEL FOR VISUALIZATION")
    for chain in analysis.iter_chains():
        console.stdoutln("CHAIN %s" % (chain.chain_id))
        for cpartition in chain.partition_collection.iter_chain_partitions():
            console.stdoutln("TLS GROUPS: %d" % (cpartition.num_tls_segments()))
            for tls in cpartition.iter_tls_segments():
                tls.fit_to_chain(cpartition.chain)
    console.endln()


def SumperimposeHomologousStructure(analysis):
    """
    """
    import structcmp

    target_struct = FileIO.LoadStructure(fil = analysis.struct2_file_path)
    target_chain = target_struct.get_chain(analysis.struct2_chain_id)

    if target_chain is None:
        console.stderrln(
            "UNABLE TO LOAD TARGET STRUCTURE/CHAIN: %s:%s" % (
            target_struct, target_chain))
        return

    analysis.target_chain = target_chain

    for chain in analysis.iter_chains():
        console.endln()
        console.kvformat("Superimposing Chain", chain.chain_id)
        hyp = structcmp.TLSConformationPredctionHypothosis(chain, target_chain)
        for ntls, cpartition in chain.partition_collection.iter_ntls_chain_partitions():
            console.endln()
            console.stdoutln("Number of TLS Segments........: %d" % (ntls))
            hyp.add_conformation_prediction_to_chain_partition(cpartition)
