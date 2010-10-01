## TLS Motion Determination (TLSMD)
## Copyright 2002-2010 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
##
## Called by:
##	tlsmdlib/html.py
##		class HTMLReport()
##			init_colors()
##			write_cwd()
##			raster3d_render_tls_graph_path()
##			write_multi_chain_alignment()
##			html_multi_chain_alignment()
##			html_refinement_prep()

## Python modules
import sys
import math
import numpy
import os
import time

## pymmlib
from mmLib import Constants, FileIO, TLS, Structure

## TLSMD
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

import signal
def SIGUSR1_handler(signum, frame):
    ## EAM FIXME - This probably should be defined at a higher level and be made
    ## available throughout the package
    raise RuntimeError, 'Caught external SIGUSR1'

def TLSMD_Main(struct_file_path  = None,
               sel_chain_ids     = None,
               html_report_dir   = None):

    console.cpu_time_stdoutln("->Starting TLSMD analysis: %s" % time.clock())

    ## create the analysis processor and load the structure, select chains
    analysis = TLSMDAnalysis(
        struct_file_path    = struct_file_path,
        sel_chain_ids       = sel_chain_ids,
        struct2_file_path   = conf.globalconf.target_struct_path,
        struct2_chain_id    = conf.globalconf.target_struct_chain_id)
    console.cpu_time_stdoutln("->LoadStructure: %s" % time.clock())

    IndependentTLSSegmentOptimization(analysis)
    console.cpu_time_stdoutln("->IndependentTLSSegmentOptimization: %s" % (
        time.clock()))

    RecombineIndependentTLSSegments(analysis)
    console.cpu_time_stdoutln("->RecombineIndependentTLSSegments: %s" % (
        time.clock()))

    if analysis.struct2_file_path is not None and \
       analysis.struct2_chain_id is not None:
        console.stdoutln("SUPERIMPOSING HOMOLOGOUS STRUCTURE")
        SuperimposeHomologousStructure(analysis)

    if html_report_dir is not None and analysis.num_chains() > 0:
        FitConstrainedTLSModel(analysis)

        ## generate summary page. Allows the user to see some information
        ## before the analysis is fully complete
        console.stdoutln("GENERATING SUMMARY PAGE")
        summary = html.HTMLSummaryReport(analysis)
        summary.write_summary(html_report_dir)
        console.cpu_time_stdoutln("->HTMLSummaryReport: %s" % time.clock())

        if not conf.globalconf.skip_html:
            ## generate in-depth HTML report pages
            console.stdoutln("GENERATING ANALYSIS/REPORT HTML PAGES")
            report = html.HTMLReport(analysis)
            report.write(html_report_dir)

        ## TODO: Figure out a way to gzip all TXT files once TLSMD is done
        ## with them. 2010-09-24
        #import fnmatch
        #import gzip
        #abs_report_path = "%s/%s/%s" % (
        #    conf.TLSMD_WORK_DIR, conf.globalconf.job_id, html_report_dir)
        #console.stdoutln("GZIP: Attempting to gzip all TXT files in: %s" % abs_report_path)
        #for path, dirs, files in os.walk(abs_report_path):
        #    console.stdoutln("GZIP: report_dir = %s" % abs_report_path)
        #    for f in files:
        #        if fnmatch.fnmatch(f, '*.txt'):
        #            console.stdoutln("GZIP: found a TXT file: %s" % f)
        #            ## If we find a TXT file, gzip it
        #            f_in = open('%s/%s' % (abs_report_path, f), 'rb')
        #            f_out = gzip.open('%s/%s.gz' % (abs_report_path, f), 'wb')
        #            f_out.writelines(f_in)
        #            f_out.close()
        #            f_in.close()
        #            os.remove(f)

    console.cpu_time_stdoutln("->TLSMD::Total_time: %s" % time.clock())


class TLSMDAnalysis(object):
    """Central object for a whole-structure TLS analysis.
    """
    def __init__(self,
                 struct             = None,
                 struct_file_path   = None,
                 struct_file_object = None,
                 struct2_file_path  = None,
                 struct2_chain_id   = None,
                 sel_chain_ids      = None):

        conf.globalconf.prnt()

        self.struct2_file_path = struct2_file_path
        self.struct2_chain_id = struct2_chain_id
        self.target_chain = None

        if isinstance(sel_chain_ids, str):
            self.sel_chain_ids = sel_chain_ids.split(",")
        elif isinstance(sel_chain_ids, list):
            self.sel_chain_ids = sel_chain_ids
        else:
            self.sel_chain_ids = None

        if struct_file_path is not None:
            self.struct = LoadStructure(struct_file_path)
        elif struct_file_object is not None:
            self.struct = LoadStructure(struct_file_object)
        elif struct is not None:
            self.struct = struct

        self.struct_id = self.struct.structure_id
        self.chains = None

        self.select_chains()
        self.prnt_settings()

    def prnt_settings(self):
        chain_ids = []
        for chain in self.chains:
            chain_ids.append(chain.chain_id)
        cids = ",".join(chain_ids)

        console.debug_stdoutln(">tlsmd_analysis->TLSMDAnalysis()")
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

            ## count the number of amino acid and/or nucleic acid residues in
            ## the chain and skip those that are too small
            ## also, warn the user if a chain contains non-standard residue
            ## types
            naa = chain.count_amino_acids()
            nna = chain.count_nucleic_acids()
            non = chain.count_non_standard_residues()

            if nna == 0 and (naa > 0 and naa < conf.MIN_AMINO_PER_CHAIN):
                console.kvformat("SKIPPING SMALL AMINO ACID CHAIN", 
                    chain.chain_id)
                continue

            if naa == 0 and (nna > 0 and nna < conf.MIN_NUCLEIC_PER_CHAIN):
                console.kvformat("SKIPPING SMALL NUCLEIC ACID CHAIN", 
                    chain.chain_id)
                continue

            if naa > nna and nna > 0:
                console.kvformat("CHAIN WITH DIFFERENT RESIDUE TYPES",
                    chain.chain_id)

            if non > 0:
                console.kvformat("CHAIN WITH NON-STANDARD RESIDUE TYPES",
                    chain.chain_id)

            num_frags = max(naa, nna)

            segment = ConstructSegmentForAnalysis(chain)
            if segment == "":
                console.kvformat("SKIPPING CHAIN WITH NO AA/NA RESIDUES",
                    chain.chain_id)
                continue

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


def LoadStructure(struct_source):
    """Loads Structure, chooses a unique struct_id string.
    Also, search the REMARK records for TLS group records. If they
    are found, then add the TLS group ADP magnitude to the B facors of
    the ATOM records.
    """
    console.debug_stdoutln(">tlsmd_analysis->LoadStructure()")

    ## determine the argument type
    if isinstance(struct_source, str):
        file_path = struct_source
        console.kvformat("LOADING STRUCTURE", file_path)
        fobj = open(file_path, "r")
    elif hasattr(struct_source, "__iter__") and hasattr(struct_source, "seek"):
        console.kvformat("LOADING STRUCTURE", str(struct_source))
        fobj = struct_source
    else:
        raise ValueError

    ## TODO: Is this the same load struct used during the sanity checks? 2009-05-28
    ## load struct
    struct = FileIO.LoadStructure(file = fobj, distance_bonds = True)
    job_dir = str(os.path.dirname(str(struct_source)))

    console.kvformat("HEADER", struct.header)
    console.kvformat("TITLE", struct.title)
    console.kvformat("EXPERIMENTAL METHOD", struct.experimental_method)
    console.kvformat("PATH", job_dir)

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

    ## return to the beginning of the file and read the REMARK/TLS records
    fobj.seek(0)
    tls_file.load(fobj)

    if len(tls_file.tls_desc_list) > 0:
        console.stdoutln("ADDING TLS GROUP Bequiv TO ATOM TEMPERATURE FACTORS")
        console.stdoutln("    NUM TLS GROUPS: %d" % (
            len(tls_file.tls_desc_list)))

        ## assume REFMAC5 groups where Utotal = Utls + Biso(temp_factor)
        for tls_desc in tls_file.tls_desc_list:
            try:
                tls_group = tls_desc.construct_tls_group_with_atoms(struct)
            except:
                print console.formatExceptionInfo()

            console.stdoutln("    TLS GROUP: %s" % (tls_group.name))
            for atm, Utls in tls_group.iter_atm_Utls():
                bresi = atm.temp_factor
                atm.temp_factor = bresi + (Constants.U2B * numpy.trace(Utls) / 3.0)
                atm.U = (Constants.B2U * bresi * numpy.identity(3, float)) + Utls

        console.endln()

    return struct


def ConstructSegmentForAnalysis(raw_chain):
    """Returns a list of Segment instance from the Chain instance which is 
    properly modified for use in the this application.
    """
    console.debug_stdoutln(">tlsmd_analysis->ConstructSegmentForAnalysis(chain %s)" % (
        raw_chain.chain_id))

    ## NOTE: raw_chain = "Chain(1:A, Res(MET,1,A)...Res(VAL,50,A))"
    ## Sets that atm.include attribute for each atom in the chains
    ## being analyzed by tlsmd
    for atm in raw_chain.iter_all_atoms():
        #atm.include = atom_selection.calc_include_atom(atm)
        atm.include = atom_selection.calc_include_atom(atm, reject_messages = True)

    ## ok, use the chain but use a segment and cut off
    ## any leading and trailing non-amino acid residues
    ## do not include a fragment with no included atoms
    naa = nna = ota = 0
    naa = raw_chain.count_amino_acids()
    nna = raw_chain.count_nucleic_acids()
    ota = raw_chain.count_fragments()

    if naa > nna:
        ## Probably a protein with (possibly) some nucleic acids.
        iter_residues = raw_chain.iter_amino_acids()
    elif nna > naa:
        ## Probably a nucleic acid with (possibly) some amino acids.
        iter_residues = raw_chain.iter_nucleic_acids()

    if naa == 0 and nna == 0 and ota > 0:
        ## This chain does not have any amino or nucleic acids, so skip.
        return ""

    segment = Structure.Segment(chain_id = raw_chain.chain_id)

    for frag in iter_residues:
        for atm in frag.iter_all_atoms():
            if atm.include:
                segment.add_fragment(frag)
                break

    ## apply data smooth if desired (default is "0")
    if conf.globalconf.adp_smoothing > 0:
        adp_smoothing.IsotropicADPDataSmoother(segment, 
                                               conf.globalconf.adp_smoothing)

    ## this is useful: for each fragment in the minimization
    ## set an attribute for its index position
    for i, frag in enumerate(segment.iter_fragments()):
        ## NOTE (by Christoph):
        ## Example output:
        ## 0 : Res(MET,1,A)
        ## 1 : Res(ILE,2,A)
        ## ...
        frag.ifrag = i

    ## create a TLSModelAnalyzer instance for the chain, and attach the
    ## instance to the chain for use by the rest of the program
    segment.tls_analyzer = tlsmdmodule.TLSModelAnalyzer()
    xlist = atom_selection.chain_to_xmlrpc_list(segment.iter_all_atoms())
    segment.tls_analyzer.set_xmlrpc_chain(xlist)

    ## INPUT : raw_chain = "Chain(1:A, Res(MET,1,A)...Res(VAL,50,A))"
    ## OUTPUT: segment   = "Segment(1:A, Res(MET,1,A)...Res(VAL,50,A))"
    return segment


def IndependentTLSSegmentOptimization(analysis):
    """Performs the TLS graph minimization on all TLSGraphs.
    """
    for chain in analysis.chains:
        isopt = independent_segment_opt.ISOptimization(
            chain,
            conf.globalconf.min_subsegment_size,
            conf.globalconf.nparts)

        ## TODO: Divide this into two CPU times, 2009-12-10
        #console.stdoutln("CPU_TIME ->ISOptResidualGraph: %s" % time.clock())

        isopt.run_minimization()
        if not isopt.minimized:
            continue

        console.endln()
        console.stdoutln("="*79)
        console.debug_stdoutln(">tlsmd_analysis->IndependentTLSSegmentOptimization()")
        console.stdoutln("MINIMIZING CHAIN %s" % (chain))
        isopt.prnt_detailed_paths()

        chain.partition_collection = isopt.construct_partition_collection(conf.globalconf.nparts)
        chain.partition_collection.struct = analysis.struct

def RecombineIndependentTLSSegments(analysis):
    console.endln()
    console.debug_stdoutln(">tlsmd_analysis->RecombineIndependentTLSSegments()")
    console.stdoutln("TLS SEGMENT RECOMBINATION")
    for chain in analysis.chains:
        ## E.g., chain="Segment(1:A, Res(ILE,16,A)...Res(SER,116,A))"
        cpartition_recombination.ChainPartitionRecombinationOptimization(chain)

def FitConstrainedTLSModel(analysis):
    """Calculates constrained TLS model for visualization.
    """
    console.endln()
    console.debug_stdoutln(">tlsmd_analysis->FitConstrainedTLSModel()")
    console.stdoutln("CALCULATING CONSTRAINED TLS MODEL FOR VISUALIZATION")

    ## EAM Feb 2008 User job was getting stuck in fit_to_chain()
    ## Obviously it would be nice to fix the actual error, but at least we would
    ## like to be able to give it a swift non-fatal kick by sending SIGUSR1
    signal.signal(signal.SIGUSR1, SIGUSR1_handler)

    ## Progress tracking 
    ##    - assume this portion of the run occupies 0.1 -> 0.5 of the total time
    progress = 0.1

    for chain in analysis.iter_chains():
        console.stdoutln("CHAIN %s" % (chain.chain_id))
        for cpartition in chain.partition_collection.iter_chain_partitions():
            ## cpartition.chain = "Segment(1:A, Res(MET,1,A)...Res(VAL,50,A))"

            console.stdoutln("TLS GROUPS: %d" % (cpartition.num_tls_segments()))

            for tls in cpartition.iter_tls_segments():
                try:
                    tls.fit_to_chain(cpartition.chain)

                    ## TODO: Write out data for residual plots.
                    #gp = gnuplots.LSQR_vs_TLS_Segments_Pre_Plot(cpartition.chain)
                    #console.stdoutln("FIT_TO_CHAIN_PATH: %s" % analysis.struct2_file_path)

                except (RuntimeError, numpy.linalg.linalg.LinAlgError), e:
                    msg  = "            Runtime error for [%s]: %s, " % (
                        tls, e)
                    msg += "trying to continue..."
                    console.stdoutln(msg)
                    print console.formatExceptionInfo()
                    pass

        ## Track progress
        progress += 0.4/analysis.num_chains()
        progress_report = open("progress","w+")
        print >> progress_report, progress
        ## progress_report.write(progress)
        progress_report.close()

    console.cpu_time_stdoutln("->FitConstrainedTLSModel: %s" % time.clock())
    console.endln()


def SuperimposeHomologousStructure(analysis):
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
