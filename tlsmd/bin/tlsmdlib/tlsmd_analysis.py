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
import datafile
import hcsssp
import tls_calcs
import fit_engine
import tlsmdmodule

def calc_include_atom(atm, reject_messages = False):
    """Filter out atoms from the model which will cause problems or
    cont contribute to the TLS analysis.
    """
    if atm.position == None:
        return False

    if atm.occupancy < 0.1:
        if reject_messages == True:
            print "calc_include_atom(%s): rejected because of low occupancy" % (atm)
	return False
    
    if numpy.trace(atm.get_U()) <= const.TSMALL:
        if reject_messages == True:
            print "calc_include_atom(%s): rejected because of small Uiso magnitude " % (atm)
        return False

    elif conf.globalconf.include_atoms == "MAINCHAIN":
        if atm.name not in const.MAINCHAIN_ATOMS:
            if reject_messages == True:
                print "calc_include_atom(%s): rejected non-mainchain atom" % (atm)
            return False
    
    return True


def calc_atom_weight(atm):
    """Weight the least-squares fit according to this function.
    """
    assert atm.occupancy >= 0.0 and atm.occupancy <= 1.0
    return atm.occupancy


def calc_num_subsegments(n, m):
    """Calculates the number of possible subsegment for the chain of length n and minimum
    subsegment length m.
    """
    m = m-1
    return (n*(n+1))/2 - (m*n - (m*(m+1))/2 + m)


def iter_ij(num_vertex, min_len):
    """Iterates over the i,j vertex indexes defining the edges
    to be built for the graph.  num_vertex gives the number
    of consecutive vertices to create, min_span is the minimum
    number of residues (fragments) a edge should span.
    """
    for i in xrange(num_vertex):
        for j in xrange(i + min_len, num_vertex):
            yield i, j


def iter_chain_subsegment_descs(chain, min_len):
    """Iterate over all possible subsegments of the given Chain object
    with a minimum size of min_span fragments.  The segments are yielded
    as Python dictionaries containing a description of the subsegment.
    """
    num_vertex = len(chain) + 1    
    for vertex_i, vertex_j in iter_ij(num_vertex, min_len):
        frag_id1 = chain.fragment_list[vertex_i].fragment_id
        frag_id2 = chain.fragment_list[vertex_j-1].fragment_id
        yield frag_id1, frag_id2, vertex_i, vertex_j


def chain_to_xmlrpc_list(chain):
    """Converts the Atoms of a Chain/Segment to a list of dictionaries 
    for transfer over xmlrpc.  Only the information required to fit
    TLS groups to the Atoms is set in the Atom dictionaries to
    reduce traffic over the xmlrpc calls.
    """
    xmlrpc_chain = []

    for atm in chain.iter_all_atoms():
        if atm.include is False:
            continue

        atm_desc = {}
        xmlrpc_chain.append(atm_desc)

        atm_desc["name"] = atm.name
        atm_desc["frag_id"] = atm.fragment_id

        frag = atm.get_fragment()
        atm_desc["ifrag"] = chain.index(frag)

        atm_desc["x"] = atm.position[0]
        atm_desc["y"] = atm.position[1]
        atm_desc["z"] = atm.position[2]

        atm_desc["u_iso"] = Constants.B2U * atm.temp_factor

        U = atm.get_U()
        atm_desc["u11"] = U[0,0]
        atm_desc["u22"] = U[1,1]
        atm_desc["u33"] = U[2,2]
        atm_desc["u12"] = U[0,1]
        atm_desc["u13"] = U[0,2]
        atm_desc["u23"] = U[1,2]

        ## calculate weight
        atm_desc["weight"] = calc_atom_weight(atm)

    return xmlrpc_chain


def IsoADPDataSmoother(chain):
    """Experimental data smoothing of temperature factors
    """
    tls_analyzer = tlsmdmodule.TLSModelAnalyzer()
    xlist = chain_to_xmlrpc_list(chain)
    tls_analyzer.set_xmlrpc_chain(xlist)

    num_frags = len(chain)

    smooth_uiso = dict()
    num_smooth = 3
    ifrag_start = num_smooth
    ifrag_end = num_frags - num_smooth - 1

    for ifrag in xrange(ifrag_start, ifrag_end + 1):
        smooth_frag = chain[ifrag]
        frag1 = chain[ifrag - num_smooth]
        frag2 = chain[ifrag + num_smooth]

        tlsdict = tls_analyzer.isotropic_fit_segment(frag1.fragment_id, frag2.fragment_id)
        IT, IL, IS, IOrigin = tls_calcs.isotlsdict2tensors(tlsdict)

        for atm, uiso in TLS.iter_itls_uiso(smooth_frag.iter_all_atoms(), IT, IL, IS, IOrigin):
            smooth_uiso[atm] = uiso
        
        if ifrag == ifrag_start:
            for i in range(ifrag_start):
                smooth_frag = chain[i]
                for atm, uiso in TLS.iter_itls_uiso(smooth_frag.iter_all_atoms(), IT, IL, IS, IOrigin):
                    smooth_uiso[atm] = uiso
        elif ifrag == ifrag_end:
            for i in range(ifrag_end + 1, num_frags):
                smooth_frag = chain[i]
                for atm, uiso in TLS.iter_itls_uiso(smooth_frag.iter_all_atoms(), IT, IL, IS, IOrigin):
                    smooth_uiso[atm] = uiso

    for atm, uiso in smooth_uiso.iteritems():
        atm.temp_factor = Constants.U2B * uiso
        atm.U = numpy.identity(3, float) * uiso
    

class TLSChainProcessor(object):
    """Fits all possible residue subsegments of the given chain object
    with TLS parameters and stores the model parameters in the tlsmdfile
    database file.
    """
    def __init__(self, tlsmdfile, chain, min_subsegment_len):
        self.tlsmdfile = tlsmdfile
        self.chain = chain

        self.min_subsegment_len = min_subsegment_len
        
        self.fit_engine = fit_engine.NewTLSGraphChain(conf.globalconf.tls_model)
        self.fit_engine.set_xmlrpc_chain(chain_to_xmlrpc_list(chain))

        self.num_subsegments = 0
        self.total_num_subsegments = calc_num_subsegments(
            chain.count_fragments(), self.min_subsegment_len)

    def prnt_percent_complete(self, p):
        print "(%10d/%10d) %2d%% Complete" % (self.num_subsegments, self.total_num_subsegments, p)

    def process_chain(self):
        print "PROCESSING CHAIN: ", self.chain.chain_id

        tlsmdfile = self.tlsmdfile

        pcomplete = 0
        pcomplete_old = 0

        chain_id = self.chain.chain_id
        
        for frag_id1, frag_id2, i, j in iter_chain_subsegment_descs(
            self.chain, self.min_subsegment_len):

            if tlsmdfile.grh_get_tls_record(chain_id, frag_id1, frag_id2) == None:
                fit_info = self.fit_engine.lsq_fit_segment(frag_id1, frag_id2)
                fit_info["method"]  = "TLS"
                fit_info["chain_id"] = chain_id
                fit_info["frag_id1"] = frag_id1
                fit_info["frag_id2"] = frag_id2
                tlsmdfile.grh_append_tls_record(fit_info)
            
            self.num_subsegments += 1
            pcomplete = round(100.0 * self.num_subsegments / self.total_num_subsegments)
            if pcomplete != pcomplete_old:
                self.prnt_percent_complete(pcomplete)
                pcomplete_old = pcomplete

        print
        print "NUMBER OF CHAIN SUBSEGMENTS FIT WITH TLS PARAMETERS: (%d/%d)" % (
            self.num_subsegments, self.total_num_subsegments)


class TLSSegment(object):
    """Information on a TLS rigid body segment of a protein chain.
    """
    def __init__(self, **args):
        self.chain_id = args["chain_id"]
        self.frag_id1 = args["frag_id1"]
        self.frag_id2 = args["frag_id2"]
        self.lsq_residual = args["lsq_residual"]
        self.method = args["method"]
        self.num_atoms = args["num_atoms"]

        ## added by TLSChainMinimizer
        self.tls_group = None
        self.tls_info = None
        self.itls_info = None
        self.model_tls_info = None
        self.model = None
        self.segment = None
        self.rmsd_b = None

        ## added by HTML generation code
        self.color = None

        ## added by structure comparison
        self.rmsd_pre_alignment = None
        self.sresult = None
        self.superposition_vscrew = None


class ChainPartition(object):
    """Collection of TLSSegment objects describing one multi-TLS
    group partitioning of a protein chain.
    """
    def __init__(self, chain, ntls_constraint):
        self.chain           = chain
        self.ntls_constraint = ntls_constraint
        self.ntls            = ntls_constraint
        self.tls_list        = []
        self.residual        = 0.0

    def add_tls_segment(self, tls):
        """The argument tls is a dictionary containing a bunch of great information
        about the TLS segment which is part of the optimal partitioning.
        """
        assert isinstance(tls, TLSSegment)
        self.tls_list.append(tls)
        self.residual += tls.lsq_residual

    def normalize_residual(self):
        """Devide the residual by the number of residues in the chain and
        convert from A^2 units to B units.
        """
        nres = len(self.chain)
        self.residual = Constants.U2B * math.sqrt(self.residual / nres)

    def is_valid(self):
        """Return True if the optimization is valid; otherwise, return False.
        """
        return len(self.tls_list) > 0

    def num_tls_segments(self):
        return len(self.tls_list)

    def iter_tls_segments(self):
        return iter(self.tls_list)

    def enumerate_tls_segments(self):
        i = 0
        for tls in self.iter_tls_segments():
            yield i, tls
            i += 1


class ChainPartitionCollection(object):
    """Contains all the ChainPartition objects for a chain.
    """
    def __init__(self, struct, struct_file_path, chain, chain_optimizer):
        self.struct = struct
        self.struct_file_path = struct_file_path
        self.chain = chain
        self.chain_optimizer = chain_optimizer
        
        self.chain_id = chain.chain_id
        self.max_ntls = conf.globalconf.nparts

        self.ntls_chain_partition_list = []
        
        for ntls in xrange(1, self.max_ntls + 1):
            cpartition = self.chain_optimizer.calc_chain_partition(ntls)
            if cpartition == None:
                continue
            self.ntls_chain_partition_list.append((ntls, cpartition))

    def num_chain_partitions(self):
        return len(self.ntls_chain_partition_list)

    def iter_ntls_chain_partitions(self):
        return iter(self.ntls_chain_partition_list)

    def iter_chain_partitions(self):
        for ntls, cpartition in self.iter_ntls_chain_partitions():
            yield cpartition

    def iter_ntls(self):
        for ntls, cpartition in self.iter_ntls_chain_partitions():
            yield ntls

    def get_chain_partition(self, find_ntls):
        for ntls, cpartition in self.ntls_chain_partition_list:
            if ntls == find_ntls:
                return cpartition
        return None


class TLSChainMinimizer(hcsssp.HCSSSP):
    """Finds the minimal TLS description of a given Chain object using
    the HCSSSP global optimization algorithm and a constraint on the number
    of TLS which can be used for the minimization.
    """
    def __init__(self, tlsmdfile, chain, min_subsegment_len, nparts):
        hcsssp.HCSSSP.__init__(self)

        self.tlsmdfile = tlsmdfile
        self.min_subsegment_len = min_subsegment_len
        self.nparts = nparts

        self.tls_cache = {}

        self.minimized = False
        self.D         = None
        self.P         = None
        self.T         = None

        ## source and destination fragment IDs; these are the
        ## begining and ending fragments(residues) of the chain
        self.frag_id_src  = None
        self.frag_id_dest = None

        if self.frag_id_src == None and self.frag_id_dest == None:
            self.chain = chain
        else:
            self.chain = chain[self.frag_id_src:self.frag_id_dest]

        ## tls analyzer is necessary for re-fitting the partitions
        ## chosen by the optimization (minimization)
        self.tls_analyzer = tlsmdmodule.TLSModelAnalyzer()
        xlist = chain_to_xmlrpc_list(self.chain)
        self.tls_analyzer.set_xmlrpc_chain(xlist)

        ## this is useful: for each fragment in the minimization
        ## set a attribute for its index position
        ichain = 0
        for frag in self.chain.iter_fragments():
            frag.ichain = ichain
            ichain = ichain + 1
        
        ## Initialize the vertex list based on the number of of
        ## residues in Chain.
        self.num_vertex = len(self.chain) + 1

    def run_minimization(self):
        """Run the HCSSSP minimization on the self.V,self.E graph, resulting
        in the creation of the self.D, self.P, and self.T arrays which
        contain 
        """
        ## build the vertex labels to reflect the protein structure
        ## the graph spans
        vertices = []
        for i in xrange(self.num_vertex):

            ## add the vertex label for i at Vi
            if i ==0 :
                vertex_label = "N-TERM"

            elif i == self.num_vertex - 1:
                vertex_label = "C-TERM"

            else:
                vertex_label = "%s{%s:%s}" % (self.chain.chain_id,
                                              self.chain[i-1].fragment_id,
                                              self.chain[i].fragment_id)

            vertex_label = "V%d[%s]" % (i, vertex_label)
            vertices.append(vertex_label)

        ## now build edges for the graph with weights given by the LSQ
        ## residual of TLS group fits
        grh_get_tls_record = self.tlsmdfile.grh_get_tls_record
        
        edges = []
        for frag_id1, frag_id2, i, j in iter_chain_subsegment_descs(
            self.chain, self.min_subsegment_len):

            tlsdict = grh_get_tls_record(self.chain.chain_id, frag_id1, frag_id2)

            assert frag_id1 == tlsdict["frag_id1"]
            assert frag_id2 == tlsdict["frag_id2"]
            
            if tlsdict == None:
                print "[ERROR] no TLS group %s{%s..%s}" % (self.chain.chain_id, frag_id1, frag_id2)
                raise SystemExit

            if tlsdict.has_key("error") == True:
                continue

            if tlsdict.has_key("lsq_residual") == False:
                print "[ERROR] no lsq_residual! %s{%s..%s}" % (self.chain.chain_id, frag_id1, frag_id2)
                raise SystemExit

            cost = tlsdict["lsq_residual"]
            frag_range = (frag_id1, frag_id2)
            edge = (i, j, cost, frag_range)
            edges.append(edge)

        ## perform the minimization
        misc.start_timing()

        if len(edges) > 0:
            print "run_minimization(chain_id=%s): HCSSSP Minimizing..." % (self.chain.chain_id)
        
            D, P, T = self.HCSSSP_minimize(vertices, edges, self.nparts)

            self.minimized = True
            self.V = vertices
            self.D = D
            self.P = P
            self.T = T
        else:
            print "run_minimization(chain_id=%s): Unable to minimize" % (self.chain.chain_id)
            self.minimized = False

        ## free memory taken up from edges
        edges = None
        gc.collect()

        print "run_minimization(): ", misc.end_timing()

    def calc_chain_partition(self, nparts):
        """Return a ChainPartition() object containing the optimal
        TLS description of self.chain using num_tls_segments.
        """
        if not self.minimized:
            return None

        print "Re-Fitting TLS Parameters of Optimized Chain %s Partitioned using %d TLS Groups" % (
            self.chain.chain_id, nparts)

        cpartition = ChainPartition(self.chain, nparts)
        
        partition_num = 0
        
        for hi, hj, edge in self.HCSSSP_path_iter(self.V, self.D, self.P, self.T, nparts):
            if edge == None:
                continue
            
            i, j, cost, frag_range = edge
            
            ## check if the edge is a bypass-edge type
            if len(frag_range) == 2:
                print "    Fitting Chain Segment %s-%s" % (frag_range[0], frag_range[1])
                tls = self.__calc_tls_record_from_edge(edge)
                tls.partition_num = partition_num
                partition_num += 1
                cpartition.add_tls_segment(tls)

        cpartition.normalize_residual()
        return cpartition

    def __calc_nonlinear_fit(self, tls):
        """Use the non-linear TLS model to calculate tensor values.
        """
        tls_group = tls.tls_group

        ## anisotropic model
        tlsdict = self.tls_analyzer.constrained_anisotropic_fit_segment(tls.frag_id1, tls.frag_id2)
        T, L, S, origin = tls_calcs.tlsdict2tensors(tlsdict)
        tls_group.T = T
        tls_group.L = L
        tls_group.S = S
        tls_group.origin = origin

        ## isotropic model
        itlsdict = self.tls_analyzer.constrained_isotropic_fit_segment(tls.frag_id1, tls.frag_id2)
        IT, IL, IS, IOrigin = tls_calcs.isotlsdict2tensors(itlsdict)
        tls_group.itls_T = IT
        tls_group.itls_L = IL
        tls_group.itls_S = IS

        assert numpy.allclose(tls_group.origin, IOrigin)

    def __calc_tls_record_from_edge(self, edge):
        """Independently calculate the TLS parameters for the segment
        to verify correctness (internal check).
        """
        i, j, cost, frag_range = edge

        frag_id1, frag_id2 = frag_range

        tlsdict = self.tlsmdfile.grh_get_tls_record(self.chain.chain_id, frag_id1, frag_id2)
        tls = TLSSegment(**tlsdict)
        tls.segment = self.chain[frag_id1:frag_id2]

        tls.tls_group = TLS.TLSGroup()
        for atm in tls.segment.iter_all_atoms():
            if atm.include == False:
                continue
            tls.tls_group.append(atm)

        ## use the constrained TLS model to calculate tensor values
        cache_key = (frag_id1, frag_id2)
        if self.tls_cache.has_key(cache_key):
            tlscache = self.tls_cache[cache_key]
            tls_group_cache = tlscache.tls_group
            tls.tls_group.origin = tls_group_cache.origin.copy()
            tls.tls_group.T = tls_group_cache.T.copy()
            tls.tls_group.L = tls_group_cache.L.copy()
            tls.tls_group.S = tls_group_cache.S.copy()
            tls.tls_group.itls_T = tls_group_cache.itls_T
            tls.tls_group.itls_L = tls_group_cache.itls_L.copy()
            tls.tls_group.itls_S = tls_group_cache.itls_S.copy()
        else:
            self.__calc_nonlinear_fit(tls)
            self.tls_cache[cache_key] = tls
        
        ## helpful additions
        tls_info  = tls.tls_group.calc_tls_info()
        itls_info = TLS.calc_itls_center_of_reaction(
            tls.tls_group.itls_T,
            tls.tls_group.itls_L,
            tls.tls_group.itls_S,
            tls.tls_group.origin)

        tls.tls_info = tls_info
        tls.itls_info = itls_info

        if conf.globalconf.tls_model in ["ISOT", "NLISOT"]:
            tls.tls_group.model = "ISOT"
            tls.model_tls_info = itls_info
        elif conf.globalconf.tls_model in ["ANISO", "NLANISO"]:
            tls.tls_group.model = "ANISO"
            tls.model_tls_info = tls_info

        tls.rmsd_b = tls_calcs.calc_rmsd_tls_biso(tls.tls_group)

        return tls

    def prnt_detailed_paths(self):
        """Debug
        """
        hops = self.nparts
        
        if not self.minimized:
            return
        dest_j = len(self.V)-1

        for h in xrange(1, hops+1):
            print
            print "MINIMIZATON VERTEX PATH FOR %d SEGMENTS" % (h)
            print "NODE LABEL              HOPS      COST      PREVIOUS NODE          EDGE"
            self.__detailed_path(self.V, self.D, self.P, self.T, h)

    def __detailed_path(self, V, D, P, T, hop_constraint):
        """Print out the path from the source vertex (vertex 0) to
        the destination vertex (end vertex) given the hop_constraint.
        """
        num_vertex = len(D[0])
        
        ## start at the destination vertex
        curr_v = num_vertex - 1
        h      = hop_constraint
        
        while curr_v >= 0:
            prev_vertex  = P[h,curr_v]
            vertex_label = V[curr_v].ljust(20)

            if prev_vertex < 0:
                prev_vertex_label = "".ljust(20)
            else:
                prev_vertex_label = V[prev_vertex].ljust(20)

            edge = T[h][curr_v]

            if edge is not None:
                i, j, cost, frag_range = edge
                wr = cost / (j - i)
                edge_label = "(%3d,%3d,%6.3f,%s) %6.3f" % (i, j, cost, frag_range, wr)
            else:
                edge_label = ""
                 
            print "%s   %3d     %10.4f   %s   %s" % (
                vertex_label, h, D[h,curr_v], prev_vertex_label, edge_label)

            curr_v = prev_vertex
            h -= 1


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

        if sel_chain_ids!=None:
            self.sel_chain_ids = sel_chain_ids.split(",")
        else:
            self.sel_chain_ids = None

        self.tlsdb_file      = tlsdb_file
        self.tlsdb_complete  = tlsdb_complete

        self.struct          = None
        self.struct_id       = None
        self.chains          = None
        self.tlsmdfile       = None

    def run_optimization(self):
        """Run the TLSMD optimization on the structure.
        """
        self.load_struct()

        ## auto name of tlsdb file then open
        if self.tlsdb_file==None:
            self.tlsdb_file = "%s_%s_%s.db" % (
                self.struct_id, conf.globalconf.tls_model,
                conf.globalconf.weight_model)

        self.tlsmdfile = datafile.TLSMDFile(self.tlsdb_file)

        ## select chains for analysis
        self.select_chains()

        ## print these settings
        self.prnt_settings()

        self.set_atom_include_flags()

        if not self.tlsdb_complete:
            self.construct_tls_segment_database()
        
        self.calc_chain_minimization()

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
            if self.sel_chain_ids!=None:
                if chain.chain_id not in self.sel_chain_ids:
                    continue

            ## count the number of amino acid residues in the chain
            if chain.count_amino_acids()<10:
                continue

            ## ok, use the chain but use a segment and cut off
            ## any leading and trailing non-amino acid residues
            frag_id1 = None
            for aa in chain.iter_amino_acids():
                if frag_id1==None:
                    frag_id1 = aa.fragment_id
                frag_id2 = aa.fragment_id
                
            segment = chain[frag_id1:frag_id2]
            segments.append(segment)
            segment.struct = self.struct
        
        self.chains = segments
        
    def set_atom_include_flags(self):
        """Sets that atm.include attribute for each atom in the chains
	being analyzed by tlsmd.
	"""
	for chain in self.chains:
	    for atm in chain.iter_all_atoms():
                atm.include = calc_include_atom(atm)
            IsoADPDataSmoother(chain)
            
	
    def construct_tls_segment_database(self):
        """Calculates the TLSGraph for each chain in self.chains, optionally
        loading pre-computed graphs from the graph_file.  Any chains
        which need TLSGraphs computed will be strored in the graph_file.
        """
        for chain in self.chains:
            misc.begin_chain_timing(chain.chain_id)

            chain_processor = TLSChainProcessor(
                self.tlsmdfile, chain, conf.globalconf.min_subsegment_size)

            print "BUILDING TLS SEGMENT DATABASE FOR %s" % (chain)
            chain_processor.process_chain()

            misc.end_chain_timing(chain.chain_id)

    def calc_chain_minimization(self):
        """Performs the TLS graph minimization on all TLSGraphs.
        """
        for chain in self.chains:
            tls_chain_minimizer = TLSChainMinimizer(
                self.tlsmdfile,
                chain,
                conf.globalconf.min_subsegment_size,
                conf.globalconf.nparts)
            
            tls_chain_minimizer.run_minimization()
            if not tls_chain_minimizer.minimized:
                continue

            print
            print "="*79
            print "MINIMIZING CHAIN %s" % (chain)
            tls_chain_minimizer.prnt_detailed_paths()

            chain.partition_collection = ChainPartitionCollection(
                self.struct,
                self.struct_file_path,
                chain,
                tls_chain_minimizer)

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
