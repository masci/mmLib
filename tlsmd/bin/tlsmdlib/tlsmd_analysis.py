## TLS Motion Determination (TLSMD)
## Copyright 2002 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import sys
import math
import numpy

from mmLib import Constants, FileLoader
from mmLib.Extensions import TLS

import misc
import const
import conf
import datafile
import hcsssp
import tls_calcs
import fit_engine


def calc_include_atom(atm, reject_messages = False):
    """Filter out atoms from the model which will cause problems or
    cont contribute to the TLS analysis.
    """
    if atm.occupancy<0.1:
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
    for i in range(num_vertex):
        for j in range(i + min_len, num_vertex):
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
        if atm.include == False:
            continue

        atm_desc = {}
        xmlrpc_chain.append(atm_desc)

        atm_desc["name"]    = atm.name
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


class TLSChainProcessor(object):
    """Uses a TLSGridServerPool object to fit all possible TLS subsegment
    to a given Chain object.  The resulting TLS groups are stored in the
    analysis.datafile.grh_*() file.
    """
    
    def __init__(self, analysis, chain, min_subsegment_len):
        self.analysis = analysis
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
        """Fits TLS groups to all possible subsegments of the given Chain,
        storing the results in the analysis.grh_*() file.  
        """
        print "PROCESSING CHAIN: ", self.chain.chain_id

        tlsmdfile = self.analysis.tlsmdfile

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
            if pcomplete!=pcomplete_old:
                self.prnt_percent_complete(pcomplete)
                pcomplete_old = pcomplete

        print
        print "NUMBER OF CHAIN SUBSEGMENTS FIT WITH TLS PARAMETERS: (%d/%d)" % (
            self.num_subsegments, self.total_num_subsegments)


class TLSOptimization(object):
    """Collection object containing one multi-TLS group partition of a protein chain.
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
        self.tls_list.append(tls)
        self.residual += tls["lsq_residual"]

    def normalize_residual(self):
        """Devide the residual by the number of residues in the chain and
        convert from A^2 units to B units.
        """
        nres = len(self.chain)
        self.residual = Constants.U2B * math.sqrt(self.residual / nres)

    def is_valid(self):
        """Return True if the optimization is valid; otherwise, return False.
        """
        return len(self.tls_list)>0


class TLSChainMinimizer(hcsssp.HCSSSP):
    """Finds the minimal TLS description of a given Chain object using
    the HCSSSP global optimization algorithm and a constraint on the number
    of TLS which can be used for the minimization.
    """
    def __init__(self, analysis, chain, min_subsegment_len, nparts):
        hcsssp.HCSSSP.__init__(self)

        self.analysis = analysis
        self.min_subsegment_len = min_subsegment_len
        self.nparts = nparts

        self.minimized = False
        self.D         = None
        self.P         = None
        self.T         = None

        ## source and destination fragment IDs; these are the
        ## begining and ending fragments(residues) of the chain
        self.frag_id_src  = None
        self.frag_id_dest = None

        if self.frag_id_src==None and self.frag_id_dest==None:
            self.chain = chain
        else:
            self.chain = chain[self.frag_id_src:self.frag_id_dest]

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
        for i in range(self.num_vertex):

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
        grh_get_tls_record = self.analysis.tlsmdfile.grh_get_tls_record
        
        edges = []
        for frag_id1, frag_id2, i, j in iter_chain_subsegment_descs(
            self.chain, self.min_subsegment_len):

            tls = grh_get_tls_record(self.chain.chain_id, frag_id1, frag_id2)

            assert frag_id1 == tls["frag_id1"]
            assert frag_id2 == tls["frag_id2"]
            
            if tls==None:
                print "[ERROR] no TLS group %s{%s..%s}" % (self.chain.chain_id, frag_id1, frag_id2)
                sys.exit(-1)

            if tls.has_key("error") == True:
                continue

            if tls.has_key("lsq_residual") == False:
                print "[ERROR] no lsq_residual! %s{%s..%s}" % (self.chain.chain_id, frag_id1, frag_id2)
                sys.exit(-1)

            cost = tls["lsq_residual"]
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
        import gc
        gc.collect()

        print "run_minimization(): ", misc.end_timing()

    def calc_tls_optimization(self, nparts):
        """Return a TLSOptimization() object containing the optimal
        TLS description of self.chain using num_tls_segments.
        """
        if not self.minimized:
            return None

        print "Re-Fitting TLS Parameters of Optimized Chain %s Partitioned using %d TLS Groups" % (
            self.chain.chain_id, nparts)

        tlsopt = TLSOptimization(self.chain, nparts)
        
        for hi, hj, edge in self.HCSSSP_path_iter(self.V, self.D, self.P, self.T, nparts):
            if edge==None: continue
            
            i, j, cost, frag_range = edge
            
            ## check if the edge is a bypass-edge type
            if len(frag_range)==2:
                print "    Fitting Chain Segment %s-%s" % (frag_range[0], frag_range[1])
                tls = self.__calc_tls_record_from_edge(edge)
                tlsopt.add_tls_segment(tls)

        tlsopt.normalize_residual()

        return tlsopt

    def __calc_nonlinear_fit(self, tls_group, segment):
        """Use the non-linear TLS model to calculate tensor values.
        """
        import nonlineartls
        nltls = nonlineartls.NLTLSModel()

        xlist = chain_to_xmlrpc_list(segment)
        nltls.set_xmlrpc_chain(xlist)

        ## anisotropic model
        tls = nltls.anisotropic_fit_segment(0, len(xlist)-1)

        tls_group.origin = numpy.array([tls["x"], tls["y"], tls["z"]], float)

        tls_group.T = numpy.array(
            [ [tls["t11"], tls["t12"], tls["t13"]],
              [tls["t12"], tls["t22"], tls["t23"]],
              [tls["t13"], tls["t23"], tls["t33"]] ], float)
        
        tls_group.L = numpy.array(
            [ [tls["l11"], tls["l12"], tls["l13"]],
              [tls["l12"], tls["l22"], tls["l23"]],
              [tls["l13"], tls["l23"], tls["l33"]] ], float)
        
        s11, s22, s33 = TLS.calc_s11_s22_s33(tls["s2211"], tls["s1133"]) 
        
        tls_group.S = numpy.array(
            [ [       s11, tls["s12"], tls["s13"]],
              [tls["s21"],        s22, tls["s23"]],
              [tls["s31"], tls["s32"],       s33] ], float)

        ## isotropic model
        itls = nltls.isotropic_fit_segment(0, len(xlist)-1)
        
        tls_group.itls_T = itls["it"]
        
        tls_group.itls_L = numpy.array(
            [ [itls["il11"], itls["il12"], itls["il13"]],
              [itls["il12"], itls["il22"], itls["il23"]],
              [itls["il13"], itls["il23"], itls["il33"]] ], float)

        tls_group.itls_S = numpy.array([itls["is1"], itls["is2"], itls["is3"]], float)

    def __calc_tls_record_from_edge(self, edge):
        """Independently calculate the TLS parameters for the segment
        to verify correctness (internal check).
        """
        i, j, cost, frag_range = edge

        ## retrieve from cache if possible
        frag_id1, frag_id2 = frag_range

        tls = self.analysis.tlsmdfile.grh_get_tls_record(self.chain.chain_id, frag_id1, frag_id2)
        segment = self.chain[frag_id1:frag_id2]

        ## create TLSGroup
        tls_group = TLS.TLSGroup()

        ## add atoms to the group
        for atm in segment.iter_all_atoms():
            if atm.include == False: continue
            tls_group.append(atm)

        ## use the constrained TLS model to calculate tensor values
        self.__calc_nonlinear_fit(tls_group, segment)
        
        ## helpful additions
        tls_info  = tls_group.calc_tls_info()
        itls_info = TLS.calc_itls_center_of_reaction(
            tls_group.itls_T, tls_group.itls_L, tls_group.itls_S, tls_group.origin)

        tls["tls_info"]  = tls_info
        tls["itls_info"] = itls_info

        if conf.globalconf.tls_model in ["ISOT", "NLISOT"]:
            tls_group.model = "ISOT"
            tls["model_tls_info"] = itls_info
        elif conf.globalconf.tls_model in ["ANISO", "NLANISO"]:
            tls_group.model = "ANISO"
            tls["model_tls_info"] = tls_info

        tls["tls_group"]            = tls_group
        tls["segment"]              = segment
        tls["rmsd_b"]               = tls_calcs.calc_rmsd_tls_biso(tls_group)

        return tls

    def prnt_detailed_paths(self):
        """Debug
        """
        hops = self.nparts
        
        if not self.minimized:
            return
        dest_j = len(self.V)-1

        for h in range(1, hops+1):
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
        
        while curr_v>=0:
            prev_vertex  = P[h,curr_v]
            vertex_label = V[curr_v].ljust(20)

            if prev_vertex<0:
                prev_vertex_label = "".ljust(20)
            else:
                prev_vertex_label = V[prev_vertex].ljust(20)

            edge = T[h][curr_v]

            if edge!=None:
                i, j, cost, frag_range = edge
                wr = cost / (j - i)
                edge_label = "(%3d,%3d,%6.3f,%s) %6.3f" % (i, j, cost, frag_range, wr)
            else:
                edge_label = ""
                 
            print "%s   %3d     %10.4f   %s   %s" % (
                vertex_label, h, D[h,curr_v], prev_vertex_label, edge_label)

            curr_v = prev_vertex
            h -= 1

    def hinge_plot(self):
	"""Spiffy new hinge-prediction algorithm.
	"""
        import lineartls
        import nonlineartls

        fil = open("hinge_chain_%s.txt" % (self.chain.chain_id), "w")

        xchain = fit_engine.XChain(chain_to_xmlrpc_list(self.chain))

        tls_model = lineartls.LinearTLSModel()
        #tls_model = nonlineartls.NLTLSModel()

        tls_model.set_xmlrpc_chain(xchain.xmlrpc_chain)

        grh_get_tls_record = self.analysis.tlsmdfile.grh_get_tls_record
        chain_id = self.chain.chain_id

        win = 12

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
                chain_id, frag_id1a, frag_id2a, frag_id1b, frag_id2b)
            
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


class TLSMDAnalysis(object):
    """Central object for a whole-structure TLS analysis.
    """
    def __init__(self,
                 struct_path    = None,
                 sel_chain_ids  = None,
                 tlsdb_file     = None,
                 tlsdb_complete = False):

        conf.globalconf.prnt()

        self.struct_path     = struct_path
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
            self.tlsdb_file = "%s_%s_%s.db" % (self.struct_id, conf.globalconf.tls_model,
                                               conf.globalconf.weight_model)
        self.tlsmdfile = datafile.TLSMDFile(self.tlsdb_file)

        ## select chains for analysis
        self.select_chains()

        ## print these settings
        self.prnt_settings()

        self.set_atom_include_flags()

        if not self.tlsdb_complete:
            self.calc_tls_segments()
        
        self.calc_chain_minimization()

    def prnt_settings(self):
        chain_ids = []
        for chain in self.chains:
            chain_ids.append(chain.chain_id)
        cids = ",".join(chain_ids)
        
        print "TLSMD ANALYSIS SETTINGS"
        print "    STRUCTURE FILE.....................: %s" % (self.struct_path)
        print "    STRUCTURE ID.......................: %s" % (self.struct_id)
        print "    CHAIN IDs SELECTED FOR ANALYSIS....: %s" % (cids)
        print "    DATABASE FILE PATH.................: %s" % (self.tlsdb_file)
        print
        
    def load_struct(self):
        """Loads Structure, chooses a unique struct_id string.
        """
        print "LOADING STRUCTURE"
        print "    PATH: %s" % (self.struct_path)

        ## load struct
        self.struct = FileLoader.LoadStructure(
            fil = self.struct_path,
            build_properties = ("library_bonds","distance_bonds"))

        ## set the structure ID
        if conf.globalconf.struct_id != None:
            struct_id = conf.globalconf.struct_id
        else:
            struct_id = self.struct.structure_id
            conf.globalconf.struct_id = struct_id

        self.struct.structure_id = struct_id
        self.struct_id = struct_id

        print "    STRUCT ID: %s" % (self.struct_id)
        print

        ## if there are REFMAC5 TLS groups in the REMARK records of
        ## the PDB file, then add those in
        tls_file = TLS.TLSFile()
        tls_file.set_file_format(TLS.TLSFileFormatPDB())

        fil = open(self.struct_path, "r")
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
        
        self.chains = segments
        
    def set_atom_include_flags(self):
        """Sets that atm.include attribute for each atom in the chains
	being analyzed by tlsmd.
	"""
	for chain in self.chains:
	    for atm in chain.iter_all_atoms():
                atm.include = calc_include_atom(atm)
	
    def calc_tls_segments(self):
        """Calculates the TLSGraph for each chain in self.chains, optionally
        loading pre-computed graphs from the graph_file.  Any chains
        which need TLSGraphs computed will be strored in the graph_file.
        """
        for chain in self.chains:
            misc.begin_chain_timing(chain.chain_id)

            chain_processor = TLSChainProcessor(self, chain, conf.globalconf.min_subsegment_size)

            print "BUILDING TLS SEGMENT DATABASE FOR %s" % (chain)
            chain_processor.process_chain()

            misc.end_chain_timing(chain.chain_id)

    def calc_chain_minimization(self):
        """Performs the TLS graph minimization on all TLSGraphs.
        """
        for chain in self.chains:
            chain.tls_chain_minimizer = TLSChainMinimizer(
                self, chain, conf.globalconf.min_subsegment_size, conf.globalconf.nparts)
            
            chain.tls_chain_minimizer.run_minimization()
            if not chain.tls_chain_minimizer.minimized: continue

            print
            print "="*79
            print "MINIMIZING CHAIN %s" % (chain)
            chain.tls_chain_minimizer.prnt_detailed_paths()
