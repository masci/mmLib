## TLS Minimized Domains (TLSMD)
## Copyright 2002 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

###############################################################################
## Grid/TLSGraph Calculation System
#

import sys
import copy
import string

from misc       import *
from fit_engine import *
from datafile   import TLSMDFile
from hcsssp     import HCSSSP


###############################################################################
## Hacked Global Options
##

## Minimum number of atoms for different models
MIN_ISOTROPIC_ATOMS   = 40
MIN_ANISOTROPIC_ATOMS = 10

## ANISO/ISO TLS Model: set externally!
TLS_MODEL = "HYBRID"

## LSQ weighting
WEIGHT_MODEL = "UNIT"

## atom selection options
INCLUDE_ATOMS = "ALL" 

## minimum span of residues for TLS subsegments
MIN_SUBSEGMENT_SIZE = 4

def print_globals():
    print "TLSMD GLOBAL OPTIONS"
    print "    TLS_MODEL ==================: %s" % (TLS_MODEL)
    print "    MIN_SUBSEGMENT_SIZE --------: %d" % (MIN_SUBSEGMENT_SIZE)
    print "    WEIGHT_MODEL ===============: %s" % (WEIGHT_MODEL)
    print "    INCLUDE ATOMS --------------: %s" % (INCLUDE_ATOMS)
    print

def set_globals():
    """Takes settings from misc.GLOBALS dictionary and applies them to
    the globals in the current module.
    """
    global TLS_MODEL
    global WEIGHT_MODEL
    global INCLUDE_ATOMS
    global MIN_SUBSEGMENT_SIZE

    TLS_MODEL     = GLOBALS["TLS_MODEL"]
    WEIGHT_MODEL  = GLOBALS["WEIGHT_MODEL"]
    INCLUDE_ATOMS = GLOBALS["INCLUDE_ATOMS"]

    assert TLS_MODEL     in ["ANISO", "HYBRID", "NLANISO", "NLISOT"]
    assert WEIGHT_MODEL  in ["UNIT", "IUISO"]
    assert INCLUDE_ATOMS in ["ALL", "MAINCHAIN", "CA"]

    if INCLUDE_ATOMS=="ALL":
        MIN_SUBSEGMENT_SIZE = 8
    elif INCLUDE_ATOMS=="MAINCHAIN":
        MIN_SUBSEGMENT_SIZE = 8
    elif INCLUDE_ATOMS=="CA":
        MIN_SUBSEGMENT_SIZE = 20

    GLOBALS["MIN_SUBSEGMENT_SIZE"] = MIN_SUBSEGMENT_SIZE

###############################################################################
## Atom Selection/Weighting Functions
##

def calc_include_atom(atm, reject_messages=False):
    """Filter out atoms from the model which will cause problems or
    cont contribute to the TLS analysis.
    """
    if atm.occupancy<0.1:
        if reject_messages==True:
            print "calc_include_atom(%s): rejected because of low occupancy" % (atm)
	return False
    
    if trace(atm.get_U())<=TSMALL:
        if reject_messages==True:
            print "calc_include_atom(%s): rejected because of small Uiso magnitude " % (atm)
        return False

    if INCLUDE_ATOMS=="ALL":
        return True

    elif INCLUDE_ATOMS=="MAINCHAIN":
        if atm.name not in ["N", "CA", "C", "O"]:
            if reject_messages==True:
                print "calc_include_atom(%s): rejected non-mainchain atom" % (atm)
            return False

    elif INCLUDE_ATOMS=="CA":
        if atm.name!="CA":
            if reject_messages==True:
                print "calc_include_atom(%s): rejected non-CA atom" % (atm)
            return False
        
    return True

def calc_atom_weight(atm):
    """Weight the least-squares fit according to this function.
    """
    weight = atm.occupancy

    if WEIGHT_MODEL=="IUISO":
        weight = weight * (1.0 / (B2U * atm.temp_factor))

    return weight

def iter_ij(num_vertex, min_span):
    """Iterates over the i,j vertex indexes defining the edges
    to be built for the graph.  num_vertex gives the number
    of consecutive vertixes to create, min_span is the minimum
    number of residues (fragments) a edge should span.
    """
    for i in range(num_vertex):
        for j in range(i+min_span, num_vertex):
            yield i, j

def iter_chain_subsegment_descs(chain, min_span):
    """Iterate over all possible subsegments of the given Chain object
    with a minimum size of min_span fragments.  The segments are yielded
    as Python dictionaries containing a description of the subsegment.
    """
    num_vertex = len(chain) + 1
    
    for vertex_i, vertex_j in iter_ij(num_vertex, min_span):
        ## begin message
        msg                 = {}
        msg["chain"]        = chain
        msg["chain_id"]     = chain.chain_id

        ## vertex span of edge info
        msg["vertex_i"]     = vertex_i
        msg["vertex_j"]     = vertex_j

        ## the definition of the subsegment
        msg["frag_id1"]     = chain[vertex_i].fragment_id
        msg["frag_id2"]     = chain[vertex_j-1].fragment_id

        yield msg


###############################################################################
## XMLRPC Client Threads which dispach jobs to the chain graphers
##

class TLSGridClientThread(Thread):
    """One instance of this class connects to a single xmlrpc server
    on the compute grid.controlls a single compute
    server on the grid.  It uses blocking xmlrpc calls within its own
    thread.
    """
    def __init__(self, server_url, request_queue, result_queue):
        Thread.__init__(self)
        self.setDaemon(True)

        ## set this to True to monitor xmlrpc messages
        self.debug         = False
        
        self.server_url    = server_url
        self.request_queue = request_queue
        self.result_queue  = result_queue

    def get_ident(self):
        if self.server_url==None:
            return "THREAD %d" % (thread.get_ident())
        return "%s THREAD %d" % (self.server_url, thread.get_ident())

    def run(self):
        """The main() for this thread. 
        """
        ## use either a in-process server, or remote
        ## server over xmlrpc
        if self.server_url!=None:
            self.compute_server = self.xmlrpc_connect()
        else:
            self.compute_server = NewTLSGraphChain(TLS_MODEL)

        self.run_grid_client()

    def xmlrpc_connect(self):
        """Connect to the xmlrpc server and return a connection object.
        """
        server = xmlrpclib.ServerProxy(self.server_url)
        server.set_tls_model(TLS_MODEL, WEIGHT_MODEL)
        return server

    def run_grid_client(self):
        """Enters a loop to marshal the request/response messages
        between this server client thread and the controller thread.
        """
        cur_chain_id = None
        
        while True:
            ## get a message requesting the LSQ fit of a TLS
            ## Chain segment (Python dictionary)
            msg = self.request_queue.get()

            ## check the xmlrpc_chain_id of the request to see if
            ## the TLSGridServer needs to be sent the working data
            ## for the chain (if it is a new chain)
            if msg["chain_id"]!=cur_chain_id:
                self.call_set_xmlrpc_chain(msg)
                cur_chain_id = msg["chain_id"]

            ## perform the LSQ fit
            msg["server_url"] = self.server_url
            return_msg = self.call_lsq_fit_segment(msg)

            ## place result on the out_queue
            self.result_queue.put(return_msg)

    def call_set_xmlrpc_chain(self, msg):
        """Calls the xmlrpc server's set_xmlrpc_chain() method.
        """
        if self.debug==True:
            print "TLSGridServerThread.call_set_xmlrpc_chain(%s)" % (
                self.get_ident())

        xmlrpc_chain = msg["xmlrpc_chain"]
        self.compute_server.set_xmlrpc_chain(xmlrpc_chain)

    def call_lsq_fit_segment(self, msg):
        """Calls the xmlrpc server's lsq_fit_segment() method.
        """
        if self.debug==True:
            print "TLSGridServerThread.call_lsq_fit_segment(%s)" % (
                self.get_ident())
        
	frag_id1 = msg["frag_id1"]
	frag_id2 = msg["frag_id2"]

        fit_info = self.compute_server.lsq_fit_segment(
            frag_id1, frag_id2)

        ## merge the fit_info data with the original message
        for key in fit_info.keys():
            if msg.has_key(key):
                assert msg[key]==fit_info[key]
            msg[key] = fit_info[key]

        return msg


class TLSGridServerPool(object):
    """Launches a pool of threads using the TLSGridClientThread() class as
    the controlling object for each of these threads.  Each thread in the
    pool waits for jobs placed in this object's self.requeust_queue, and
    the results are placed by these threads back into self.result_queue.
    """
    def __init__(self):
        self.client_thread_list = []
        self.request_queue      = Queue(1)
        self.result_queue       = Queue(1)

    def launch_client_thread(self, server_url):
        """Launches one TLSGridServerThread.
        """
        client_thread = TLSGridClientThread(
            server_url, self.request_queue, self.result_queue)
        
        self.client_thread_list.append(client_thread)
        client_thread.start()

    def iter_segment_processor(self, iter_segment_producer):
        """Consume the seg_info dictionaries produced by the seg_producer
        iterator, and run them concurrently on all avilible
        TLSGridServerThread threads.
        Yield back the calculation results in the form of a fit_info
        dict for each calcuation request.
        """
        open_requests = 0

        ## make all requests
        for seg_info in iter_segment_producer:
            self.request_queue.put(seg_info)
            open_requests += 1

            while not self.result_queue.empty():
                yield self.result_queue.get()
                open_requests -= 1
        
        ## wait for any remaining requests
        while open_requests>0:
            yield self.result_queue.get()
            open_requests -= 1


class TLSChainProcessor(object):
    """Uses a TLSGridServerPool object to fit all possible TLS subsegment
    to a given Chain object.  The resulting TLS groups are stored in the
    analysis.datafile.grh_*() file.
    """
    def __init__(self, server_pool):
        self.server_pool = server_pool

    def convert_chain_for_xmlrpc(self, chain):
        """Converts the Atoms of a Chain/Segment to a list of dictionaries 
        for transfer over xmlrpc.  Only the information required to fit
        TLS groups to the Atoms is set in the Atom dictionaries to
        reduce traffic over the xmlrpc calls.
        """
        print "convert_chain_for_xmlrpc(chain_id=%s)" % (chain.chain_id)
        
        xmlrpc_chain = []

        for atm in chain.iter_all_atoms():

            if calc_include_atom(atm, True)==False:
                continue

            atm_desc = {}
            xmlrpc_chain.append(atm_desc)

            atm_desc["name"]    = atm.name
            atm_desc["frag_id"] = atm.fragment_id

            atm_desc["x"] = atm.position[0]
            atm_desc["y"] = atm.position[1]
            atm_desc["z"] = atm.position[2]

            atm_desc["u_iso"] = B2U * atm.temp_factor

            U = atm.get_U()
            atm_desc["u11"] = U[0,0]
            atm_desc["u22"] = U[1,1]
            atm_desc["u33"] = U[2,2]
            atm_desc["u12"] = U[0,1]
            atm_desc["u13"] = U[0,2]
            atm_desc["u23"] = U[1,2]

            ## calculate weight
            atm_desc["w"] = calc_atom_weight(atm)

        return xmlrpc_chain

    def iter_lsq_fit_messages(self, analysis, chain, min_span, xmlrpc_chain):
        """Iterate over (and create) the message dictionaries sent
        to consumer threads which are the the graph edge definitions
        of the TLS segments we need LSQ TLS fit data on.
        """
        for msg in iter_chain_subsegment_descs(chain, min_span):
            if analysis.tlsmdfile.grh_get_tls_record(
                msg["chain_id"], msg["frag_id1"], msg["frag_id2"])!=None:
                continue

            msg["method"] = "TLS"
            msg["xmlrpc_chain"] = xmlrpc_chain
            yield msg

    def process_chain(self, analysis, chain, min_span):
        """Fits TLS groups to all possible subsegments of the given Chain,
        storing the results in the analysis.grh_*() file.  
        """
        print "PROCESSING CHAIN: ",chain.chain_id

        ## convert the Chain object to a list of atoms
        xmlrpc_chain = self.convert_chain_for_xmlrpc(chain)

        ## create a iterator for the protein chain which yields
        ## a series of messages to call lsq_tls_fit on all the
        ## segments
        ##
        ## XXX: come up with a unique ID for the structure which
        ##      is included in the messages, so the return messages
        ##      from the pool can be checked
        iter_msg = self.iter_lsq_fit_messages(
            analysis, chain, min_span, xmlrpc_chain)

        iter_segment_fit_info = self.server_pool.iter_segment_processor(
            iter_msg)

        ## iterate over the fit_info dictionaries which hold the
        ## information on the TLS fits comming back from the grid servers
        num_edges          = 0
        num_rejected_edges = 0
        num_accepted_edges = 0
        
        for fit_info in iter_segment_fit_info:
            assert fit_info["chain_id"]==chain.chain_id

            num_edges += 1

            ## remove the reference to the xmlrpc_chain -- it's huge
            del fit_info["xmlrpc_chain"]

            ## save the fit_info record the graph state file
            analysis.tlsmdfile.grh_append_tls_record(fit_info)
            
            if fit_info.has_key("error"):
                num_rejected_edges += 1
                print "process_chain(chain_id=%s, frag_id={%s..%s}, "\
                      "discard=%s)" % (
                    chain.chain_id,
                    fit_info["frag_id1"], fit_info["frag_id2"],
                    fit_info["error"])
    
            else:
                num_accepted_edges += 1
                print "process_chain(chain_id=%s, frag_id={%s..%s}, "\
                      "lsqr=%6.4f)" % (
                    chain.chain_id,
                    fit_info["frag_id1"], fit_info["frag_id2"],
                    fit_info["lsq_residual"])

        print
        print "NUM POSSIBLE  TLS GROUPS: %d" % (num_edges)
        print "NUM ACCEPTED  TLS GROUPS: %d" % (num_accepted_edges)
        print "NUM DISCARDED TLS GROUPS: %d" % (num_rejected_edges)


class TLSOptimization(object):
    """Collection object containing one TLS description of a protein chain.
    """
    def __init__(self, chain, ntls_constraint):
        self.chain           = chain
        self.ntls_constraint = ntls_constraint
        self.ntls            = ntls_constraint
        self.tls_list        = []
        self.residual        = 0.0

    def is_valid(self):
        """Return True if the optimization is valid; otherwise, return False.
        """
        return len(self.tls_list)>0


class TLSChainMinimizer(HCSSSP):
    """Finds the minimal TLS description of a given Chain object using
    the HCSSP global optimization algorithm and a constraint on the number
    of TLS which can be used for the minimization.
    """
    def __init__(self, analysis, chain, min_span):
        self.analysis   = analysis
        self.min_span   = min_span

        self.minimized = False
        self.D         = None
        self.P         = None
        self.T         = None

        ## source and destination fragment IDs; these are the
        ## begining and ending fragments(residues) of the chain
        self.frag_id_src  = None
        self.frag_id_dest = None

        print "TLSChainMinimizer(chain_id=%s, range={%s..%s})" % (
            chain.chain_id, self.frag_id_src, self.frag_id_dest)

        if self.frag_id_src==None and self.frag_id_dest==None:
            self.chain = chain
        else:
            self.chain = chain[self.frag_id_src:self.frag_id_dest]
        
        ## Initalize the vertex list based on the number of of
        ## residues in Chain.
        self.num_vertex = len(self.chain) + 1

        ## calculate the mean number of atoms per residue in the
	## chain
	num_atoms = 0
	for atm in self.chain.iter_atoms():
            if calc_include_atom(atm):
                num_atoms += 1
	self.mean_res_atoms = round(float(num_atoms) / len(self.chain)) 

        ## calculate the minimum temperature factor in the chain
        ## XXX: hack
        self.min_temp_factor = None
	self.max_temp_factor = None
        for atm in self.chain.iter_all_atoms():
            if self.min_temp_factor==None:
                self.min_temp_factor = atm.temp_factor
		self.max_temp_factor = atm.temp_factor
                continue
            self.min_temp_factor = min(self.min_temp_factor, atm.temp_factor)
            self.max_temp_factor = max(self.max_temp_factor, atm.temp_factor)
        
    def run_minimization(self, max_tls_segments=20):
        """Run the HCSSSP minimization on the self.V,self.E graph, resulting
        in the creation of the self.D, self.P, and self.T arrays which
        contain 
        """
        ## build the vertex labels to reflect the protein structure
        ## the graph spans
        V = []
        for i in range(self.num_vertex):

            ## add the vertex label for i at Vi
            if i==0:
                vertex_label = "N-TERM"

            elif i==self.num_vertex - 1:
                vertex_label = "C-TERM"

            else:
                vertex_label = "%s{%s:%s}" % (
                    self.chain.chain_id,
                    self.chain[i-1].fragment_id,
                    self.chain[i].fragment_id)

            vertex_label = "V%d[%s]" % (i, vertex_label)
            V.append(vertex_label)

        self.V = V

        ## now build edges for the graph with weights given by the LSQ
        ## residual of TLS group fits
        grh_get_tls_record = self.analysis.tlsmdfile.grh_get_tls_record
        
        E = []
        for msg in iter_chain_subsegment_descs(self.chain, self.min_span):
            tls = grh_get_tls_record(
                self.chain.chain_id,
                msg["frag_id1"],
                msg["frag_id2"])
            
            if tls==None:
                print "build_graph(): [ERROR] no TLS group %s{%s..%s}" % (
                    self.chain.chain_id,
                    msg["frag_id1"],
                    msg["frag_id2"])

                sys.exit(-1)

            ## filter out the bad TLS segments
            if self.__minimization_filter(tls)==False:
                if True:
                    i, j = msg["vertex_i"], msg["vertex_j"]
		    num_res = j - i
                    if tls.has_key("num_atoms"):
                        weight = tls["num_atoms"] * (B2U * self.max_temp_factor)**2
                    else:
                        weight = num_res * self.mean_res_atoms * (B2U * self.max_temp_factor)**2
                    edge = (i, j, weight, (None, ))
                    E.append(edge)
                continue

            weight = tls["lsq_residual"]
            frag_range = (tls["frag_id1"], tls["frag_id2"])
                
            edge = (msg["vertex_i"], msg["vertex_j"], weight, frag_range)
            E.append(edge)

        self.E = E

        ## fill in any un-reachable gaps in the structure by adding
        ## fake 0.0 cost edges where needed
        ## XXX: FIXME

        ## perform the minimization
        start_timing()

        if len(self.E)>0:
            print "run_minimization(chain_id=%s): HCSSSP Minimizing..." % (
                self.chain.chain_id)
        
            D, P, T = self.HCSSSP_minimize(self.V, self.E, max_tls_segments)

            self.minimized = True
            self.D = D
            self.P = P
            self.T = T
        else:
            print "run_minimization(chain_id=%s): Unable to minimize" % (
                self.chain.chain_id)
            self.minimized = False

        ## free memory taken up from edges
        E      = None
        self.E = None
        import gc
        gc.collect()

        print "run_minimization(): ",end_timing()

    def __minimization_filter(self, tls):
        """Returns False if the tls group descibed by the tls database
        record should not be included in the minimization.
        """
        if tls.has_key("error"):
            return False
    
        T,L,S,O = self.__TLSO(tls)
        cdict = temp_calc_cor(T, L, S, O)

        ## sanity checks on rT
        rT = cdict["rT'"]
        evals = eigenvalues(rT)
        
        min_rT = U2B * min(evals)
        if min_rT < (0.1 * self.min_temp_factor):
            print "[FILTER] small rT chain=%s frag_rng={%s..%s}" % (
                tls["chain_id"], tls["frag_id1"], tls["frag_id2"])
            return False

        aniso = min(evals)/max(evals)
        if aniso<0.05:
            print "[FILTER] anisotropic rT chain=%s frag_rng={%s..%s}" % (
                tls["chain_id"], tls["frag_id1"], tls["frag_id2"])
            return False

        ## sanity checks on rho
        max_len = 50.0
        if length(cdict["L1_rho"])>max_len or \
           length(cdict["L2_rho"])>max_len or \
           length(cdict["L3_rho"])>max_len:
            print "[FILTER] length(rho) chain=%s frag_rng={%s..%s}" % (
                tls["chain_id"], tls["frag_id1"], tls["frag_id2"])
            return False
           
        return True

    def calc_tls_optimization(self, ntls_constraint):
        """Return a TLSOptimization() object containing the optimal
        TLS description of self.chain using num_tls_segments.
        """
        if not self.minimized:
            return None

        tlsopt = TLSOptimization(self.chain, ntls_constraint)
        
        for hi, hj, edge in self.HCSSSP_path_iter(
            self.V, self.D, self.P, self.T, ntls_constraint):

            if edge==None:
                continue
            
            ## check if the edge is a bypass-edge type
            i, j, weight, frag_range = edge

            tlsopt.residual += weight
            
            if len(frag_range)==2:
                tls = self.__calc_tls_record_from_edge(edge)
                tlsopt.tls_list.append(tls)

        return tlsopt

    def __calc_tls_record_from_edge(self, edge):
        """Independently calculate the TLS parameters for the segment
        to verify correctness (internal check).
        """
        i, j, weight, frag_range = edge

        ## retrieve from cache if possible
        frag_id1, frag_id2 = frag_range

        tls = self.analysis.tlsmdfile.grh_get_tls_record(
            self.chain.chain_id, frag_id1, frag_id2)

        segment  = self.chain[frag_id1:frag_id2]

        print "calc_tls_record_from_edge("\
              "chain_id=%s frag_id={%s..%s})" % (
            self.chain.chain_id, frag_id1, frag_id2)

        ## create TLSGroup
        tls_group = TLSGroup()

        ## add atoms to the group
        for atm in segment.iter_all_atoms():
            if calc_include_atom(atm)==True:
                tls_group.append(atm)

        ## take the TLS group tensors from the database set
        ## the TLSGroup object with them
        self.__add_tensors(tls_group, tls)

        ## helpful additions
        tls_info                    = tls_group.calc_tls_info()
        tls["tls_group"]            = tls_group
        tls["tls_info"]             = tls_info
        tls["segment"]              = segment
        tls["lsq_residual_per_res"] = tls["lsq_residual"] / (len(segment))

        return tls

    def __add_tensors(self, tls_group, tls):
        """Adds the TLS tensors from the tls database dictionary to the
        tls_group.
        """
        T,L,S,O = self.__TLSO(tls)
        tls_group.T      = T
        tls_group.L      = L
        tls_group.S      = S
        tls_group.origin = O
        
    def __TLSO(self, tls):
        """Builds Numeric Array tensors from the tls database dictionary.
        """
        O = array(
            [tls["cor_x"], tls["cor_y"], tls["cor_z"]], Float)

        T = array(
            [ [tls["t11"], tls["t12"], tls["t13"]],
              [tls["t12"], tls["t22"], tls["t23"]],
              [tls["t13"], tls["t23"], tls["t33"]] ], Float)
        
        L = array(
            [ [tls["l11"], tls["l12"], tls["l13"]],
              [tls["l12"], tls["l22"], tls["l23"]],
              [tls["l13"], tls["l23"], tls["l33"]] ], Float)
        
        s11, s22, s33 = calc_s11_s22_s33(tls["s2211"], tls["s1133"]) 
        
        S = array(
            [ [       s11, tls["s12"], tls["s13"]],
              [tls["s21"],        s22, tls["s23"]],
              [tls["s31"], tls["s32"],       s33] ], Float)

        return T,L,S,O
        
    def __add_bypass_edges(self):
        """Modify the graph before minimization to include the possibility
        of short segments of the protein which are not described well
        by the TLS model.
        """
        for i in range(self.num_vertex):
            for j in range(i, self.num_vertex):

                tls = {"method": "BYPASS"}

                weight = (j-i) * 100.0
                edge = (i, j, weight, tls)
                self.E.append(edge)
        
    def __add_bypass_edges_old(self):
        """Modify the graph before minimization to include the possibility
        of short segments of the protein which are not described well
        by the TLS model.
        """
        ## compute the mean weight/residue for all edges
        nres = 0
        mean = 0.0
        edge_mean_list = []

        for i,j,weight,tls in self.E:
            ## weight/residue for this edge
            edge_mean = weight / (j-i)
            edge_mean_list.append(edge_mean)

            ## for overall mean weight/residue
            nres += j-i
            mean += weight

        mean = mean / float(nres)

        ## write out histogram file
        num_bins = 100
        bins     = [0 for i in range(num_bins)]

        min_mean = min(edge_mean_list)
        max_mean = max(edge_mean_list)

        bin_range = max_mean - min_mean
        bin_width = bin_range / num_bins

        ## calculate the histogram if it's worth it to do so
        if bin_width>0.0:
            for edge_mean in edge_mean_list:
                bin_i = int((edge_mean - min_mean) / bin_width)

                try:
                    bins[bin_i] += 1
                except IndexError:
                    pass

            ## bin calculations for the mean value
            mean_bin_i = int( (mean - min_mean) / bin_width ) 

            hist = open("mean_weight_per_residue_histogram.dat", "w")

            hist.write("## min mean weight/residue:    %f\n" % (min_mean))
            hist.write("## max mean weight/residue:    %f\n" % (max_mean))
            hist.write("## mean weight/residue:        %f\n" % (mean))
            hist.write("## bin of mean weight/residue: %d\n" % (mean_bin_i)) 

            for i in range(len(bins)):
                hist.write("%10d   %10d\n" % (i+1, bins[i]))

            hist.close()

        ## create bypass edges with a weight equal to
        ## the average
        max_bypass_len = 20
    
        for i in range(self.num_vertex - max_bypass_len):
            for j in range(1, max_bypass_len+1):

                tls = {"method": "bypass"}

                weight = j * 100.0
                edge = (i, i+j, weight, tls)
                self.E.append(edge)
        
    def prnt_detailed_paths(self, hops=20):
        """Debug
        """
        if not self.minimized:
            return
        dest_j = len(self.V)-1

        for h in range(1, hops+1):
            print
            print "MINIMIZATON VERTEX PATH FOR %d SEGMENTS" % (h)
            print "NODE LABEL              HOPS      COST"\
                  "      PREVIOUS NODE          EDGE"

            self.__detailed_path(self.V, self.D, self.P, self.T, h)

    def __detailed_path(self, V, D, P, T, hop_constraint):
        """Print out the path from the source vertex (vertex 0) to
        the destination vertex (end vertex) given the hop_constarint.
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
                i, j, weight, frag_range = edge
                wr = weight / (j - i)
                edge_label = "(%3d,%3d,%6.3f,%s) %6.3f" % (
                    i, j, weight, frag_range, wr)
            else:
                edge_label = ""
                 
            print "%s   %3d     %10.4f   %s   %s" % (
                vertex_label, h, D[h,curr_v], prev_vertex_label,
                edge_label)

            curr_v = prev_vertex
            h -= 1


class TLSMDAnalysis(object):
    """Central object for a whole-structure TLS analysis.
    """
    def __init__(self,
                 struct_path    = None,
                 sel_chain_ids  = None,
                 tlsdb_file     = None,
                 tlsdb_complete = False,
                 gridconf_file  = None,
                 num_threads    = 1):

        set_globals()
        print_globals()

        self.struct_path     = struct_path
        if sel_chain_ids!=None:
            self.sel_chain_ids = sel_chain_ids.split(",")
        else:
            self.sel_chain_ids = None
        self.tlsdb_file      = tlsdb_file
        self.tlsdb_complete  = tlsdb_complete
        self.gridconf_file   = gridconf_file
        self.num_threads     = num_threads

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
            self.tlsdb_file = "%s.db" % (self.struct_id)
        self.tlsmdfile = TLSMDFile(self.tlsdb_file)

        ## select chains for analysis
        self.select_chains()

        ## print these settings
        self.prnt_settings()

        if not self.tlsdb_complete:
            self.calc_tls_segments()
        
        self.calc_chain_minimization()

    def prnt_settings(self):
        chain_ids = []
        for chain in self.chains:
            chain_ids.append(chain.chain_id)
        cids = string.join(chain_ids, ",")
        
        print "TLSMD ANALYSIS SETTINGS"
        print "    STRUCTURE FILE ================: %s"%(self.struct_path)
        print "    STRUCTURE ID ------------------: %s"%(self.struct_id)
        print "    CHAIN IDs SELECTED FOR ANALYSIS: %s"%(cids)
        print "    DATABASE FILE PATH ------------: %s"%(self.tlsdb_file)
        print "    GRID SERVER CONFIG FILE =======: %s"%(self.gridconf_file)
        print
        
    def load_struct(self):
        """Loads Structure, chooses a unique struct_id string.
        """
        print "LOADING STRUCTURE"
        print "    PATH: %s" % (self.struct_path)

        ## load struct
        self.struct = LoadStructure(
            fil = self.struct_path,
            build_properties = ("library_bonds","distance_bonds"))

        ## set the structure ID
        if GLOBALS.has_key("STRUCT_ID"):
            struct_id = GLOBALS["STRUCT_ID"]
        else:
            struct_id = self.struct.structure_id

        self.struct.structure_id = struct_id
        self.struct_id = struct_id

        print "    STRUCT ID: %s" % (self.struct_id)
        print

        ## if there are REFMAC5 TLS groups in the REMARK records of
        ## the PDB file, then add those in
        tls_file = TLSFile()
        tls_file.set_file_format(TLSFileFormatPDB())

        fil = open(self.struct_path, "r")
        tls_file.load(fil)

        if len(tls_file.tls_desc_list)>0:
            print "MERGING TLS GROUPS"
            print "    NUM TLS GROUPS: %d" % (len(tls_file.tls_desc_list))

            ## assume REFMAC5 groups where Utotal = Utls + Biso(temp_factor)
            for tls_desc in tls_file.tls_desc_list:
                tls_group = tls_desc.construct_tls_group_with_atoms(
                    self.struct)

                print "    TLS GROUP: %s" % (tls_group.name)
                
                for atm, Utls in tls_group.iter_atm_Utls():
                    bresi = atm.temp_factor
                    atm.temp_factor = bresi + (U2B * trace(Utls) / 3.0)
                    atm.U = (B2U * bresi * identity(3, Float)) + Utls
            
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
            if chain.count_amino_acids()<MIN_SUBSEGMENT_SIZE:
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
        
    def calc_tls_segments(self):
        """Calculates the TLSGraph for each chain in self.chains, optionally
        loading pre-computed graphs from the graph_file.  Any chains
        which need TLSGraphs computed will be strored in the graph_file.
        """
        chain_processor = self.launch_chain_processor()
        
        for chain in self.chains:
            begin_chain_timing(chain.chain_id)

            print "BUILDING TLS DATABASE FOR %s" % (chain)
            chain_processor.process_chain(self, chain, MIN_SUBSEGMENT_SIZE)

            end_chain_timing(chain.chain_id)

    def launch_chain_processor(self):
        """Starts up a in-process or grid server pool for graphing
        this structure.
        """
        ## setup the server pool for local in-process running (possibly
        ## using threads), or using a compute server grid
        server_pool = TLSGridServerPool()

        if self.gridconf_file!=None:
            ## open the grid config file and read the URLs of the compute
            ## server and add them to the server pool
            fil = open(self.gridconf_file, "r")

            ## each line is a comment starting with # or the URL of a
            ## grid server
            for ln in fil.readlines():
                server_url = ln.strip()
                if server_url.startswith("#"):
                    continue
                server_pool.launch_client_thread(server_url)

        else:
            ## launch local processing threads -- this is no good; they
            ## only run on one processor!
            for x in range(self.num_threads):
                server_pool.launch_client_thread(None)

        ## run the job to compute the full TLS Graph using the
        ## server_pool
        chain_processor = TLSChainProcessor(server_pool)
        return chain_processor

    def calc_chain_minimization(self):
        """Performs the TLS graph minimization on all TLSGraphs.
        """
        for chain in self.chains:
            chain.tls_chain_minimizer = TLSChainMinimizer(
                self, chain, MIN_SUBSEGMENT_SIZE)
            
            chain.tls_chain_minimizer.run_minimization()
            if not chain.tls_chain_minimizer.minimized:
                continue

            print
            print "="*79
            print "MINIMIZING CHAIN %s" % (chain)
            chain.tls_chain_minimizer.prnt_detailed_paths()
