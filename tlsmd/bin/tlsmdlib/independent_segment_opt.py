## TLS Motion Determination (TLSMD)
## Copyright 2002-2006 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import math
import numpy
import gc

from mmLib import Constants, TLS

import misc
import const
import conf
import console
import atom_selection
import hcsssp
import tls_calcs
import tlsmdmodule
import opt_containers


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
    frag_ids = []
    for frag in chain.iter_fragments():
        frag_ids.append(frag.fragment_id)

    num_vertex = len(frag_ids) + 1    
    for vertex_i, vertex_j in iter_ij(num_vertex, min_len):
        frag_id1 = frag_ids[vertex_i]
        frag_id2 = frag_ids[vertex_j-1]
        yield frag_id1, frag_id2, vertex_i, vertex_j

class ISOptimization(hcsssp.HCSSSP):
    """Finds the minimal TLS description of a given Chain instance using
    the HCSSSP global optimization algorithm and a constraint on the number
    of TLS which can be used for the optimization.
    """
    def __init__(self, chain, min_subsegment_len, nparts):
        hcsssp.HCSSSP.__init__(self)

        self.chain = chain
        self.min_subsegment_len = min_subsegment_len
        self.nparts = nparts

        self.minimized = False
        self.D = None
        self.P = None
        self.T = None

    def get_fit_method(self, chain):
        fit_method = None
        if conf.globalconf.tls_model == "ISOT":
            fit_method = chain.tls_analyzer.isotropic_fit_segment
        elif conf.globalconf.tls_model == "ANISO":
            fit_method = chain.tls_analyzer.anisotropic_fit_segment
        elif conf.globalconf.tls_model=="NLISOT":
            fit_method = chain.tls_analyzer.constrained_isotropic_fit_segment
        elif conf.globalconf.tls_model=="NLANISO":
            fit_method = chain.tls_analyzer.constrained_anisotropic_fit_segment
        return fit_method

    def run_minimization(self):
        """Run the HCSSSP minimization on the self.V,self.E graph, resulting
        in the creation of the self.D, self.P, and self.T arrays which
        contain 
        """
        chain = self.chain
        chain_id = self.chain.chain_id
        min_subsegment_len = self.min_subsegment_len
        num_vertex = len(chain) + 1

        ## choose the TLS Model to fit for the chain
        fit_method = self.get_fit_method(chain)

        ## build the vertex labels to reflect the protein structure
        ## the graph spans
        vertices = []
        for i in xrange(num_vertex):
            ## add the vertex label for i at Vi
            if i == 0 :
                vertex_label = "N-TERM"
            elif i == num_vertex - 1:
                vertex_label = "C-TERM"
            else:
                vertex_label = "%s{%s:%s}" % (chain_id, chain[i-1].fragment_id, chain[i].fragment_id)
            vertex_label = "V%d[%s]" % (i, vertex_label)
            vertices.append(vertex_label)

        ## fit chain segments with TLS model and build residual graph to minimize
        total_num_subsegments = calc_num_subsegments(chain.count_fragments(), min_subsegment_len)
        num_subsegments = 0
        pcomplete = 0
        pcomplete_old = 0
        edges = []
        for frag_id1, frag_id2, i, j in iter_chain_subsegment_descs(chain, min_subsegment_len):
            tlsdict = fit_method(frag_id1, frag_id2)

            num_subsegments += 1
            pcomplete = round(100.0 * num_subsegments / total_num_subsegments)
            if pcomplete != pcomplete_old:
                console.stdoutln("(%10d/%10d) %2d%% Complete" % (num_subsegments, total_num_subsegments, pcomplete))
                pcomplete_old = pcomplete

            if tlsdict == None:
                console.stderrln("no TLS group %s{%s..%s}" % (chain_id, frag_id1, frag_id2))
                raise SystemExit
            if tlsdict.has_key("error") is True:
                continue
            if not tlsdict.has_key("residual"):
                console.stderrln("no residual! %s{%s..%s}" % (chain_id, frag_id1, frag_id2))
                raise SystemExit

            residual = tlsdict["residual"]
            num_atoms = tlsdict["num_atoms"]
            num_residues = tlsdict["num_residues"]
            msd = residual / num_residues
            rmsd = math.sqrt(msd)
            rmsd_b = rmsd * Constants.U2B
            chi2 = msd * num_atoms

            if num_atoms < 40:
                continue

            cost = residual

            frag_range = (frag_id1, frag_id2)
            edge = (i, j, cost, frag_range, tlsdict)
            edges.append(edge)

        ## perform the minimization
        if len(edges) > 0:
            console.stdoutln("run_minimization(chain_id=%s): HCSSSP Minimizing..." % (chain_id))
        
            D, P, T = self.HCSSSP_minimize(vertices, edges, self.nparts)

            self.minimized = True
            self.V = vertices
            self.D = D
            self.P = P
            self.T = T
        else:
            console.stdoutln("run_minimization(chain_id=%s): Unable to minimize" % (chain_id))
            self.minimized = False

        ## free memory taken up from edges
        edges = None
        gc.collect()

    def construct_tls_segment(self, edge):
        """Returns a instance of TLSSegment fully constructed for
        self.chain and the fragment range given in edge.
        """
        i, j, cost, frag_range, tlsdict = edge
        tls = opt_containers.TLSSegment(chain_id = self.chain.chain_id,
                                        segment_ranges = [frag_range],
                                        method = "TLS",
                                        residual = cost,
                                        num_atoms = tlsdict["num_atoms"],
                                        num_residues = tlsdict["num_residues"])
        return tls

    def construct_chain_partition(self, nparts):
        """Return a ChainPartition instance containing the optimal
        TLS description of self.chain using num_tls_segments.
        """
        if not self.minimized:
            return None

        cpartition = opt_containers.ChainPartition(self.chain, nparts)
        
        for hi, hj, edge in self.HCSSSP_path_iter(self.V, self.D, self.P, self.T, nparts):
            if edge is None:
                continue
            i, j, cost, frag_range, tlsdict = edge
            
            ## check if the edge is a bypass-edge type
            if len(frag_range) == 2:
                tls_segment = self.construct_tls_segment(edge)
                cpartition.add_tls_segment(tls_segment)

        return cpartition

    def construct_partition_collection(self, nparts_max = None):
        """Returns a ChainPartitionCollection instance containing
        the optimized ChainPartition instances using from 1 to
        nparts_max TLS groups.
        """
        if nparts_max is None:
            nparts_max = self.nparts

        partition_collection = opt_containers.ChainPartitionCollection(self.chain)
        for nparts in xrange(1, nparts_max + 1):
            cpartition = self.construct_chain_partition(nparts)
            if cpartition is not None:
                partition_collection.insert_chain_partition(cpartition)

        return partition_collection

    def prnt_detailed_paths(self):
        """Debug
        """
        hops = self.nparts
        
        if not self.minimized:
            return
        dest_j = len(self.V) - 1

        for h in xrange(1, hops + 1):
            console.endln()
            console.stdoutln("MINIMIZATON VERTEX PATH FOR %d SEGMENTS" % (h))
            console.stdoutln("NODE LABEL              HOPS      COST      PREVIOUS NODE          EDGE")
            self.__detailed_path(self.V, self.D, self.P, self.T, h)

    def __detailed_path(self, V, D, P, T, hop_constraint):
        """Print out the path from the source vertex (vertex 0) to
        the destination vertex (end vertex) given the hop_constraint.
        """
        num_vertex = len(D[0])
        
        ## start at the destination vertex
        curr_v = num_vertex - 1
        h = hop_constraint
        
        while curr_v >= 0:
            prev_vertex  = P[h,curr_v]
            vertex_label = V[curr_v].ljust(20)

            if prev_vertex < 0:
                prev_vertex_label = "".ljust(20)
            else:
                prev_vertex_label = V[prev_vertex].ljust(20)

            edge = T[h][curr_v]

            if edge is not None:
                i, j, cost, frag_range, tlsdict = edge
                wr = cost / (j - i)
                edge_label = "(%3d,%3d,%6.3f,%s) %6.3f" % (i, j, cost, frag_range, wr)
            else:
                edge_label = ""

            console.stdoutln(
                "%s   %3d     %10.4f   %s   %s" % (
                vertex_label, h, D[h,curr_v], prev_vertex_label, edge_label))

            curr_v = prev_vertex
            h -= 1
