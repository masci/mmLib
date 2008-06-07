## TLS Motion Determination (TLSMD)
## Copyright 2002-2008 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import math
import numpy

from mmLib import Constants, Structure

import const
import conf
import console
import tree
import opt_containers


def segment_range_cmp(segrange1, segrange2):
    return Structure.fragment_id_cmp(segrange1[0], segrange2[0])


def join_segment_ranges(chain, segrange1, segrange2):
    tmp_segrange = segrange1 + segrange2
    tmp_segrange.sort(segment_range_cmp)

    if len(tmp_segrange) == 0:
        return []
    
    segrange = []
    frag1 = None
    frag2 = None
    
    for fid1, fid2 in tmp_segrange:
        xfrag1 = chain[fid1]
        xfrag2 = chain[fid2]
        if frag2 is not None:
            if frag2.ifrag == (xfrag1.ifrag - 1):
                frag2 = xfrag2
                continue
            segrange.append((frag1.fragment_id, frag2.fragment_id))
        frag1 = xfrag1
        frag2 = xfrag2
    
    segrange.append((frag1.fragment_id, frag2.fragment_id))
    return segrange
    

def recombination2_iter(nparts):
    for part1 in xrange(0, nparts):
        for part2 in xrange(part1 + 1, nparts):
            yield part1, part2


def JoinTLSSegments(tls1, tls2, chain, tlsdict):
    """Returns a new TLSSegment object which is the combination
    of TLSSegment objects tls1 and tls2.  The argument tlsdict is
    the results from a tls_analyzer fit of the combined tls1 and tls2
    atoms, which should already have been performed.
    """

    assert tls1.chain_id == tls2.chain_id
    segment_ranges = join_segment_ranges(chain, tls1.segment_ranges, tls2.segment_ranges)
    tls = opt_containers.TLSSegment(chain_id = tls1.chain_id,
                                    segment_ranges = segment_ranges,
                                    residual = tlsdict["residual"],
                                    method = tls1.method,
                                    num_atoms = tlsdict["num_atoms"],
                                    num_residues = tlsdict["num_residues"])
    return tls


def ChainPartitionRecombination(cpartition, num_return = 1):
    """Returns a new TLSSegment object which is the best
    combination of any two TLSSegment instances in cpartition.
    """

    chain = cpartition.chain 		## E.g., "Segment(1:A, Res(ILE,16,A)...Res(SER,116,A))"
    tls_list = cpartition.tls_list	## 2D object
    tls_analyzer = chain.tls_analyzer	## 1D object
    nparts = len(tls_list)		## integer

    ## set the diagonal to the rmsd_b of the individual groups
    rmsd_b_mtx = numpy.zeros((nparts, nparts), float)
    for i, tls in enumerate(tls_list):
	## E.g., "i=0, tls=A:16-75; 82-94; 101-116"
        rmsd_b_mtx[i,i] = tls.residual_rmsd_b()

    recombination_list = []

    for i, j in recombination2_iter(nparts):
	## NOTE (by Christoph):
	## i and j are integer values and each combination is considered.
	## i starts at 0 and j starts at 1 even though it is a MxN matrix (not sure why?)
	## E.g., if i=[0,1,2] and j=[1,2,3], then:
	##	i=0, j=1
	##	i=0, j=2
	##	i=0, j=3
	##	i=1, j=2
	##	i=1, j=3
	##	i=2, j=3
	## which is the non-redundant values in a square matrix, e.g.,
	##	    1 2 3=j
	##	i=0 x x x
	##	  1   x x
	##	  2     x

        tls1 = tls_list[i]
        tls2 = tls_list[j]

        segment_ranges = tls1.segment_ranges + tls2.segment_ranges
        segment_ranges.sort(segment_range_cmp)
        tlsdict = tls_analyzer.isotropic_fit(segment_ranges)

        residual = tlsdict["residual"]
        num_residues = tlsdict["num_residues"]
        msd = residual / num_residues
        rmsd_b = Constants.U2B * math.sqrt(msd)

	## NOTE (by Christoph):
	## rmsd_b is the off-diagonal values found in the 'xxxx_CHAINa_NTLSn_RECOMBINATION.txt'
	## files, where xxxx = pdbID, a = chainID, n = TLS segment
        rmsd_b_mtx[i, j] = rmsd_b
        rmsd_b_mtx[j, i] = rmsd_b
        
        residual_delta = (residual - (tls1.residual() + tls2.residual())) / num_residues
        recombination_list.append((residual_delta, tls1, tls2, tlsdict))

    cpartition.rmsd_b_mtx = rmsd_b_mtx

    recombination_list.sort()
    lenreco=len(recombination_list) ## DEBUG
    return_list = []
    for i, reco in enumerate(recombination_list):
	## NOTE (by Christoph):
	## For some reason, i always equals zero?
	## It shouldn't get past the break
        if i >= num_return:
            break

        residual_delta, tls1, tls2, tlsdict = reco
	## NOTE (by Christoph):
	## float(residual_delta),
	## tls1=A:16-75,  # <-example
	## tls2=A:76-116, # <-example
	## tls1+tls2 => (A:16-75; A:76-116) # <- they produce these types of "TLS SEGMENT RECOMBINATION" in 'log.txt'
	## tlsdict={	'residual':foat, 
	##		'is3':float, 
	##		'il11':float, 
	##		'il12':float, 
	##		'il13':float, 
	##		'it':float, 
	##		'num_residues':int, 
	##		'il23':float, 
	##		'il22':float, 
	##		'il33':float, 
	##		'is2':float, 
	##		'y':float, 
	##		'x':float,
	##		'z':float, 
	##		'is1':float, 
	##		'num_atoms':int}
	#console.stdoutln(">>>i=%s, len(reco)=%s, residual_delta=%s, tls1=%s, tls2=%s"%(i,lenreco,residual_delta, tls1, tls2)) ## DEBUG
        tls12 = JoinTLSSegments(tls1, tls2, cpartition.chain, tlsdict)

        combined_cp = opt_containers.ChainPartition(cpartition.chain, cpartition.num_tls_segments() - 1)
        return_list.append(combined_cp)
        combined_cp.residual_delta = residual_delta

        for tls in tls_list:
	    ## NOTE (by Christoph):
	    ## In 'log.txt', an example line might look like:
	    ## "(A:16-26)(A:27-61; 70-75)(A:62-69)(A:76-94)(A:95-105)(A:106-116)  5.10( 5.09) 6"
	    ## For fields having double segment, e.g., "(A:27-61; 70-75)":
            if tls == tls1:
		## Creates the _first_ part of, e.g., "(A:27-61; 70-75)". I.e., "A:27-61"
                combined_cp.add_tls_segment(tls12)
                continue
            if tls == tls2:
		## Creates the _last_ part of, e.g., "(A:27-61; 70-75)". I.e., "70-75"
                continue
	    ## For fields _not_ having double segment, e.g., "(A:16-26)":
            combined_cp.add_tls_segment(tls.copy())

    return return_list


def ExtendRecombinationTree(ptree, depth, width):

    if depth == 0:
        return
    if ptree.empty():
        combined_list = ChainPartitionRecombination(ptree.cp, width)
        for combined_cp in combined_list:
            child = tree.Tree()
            child.cp = combined_cp
            ptree.append(child)
            ExtendRecombinationTree(child, depth - 1, width)
    else:
        for child in ptree.iter_children():
            ExtendRecombinationTree(child, depth - 1, width)


def ChainPartitionRecombinationOptimization(chain):
    console.debug_stdoutln(">Entering: cpartition_recombination.py->ChainPartitionRecombinationOptimization()...") ## DEBUG

    visited = {}
    
    ntls_best = {}
    orig_best = {}
    for ntls, cpartition in chain.partition_collection.iter_ntls_chain_partitions():
        ntls_best[ntls] = cpartition
        orig_best[ntls] = cpartition

    for ntls, cpartition in chain.partition_collection.iter_ntls_chain_partitions():
        if ntls < 2:
            continue
        
        console.stdoutln("%d INTO %d TO 2" % (ntls, ntls - 1)) ## LOGLINE

        search_width = 1
        search_depth = 1

        proot = tree.Tree()
        proot.cp = cpartition

        while True:
            ExtendRecombinationTree(proot, search_depth, search_width)

            max_depth = proot.depth()
            best_at_depth = None

            for depth, ptree in proot.iter_depth_first():
                if depth == max_depth:
                    if best_at_depth is None or ptree.cp.residual_delta < best_at_depth.cp.residual_delta:
                        best_at_depth = ptree

                tmp_ntls = ptree.cp.num_tls_segments()

                if not visited.has_key(ptree):
                    visited[ptree] = True

                    if ntls_best.has_key(tmp_ntls):
                        best_rmsd = "%5.2f" % (ntls_best[tmp_ntls].rmsd_b())
                    else:
                        best_rmsd = "-----"

		    ## E.g., "(A:2-11; 18-23; 39-52; 58-95)(A:12-17; 24-38; 53-57)  5.42( 4.98) 2"
                    console.stdoutln("%s %5.2f(%s) %d" % (ptree.cp, ptree.cp.rmsd_b(), best_rmsd, tmp_ntls)) ## LOGLINE

                if ntls_best.has_key(tmp_ntls):
                    if ptree.cp.residual() < ntls_best[tmp_ntls].residual():
                        ntls_best[tmp_ntls] = ptree.cp
                else:
                    ntls_best[tmp_ntls] = ptree.cp

            ptree = best_at_depth
            for i in xrange(max_depth - 1):
                ptree = ptree.parent()

            proot = ptree
            if search_depth > max_depth:
                break
            
    ## insert replacement ChainPartitions
    if conf.globalconf.recombination:
        for ntls, cpartition in ntls_best.iteritems():
            cp = chain.partition_collection.get_chain_partition(ntls)
            if cp != cpartition:
                chain.partition_collection.insert_chain_partition(cpartition)

