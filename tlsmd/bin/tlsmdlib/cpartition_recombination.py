## TLS Motion Determination (TLSMD)
## Copyright 2002-2006 by TLSMD Development Group (see AUTHORS file)
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

## TODO Allow multiple chains. Christoph Champ, 2008-03-20
#def ChainPartitionRecombination(cpartition, num_return = 1):


def ChainPartitionRecombination(cpartition, num_return = 1):
    """Returns a new TLSSegment object which is the best
    combination of any two TLSSegment instances in cpartition.
    """
    chain = cpartition.chain
    tls_list = cpartition.tls_list
    tls_analyzer = chain.tls_analyzer
    nparts = len(tls_list)

    ## set the diagonal to the rmsd_b of the individual groups
    rmsd_b_mtx = numpy.zeros((nparts, nparts), float)
    for i, tls in enumerate(tls_list):
        rmsd_b_mtx[i,i] = tls.residual_rmsd_b()

    recombination_list = []

    for i, j in recombination2_iter(nparts):
        tls1 = tls_list[i]
        tls2 = tls_list[j]

        segment_ranges = tls1.segment_ranges + tls2.segment_ranges
        segment_ranges.sort(segment_range_cmp)
        tlsdict = tls_analyzer.isotropic_fit(segment_ranges)

        residual = tlsdict["residual"]
        num_residues = tlsdict["num_residues"]
        msd = residual / num_residues
        rmsd_b = Constants.U2B * math.sqrt(msd)

        rmsd_b_mtx[i, j] = rmsd_b
        rmsd_b_mtx[j, i] = rmsd_b
        
        residual_delta = (residual - (tls1.residual() + tls2.residual())) / num_residues
        recombination_list.append((residual_delta, tls1, tls2, tlsdict))

    cpartition.rmsd_b_mtx = rmsd_b_mtx

    recombination_list.sort()
    return_list = []
    for i, reco in enumerate(recombination_list):
        if i >= num_return:
            break
        residual_delta, tls1, tls2, tlsdict = reco
        tls12 = JoinTLSSegments(tls1, tls2, cpartition.chain, tlsdict)

        combined_cp = opt_containers.ChainPartition(cpartition.chain, cpartition.num_tls_segments() - 1)
        return_list.append(combined_cp)
        combined_cp.residual_delta = residual_delta

        for tls in tls_list:
            if tls == tls1:
                combined_cp.add_tls_segment(tls12)
                continue
            if tls == tls2:
                continue
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
    visited = {}
    
    ntls_best = {}
    orig_best = {}
    for ntls, cpartition in chain.partition_collection.iter_ntls_chain_partitions():
        ntls_best[ntls] = cpartition
        orig_best[ntls] = cpartition

    for ntls, cpartition in chain.partition_collection.iter_ntls_chain_partitions():
        if ntls < 2:
            continue
        
        console.stdoutln("%d INTO %d TO 2" % (ntls, ntls - 1))

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

                    console.stdoutln("%s %5.2f(%s) %d" % (ptree.cp, ptree.cp.rmsd_b(), best_rmsd, tmp_ntls))

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

