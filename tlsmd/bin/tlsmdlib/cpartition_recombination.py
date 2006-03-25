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
            if frag2.ichain == (xfrag1.ichain - 1):
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
    chain = cpartition.chain
    tls_list = cpartition.tls_list
    tls_analyzer = chain.tls_analyzer
    nparts = len(tls_list)

    recombination_list = []

    for part1, part2 in recombination2_iter(nparts):
        tls1 = tls_list[part1]
        tls2 = tls_list[part2]

        segment_ranges = tls1.segment_ranges + tls2.segment_ranges
        segment_ranges.sort(segment_range_cmp)
        tlsdict = tls_analyzer.isotropic_fit(segment_ranges)
        chain_residual = tlsdict["residual"]
        
        for tls in tls_list:
            if tls != tls1 and tls != tls2:
                chain_residual += tls.residual()

        recombination_list.append((chain_residual, tls1, tls2, tlsdict))

    recombination_list.sort()
    return_list = []
    for i, reco in enumerate(recombination_list):
        if i >= num_return:
            break
        chain_residual, tls1, tls2, tlsdict = reco
        tls12 = JoinTLSSegments(tls1, tls2, cpartition.chain, tlsdict)

        combined_cp = opt_containers.ChainPartition(cpartition.chain, cpartition.num_tls_segments() - 1)
        return_list.append(combined_cp)

        for tls in tls_list:
            if tls == tls1:
                combined_cp.add_tls_segment(tls12)
                continue
            if tls == tls2:
                continue
            combined_cp.add_tls_segment(tls.copy())

        assert numpy.allclose(chain_residual, combined_cp.residual())

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
    ntls_best = {}
    orig_best = {}
    for ntls, cpartition in chain.partition_collection.iter_ntls_chain_partitions():
        ntls_best[ntls] = cpartition
        orig_best[ntls] = cpartition

    for ntls, cpartition in chain.partition_collection.iter_ntls_chain_partitions():
        if ntls < 2:
            continue
        
        print "%d INTO %d TO 2" % (ntls, ntls - 1)

        search_width = 3
        search_depth = 3

        proot = tree.Tree()
        proot.cp = cpartition

        while True:
            ExtendRecombinationTree(proot, search_depth, search_width)

            max_depth = proot.depth()
            best_at_depth = None

            for depth, ptree in proot.iter_depth_first():
                if depth == max_depth:
                    if best_at_depth is None or ptree.cp.residual() < best_at_depth.cp.residual():
                        best_at_depth = ptree

                tmp_ntls = ptree.cp.num_tls_segments()
                if (ptree.cp.rmsd_b() + 0.05) < orig_best[tmp_ntls].rmsd_b():
                    if ptree.cp.residual() < ntls_best[tmp_ntls].residual():
                        print "%s %5.2f(%5.2f) %d" % (ptree.cp,
                                                      ptree.cp.rmsd_b(),
                                                      ntls_best[tmp_ntls].rmsd_b(),
                                                      tmp_ntls)
                        ntls_best[tmp_ntls] = ptree.cp
                    

            ptree = best_at_depth
            for i in xrange(max_depth - 1):
                ptree = ptree.parent()

            proot = ptree
            if search_depth > max_depth:
                break

    ## insert replacement ChainPartitions
    for ntls, cpartition in ntls_best.iteritems():
        cp = chain.partition_collection.get_chain_partition(ntls)
        if cp != cpartition:
            chain.partition_collection.insert_chain_partition(cpartition)



