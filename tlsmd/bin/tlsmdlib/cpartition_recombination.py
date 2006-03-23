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
import opt_containers


def segment_range_cmp(range1, range2):
    return Structure.fragment_id_cmp(range1[0], range2[0])

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

    segment_ranges = tls1.segment_ranges + tls2.segment_ranges
    segment_ranges.sort(segment_range_cmp)
    
    tls = opt_containers.TLSSegment(chain_id = tls1.chain_id,
                                    segment_ranges = segment_ranges,
                                    residual = tlsdict["residual"],
                                    method = tls1.method,
                                    num_atoms = tlsdict["num_atoms"])
    tls.fit_to_chain(chain)
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
    for i in range(num_return):
        chain_residual, tls1, tls2, tlsdict = recombination_list[i]
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


def ChainPartitionRecombinationOptimization(chain):
    ntls_best = {}
    for ntls, cpartition in chain.partition_collection.iter_ntls_chain_partitions():
        ntls_best[ntls] = cpartition

    for ntls, cpartition in chain.partition_collection.iter_ntls_chain_partitions():
        if ntls < 2:
            continue
        
        print "%d INTO %d TO 2" % (ntls, ntls - 1)

        tmp_cp = cpartition
        while tmp_cp.num_tls_segments() > 1:

            combined_list = ChainPartitionRecombination(tmp_cp)
            for combined_cp in combined_list:
                tmp_ntls = combined_cp.num_tls_segments()

                print "%s %8.6f(%8.6f) %d" % (combined_cp,
                                              combined_cp.residual(),
                                              ntls_best[tmp_ntls].residual(),
                                              tmp_ntls)

                if combined_cp.residual() < ntls_best[tmp_ntls].residual():
                    ntls_best[tmp_ntls] = combined_cp

                tmp_cp = combined_cp

    ## insert replacement ChainPartitions
    for ntls, cpartition in ntls_best.iteritems():
        cp = chain.partition_collection.get_chain_partition(ntls)
        if cp != cpartition:
            chain.partition_collection.insert_chain_partition(cpartition)



