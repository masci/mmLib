## TLS Motion Determination (TLSMD)
## Copyright 2002-2006 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import math
import numpy

from mmLib import Constants

import const
import conf
import opt_containers


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
    tls = TLSSegment(chain_id = tls1.chain_id,
                     segment_ranges = tls1.segment_ranges + tls2.segment_ranges,
                     residual = tlsdict["residual"],
                     method = tls1.method,
                     num_atoms = tlsdict["num_atoms"])
    tls.reset(chain)
    return tls


def ChainPartitionRecombination(cpartition):
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
        tlsdict = tls_analyzer.isotropic_fit(segment_ranges)
        chain_residual = tlsdict["residual"]
        
        for itls, tls in enumerate(tls_list):
            if itls == part1 or itls == part2:
                continue
            chain_residual += tls.residual()

        recombination_list.append((chain_residual, tls1, tls2, tlsdict))

    recombination_list.sort()
    chain_residual, tls1, tls2, tlsdict = recombination_list[0]
    tls12 = JoinTLSSegments(tls1, tls2, cpartition.chain, tlsdict)

    combined_cp = opt_container.ChainPartition(cpartition.chain, cpartition.ntls - 1)
    for tls in tls_list:
        if tls == tls1:
            combined_cp.add_tls_segment(tls12)
            continue
        if tls == tls2:
            continue
        combined_cp.add_tls_segment(tls)

    assert numpy.allclose(chain_residual, combined_cp.residual())
    return combined_cp


def ChainPartitionRecombinationOptimization(chain):
    improvements = []
    prev_cpartition = None
    for ntls, cpartition in chain.partition_collection.iter_ntls_chain_partitions():
        if prev_cpartition is None:
            prev_cpartition = cpartition
            continue

        print
        print "%d->%d Segment Recombination beat(%8.6f)" % (ntls, ntls-1, prev_cpartition.residual())
        combined_cp = ChainPartitionRecombination(cpartition)
        print "%s %8.6f" % (combined_cp, combined_cp.residual())

        if combined_cp.residual() < prev_cpartition.residual():
            improvements.append(combined_cp)

        prev_cpartition = cpartition

    ## insert replacement ChainPartitions
    for cpartition in improvements:
        chain.partition_collection.insert_chain_partition(cpartition)



