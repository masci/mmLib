## TLS Motion Determination (TLSMD)
## Copyright 2002-2006 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import math
import numpy

from mmLib import Constants, Structure, FileIO

import const
import conf
import tree
import opt_containers


class Partition(object):
    def __init__(self, cpartition, lhs_tls, lhs_range, rhs_tls, rhs_range):
        self.cpartition = cpartition
        self.lhs_tls = lhs_tls
        self.lhs_range_idx = lhs_tls.segment_ranges.index(lhs_range)
        self.rhs_tls = rhs_tls
        self.rhs_range_idx = rhs_tls.segment_ranges.index(rhs_range)
        self.min_num_frags = 2

    def __repr__(self):
        return str(self)

    def __str__(self):
        lhs_frag_id1, lhs_frag_id2 = self.lhs_range()
        rhs_frag_id1, rhs_frag_id2 = self.rhs_range()
        return "%s-%s:%s-%s" % (
            lhs_frag_id1, lhs_frag_id2,
            rhs_frag_id1, rhs_frag_id2)

    def lhs_range(self):
        """The fragment range tuple of the left hand side segment.
        """
        return self.lhs_tls.segment_ranges[self.lhs_range_idx]

    def set_lhs_range(self, range):
        self.lhs_tls.segment_ranges[self.lhs_range_idx] = range

    def rhs_range(self):
        """The fragment range tuple of the right hand side segment.
        """
        return self.rhs_tls.segment_ranges[self.rhs_range_idx]

    def set_rhs_range(self, range):
        self.rhs_tls.segment_ranges[self.rhs_range_idx] = range

    def move_left(self, num_frags = 1):
        """Move the partition point num_frags to the left.
        """
        chain = self.cpartition.chain
        
        lhs_frag_id1, lhs_frag_id2 = self.lhs_range()
        lhs_frag1 = chain.get_fragment(lhs_frag_id1)
        lhs_frag2 = chain.get_fragment(lhs_frag_id2)
        lhs_nfrag = lhs_frag2.ifrag - lhs_frag1.ifrag + 1
        if lhs_nfrag < self.min_num_frags:
            raise ValueError
        new_lhs_frag2 = chain[lhs_frag2.ifrag - num_frags]
        new_lhs_range = (lhs_frag_id1, new_lhs_frag2.fragment_id)

        rhs_frag_id1, rhs_frag_id2 = self.rhs_range()
        rhs_frag1 = chain.get_fragment(rhs_frag_id1)
        rhs_frag2 = chain.get_fragment(rhs_frag_id2)
        rhs_nfrag = rhs_frag2.ifrag - rhs_frag1.ifrag + 1
        if rhs_nfrag < self.min_num_frags:
            raise ValueError
        new_rhs_frag1 = chain[rhs_frag1.ifrag - num_frags]
        new_rhs_range = (new_rhs_frag1.fragment_id, rhs_frag_id2)

        self.set_lhs_range(new_lhs_range)
        self.set_rhs_range(new_rhs_range)

    def move_right(self, num_frags = 1):
        """Move the partition point num_frgs to the right.
        """
        self.move_left(-num_frags)


def ChainPartitionList(cpartition):
    segment_range_list = []
    for tls in cpartition.iter_tls_segments():
        for range in tls.iter_segment_ranges():
            segment_range_list.append((range, tls))

    def srlcmp(item1, item2):
        return Structure.fragment_id_cmp(item1[0][0], item2[0][0])
    segment_range_list.sort(srlcmp)

    partitions = []
    prev_sr = None
    for sr in segment_range_list:
        if prev_sr is None:
            prev_sr = sr
            continue
        prev_range, prev_tls = prev_sr
        range, tls = sr
        partition = Partition(cpartition, prev_tls, prev_range, tls, range)
        partitions.append(partition)
        prev_sr = sr
        
    return partitions

    
def RefineChainPartitionPositions(cpartition):
    """Refines positions of the partition points of the
    given ChainPartition instance.
    """
    partitions = ChainPartitionList(cpartition)
    print partitions


## testing
def testmain():
    import tlsmd_analysis
    
    struct = FileIO.LoadStructure(fil="/home/jpaint/PDB/8rxn.pdb")
    chain = tlsmd_analysis.ConstructSegmentForAnalysis(struct.get_chain("A"))
    
    cpartition = opt_containers.ChainPartition(chain, 3)

    groups = [
        [("1", "15"), ("35", "52")],
        [("16", "21")],
        [("22", "34")] ]

    for segment_ranges in groups:
        tls = opt_containers.TLSSegment(segment_ranges)
        cpartition.add_tls_segment(tls)
    
    partitions = ChainPartitionList(cpartition)
    print partitions

    print "Move 0 left"
    for x in range(3):
        partitions[0].move_left()
        print partitions

    print "Move 1 right"
    for x in range(5):
        partitions[1].move_right()
        print partitions

    

if __name__ == "__main__":
    testmain()
