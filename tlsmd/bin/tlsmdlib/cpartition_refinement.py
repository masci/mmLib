## TLS Motion Determination (TLSMD)
## Copyright 2002-2009 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

## Python modules
import math
import numpy

## Pymmlib
from mmLib import Constants, Structure, FileIO

## TLSMD
import const
import conf
import console
import opt_containers


class Partition(object):
    """A Partition instance represents the interface of two adjacent
    TLS segments.
    """

    def __init__(self, cpartition, lhs_tls, lhs_range, rhs_tls, rhs_range, auto_fit_residual = True):
        assert isinstance(cpartition, opt_containers.ChainPartition)
        assert isinstance(lhs_tls, opt_containers.TLSSegment)
        assert isinstance(rhs_tls, opt_containers.TLSSegment)

        self.cpartition = cpartition
        self.lhs_tls = lhs_tls
        self.lhs_range_idx = lhs_tls.segment_ranges.index(lhs_range)
        self.rhs_tls = rhs_tls
        self.rhs_range_idx = rhs_tls.segment_ranges.index(rhs_range)
        self.min_num_frags = 4
        self.auto_fit_residual = auto_fit_residual

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

    def move(self, num_frags):
        """Move the partition point num_frags.
        """
        chain = self.cpartition.chain

        lhs_frag_id1, lhs_frag_id2 = self.lhs_range()
        rhs_frag_id1, rhs_frag_id2 = self.rhs_range()

        lhs_frag1 = chain.get_fragment(lhs_frag_id1)
        lhs_frag2 = chain.get_fragment(lhs_frag_id2)
        rhs_frag1 = chain.get_fragment(rhs_frag_id1)
        rhs_frag2 = chain.get_fragment(rhs_frag_id2)

        new_lhs_ifrag2 = lhs_frag2.ifrag + num_frags
        new_rhs_ifrag1 = rhs_frag1.ifrag + num_frags

        if new_lhs_ifrag2 < (lhs_frag1.ifrag + (self.min_num_frags-1)) or \
           new_rhs_ifrag1 > (rhs_frag2.ifrag - (self.min_num_frags-1)):
            raise ValueError

        new_lhs_frag2 = chain[new_lhs_ifrag2]
        new_lhs_range = (lhs_frag_id1, new_lhs_frag2.fragment_id)

        new_rhs_frag1 = chain[new_rhs_ifrag1]
        new_rhs_range = (new_rhs_frag1.fragment_id, rhs_frag_id2)

        self.set_lhs_range(new_lhs_range)
        self.set_rhs_range(new_rhs_range)

        if self.auto_fit_residual:
            self.fit_residual()

    def move_left(self, num_frags = 1):
        """Move the partition point num_frags to the left.
        """
        self.move(-num_frags)

    def move_right(self, num_frags = 1):
        """Move the partition point num_frgs to the right.
        """
        self.move(num_frags)

    def fit_residual(self):
        chain = self.cpartition.chain
        self.lhs_tls.fit_residual(chain)
        self.rhs_tls.fit_residual(chain)


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
    console.stdoutln(partitions)


## testing
def testmain():
    import tlsmd_analysis

    struct = FileIO.LoadStructure(
        file = "/home/tlsmd/public_html/examples/1KP8/ANALYSIS/1KP8.pdb")
    chain = tlsmd_analysis.ConstructSegmentForAnalysis(struct.get_chain("A"))

    cpartition = opt_containers.ChainPartition(chain, 3)
    groups = [
        [("2", "135"), ("411", "525")],
        [("136", "190"), ("375", "410")],
        [("191", "374")]]

    for segment_ranges in groups:
        tls = opt_containers.TLSSegment(segment_ranges)
        cpartition.add_tls_segment(tls)

    cpartition.fit_residual()

    partitions = ChainPartitionList(cpartition)
    print partitions,  cpartition.residual()

    for part in partitions:
        print "PARTITION ",part
        part.move_left(10)
        for x in range(20):
            try:
                part.move_right()
            except ValueError:
                print "reached limit"
            print str(partitions).replace(" ",""), cpartition.residual()
        part.move_left(10)

    return

    for x in range(60):
        try:
            partitions[2].move_right()
        except ValueError:
            print "reached limit"
        print str(partitions).replace(" ",""), cpartition.residual()

    return

    print "Move 1 right"
    for x in range(15):
        try:
            partitions[0].move_right()
        except ValueError:
            print "reached limit"
        print partitions, cpartition.residual()


if __name__ == "__main__":
    testmain()
