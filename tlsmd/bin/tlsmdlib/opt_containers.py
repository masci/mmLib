## TLS Motion Determination (TLSMD)
## Copyright 2002-2006 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import math
import copy
import numpy

from mmLib import Constants, Structure, TLS

import conf
import tls_calcs


class TLSSegment(object):
    """Information on a TLS rigid body segment of a protein chain.
    """
    def __init__(self, segment_ranges = [], **args):
        self.segment_ranges = segment_ranges
        self.chain_id = args.get("chain_id")
        self.method = args.get("method")
        self.__residual = args.get("residual")
        self.__num_atoms = args.get("num_atoms")
        self.__num_residues = args.get("num_residues")
        
        ## added by a call to the method fit_to_chain()
        self.tls_group = None
        self.tls_info = None
        self.itls_info = None
        self.model_tls_info = None
        self.model = None
        self.rmsd_b = None
        self.segments = []

        ## added by HTML generation code
        self.color = None

        ## added by structure comparison
        self.rmsd_pre_alignment = None
        self.sresult = None
        self.superposition_vscrew = None

    def __str__(self):
        return "%s:%s" % (self.chain_id, self.display_label())

    def copy(self):
        return copy.copy(self)

    def display_label(self):
        l = []
        for frag_id1, frag_id2 in self.segment_ranges:
            l.append("%s-%s" % (frag_id1, frag_id2))
        return "; ".join(l)

    def filename_label(self):
        l = []
        for frag_id1, frag_id2 in self.segment_ranges:
            l.append("%s%s_%s" % (self.chain_id, frag_id1, frag_id2))
        return "_".join(l)

    def jmol_select(self):
        l = []
        for frag_id1, frag_id2 in self.segment_ranges:
            l.append("%s-%s:%s" % (frag_id1, frag_id2, self.chain_id))
        return ",".join(l)

    def set_color(self, color):
        assert self.color is None
        self.color = color

    def num_segment_ranges(self):
        return len(self.segment_ranges)

    def iter_segment_ranges(self):
        return iter(self.segment_ranges)

    def iter_segments(self):
        return iter(self.segments)

    def num_residues(self):
        return self.__num_residues

    def num_residues_segment(self):
        n = 0
        for segment in self.segments:
            n += len(segment)
        return n

    def iter_fragments(self):
        for segment in self.segments:
            for frag in segment.iter_fragments():
                yield frag

    def num_atoms(self):
        return self.__num_atoms

    def num_atoms_segment(self):
        n = 0
        for segment in self.segments:
            n += segment.count_all_atoms()
        return n

    def iter_atoms(self):
        for segment in self.segments:
            for atm in segment.iter_all_atoms():
                yield atm

    def mean_b(self):
        return self.tls_info["exp_mean_temp_factor"]

    def mean_anisotropy(self):
        return self.tls_info["exp_mean_anisotropy"]

    def tls_mean_b(self):
        return self.tls_info["tls_mean_temp_factor"]

    def tls_mean_anisotropy(self):
        return self.tls_info["tls_mean_anisotropy"]

    def chi2(self):
        return self.__residual / self.__num_residues
    
    def residual(self):
        return self.__residual

    def fit_residual(self, chain):
        tlsdict = chain.tls_analyizer.isotropic_fit(self.segment_ranges)
        if tlsdict:
            self.__residual = tlsdict["residual"]
            self.__num_atoms = tlsdict["num_atoms"]
            return self.__residual
        return None

    def fit_to_chain(self, chain):
        """Re-sets all derived information in the TLSSegment.
        """
        ## cut segment from chain using segment ranges
        segments = []
        for frag_id1, frag_id2 in self.segment_ranges:
            segments.append(chain[frag_id1:frag_id2])
        self.segments = segments

        ## put all atoms in the segment into a new TLSGroup instance
        tls_group = TLS.TLSGroup()
        for segment in self.segments:
            for atm in segment.iter_all_atoms():
                if atm.include is True:
                    tls_group.append(atm)
        self.tls_group = tls_group

        if len(self.tls_group) != self.num_atoms():
            print "fit_to_chain: EEK! (%s) len(self.tls_group)=%d != self.num_atoms()=%d" % (
                self, len(self.tls_group), self.num_atoms())
            raise SystemExit

        ## fit the TLS group parameters
        self.fit_tls_parameters(chain)
        
        ## helpful additions
        tls_info  = self.tls_group.calc_tls_info()
        itls_info = TLS.calc_itls_center_of_reaction(
            self.tls_group.itls_T,
            self.tls_group.itls_L,
            self.tls_group.itls_S,
            self.tls_group.origin)

        self.tls_info = tls_info
        self.itls_info = itls_info

        if conf.globalconf.tls_model in ["ISOT", "NLISOT"]:
            self.tls_group.model = "ISOT"
            self.model_tls_info = itls_info
        elif conf.globalconf.tls_model in ["ANISO", "NLANISO"]:
            self.tls_group.model = "ANISO"
            self.model_tls_info = tls_info

        self.rmsd_b = tls_calcs.calc_rmsd_tls_biso(self.tls_group)

    def fit_tls_parameters(self, chain):
        """Use the non-linear TLS model to calculate tensor values.
        """
        tls_group = self.tls_group
        
        ## anisotropic model
        tlsdict = chain.tls_analyzer.constrained_anisotropic_fit(self.segment_ranges)
        T, L, S, origin = tls_calcs.tlsdict2tensors(tlsdict)
        tls_group.T = T
        tls_group.L = L
        tls_group.S = S
        tls_group.origin = origin

        ## isotropic model
        itlsdict = chain.tls_analyzer.constrained_isotropic_fit(self.segment_ranges)
        IT, IL, IS, IOrigin = tls_calcs.isotlsdict2tensors(itlsdict)
        tls_group.itls_T = IT
        tls_group.itls_L = IL
        tls_group.itls_S = IS

        assert numpy.allclose(tls_group.origin, IOrigin)


class ChainPartition(object):
    """Collection of TLSSegment objects describing one multi-TLS
    group partitioning of a protein chain.
    """
    def __init__(self, chain, ntls):
        self.chain = chain
        self.__ntls = ntls
        self.tls_list = []

    def __str__(self):
        return "".join(["(%s)" % (str(x)) for x in self.tls_list])

    def first_frag_id(self):
        return self.chain[0].fragment_id

    def last_frag_id(self):
        return self.chain[-1].fragment_id

    def add_tls_segment(self, tls):
        """The argument tls is a dictionary containing a bunch of great information
        about the TLS segment which is part of the optimal partitioning.
        """
        assert isinstance(tls, TLSSegment)
        self.tls_list.append(tls)

    def residual(self):
        """Calculate the total residual from the component TLS groups.
        """
        residual = 0.0
        for tls in self.iter_tls_segments():
            residual += tls.residual()
        return residual

    def rmsd_b(self):
        """Calculate the RMSD of the multi-TLS group prediction for the
        entire chain.
        """
        num_residues = 0
        msd_sum = 0.0
        for tls in self.iter_tls_segments():
            num_residues += tls.num_residues()
            msd_sum += tls.residual()
        return Constants.U2B * math.sqrt(msd_sum / num_residues)

    def chi2(self):
        chi2 = 0.0
        for tls in self.iter_tls_segments():
            chi2 += tls.chi2()
        return chi2

    def is_valid(self):
        """Return True if the optimization is valid; otherwise, return False.
        """
        return len(self.tls_list) > 0

    def num_tls_segments(self):
        return len(self.tls_list)

    def iter_tls_segments(self):
        return iter(self.tls_list)



class ChainPartitionCollection(object):
    """Contains all the ChainPartition objects for a chain.
    """
    def __init__(self, chain):
        self.chain = chain
        self.chain_id = chain.chain_id
        self.ntls_chain_partition_list = []

    def insert_chain_partition(self, cpartition):
        """Inserts a new ChainPartitin instance into the collection.
        Removes any ChainPartition currently in the collection with
        the same number of TLS groups.
        """
        assert isinstance(cpartition, ChainPartition)
        for i, (ntls, cpart) in enumerate(self.ntls_chain_partition_list):
            if ntls == cpartition.num_tls_segments():
                if cpart == cpartition:
                    return
                del self.ntls_chain_partition_list[i]
                break
        self.ntls_chain_partition_list.append((cpartition.num_tls_segments(), cpartition))
        self.ntls_chain_partition_list.sort()

    def num_chain_partitions(self):
        return len(self.ntls_chain_partition_list)

    def iter_ntls_chain_partitions(self):
        return iter(self.ntls_chain_partition_list)

    def iter_chain_partitions(self):
        for ntls, cpartition in self.iter_ntls_chain_partitions():
            yield cpartition

    def iter_ntls(self):
        for ntls, cpartition in self.iter_ntls_chain_partitions():
            yield ntls

    def get_chain_partition(self, find_ntls):
        for ntls, cpartition in self.ntls_chain_partition_list:
            if ntls == find_ntls:
                return cpartition
        return None

    def max_ntls(self):
        return max(self.iter_ntls())

    def min_ntls(self):
        return min(self.iter_ntls())
    
