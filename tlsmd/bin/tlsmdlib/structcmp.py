## TLS Motion Determination (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import sys
import math
import numpy

from Bio import pairwise2

from mmLib import Constants, AtomMath, Structure, FileLoader, Superposition
from mmLib.Extensions import TLS

SUPER_ATOMS  = ["N","CA","C"]


def calc_angle(a, b):
    cos_ab = numpy.dot(a, b)/ (AtomMath.length(a) * AtomMath.length(b))
    return Constants.RAD2DEG * abs(math.acos(cos_ab))


def calc_directional_overlap(a, b):
    cos_ab = numpy.dot(a, b) / (AtomMath.length(a) * AtomMath.length(b))
    if cos_ab < 0.0:
        return 0.0
    return abs(cos_ab)


def align_chains(src_chn, dst_chn):
    """Adds a .equiv attribute to each fragment of the chain
    referencing the equivalent fragment in the other chain.
    """
    print "Chain Alignment"
    srcseq = src_chn.calc_sequence_one_letter_code()
    print "Length of Source Chain: %d" % (len(srcseq))
    dstseq = dst_chn.calc_sequence_one_letter_code()
    print "Length of Destination Chain: %d" % (len(dstseq))
    align = pairwise2.align.globalxs(srcseq, dstseq, -0.25, -0.125)
    print align[0]

    srcseq_align = align[0][0]
    dstseq_align = align[0][1]

    srciter = src_chn.iter_fragments()
    dstiter = dst_chn.iter_fragments()

    srcdst_equiv = {}
    dstsrc_equiv = {}
    
    for i in range(len(srcseq_align)):
        src_frag = None
        dst_frag = None

        if srcseq_align[i] != '-':
            src_frag = srciter.next()

        if dstseq_align[i] != '-':
            dst_frag = dstiter.next()

        if src_frag and dst_frag:
            srcdst_equiv[src_frag] = dst_frag
            dstsrc_equiv[dst_frag] = src_frag

    return srcdst_equiv, dstsrc_equiv


def SuperimposeChains(source_chain, target_chain, srcdst_equiv, atom_names = ["CA"]):
    """Superimpose two homologus protein chains.
    """
    alist = []

    for frag1 in source_chain.iter_fragments():
        try:
            frag2 = srcdst_equiv[frag1]
        except KeyError:
            continue

        for name in atom_names:
            atm1 = frag1.get_atom(name)
            atm2 = frag2.get_atom(name)
            if atm1 == None or atm2 == None:
                continue
            alist.append((atm1,atm2))

    return Superposition.SuperimposeAtomsOutlierRejection(alist, 0.5)

        
class TLSConformationPredctionHypothosis(object):
    """Calculate the directional overlap of the TLS displacements in
    struct_file/tls_file with the displacements in target_struct_file.
    """    
    def __init__(self, chain, cpartition, target_chain):
        self.chain = chain
        self.cpartition = cpartition
        self.target_chain = target_chain

    def test_hypososis(self):
        self.push_target_struct_atom_positions()
        
        ## align the overall structures
        self.align_source_target_chains()

        ## calculate CA superposition of the segments
        for tls in self.cpartition.iter_tls_segments():
            self.calc_superposition(tls)

        self.pop_target_atom_positions()

    def push_target_struct_atom_positions(self):
        for atm in self.target_struct.iter_all_atoms():
            atm.orig_position = atm.position

    def pop_target_struct_atom_positions(self):
        for atm in self.target_struct.iter_all_atoms():
            atm.position = atm.orig_position
            del atm.orig_position

    def align_source_target_chains(self):
        """Performs a sequence alginment folled by a structure alignment
        of the target chain to the source chain.  The coordinates of the
        target chain is altered.
        """
        schn = self.chain
        tchn = self.target_chain

        srcdst_equiv, dstsrc_equiv = align_chains(tchn, schn)
        self.srctgt_equiv = dstsrc_equiv

        sresult = SuperimposeChains(tchn, schn, dstsrc_equiv, SUPER_ATOMS)
        print "Structure Superposition RMSD: %6.2f" % (sresult.rmsd)

        for atm in self.target_struct.iter_all_atoms():
            pos =  numpy.matrixmultiply(sresult.R, atm.position - sresult.src_origin)
            atm.align_position = pos + sresult.dst_origin

    def calc_tls_segment_superposition(self, tls):
        al = []
        msd = 0.0

        try:
            segment = tls["segment"]
        except KeyError:
            continue

        for frag1 in segment.iter_fragments():
            if frag1.equiv == None:
                continue
            
            frag2 = frag1.equiv

            for name in SUPER_ATOMS:
                atm1 = frag1.get_atom(name)
                atm2 = frag2.get_atom(name)
                if atm1 == None or atm2 == None:
                    continue

                try:
                    assert atm1.res_name == atm2.res_name
                except AssertionError:
                    print "EEK! %s::%s != %s::%s" % (atm1.fragment_id, atm1.res_name,
                                                     atm2.fragment_id, atm2.res_name)
                al.append((atm1,atm2))
                d = atm1.position - atm2.position
                msd += numpy.dot(d,d)
                    
        rmsd_pre_alignment = math.sqrt(msd / len(al))
        tls_group.super = Superposition.SuperpositionAtoms(al)

        print "TLS Group::%20s  Num Atoms::%4d  RMSD PRE ALIGN::%6.2f  RMSD::%6.2f" % (
            tls_group.name, len(al), rmsd_pre_alignment, tls_group.super.rmsd)

        ## screw displacement vector
        Q = tls_group.super.Q
        vscrew = AtomMath.normalize(numpy.array([Q[1],Q[2],Q[3]], float))

        ## fit the isotropic TLS model to the group
        fit_itls_group(tls_group)
        evals, evecs = numpy.linalg.eigenvectors(tls_group.itls_L)

        for i in range(3):
            eval = evals[i]
            evec = evecs[i]
            
            lname = "L%d_eigen_val" % (i)

            if (eval * Constants.RAD2DEG2) < 1.0: continue

            ang = min(calc_angle(evec, vscrew), calc_angle(-evec, vscrew))

            print "%s  magnitude::%6.2f  vector angle::%6.2f" % (
                lname, eval*Constants.RAD2DEG2, ang)
            

