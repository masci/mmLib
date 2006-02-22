## TLS Motion Determination (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import sys
import math
import numpy

try:
    from Bio import pairwise2
except ImportError:
    print "You need to install BioPython to use this feature"
    sys.exit(-1)

from mmLib import Constants, AtomMath, Structure, Superposition, TLS

SUPER_ATOMS  = ["N","CA","C"]


def calc_angle(a, b):
    cos_ab = numpy.dot(a, b)/ (AtomMath.length(a) * AtomMath.length(b))
    return Constants.RAD2DEG * abs(math.acos(cos_ab))


def calc_directional_overlap(a, b):
    cos_ab = numpy.dot(a, b) / (AtomMath.length(a) * AtomMath.length(b))
    if cos_ab < 0.0:
        return 0.0
    return abs(cos_ab)


def align_chains(chain1, chain2):
    """Adds a .equiv attribute to each fragment of the chain
    referencing the equivalent fragment in the other chain.
    """
    print "Chain Alignment"
    seq1 = chain1.calc_sequence_one_letter_code()
    print "Length of Chain 1: %d" % (len(seq1))
    seq2 = chain2.calc_sequence_one_letter_code()
    print "Length of Chain 2: %d" % (len(seq2))
    align = pairwise2.align.globalxs(seq1, seq2, -0.5, -0.125)
    print align[0]

    seq1_align = align[0][0]
    seq2_align = align[0][1]

    iter1 = chain1.iter_fragments()
    iter2 = chain2.iter_fragments()

    chain1_equiv = {}
    chain2_equiv = {}
    
    for i in range(len(seq1_align)):
        frag1 = None
        frag2 = None

        if seq1_align[i] != '-':
            frag1 = iter1.next()

        if seq2_align[i] != '-':
            frag2 = iter2.next()

        if frag1 and frag2:
            chain1_equiv[frag1] = frag2
            chain2_equiv[frag2] = frag1

    return align[0], chain1_equiv, chain2_equiv


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
    #return Superposition.SuperimposeAtoms(alist)

        
class TLSConformationPredctionHypothosis(object):
    """Calculate the directional overlap of the TLS displacements in
    struct_file/tls_file with the displacements in target_struct_file.
    """    
    def __init__(self, chain, target_chain):
        self.chain = chain
        self.target_chain = target_chain
        self.alignment_score = self.align_source_target_chains()

    def add_conformation_prediction_to_chain_partition(self, cpartition):
        for tls in cpartition.iter_tls_segments():
            self.calc_superposition(tls)

    def align_source_target_chains(self):
        """Performs a sequence alginment folled by a structure alignment
        of the target chain to the source chain.  The coordinates of the
        target chain is altered.
        """
        alignment_score, chain1_equiv, chain2_equiv, = align_chains(self.chain, self.target_chain)
        self.srctgt_equiv = chain1_equiv

        sresult = SuperimposeChains(self.target_chain, self.chain, chain2_equiv, ["CA"])
        self.chain.target_chain_sresult = sresult

        print "Structure Superposition RMSD: %6.2f" % (sresult.rmsd)

        for atm in self.target_chain.iter_all_atoms():
            atm.align_position = sresult.transform(atm.position)

        ## residue type mismatches in the sequence alignment of the fragments
        for frag1 in self.chain.iter_fragments():
            try:
                frag2 =  self.srctgt_equiv[frag1]
            except KeyError:
                continue
            if frag1.res_name != frag2.res_name:
                print "EEK! %s::%s != %s::%s" % (
                    frag1.fragment_id, frag1.res_name,
                    frag2.fragment_id, frag2.res_name)

        return alignment_score

    def calc_superposition(self, tls):
        plist = []
        msd = 0.0
        segment = tls.segment

        for frag1 in segment.iter_fragments():
            try:
                frag2 = self.srctgt_equiv[frag1]
            except KeyError:
                continue
            for name in SUPER_ATOMS:
                atm1 = frag1.get_atom(name)
                atm2 = frag2.get_atom(name)
                if atm1 == None or atm2 == None:
                    print "EEK! No Equivalent Atom ", atm1
                    continue
                plist.append((atm1.position, atm2.align_position))
                d = atm1.position - atm2.align_position
                msd += numpy.dot(d,d)
                    
        rmsd_pre_alignment = math.sqrt(msd / len(plist))
        tls.rmsd_pre_alignment = rmsd_pre_alignment

        sresult = Superposition.SuperimposePositions(plist)
        tls.sresult = sresult

        rotation = math.degrees(2.0 * math.acos(sresult.Q[0]))
        if rotation > 180.0:
            rotation = 360.0 - rotation

        fragstr = "%s:%s-%s" % (self.chain.chain_id, tls.frag_id1, tls.frag_id2)
        print "TLS Group::%20s  Num Atoms::%4d  RMSD PRE ALIGN::%6.2f  RMSD::%6.2f  TRANSORM ROTATION::%6.2f" % (
            fragstr, len(plist), rmsd_pre_alignment, sresult.rmsd, rotation)

        ## screw displacement vector
        vscrew = AtomMath.normalize(numpy.array([sresult.Q[1],sresult.Q[2],sresult.Q[3]], float))
        print "superposition rotation vector: ",vscrew
        tls.superposition_vscrew = vscrew * rotation

        ## fit the isotropic TLS model to the group
        evals, evecs = numpy.linalg.eigenvectors(tls.tls_group.itls_L)

        for i in range(3):
            eval = evals[i]
            evec = evecs[i]
            
            lname = "L%d_eigen_val" % (i)

            if (eval * Constants.RAD2DEG2) < 1.0:
                continue

            ang = min(calc_angle(evec, vscrew), calc_angle(-evec, vscrew))

            print "%s  magnitude::%6.2f  vector angle::%6.2f" % (
                lname, eval*Constants.RAD2DEG2, ang)
            

