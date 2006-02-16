
import sys
import getopt
import copy
import string
import math
import numpy

from Bio import pairwise2

from mmLib import Constants, AtomMath, Structure, FileLoader, Superposition
from mmLib.Extensions import TLS

import nonlineartls


SUPER_ATOMS  = ["N","CA","C"]


def iter_atoms_to_dict(atom_iter):
    """Converts the Atoms of a Chain/Segment to a list of dictionaries 
    for transfer over xmlrpc.  Only the information required to fit
    TLS groups to the Atoms is set in the Atom dictionaries to
    reduce traffic over the xmlrpc calls.
    """
    for atm in atom_iter:
        if atm.name not in SUPER_ATOMS: continue
        
        atm_desc = {}

        atm_desc["name"]    = atm.name
        atm_desc["frag_id"] = atm.fragment_id

        atm_desc["x"] = atm.position[0]
        atm_desc["y"] = atm.position[1]
        atm_desc["z"] = atm.position[2]

        atm_desc["u_iso"] = Constants.B2U * atm.temp_factor

        U = atm.get_U()
        atm_desc["u11"] = U[0,0]
        atm_desc["u22"] = U[1,1]
        atm_desc["u33"] = U[2,2]
        atm_desc["u12"] = U[0,1]
        atm_desc["u13"] = U[0,2]
        atm_desc["u23"] = U[1,2]

        ## calculate weight
        atm_desc["weight"] = 1.0

        atm_desc["ifrag"] = 1

        yield atm_desc



def fit_itls_group(tls_group):
    alist = list(iter_atoms_to_dict(iter(tls_group)))
    nltls = nonlineartls.NLTLSModel()
    nltls.set_xmlrpc_chain(alist)
    itls = nltls.isotropic_fit_segment(0, len(alist)-1)

    tls_group.origin = numpy.array(
        [itls["x"], itls["y"], itls["z"]], float)
    
    tls_group.itls_T = itls["it"]
    tls_group.itls_L = numpy.array(
        [ [itls["il11"], itls["il12"], itls["il13"]],
          [itls["il12"], itls["il22"], itls["il23"]],
          [itls["il13"], itls["il23"], itls["il33"]] ], float)

    tls_group.itls_S = numpy.array(
        [itls["is1"], itls["is2"], itls["is3"]], float)


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
    
    for i in range(len(srcseq_align)):
        src_frag = None
        dst_frag = None

        if srcseq_align[i] != '-':
            src_frag = srciter.next()

        if dstseq_align[i] != '-':
            dst_frag = dstiter.next()

        if src_frag and dst_frag:
            src_frag.equiv = dst_frag
            dst_frag.equiv = src_frag

        if src_frag and not dst_frag:
            src_frag.equiv = None

        if dst_frag and not src_frag:
            dst_frag.equiv = None


def SuperimposeChains(source_chain, target_chain, atom_names = ["CA"]):
    """Superimpose two homologus protein chains.
    """
    alist = []

    for frag1 in source_chain.iter_fragments():
        if frag1.equiv == None: continue
        frag2 = frag1.equiv

        for name in atom_names:
            atm1 = frag1.get_atom(name)
            atm2 = frag2.get_atom(name)
            if atm1 == None or atm2 == None: continue
            alist.append((atm1,atm2))

    return Superposition.SuperimposeAtomsOutlierRejection(alist, 0.5)

        
class TLSConformationPredctionHypothosis(object):
    """Calculate the directional overlap of the TLS displacements in
    struct_file/tls_file with the displacements in target_struct_file.
    """    
    def __init__(self):
        self.source_chain = None
        self.source_tls_group_list = None
        self.target_chain = None

        ## align the overall structures
        self.align_structures()

        ## load the TLS description
        tls_file = TLS.TLSFile()
        tls_file.set_file_format(TLS.TLSFileFormatTLSOUT())
        tls_file.load(open(self.tls_file, "r"))

        self.tls_group_list = []

        for tls_desc in tls_file.tls_desc_list:
            tls_group = tls_desc.construct_tls_group_with_atoms(self.struct)
            self.tls_group_list.append(tls_group)
            
            tls_group.tls_desc = tls_desc
            tls_group.tls_info = tls_group.calc_tls_info()

        ## calculate CA superposition of the segments
        for tls_group in self.tls_group_list:
            self.calc_superposition(tls_group)

    def load_source_chain(self, struct_path, chain_id):
        struct = FileLoader.LoadStructure(fil = struct_path)
        self.set_source_chain(struct.get_chain(chain_id))
    
    def load_target_chain(self, struct_path, chain_id):
        struct = FileLoader.LoadStructure(fil = struct_path)
        self.set_target_chain(struct.get_chain(chain_id))

    def load_tls_file(self, tls_file_path):
        """Load the TLS groups from a given file.
        """
        tls_file = TLS.TLSFile()
        tls_file.set_file_format(TLS.TLSFileFormatTLSOUT())
        tls_file.load(open(self.tls_file, "r"))

        tls_group_list = []

        for tls_desc in tls_file.tls_desc_list:
            tls_group = tls_desc.construct_tls_group_with_atoms(self.struct)
            tls_group_list.append(tls_group)
            
            tls_group.tls_desc = tls_desc
            tls_group.tls_info = tls_group.calc_tls_info()

        self.set_tls_group_list(tls_group_list)

    def set_source_chain(self, chain):
        self.source_chain = chain

    def set_target_chain(self, chain):
        self.target_chain = chain

    def set_tls_group_list(self, tls_group_list):
        self.tls_group_list = tls_group_list



    def test_hypososis(self):
        ## align the overall structures
        self.align_source_target_chains()

        ## calculate CA superposition of the segments
        for tls_group in self.tls_group_list:
            self.calc_superposition(tls_group)

    def align_source_target_chains(self):
        """Performs a sequence alginment folled by a structure alignment
        of the target chain to the source chain.  The coordinates of the
        target chain is altered.
        """
        src_chn = self.source_chain
        trg_chn = self.target_chain

        align_chains(src_chn, trg_chn)

        sresult = SuperimposeChains(trg_chn, src_chn, SUPER_ATOMS)
        print "Structure Superposition RMSD: %6.2f" % (sresult.rmsd)

        for atm in self.target_struct.iter_all_atoms():
            pos =  numpy.matrixmultiply(sresult.R, atm.position - sresult.src_origin)
            atm.align_position = pos + sresult.dst_origin

    def calc_superposition(self, tls_group):
        al = []
        msd = 0.0
        
        for (chain_id1, frag_id1, chain_id2, frag_id2, sel) in tls_group.tls_desc.range_list:
            chain1 = self.struct.get_chain(chain_id1)
            if chain1 == None:
                continue

            try:
                segment = chain1[frag_id1:frag_id2]
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
            

def main():
    schain = SOURCE_CHAIN
    tchain = TARGET_CHAIN
    
    try:
        (opts, args) = getopt.getopt(sys.argv[1:], "s:t:")
    except getopt.GetoptError:
        sys.exit(1)

    for key, val in opts:
        if key == "-s": schian = val
        if key == "-t": tchain = val

    scalc = TLS_Dot(args[0], args[1], args[2], schain, tchain)


if __name__=="__main__":
    main()
