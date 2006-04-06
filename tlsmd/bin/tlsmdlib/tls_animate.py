## TLS Motion Determination (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import copy
import string
import math
import numpy

from mmLib import Structure, FileIO, Gaussian, AtomMath, TLS

import conf

TWO_PI = 2.0 * math.pi


def iter_fragment_range(chain, frag_id1, frag_id2):
    for frag in chain.iter_fragments():
        if Structure.fragment_id_lt(frag.fragment_id, frag_id1):
            continue
        if Structure.fragment_id_gt(frag.fragment_id, frag_id2):
            break
        yield frag


class TLSAnimateFailure(Exception):
    pass


class TLSAnimate(object):
    """Create a multi-model PDB file which each model a frame of a TLS
    animation.
    """
    
    def __init__(self, chain, cpartition):
        ## copy and get the anisotropic ADPs out of the structure
        self.struct = self.copy_struct(chain.struct, chain)
        
        self.chain = chain
        self.cpartition = cpartition

        self.L1_chain = None
        self.L2_chain = None
        self.L3_chain = None

        self.construct_chain_copies()

    def copy_struct(self, struct, chain):
        """Make a copy of the argument structure to use for generating
        the animation.  Only copy the chain specified in chain_id.
        """
        cp_struct = Structure.Structure(structure_id = struct.structure_id)
        chain_id = chain.chain_id

        for chain in struct.iter_chains():

            if chain.chain_id == chain_id or chain.count_fragments() < 200:
                include_chain = True
	    else:
		include_chain = False
		
            for frag in chain.iter_fragments():

                ## skip all waters
                if frag.is_water() == True:
                    continue

                elif frag.is_standard_residue() == True and include_chain:
                    for atm in frag.iter_atoms():

                        if atm.name != "CA":
                            continue
                        
                        cp_atom = Structure.Atom(
                            chain_id    = atm.chain_id,
                            fragment_id = atm.fragment_id,
                            res_name    = atm.res_name,
                            element     = atm.element,
                            name        = atm.name,
                            temp_factor = atm.temp_factor,
                            occupancy   = atm.occupancy,
                            position    = atm.position.copy())
                        cp_struct.add_atom(cp_atom, True)

                elif frag.is_standard_residue() == False:
                    for atm in frag.iter_atoms():
                        cp_atom = Structure.Atom(
                            chain_id    = atm.chain_id,
                            fragment_id = atm.fragment_id,
                            res_name    = atm.res_name,
                            element     = atm.element,
                            name        = atm.name,
                            temp_factor = atm.temp_factor,
                            occupancy   = atm.occupancy,
                            position    = atm.position.copy())
                        cp_struct.add_atom(cp_atom, True)

        cp_struct.sort()
        return cp_struct
                
    def construct_animation(self, filename):
        """Save the animated structure to the given filename.
        """
        ##self.phase_assignment()
        for phase in (0.5, 1.0, 0.5, 0.0, -0.5, -1.0, -0.5):
            self.construct_frame(phase)
        FileIO.SaveStructure(struct = self.struct, fil = filename)

    def next_model_id(self):
        """Return the next availible model_id in self.struct
        """
        model_id = 0
        while True:
            model_id += 1            
            model = self.struct.get_model(model_id)
            if model == None:
                return model_id
            
    def next_chain_id(self):
        """Return the next availible chain_id in self.astruct
        """
        for chain_id in string.uppercase:
            chain = self.struct.get_chain(chain_id)
            if chain == None:
                return chain_id
        raise TLSAnimateFailure()

    def construct_chain_copies(self):
        """Create two other copies 
        """
        chain_id = self.chain.chain_id        
        self.L1_chain = self.struct.get_chain(chain_id)

        self.L2_chain = copy.deepcopy(self.L1_chain)
        self.L2_chain.set_chain_id(self.next_chain_id())
        self.struct.add_chain(self.L2_chain, True)

        self.L3_chain = copy.deepcopy(self.L1_chain)
        self.L3_chain.set_chain_id(self.next_chain_id())
        self.struct.add_chain(self.L3_chain, True)
        
    def construct_frame(self, phase):
        """Create a new model in self.struct with the TLS displacements
        caused from the three screw dispacement axes displaced by the
        sin(phase).
        """
        ## copy the original model and add it to the structure
        model1 = self.struct.get_model(1)
        
        model = copy.deepcopy(model1)
        model.set_model_id(self.next_model_id())
        self.struct.add_model(model, True)

        ## now displace the new model according to the tls group
        ## screw axes and given phase
        for tls in self.cpartition.iter_tls_segments():
            self.displace_model(model, tls, phase)

    def displace_model(self, model, tls, phase):
        """Displace the given model by the tls. 
        """
        tls_info  = tls.model_tls_info
        cor       = tls_info["COR"]

        for n, Lx_rmsd, Lx_vec, Lx_rho, Lx_pitch in [
            (1, "L1_rmsd", "L1_eigen_vec", "L1_rho", "L1_pitch"),
            (2, "L2_rmsd", "L2_eigen_vec", "L2_rho", "L2_pitch"),
            (3, "L3_rmsd", "L3_eigen_vec", "L3_rho", "L3_pitch") ]:

            Lrmsd  = tls_info[Lx_rmsd]
            Lvec   = tls_info[Lx_vec]
            Lrho   = tls_info[Lx_rho]
            Lpitch = tls_info[Lx_pitch]

            ## pre-calculations for screw displacement
            Lorigin = cor + Lrho
            Lrot = Gaussian.GAUSS3C[conf.ADP_PROB] * Lrmsd * phase
            D = AtomMath.dmatrixu(Lvec, Lrot)
            d_screw = (Lrot * Lpitch) * Lvec
            
            if n==1:
                chain_id = self.L1_chain.chain_id
            elif n==2:
                chain_id = self.L2_chain.chain_id
            elif n==3:
                chain_id = self.L3_chain.chain_id

            chain = model.get_chain(chain_id)

            for frag_id1, frag_id2 in tls.iter_segment_ranges():
                for frag in iter_fragment_range(chain, frag_id1, frag_id2):
                    for atm in frag.iter_atoms():
                        d = numpy.matrixmultiply(D, atm.position - Lorigin) + d_screw
                        atm.position += d
