## TLS Minimized Domains (TLSMD)
## Copyright 2002 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import copy
import string
import math

from mmLib.Structure      import *
from mmLib.FileLoader     import *
from mmLib.Extensions.TLS import *

TWO_PI = 2.0 * math.pi
ADP_PROB = 85


class ScrewDisplacer(object):
    def __init__(self, cor, Lval, Lvec, Lrho, Lpitch, phase):
        self.Lorigin = cor + Lrho

        Lrot = GAUSS3C[ADP_PROB]*math.sqrt(Lval)*math.sin(TWO_PI*phase)
        self.D = dmatrixu(Lvec, Lrot)

        self.d_screw = (Lrot * Lpitch) * Lvec

    def calc_d(self, position):
        d_rotation = matrixmultiply(self.D, position - self.Lorigin)
        return d_rotation + self.d_screw
        

def iter_fragment_range(chain, frag_id1, frag_id2):
    for frag in chain.iter_fragments():
        if fragment_id_lt(frag.fragment_id, frag_id1):
            continue
        if fragment_id_gt(frag.fragment_id, frag_id2):
            break
        yield frag


class TLSAnimateFailure(Exception):
    pass


class TLSAnimateAltLoc(object):
    """Create a multi-model PDB file which each model a frame of a TLS
    animation.
    """
    
    def __init__(self, struct, chainopt, tlsopt):

        ## copy and get the anisotropic ADPs out of the structure
        self.struct = self.copy_struct(struct, chainopt["chain_id"])
        
        self.chainopt = chainopt
        self.tlsopt   = tlsopt

        self.L1_chain = None
        self.L2_chain = None
        self.L3_chain = None

        #self.construct_chain_copies()

    def copy_struct(self, struct, chain_id):
        """Make a copy of the argument structure to use for generating
        the animation.  Only copy the chain specified in chain_id.
        """
        cp_struct = Structure(structure_id = struct.structure_id)

        for chain in struct.iter_chains():
            if chain.chain_id==chain_id:
                for atm in chain.iter_atoms():
                    cp_atm = Atom(
                        chain_id    = atm.chain_id,
                        fragment_id = atm.fragment_id,
                        res_name    = atm.res_name,
                        element     = atm.element,
                        name        = atm.name,
                        alt_loc     = "A",
                        position    = atm.position.copy())
                    cp_struct.add_atom(cp_atm, True)
                    cp_atm = Atom(
                        chain_id    = atm.chain_id,
                        fragment_id = atm.fragment_id,
                        res_name    = atm.res_name,
                        element     = atm.element,
                        name        = atm.name,
                        alt_loc     = "B",
                        position    = atm.position.copy())
                    cp_struct.add_atom(cp_atm, True)
                    cp_atm = Atom(
                        chain_id    = atm.chain_id,
                        fragment_id = atm.fragment_id,
                        res_name    = atm.res_name,
                        element     = atm.element,
                        name        = atm.name,
                        alt_loc     = "C",
                        position    = atm.position.copy())
                    cp_struct.add_atom(cp_atm, True)
            else:
                for atm in chain.iter_atoms():
                    cp_atom = Atom(
                        chain_id    = atm.chain_id,
                        fragment_id = atm.fragment_id,
                        res_name    = atm.res_name,
                        element     = atm.element,
                        name        = atm.name,
                        position    = atm.position.copy())
                    cp_struct.add_atom(cp_atom, True)


        cp_struct.sort()
        return cp_struct
                
    def construct_animation(self, filename):
        """Save the animated structure to the given filename.
        """
        for phase in (0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875):
            self.construct_frame(phase)
        SaveStructure(struct=self.struct, fil=filename)

    def next_model_id(self):
        """Return the next availible model_id in self.struct
        """
        model_id = 0
        while True:
            model_id += 1            
            model = self.struct.get_model(model_id)
            if model==None:
                return model_id
            
    def next_chain_id(self):
        """Return the next availible chain_id in self.astruct
        """
        for chain_id in string.uppercase:
            chain = self.struct.get_chain(chain_id)
            if chain==None:
                return chain_id
        raise TLSAnimateFailure()

    def construct_chain_copies(self):
        """Create two other copies 
        """
        chain_id = self.chainopt["chain_id"]        
        self.L1_chain = self.struct.get_chain(chain_id)

        self.L2_chain = copy.deepcopy(self.L1_chain)
        self.L2_chain.set_chain_id(self.next_chain_id())
        self.struct.add_chain(self.L2_chain)

        self.L3_chain = copy.deepcopy(self.L1_chain)
        self.L3_chain.set_chain_id(self.next_chain_id())
        self.struct.add_chain(self.L3_chain)
        
    def construct_frame(self, phase):
        """Create a new model in self.struct with the TLS displacements
        caused from the three screw dispacement axes displaced by the
        sin(phase).
        """
        ## copy the original model and add it to the structure
        model1 = self.struct.get_model(1)
        
        model = copy.deepcopy(model1)
        model.set_model_id(self.next_model_id())
        self.struct.add_model(model)

        ## now displace the new model according to the tls group
        ## screw axes and given phase
        for tls in self.tlsopt.tls_list:
            self.displace_model(model, tls, phase)

    def displace_model(self, model, tls, phase):
        """Displace the given model by the tls. 
        """
        chain_id = tls["chain_id"]
        chain    = model.get_chain(chain_id)
        
        tls_group = tls["tls_group"]
        tls_info  = tls["tls_info"]
        cor       = tls_info["COR"]

        frag_id1  = tls["frag_id1"]
        frag_id2  = tls["frag_id2"]

        for alt_loc, Lx_val, Lx_vec, Lx_rho, Lx_pitch in [
            ("A", "L1_eigen_val", "L1_eigen_vec", "L1_rho", "L1_pitch"),
            ("B", "L2_eigen_val", "L2_eigen_vec", "L2_rho", "L2_pitch"),
            ("C", "L3_eigen_val", "L3_eigen_vec", "L3_rho", "L3_pitch") ]:

            Lval   = tls_info[Lx_val]
            Lvec   = tls_info[Lx_vec]
            Lrho   = tls_info[Lx_rho]
            Lpitch = tls_info[Lx_pitch]

            ## pre-calculations for screw displacement
            Lorigin = cor + Lrho
            Lrot = GAUSS3C[ADP_PROB]*math.sqrt(Lval)*math.sin(TWO_PI*phase)
            D = dmatrixu(Lvec, Lrot)
            d_screw = (Lrot * Lpitch) * Lvec
            
            self.struct.set_default_alt_loc(alt_loc)

            for frag in iter_fragment_range(chain, frag_id1, frag_id2):
                for atm in frag.iter_atoms():
                    d = matrixmultiply(D, atm.position - Lorigin) + d_screw
                    atm.position += d


class TLSAnimate(object):
    """Create a multi-model PDB file which each model a frame of a TLS
    animation.
    """
    
    def __init__(self, struct, chainopt, tlsopt):
        ## copy and get the anisotropic ADPs out of the structure
        self.struct = self.copy_struct(struct, chainopt)
        
        self.chainopt = chainopt
        self.tlsopt   = tlsopt

        self.L1_chain = None
        self.L2_chain = None
        self.L3_chain = None

        self.construct_chain_copies()

    def copy_struct(self, struct, chainopt):
        """Make a copy of the argument structure to use for generating
        the animation.  Only copy the chain specified in chain_id.
        """
        cp_struct = Structure(structure_id = struct.structure_id)
        chain_id = chainopt["chain_id"]

        for chain in struct.iter_chains():

            if chain.chain_id==chain_id or chain.count_fragments()<200:
                include_chain = True
	    else:
		include_chain = False
		
            for frag in chain.iter_fragments():

                ## skip all waters
                if frag.is_water()==True:
                    continue

                elif frag.is_amino_acid()==True and include_chain:
                    for atm in frag.iter_atoms():

                        ## we only want main-chain atoms necessary for
                        ## generating cartoon/trace
                        if atm.name not in ["N","CA","C","O"]:
                            continue
                        
                        cp_atom = Atom(
                            chain_id    = atm.chain_id,
                            fragment_id = atm.fragment_id,
                            res_name    = atm.res_name,
                            element     = atm.element,
                            name        = atm.name,
                            temp_factor = atm.temp_factor,
                            occupancy   = atm.occupancy,
                            position    = atm.position.copy())
                        cp_struct.add_atom(cp_atom, True)

                elif frag.is_standard_residue()==False:
                    for atm in frag.iter_atoms():
                        cp_atom = Atom(
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
        for phase in (0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875):
            self.construct_frame(phase)
        SaveStructure(struct=self.struct, fil=filename)

    def next_model_id(self):
        """Return the next availible model_id in self.struct
        """
        model_id = 0
        while True:
            model_id += 1            
            model = self.struct.get_model(model_id)
            if model==None:
                return model_id
            
    def next_chain_id(self):
        """Return the next availible chain_id in self.astruct
        """
        for chain_id in string.uppercase:
            chain = self.struct.get_chain(chain_id)
            if chain==None:
                return chain_id
        raise TLSAnimateFailure()

    def construct_chain_copies(self):
        """Create two other copies 
        """
        chain_id = self.chainopt["chain_id"]        
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
        for tls in self.tlsopt.tls_list:
            self.displace_model(model, tls, phase)

    def displace_model(self, model, tls, phase):
        """Displace the given model by the tls. 
        """
        tls_info  = tls["tls_info"]
        cor       = tls_info["COR"]

        frag_id1  = tls["frag_id1"]
        frag_id2  = tls["frag_id2"]

        for n, Lx_val, Lx_vec, Lx_rho, Lx_pitch in [
            (1, "L1_eigen_val", "L1_eigen_vec", "L1_rho", "L1_pitch"),
            (2, "L2_eigen_val", "L2_eigen_vec", "L2_rho", "L2_pitch"),
            (3, "L3_eigen_val", "L3_eigen_vec", "L3_rho", "L3_pitch") ]:

            Lval   = tls_info[Lx_val]
            Lvec   = tls_info[Lx_vec]
            Lrho   = tls_info[Lx_rho]
            Lpitch = tls_info[Lx_pitch]

            ## pre-calculations for screw displacement
            Lorigin = cor + Lrho
            Lrot = GAUSS3C[ADP_PROB]*math.sqrt(Lval)*math.sin(TWO_PI*phase)
            D = dmatrixu(Lvec, Lrot)
            d_screw = (Lrot * Lpitch) * Lvec
            
            if n==1:
                chain_id = self.L1_chain.chain_id
            elif n==2:
                chain_id = self.L2_chain.chain_id
            elif n==3:
                chain_id = self.L3_chain.chain_id

            chain = model.get_chain(chain_id)

            for frag in iter_fragment_range(chain, frag_id1, frag_id2):
                for atm in frag.iter_atoms():
                    d = matrixmultiply(D, atm.position - Lorigin) + d_screw
                    atm.position += d

    def phase_assignment(self):
        """
        """
        print "CHAIN %s NTLS %d" % (
            self.chainopt["chain_id"],
            self.tlsopt.ntls)
        
        tls_list = self.tlsopt.tls_list
        
        n = len(tls_list)
        for i in range(n-1):
            tls1 = tls_list[i]
            tls2 = tls_list[i+1]

            print "  TLS %d - TLS %d" % (i, i+1)

            self.relative_phase_assignment(tls1, tls2)

        print
        print
        
    def relative_phase_assignment(self, tls1, tls2):
        """Make a rough guess at phase assignment based on the bond length
        of the C-N bonds which connect adjecent TLS groups.
        """
        segment1 = tls1["segment"]
        segment2 = tls2["segment"]
        
        frag1 = segment1[-1]
        frag2 = segment2[0]

        atmC = frag1.get_atom("C")
        atmN = frag2.get_atom("N")

        if atmC==None or atmN==None:
            return None

        bond = atmC.get_bond(atmN)
        if bond==None:
            return None

        target_len = bond.calc_length()
        print "  TARGET LEN %6.4f (%6.4f - %6.4f)" % (
            target_len, target_len - 0.02, target_len + 0.02)

        six1 = self.calc_six(tls1, atmC.position)
        six2 = self.calc_six(tls2, atmN.position)

        ## is this matrix symmetric???
        M = zeros((6,6), Float)
        for i in range(6):
            for j in range(6):
                M[i,j] = length(six1[i] - six2[j])

        print M

    def calc_six(self, tls, position):
        """Find the six displacement positions of a atom at the argument
        position and return them in a 6-element list.  See code for order.
        """
        tls_info = tls["tls_info"]
        cor = tls_info["COR"]

        dlist = []

        for n, Lx_val, Lx_vec, Lx_rho, Lx_pitch in [
            (1, "L1_eigen_val", "L1_eigen_vec", "L1_rho", "L1_pitch"),
            (2, "L2_eigen_val", "L2_eigen_vec", "L2_rho", "L2_pitch"),
            (3, "L3_eigen_val", "L3_eigen_vec", "L3_rho", "L3_pitch") ]:

            for phase in [1.0, -1.0]:

                Lval   = tls_info[Lx_val]
                Lvec   = tls_info[Lx_vec]
                Lrho   = tls_info[Lx_rho]
                Lpitch = tls_info[Lx_pitch]

                ## pre-calculations for screw displacement
                Lorigin = cor + Lrho
                Lrot = GAUSS3C[ADP_PROB] * math.sqrt(Lval) * phase
                D = dmatrixu(Lvec, Lrot)
                d_screw = (Lrot * Lpitch) * Lvec
                d = matrixmultiply(D, position - Lorigin) + d_screw
                dposition = d + position

                dlist.append(dposition)

        return dlist
