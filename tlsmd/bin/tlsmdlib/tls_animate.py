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


ADP_PROB = 85


def calc_screw_d(cor, Lval, Lvec, Lrho, Lpitch, position, phase):
    """Returns the amount of rotational displacement from L
    for a atom at the given position.
    """
    assert Lval>=0.0
    
    Lrot = GAUSS3C[ADP_PROB] * math.sqrt(Lval) * math.sin(2.0*math.pi*phase)
    Lorigin = cor + Lrho
    D = dmatrixu(Lvec, Lrot)

    drot = matrixmultiply(D, position - Lorigin)
    dscw = (Lrot * Lpitch) * Lvec
    
    return drot + dscw


class TLSAnimateFailure(Exception):
    pass


class TLSAnimate(object):
    """Create a multi-model PDB file which each model a frame of a TLS
    animation.
    """
    
    def __init__(self, struct, chainopt, tlsopt):

        ## copy and get the anisotropic ADPs out of the structure
        self.struct   = copy.deepcopy(struct)
        for atm in self.struct.iter_all_atoms():
            atm.U = None
        
        self.chainopt = chainopt
        self.tlsopt   = tlsopt

        self.L1_chain = None
        self.L2_chain = None
        self.L3_chain = None

        self.construct_chain_copies()

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
        tls_group = tls["tls_group"]
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

            if n==1:
                chain_id = self.L1_chain.chain_id
            elif n==2:
                chain_id = self.L2_chain.chain_id
            elif n==3:
                chain_id = self.L3_chain.chain_id

            chain = model.get_chain(chain_id)
            segment = chain[frag_id1:frag_id2]

            for atm in segment.iter_all_atoms():
                d = calc_screw_d(
                    cor, Lval, Lvec, Lrho, Lpitch, atm.position, phase)
                atm.position += d
