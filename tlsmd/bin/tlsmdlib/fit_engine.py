## TLS Minimized Domains (TLSMD)
## Copyright 2002 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import sys
import copy
import string

from mmLib import Structure

import conf
import lineartls
import nonlineartls



###############################################################################
## Utility Class
##

class XChain(object):
    """
    """
    def __init__(self, xmlrpc_chain):
        self.xmlrpc_chain = xmlrpc_chain
        self.num_atoms = len(xmlrpc_chain)

        ## construct a list of fragment IDs
        last_frag_id = None
        frag_id_list = []
        for atm_desc in xmlrpc_chain:
            frag_id = atm_desc["frag_id"]
            if frag_id!=last_frag_id:
                frag_id_list.append(frag_id)
                last_frag_id = frag_id

        self.frag_id_list = frag_id_list
        self.num_frags = len(frag_id_list)

        ## construct a 1:1 list of frag_id_list with istart/iend indexes
        istart_list = [None for x in range(len(frag_id_list))]
        iend_list = [None for x in range(len(frag_id_list))]

        for i in range(len(frag_id_list)):
            frag_id = frag_id_list[i]

            istart, iend = self.__find_istart_iend(frag_id, frag_id)
            istart_list[i] = istart
            iend_list[i] = iend

            assert xmlrpc_chain[istart]["frag_id"]==frag_id
            assert xmlrpc_chain[iend]["frag_id"]==frag_id

        self.istart_list = istart_list
        self.iend_list = iend_list

    def get_istart(self, frag_id):
        try:
            i = self.frag_id_list.index(frag_id)
        except ValueError:
            return None
        return self.istart_list[i]

    def get_iend(self, frag_id):
        try:
            i = self.frag_id_list.index(frag_id)
        except ValueError:
            return None
        return self.iend_list[i]
    
    def __find_istart_iend(self, frag_id1, frag_id2):
        istart = None
        iend   = None
        state  = "find_istart"

        for icur in range(len(self.xmlrpc_chain)):
            if state=="find_istart":
                if Structure.fragment_id_ge(self.xmlrpc_chain[icur]["frag_id"], frag_id1):
                    state  = "find_iend"
                    istart = icur
            elif state=="find_iend":
                if Structure.fragment_id_gt(self.xmlrpc_chain[icur]["frag_id"], frag_id2):
                    iend = icur - 1
                    break

        if istart==None:
            return None, None

        if iend==None:
            iend = len(self.xmlrpc_chain) - 1

        return istart, iend


###############################################################################
## Unconstrained Linear TLS Parameter Fitting Engines
##

class TLSGraphChain(object):
    """Stores a protein chain of atoms in the tuple list xmlrpc_chain along
    with the residue(fragment) indexes and U ADP tensor of each atom.
    Subsegments of the protein chain then can be fit with TLS tensors.
    """
    def __init__(self, name):
        self.name = name
        self.tls_model = None
        self.model_parameters = 0
        self.verbose = conf.globalconf.verbose
    
    def set_xmlrpc_chain(self, xmlrpc_chain):
        self.xchain = XChain(xmlrpc_chain)
        self.tls_model.set_xmlrpc_chain(xmlrpc_chain)
        print "[%s] set_xmlrpc_chain(num_atoms=%d)" % (self.name, len(xmlrpc_chain))
        return True

    def lsq_fit_segment(self, frag_id1, frag_id2):
        """Performs a LSQ fit of TLS parameters for the protein segment
        starting with fragment index ifrag_start to (and including) the
        fragment ifrag_end.
        """
        ## all return values here
        fit_info = {}
        
        ## calculate the start/end indexes of the start fragment
        ## and end fragment so the A matrix and b vector can be sliced
        ## in the correct places
        istart = self.xchain.get_istart(frag_id1)
        iend = self.xchain.get_iend(frag_id2)
        if istart==None or iend==None:
            fit_info["error"] = "No Atoms In Segment"
            return fit_info

        ## are there enough atoms in this chain segment
        num_atoms = iend - istart + 1
        fit_info["num_atoms"] = num_atoms
        if num_atoms<self.model_parameters:
            fit_info["error"] = "data/parameter raito = %d/%d less than 1.0" % (num_atoms, self.model_parameters)
            return fit_info

        ## perform the LSQR fit
        fdict = self.fit_segment(istart, iend)

        ## return information
        fit_info["lsq_residual"] = fdict["lsq_residual"]

        if self.verbose:
            print "[%s] lsq_fit_segment(frag_id={%s..%s}, num_atoms=%d, lsqr=%8.6f, gammaq=%6.4f)" % (
                self.name, frag_id1, frag_id2, fit_info["num_atoms"], fit_info["lsq_residual"], fdict["gammaq"])

        return fit_info


class TLSGraphChainLinear(TLSGraphChain):
    """Graph the chain using the linear TLS parameter fit engine. 
    """
    def __init__(self, name):
        TLSGraphChain.__init__(self, name)        
        self.tls_model = lineartls.LinearTLSModel()

class TLSGraphChainLinearIsotropic(TLSGraphChainLinear):
    """Linear fit of TLS parameters to isotropically refined ADPs.
    """
    def __init__(self):
        TLSGraphChainLinear.__init__(self, "ISOT")
        self.model_parameters = 10

    def set_xmlrpc_chain(self, xmlrpc_chain):
        TLSGraphChainLinear.set_xmlrpc_chain(self, xmlrpc_chain)

    def fit_segment(self, istart, iend):
        fdict = self.tls_model.isotropic_fit_segment(istart, iend)
        fdict["lsq_residual"] = fdict["ilsqr"]
        return fdict

class TLSGraphChainLinearAnisotropic(TLSGraphChainLinear):
    """Linear fit of TLS parameters to anisotropically refined ADPs.
    """
    def __init__(self):
        TLSGraphChainLinear.__init__(self, "ANISO")
        self.model_parameters = 20
        
    def fit_segment(self, istart, iend):
        fdict = self.tls_model.anisotropic_fit_segment(istart, iend)
        fdict["lsq_residual"] = fdict["alsqr"]
        return fdict


###############################################################################
## Constrained Non-Linear TLS Parameter Fitting Engines
##

class TLSGraphChainNonlinear(TLSGraphChain):
    """Abstract class.
    """
    def __init__(self, name):
        TLSGraphChain.__init__(self, name)
        self.tls_model = nonlineartls.NLTLSModel()
    
class TLSGraphChainNonlinearIsotropic(TLSGraphChainNonlinear):
    """Nonlinear fit of TLS parameters to isotropically refined ADPs.
    """
    def __init__(self):
        TLSGraphChainNonlinear.__init__(self, "NLISOT")
        self.model_parameters = 10
        
    def fit_segment(self, istart, iend):
        fdict = self.tls_model.isotropic_fit_segment(istart, iend)
        fdict["lsq_residual"] = fdict["ilsqr"]
        return fdict

class TLSGraphChainNonlinearAnisotropic(TLSGraphChainNonlinear):
    """Nonlinear fit of TLS parameters to anisotropically refined ADPs.
    """
    def __init__(self):
        TLSGraphChainNonlinear.__init__(self, "NLANISO")
        self.model_parameters = 20

    def fit_segment(self, istart, iend):
        fdict = self.tls_model.anisotropic_fit_segment(istart, iend)
        fdict["lsq_residual"] = fdict["alsqr"]
        return fdict


###############################################################################
## interface...
##

def NewTLSGraphChain(tls_model):
    """Generate and return the proper TLSGraphChain subclass for the
    requested TLS model.
    """
    if tls_model=="ISOT":
        return TLSGraphChainLinearIsotropic()

    elif tls_model=="ANISO":
        return TLSGraphChainLinearAnisotropic()

    elif tls_model=="NLISOT":
        return TLSGraphChainNonlinearIsotropic()

    elif tls_model=="NLANISO":
        return TLSGraphChainNonlinearAnisotropic()

    raise Exception()

