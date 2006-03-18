## TLS Motion Determination (TLSMD)
## Copyright 2002 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

from mmLib import Structure

import conf
import tlsmdmodule


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
        self.tls_model.set_xmlrpc_chain(xmlrpc_chain)
        print "[%s] set_xmlrpc_chain(num_atoms=%d)" % (self.name, len(xmlrpc_chain))
        return True

    def lsq_fit_segment(self, frag_id1, frag_id2):
        """Performs a LSQ fit of TLS parameters for the protein segment
        starting with fragment index ifrag_start to (and including) the
        fragment ifrag_end.
        """
        ## perform the LSQR fit
        fdict = self.fit_segment(frag_id1, frag_id2)
        fdict["lsq_residual"] = fdict["residual"]
        return fdict

    def fit_segment(self, frag_id1, iend):
        raise Exception()
    

class TLSGraphChainLinear(TLSGraphChain):
    """Graph the chain using the linear TLS parameter fit engine. 
    """
    def __init__(self, name):
        TLSGraphChain.__init__(self, name)        
        self.tls_model = tlsmdmodule.TLSModelAnalyzer()

class TLSGraphChainLinearIsotropic(TLSGraphChainLinear):
    """Linear fit of TLS parameters to isotropically refined ADPs.
    """
    def __init__(self):
        TLSGraphChainLinear.__init__(self, "ISOT")
        self.model_parameters = 10

    def set_xmlrpc_chain(self, xmlrpc_chain):
        TLSGraphChainLinear.set_xmlrpc_chain(self, xmlrpc_chain)

    def fit_segment(self, frag_id1, frag_id2):
        fdict = self.tls_model.isotropic_fit_segment(frag_id1, frag_id2)
        return fdict

class TLSGraphChainLinearAnisotropic(TLSGraphChainLinear):
    """Linear fit of TLS parameters to anisotropically refined ADPs.
    """
    def __init__(self):
        TLSGraphChainLinear.__init__(self, "ANISO")
        self.model_parameters = 20
        
    def fit_segment(self, frag_id1, frag_id2):
        fdict = self.tls_model.anisotropic_fit_segment(frag_id1, frag_id2)
        return fdict


###############################################################################
## Constrained Non-Linear TLS Parameter Fitting Engines
##

class TLSGraphChainNonlinear(TLSGraphChain):
    """Abstract class.
    """
    def __init__(self, name):
        TLSGraphChain.__init__(self, name)
        self.tls_model = tlsmdmodule.TLSModelAnalyzer()
    
class TLSGraphChainNonlinearIsotropic(TLSGraphChainNonlinear):
    """Nonlinear fit of TLS parameters to isotropically refined ADPs.
    """
    def __init__(self):
        TLSGraphChainNonlinear.__init__(self, "NLISOT")
        self.model_parameters = 10
        
    def fit_segment(self, frag_id1, frag_id2):
        fdict = self.tls_model.constrained_isotropic_fit_segment(frag_id1, frag_id2)
        return fdict

class TLSGraphChainNonlinearAnisotropic(TLSGraphChainNonlinear):
    """Nonlinear fit of TLS parameters to anisotropically refined ADPs.
    """
    def __init__(self):
        TLSGraphChainNonlinear.__init__(self, "NLANISO")
        self.model_parameters = 20

    def fit_segment(self, frag_id1, frag_id2):
        fdict = self.tls_model.constrained_anisotropic_fit_segment(frag_id1, frag_id2)
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

