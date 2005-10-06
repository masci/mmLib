## TLS Minimized Domains (TLSMD)
## Copyright 2002 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import sys
import copy
import string

from misc     import *
from datafile import TLSMDFile
from hcsssp   import HCSSSP

## C/FORTRAN TLS Fitting Engines
import lineartls
import nonlineartls


###############################################################################
## Unconsrained Linear TLS Parameter Fitting Engines
##

class TLSGraphChain(object):
    """Stores a protein chain of atoms in the tuple list xmlrpc_chain along
    with the residue(fragment) indexes and U ADP tensor of each atom.
    Subsegments of the protein chain then can be fit with TLS tensors.
    """
    pass

class TLSGraphChainLinear(TLSGraphChain):
    """Graph the chain using the linear TLS parameter fit engine. 
    """
    def __init__(self, name):
        TLSGraphChain.__init__(self)
        
        self.name         = name
        self.xmlrpc_chain = None
        self.f            = None
        self.tls_model    = lineartls.LinearTLSModel()
    
    def set_xmlrpc_chain(self, xmlrpc_chain):
        self.xmlrpc_chain = xmlrpc_chain

        ## create a Numeric Array vector for each atom position (faster)
        frag_id_list = []
        
        for atm_desc in xmlrpc_chain:
            frag_id_list.append(atm_desc["frag_id"])

        ## set the class frag list
        self.f = frag_id_list

        ## set the C solver
        self.tls_model.set_xmlrpc_chain(xmlrpc_chain)

        print "[%s] set_xmlrpc_chain(num_atoms=%d)" % (self.name, len(frag_id_list))
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
        ## in the correct placees
	istart = None
        iend   = None
        state  = "find_istart"

        for icur in range(len(self.f)):
            if state=="find_istart":
                if fragment_id_ge(self.f[icur], frag_id1):
                    state  = "find_iend"
                    istart = icur
            elif state=="find_iend":
                if fragment_id_gt(self.f[icur], frag_id2):
                    iend = icur - 1
                    break

        if istart==None:
            fit_info["error"] = "No Atoms In Segment"
            return fit_info

	if iend==None:
	    iend = len(self.f) - 1

        ## are there enough atoms in this chain segment
        num_atoms = iend - istart + 1
        fit_info["num_atoms"] = num_atoms
        if num_atoms<20:
            fit_info["error"] = "data/parameter raito = %d/20 less than 1.0" % (num_atoms)
            return fit_info

        ## perform the LSQR fit
        fdict = self.fit_segment(istart, iend)

        ## return information
        fit_info["lsq_residual"] = fdict["lsq_residual"]

        fit_info["cor_x"] = fdict["x"]
        fit_info["cor_y"] = fdict["y"]
        fit_info["cor_z"] = fdict["z"]

        for key in ["t11","t22","t33","t12","t13","t23",
                    "l11","l22","l33","l12","l13","l23",
                    "s2211", "s1133", "s12","s13","s23","s21","s31","s32"]:
            fit_info[key] = fdict[key]

        print "[%s] lsq_fit_segment(frag_id={%s..%s}, num_atoms=%d, lsqr=%8.6f)" % (
            self.name, frag_id1, frag_id2, fit_info["num_atoms"], fit_info["lsq_residual"])

        return fit_info


class TLSGraphChainLinearIsotropic(TLSGraphChainLinear):
    """Linear fit of TLS parameters to isotropically refined ADPs.
    """
    def __init__(self):
        TLSGraphChainLinear.__init__(self, "ISOT")

    def fit_segment(self, istart, iend):
        fdict = self.tls_model.isotropic_fit_segment(istart, iend)
        fdict["lsq_residual"] = fdict["ilsqr"]
        return fdict


class TLSGraphChainLinearAnisotropic(TLSGraphChainLinear):
    """Linear fit of TLS parameters to anisotropically refined ADPs.
    """
    def __init__(self):
        TLSGraphChainLinear.__init__(self, "ANISO")

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
        TLSGraphChain.__init__(self)
        
        self.name         = name
        self.xmlrpc_chain = None
        self.f            = None

        ## the non-linear fit engine is implemented in C/FORTRAN
        import nonlineartls
        self.tls_model = nonlineartls.NLTLSModel()
    
    def set_xmlrpc_chain(self, xmlrpc_chain):
        self.xmlrpc_chain = xmlrpc_chain

        ## create a Numeric Array vector for each atom position (faster)
        frag_id_list = []
        
        for atm_desc in xmlrpc_chain:
            frag_id_list.append(atm_desc["frag_id"])

        ## set the class frag list
        self.f = frag_id_list

        ## set the C solver
        self.tls_model.set_xmlrpc_chain(xmlrpc_chain)

        print "[%s] set_xmlrpc_chain(num_atoms=%d)" % (self.name, len(frag_id_list))
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
        ## in the correct placees
	istart = None
        iend   = None
        state  = "find_istart"

        for icur in range(len(self.f)):
            if state=="find_istart":
                if fragment_id_ge(self.f[icur], frag_id1):
                    state  = "find_iend"
                    istart = icur
            elif state=="find_iend":
                if fragment_id_gt(self.f[icur], frag_id2):
                    iend = icur - 1
                    break

        if istart==None:
            fit_info["error"] = "No Atoms In Segment"
            return fit_info

	if iend==None:
	    iend = len(self.f) - 1

        ## are there enough atoms in this chain segment
        num_atoms = iend - istart + 1
        fit_info["num_atoms"] = num_atoms
        if num_atoms<20:
            fit_info["error"] = "data/parameter raito = %d/20 less than 1.0" % (num_atoms)
            return fit_info

        ## perform the LSQR fit
        fdict = self.fit_segment(istart, iend)
        
        ## return information
        fit_info["lsq_residual"] = fdict["lsq_residual"]

        fit_info["cor_x"] = fdict["x"]
        fit_info["cor_y"] = fdict["y"]
        fit_info["cor_z"] = fdict["z"]

        for key in ["t11","t22","t33","t12","t13","t23",
                    "l11","l22","l33","l12","l13","l23",
                    "s2211", "s1133", "s12","s13","s23","s21","s31","s32"]:
            fit_info[key] = fdict[key]

        print "[%s] lsq_fit_segment(frag_id={%s..%s}, num_atoms=%d, lsqr=%6.4f)" % (
            self.name, frag_id1, frag_id2, fit_info["num_atoms"], fit_info["lsq_residual"])

        return fit_info

    def fit_segment(self, istart, iend):
        """Implement me.
        """
        pass

    
class TLSGraphChainNonlinearIsotropic(TLSGraphChainNonlinear):
    """Nonlinear fit of TLS parameters to isotropically refined ADPs.
    """
    def __init__(self):
        TLSGraphChainNonlinear.__init__(self, "NLISOT")

    def fit_segment(self, istart, iend):
        fdict = self.tls_model.isotropic_fit_segment(istart, iend)
        fdict["lsq_residual"] = fdict["ilsqr"]
        return fdict


class TLSGraphChainNonlinearAnisotropic(TLSGraphChainNonlinear):
    """Nonlinear fit of TLS parameters to anisotropically refined ADPs.
    """
    def __init__(self):
        TLSGraphChainNonlinear.__init__(self, "NLANISO")

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
    if tls_model=="HYBRID" or tls_model=="ISOT":
        return TLSGraphChainLinearIsotropic()

    elif tls_model=="ANISO":
        return TLSGraphChainLinearAnisotropic()

    elif tls_model=="NLISOT":
        return TLSGraphChainNonlinearIsotropic()

    elif tls_model=="NLANISO":
        return TLSGraphChainNonlinearAnisotropic()

    raise Exception()


class TLSGraphChainXMLRPCServer(TLSGraphChain):
    """Runs this object as a xmlrpc server servicing requests for
    functions: set_xmlrpc_chain, and lsq_fit_segment.
    """
    def __init__(self):
        self.proxy = None

    def set_tls_model(self, tls_model, weight_model):
        print "TLSGraphChainXMLRPCServer.set_tls_model(%s, %s)" % (tls_model, weight_model)

        ## select proper TLS model
        self.proxy = NewTLSGraphChain(tls_model)
        
        return True

    def set_xmlrpc_chain(self, xmlrpc_chain):
        return self.proxy.set_xmlrpc_chain(xmlrpc_chain)

    def lsq_fit_segment(self, frag_id1, frag_id2):
        return self.proxy.lsq_fit_segment(frag_id1, frag_id2)

    def run_server(self, host_port):
        xmlrpc_server = SimpleXMLRPCServer.SimpleXMLRPCServer(
            host_port,
            SimpleXMLRPCServer.SimpleXMLRPCRequestHandler,
            False)
        
        xmlrpc_server.register_function(self.set_tls_model,    "set_tls_model")
        xmlrpc_server.register_function(self.set_xmlrpc_chain, "set_xmlrpc_chain")
        xmlrpc_server.register_function(self.lsq_fit_segment,   "lsq_fit_segment")
        
        xmlrpc_server.serve_forever()

