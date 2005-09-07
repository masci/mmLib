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
try:
    import lineartls
except ImportError:
    USE_LINEARTLS = False
else:
    USE_LINEARTLS = True


try:
    import nonlineartls
except ImportError:
    USE_NONLINEARTLS = False
else:
    USE_NONLINEARTLS = True


###############################################################################
## Unconsrained Linear TLS Parameter Fitting Engines
##

class TLSGraphChain(object):
    """Stores a protein chain of atoms in the tuple list xmlrpc_chain along
    with the residue(fragment) indexes and U ADP tensor of each atom.
    Subsegments of the protein chain then can be fit with TLS tensors.
    """
    pass


class TLSGraphChainHybrid(TLSGraphChain):
    """Graph the chain with a hybrid isotropic/anisotropic TLS model. 
    This is the slow Python version.
    """
    def __init__(self):
        TLSGraphChain.__init__(self)
        
        self.xmlrpc_chain = None
        self.f            = None
    
    def set_xmlrpc_chain(self, xmlrpc_chain):
        """Sets a new TLS chain. which generates new self.A, self.b,
        and self.f arrays.
        """
        self.xmlrpc_chain = xmlrpc_chain

        ## create a Numeric Array vector for each atom position (faster)
        frag_id_list = []
        
        for atm_desc in xmlrpc_chain:
            frag_id_list.append(atm_desc["frag_id"])
            
            atm_desc["position"] = array((atm_desc["x"], atm_desc["y"], atm_desc["z"]), Float)
            atm_desc["sqrt_w"] = math.sqrt(atm_desc["w"])

        ## set the class frag list
        self.f = frag_id_list

        print "[HYBRID] set_xmlrpc_chain(num_atoms=%d)" % (len(frag_id_list))
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

        ## ensure enough data points to solve for parameters
        if num_atoms<20:
            fit_info["error"] = "data/parameter raito = %d/20 less than 1.0" % (num_atoms)
            return fit_info

        xmlrpc_chain = self.xmlrpc_chain

        ## CALCULATE CENTROID
        centroid = zeros(3, Float)
        ia = istart - 1
        while ia<iend:
            ia += 1
            centroid += xmlrpc_chain[ia]["position"]
        centroid /= num_atoms


        ## SOLVE ISOTROPIC AND ANISOTROPIC TLS MODELs
        A_ISOW = zeros((num_atoms, 13), Float)
        B_ISOW = zeros(num_atoms, Float)

        A_ANISOW = zeros((num_atoms*6, 20), Float)
        B_ANISOW = zeros(num_atoms*6, Float)
        
        i = -1
        ia = istart - 1
        while ia<iend:
            i += 1
            ia += 1
            atm_desc = xmlrpc_chain[ia]

            ## calculate atom position relative to the centroid
            x, y, z = atm_desc["position"] - centroid
            
            ## w is actually w^2
            w = atm_desc["sqrt_w"]

            ## set the A Matrix, B vector
            set_TLSiso_A(A_ISOW, i, 0, x, y, z, w)
            set_TLSiso_b(B_ISOW, i, atm_desc["u_iso"], w)

            set_TLS_A(A_ANISOW, i*6, 0, x, y, z, w)
            set_TLS_b(B_ANISOW, i*6,
                      atm_desc["u11"], atm_desc["u22"], atm_desc["u33"],
                      atm_desc["u12"], atm_desc["u13"], atm_desc["u23"], w)
            
        X_ISO  = solve_TLS_Ab(A_ISOW, B_ISOW)
        U_ISOW = matrixmultiply(A_ISOW, X_ISO)

        ## calculate the lsq residual from the isotropic model
        D_ISOW = U_ISOW - B_ISOW
        fit_info["lsq_residual"] = dot(D_ISOW, D_ISOW)
        
        ## everything else comes from the anisotropic model
        X_ANISO = solve_TLS_Ab(A_ANISOW, B_ANISOW)

        T11, T22, T33, T12, T13, T23, L11, L22, L33, L12, L13, L23, S1133, S2211, S12, S13, S23, S21, S31, S32 = (
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)

        T = array([ [ X_ANISO[T11], X_ANISO[T12], X_ANISO[T13] ],
                    [ X_ANISO[T12], X_ANISO[T22], X_ANISO[T23] ],
                    [ X_ANISO[T13], X_ANISO[T23], X_ANISO[T33] ] ], Float)

        L = array([ [ X_ANISO[L11], X_ANISO[L12], X_ANISO[L13] ],
                    [ X_ANISO[L12], X_ANISO[L22], X_ANISO[L23] ],
                    [ X_ANISO[L13], X_ANISO[L23], X_ANISO[L33] ] ], Float)

        s11, s22, s33 = calc_s11_s22_s33(X_ANISO[S2211], X_ANISO[S1133])
	
        S = array([ [          s11, X_ANISO[S12], X_ANISO[S13] ],
                    [ X_ANISO[S21],          s22, X_ANISO[S23] ],
                    [ X_ANISO[S31], X_ANISO[S32],        s33 ] ], Float)

        ## caculate the tensors shifted to the center of reaction
        cor_info = calc_TLS_center_of_reaction(T, L, S, centroid)

        cor    = cor_info["COR"]
        T_cor  = cor_info["T'"]
        L_cor  = cor_info["L'"]
        S_cor  = cor_info["S'"]
        T_red  = cor_info["rT'"]

        ## return information
        fit_info["cor_x"] = cor[0]
        fit_info["cor_y"] = cor[1]
        fit_info["cor_z"] = cor[2]

        fit_info["t11"] = T_cor[0,0]
        fit_info["t22"] = T_cor[1,1]
        fit_info["t33"] = T_cor[2,2]
        fit_info["t12"] = T_cor[0,1]
        fit_info["t13"] = T_cor[0,2]
        fit_info["t23"] = T_cor[1,2]

        fit_info["l11"] = L_cor[0,0]
        fit_info["l22"] = L_cor[1,1]
        fit_info["l33"] = L_cor[2,2]
        fit_info["l12"] = L_cor[0,1]
        fit_info["l13"] = L_cor[0,2]
        fit_info["l23"] = L_cor[1,2]

        fit_info["s2211"] = S_cor[1,1] - S_cor[0,0]
        fit_info["s1133"] = S_cor[0,0] - S_cor[2,2]
        fit_info["s12"]   = S_cor[0,1]
        fit_info["s13"]   = S_cor[0,2]
        fit_info["s23"]   = S_cor[1,2]
        fit_info["s21"]   = S_cor[1,0]
        fit_info["s31"]   = S_cor[2,0]
        fit_info["s32"]   = S_cor[2,1]

        if min(eigenvalues(T_red))<=TSMALL:
            errx = "Invalid Tr Eigenvalue"
            fit_info["error"] = errx

        elif min(eigenvalues(L_cor))<=LSMALL:
            errx = "Invalid L Eigenvalue"
            fit_info["error"] = errx

        if fit_info.has_key("error"):
            print "[HYBRID] lsq_fit_segment(frag_id={%s..%s}, num_atoms=%d, lsqr=%6.4f, discard=%s)" % (
                frag_id1, frag_id2, fit_info["num_atoms"], fit_info["lsq_residual"], fit_info["error"]) 
        else:
            print "[HYBRID] lsq_fit_segment(frag_id={%s..%s}, num_atoms=%d, lsqr=%6.4f)" % (
                frag_id1, frag_id2, fit_info["num_atoms"], fit_info["lsq_residual"])

        return fit_info


class TLSGraphChainFastHybrid(TLSGraphChain):
    """Graph the chain with a hybrid isotropic/anisotropic TLS model. 
    """
    def __init__(self):
        TLSGraphChain.__init__(self)
        
        self.xmlrpc_chain = None
        self.f            = None
        self.itls_model   = lineartls.ITLSModel()
    
    def set_xmlrpc_chain(self, xmlrpc_chain):
        self.xmlrpc_chain = xmlrpc_chain

        ## create a Numeric Array vector for each atom position (faster)
        frag_id_list = []
        
        for atm_desc in xmlrpc_chain:
            frag_id_list.append(atm_desc["frag_id"])
            atm_desc["sqrt_w"] = math.sqrt(atm_desc["w"])

        ## set the class frag list
        self.f = frag_id_list

        ## set the C solver
        self.itls_model.set_xmlrpc_chain(xmlrpc_chain)

        print "[FAST-HYBRID] set_xmlrpc_chain(num_atoms=%d)" % (len(frag_id_list))
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
        fdict = self.itls_model.fit_segment(istart, iend)

        fit_info["lsq_residual"] = fdict["ilsqr"]

        T = array([ [ fdict["t11"], fdict["t12"], fdict["t13"] ],
                    [ fdict["t12"], fdict["t22"], fdict["t23"] ],
                    [ fdict["t13"], fdict["t23"], fdict["t33"] ] ], Float)

        L = array([ [ fdict["l11"], fdict["l12"], fdict["l13"] ],
                    [ fdict["l12"], fdict["l22"], fdict["l23"] ],
                    [ fdict["l13"], fdict["l23"], fdict["l33"] ] ], Float)
  
        s11, s22, s33 = calc_s11_s22_s33(fdict["s2211"], fdict["s1133"])
	
        S = array([ [          s11, fdict["s12"], fdict["s13"] ],
                    [ fdict["s21"],          s22, fdict["s23"] ],
                    [ fdict["s31"], fdict["s32"],        s33 ] ], Float)

        centroid = array([ fdict["x"], fdict["y"], fdict["z"] ], Float)

        ## caculate the tensors shifted to the center of reaction
        cor_info = calc_TLS_center_of_reaction(T, L, S, centroid)

        cor    = cor_info["COR"]
        T_cor  = cor_info["T'"]
        L_cor  = cor_info["L'"]
        S_cor  = cor_info["S'"]
        T_red  = cor_info["rT'"]

        ## return information
        fit_info["cor_x"] = cor[0]
        fit_info["cor_y"] = cor[1]
        fit_info["cor_z"] = cor[2]

        fit_info["t11"] = T_cor[0,0]
        fit_info["t22"] = T_cor[1,1]
        fit_info["t33"] = T_cor[2,2]
        fit_info["t12"] = T_cor[0,1]
        fit_info["t13"] = T_cor[0,2]
        fit_info["t23"] = T_cor[1,2]

        fit_info["l11"] = L_cor[0,0]
        fit_info["l22"] = L_cor[1,1]
        fit_info["l33"] = L_cor[2,2]
        fit_info["l12"] = L_cor[0,1]
        fit_info["l13"] = L_cor[0,2]
        fit_info["l23"] = L_cor[1,2]

        fit_info["s2211"] = S_cor[1,1] - S_cor[0,0]
        fit_info["s1133"] = S_cor[0,0] - S_cor[2,2]
        fit_info["s12"]   = S_cor[0,1]
        fit_info["s13"]   = S_cor[0,2]
        fit_info["s23"]   = S_cor[1,2]
        fit_info["s21"]   = S_cor[1,0]
        fit_info["s31"]   = S_cor[2,0]
        fit_info["s32"]   = S_cor[2,1]

        if min(eigenvalues(T_red))<=TSMALL:
            errx = "Invalid Tr Eigenvalue"
            fit_info["error"] = errx

        elif min(eigenvalues(L_cor))<=LSMALL:
            errx = "Invalid L Eigenvalue"
            fit_info["error"] = errx

        if fit_info.has_key("error"):
            print "[FAST-HYBRID] lsq_fit_segment(frag_id={%s..%s}, num_atoms=%d, lsqr=%6.4f, discard=%s)" % (
                frag_id1, frag_id2, fit_info["num_atoms"], fit_info["lsq_residual"], fit_info["error"]) 
        else:
            print "[FAST-HYBRID] lsq_fit_segment(frag_id={%s..%s}, num_atoms=%d, lsqr=%6.4f)" % (
                frag_id1, frag_id2, fit_info["num_atoms"], fit_info["lsq_residual"])

        return fit_info


class TLSGraphChainAnisotropic(TLSGraphChain):
    """Graph the chain using the anisotropic TLS model.
    """
    def __init__(self):
        TLSGraphChain.__init__(self)
        
        self.xmlrpc_chain = None
        self.A            = None
        self.b            = None
        self.o            = None
        self.f            = None
        
    def set_xmlrpc_chain(self, xmlrpc_chain):
        """Sets a new TLS chain. which generates new self.A, self.b,
        and self.f arrays.
        """
        self.xmlrpc_chain = xmlrpc_chain

        ## calculate origin using chain centroid
        n = 0
        centroid_sum = zeros(3, Float)
        for atm_desc in xmlrpc_chain:
            n += 1
            
            centroid_sum[0] += atm_desc["x"]
            centroid_sum[1] += atm_desc["y"]
            centroid_sum[2] += atm_desc["z"]

        centroid = centroid_sum / float(n)
        self.o = centroid

        A, b, f = self.generate_Abd_aniso(xmlrpc_chain, self.o)

        self.A  = A
        self.b  = b
	self.f  = f

        print "[ANISO] set_xmlrpc_chain(A=%s, b=%s, f=%d)" % (shape(self.A), shape(self.b), len(f))

        return True

    def generate_Abd_aniso(self, xmlrpc_chain, origin):
        """Return the 3-tuple (A, b, d).  The A matrix is the TLS position
        dependent coefficents(6 rows per atom), the b vector is the
        experimental U ADP values from xmlrpc_chain, and the d vector is the
        fragment number each row belongs to.
        """
        A  = zeros((len(xmlrpc_chain)*6, 20), Float)
        b  = zeros(len(xmlrpc_chain)*6,  Float)
        f  = []

        i = -1
        for atm_desc in xmlrpc_chain:
            i += 1
            iU11 = i * 6

            ## w is actually w^2
            sqrt_w  = math.sqrt(atm_desc["w"])

            ## keep a list containing a 1:1 mapping
            ## of A,b rows to f frag_ids
            frag_id = atm_desc["frag_id"]
	    for k in range(6):
	        f.append(frag_id)

            ## set x, y, z as the vector components from the TLS origin
            x = atm_desc["x"] - origin[0]
            y = atm_desc["y"] - origin[1]
            z = atm_desc["z"] - origin[2]

            ## set the b vector
            set_TLS_b(b, iU11,
                      atm_desc["u11"], atm_desc["u22"], atm_desc["u33"],
                      atm_desc["u12"], atm_desc["u13"], atm_desc["u23"], sqrt_w)

            ## set the A matrix
            set_TLS_A(A, iU11, 0, x, y, z, sqrt_w)

        return A, b, f

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
                    iend = icur
                    break

        if istart==None:
            fit_info["error"] = "No Atoms In Segment"
            return fit_info
	    
	if iend==None:
	    iend = len(self.f)

        ## now slice A and b
        Aw = self.A[istart:iend,:]
        Bw = self.b[istart:iend]

        ## deal with the no-atoms case
        if len(Bw)==1:
            fit_info["error"] = "No suitable atoms in segment"
            return fit_info
       
        ## double-check six rows per atom
        try:
            assert len(Bw)%6==0
        except AssertionError, err:
            print "AssertionError: len(%d)%%6==0" % (len(Bw))
            raise

        ## calculate number of atoms
        num_atoms = len(Bw)/6
	fit_info["num_atoms"] = len(Bw)/6

        ## check if there are enough atoms
        if num_atoms<4:
            fit_info["error"] = "data/parameter raito = %d/20 less than 1.0" % (num_atoms*6)
            return fit_info

        ## LSQ Fit
        X = solve_TLS_Ab(Aw, Bw)

        ## calculate the weighted and unweighted TLS-predicted Uij values
        Uw = matrixmultiply(Aw, X)

        ## calculate the lsq residual since silly Numeric Python won't
        ## do it for us
        Dw = Uw - Bw
        fit_info["lsq_residual"] = dot(Dw, Dw)

	## shift TLS tensors to the center of reaction
        T11, T22, T33, T12, T13, T23, L11, L22, L33, L12, L13, L23, S1133, S2211, S12, S13, S23, S21, S31, S32 = (
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)

        T = array([ [ X[T11], X[T12], X[T13] ],
                    [ X[T12], X[T22], X[T23] ],
                    [ X[T13], X[T23], X[T33] ] ], Float)

        L = array([ [ X[L11], X[L12], X[L13] ],
                    [ X[L12], X[L22], X[L23] ],
                    [ X[L13], X[L23], X[L33] ] ], Float)

        s11, s22, s33 = calc_s11_s22_s33(X[S2211], X[S1133])
	
        S = array([ [    s11, X[S12], X[S13] ],
                    [ X[S21],    s22, X[S23] ],
                    [ X[S31], X[S32],    s33 ] ], Float)

        ## caculate the tensors shifted to the center of reaction
        cor_info = calc_TLS_center_of_reaction(T, L, S, self.o)

        cor    = cor_info["COR"]
        T_cor  = cor_info["T'"]
        L_cor  = cor_info["L'"]
        S_cor  = cor_info["S'"]
        T_red  = cor_info["rT'"]
        
        ## return information
        fit_info["cor_x"] = cor[0]
        fit_info["cor_y"] = cor[1]
        fit_info["cor_z"] = cor[2]

        fit_info["t11"] = T_cor[0,0]
        fit_info["t22"] = T_cor[1,1]
        fit_info["t33"] = T_cor[2,2]
        fit_info["t12"] = T_cor[0,1]
        fit_info["t13"] = T_cor[0,2]
        fit_info["t23"] = T_cor[1,2]

        fit_info["l11"] = L_cor[0,0]
        fit_info["l22"] = L_cor[1,1]
        fit_info["l33"] = L_cor[2,2]
        fit_info["l12"] = L_cor[0,1]
        fit_info["l13"] = L_cor[0,2]
        fit_info["l23"] = L_cor[1,2]

        fit_info["s2211"] = S_cor[1,1] - S_cor[0,0]
        fit_info["s1133"] = S_cor[0,0] - S_cor[2,2]
        fit_info["s12"]   = S_cor[0,1]
        fit_info["s13"]   = S_cor[0,2]
        fit_info["s23"]   = S_cor[1,2]
        fit_info["s21"]   = S_cor[1,0]
        fit_info["s31"]   = S_cor[2,0]
        fit_info["s32"]   = S_cor[2,1]

        if min(eigenvalues(T_red))<=TSMALL:
            errx = "Invalid Tr Eigenvalue"
            fit_info["error"] = errx

        elif min(eigenvalues(L_cor))<=LSMALL:
            errx = "Invalid L Eigenvalue"
            fit_info["error"] = errx

        if fit_info.has_key("error"):
            print "[ANISO] lsq_fit_segment(frag_id={%s..%s}, num_atoms=%d, lsqr=%6.4f, discard=%s)" % (
                frag_id1, frag_id2, fit_info["num_atoms"], fit_info["lsq_residual"], fit_info["error"]) 
        else:
            print "[ANISO] lsq_fit_segment(frag_id={%s..%s}, num_atoms=%d, lsqr=%6.4f)" % (
                frag_id1, frag_id2, fit_info["num_atoms"], fit_info["lsq_residual"]) 

        return fit_info


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
            atm_desc["sqrt_w"] = math.sqrt(atm_desc["w"])

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
        print "info = %d" % (fdict["info"])
        
        fit_info["lsq_residual"] = fdict["lsq_residual"]

        T = array([ [ fdict["t11"], fdict["t12"], fdict["t13"] ],
                    [ fdict["t12"], fdict["t22"], fdict["t23"] ],
                    [ fdict["t13"], fdict["t23"], fdict["t33"] ] ], Float)

        L = array([ [ fdict["l11"], fdict["l12"], fdict["l13"] ],
                    [ fdict["l12"], fdict["l22"], fdict["l23"] ],
                    [ fdict["l13"], fdict["l23"], fdict["l33"] ] ], Float)
  
        s11, s22, s33 = calc_s11_s22_s33(fdict["s2211"], fdict["s1133"])
	
        S = array([ [          s11, fdict["s12"], fdict["s13"] ],
                    [ fdict["s21"],          s22, fdict["s23"] ],
                    [ fdict["s31"], fdict["s32"],        s33 ] ], Float)

        centroid = array([ fdict["x"], fdict["y"], fdict["z"] ], Float)

        levals = eigenvalues(L) * RAD2DEG2
        tevals = eigenvalues(T) * U2B
        print "T: %8.4f %8.4f %8.4f" % (tevals[0],tevals[1],tevals[2])
        print "L: %8.4f %8.4f %8.4f" % (levals[0],levals[1],levals[2])

        ## caculate the tensors shifted to the center of reaction
        cor_info = calc_TLS_center_of_reaction(T, L, S, centroid)

        cor    = cor_info["COR"]
        T_cor  = cor_info["T'"]
        L_cor  = cor_info["L'"]
        S_cor  = cor_info["S'"]
        T_red  = cor_info["rT'"]

        ## return information
        fit_info["cor_x"] = cor[0]
        fit_info["cor_y"] = cor[1]
        fit_info["cor_z"] = cor[2]

        fit_info["t11"] = T_cor[0,0]
        fit_info["t22"] = T_cor[1,1]
        fit_info["t33"] = T_cor[2,2]
        fit_info["t12"] = T_cor[0,1]
        fit_info["t13"] = T_cor[0,2]
        fit_info["t23"] = T_cor[1,2]

        fit_info["l11"] = L_cor[0,0]
        fit_info["l22"] = L_cor[1,1]
        fit_info["l33"] = L_cor[2,2]
        fit_info["l12"] = L_cor[0,1]
        fit_info["l13"] = L_cor[0,2]
        fit_info["l23"] = L_cor[1,2]

        fit_info["s2211"] = S_cor[1,1] - S_cor[0,0]
        fit_info["s1133"] = S_cor[0,0] - S_cor[2,2]
        fit_info["s12"]   = S_cor[0,1]
        fit_info["s13"]   = S_cor[0,2]
        fit_info["s23"]   = S_cor[1,2]
        fit_info["s21"]   = S_cor[1,0]
        fit_info["s31"]   = S_cor[2,0]
        fit_info["s32"]   = S_cor[2,1]


        if fdict["info"] not in [1, 2, 3]:
            errx = "info = %d" % (fdict["info"])
            fit_info["error"] = errx

        elif min(eigenvalues(T_red))<=TSMALL:
            errx = "Invalid Tr Eigenvalue"
            fit_info["error"] = errx

        elif min(eigenvalues(L_cor))<0.0:
            errx = "Invalid L Eigenvalue"
            fit_info["error"] = errx

        if fit_info.has_key("error"):
            print "[%s] lsq_fit_segment(frag_id={%s..%s}, num_atoms=%d, lsqr=%6.4f, discard=%s)" % (
                self.name, frag_id1, frag_id2, fit_info["num_atoms"], fit_info["lsq_residual"], fit_info["error"]) 
        else:
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
    if tls_model=="HYBRID":
        if USE_LINEARTLS==True:
            return TLSGraphChainFastHybrid()
        else:
            return TLSGraphChainHybrid()

    elif tls_model=="ANISO":
        return TLSGraphChainAnisotropic()

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

