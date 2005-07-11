## TLS Minimized Domains (TLSMD)
## Copyright 2002 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

###############################################################################
## Grid/TLSGraph Calculation System
##
import sys
import copy
import string

from misc     import *
from datafile import TLSMDFile
from hcsssp   import HCSSSP

## spiffy new C tlsmdmodule for MINPACK minimization
#sys.path.append("/home/jpaint/tlsmd/src")
#import tlsmdmodule
USE_TLSMDMODULE = False

###############################################################################
## Hacked Global Options
##

## ANISO/ISO TLS Model: set externally!
TLS_MODEL = "HYBRID"

## LSQ weighting
WEIGHT_MODEL = "UNIT"

## atom selection options
INCLUDE_ATOMS = "ALL" 

## minimum span of residues for TLS subsegments
MIN_SUBSEGMENT_SIZE = 3

## use Uiso residual
USE_UISO_RESIDUAL = True

def print_globals():
    print "TLSMD GLOBAL OPTIONS"
    print "    TLS_MODEL ==================: %s" % (TLS_MODEL)
    print "    MIN_SUBSEGMENT_SIZE --------: %d" % (MIN_SUBSEGMENT_SIZE)
    print "    USE_UISO_RESIDUAL ==========: %s" % (USE_UISO_RESIDUAL)
    print "    WEIGHT_MODEL ===============: %s" % (WEIGHT_MODEL)
    print "    INCLUDE ATOMS --------------: %s" % (INCLUDE_ATOMS)
    print

def set_globals():
    """Takes settings from misc.GLOBALS dictionary and applies them to
    the globals in the current module.
    """
    global TLS_MODEL
    global WEIGHT_MODEL
    global INCLUDE_ATOMS
    global MIN_SUBSEGMENT_SIZE

    TLS_MODEL     = GLOBALS["TLS_MODEL"]
    WEIGHT_MODEL  = GLOBALS["WEIGHT_MODEL"]
    INCLUDE_ATOMS = GLOBALS["INCLUDE_ATOMS"]

    assert TLS_MODEL     in ["ANISO", "HYBRID"]
    assert WEIGHT_MODEL  in ["UNIT", "IUISO"]
    assert INCLUDE_ATOMS in ["ALL", "MAINCHAIN", "CA"]

    if INCLUDE_ATOMS=="ALL":
        MIN_SUBSEGMENT_SIZE = 3
    elif INCLUDE_ATOMS=="MAINCHAIN":
        MIN_SUBSEGMENT_SIZE = 7
    elif INCLUDE_ATOMS=="CA":
        MIN_SUBSEGMENT_SIZE = 20

    GLOBALS["MIN_SUBSEGMENT_SIZE"] = MIN_SUBSEGMENT_SIZE

###############################################################################
## Atom Selection/Weighting Functions
##

def calc_include_atom(atm):
    """Filter out atoms from the model which will cause problems or
    cont contribute to the TLS analysis.
    """
    if trace(atm.get_U())<=TSMALL:
        return False

    if INCLUDE_ATOMS=="ALL":
        return True

    elif INCLUDE_ATOMS=="MAINCHAIN":
        if atm.name in ["N", "CA", "C"]:
            return True
        else:
            return False

    elif INCLUDE_ATOMS=="CA":
        if atm.name=="CA":
            return True
        else:
            return False
        
    return True

def calc_atom_weight(atm):
    """Weight the least-squares fit according to this function.
    """
    weight = atm.occupancy

    if WEIGHT_MODEL=="IUISO":
        weight = weight * (1.0 / (B2U * atm.temp_factor))

    return weight

def iter_ij(num_vertex, min_span):
    """Iterates over the i,j vertex indexes defining the edges
    to be built for the graph.  num_vertex gives the number
    of consecutive vertixes to create, min_span is the minimum
    number of residues (fragments) a edge should span.
    """
    for i in range(num_vertex):
        for j in range(i+min_span, num_vertex):
            yield i, j

def iter_chain_subsegment_descs(chain, min_span):
    """Iterate over all possible subsegments of the given Chain object
    with a minimum size of min_span fragments.  The segments are yielded
    as Python dictionaries containing a description of the subsegment.
    """
    num_vertex = len(chain) + 1
    
    for vertex_i, vertex_j in iter_ij(num_vertex, min_span):
        ## begin message
        msg                 = {}
        msg["chain"]        = chain
        msg["chain_id"]     = chain.chain_id

        ## vertex span of edge info
        msg["vertex_i"]     = vertex_i
        msg["vertex_j"]     = vertex_j

        ## the definition of the subsegment
        msg["frag_id1"]     = chain[vertex_i].fragment_id
        msg["frag_id2"]     = chain[vertex_j-1].fragment_id

        yield msg


###############################################################################
## Chain Graphing Classes (also can be run as XMLRPC servers)
##

class TLSGraphChain(object):
    """Stores a protein chain of atoms in the tuple list xmlrpc_chain along
    with the residue(fragment) indexes and U ADP tensor of each atom.
    Subsegments of the protein chain then can be fit with TLS tensors.
    """
    pass


class TLSGraphChainHybrid(TLSGraphChain):
    """Graph the chain with a hybrid isotropic/anisotropic TLS model. 
    """
    def __init__(self):
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
            
            atm_desc["position"] = array(
                (atm_desc["x"], atm_desc["y"], atm_desc["z"]), Float)
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
                
	if iend==None:
	    iend = len(self.f) - 1

        ## are there enough atoms in this chain segment
        num_atoms = iend - istart + 1
        fit_info["num_atoms"] = num_atoms
        if num_atoms<13:
            fit_info["error"] = "%d Atoms In Segment" % (num_atoms)
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

            ## uiso
            u_iso = atm_desc["u_iso"]

            ## set the A Matrix, B vector
            set_TLSiso_A(A_ISOW, i, 0, x, y, z, w)
            set_TLSiso_b(B_ISOW, i, u_iso, w)

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
        U_ISOW  = matrixmultiply(A_ANISOW, X_ANISO)

        T11, T22, T33, T12, T13, T23, L11, L22, L33, L12, L13, L23, \
             S1133, S2211, S12, S13, S23, S21, S31, S32 = (
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
            print "[HYBRID] lsq_fit_segment("\
                  "frag_id={%s..%s}, "\
                  "num_atoms=%d, lsqr=%6.4f, discard=%s)" % (
                frag_id1, frag_id2,
                fit_info["num_atoms"], fit_info["lsq_residual"],
                fit_info["error"]) 
        else:
            print "[HYBRID] lsq_fit_segment("\
                  "frag_id={%s..%s}, "\
                  "num_atoms=%d, lsqr=%6.4f)" % (
                frag_id1, frag_id2,
                fit_info["num_atoms"], fit_info["lsq_residual"])

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

        print "[ANISO] set_xmlrpc_chain(A=%s, b=%s, f=%d)" % (
            shape(self.A), shape(self.b), len(f))

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
            w  = math.sqrt(atm_desc["w"])

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
                      atm_desc["u12"], atm_desc["u13"], atm_desc["u23"], w)

            ## set the A matrix
            set_TLS_A(A, iU11, 0, x, y, z, w)

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
        if num_atoms<5:
            fit_info["error"] = "%d Atoms In Segment" % (num_atoms)
            return fit_info

        ## LSQ Fit
        X = solve_TLS_Ab(Aw, Bw)

        ## calculate the weighted and unweighted TLS-predicted Uij values
        Uw = matrixmultiply(Aw, X)

        ## calculate the lsq residual since silly Numeric Python won't
        ## do it for us
        Dw = Uw - Bw

        if USE_UISO_RESIDUAL:
            ## this residual attempts not to penalize highly anisotropic
            ## TLS groups by returning a residual aginst Uiso/Utlsiso
            u_iso_residual = 0.0
            
            for i in range(num_atoms):
                iU11 = i * 6

                u_iso     = (Bw[iU11] + Bw[iU11+1] + Bw[iU11+2])/3.0
                u_iso_tls = (Uw[iU11] + Uw[iU11+1] + Uw[iU11+2])/3.0

                u_iso_residual += (u_iso_tls - u_iso)**2
                
            fit_info["lsq_residual"] = u_iso_residual
            
        else:
            fit_info["lsq_residual"] = dot(Dw, Dw)

	## shift TLS tensors to the center of reaction
        T11, T22, T33, T12, T13, T23, L11, L22, L33, L12, L13, L23, \
             S1133, S2211, S12, S13, S23, S21, S31, S32 = (
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
            print "[ANISO] lsq_fit_segment("\
                  "frag_id={%s..%s}, "\
                  "num_atoms=%d, lsqr=%6.4f, discard=%s)" % (
                frag_id1, frag_id2,
                fit_info["num_atoms"], fit_info["lsq_residual"],
                fit_info["error"]) 
        else:
            print "[ANISO] lsq_fit_segment("\
                  "frag_id={%s..%s}, "\
                  "num_atoms=%d, lsqr=%6.4f)" % (
                frag_id1, frag_id2,
                fit_info["num_atoms"], fit_info["lsq_residual"]) 

        return fit_info


class PluginAtom(object):
    def __init__(self, position, temp_factor, U):
        self.position = position
        self.temp_factor = temp_factor
        self.U = U


class TLSGraphChainPlugin(TLSGraphChain):
    """Graph the chain with a hybrid isotropic/anisotropic TLS model. 
    """
    def __init__(self):
        self.patom_list  = None
        self.weight_dict = None
        self.f           = None
    
    def set_xmlrpc_chain(self, xmlrpc_chain):
        """Sets a new TLS chain. which generates new self.A, self.b,
        and self.f arrays.
        """
        ## create a Numeric Array vector for each atom position (faster)
        patom_list   = []
        weight_dict  = {}
        frag_id_list = []
        
        for atm_desc in xmlrpc_chain:
            patm = PluginAtom(
                array((atm_desc["x"], atm_desc["y"], atm_desc["z"]), Float),
                U2B * atm_desc["u_iso"],
                array(((atm_desc["u11"], atm_desc["u12"], atm_desc["u13"]),
                       (atm_desc["u12"], atm_desc["u22"], atm_desc["u23"]),
                       (atm_desc["u13"], atm_desc["u23"], atm_desc["u33"])),
                      Float))
            patom_list.append(patm)
            
            weight_dict[patm] = atm_desc["w"]
            frag_id_list.append(atm_desc["frag_id"])
            

        ## set the class frag list
        self.patom_list  = patom_list
        self.weight_dict = weight_dict
        self.f           = frag_id_list

        print "[PLUGIN] set_xmlrpc_chain(num_atoms=%d)" % (len(patom_list))
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
                
	if iend==None:
	    iend = len(self.f) - 1

        ## are there enough atoms in this chain segment
        num_atoms = iend - istart + 1
        fit_info["num_atoms"] = num_atoms
        if num_atoms<13:
            fit_info["error"] = "%d Atoms In Segment" % (num_atoms)
            return fit_info

        patom_list = self.patom_list[istart:iend+1]

        ## CALCULATE CENTROID
        centroid = zeros(3, Float)
        for patm in patom_list:
            centroid += patm.position
        centroid /= num_atoms

        ## SOLVE MODEL
        rdict = calc_TLS_least_squares_fit_for_iso(
            patom_list, centroid, self.weight_dict)

        fit_info["lsq_residual"] = rdict["lsq_residual"]
        
        ## caculate the tensors shifted to the center of reaction
        T = rdict["T"]
        L = rdict["L"]
        S = rdict["S"]
        
        cor_info = calc_TLS_center_of_reaction(T, L, S, centroid)

        cor    = cor_info["COR"]
        T_cor  = cor_info["T'"]
        L_cor  = cor_info["L'"]
        S_cor  = cor_info["S'"]
        T_red  = cor_info["rT'"]

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

        if fit_info.has_key("error"):
            print "[PLUGIN] lsq_fit_segment("\
                  "frag_id={%s..%s}, "\
                  "num_atoms=%d, lsqr=%6.4f, error=%s)" % (
                frag_id1, frag_id2,
                fit_info["num_atoms"], fit_info["lsq_residual"],
                fit_info["error"]) 
        else:
            print "[PLUGIN] lsq_fit_segment("\
                  "frag_id={%s..%s}, "\
                  "num_atoms=%d, lsqr=%6.4f)" % (
                frag_id1, frag_id2,
                fit_info["num_atoms"], fit_info["lsq_residual"]) 

        return fit_info


def NewTLSGraphChain0(tls_model):
    """Generate and return the proper TLSGraphChain subclass for the
    requested TLS model.
    """
    if tls_model=="HYBRID":
        return TLSGraphChainHybrid()
    if tls_model=="ANISO":
        return TLSGraphChainAnisotropic()
    if tls_model=="PLUGIN":
        return TLSGraphChainPlugin()
    raise Exception()


def NewTLSGraphChain():
    """Return the proper TLSGraphChain subclass for the default tls_model.
    """
    return NewTLSGraphChain0(TLS_MODEL)


class TLSGraphChainXMLRPCServer(TLSGraphChain):
    """Runs this object as a xmlrpc server servicing requests for
    functions: set_xmlrpc_chain, and lsq_fit_segment.
    """
    def __init__(self):
        self.proxy = None

    def set_tls_model(self, tls_model, weight_model):
        print "TLSGraphChainXMLRPCServer.set_tls_model(%s, %s)" % (
            tls_model, weight_model)

        ## select proper TLS model
        self.proxy = NewTLSGraphChain0(tls_model)

        ## set weighting scheme
        global WEIGHT_MODEL
        WEIGHT_MODEL = weight_model
        
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
        
        xmlrpc_server.register_function(
            self.set_tls_model, "set_tls_model")

        xmlrpc_server.register_function(
            self.set_xmlrpc_chain, "set_xmlrpc_chain")
        
        xmlrpc_server.register_function(
            self.lsq_fit_segment, "lsq_fit_segment")
        
        xmlrpc_server.serve_forever()


###############################################################################
## XMLRPC Client Threads which dispach jobs to the chain graphers
##

class TLSGridClientThread(Thread):
    """One instance of this class connects to a single xmlrpc server
    on the compute grid.controlls a single compute
    server on the grid.  It uses blocking xmlrpc calls within its own
    thread.
    """
    def __init__(self, server_url, request_queue, result_queue):
        Thread.__init__(self)
        self.setDaemon(True)

        ## set this to True to monitor xmlrpc messages
        self.debug         = False
        
        self.server_url    = server_url
        self.request_queue = request_queue
        self.result_queue  = result_queue

    def get_ident(self):
        if self.server_url==None:
            return "THREAD %d" % (thread.get_ident())
        return "%s THREAD %d" % (self.server_url, thread.get_ident())

    def run(self):
        """The main() for this thread. 
        """
        ## use either a in-process server, or remote
        ## server over xmlrpc
        if self.server_url!=None:
            self.compute_server = self.xmlrpc_connect()
        else:
            self.compute_server = NewTLSGraphChain()

        self.run_grid_client()

    def xmlrpc_connect(self):
        """Connect to the xmlrpc server and return a connection object.
        """
        server = xmlrpclib.ServerProxy(self.server_url)
        server.set_tls_model(TLS_MODEL, WEIGHT_MODEL)
        return server

    def run_grid_client(self):
        """Enters a loop to marshal the request/response messages
        between this server client thread and the controller thread.
        """
        cur_chain_id = None
        
        while True:
            ## get a message requesting the LSQ fit of a TLS
            ## Chain segment (Python dictionary)
            msg = self.request_queue.get()

            ## check the xmlrpc_chain_id of the request to see if
            ## the TLSGridServer needs to be sent the working data
            ## for the chain (if it is a new chain)
            if msg["chain_id"]!=cur_chain_id:
                self.call_set_xmlrpc_chain(msg)
                cur_chain_id = msg["chain_id"]

            ## perform the LSQ fit
            msg["server_url"] = self.server_url
            return_msg = self.call_lsq_fit_segment(msg)

            ## place result on the out_queue
            self.result_queue.put(return_msg)

    def call_set_xmlrpc_chain(self, msg):
        """Calls the xmlrpc server's set_xmlrpc_chain() method.
        """
        if self.debug==True:
            print "TLSGridServerThread.call_set_xmlrpc_chain(%s)" % (
                self.get_ident())

        xmlrpc_chain = msg["xmlrpc_chain"]
        self.compute_server.set_xmlrpc_chain(xmlrpc_chain)

    def call_lsq_fit_segment(self, msg):
        """Calls the xmlrpc server's lsq_fit_segment() method.
        """
        if self.debug==True:
            print "TLSGridServerThread.call_lsq_fit_segment(%s)" % (
                self.get_ident())
        
	frag_id1 = msg["frag_id1"]
	frag_id2 = msg["frag_id2"]

        fit_info = self.compute_server.lsq_fit_segment(
            frag_id1, frag_id2)

        ## merge the fit_info data with the original message
        for key in fit_info.keys():
            if msg.has_key(key):
                assert msg[key]==fit_info[key]
            msg[key] = fit_info[key]

        return msg


class TLSGridServerPool(object):
    """Launches a pool of threads using the TLSGridClientThread() class as
    the controlling object for each of these threads.  Each thread in the
    pool waits for jobs placed in this object's self.requeust_queue, and
    the results are placed by these threads back into self.result_queue.
    """
    def __init__(self):
        self.client_thread_list = []
        self.request_queue      = Queue(1)
        self.result_queue       = Queue(1)

    def launch_client_thread(self, server_url):
        """Launches one TLSGridServerThread.
        """
        client_thread = TLSGridClientThread(
            server_url, self.request_queue, self.result_queue)
        
        self.client_thread_list.append(client_thread)
        client_thread.start()

    def iter_segment_processor(self, iter_segment_producer):
        """Consume the seg_info dictionaries produced by the seg_producer
        iterator, and run them concurrently on all avilible
        TLSGridServerThread threads.
        Yield back the calculation results in the form of a fit_info
        dict for each calcuation request.
        """
        open_requests = 0

        ## make all requests
        for seg_info in iter_segment_producer:
            self.request_queue.put(seg_info)
            open_requests += 1

            while not self.result_queue.empty():
                yield self.result_queue.get()
                open_requests -= 1
        
        ## wait for any remaining requests
        while open_requests>0:
            yield self.result_queue.get()
            open_requests -= 1


class TLSChainProcessor(object):
    """Uses a TLSGridServerPool object to fit all possible TLS subsegment
    to a given Chain object.  The resulting TLS groups are stored in the
    analysis.datafile.grh_*() file.
    """
    def __init__(self, server_pool):
        self.server_pool = server_pool

    def convert_chain_for_xmlrpc(self, chain):
        """Converts the Atoms of a Chain/Segment to a list of dictionaries 
        for transfer over xmlrpc.  Only the information required to fit
        TLS groups to the Atoms is set in the Atom dictionaries to
        reduce traffic over the xmlrpc calls.
        """
        print "convert_chain_for_xmlrpc(chain_id=%s)" % (chain.chain_id)
        
        xmlrpc_chain = []

        for atm in chain.iter_all_atoms():

            if not calc_include_atom(atm):
                print "calc_include_atom(%s)==False" % (atm)
                continue

            atm_desc = {}
            xmlrpc_chain.append(atm_desc)

            atm_desc["name"]    = atm.name
            atm_desc["frag_id"] = atm.fragment_id

            atm_desc["x"] = atm.position[0]
            atm_desc["y"] = atm.position[1]
            atm_desc["z"] = atm.position[2]

            atm_desc["u_iso"] = B2U * atm.temp_factor

            U = atm.get_U()
            atm_desc["u11"] = U[0,0]
            atm_desc["u22"] = U[1,1]
            atm_desc["u33"] = U[2,2]
            atm_desc["u12"] = U[0,1]
            atm_desc["u13"] = U[0,2]
            atm_desc["u23"] = U[1,2]

            ## calculate weight
            atm_desc["w"] = calc_atom_weight(atm)

        return xmlrpc_chain

    def iter_lsq_fit_messages(self, analysis, chain, min_span, xmlrpc_chain):
        """Iterate over (and create) the message dictionaries sent
        to consumer threads which are the the graph edge definitions
        of the TLS segments we need LSQ TLS fit data on.
        """
        for msg in iter_chain_subsegment_descs(chain, min_span):
            if analysis.tlsmdfile.grh_get_tls_record(
                msg["chain_id"], msg["frag_id1"], msg["frag_id2"])!=None:
                continue

            msg["method"] = "TLS"
            msg["xmlrpc_chain"] = xmlrpc_chain
            yield msg
        
    def process_chain(self, analysis, chain, min_span):
        """Fits TLS groups to all possible subsegments of the given Chain,
        storing the results in the analysis.grh_*() file.  
        """
        print "PROCESSING CHAIN: ",chain.chain_id

        ## convert the Chain object to a list of atoms
        xmlrpc_chain = self.convert_chain_for_xmlrpc(chain)

        ## create a iterator for the protein chain which yields
        ## a series of messages to call lsq_tls_fit on all the
        ## segments
        ##
        ## XXX: come up with a unique ID for the structure which
        ##      is included in the messages, so the return messages
        ##      from the pool can be checked
        iter_msg = self.iter_lsq_fit_messages(
            analysis, chain, min_span, xmlrpc_chain)

        iter_segment_fit_info = self.server_pool.iter_segment_processor(
            iter_msg)

        ## iterate over the fit_info dictionaries which hold the
        ## information on the TLS fits comming back from the grid servers
        num_edges          = 0
        num_rejected_edges = 0
        num_accepted_edges = 0
        
        for fit_info in iter_segment_fit_info:
            assert fit_info["chain_id"]==chain.chain_id

            num_edges += 1

            ## remove the reference to the xmlrpc_chain -- it's huge
            del fit_info["xmlrpc_chain"]

            ## save the fit_info record the graph state file
            analysis.tlsmdfile.grh_append_tls_record(fit_info)
            
            if fit_info.has_key("error"):
                num_rejected_edges += 1
                print "process_chain(chain_id=%s, frag_id={%s..%s}, "\
                      "discard=%s)" % (
                    chain.chain_id,
                    fit_info["frag_id1"], fit_info["frag_id2"],
                    fit_info["error"])
    
            else:
                num_accepted_edges += 1
                print "process_chain(chain_id=%s, frag_id={%s..%s}, "\
                      "lsqr=%6.4f)" % (
                    chain.chain_id,
                    fit_info["frag_id1"], fit_info["frag_id2"],
                    fit_info["lsq_residual"])

        print
        print "NUM POSSIBLE  TLS GROUPS: %d" % (num_edges)
        print "NUM ACCEPTED  TLS GROUPS: %d" % (num_accepted_edges)
        print "NUM DISCARDED TLS GROUPS: %d" % (num_rejected_edges)


class TLSOptimization(object):
    """Collection object containing one TLS description of a protein chain.
    """
    def __init__(self, chain, ntls_constraint):
        self.chain           = chain
        self.ntls_constraint = ntls_constraint
        self.ntls            = ntls_constraint
        self.tls_list        = []
        self.residual        = 0.0

    def is_valid(self):
        """Return True if the optimization is valid; otherwise, return False.
        """
        return len(self.tls_list)>0

    def add_tls_record(self, tls):
        """Adds a tls informatio dictionary.
        """
        self.tls_list.append(tls)
        self.residual += tls["lsq_residual"]


class TLSChainMinimizer(HCSSSP):
    """Finds the minimal TLS description of a given Chain object using
    the HCSSP global optimization algorithm and a constraint on the number
    of TLS which can be used for the minimization.
    """
    def __init__(self, analysis, chain, min_span):
        self.analysis   = analysis
        self.min_span   = min_span

        self.minimized = False
        self.D         = None
        self.P         = None
        self.T         = None

        ## source and destination fragment IDs; these are the
        ## begining and ending fragments(residues) of the chain
        self.frag_id_src  = None
        self.frag_id_dest = None

        print "TLSChainMinimizer(chain_id=%s, range={%s..%s})" % (
            chain.chain_id, self.frag_id_src, self.frag_id_dest)

        if self.frag_id_src==None and self.frag_id_dest==None:
            self.chain = chain
        else:
            self.chain = chain[self.frag_id_src:self.frag_id_dest]
        
        ## Initalize the vertex list based on the number of of
        ## residues in Chain.
        self.num_vertex = len(self.chain) + 1

        ## calculate the minimum temperature factor in the chain
        ## XXX: hack
        self.min_temp_factor = None
        for atm in self.chain.iter_atoms():
            if self.min_temp_factor==None:
                self.min_temp_factor = atm.temp_factor
                continue
            self.min_temp_factor = min(self.min_temp_factor, atm.temp_factor)
        
    def run_minimization(self, max_tls_segments=20):
        """Run the HCSSSP minimization on the self.V,self.E graph, resulting
        in the creation of the self.D, self.P, and self.T arrays which
        contain 
        """
        ## build the vertex labels to reflect the protein structure
        ## the graph spans
        V = []
        for i in range(self.num_vertex):
            
            ## add the vertex label for i at Vi
            if i==0:
                vertex_label = "N-TERM"

            elif i==self.num_vertex - 1:
                vertex_label = "C-TERM"

            else:
                vertex_label = "%s{%s:%s}" % (
                    self.chain.chain_id,
                    self.chain[i-1].fragment_id,
                    self.chain[i].fragment_id)

            vertex_label = "V%d[%s]" % (i, vertex_label)
            V.append(vertex_label)

        self.V = V

        ## now build edges for the graph with weights given by the LSQ
        ## residual of TLS group fits
        grh_get_tls_record = self.analysis.tlsmdfile.grh_get_tls_record
        
        E = []
        for msg in iter_chain_subsegment_descs(self.chain, self.min_span):
            tls = grh_get_tls_record(
                self.chain.chain_id,
                msg["frag_id1"],
                msg["frag_id2"])
            
            if tls==None:
                print "build_graph(): [ERROR] no TLS group %s{%s..%s}" % (
                    self.chain.chain_id,
                    msg["frag_id1"],
                    msg["frag_id2"])

                sys.exit(-1)

            ## filter out the bad TLS segments
            if self.__minimization_filter(tls)==False:
                continue

            weight = tls["lsq_residual"]
            frag_range = (tls["frag_id1"], tls["frag_id2"])
                
            edge = (msg["vertex_i"], msg["vertex_j"], weight, frag_range)
            E.append(edge)

        self.E = E


        ## fill in any un-reachable gaps in the structure by adding
        ## fake 0.0 cost edges where needed
        
        

        ## perform the minimization
        start_timing()

        if len(self.E)>0:
            print "run_minimization(chain_id=%s): HCSSSP Minimizing..." % (
                self.chain.chain_id)
        
            D, P, T = self.HCSSSP_minimize(self.V, self.E, max_tls_segments)

            self.minimized = True
            self.D = D
            self.P = P
            self.T = T
        else:
            print "run_minimization(chain_id=%s): Unable to minimize" % (
                self.chain.chain_id)
            self.minimized = False

        ## free memory taken up from edges
        E      = None
        self.E = None
        import gc
        gc.collect()

        print "run_minimization(): ",end_timing()

    def __minimization_filter(self, tls):
        """Returns False if the tls group descibed by the tls database
        record should not be included in the minimization.
        """
        if tls.has_key("error"):
            return False
        
        T,L,S,O = self.__TLSO(tls)
        cdict = calc_TLS_center_of_reaction(T, L, S, O)

        ## sanity checks on rT
        rT = cdict["rT'"]
        evals = eigenvalues(rT)
        
        min_rT = U2B * min(evals)
        if min_rT < (0.1 * self.min_temp_factor):
            print "[FILTER] small rT chain=%s frag_rng={%s..%s}" % (
                tls["chain_id"], tls["frag_id1"], tls["frag_id2"])
            return False

        aniso = min(evals)/max(evals)
        if aniso<0.05:
            print "[FILTER] anisotropic rT chain=%s frag_rng={%s..%s}" % (
                tls["chain_id"], tls["frag_id1"], tls["frag_id2"])
            return False

        ## sanity checks on rho
        max_len = 30.0
        if length(cdict["L1_rho"])>max_len or \
           length(cdict["L2_rho"])>max_len or \
           length(cdict["L3_rho"])>max_len:
            print "[FILTER] length(rho) chain=%s frag_rng={%s..%s}" % (
                tls["chain_id"], tls["frag_id1"], tls["frag_id2"])
            return False
           
        return True

    def calc_tls_optimization(self, ntls_constraint):
        """Return a TLSOptimization() object containing the optimal
        TLS description of self.chain using num_tls_segments.
        """
        if not self.minimized:
            return None

        tlsopt = TLSOptimization(self.chain, ntls_constraint)
        
        for hi, hj, edge in self.HCSSSP_path_iter(
            self.V, self.D, self.P, self.T, ntls_constraint):

            if edge==None:
                continue
            
            ## check if the edge is a bypass-edge type
            i, j, weight, frag_range = edge
            
            if len(frag_range)==2:
                tls = self.__calc_tls_record_from_edge(edge)
                tlsopt.add_tls_record(tls)

        return tlsopt

    def __calc_tls_record_from_edge(self, edge):
        """Independently calculate the TLS parameters for the segment
        to verify correctness (internal check).
        """
        i, j, weight, frag_range = edge

        ## retrieve from cache if possible
        frag_id1, frag_id2 = frag_range

        tls = self.analysis.tlsmdfile.grh_get_tls_record(
            self.chain.chain_id, frag_id1, frag_id2)

        frag1    = self.chain[frag_id1]
        frag2    = self.chain[frag_id2]

        segment  = self.chain[frag_id1:frag_id2]

        print "calc_tls_record_from_edge("\
              "chain_id=%s frag_id={%s..%s})" % (
            self.chain.chain_id, frag_id1, frag_id2)

        ## create TLSGroup
        tls_group = TLSGroup()

        ## add atoms to the group
        for atm in segment.iter_atoms():
            if calc_include_atom(atm):
                tls_group.append(atm)

        ## take the TLS group tensors from the database set
        ## the TLSGroup object with them
        self.__add_tensors(tls_group, tls)

        ## helpful additions
        tls_info                    = tls_group.calc_tls_info()
        tls["tls_group"]            = tls_group
        tls["tls_info"]             = tls_info
        tls["segment"]              = segment
        tls["lsq_residual_per_res"] = tls["lsq_residual"] / (len(segment))

        return tls

    def __add_tensors(self, tls_group, tls):
        """Adds the TLS tensors from the tls database dictionary to the
        tls_group.
        """
        T,L,S,O = self.__TLSO(tls)
        tls_group.T      = T
        tls_group.L      = L
        tls_group.S      = S
        tls_group.origin = O
        
    def __TLSO(self, tls):
        """Builds Numeric Array tensors from the tls database dictionary.
        """
        O = array(
            [tls["cor_x"], tls["cor_y"], tls["cor_z"]], Float)

        T = array(
            [ [tls["t11"], tls["t12"], tls["t13"]],
              [tls["t12"], tls["t22"], tls["t23"]],
              [tls["t13"], tls["t23"], tls["t33"]] ], Float)
        
        L = array(
            [ [tls["l11"], tls["l12"], tls["l13"]],
              [tls["l12"], tls["l22"], tls["l23"]],
              [tls["l13"], tls["l23"], tls["l33"]] ], Float)
        
        s11, s22, s33 = calc_s11_s22_s33(tls["s2211"], tls["s1133"]) 
        
        S = array(
            [ [       s11, tls["s12"], tls["s13"]],
              [tls["s21"],        s22, tls["s23"]],
              [tls["s31"], tls["s32"],       s33] ], Float)

        return T,L,S,O
        
    def __add_bypass_edges(self):
        """Modify the graph before minimization to include the possibility
        of short segments of the protein which are not described well
        by the TLS model.
        """
        for i in range(self.num_vertex):
            for j in range(i, self.num_vertex):


                tls = {"method": "BYPASS"}

                weight = (j-i) * 100.0
                edge = (i, j, weight, tls)
                self.E.append(edge)
        
    def __add_bypass_edges_old(self):
        """Modify the graph before minimization to include the possibility
        of short segments of the protein which are not described well
        by the TLS model.
        """
        ## compute the mean weight/residue for all edges
        nres = 0
        mean = 0.0
        edge_mean_list = []

        for i,j,weight,tls in self.E:
            ## weight/residue for this edge
            edge_mean = weight / (j-i)
            edge_mean_list.append(edge_mean)

            ## for overall mean weight/residue
            nres += j-i
            mean += weight

        mean = mean / float(nres)

        ## write out histogram file
        num_bins = 100
        bins     = [0 for i in range(num_bins)]

        min_mean = min(edge_mean_list)
        max_mean = max(edge_mean_list)

        bin_range = max_mean - min_mean
        bin_width = bin_range / num_bins

        ## calculate the histogram if it's worth it to do so
        if bin_width>0.0:
            for edge_mean in edge_mean_list:
                bin_i = int((edge_mean - min_mean) / bin_width)

                try:
                    bins[bin_i] += 1
                except IndexError:
                    pass

            ## bin calculations for the mean value
            mean_bin_i = int( (mean - min_mean) / bin_width ) 

            hist = open("mean_weight_per_residue_histogram.dat", "w")

            hist.write("## min mean weight/residue:    %f\n" % (min_mean))
            hist.write("## max mean weight/residue:    %f\n" % (max_mean))
            hist.write("## mean weight/residue:        %f\n" % (mean))
            hist.write("## bin of mean weight/residue: %d\n" % (mean_bin_i)) 

            for i in range(len(bins)):
                hist.write("%10d   %10d\n" % (i+1, bins[i]))

            hist.close()

        ## create bypass edges with a weight equal to
        ## the average
        max_bypass_len = 20
    
        for i in range(self.num_vertex - max_bypass_len):
            for j in range(1, max_bypass_len+1):

                tls = {"method": "bypass"}

                weight = j * 100.0
                edge = (i, i+j, weight, tls)
                self.E.append(edge)
        
    def prnt_detailed_paths(self, hops=20):
        """Debug
        """
        if not self.minimized:
            return
        dest_j = len(self.V)-1

        for h in range(1, hops+1):
            print
            print "MINIMIZATON VERTEX PATH FOR %d SEGMENTS" % (h)
            print "NODE LABEL              HOPS      COST"\
                  "      PREVIOUS NODE          EDGE"

            self.__detailed_path(self.V, self.D, self.P, self.T, h)

    def __detailed_path(self, V, D, P, T, hop_constraint):
        """Print out the path from the source vertex (vertex 0) to
        the destination vertex (end vertex) given the hop_constarint.
        """
        num_vertex = len(D[0])
        
        ## start at the destination vertex
        curr_v = num_vertex - 1
        h      = hop_constraint
        
        while curr_v>=0:
            prev_vertex  = P[h,curr_v]
            vertex_label = V[curr_v].ljust(20)

            if prev_vertex<0:
                prev_vertex_label = "".ljust(20)
            else:
                prev_vertex_label = V[prev_vertex].ljust(20)

            edge = T[h][curr_v]

            if edge!=None:
                i, j, weight, frag_range = edge
                wr = weight / (j - i)
                edge_label = "(%3d,%3d,%6.3f,%s) %6.3f" % (
                    i, j, weight, frag_range, wr)
            else:
                edge_label = ""
                 
            print "%s   %3d     %10.4f   %s   %s" % (
                vertex_label, h, D[h,curr_v], prev_vertex_label,
                edge_label)

            curr_v = prev_vertex
            h -= 1


class TLSMDAnalysis(object):
    """Central object for a whole-structure TLS analysis.
    """
    def __init__(self,
                 struct_path    = None,
                 sel_chain_ids  = None,
                 tlsdb_file     = None,
                 tlsdb_complete = False,
                 gridconf_file  = None,
                 num_threads    = 1):

        set_globals()
        print_globals()

        self.struct_path     = struct_path
        if sel_chain_ids!=None:
            self.sel_chain_ids = sel_chain_ids.split(",")
        else:
            self.sel_chain_ids = None
        self.tlsdb_file      = tlsdb_file
        self.tlsdb_complete  = tlsdb_complete
        self.gridconf_file   = gridconf_file
        self.num_threads     = num_threads

        self.struct          = None
        self.struct_id       = None
        self.chains          = None
        self.tlsmdfile       = None

    def run_optimization(self):
        """Run the TLSMD optimization on the structure.
        """
        self.load_struct()

        ## auto name of tlsdb file then open
        if self.tlsdb_file==None:
            self.tlsdb_file = "%s.db" % (self.struct_id)
        self.tlsmdfile = TLSMDFile(self.tlsdb_file)

        ## select chains for analysis
        self.select_chains()

        ## print these settings
        self.prnt_settings()

        if not self.tlsdb_complete:
            self.calc_tls_segments()
        
        self.calc_chain_minimization()

    def prnt_settings(self):
        chain_ids = []
        for chain in self.chains:
            chain_ids.append(chain.chain_id)
        cids = string.join(chain_ids, ",")
        
        print "TLSMD ANALYSIS SETTINGS"
        print "    STRUCTURE FILE ================: %s"%(self.struct_path)
        print "    STRUCTURE ID ------------------: %s"%(self.struct_id)
        print "    CHAIN IDs SELECTED FOR ANALYSIS: %s"%(cids)
        print "    DATABASE FILE PATH ------------: %s"%(self.tlsdb_file)
        print "    GRID SERVER CONFIG FILE =======: %s"%(self.gridconf_file)
        print
        
    def load_struct(self):
        """Loads Structure, chooses a unique struct_id string.
        """
        print "LOADING STRUCTURE"
        print "    PATH: %s" % (self.struct_path)

        ## load struct
        self.struct = LoadStructure(
            fil = self.struct_path,
            build_properties = ("library_bonds","distance_bonds"))

        ## set the structure ID
        if GLOBALS.has_key("STRUCT_ID"):
            struct_id = GLOBALS["STRUCT_ID"]
        else:
            struct_id = self.struct.structure_id

        if struct_id==None or struct_id=="XXXX":
            directory, filename = os.path.split(self.struct_path)
            basename, ext       = os.path.splitext(filename)
            struct_id           = basename

        self.struct_id = struct_id

        print "    STRUCT ID: %s" % (self.struct_id)
        print

        ## if there are REFMAC5 TLS groups in the REMARK records of
        ## the PDB file, then add those in
        tls_file = TLSFile()
        tls_file.set_file_format(TLSFileFormatPDB())

        fil = open(self.struct_path, "r")
        tls_file.load(fil)

        if len(tls_file.tls_desc_list)>0:
            print "MERGING TLS GROUPS"
            print "    NUM TLS GROUPS: %d" % (len(tls_file.tls_desc_list))

            ## assume REFMAC5 groups where Utotal = Utls + Biso(temp_factor)
            for tls_desc in tls_file.tls_desc_list:
                tls_group = tls_desc.construct_tls_group_with_atoms(
                    self.struct)

                print "    TLS GROUP: %s" % (tls_group.name)
                
                for atm, Utls in tls_group.iter_atm_Utls():
                    bresi = atm.temp_factor
                    atm.temp_factor = bresi + (U2B * trace(Utls) / 3.0)
                    atm.U = (B2U * bresi * identity(3, Float)) + Utls
            
            print

    def select_chains(self):
        """Selects chains for analysis.
        """
        ## select viable chains for TLS analysis
        segments = []
        
        for chain in self.struct.iter_chains():
            ## if self.sel_chain_ids is set, then only use those
            ## selected chain ids
            if self.sel_chain_ids!=None:
                if chain.chain_id not in self.sel_chain_ids:
                    continue

            ## count the number of amino acid residues in the chain
            if chain.count_amino_acids()<MIN_SUBSEGMENT_SIZE:
                continue

            ## ok, use the chain but use a segment and cut off
            ## any leading and trailing non-amino acid residues
            frag_id1 = None
            for aa in chain.iter_amino_acids():
                if frag_id1==None:
                    frag_id1 = aa.fragment_id
                frag_id2 = aa.fragment_id
                
            segment = chain[frag_id1:frag_id2]
            segments.append(segment)
        
        self.chains = segments
        
    def calc_tls_segments(self):
        """Calculates the TLSGraph for each chain in self.chains, optionally
        loading pre-computed graphs from the graph_file.  Any chains
        which need TLSGraphs computed will be strored in the graph_file.
        """
        chain_processor = self.launch_chain_processor()
        
        for chain in self.chains:
            begin_chain_timing(chain.chain_id)

            print "BUILDING TLS DATABASE FOR %s" % (chain)
            chain_processor.process_chain(self, chain, MIN_SUBSEGMENT_SIZE)

            end_chain_timing(chain.chain_id)

    def launch_chain_processor(self):
        """Starts up a in-process or grid server pool for graphing
        this structure.
        """
        ## setup the server pool for local in-process running (possibly
        ## using threads), or using a compute server grid
        server_pool = TLSGridServerPool()

        if self.gridconf_file!=None:
            ## open the grid config file and read the URLs of the compute
            ## server and add them to the server pool
            fil = open(self.gridconf_file, "r")

            ## each line is a comment starting with # or the URL of a
            ## grid server
            for ln in fil.readlines():
                server_url = ln.strip()
                if server_url.startswith("#"):
                    continue
                server_pool.launch_client_thread(server_url)

        else:
            ## launch local processing threads -- this is no good; they
            ## only run on one processor!
            for x in range(self.num_threads):
                server_pool.launch_client_thread(None)

        ## run the job to compute the full TLS Graph using the
        ## server_pool
        chain_processor = TLSChainProcessor(server_pool)
        return chain_processor

    def calc_chain_minimization(self):
        """Performs the TLS graph minimization on all TLSGraphs.
        """
        for chain in self.chains:
            chain.tls_chain_minimizer = TLSChainMinimizer(
                self, chain, MIN_SUBSEGMENT_SIZE)
            
            chain.tls_chain_minimizer.run_minimization()
            if not chain.tls_chain_minimizer.minimized:
                continue

            print
            print "="*79
            print "MINIMIZING CHAIN %s" % (chain)
            chain.tls_chain_minimizer.prnt_detailed_paths()
