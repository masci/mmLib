## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Utility classes for loading, manipulating, and analyzing TLS parameters.
"""
from __future__ import generators
import re
import sys
import string
import fpformat

from mmLib.mmTypes   import *
from mmLib.Structure import *


## regular expressions used for extracting the TLS groups from the REFMAC
## TLS tensor file format
refmac_tls_regex_dict = {
    "group": re.compile(
    "(?:^TLS\s*$)|(?:^TLS\s+(.*)$)"),
    
    "range": re.compile(
    "^RANGE\s+[']([A-Z])\s*([0-9A-Z.]+)\s*[']\s+"\
    "[']([A-Z])\s*([0-9A-Z.]+)\s*[']\s*(\w*).*$"),

    "origin": re.compile(
    "^ORIGIN\s+(\S+)\s+(\S+)\s+(\S+).*$"),
    
    "T": re.compile(
    "^T\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+).*$"),
    
    "L": re.compile(
    "^L\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+).*$"),

    "S": re.compile(
    "^S\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+).*$")
    }


## regular expressions used for extracting the TLS groups from the REMARK
## fields of the PDB file output by REFMAC
refmac_pdb_regex_dict = {
    "group": re.compile(
    "\s*TLS GROUP :\s+(\d+)\s*$"),

    "range": re.compile(
    "\s*RESIDUE RANGE :\s+(\w)\s+(\w+)\s+(\w)\s+(\w+)\s*$"),

    "origin": re.compile(
    "\s*ORIGIN\s+FOR\s+THE\s+GROUP\s+[(]A[)]:\s+(\S+)\s+(\S+)\s+(\S+)\s*$"),

    "t11_t22": re.compile(
    "\s*T11:\s*(\S+)\s+T22:\s*(\S+)\s*$"),

    "t33_t12": re.compile(
    "\s*T33:\s*(\S+)\s+T12:\s*(\S+)\s*$"),

    "t13_t23": re.compile(
    "\s*T13:\s*(\S+)\s+T23:\s*(\S+)\s*$"),

    "l11_l22": re.compile(
    "\s*L11:\s*(\S+)\s+L22:\s*(\S+)\s*$"),

    "l33_l12": re.compile(
    "\s*L33:\s*(\S+)\s+L12:\s*(\S+)\s*$"),

    "l13_l23": re.compile(
    "\s*L13:\s*(\S+)\s+L23:\s*(\S+)\s*$"),

    "s11_s12_s13": re.compile(
    "\s*S11:\s*(\S+)\s+S12:\s*(\S+)\s+S13:\s*(\S+)\s*$"),

    "s21_s22_s23": re.compile(
    "\s*S21:\s*(\S+)\s+S22:\s*(\S+)\s+S23:\s*(\S+)\s*$"),

    "s31_s32_s33": re.compile(
    "\s*S31:\s*(\S+)\s+S32:\s*(\S+)\s+S33:\s*(\S+)\s*$")
    }


class TLSInfo(object):
    """Contains information on one TLS group
    """
    def __init__(self):
        self.name       = ""     # text name (maybe numeric label?)
        self.origin     = None   # (x, y, z)
        self.range_list = []     # (chain1, res1, chain2, res2, selection)
        self.T          = None   # (t11, t22, t33, t12, t13, t23)
        self.L          = None   # (l11, l22, l33, l12, l13, l23)
        self.S          = None   # (s2211, s1133, s12, s13, s23, s21, s31, s32)

    def __str__(self):
        return self.refmac_description()

    def create_tls_name(self):
        """Creates a name for the TLS group using the selected
        residue ranges.
        """
        listx = []
        for (chain_id1, frag_id1, chain_id2, frag_id2, sel) in self.range_list:
            listx.append("%s%s-%s%s %s" % (
                chain_id1, frag_id1, chain_id2, frag_id2, sel))
        return string.join(listx,';')
        
    def refmac_description(self):
        """Return a string describing the TLS group in REFMAC/CCP4 format.
        """
        def ft8(tx):
            strx = ""
            for x in tx:
                strx += fpformat.fix(x, 4).rjust(8)
            return strx
        
        listx = []

        if self.name != "":
            listx.append("TLS %s" % (self.name))
        else:
            listx.append("TLS")

        for (c1, f1, c2, f2, sel) in self.range_list:
            try:
                if f1[-1] in string.digits: f1 += "."
            except IndexError:
                f1 = "."
            try:
                if f2[-1] in string.digits: f2 += "."
            except IndexError:
                f2 = "."
            listx.append("RANGE  '%s%s' '%s%s' %s" % (
                c1, f1.rjust(5), c2, f2.rjust(5), sel))

        if self.origin != None:
            listx.append("ORIGIN %s" % ft8(self.origin))
        if self.T != None:
            listx.append("T   %s" % ft8(self.T))
        if self.L != None:
            listx.append("L   %s" % ft8(self.L))
        if self.S != None:
            listx.append("S   %s" % ft8(self.S))

        return string.join(listx, "\n")

    def make_tls_group(self, struct):
        """Returns a TLSGroup containing the appropriate atoms from the
        argument Structure object, and with the origin, T, L, S tensors
        set.
        """
        tls = TLSGroup()

        if self.name=="":
            tls.name = self.create_tls_name()
        else:
            tls.name = self.name

        tls.set_origin(*self.origin)
        tls.set_T(*self.T)
        tls.set_L(*self.L)
        tls.set_S(*self.S)

        for (chain_id1, frag_id1, chain_id2, frag_id2, sel) in self.range_list:

            chain1 = struct.get_chain(chain_id1)
            if chain1 == None:
                print "[ERROR] no chain id: %s" % (chain_id1)
                sys.exit(1)

            try:
                seg = chain1[frag_id1:frag_id2]
            except KeyError:
                print "[ERROR] unable to find segment: %s-%s" % (
                    frag_id1, frag_id2)
                sys.exit(1)

            for atm in seg.iter_atoms():
                tls.append(atm)
            
        return tls


class TLSInfoList(list):
    """This class reads and writes TLS information stored in the same
    format as REFMAC from CCP4 >= 4.1.0.  The TLS groups are stored as a list
    of TLSInfo classes.
    """
    def __str__(self):
        return self.refmac_description()

    def refmac_description(self):
        """Return a string describing the TLS groups in REFMAC/CCP4 format.
        """
        listx = ["REFMAC"]
        for tls_info in self:
            listx.append(str(tls_info))
        return string.join(listx, "\n")
        
    def load_refmac_tlsout_file(self, fil):
        """Read the TLS information from file object fil, and store that
        inforamtion in the class instance variables.
        """
        tls_info = None
        
        for ln in fil.readlines():
            ln = ln.rstrip()

            ## search all regular expressions for a match
            for (re_key, re_tls) in refmac_tls_regex_dict.items():
                mx = re_tls.match(ln)
                if mx != None:
                    break

            ## no match was found for the line
            if mx == None:
                continue
            
            ## do not allow a match if tls_info == None, because then
            ## somehow we've missed the TLS group begin line
            if tls_info == None and re_key != "group":
                print "[ERROR] TLS group info without TLS group begin line"
                sys.exit(1)

            if re_key == "group":
                tls_info = TLSInfo()
                self.append(tls_info)
                if mx.group(1) != None:
                    tls_info.name = mx.group(1)

            elif re_key == "origin":
                (x, y, z) = mx.groups()
                tls_info.origin = (float(x), float(y), float(z))

            elif re_key == "range":
                (c1, f1, c2, f2, sel) = mx.groups()

                try:
                    if f1[-1] == ".": f1 = f1[:-1]
                except IndexError:
                    pass
                
                try:
                    if f2[-1] == ".": f2 = f2[:-1]
                except IndexError:
                    pass

                tls_info.range_list.append((c1, f1, c2, f2, sel))

            elif re_key == "T":
                ## REFMAC ORDER: t11 t22 t33 t23 t13 t12
                (t11, t22, t33, t23, t13, t12) = mx.groups()
                tls_info.T = (float(t11),
                              float(t22),
                              float(t33),
                              float(t12),
                              float(t13),
                              float(t23))

            elif re_key == "L":
                ## REFMAC ORDER: t11 t22 t33 t23 t13 t12
                (l11, l22, l33, l23, l13, l12) = mx.groups()
                tls_info.L = (float(l11),
                              float(l22),
                              float(l33),
                              float(l12),
                              float(l13),
                              float(l23))

            elif re_key == "S":
                ## REFMAC ORDER: s2211 s1133 s23 s31 s12 s32 s13 s21
                (s2211, s1133, s23, s31, s12, s32, s13, s21) = mx.groups()
                tls_info.S = (float(s2211),
                              float(s1133),
                              float(s12),
                              float(s13),
                              float(s23),
                              float(s21),
                              float(s31),
                              float(s32))

    def process_REMARK(self, rec):
        """Callback for the PDBFile parser for REMARK records.  If the
        REMARK records contain TLS group information, then it is
        extracted and added to the TLSGroups list.
        """
        try:
            if rec["remarkNum"] != 3:
                return
        except KeyError:
            return

        try:
            text = rec["text"]
        except KeyError:
            return

        try:
            tls_info = self[-1]
        except IndexError:
            tls_info = None

        ## attempt to extract information from the text
        for (re_key, re_tls) in refmac_pdb_regex_dict.items():
            mx = re_tls.match(text)
            if mx != None:
                break

        ## no match
        if mx == None:
            return

        if re_key == "group":
            tls_info = TLSInfo()
            self.append(tls_info)

        elif re_key == "origin":
            (x, y, z) = mx.groups()
            tls_info.origin = (float(x), float(y), float(z))

        elif re_key == "range":
            (c1, f1, c2, f2) = mx.groups()
            sel = ""
            tls_info.range_list.append((c1, f1, c2, f2, sel))

        elif re_key == "t11_t22":
            (t11, t22) = mx.groups()
            if tls_info.T == None:
                T = (None, None, None, None, None, None)
            else:
                T = tls_info.T
            tls_info.T = (float(t11), float(t22), T[2], T[3], T[4], T[5])

        elif re_key == "t33_t12":
            (t33, t12) = mx.groups()
            if tls_info.T == None:
                T = (None, None, None, None, None, None)
            else:
                T = tls_info.T
            tls_info.T = (T[0], T[1], float(t33), float(t12), T[4], T[5])
            
        elif re_key == "t13_t23":
            (t13, t23) = mx.groups()
            if tls_info.T == None:
                T = (None, None, None, None, None, None)
            else:
                T = tls_info.T
            tls_info.T = (T[0], T[1], T[2], T[3], float(t13), float(t23))
            
        elif re_key == "l11_l22":
            (l11, l22) = mx.groups()
            if tls_info.L == None:
                L = (None, None, None, None, None, None)
            else:
                L = tls_info.L
            tls_info.L = (float(l11),
                          float(l22),
                          L[2], L[3], L[4], L[5])

        elif re_key == "l33_l12":
            (l33, l12) = mx.groups()
            if tls_info.L == None:
                L = (None, None, None, None, None, None)
            else:
                L = tls_info.L
            tls_info.L = (L[0], L[1],
                          float(l33), float(l12),
                          L[4], L[5])
            
        elif re_key == "l13_l23":
            (l13, l23) = mx.groups()
            if tls_info.L == None:
                L = (None, None, None, None, None, None)
            else:
                L = tls_info.L
            tls_info.L = (L[0], L[1], L[2], L[3],
                          float(l13), float(l23))

        elif re_key == "s11_s12_s13":
            (s11, s12, s13) = mx.groups()
            if tls_info.S == None:
                S = (None, None, None, None, None, None, None, None)
            else:
                S = tls_info.S
                
            #<S22 - S11> <S11 - S33> <S12> <S13> <S23> <S21> <S31> <S32>   

            tls_info.s11 = float(s11)

            if hasattr(tls_info, "s22"):
                s2211 = tls_info.s22 - tls_info.s11
            else:
                s2211 = None

            if hasattr(tls_info, "s33"):
                s1133 = tls_info.s11 - tls_info.s33
            else:
                s1133 = None

            tls_info.S = (
                s2211, s1133,
                float(s12), float(s13),
                S[4], S[5], S[6], S[7])

        elif re_key == "s21_s22_s23":
            (s21, s22, s23) = mx.groups()
            if tls_info.S == None:
                S = (None, None, None, None, None, None, None, None)
            else:
                S = tls_info.S
                
            tls_info.s22 = float(s22)

            if hasattr(tls_info, "s11"):
                s2211 = tls_info.s22 - tls_info.s11
            else:
                s2211 = None

            tls_info.S = (
                s2211, S[1], S[2], S[3],
                float(s23), float(s21),
                S[6], S[7])

        elif re_key == "s31_s32_s33":
            (s31, s32, s33) = mx.groups()
            if tls_info.S == None:
                S = (None, None, None, None, None, None, None, None)
            else:
                S = tls_info.S
                
            tls_info.s33 = float(s33)

            if hasattr(tls_info, "s11"):
                s1133 = tls_info.s11 - tls_info.s33
            else:
                s1133 = None

            tls_info.S = (
                S[0], s1133, S[2], S[3], S[4], S[5],
                float(s31), float(s32))


class TLSGroup(AtomList):
    """A subclass of AtomList implementing methods for performing TLS
    calculations on the contained Atom instances.
    """
    def __init__(self, *args):
        AtomList.__init__(self, *args)

        self.name           = "" 
        self.origin         = zeros(3, Float)
        self.T              = array([[0.0, 0.0, 0.0],
                                     [0.0, 0.0, 0.0],
                                     [0.0, 0.0, 0.0]], Float)
        self.L              = array([[0.0, 0.0, 0.0],
                                     [0.0, 0.0, 0.0],
                                     [0.0, 0.0, 0.0]], Float)
        self.S              = array([[0.0, 0.0, 0.0],
                                     [0.0, 0.0, 0.0],
                                     [0.0, 0.0, 0.0]], Float)

    def __str__(self):
        tstr  = "TLS %s\n" % (self.name)

        tstr += "ORIGIN  %f %f %f\n" % (
            self.origin[0], self.origin[1], self.origin[2])

        tstr += "T %f %f %f %f %f %f\n" % (
            self.T[0,0], self.T[1,1], self.T[2,2], self.T[0,1], self.T[0,2],
            self.T[1,2])

        tstr += "L %f %f %f %f %f %f\n" % (
            self.L[0,0], self.L[1,1], self.L[2,2], self.L[0,1], self.L[0,2],
            self.L[1,2])

        tstr += "S %f %f %f %f %f %f %f %f\n" % (
            self.S[1,1]-self.S[0,0], self.S[0,0]-self.S[2,2], self.S[0,1],
            self.S[0,2], self.S[1,2], self.S[1,0], self.S[2,0], self.S[2,1])

        return tstr

    def set_origin(self, x, y, z):
        """Sets the x, y, z components of the TLS origin vector.
        """
        self.origin = array([x, y, z])

    def set_T(self, t11, t22, t33, t12, t13, t23):
        """Sets the components of the symmetric T tensor.
        """
        self.T = array([[t11, t12, t13],
                        [t12, t22, t23],
                        [t13, t23, t33]])

    def set_L(self, l11, l22, l33, l12, l13, l23):
        """Sets the components of the symmetric L tensor from values
        assumed to be in degrees^2.
        """
        self.L = array([[l11, l12, l13],
                        [l12, l22, l23],
                        [l13, l23, l33]]) * deg2rad2

    def set_S(self, s2211, s1133, s12, s13, s23, s21, s31, s32):
        """Sets the componets of the asymmetric S tenssor.  The trace
        of the S tensor is set with the standard convention of
        the Trace(S) = 0.
        """
        s22 = 2.0*(s2211)/3.0 + s1133/3.0
        s11 = s22 - s2211
        s33 = s11 - s1133

        self.S = array([[s11, s12, s13],
                        [s21, s22, s23],
                        [s31, s32, s33]]) * deg2rad

    def calc_TLS_least_squares_fit(self):
        """Perform a least-squares fit of the atoms contained in self
        to the three TLS tensors: self.T, self.L, and self.S using the
        origin given by self.origin.
        """

        ## use label indexing to avoid confusion!
        T11, T22, T33, T12, T13, T23, L11, L22, L33, L12, L13, L23, \
        S11, S22, S33, S12, S13, S23, S21, S31, S32 = (
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)

        A = zeros((len(self)*6, 21), Float)
        b = zeros(len(self) * 6,  Float)

        i = -1
        for atm in self:
            i += 1

            ## set x, y, z as the vector components from the TLS origin
            x, y, z = atm.position - self.origin

            ## indecies of the components of U
            u11 =   i * 6
            u22 = u11 + 1
            u33 = u11 + 2
            u12 = u11 + 3
            u13 = u11 + 4
            u23 = u11 + 5

            ## set the B vector
            Uatm   = atm.get_U()            
            b[u11] = Uatm[0,0]
            b[u22] = Uatm[1,1]
            b[u33] = Uatm[2,2]
            b[u12] = Uatm[0,1]
            b[u13] = Uatm[0,2]
            b[u23] = Uatm[1,2]

            ## C Matrix
            xx = x*x
            yy = y*y
            zz = z*z

            xy = x*y
            xz = x*z
            yz = y*z

            A[u11, T11] = 1.0
            A[u11, L22] = zz
            A[u11, L33] = yy
            A[u11, L23] = -2.0*yz
            A[u11, S31] = -2.0*y
            A[u11, S21] = 2.0*z

            A[u22, T22] = 1.0
            A[u22, L11] = zz
            A[u22, L33] = xx
            A[u22, L13] = -2.0*xz
            A[u22, S12] = -2.0*z
            A[u22, S32] = 2.0*x

            A[u33, T33] = 1.0
            A[u33, L11] = yy
            A[u33, L22] = xx
            A[u33, L12] = -2.0*xy
            A[u33, S23] = -2.0*x
            A[u33, S13] = 2.0*y

            A[u12, T12] = 1.0
            A[u12, L33] = -xy
            A[u12, L23] = xz
            A[u12, L13] = yz
            A[u12, L12] = -zz
            A[u12, S11] = -z
            A[u12, S22] = z
            A[u12, S31] = x
            A[u12, S32] = -y

            A[u13, T13] = 1.0
            A[u13, L22] = -xz
            A[u13, L23] = xy
            A[u13, L13] = -yy
            A[u13, L12] = yz
            A[u13, S11] = y
            A[u13, S33] = -y
            A[u13, S23] = z
            A[u13, S21] = -x

            A[u23, T23] = 1.0
            A[u23, L11] = -yz
            A[u23, L23] = -xx
            A[u23, L13] = xy
            A[u23, L12] = xz
            A[u23, S22] = -x
            A[u23, S33] = x
            A[u23, S12] = y
            A[u23, S13] = -z

        (C, resids, rank, s) = linear_least_squares(A, b)

        ## verify the trase of S is zero
        assert allclose(C[S11]+C[S22]+C[S33], 0.0)

        self.T = array([ [ C[T11], C[T12], C[T13] ],
                         [ C[T12], C[T22], C[T23] ],
                         [ C[T13], C[T23], C[T33] ] ], Float)

        self.L = array([ [ C[L11], C[L12], C[L13] ],
                         [ C[L12], C[L22], C[L23] ],
                         [ C[L13], C[L23], C[L33] ] ], Float)

        self.S = array([ [ C[S11], C[S12], C[S13] ],
                         [ C[S21], C[S22], C[S23] ],
                         [ C[S31], C[S32], C[S33] ] ], Float)

    def calc_Utls(self, T, L, S, vec):
        """Returns the calculated value for the anisotropic U tensor for
        the given position.
        """
        x, y, z = vec

        xx = x*x
        yy = y*y
        zz = z*z

        xy = x*y
        yz = y*z
        xz = x*z

        u11 = T[0,0] \
              + L[1,1]*zz + L[2,2]*yy - 2.0*L[1,2]*yz \
              + 2.0*S[1,0]*z - 2.0*S[2,0]*y

        u22 = T[1,1] \
              + L[0,0]*zz + L[2,2]*xx - 2.0*L[2,0]*xz \
              - 2.0*S[0,1]*z + 2.0*S[2,1]*x

        u33 = T[2,2] \
              + L[0,0]*yy + L[1,1]*xx - 2.0*L[0,1]*xy \
              - 2.0*S[1,2]*x + 2.0*S[0,2]*y

        u12 = T[0,1] \
              - L[2,2]*xy + L[1,2]*xz + L[2,0]*yz - L[0,1]*zz \
              - S[0,0]*z + S[1,1]*z + S[2,0]*x - S[2,1]*y

        u13 = T[0,2] \
              - L[1,1]*xz + L[1,2]*xy - L[2,0]*yy + L[0,1]*yz \
              + S[0,0]*y - S[2,2]*y + S[1,2]*z - S[1,0]*x

        u23 = T[1,2] \
              - L[0,0]*yz - L[1,2]*xx + L[2,0]*xy + L[0,1]*xz \
              - S[1,1]*x + S[2,2]*x + S[0,1]*y - S[0,2]*z

        return array([[u11, u12, u13],
                      [u12, u22, u23],
                      [u13, u23, u33]])

    def iter_atm_Utls(self, T = None, L = None, S = None, o = None):
        """Iterates all the atoms in the TLS object, returning the 2-tuple
        (atm, U) where U is the calcuated U value from the current values
        of the TLS object's T,L,S, tensors and origin.
        """
        T      = T or self.T
        L      = L or self.L
        S      = S or self.S
        origin = o or self.origin
        
        for atm in self:
            Utls = self.calc_Utls(T, L, S, atm.position - origin)
            yield atm, Utls

    def calc_R(self):
        """Calculate the R factor of U vs. Utls.
        """
        Rn = 0.0
        Rd = 0.0

        for atm, Utls in self.iter_atm_Utls():
            U = atm.get_U()

            Rn += abs(U[0,0] - Utls[0,0])
            Rn += abs(U[1,1] - Utls[1,1])
            Rn += abs(U[2,2] - Utls[2,2])
            Rn += abs(U[0,1] - Utls[0,1])
            Rn += abs(U[0,2] - Utls[0,2])
            Rn += abs(U[1,2] - Utls[1,2])

            Rd += abs(U[0,0])
            Rd += abs(U[1,1])
            Rd += abs(U[2,2])
            Rd += abs(U[0,1])
            Rd += abs(U[0,2])
            Rd += abs(U[1,2])
            
        return Rn / Rd

    def calc_mean_S(self):
        """Calculates the mean Suij of U vs. Utls.
        """
        num      = 0
        mean_S   = 0.0
        MSD      = 0.0

        
        S_dict = {}
        for atm, Utls in self.iter_atm_Utls():
            Uatm = atm.get_U()
            
            try:
                atm_S = calc_Suij(Uatm, Utls)
            except ValueError:
                print "bad Suij"
            else:
                S_dict[atm] = atm_S
                mean_S += atm_S
                num += 1

        mean_S = mean_S / float(num)

        for S in S_dict.values():
            MSD += (mean_S - S)**2

        MSD = MSD / float(num)

        return mean_S, math.sqrt(MSD)
        
    def calc_mean_DP2(self):
        """Calculates the mean dP2 and standard deviation of U vs. Utls
        for every atom in the TLS group.
        Returns the 2-tuple. (mean DP2, rmsd)
        """
        num      = 0
        mean_DP2 = 0.0
        MSD      = 0.0
        
        DP2_dict = {}
        for atm, Utls in self.iter_atm_Utls():
            Uatm = atm.get_U()

            try:
                atm_DP2 = calc_DP2uij(Uatm, Utls)
            except ValueError:
                print "bad DP2uij"
            else:
                DP2_dict[atm] = atm_DP2 
                mean_DP2 += atm_DP2
                num += 1

        mean_DP2 = mean_DP2 / float(num)

        for DP2 in DP2_dict.values():
            MSD += (mean_DP2 - DP2)**2
        
        MSD = MSD / float(num)

        return mean_DP2, math.sqrt(MSD)
        
    def calc_mean_DP2N(self):
        """Calculates the mean DP2uij of the normalized U vs Utls.
        """
        B   = 10.0
        U2B = 8.0*math.pi*math.pi
        B2U = 1.0/U2B
        
        num       = 0
        mean_DP2N = 0.0
        MSD       = 0.0

        DP2N_dict = {}
        
        for atm, Utls in self.iter_atm_Utls():

            ## normalize the trace of Utls
            eval_U, evec_U = eigenvectors(Utls)
            evec_Ui        = inverse(evec_U)
            
            U2 = matrixmultiply(
                evec_U,
                matrixmultiply(Utls, transpose(evec_U)))
            
            U2    = U2 * U2B
            scale = B / trace(U2)
            U2    = scale * U2

            assert allclose(trace(U2), B)

            U2    = U2 * B2U
            UNtls = matrixmultiply(
                evec_Ui,
                matrixmultiply(U2, transpose(evec_Ui)))
            
            ## normalize the trace of atm.U
            U              = atm.get_U().copy()
            eval_U, evec_U = eigenvectors(U)
            evec_Ui        = inverse(evec_U)
            
            U2 = matrixmultiply(
                evec_U,
                matrixmultiply(U, transpose(evec_U)))

            U2    = U2 * U2B
            scale = B / trace(U2)
            U2    = scale * U2

            assert allclose(trace(U2), B)
            
            U2 = U2 * B2U
            UN = matrixmultiply(
                evec_Ui,
                matrixmultiply(U2, transpose(evec_Ui)))

            assert allclose(trace(UN), trace(UNtls))
            
            try:
                atm_DP2N = calc_DP2uij(UN, UNtls)
            except ValueError:
                pass
            else:
                num       += 1
                mean_DP2N += atm_DP2N
                DP2N_dict[atm] = atm_DP2N

        mean_DP2N = mean_DP2N / float(num)

        for DP2N in DP2N_dict.values():
            MSD += (mean_DP2N - DP2N)**2
        
        MSD = MSD / float(num)

        return mean_DP2N, math.sqrt(MSD)

    def shift_COR(self):
        """Shift the TLS group to the center of reaction.
        """
        calcs       = self.calc_COR()
        
        self.origin = calcs["COR"].copy()
        self.T      = calcs["T'"].copy()
        self.L      = calcs["L'"].copy()
        self.S      = calcs["S'"].copy()

        return calcs
        
    def calc_COR(self):
        """Calculate new tensors based on the center for reaction.
        This method returns a dictionary of the calculations:

        T^: T tensor in the coordinate system of L
        L^: L tensor in the coordinate system of L
        S^: S tensor in the coordinate system of L

        COR: Center of Reaction

        T',S',L': T,L,S tensors in origonal coordinate system
                  with the origin shifted to the center of reaction.
        """
        calcs = {}

        (eval_L, evec_L) = eigenvectors(self.L)
        evec_L = transpose(evec_L)
        
        ## carrot-L tensor (tensor WRT principal axes of L)
        cL      = zeros([3,3], Float)
        cL[0,0] = eval_L[0]
        cL[1,1] = eval_L[1]
        cL[2,2] = eval_L[2]

        calcs["L^"] = cL

        ## carrot-T tensor (T tensor WRT principal axes of L)
        cT = matrixmultiply(
            matrixmultiply(transpose(evec_L), self.T), evec_L)

        calcs["T^"] = cT

        ## carrot-S tensor (S tensor WRT principal axes of L)
        cS = matrixmultiply(
            matrixmultiply(transpose(evec_L), self.S), evec_L)

        ## correct for left-handed libration eigenvectors
        det = determinant(evec_L)
        if int(det) != 1:
            cS = -cS

        calcs["S^"] = cS
        
        ## ^rho: the origin-shift vector in the coordinate system of L
        small = 0.2 * deg2rad2

        cL1122 = cL[1,1] + cL[2,2]
        cL2200 = cL[2,2] + cL[0,0]
        cL0011 = cL[0,0] + cL[1,1]

        if cL1122>=small:
            crho0 = (cS[1,2]-cS[2,1]) / cL1122
        else:
            crho0 = 0.0

        if cL2200>=small:
            crho1 = (cS[2,0]-cS[0,2]) / cL2200
        else:
            crho1 = 0.0

        if cL0011>=small:
            crho2 = (cS[0,1]-cS[1,0]) / cL0011
        else:
            crho2 = 0.0

        crho = array([crho0, crho1, crho2], Float)

        calcs["RHO^"] = crho
        
        ## rho: the origin-shift vector in orthogonal coordinates
        rho = matrixmultiply(evec_L, crho)
        
        calcs["RHO"] = rho
        calcs["COR"] = array(self.origin, Float) + rho

        ## set up the origin shift matrix PRHO WRT orthogonal axes
        PRHO = array([ [    0.0,  rho[2], -rho[1]],
                       [-rho[2],     0.0,  rho[0]],
                       [ rho[1], -rho[0],     0.0] ], Float)

        ## set up the origin shift matrix cPRHO WRT libration axes
        cPRHO = array([ [    0.0,  crho[2], -crho[1]],
                        [-crho[2],     0.0,  crho[0]],
                        [ crho[1], -crho[0],     0.0] ], Float)

        ## calculate tranpose of cPRHO, ans cS
        cSt = transpose(cS)
        cPRHOt = transpose(cPRHO)

        ## calculate S'^ = S^ + L^*pRHOt
        cSp = cS + matrixmultiply(cL, cPRHOt)
        calcs["S'^"] = cSp

        ## L'^ = L^ = cL
        calcs["L'^"] = cL

        ## calculate T'^ = cT + cPRHO*S^ + cSt*cPRHOt + cPRHO*cL*cPRHOt *
        cTp = cT + \
              matrixmultiply(cPRHO, cS) + \
              matrixmultiply(cSt, cPRHOt) + \
              matrixmultiply(matrixmultiply(cPRHO, cL), cPRHOt)
        calcs["T'^"] = cTp

        ## transpose of PRHO and S
        PRHOt = transpose(PRHO)
        St = transpose(self.S)

        ## calculate S' = S + L*PRHOt
        Sp = self.S + matrixmultiply(self.L, PRHOt)
        calcs["S'"] = Sp

        ## calculate T' = T + PRHO*S + St*PRHOT + PRHO*L*PRHOt
        Tp = self.T + \
             matrixmultiply(PRHO, self.S) + \
             matrixmultiply(St, PRHOt) + \
             matrixmultiply(matrixmultiply(PRHO, self.L), PRHOt)
        calcs["T'"] = Tp

        ## L' is just L
        calcs["L'"] = self.L.copy()


        ### Verify that the the math we've just done is correct by
        ### comparing the original TLS calculated U tensors with the
        ### U tensors calculated from the COR
        T_cor = calcs["T'"].copy()
        L_cor = calcs["L'"].copy()
        S_cor = calcs["S'"].copy()

        for atm in self:
            x  = atm.position - self.origin
            xp = atm.position - calcs["COR"]

            Utls     = self.calc_Utls(self.T, self.L, self.S, x)
            Utls_cor = self.calc_Utls(T_cor, L_cor, S_cor, xp)

            assert allclose(Utls, Utls_cor)
        ###


        ## now calculate the TLS motion description using 3 non
        ## intersecting screw axes, with one

        ## libration axis 1 shift in the L coordinate system
        
        if cL[0,0]>=small:
            cL1rho = array([0.0, -cSp[0,2]/cL[0,0], cSp[0,1]/cL[0,0]], Float)
        else:
            cL1rho = zeros(3, Float)

        ## libration axis 2 shift in the L coordinate system
        if cL[1,1]>=small:
            cL2rho = array([cSp[1,2]/cL[1,1], 0.0, -cSp[1,0]/cL[1,1]], Float)
        else:
            cL2rho = zeros(3, Float)

        ## libration axis 2 shift in the L coordinate system
        if cL[2,2]>=small:
            cL3rho = array([-cSp[2,1]/cL[2,2], cSp[2,0]/cL[2,2], 0.0], Float)
            
        else:
            cL3rho = zeros(3, Float)

        ## libration axes shifts in the origional orthogonal
        ## coordinate system
        calcs["L1_rho"] = matrixmultiply(evec_L, cL1rho)
        calcs["L2_rho"] = matrixmultiply(evec_L, cL2rho)
        calcs["L3_rho"] = matrixmultiply(evec_L, cL3rho)

        ## calculate screw pitches (A*R / R*R) = (A/R)
        if cL[0,0]>=small:
            calcs["L1_pitch"] = cS[0,0]/cL[0,0]
        else:
            calcs["L1_pitch"] = 0.0
            
        if cL[1,1]>=small:
            calcs["L2_pitch"] = cS[1,1]/cL[1,1]
        else:
            calcs["L2_pitch"] = 0.0

        if cL[2,2]>=small:
            calcs["L3_pitch"] = cS[2,2]/cL[2,2]
        else:
            calcs["L3_pitch"] = 0.0

        ## now calculate the reduction in T for the screw rotation axes
        cTred = cT.copy()

        for i in range(3):
            for k in range(3):
                if i==k:
                    continue
                cTred[i,i] -= (cS[k,i]**2) / cL[k,k]

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    if j==i:
                        continue
                    cTred[i,j] -= (cS[k,i]*cS[k,j]) / cL[k,k]
                

        calcs["rT'"] = matrixmultiply(transpose(evec_L),
                                      matrixmultiply(cTred, evec_L))

        return calcs

    def check_positive_eigenvalues(self):
        """Returns True if the eigenvalues of the T and L tensors are all
        positive, otherwise returns False.
        """
        return min(eigenvalues(self.L))>=0 and min(eigenvalues(self.T))>=0.0

    def write(self, out = sys.stdout):
        """Write a nicely formatted tensor description.
        """
        ## the TLS tensors are in radians
        ## convert from radians to degrees
        T = self.T
        L = self.L * rad2deg2
        S = self.S * rad2deg

        (eval_t, evec_t) = eigenvectors(T)
        (eval_l, evec_l) = eigenvectors(L)
        
        listx = [
            "Tensor: T(A)",
            str(T),
            "Eigenvectors",
            str(evec_t),
            "Eigenvalues",
            str(eval_t),
            "",
            "Tensor: L(deg*deg)",
            str(L),
            "Eigenvectors",
            str(evec_l),
            "Eigenvalues",
            str(eval_l),
            "",
            "Tensor: S(deg*A)",
            str(S),
            ""
            ]

        out.write(string.join(listx, "\n"))


class TLSStructureAnalysis(object):
    """
    """
    def __init__(self, struct):
        self.struct = struct

    def iter_segments(self, chain, seg_len):
        """This iteratar yields a series of AtomList objects.  The first
        AtomList is all the atoms in the chain's first seg_len residues.
        The AtomLists yielded after the first one are formed by moving the
        starting residue up by one, and including all the atoms in the next
        seg_len residues.  The iteration terminates when there are are
        less than seg_len residues left in the chain from the iterators
        current position.
        """
        res_segment = []
        
        for res in chain.iter_amino_acids():
            res_segment.append(res)

            if len(res_segment)<seg_len:
                continue

            if len(res_segment)>seg_len:
                res_segment = res_segment[1:]

            atom_list = AtomList()
            for resx in res_segment:
                for atm in resx.iter_atoms():
                    atom_list.append(atm)

            yield res_segment[:], atom_list

    def fit_TLS_segments(self, **args):
        """Run the algorithm to fit TLS parameters to segments of the
        structure.  This method has many options, which are outlined in
        the source code for the method.  This returns a list of dictionaries
        containing statistics on each of the fit TLS groups, the residues
        involved, and the TLS object itself.
        """
        
        ## arguments
        origin                  = args.get("origin_of_calc")
        residue_width           = args.get("residue_width", 6)
        use_side_chains         = args.get("use_side_chains", True)
        filter_neg_eigen_values = args.get("filter_neg_eigen_values", True)
        include_hydrogens       = args.get("include_hydrogens", False)
        include_frac_occupancy  = args.get("include_frac_occupancy", False)
        include_single_bond     = args.get("include_single_bond", True)
        
        ## list of all TLS groups
        stats_list = []

        for chain in self.struct.iter_chains():
            for res_segment, seg_atom_list in self.iter_segments(
                chain, residue_width):

                stats             = {}
                stats["tls"]      = TLSGroup()
                stats["residues"] = res_segment
                stats["name"]     = "%s-%s" % (res_segment[0].fragment_id,
                                               res_segment[-1].fragment_id)

                tls      = stats["tls"]
                tls.name = stats["name"]

                ## add atoms into the TLSGroup
                ## filter the atoms going into the TLS group                
                for atm in seg_atom_list:
                    ## omit atoms which are not at full occupancy
                    if atm.occupancy<1.0:
                        continue
                    ## omit hydrogens
                    if include_hydrogens==False and atm.element=="H":
                        continue
                    ## omit side chain atoms
                    if use_side_chains==False:
                        if atm.name not in ("C", "N", "CA", "O"):
                            continue
                    ## omit atoms with a single bond 
                    if include_single_bond==False:
                        if len(atm.bond_list)<=1:
                            continue
                        
                    tls.append(atm)

                ## skip if there are not enough atoms for the TLS parameters
                ## the least squares fit of the TLS parameters requires at
                ## least 21 paramters (6 per atom)
                if len(tls)<4:
                    continue

                ## set the origin of the TLS group to the centroid, and also
                ## save it under calc_origin because origin will be overwritten
                ## using the COR after the least squares fit
                if origin!=None:
                    tls.origin = origin.copy()
                else:
                    tls.origin = tls.calc_centroid()

                stats["calc_origin"] = tls.origin

                ## calculate tensors and print
                tls.calc_TLS_least_squares_fit()

                ## shift the TLSGroup tensors and origin to the Center of
                ## Reaction
                calc = tls.shift_COR()

                ## negitive eigen values mean the TLS group cannot be
                ## physically intrepreted
                if filter_neg_eigen_values==True:
                    if min(eigenvalues(tls.L))<=0.0:
                        continue
                    elif min(eigenvalues(tls.T))<=0.0:
                        continue
                    elif  min(eigenvalues(calc["rT'"]))<=0.0:
                        continue

                ## this TLS group passes all our tests -- add it to the
                ## stats list
                stats_list.append(stats)

                stats["R"]                              = tls.calc_R()
                stats["mean_DP2"], stats["sigma_DP2"]   = tls.calc_mean_DP2()
                stats["mean_DP2N"], stats["sigma_DP2N"] = tls.calc_mean_DP2N()
                stats["mean_S"], stats["sigma_S"]       = tls.calc_mean_S()
                stats["num_atoms"]                      = len(tls)

        return stats_list


    def iter_chain_split(self, chain, num_splits, min_width):

        residues    = len(chain)
        split_point = min_width

        while split_point<(residues - min_width):
            seg1 = chain[:split_point]
            seg2 = chain[split_point+1:]

            if num_splits>1:

                for slist in self.iter_chain_split(
                    seg2, num_splits-1, min_width):

                    slist.insert(0, seg1)
                    yield slist

            else:
                yield [seg1, seg2]
                

    def split_TLS(self, **args):

        num_splits        = 5
        min_residue_width = 3


        
        

        

    def fit_common_TLS(self, **args):
        """Run the algorithm to fit TLS parameters to segments of the
        structure.  This method has many options, which are outlined in
        the source code for the method.  This returns a list of dictionaries
        containing statistics on each of the fit TLS groups, the residues
        involved, and the TLS object itself.
        """
        residue_width = args["residue_width"]
        
        ## calculate the centroid for all the TLS groiup fits
        al = AtomList()
        for atm in self.struct.iter_atoms():
            al.append(atm)
        origin = al.calc_centroid()

        ## list of all TLS groups
        stats_list = []
        for chain in self.struct.iter_chains():
            for res_segment, seg_atom_list in self.iter_segments(
                chain, residue_width):

                stats             = {}
                stats["tls"]      = TLSGroup()
                stats["residues"] = res_segment
                stats["name"]     = "%s-%s" % (res_segment[0].fragment_id,
                                               res_segment[-1].fragment_id)

                tls        = stats["tls"]
                tls.name   = stats["name"]
                tls.origin = origin.copy()
                
                for atm in seg_atom_list:
                    if atm.occupancy<1.0:
                        continue
                    tls.append(atm)

                ## calculate tensors and print
                tls.calc_TLS_least_squares_fit()

                ## negitive eigen values mean the TLS group cannot be
                ## physically intrepreted
                if min(eigenvalues(tls.L))<=0.0:
                    continue
                elif min(eigenvalues(tls.T))<=0.0:
                    continue

                ## this TLS group passes all our tests -- add it to the
                ## stats list
                stats_list.append(stats)

        ## find the common part of the TLS groups
        tls_common = TLSGroup()
        tls_common.origin = origin.copy()

        tls0         = stats_list[0]["tls"]
        tls_common.T = tls0.T.copy()
        tls_common.L = tls0.L.copy()
        tls_common.S = tls0.S.copy()

        for stats in stats_list[1:]:
            tls_s = stats["tls"]
            
            for C, S in ( (tls_common.T, tls_s.T),
                          (tls_common.L, tls_s.L),
                          (tls_common.S, tls_s.S) ):
                
                for i in (0,1,2):
                    for j in (0,1,2):

                        C[i,j] = min(C[i,j], S[i,j]) 

        print tls_common.T
        print tls_common.L
        print tls_common.S

        ## subtract TLS_common and shift all TLS groups to the COR
        for stats in stats_list:
            tls   = stats["tls"]

            tls.T = tls.T - tls_common.T
            tls.L = tls.L - tls_common.L
            tls.S = tls.S - tls_common.S

            calcs = tls.shift_COR()
            
            stats["R"]                              = tls.calc_R()
            stats["mean_DP2"], stats["sigma_DP2"]   = tls.calc_mean_DP2()
            stats["mean_DP2N"], stats["sigma_DP2N"] = tls.calc_mean_DP2N()
            stats["mean_S"], stats["sigma_S"]       = tls.calc_mean_S()
            stats["num_atoms"]                      = len(tls)
        
        return stats_list

    def iter_splits(self, length, min_seg, segments):
        """
        """
        pass

    def brute_fit(self, **args):
        """
        """
        num_tls_segments = args.get("num_tls_segments", 2)

    def join_tls_segments(self, stats_list):
        """
        Stage 1: fit with small segment length.
        Stage 2: examine segments, combine segment runs
        Stage 3: determine the longest segment run
        """
        pass


## <testing>
if __name__ == "__main__":
    print "==============================================="
    print "TEST CASE 1: TLS Class"
    print

    tls = TLS()
    tls.name = "All protein"
    tls.set_origin(18.885, 49.302, 13.315)
    tls.set_T(0.0263, 0.0561, 0.0048, -0.0128, 0.0065, -0.0157)
    tls.set_L(0.9730, 5.1496, 0.8488,  0.2151,-0.1296,  0.0815)
    tls.set_S(0.0007, 0.0281, 0.0336, -0.0446,-0.2288, -0.0551,
              0.0487, 0.0163)

    print tls

    print "eigenvalues(T)"
    print eigenvalues(tls.T)
    print "eigenvalues(L)"
    print eigenvalues(tls.L)

    print "==============================================="

    print 

    print "==============================================="
    print "TEST CASE 2: TLS/REFMAC file loading"

    import sys

    fil = open(sys.argv[1], "r")
    tls_list = LoadTLSGroups(fil)
    
    for tls in tls_list:
        print "----- TLS Group #%d ---" % (tls_list.index(tls))
        print tls

        print "eigenvalues(T)"
        print eigenvalues(tls.T)
        print "eigenvalues(L)"
        print eigenvalues(tls.L)

        print "-----------------------"

    print "==============================================="
## </testing>
