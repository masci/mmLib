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

from mmLib.mmTypes import *
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


class TLSGroupInfo:
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
        """Return a string describing the TLS group in REFMAC/CCP4 format.
        """
        def ft8(tx):
            strx = ""
            for x in tx:
                strx += fpformat.fix(x, 4).rjust(8)
            return strx
        
        listx = []

        if len(self.name) > 0:
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

        listx.append("ORIGIN %s" % ft8(self.origin))
        listx.append("T   %s" % ft8(self.T))
        listx.append("L   %s" % ft8(self.L))
        listx.append("S   %s" % ft8(self.S))

        return string.join(listx, "\n")

    def make_tls_group(self, struct):
        """Returns a TLSGroup containing the appropriate atoms from the
        argument Structure object, and with the origin, T, L, S tensors
        set.
        """
        tls = TLSGroup()
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


class TLSGroupFile:
    """This class reads and writes TLS information stored in the same
    format as REFMAC from CCP4 >= 4.1.0.  The TLS groups are stored as a list
    of TLSGroupInfo classes.
    """
    def __init__(self):
        self.tls_info_list = []
        self.state = {}

    def __str__(self):
        """Return a string describing the TLS groups in REFMAC/CCP4 format.
        """
        listx = []
        
        for tls_info in self.tls_info_list:
            listx.append(str(tls_info))

        return string.join(listx, "\n\n")
        
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
                tls_info = TLSGroupInfo()
                self.tls_info_list.append(tls_info)
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
                (t11, t22, t33, t12, t13, t23) = mx.groups()
                tls_info.T = (float(t11), float(t22), float(t33),
                              float(t12), float(t13), float(t23))

            elif re_key == "L":
                (l11, l22, l33, l12, l13, l23) = mx.groups()
                tls_info.L = (float(l11),
                              float(l22),
                              float(l33),
                              float(l12),
                              float(l13),
                              float(l23))

            elif re_key == "S":
                (s2211, s1133, s12, s13, s23, s21, s31, s32) = mx.groups()
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
            tls_info = self.tls_info_list[-1]
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
            tls_info = TLSGroupInfo()
            self.tls_info_list.append(tls_info)

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
    def __init__(self):
        AtomList.__init__(self)

        self.name           = "" 
        self.origin         = Vector(0.0, 0.0, 0.0)
        self.T              = array([[0.0, 0.0, 0.0],
                                     [0.0, 0.0, 0.0],
                                     [0.0, 0.0, 0.0]])
        self.L              = array([[0.0, 0.0, 0.0],
                                     [0.0, 0.0, 0.0],
                                     [0.0, 0.0, 0.0]])
        self.S              = array([[0.0, 0.0, 0.0],
                                     [0.0, 0.0, 0.0],
                                     [0.0, 0.0, 0.0]])

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
        self.origin = Vector(x, y, z)

    def set_T(self, t11, t22, t33, t12, t13, t23):
        """Sets the components of the symmetric T tensor.
        """
        self.T[0,0] = t11
        self.T[1,1] = t22
        self.T[2,2] = t33
        self.T[0,1] = t12
        self.T[1,0] = t12
        self.T[0,2] = t13
        self.T[2,0] = t13
        self.T[1,2] = t23
        self.T[2,1] = t23

    def set_L(self, l11, l22, l33, l12, l13, l23):
        """Sets the components of the symmetric L tensor from values
        assumed to be in degrees^2.
        """
        self.L[0,0] = l11 
        self.L[1,1] = l22
        self.L[2,2] = l33
        self.L[0,1] = l12
        self.L[1,0] = l12
        self.L[0,2] = l13
        self.L[2,0] = l13
        self.L[1,2] = l23
        self.L[2,1] = l23

        self.L = self.L * deg2rad2

    def set_S(self, s2211, s1133, s12, s13, s23, s21, s31, s32):
        """Sets the componets of the asymmetric S tenssor.  The trace
        of the S tensor is set with the standard convention of
        the Trace(S) = 0.
        """
        s22 = 2.0*(s2211)/3.0 + s1133/3.0
        s11 = s22 - s2211
        s33 = s11 - s1133

        self.S[0,0] = s11
        self.S[1,1] = s22
        self.S[2,2] = s33
        self.S[0,1] = s12
        self.S[0,2] = s13
        self.S[1,0] = s21
        self.S[1,2] = s23
        self.S[2,0] = s31
        self.S[2,1] = s32

        self.S = self.S * deg2rad

    def calc_tls_tensors(self):
        """Perform a least-squares fit of the atoms contained in self
        to the three TLS tensors: self.T, self.L, and self.S using the
        origin given by self.origin.
        """
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
            if atm.U == None:
                b[u11] = atm.temp_factor / (24.0 * math.pi * math.pi)
                b[u22] = b[u11]
                b[u33] = b[u11]
            else:
                b[u11] = atm.U[0,0]
                b[u22] = atm.U[1,1]
                b[u33] = atm.U[2,2]
                b[u12] = atm.U[0,1]
                b[u13] = atm.U[0,2]
                b[u23] = atm.U[1,2]

            ## C Matrix
            xx = x*x
            yy = y*y
            zz = z*z

            xy = x*y
            xz = x*z
            yz = y*z

            A[u11,  0] = 1.0
            A[u11, 10] = 2.0*z
            A[u11, 14] = zz
            A[u11, 15] = -2.0*y
            A[u11, 19] = -2.0*yz
            A[u11, 20] = yy

            A[u12,  1] = 1.0
            A[u12,  6] = -z
            A[u12,  7] = -2.0*y
            A[u12, 11] = z
            A[u12, 13] = -zz
            A[u12, 15] = x
            A[u12, 16] = -y
            A[u12, 18] = yz
            A[u12, 19] = xz
            A[u12, 20] = -xy

            A[u22,  2] = 1.0
            A[u22,  9] = zz
            A[u22, 16] = 2.0*x
            A[u22, 18] = -2.0*xz
            A[u22, 20] = xx

            A[u13,  3] = 1.0
            A[u13,  6] = y
            A[u13, 10] = -x
            A[u13, 12] = z
            A[u13, 13] = yz
            A[u13, 14] = -xz
            A[u13, 17] = -y
            A[u13, 18] = -yy
            A[u13, 19] = xy

            A[u23,  4] = 1.0
            A[u23,  7] = y
            A[u23,  8] = -z
            A[u23,  9] = -yz
            A[u23, 11] = -x
            A[u23, 13] = xz
            A[u23, 17] = x
            A[u23, 18] = xy
            A[u23, 19] = -xx

            A[u33,  5] = 1.0
            A[u33,  8] = 2.0*y
            A[u33,  9] = yy
            A[u33, 12] = -2.0*x
            A[u33, 13] = -2.0*xy
            A[u33, 14] = xx

        (C, resids, rank, s) = linear_least_squares(A, b)

        self.T = array([ [ C[ 0], C[ 1], C[ 3] ],
                         [ C[ 1], C[ 2], C[ 4] ],
                         [ C[ 3], C[ 4], C[ 5] ] ])

        self.L = array([ [ C[ 9], C[13], C[18] ],
                         [ C[13], C[14], C[19] ],
                         [ C[18], C[19], C[20] ] ])

        self.S = array([ [ C[ 6], C[ 7], C[ 8] ],
                         [ C[10], C[11], C[12] ],
                         [ C[15], C[16], C[17] ] ])

    def iter_atm_Ucalc(self):
        """Iterates all the atoms in the TLS object, returning the 2-tuple
        (atm, U) where U is the calcuated U value from the current values
        of the TLS object's T,L,S, tensors and origin.
        """
        T = self.T
        L = self.L
        S = self.S
        O = self.origin
        
        for atm in self:
            x, y, z = atm.position - O

            xx = x*x
            yy = y*y
            zz = z*z

            xy = x*y
            yz = y*z
            xz = x*z

            u00 = (         T[0,0]
                    + 2.0 * S[1,0] * z
                    +       L[1,1] * zz 
                    - 2.0 * S[2,0] * y 
                    - 2.0 * L[2,1] * yz
                    +       L[2,2] * yy )
            u10 = (         T[1,0]
                    -       S[0,0] * z
                    - 2.0 * S[0,1] * y
                    +       S[1,1] * z
                    -       L[1,0] * zz
                    +       S[2,0] * x
                    -       S[2,1] * y
                    +       L[2,0] * yz
                    +       L[2,1] * xz
                    -       L[2,2] * xy )
            u11 = (         T[1,1]
                    +       L[0,0] * zz
                    + 2.0 * S[2,1] * x
                    - 2.0 * L[2,0] * xz
                    +       L[2,2] * xx )
            u20 = (         T[2,0]
                    +       S[0,0] * y
                    -       S[1,0] * x
                    +       S[1,2] * z
                    +       L[1,0] * yz
                    -       L[1,1] * xz
                    -       S[2,2] * y
                    -       L[2,0] * yy
                    +       L[2,1] * xy )
            u21 = (         T[2,1]
                    +       S[0,1] * y
                    -       S[0,2] * z
                    -       L[0,0] * yz
                    -       S[1,1] * x
                    +       L[1,0] * xz
                    +       S[2,2] * x 
                    +       L[2,0] * xy
                    -       L[2,1] * xx )
            u22 = (         T[2,2]
                    + 2.0 * S[0,2] * y
                    +       L[0,0] * yy
                    - 2.0 * S[1,2] * x
                    - 2.0 * L[1,0] * xy
                    +       L[1,1] * xx )

            U = array ([[u00, u10, u20],
                        [u10, u11, u21],
                        [u20, u21, u22]])

            yield (atm, U)

    def calc_R(self):
        """Calculate the R factor of U vs. Ucalc.
        """
        Rn = 0.0
        Rd = 0.0

        for (atm, Ucalc) in self.iter_atm_Ucalc():
            if atm.U == None:
                Uobs = identity(3) * (atm.temp_factor/(24.0*math.pi*math.pi))
            else:
                Uobs = atm.U

            for i in range(3):
                for j in range(3):
                    Rn += abs(Uobs[i,j] - Ucalc[i,j])
                    Rd += abs(Uobs[i,j])

        return Rn / Rd

    def calc_COR(self):
        """Calculate new tensors based on the center for reaction.
        This method returns a dictionary of the calculations:

        T^: T tensor in the coordinate system of L
        L^: L tensor in the coordinate system of L
        S^: S tensor in the coordinate system of L
        
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
        crho = array([ (cS[1,2] - cS[2,1]) / (cL[1,1] + cL[2,2]),
                       (cS[2,0] - cS[0,2]) / (cL[2,2] + cL[0,0]),
                       (cS[0,1] - cS[1,0]) / (cL[0,0] + cL[1,1]) ])

        calcs["RHO^"] = crho
        
        ## rho: the origin-shift vector in orthogonal coordinates
        rho = matrixmultiply(evec_L, crho)
        
        calcs["RHO"] = rho
        calcs["COR"] = array(self.origin) + rho

        ## set up the origin shift matrix PRHO WRT orthogonal axes
        PRHO = array([ [    0.0,  rho[2], -rho[1]],
                       [-rho[2],     0.0,  rho[0]],
                       [ rho[1], -rho[0],     0.0] ])

        ## set up the origin shift matrix cPRHO WRT libration axes
        cPRHO = array([ [    0.0,  crho[2], -crho[1]],
                        [-crho[2],     0.0,  crho[0]],
                        [ crho[1], -crho[0],     0.0] ])

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

        # calculate T' = T + PRHO*S + St*PRHOT + PRHO*L*PRHOt
        Tp = self.T + \
             matrixmultiply(PRHO, self.S) + \
             matrixmultiply(St, PRHOt) + \
             matrixmultiply(matrixmultiply(PRHO, self.L), PRHOt)
        calcs["T'"] = Tp

        ## L' is just L
        calcs["L'"] = self.L

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
