## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Utility classes for loading, manipulating, and analyzing TLS parameters.
"""
import re
import sys
import string

from mmLib.mmTypes import *
from mmLib.Structure import *


class RefmacTLSFile:
    """This class reads and writes TLS information stored in the same
    format as REFMAC from CCP4 >= 4.1.0.
    """

    re_RANGE  = re.compile(
        "^RANGE\s+\'([A-Z])\s*([^\']+)\'\s+\'([A-Z])\s*([^\']+)\'\s*(\w+).*$")
    re_ORIGIN = re.compile(
        "^ORIGIN\s+(\S+)\s+(\S+)\s+(\S+).*$")
    re_T      = re.compile(
        "^T\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+).*$")
    re_L      = re.compile(
        "^L\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+).*$")
    re_S      = re.compile(
        "^S\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)"\
        "\s+(\S+).*$")

    def __init__(self):
        self.tls_list = []

    def __str__(self):
        tstr  = "TLS %s\n" % (self.name)
        for ((chain_id1, frag_id1),(chain_id2, frag_id2)) in self.range_list:
            tstr += "RANGE '%s%s' '%s%s'ALL\n" % (
                chain_id1, frag_id1.rjust(5), chain_id2, frag_id2.rjust(5))
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

    def load(self, fil):
        """Read the TLS information from file object fil, and store that
        inforamtion in the class instance variables.
        """
        self.tls_list = []

        while 1:
            ln = fil.readline()
            if not ln: break
            if not ln.strip(): continue

            if ln[:3] == "TLS":
                tls = TLS()
                residue_selection_list = []

                tls.name = ln[4:].strip()

                while 1:
                    ln = fil.readline()
                    if not ln: break

                    m = self.re_RANGE.match(ln)
                    if m:
                        (chain1, frag_id1, chain2, frag_id2, cmd)=m.groups()
                        residue_selection_list.append(
                            (chain1, frag_id1, chain2, frag_id2, cmd))
                        continue

                    m = self.re_ORIGIN.match(ln)
                    if m:
                        (x, y, z)=m.groups()
                        tls.setOrigin(float(x), float(y), float(z))
                        continue

                    m = self.re_T.match(ln)
                    if m:
                        (t11, t22, t33, t12, t13, t23)=m.groups()
                        tls.setT(float(t11), float(t22), float(t33),
                                 float(t12), float(t13), float(t23))
                        continue

                    m = self.re_L.match(ln)
                    if m:
                        (l11, l22, l33, l12, l13, l23)=m.groups()
                        tls.setL(float(l11), float(l22), float(l33),
                                 float(l12), float(l13), float(l23))
                        continue

                    m = self.re_S.match(ln)
                    if m:
                        (s2211, s1133, s12, s13, s23, s21, s31, s32)=m.groups()
                        tls.setS(float(s2211), float(s1133),
                                 float(s12), float(s13), float(s23),
                                 float(s21), float(s31), float(s32))
                        break

                tls_list.append((tls, residue_selection_list))

    def set_structure_atoms(self, struct):
        """For each TLS group in the Refmac input file, there are one
        or more RANGE records describing the residues included in the
        TLS group.  This method uses these ranges to set each TLS
        instance tls.atom_list method.
        """
        pass


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
        self.origin[0] = x
        self.origin[1] = y
        self.origin[2] = z

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
        """Sets the components of the symmetric L tensor.
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

    def set_S(self, s2211, s1133, s12, s13, s23, s21, s31, s32):
        """Sets teh componets of the asymmetric S tenssor.  The trace
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

    def calc_tls_tensors(self):
        """Perform a least-squares fit of the atoms contained in self
        to the three TLS tensors: self.T, self.L, and self.S using the
        origin given by self.origin.
        """
        A = zeros((len(self)*6, 21), Float)
        b = zeros((len(self)*6, 1),  Float)

        for i in range(len(self)):
            atm = self[i]

            ## set x, y, z as the vector components from the TLS origin
            r = atm.position - self.origin
            x = r[0]
            y = r[1]
            z = r[2]

            ## indecies of the components of U
            u11 =   i * 6
            u22 = u11 + 1
            u33 = u11 + 2
            u12 = u11 + 3
            u13 = u11 + 4
            u23 = u11 + 5

            ## set the U vector
            b[u11, 0] = atm.U[0,0]
            b[u22, 0] = atm.U[1,1]
            b[u33, 0] = atm.U[2,2]
            b[u12, 0] = atm.U[0,1]
            b[u13, 0] = atm.U[0,2]
            b[u23, 0] = atm.U[1,2]

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

        self.T = array([ [ C[ 0,0], C[ 1,0], C[ 3,0] ],
                         [ C[ 1,0], C[ 2,0], C[ 4,0] ],
                         [ C[ 3,0], C[ 4,0], C[ 5,0] ] ])

        self.L = array([ [ C[ 9,0], C[13,0], C[18,0] ],
                         [ C[13,0], C[14,0], C[19,0] ],
                         [ C[18,0], C[19,0], C[20,0] ] ])

        self.S = array([ [ C[ 6,0], C[ 7,0], C[ 8,0] ],
                         [ C[10,0], C[11,0], C[12,0] ],
                         [ C[15,0], C[16,0], C[17,0] ] ])

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
