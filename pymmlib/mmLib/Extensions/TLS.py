
import re
from Numeric import *
from LinearAlgebra import *
from Scientific.Geometry import Vector


def LoadTLSGroups(fil):
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

    tls_list = []

    while 1:
        ln = fil.readline()
        if not ln: break
        if not ln.strip(): continue

        if ln[:3] == "TLS":
            tls      = TLS()
            tls.name = ln[4:].strip()

            while 1:
                ln = fil.readline()
                if not ln: break
                
                m = re_RANGE.match(ln)
                if m:
                    (chain1, frag_id1, chain2, frag_id2, cmd) = m.groups()
                    tls.range_list.append(
                        ((chain1, frag_id1),(chain2, frag_id2)) )
                    continue
                    
                m = re_ORIGIN.match(ln)
                if m:
                    (x, y, z) = m.groups()
                    tls.setOrigin(float(x), float(y), float(z))
                    continue

                m = re_T.match(ln)
                if m:
                    (t11, t22, t33, t12, t13, t23) = m.groups()
                    tls.setT(float(t11), float(t22), float(t33),
                             float(t12), float(t13), float(t23))
                    continue

                m = re_L.match(ln)
                if m:
                    (l11, l22, l33, l12, l13, l23) = m.groups()
                    tls.setL(float(l11), float(l22), float(l33),
                             float(l12), float(l13), float(l23))
                    continue

                m = re_S.match(ln)
                if m:
                    (s2211, s1133, s12, s13, s23, s21, s31, s32) = m.groups()
                    tls.setS(float(s2211), float(s1133),
                             float(s12), float(s13), float(s23),
                             float(s21), float(s31), float(s32))
                    break

            tls_list.append(tls)

    return tls_list


class TLS:
    """Read, write, and perform calculations on TLS parameters.  The TLS
    parameter file supports REFMAC TLS format."""

    def __init__(self):
        self.name           = "" 
        self.range_list     = []
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

    def setOrigin(self, x, y, z):
        self.origin = Vector(x, y, z)
        
    def setT(self, t11, t22, t33, t12, t13, t23):
        self.T = array([[t11, t12, t13],
                        [t12, t22, t23],
                        [t13, t23, t33]])

    def setL(self, l11, l22, l33, l12, l13, l23):
        self.L = array([[l11, l12, l13],
                        [l12, l22, l23],
                        [l13, l23, l33]])

    def setS(self, s2211, s1133, s12, s13, s23, s21, s31, s32):
        ## assume: s11 + s22 + s33 = 0
        s22 = 2.0*(s2211)/3.0 + s1133/3.0
        s11 = s22 - s2211
        s33 = s11 - s1133
        self.S = array([[s11, s12, s13],
                        [s21, s22, s23],
                        [s31, s32, s33]])


if __name__ == "__main__":

    import sys

    fil = open(sys.argv[1], "r")
    tls_list = LoadTLSGroups(fil)

    
    for tls in tls_list:
        print "----- TLS Group #%d ---" % (tls_list.index(tls))
        print tls

    sys.exit(0)

    t = array([[ 0.0263, -0.0128,  0.0065],
               [-0.0128,  0.0561, -0.0157],
               [ 0.0065, -0.0157,  0.0048]])

    l = array([[ 0.9730,  0.2151, -0.1296],
               [ 0.2151,  5.1496,  0.0815],
               [-0.1296,  0.0815,  0.8488]])

    print "T", eigenvalues(t)
    print "L", eigenvalues(l)
