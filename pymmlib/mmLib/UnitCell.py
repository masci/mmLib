## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Classes for handling unit cell transformation.
"""
import math
from mmTypes import *
from AtomMath import *
from SpaceGroups import GetSpaceGroup, SymOp


class UnitCell(object):
    """Class for storing and performing calculations on unit cell
    parameters.  The constructor expects alpha, beta, and gamma to be in
    degrees, but converts them to radians.  Set angle_units = "rad" if
    the alpha, beta, and gamma are already in radians.
    """
    def __init__(self,
                 a = 1.0,
                 b = 1.0,
                 c = 1.0,
                 alpha = 90.0,
                 beta = 90.0,
                 gamma = 90.0,
                 space_group = "P1",
                 angle_units = "deg"):

        assert angle_units == "deg" or angle_units == "rad"

        self.a = a
        self.b = b
        self.c = c

        if angle_units == "deg":
            self.alpha = math.radians(alpha) 
            self.beta  = math.radians(beta)
            self.gamma = math.radians(gamma)

        elif angle_units == "rad":
            self.alpha = alpha
            self.beta  = beta
            self.gamma = gamma

        self.space_group = GetSpaceGroup(space_group)

        self.orth_to_frac = self.calc_fractionalization_matrix()
        self.frac_to_orth = self.calc_orthogonalization_matrix()

        assert allclose(self.orth_to_frac, inverse(self.frac_to_orth))
        print self

        print self.frac_to_orth


    def __str__(self):
        alpha = self.alpha * rad2deg
        beta = self.beta * rad2deg
        gamma = self.gamma * rad2deg

        return "UnitCell(a=%f, b=%f, c=%f, alpha=%f, beta=%f, gamma=%f)" % (
            self.a, self.b, self.c, alpha, beta, gamma)

    def calc_alpha_deg(self):
        """Returns the alpha angle in degrees.
        """
        return self.alpha * rad2deg
    
    def calc_beta_deg(self):
        """Returns the beta angle in degrees.
        """
        return self.beta * rad2deg

    def calc_gamma_deg(self):
        """Returns the gamma angle in degrees.
        """
        return self.gamma * rad2deg

    def calc_v(self):
        """Calculates the volume of the rhombohedrial created by the
        unit vectors a1/|a1|, a2/|a2|, a3/|a3|.
        """
        cos_alpha = math.cos(self.alpha)
        cos_beta  = math.cos(self.beta)
        cos_gamma = math.cos(self.gamma)
        return math.sqrt(1                       -
                         (cos_alpha * cos_alpha) -
                         (cos_beta  * cos_beta)  -
                         (cos_gamma * cos_gamma) +
                         (2 * cos_alpha * cos_beta * cos_gamma) )

    def calc_volume(self):
        """Calculates the volume of the unit cell.
        """
        return self.a * self.b * self.c * self.calc_v()

    def calc_reciprocal_unit_cell(self):
        """Corresponding reciprocal unit cell.
        """
        V = self.calc_volume()

        sin_alpha = math.sin(self.alpha)
        sin_beta  = math.sin(self.beta)
        sin_gamma = math.sin(self.gamma)

        cos_alpha = math.cos(self.alpha)
        cos_beta  = math.cos(self.beta)
        cos_gamma = math.cos(self.gamma)

        ra  = (self.b * self.c * sin_alpha) / V
        rb  = (self.a * self.c * sin_beta)  / V
        rc  = (self.a * self.b * sin_gamma) / V

        ralpha = math.acos(
            (cos_beta  * cos_gamma - cos_alpha) / (sin_beta  * sin_gamma))
        rbeta  = math.acos(
            (cos_alpha * cos_gamma - cos_beta)  / (sin_alpha * sin_gamma))
        rgamma = math.acos(
            (cos_alpha * cos_beta  - cos_gamma) / (sin_alpha * sin_beta))

        return UnitCell(ra, rb, rc, ralpha, rbeta, rgamma)

    def calc_orthogonalization_matrix(self):
        """Cartesian to fractional coordinates.
        """
        sin_alpha = math.sin(self.alpha)
        sin_beta  = math.sin(self.beta)
        sin_gamma = math.sin(self.gamma)

        cos_alpha = math.cos(self.alpha)
        cos_beta  = math.cos(self.beta)
        cos_gamma = math.cos(self.gamma)

        v = self.calc_v()

        f11 = self.a
        f12 = self.b * cos_gamma
        f13 = self.c * cos_beta
        f22 = self.b * sin_gamma
        f23 = (self.c * (cos_alpha - cos_beta * cos_gamma)) / (sin_gamma)
        f33 = (self.c * v) / sin_gamma

        orth_to_frac = array([ [f11, f12, f13],
                               [0.0, f22, f23],
                               [0.0, 0.0, f33] ])
        
        return orth_to_frac
        
    def calc_fractionalization_matrix(self):
        """Fractional to cartesian coordinates.
        """
        sin_alpha = math.sin(self.alpha)
        sin_beta  = math.sin(self.beta)
        sin_gamma = math.sin(self.gamma)

        cos_alpha = math.cos(self.alpha)
        cos_beta  = math.cos(self.beta)
        cos_gamma = math.cos(self.gamma)

        v = self.calc_v()

        o11 = 1.0 / self.a
        o12 = - cos_gamma / (self.a * sin_gamma)
        o13 = (cos_gamma * cos_alpha - cos_beta) / (self.a * v * sin_gamma)
        o22 = 1.0 / (self.b * sin_gamma)
        o23 = (cos_gamma * cos_beta - cos_alpha) / (self.b * v * sin_gamma)
        o33 = sin_gamma / (self.c * v)

        frac_to_orth = array([ [o11, o12, o13],
                               [0.0, o22, o23],
                               [0.0, 0.0, o33] ])

        return frac_to_orth

    def calc_orth_to_frac(self, v):
        """Calculates and returns the fractional coordinate vector of
        orthogonal vector v.
        """
        return matrixmultiply(self.orth_to_frac, v)

    def calc_frac_to_orth(self, v):
        """Calculates and returns the orthogonal coordinate vector of
        fractional vector v.
        """
        return matrixmultiply(self.frac_to_orth, v)

    def calc_symop_RT_old(self, symop):
        ## rotation matrix
        x1 = self.calc_orth_to_frac(array([1.0, 0.0, 0.0]))
        y1 = self.calc_orth_to_frac(array([0.0, 1.0, 0.0]))
        z1 = self.calc_orth_to_frac(array([0.0, 0.0, 1.0]))

        x2 = matrixmultiply(symop[0], x1)
        y2 = matrixmultiply(symop[0], y1)
        z2 = matrixmultiply(symop[0], z1)

        x3 = self.calc_frac_to_orth(x2)
        y3 = self.calc_frac_to_orth(y2)
        z3 = self.calc_frac_to_orth(z2)

        R = array([ [x3[0], y3[0], z3[0]],
                    [x3[1], y3[1], z3[1]],
                    [x3[2], y3[2], z3[2]] ])

        ## translation matrix
        T = self.calc_frac_to_orth(symop[1])

        x  = "[%6.3f %6.3f %6.3f %6.3f]\n" % (
            R[0,0], R[0,1], R[0,2], T[0])
        x += "[%6.3f %6.3f %6.3f %6.3f]\n" % (
            R[1,0], R[1,1], R[1,2], T[1])
        x += "[%6.3f %6.3f %6.3f %6.3f]\n" % (
            R[2,0], R[2,1], R[2,2], T[2])

        return R, T

    def calc_symop_RT(self, symop):
        R = matrixmultiply(
            self.frac_to_orth,
            matrixmultiply(symop[0], self.orth_to_frac))

        T = matrixmultiply(self.frac_to_orth, symop[1])

        return R, T
    
    def calc_xyz_coords(self, atom_list, symop):
        """
        """
        xyz_dict = {}
        for atm in atom_list:
            xyz_dict[atm] = symop(self.calc_orth_to_frac(atm.position))
        return xyz_dict

    def calc_cell(self, xyz):
        """Returns the cell integer 3-Tuple where the xyz fractional
        coordinates are located.
        """
        if xyz[0]<0.0:
            cx = int(xyz[0] - 1.0)
        else:
            cx = int(xyz[0] + 1.0)
            
        if xyz[1]<0.0:
            cy = int(xyz[1] - 1.0)
        else:
            cy = int(xyz[1] + 1.0)

        if xyz[2]<0.0:
            cz = int(xyz[2] - 1.0)
        else:
            cz = int(xyz[2] + 1.0)
            
        return (cx, cy, cz)

    def cell_search_iter(self):
         for i in (-1.0, 0.0, 1.0):
                for j in (-1.0, 0.0, 1.0):
                    for k in (-1.0, 0.0, 1.0):
                        yield i, j, k

    def iter_RT(self, struct):
        """
        """
        print "===================================="

        
        atom_list = list(struct.iter_atoms())

        fill_cell = (1,1,1)

        for symop in self.space_group.iter_symops():
            
            print "NEW SYMOP"
            print symop

            for i, j, k in self.cell_search_iter():
                cell_tra = array([i, j, k])

                print "## shifting translation = %s" % (str(cell_tra))
                s2 = SymOp((symop[0], symop[1] + cell_tra))

                in_target = False

                for atm in atom_list:
                    xyz       = s2(self.calc_orth_to_frac(atm.position))
                    calc_cell = self.calc_cell(xyz)

                    if calc_cell==fill_cell:
                        in_target = True
                        
                        print "## fill_cell=%s calc_cell=%s" % (str(fill_cell), str(calc_cell))
                        
                        R, T  = self.calc_symop_RT(s2)

                        print "## R,T"
                        print strRT(R, T)

                        atm_pos = matrixmultiply(R, atm.position) + T
                        print "## pos=%s xyz=%s:%s sym=%s" % (
                            str(atm.position),
                            str(xyz),
                            str(self.calc_frac_to_orth(xyz)),
                            str(atm_pos))

                        break

                    
                if in_target==True:
                    yield self.calc_symop_RT(s2)

def strRT(R, T):
    x  = "[%6.3f %6.3f %6.3f %6.3f]\n" % (
        R[0,0], R[0,1], R[0,2], T[0])
    x += "[%6.3f %6.3f %6.3f %6.3f]\n" % (
        R[1,0], R[1,1], R[1,2], T[1])
    x += "[%6.3f %6.3f %6.3f %6.3f]\n" % (
        R[2,0], R[2,1], R[2,2], T[2])
    
    return x
        
## <testing>
def main():
    print "================================================="
    print "TEST CASE #1: Triclinic unit cell"
    print

    uc = UnitCell(7.877, 7.210, 7.891, 105.563, 116.245, 79.836)

    e = array([[1.0, 0.0, 0.0],
               [0.0, 1.0, 0.0],
               [0.0, 0.0, 1.0]])


    print uc
    print "volume                   = ",uc.calc_v()
    print "cell volume              = ",uc.calc_volume()
    print "fractionalization matrix =\n",uc.calc_fractionalization_matrix()
    print "orthogonalization matrix =\n",uc.calc_orthogonalization_matrix()

    print "orth * e =\n",matrixmultiply(
        uc.calc_orthogonalization_matrix(), e)


    print "calc_frac_to_orth"
    vlist = [
        array([0.0, 0.0, 0.0]),
        array([0.5, 0.5, 0.5]),
        array([1.0, 1.0, 1.0]),
        array([-0.13614, 0.15714, -0.07165]) ]
    
    for v in vlist:
        ov = uc.calc_frac_to_orth(v)
        v2 = uc.calc_orth_to_frac(ov)
        print "----"
        print "    ",v
        print "    ",ov
        print "    ",v2
        print "----"


    print "================================================="

    print
    
    print "================================================="
    print "TEST CASE #2: Reciprocal of above unit cell "
    print

    ruc = uc.calc_reciprocal_unit_cell()
    print ruc
    print "volume      = ",ruc.calc_v()
    print "cell volume = ",ruc.calc_volume()
    
    print "================================================="


##     ## a more interesting space group
##     unitx = UnitCell(a=64.950,
##                      b=64.950,
##                      c=68.670,
##                      alpha=90.00,
##                      beta=90.00,
##                      gamma=120.00,
##                      space_group="P 32 2 1")


##     print "================================================="
##     print unitx
##     print

##     print unitx.frac_to_orth

##     for symop in unitx.space_group.iter_symops():
##         print "SYMOP:"
##         print symop
##         print "REAL:"
##         print
##         unitx.calc_rot_tra(symop)
##         print

if __name__=="__main__":
    main()

## </testing>
