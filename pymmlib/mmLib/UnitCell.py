## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Classes for handling unit cell transformation.
"""
import math
from AtomMath import *


class UnitCell(object):
    """Class for storing and performing calculations on unit cell
    parameters.
    """
    def __init__(self, a=1.0, b=1.0, c=1.0, alpha=90.0, beta=90.0, gamma=90.0):
        self.a     = a
        self.b     = b
        self.c     = c
        self.alpha = alpha
        self.beta  = beta
        self.gamma = gamma

    def __str__(self):
        return "UnitCell(a=%f, b=%f, c=%f, alpha=%f, beta=%f, gamma=%f)" % (
            self.a, self.b, self.c, self.alpha, self.beta, self.gamma)

    def calc_v(self):
        """Calculates the volume of the rhombohedrial created by the
        unit vectors a1/|a1|, a2/|a2|, a3/|a3|.
        """
        cos_alpha = cos(self.alpha)
        cos_beta  = cos(self.beta)
        cos_gamma = cos(self.gamma)
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

        sin_alpha = sin(self.alpha)
        sin_beta  = sin(self.beta)
        sin_gamma = sin(self.gamma)

        cos_alpha = cos(self.alpha)
        cos_beta  = cos(self.beta)
        cos_gamma = cos(self.gamma)

        ra  = (self.b * self.c * sin_alpha) / V
        rb  = (self.a * self.c * sin_beta)  / V
        rc  = (self.a * self.b * sin_gamma) / V

        ralpha = acos(
            (cos_beta  * cos_gamma - cos_alpha) / (sin_beta  * sin_gamma))
        rbeta  = acos(
            (cos_alpha * cos_gamma - cos_beta)  / (sin_alpha * sin_gamma))
        rgamma = acos(
            (cos_alpha * cos_beta  - cos_gamma) / (sin_alpha * sin_beta))

        return UnitCell(ra, rb, rc, ralpha, rbeta, rgamma)

    def calc_fractionalization_matrix(self):
        """Cartesian to fractional coordinates.
        """
        return transpose(self.calc_orthogonalization_matrix())
        
    def calc_orthogonalization_matrix(self):
        """Fractional to cartesian coordinates.
        """
        sin_alpha = sin(self.alpha)
        sin_beta  = sin(self.beta)
        sin_gamma = sin(self.gamma)

        cos_alpha = cos(self.alpha)
        cos_beta  = cos(self.beta)
        cos_gamma = cos(self.gamma)

        o11 = self.a
        o12 = self.b * cos_gamma
        o13 = self.c * cos_beta
        o21 = 0.0
        o22 = self.b * sin_gamma
        o23 = (self.c * (cos_alpha - cos_beta*cos_gamma)) / sin_gamma
        o31 = 0.0
        o32 = 0.0
        o33 = (self.c * self.calc_v()) / sin_gamma

        return array([[o11, o12, o13],
                      [o21, o22, o23],
                      [o31, o32, o33]])

    def calc_cartisian_unit_cell_axes(self):
        return matrixmultiply(
            self.calc_orthogonalization_matrix(),
            array([[1.0, 0.0, 0.0],
                   [0.0, 1.0, 0.0],
                   [0.0, 0.0, 1.0]]))


if __name__ == "__main__":
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
