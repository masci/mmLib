## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Please ignore this for now -- this code doesn't do anything useful yet.
"""
from __future__ import generators
import cmath

from mmTypes import *
from AtomMath import *
from UnitCell import UnitCell
from Structure import *


class StructureFactors(object):
    """Class for calculating crystallographic structure factors given
    a Structure model.  Please ignore this for now -- this code doesn't
    do anything useful yet.
    """
    def __init__(self, struct):
        self.struct = struct

    def iter_hkl(self, rh, rk, rl):
        """
        """
        h = -rh
        k = -rk
        l = -rl
        
        while h <= rh:
            while k <= rk:
                while l <= rl:
                    yield h, k, l
                    l += 1
                k += 1
                l = -rl
            h += 1
            k = -rk

    def iter_equiv_pos_frac(self, orth):
        """
        """
        unit_cell = struct.unit_cell
        xyz = unit_cell.calc_orth_to_frac(orth)

        for x, y, z in unit_cell.space_group.iter_equivalent_positions(xyz):
            if x < 0.0:
                x += 1.0
            elif x >= 1.0:
                x -= 1.0

            if y < 0.0:
                y += 1.0
            elif y >= 1.0:
                y -= 1.0

            if z < 0.0:
                z += 1.0
            elif z >= 1.0:
                z -= 1.0

            yield x, y, z

    def calc_Fcalc(self, h_range, k_range, l_range):
        """Returns a 3 dimentional hkl array of structure factors calculated
        from the negitive to positive hkl values of the given ranges.
        """
        two_pi_i = 2.0 * cmath.pi * 1.0j

        hkl = zeros((2*h_range + 1, 2*k_range + 1, 2*l_range + 1), Complex)

        for atm in self.struct.iter_atoms():
            for x, y, z in self.iter_equiv_pos_frac(atm.position):
                for h, k, l in self.iter_hkl(h_range, k_range, l_range):

                    F = cmath.exp(two_pi_i * (x*h + y*k + z*l))
                    hkl[h,k,l] += F

        return hkl

    def calc_R(self):
        """
        """
        pass


## <testing>
if __name__ == "__main__":
    from FileLoader import LoadStructure

    struct = LoadStructure(fil = sys.argv[1])
    for chain in struct.iter_chains():
        print chain

    print "Calculating Structure Factors..."
    
    sf  = StructureFactors(struct)
    hkl = sf.calc_Fcalc(10, 10, 0)

    for h, k, l in sf.iter_hkl(10, 10, 0):
        print "%4d  %4d  %4d  %40s  %25s" % (
            h, k, l,
            hkl[h,k,l],
            hkl[h,k,l] * hkl[h,k,l].conjugate())
## </testing>
