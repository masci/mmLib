## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Classes for representing biological macromolecules.
"""
from __future__ import generators
import cmath

from mmTypes import *
from AtomMath import *
from UnitCell import UnitCell
from Structure import *


class StructureFactors(object):

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

            yield xyz

    def calc(self, h_range, k_range, l_range):
        """
        """
        hkl = zeros((2*h_range + 1, 2*k_range + 1, 2*l_range + 1), Complex)

        for atm in self.struct.iter_atoms():
            for x, y, z in self.iter_equiv_pos_frac(atm.position):
                for h, k, l in self.iter_hkl(h_range, k_range, l_range):

                    F = cmath.exp(2.0 * cmath.pi * 1.0j * (x*h + y*k + z*l))
                    hkl[h,k,l] += F

        return hkl


if __name__ == "__main__":
    from FileLoader import LoadStructure

    struct = LoadStructure(fil = sys.argv[1])
    for atm in struct:
        print atm
    
    sf  = StructureFactors(struct)
    hkl = sf.calc(5, 5, 5)

    for h, k, l in sf.iter_hkl(3, 3, 0):
        print "%4d  %4d  %4d  %40s  %40s" % (
            h, k, l,
            hkl[h,k,l],
            hkl[h,k,l] * hkl[h,k,l].conjugate())


        
