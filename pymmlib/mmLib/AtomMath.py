## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

"""Mathmatical operations performed on mmLib.Strcuture.Atom objects."""

import math
from   mmTypes import *

def calculateDistance(a1, a2):
    """Returns the distance between two argument atoms."""
    if not (a1 and a2): return None

    a12 = a1.position - a2.position
    distance = abs(a12.length())
    return distance


def calculateAngle(a1, a2, a3):
    """Return the angle between the three argument atoms."""
    if not (a1 and a2 and a3): return None

    a21 = a1.position - a2.position
    a23 = a3.position - a2.position
    return a21.angle(a23)


def calculateTorsionAngle(a1, a2, a3, a4):
    """Calculates the torsion angle between the four argument atoms."""
    if not (a1 and a2 and a3 and a4): return None

    a12 = a2.position - a1.position
    a23 = a3.position - a2.position
    a34 = a4.position - a3.position

    n12 = a12.cross(a23)
    n34 = a23.cross(a34)

    n12 = n12.normal()
    n34 = n34.normal()

    cross_n12_n34 = n12.cross(n34)
    direction = cross_n12_n34 * a23
    scalar_product = n12 * n34

    if scalar_product > 1.0:
        scalar_product = 1.0
    if scalar_product < -1.0:
        scalar_product = -1.0

    angle = math.acos(scalar_product)

    if direction < 0.0:
        angle = -angle

    return angle


##
## <testing>
##
if __name__ == "__main__":
    from Structure import *

    a1 = Atom()
    a1.setPosition(Vector(0.0,0.0,0.0))

    a2 = Atom()
    a2.setPosition(Vector(1.0,0.0,0.0))

    a3 = Atom()
    a3.setPosition(Vector(0.0,1.0,0.0))

    a4 = Atom()
    a4.setPosition(Vector(0.0,1.0,-1.0))

    print calculateTorsionAngle(a1, a2, a3, a4) * 180 / math.pi
##
## </testing>
##
