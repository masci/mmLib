## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Mathmatical operations performed on mmLib.Strcuture.Atom objects.
"""
import math
from mmTypes import *

def length(u):
    """Calculates the length of u.
    """
    return sqrt(dot(u, u))

def normalize(u):
    """Returns the normalized vector along u.
    """
    return u/sqrt(dot(u, u))

def cross(u, v):
    """Cross product of u and v:
    Cross[u,v] = {-u3 v2 + u2 v3, u3 v1 - u1 v3, -u2 v1 + u1 v2}
    """
    return array([ u[1]*v[2] - u[2]*v[1],
                   u[2]*v[0] - u[0]*v[2],
                   u[0]*v[1] - u[1]*v[0] ], Float)

def rmatrix(alpha, beta, gamma):
    """Return a rotation matrix based on the Euler angles alpha,
    beta, and gamma in radians.
    """
    cosA = math.cos(alpha)
    cosB = math.cos(beta)
    cosG = math.cos(gamma)

    sinA = math.sin(alpha)
    sinB = math.sin(beta)
    sinG = math.sin(gamma)
    
    return array(
        [[cosB*cosG, cosG*sinA*sinB-cosA*sinG, cosA*cosG*sinB+sinA*sinG],
         [cosB*sinG, cosA*cosG+sinA*sinB*sinG, cosA*sinB*sinG-cosG*sinA ],
         [-sinB,     cosB*sinA,                cosA*cosB ]], Float)

def rmatrixu(u, theta):
    """Return a rotation matrix caused by a right hand rotation of theta
    radians around vector u.
    """
    u    = normalize(u)
    U    = array([ [   0.0, -u[2],  u[1] ],
                   [  u[2],   0.0, -u[0] ],
                   [ -u[1],  u[0],   0.0 ] ], Float)
        
    cosT = math.cos(theta)
    sinT = math.sin(theta)

    return identity(3, Float)*cosT + outerproduct(u,u)*(1-cosT) + U*sinT

def dmatrix(alpha, beta, gamma):
    """Returns the displacment matrix based on rotation about Euler
    angles alpha, beta, and gamma.
    """
    return rmatrix(alpha, beta, gamma) - identity(3, Float)

def dmatrixu(u, theta):
    """Return a displacement matrix caused by a right hand rotation of theta
    radians around vector u.
    """
    return rmatrixu(u, theta) - identity(3, Float)

def rmatrixz(u):
    """Return a rotation matrix which transforms the coordinate system
    such that the vector u is aligned along the z axis.
    """
    u, v, w = normalize(u)

    d = math.sqrt(u*u + v*v)

    if d!=0.0:
        Rxz = array([ [  u/d, v/d,  0.0 ],
                      [ -v/d, u/d,  0.0 ],
                      [  0.0, 0.0,  1.0 ] ], Float)
    else:
        Rxz = identity(3, Float)


    Rxz2z = array([ [   w, 0.0,    -d],
                    [ 0.0, 1.0,   0.0],
                    [   d, 0.0,     w] ], Float)

    return matrixmultiply(Rxz2z, Rxz)
	 
def calc_distance(a1, a2):
    """Returns the distance between two argument atoms.
    """
    if a1==None or a2==None:
        return None
    return length(a1.position - a2.position)

def calc_angle(a1, a2, a3):
    """Return the angle between the three argument atoms.
    """
    if a1==None or a2==None or a3==None:
        return None
    a21 = a1.position - a2.position
    a21 = a21 / (length(a21))

    a23 = a3.position - a2.position
    a23 = a23 / (length(a23))

    return arccos(dot(a21, a23))

def calc_torsion_angle(a1, a2, a3, a4):
    """Calculates the torsion angle between the four argument atoms.
    """
    if a1==None or a2==None or a3==None or a4==None:
        return None

    a12 = a2.position - a1.position
    a23 = a3.position - a2.position
    a34 = a4.position - a3.position

    n12 = cross(a12, a23)
    n34 = cross(a23, a34)

    n12 = n12 / length(n12)
    n34 = n34 / length(n34)

    cross_n12_n34  = cross(n12, n34)
    direction      = cross_n12_n34 * a23
    scalar_product = dot(n12, n34)

    if scalar_product > 1.0:
        scalar_product = 1.0
    if scalar_product < -1.0:
        scalar_product = -1.0

    angle = arccos(scalar_product)

    if direction < 0.0:
        angle = -angle

    return angle

def calc_CCuij(U, V):
    invU = inverse(U)
    invV = inverse(V)
    
    det_invU = determinant(invU)
    det_invV = determinant(invV)

    return ( math.sqrt(math.sqrt(det_invU * det_invV)) /
             math.sqrt((1.0/8.0) * determinant(invU + invV)) )

def calc_Suij(U, V):
    eqU = trace(U) / 3.0
    eqV = trace(V) / 3.0

    isoU = eqU * identity(3, Float)
    isoV = eqV * identity(3, Float)

    return ( calc_CCuij(U, (eqU/eqV)*V) /
             (calc_CCuij(U, isoU) * calc_CCuij(V, isoV)) )

def calc_DP2uij(U, V):
    invU = inverse(U)
    invV = inverse(V)

    det_invU = determinant(invU)
    det_invV = determinant(invV)

    pi3 = math.pi * math.pi * math.pi

    Pu2 = math.sqrt( det_invU / (64.0 * pi3) )
    Pv2 = math.sqrt( det_invV / (64.0 * pi3) )

    Puv = math.sqrt( (det_invU * det_invV) /
                     (8.0 * pi3 * determinant(invU + invV)) )

    return Pu2 + Pv2 - (2.0 * Puv)


### <TESTING>
if __name__ == "__main__":
    from Structure import *

    a1 = Atom(x=0.0, y=-1.0, z=0.0)
    a2 = Atom(x=0.0, y=0.0,  z=0.0)
    a3 = Atom(x=1.0, y=0.0,  z=0.0)
    a4 = Atom(x=1.0, y=1.0,  z=-1.0)

    print "a1:",a1.position
    print "calc_angle:",calc_angle(a1, a2, a3)
    print "calc_torsion_angle:",calc_torsion_angle(a1, a2, a3, a4)
### </TESTING>

