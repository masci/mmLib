## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Classes relating to space groups and symmetry.
"""
from mmTypes import *


class SymmetryOperator:

    mat_i    = identity(3)
    mat_inv  = -1.0 * identity(3)
    mat_2z   = array([[-1.0, 0.0, 0.0], [0.0,-1.0, 0.0], [0.0, 0.0, 1.0]])
    mat_3z   = array([[0.0,-1.0, 0.0], [1.0,-1.0, 0.0], [0.0, 0.0, 1.0]])
    mat_4z   = array([[0.0,-1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    mat_6z   = array([[1.0,-1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    mat_2q   = array([[0.0,-1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0,-1.0]])
    mat_2qq  = array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0,-1.0]])
    mat_3abc = array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])

    vec_0 = array([0.0, 0.0, 0.0])
    vec_a = array([0.5, 0.0, 0.0])
    vec_b = array([0.0, 0.5, 0.0])
    vec_c = array([0.0, 0.0, 0.5])
    vec_n = array([0.5, 0.5, 0.5])
    vec_u = array([0.25, 0.0, 0.0])
    vec_v = array([0.0, 0.25, 0.0])
    vec_w = array([0.0, 0.0, 0.25])
    vec_d = array([0.25, 0.25, 0.25])
    vec_A = array([0.0, 0.5, 0.5])
    vec_B = array([0.5, 0.0, 0.5])
    vec_C = array([0.5, 0.5, 0.0])
    vec_I = array([0.5, 0.5, 0.5])
    vec_R = array([0.666667, 0.333333, 0.333333])
    vec_S = array([0.333333, 0.666667, 0.333333])
    vec_T = array([0.333333, 0.333333, 0.666667])
    vec_H = array([0.666667, 0.333333, 0.0])
    vec_F1= array([0.0, 0.5, 0.5])
    vec_F2= array([0.5, 0.0, 0.5])

    def __init__(self, symop_desc):
        (symX, symY, symZ) = symop_desc.split(",")


        self.R = identity(3)
        self.t = Vector(0.0, 0.0, 0.0)

    def decode_symop(self, symop_desc):
        """Decode symmetry operator strings of the form: 1/2-X,-Y,1/2+Z
        """
        try:
            (symX, symY, symZ) = symop_desc.upper().split(",")
        except ValueError:
            return None

        def opsplit(x):
            tr = ""
            rot = ""
            for t in ["1/2", "1/4", "1/3", "2/3"]:
                if x.find(t) > -1:
                    tr = 


class SpaceGroup:
    pass
