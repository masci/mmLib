## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Custom solvers for systems of equations.
"""
import math
from mmTypes import *


class Solver(object):
    """Generic equation solving interface.
    """
    pass


class NewtonsMethod(Solver):
    """Newton's Method for solving a system of non-linear equations
    using a finite difference jacobian.  This algorithm has been modified
    to handle a overdetermined system of equations.
    """
    def __init__(self, F, N=100, h=1.0e-2, tol=1.0e-10):
        """Initalize with the list of non-linear functions to solve.  F
        is a Python list of function objects.  Each function object needs
        to implement a __call__ method which takes a single x vector
        (Numeric Python array()) as its argument.
        """
        self.F   = F
        self.N   = N
        self.h   = h
        self.tol = tol 

    def calc_F(self, x):
        """Returns the vector F(x).
        """
        F = zeros(len(self.F), Float)
        for i in range(len(self.F)):
            F[i] = self.F[i](x)
        return F
            
    def calc_finite_difference_jacobian(self, x):
        """Returns the matrix J(x) using finite difference derivitives.
        """
        J  = zeros((len(self.F), len(x)), Float)
        h  = self.h
        Eh = h * identity(len(x))

        for i in range(len(self.F)):
            f = self.F[i]
            for j in range(len(x)):
                J[i,j] = (f(x + Eh[j]) - f(x)) / h

        return J

    def solve(self, x):
        """Solves the system of nonlinear equations given the initial vector
        x.
        """
        x = x.copy()

        for n in range(self.N):
            #print "Cycle: ",n
            
            f = self.calc_F(x)
            J = self.calc_finite_difference_jacobian(x)

            (y, resids, rank, s) = linear_least_squares(J, -f)

            x += y

            if math.sqrt(dot(y, y)) < self.tol:
                break

        return x


## <testing>
if __name__ == '__main__':

    def f1(x):
        return 3.0*x[0] - math.cos(x[1]*x[2]) - 0.5
    def f2(x):
        return x[0]*x[0] - 81.0*((x[1]+0.1)**2) + math.sin(x[2]) + 1.06
    def f3(x):
        return math.exp(-x[0]*x[1]) + 20.0*x[2] + (10.0*math.pi - 3.0)/3.0

    solver = NewtonsMethod([f1, f2, f3])
    print solver.solve(array([0.1, 0.1, -0.1]))

## </testing>
