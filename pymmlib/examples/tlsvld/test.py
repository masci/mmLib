#!/usr/bin/python
import tlsvld
import scipy
import numpy
import math

U = tlsvld.array([[0.37471362,-0.07411912,-0.0379023],[-0.07411912,0.39010298,-0.03240841],[-0.0379023,-0.03240841,0.33057791]])
V = tlsvld.array([[0.38484156,-0.06933809,-0.04994903],[-0.06933809,0.41376374,-0.03275701],[-0.04994903,-0.03275701,0.33967632]])
#print U
#print V
#print "="*80
#U = scipy.array([[0.37471362,-0.07411912,-0.0379023],[-0.07411912,0.39010298,-0.03240841],[-0.0379023,-0.03240841,0.33057791]])
#V = scipy.array([[0.38484156,-0.06933809,-0.04994903],[-0.06933809,0.41376374,-0.03275701],[-0.04994903,-0.03275701,0.33967632]])
#V = scipy.array([[100.38484156,-0.06933809,-0.04994903],[-0.06933809,0.41376374,-0.03275701],[-0.04994903,-0.03275701,0.33967632]])
#print U
#print V
#print "="*80

#U = tlsvld.array([[7,7,7],[7,7,7],[7,7,7]])
#V = tlsvld.array([[4,4,4],[4,4,4],[4,4,4]])
W = tlsvld.array([[1.0,1.0,1.0],[1.0,1.0,1.0],[1.0,1.0,1.0]])
X = tlsvld.array([[4.0,4.0,4.0],[4.0,4.0,4.0],[4.0,4.0,4.0]])
print "SQD = ",tlsvld.sum_square_diff(W,X)

#print "trace(U) = ", tlsvld.trace3x3(X)
print "dtUV = ",tlsvld.diff_trace_uv(W,X)

print tlsvld.mul3x3(U,V)
print tlsvld.add3x3(U,V)
print "SUM(W) = ",tlsvld.sum3x3(X)
print tlsvld.subtract3x3(U,V)
#print tlsvld.trn3x3(U)
#print tlsvld.inv3x3(U)

#W = scipy.array([[0,0,0],[0,0,0],[0,0,0]])
#det, invU = tlsvld.inv3x3(U)
#print "PYTHON: det(U) = %s" % det
#print "PYTHON: inv(U) = %s" % invU
#print "="*80
#X = scipy.array([[1,0,0],[0,2,0],[0,0,3]])
#Y = scipy.array([[1,0,0],[0,2,0],[0,0,3]])
#print tlsvld.add3x3(X,Y)
#print tlsvld.calc_ccuij(X, Y)
#print "="*80

print "CCuij = %s" % tlsvld.calc_ccuij(U, V)

#print "EIGENS = ", tlsvld.jacobi(U,3,3, EIGENS, EVECS )
print U
print "EIGENS = ", tlsvld.eigenvalues(U)

x1 = 1.00
y1 = 2.00
z1 = 3.00
x2 = 6.00
y2 = 5.00
z2 = 4.00
d = math.sqrt( (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 )
U = tlsvld.array([[0.37471362,-0.07411912,-0.0379023],[-0.07411912,0.39010298,-0.03240841],[-0.0379023,-0.03240841,0.33057791]])
V = tlsvld.array([[0.38484156,-0.06933809,-0.04994903],[-0.06933809,0.41376374,-0.03275701],[-0.04994903,-0.03275701,0.33967632]])
print d
a = tlsvld.array([x1,y1,z1])
b = tlsvld.array([x2,y2,z2])

#     n = numpy.array([(a[0] - b[0])/d, (a[1] - b[1])/d, (a[2] - b[2])/d])
#     Un = numpy.dot(numpy.dot(n, U), numpy.transpose(n))
#     Vn = numpy.dot(numpy.dot(n, V), numpy.transpose(n))
#     return abs(Un - Vn)
n = numpy.array([(a[0] - b[0])/d, (a[1] - b[1])/d, (a[2] - b[2])/d])
Un = numpy.dot(numpy.dot(n, U), numpy.transpose(n))
Vn = numpy.dot(numpy.dot(n, V), numpy.transpose(n))
print "[Python] DOT(n, U) = ", numpy.dot(n, U)
print "[Python] d = ", d
print "[Python] Rosenfeld = ", abs(Un - Vn)
print "="*80
print "[Fortran] Rosenfeld = ", tlsvld.rosenfeld(a,b,U,V)
