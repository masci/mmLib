// Copyright 2006 by TLSMD Development Group (see AUTHORS file)
// This code is part of the TLSMD distribution and governed by
// its license.  Please see the LICENSE file that should have been
// included as part of this package.

#include <stdio.h>
#include <math.h>
#include "tls_model_nl.h"

#define PI     3.1415926535897931
#define PI2    (PI * PI)
#define PI3    (PI * PI * PI)

#define RAD2DEG  (180.0   / PI)
#define RAD2DEG2 (RAD2DEG*RAD2DEG)
#define DEG2RAD  (PI / 180.0)
#define DEG2RAD2 (DEG2RAD*DEG2RAD)

#define LSMALL (0.0001 * DEG2RAD2)

namespace TLSMD {

// prototype for MINPACK FORTRAN subroutine
typedef void (*FCN)(int*, int*, double*, double*, double*, int*, int*);
extern "C" void lmder1_(FCN, int*, int*, double*, double*, double*, int*, double*, int*, int*, double*, int*);

static ConstrainedFitIsotropicTLSModel *g_pISolver = 0;
static ConstrainedFitAnisotropicTLSModel *g_pASolver = 0;

// zero a M(m,n) double matrix
inline void
zero_dmatrix(double *M, int m, int n) {
  int sz = m * n;
  for (int i = 0; i < sz; ++i, ++M) *M = 0.0;
}

// calculate the 6x6 Jacobian of the positive definite tensor 
// x,y,z are the root-mean-squre values of the principal tensor components
// a,b,c are Alpha/Beta/Gamma convention Euler angles which orient the tensor
inline void
calc_pdtensor_derivs(double x, double y, double z, double a, double b, double c, double dT[6][6]) {
  double xx, yy, zz;
  double sa, sb, sc, ca, cb, cc;
  double sa2, sb2, sc2, ca2, cb2, cc2;
  double s2a, s2b, s2c, c2a, c2b, c2c;
  double tmp1;

  xx = x*x;
  yy = y*y;
  zz = z*z;

  sa = sin(a);
  sb = sin(b);
  sc = sin(c);

  ca = cos(a);
  cb = cos(b);
  cc = cos(c);

  sa2 = sa*sa;
  sb2 = sb*sb;
  sc2 = sc*sc;

  ca2 = ca*ca;
  cb2 = cb*cb;
  cc2 = cc*cc;

  s2a = sin(2.0*a);
  s2b = sin(2.0*b);
  s2c = sin(2.0*c);

  c2a = cos(2.0*a);
  c2b = cos(2.0*b);
  c2c = cos(2.0*c);

  /* dT11/dx */
  tmp1 = ca*cc - cb*sa*sc;
  dT[0][0] = 2.0*x*tmp1*tmp1;
  /* dT22/dx */
  tmp1 = cb*cc*sa + ca*sc;
  dT[1][0] = 2.0*x*tmp1*tmp1;
  /* dT33/dx */
  dT[2][0] = 2.0*x*sa2*sb2;
  /* dT12/dx */
  dT[3][0] = -0.25*x*(4.0*cb*c2c*s2a + (c2a*(3.0+c2b) + 2.0*sb2)*s2c); 
  /* dT13/dx */
  dT[4][0] = 2.0*x*sa*sb*(ca*cc - cb*sa*sc);
  /* dT23/dx */
  dT[5][0] = -2.0*x*sa*sb*(cb*cc*sa + ca*sc);

  /* dT11/dy */
  tmp1 = cc*sa + ca*cb*sc;
  dT[0][1] = 2.0*y*tmp1*tmp1;
  /* dT22/dy */
  tmp1 = ca*cb*cc - sa*sc;
  dT[1][1] = 2.0*y*tmp1*tmp1;
  /* dT33/dy */
  dT[2][1] = 2.0*y*ca2*sb2;
  /* dT12/dy */
  dT[3][1] = 0.25*y*(4.0*cb*c2c*s2a + (c2a*(3.0+c2b) - 2.0*sb2)*s2c);
  /* dT13/dy */
  dT[4][1] = -2.0*y*ca*sb*(cc*sa + ca*cb*sc);
  /* dT23/dy */
  dT[5][1] = 2.0*y*ca*sb*(-ca*cb*cc + sa*sc);

  /* dT11/dz */
  dT[0][2] = 2.0*z*sb2*sc2;
  /* dT22/dz */
  dT[1][2] = 2.0*z*cc2*sb2;
  /* dT33/dz */
  dT[2][2] = 2.0*z*cb2;
  /* dT12/dz */
  dT[3][2] = z*sb2*s2c;
  /* dT13/dz */
  dT[4][2] = z*s2b*sc;
  /* dT23/dz */
  dT[5][2] = z*cc*s2b;

  /* dT11/da */
  dT[0][3] = -2.0*(x-y)*(x+y)*(cc*sa + ca*cb*sc)*(ca*cc - cb*sa*sc);
  /* dT22/da */
  dT[1][3] = 2.0*(x-y)*(x+y)*(cb*cc*sa + ca*sc)*(ca*cb*cc - sa*sc);
  /* dT33/da */
  dT[2][3] = (x-y)*(x+y)*s2a*sb2;
  /* dT12/da */
  dT[3][3] = 0.25*(x-y)*(x+y)*(-4.0*c2a*cb*c2c + (3.0+c2b)*s2a*s2c);
  /* dT13/da */
  dT[4][3] = (x-y)*(x+y)*sb*(c2a*cc - cb*s2a*sc);
  /* dT23/da */
  dT[5][3] = -(x-y)*(x+y)*sb*(cb*cc*s2a + c2a*sc);

  /* dT11/db */
  dT[0][4] = 2.0*sb*sc*((x-y)*(x+y)*ca*cc*sa + cb*(zz - yy*ca2 - xx*sa2)*sc);
  /* dT22/db */
  dT[1][4] = 2.0*cc*sb*(cb*cc*(zz - yy*ca2 - xx*sa2) + (-xx + yy)*ca*sa*sc);
  /* dT33/db */
  dT[2][4] = 2.0*cb*(-zz + yy*ca2 + xx*sa2)*sb;
  /* dT12/db */
  dT[3][4] = 0.25*(2.0*(x-y)*(x+y)*c2c*s2a*sb + (-xx -yy + 2.0*zz + (x-y)*(x+y)*c2a)*s2b*s2c);
  /* dT13/db */
  dT[4][4] = (x-y)*(x+y)*ca*cb*cc*sa + 0.5*(-xx + 2.0*zz - 2.0*yy*ca2 + xx*c2a)*c2b*sc;
  /* dT23/db */
  dT[5][4] = c2b*cc*(zz - yy*ca2 - xx*sa2) + (-xx+yy)*ca*cb*sa*sc;

  /* dT11/dc */
  dT[0][5] = -(x-y)*(x+y)*cb*c2c*s2a - (xx + yy - 2.0*zz)*cc*sb2*sc - 0.25*(x-y)*(x+y)*c2a*(3.0+c2b)*s2c;
  /* dT22/dc */
  dT[1][5] =(x-y)*(x+y)*cb*c2c*s2a + (xx + yy - 2.0*zz)*cc*sb2*sc + 0.25*(x-y)*(x+y)*c2a*(3.0+c2b)*s2c;
  /* dT33/dc */
  dT[2][5] = 0.0;
  /* dT12/dc */
  dT[3][5] = 0.25*(c2c*(-(x-y)*(x+y)*c2a*(3.0+c2b) - 2.0*(xx + yy - 2.0*zz)*sb2) + 4.0*(x-y)*(x+y)*cb*s2a*s2c);
  /* dT13/dc */
  dT[4][5] = sb*(cb*cc*(zz - yy*ca2 - xx*sa2) + (-xx + yy)*ca*sa*sc);
  /* dT23/dc */
  dT[5][5] = sb*(-(x-y)*(x+y)*ca*cc*sa + cb*(-zz + yy*ca2 + xx*sa2)*sc);
}

// calculate the 3x3 Gaussian covariance tensor given
// x,y,z are the root-mean-squre values of the principal tensor components
// a,b,c are Alpha/Beta/Gamma convention Euler angles which orient the tensor
inline void
calc_pdtensor(double x, double y, double z, double a, double b, double c, double T[6]) {
  double xx, yy, zz;
  double sa, sb, sc, ca, cb, cc;
  double sa2, sb2, sc2, ca2, cb2, cc2;
  double s2a, s2b, s2c, c2a, c2b, c2c;
  double tmp1, tmp2;

  xx = x*x;
  yy = y*y;
  zz = z*z;

  sa = sin(a);
  sb = sin(b);
  sc = sin(c);

  ca = cos(a);
  cb = cos(b);
  cc = cos(c);

  sa2 = sa*sa;
  sb2 = sb*sb;
  sc2 = sc*sc;

  ca2 = ca*ca;
  cb2 = cb*cb;
  cc2 = cc*cc;

  s2a = sin(2.0*a);
  s2b = sin(2.0*b);
  s2c = sin(2.0*c);

  c2a = cos(2.0*a);
  c2b = cos(2.0*b);
  c2c = cos(2.0*c);

  /* T11 */
  tmp1 = cc*sa + ca*cb*sc;
  tmp2 = ca*cc - cb*sa*sc;
  T[0] = zz*sb2*sc2 + yy*tmp1*tmp1 + xx*tmp2*tmp2;

  /* T22 */
  tmp1 = cb*cc*sa + ca*sc;
  tmp2 = ca*cb*cc - sa*sc;
  T[1] = zz*cc2*sb2 + xx*tmp1*tmp1 + yy*tmp2*tmp2;
  
  /* T33 */
  T[2] = zz*cb2 + (yy*ca2 + xx*sa2)*sb2;

  /* T12 */
  T[3] = 0.125*(-4.0*(x-y)*(x+y)*cb*c2c*s2a + (-(x-y)*(x+y)*c2a*(3.0 + c2b) - 2.0*(xx + yy - 2.0*zz)*sb2) * s2c);
  
  /* T13 */
  T[4] = sb*((x-y)*(x+y)*ca*cc*sa + cb*(zz - yy*ca2 - xx*sa2)*sc);

  /* T23 */
  T[5] = sb*(cb*cc*(zz - yy*ca2 - xx*sa2) + (-xx + yy)*ca*sa*sc);
}

/* calculate linear anisotropic TLS parameters from non-linear TLS parameters */
inline void
calc_isotropic_tls_parameters(double NL_ITLS[NL_ITLS_NUM_PARAMS], double ITLS[ITLS_NUM_PARAMS]) {
  for (int i = 0; i < ITLS_NUM_PARAMS; ++i) {
    ITLS[i] = NL_ITLS[i];
  }

  double lx, ly, lz;
  lx = fabs(NL_ITLS[NL_ITLS_LX]);
  ly = fabs(NL_ITLS[NL_ITLS_LY]);
  lz = fabs(NL_ITLS[NL_ITLS_LZ]);

  if ((lx*lx) < LSMALL) {
    lx = 0.0;
  }
  if ((ly*ly) < LSMALL) {
    ly = 0.0;
  }
  if ((lz*lz) < LSMALL) {
    lz = 0.0;
  }

  calc_pdtensor(lx, ly, lz, NL_ITLS[NL_ITLS_LA], NL_ITLS[NL_ITLS_LB], NL_ITLS[NL_ITLS_LC], &ITLS[ITLS_L11]);
}

/* calculate linear anisotropic TLS parameters from non-linear TLS parameters */
inline void
calc_anisotropic_tls_parameters(double NL_ATLS[ATLS_NUM_PARAMS], double ATLS[ATLS_NUM_PARAMS]) {
  int i;
  double lx, ly, lz;

  for (i = 0; i < ATLS_NUM_PARAMS; i++) {
    ATLS[i] = NL_ATLS[i];
  }

  lx = fabs(NL_ATLS[NL_ATLS_LX]);
  ly = fabs(NL_ATLS[NL_ATLS_LY]);
  lz = fabs(NL_ATLS[NL_ATLS_LZ]);

  if ((lx*lx)<LSMALL) {
    lx = 0.0;
  }
  if ((ly*ly)<LSMALL) {
    ly = 0.0;
  }
  if ((lz*lz)<LSMALL) {
    lz = 0.0;
  }

  calc_pdtensor(lx, ly, lz, NL_ATLS[NL_ATLS_LA], NL_ATLS[NL_ATLS_LB], NL_ATLS[NL_ATLS_LC], &ATLS[ATLS_L11]);
}

inline void
set_isotropic_jacobian(double *A, int m, int n, int row, double x, double y, double z, double dNL[6][6]) {
#define FA(__i, __j) A[__i + (m * __j)]

  double xx, yy, zz, xy, xz, yz;
   
  xx = x*x;
  yy = y*y;
  zz = z*z;
  xy = x*y;
  xz = x*z;
  yz = y*z;

  /* T iso */
  FA(row, NL_ITLS_T) = 1.0;

  /* Uiso = 1/3((zz+yy)L11 + 1/3(xx+zz)L22 + 1/3(xx+yy)L33 -2xyL12 -2xzL13 -2yzL23) */
  FA(row, NL_ITLS_LX) = ((zz+yy)*dNL[0][0] + (xx+zz)*dNL[1][0] + (xx+yy)*dNL[2][0] - 2.0*xy*dNL[3][0] - 2.0*xz*dNL[4][0] - 2.0*yz*dNL[5][0])/3.0;
  FA(row, NL_ITLS_LY) = ((zz+yy)*dNL[0][1] + (xx+zz)*dNL[1][1] + (xx+yy)*dNL[2][1] - 2.0*xy*dNL[3][1] - 2.0*xz*dNL[4][1] - 2.0*yz*dNL[5][1])/3.0;
  FA(row, NL_ITLS_LZ) = ((zz+yy)*dNL[0][2] + (xx+zz)*dNL[1][2] + (xx+yy)*dNL[2][2] - 2.0*xy*dNL[3][2] - 2.0*xz*dNL[4][2] - 2.0*yz*dNL[5][2])/3.0;
  FA(row, NL_ITLS_LA) = ((zz+yy)*dNL[0][3] + (xx+zz)*dNL[1][3] + (xx+yy)*dNL[2][3] - 2.0*xy*dNL[3][3] - 2.0*xz*dNL[4][3] - 2.0*yz*dNL[5][3])/3.0;
  FA(row, NL_ITLS_LB) = ((zz+yy)*dNL[0][4] + (xx+zz)*dNL[1][4] + (xx+yy)*dNL[2][4] - 2.0*xy*dNL[3][4] - 2.0*xz*dNL[4][4] - 2.0*yz*dNL[5][4])/3.0;
  FA(row, NL_ITLS_LC) = ((zz+yy)*dNL[0][5] + (xx+zz)*dNL[1][5] + (xx+yy)*dNL[2][5] - 2.0*xy*dNL[3][5] - 2.0*xz*dNL[4][5] - 2.0*yz*dNL[5][5])/3.0;

  FA(row, NL_ITLS_S1) =  (2.0/3.0)*z;
  FA(row, NL_ITLS_S2) =  (2.0/3.0)*y;
  FA(row, NL_ITLS_S3) =  (2.0/3.0)*x;

#undef FA
}


// Sets the six rows of matrix A starting at A[i,j] with the anistropic
// TLS model coefficents for a atom located at t position x, y, z with
// least-squares weight w.  Matrix A is filled to coumn j+12.
// 
// The matrix A(m,n) is filled in FORTRAN-style, that is, assuming
// a column-major memory layout.  That's because, athough the
// matrix is contructed in C, the SVD subroutine from LAPACK is 
// written in FORTRAN.
inline void
set_anisotropic_jacobian(double *A, int m, int n, int row, double x, double y, double z, double dNL[6][6]) {
#define FA(__i,__j) A[__i + (m * __j)]

  int iU11, iU22, iU33, iU12, iU13, iU23;
  double xx, yy, zz, xy, xz, yz;
  
  xx = x*x;
  yy = y*y;
  zz = z*z;
  xy = x*y;
  xz = x*z;
  yz = y*z;

  /* calculate row indexes */
  iU11 = row;
  iU22 = row + 1;
  iU33 = row + 2;
  iU12 = row + 3;
  iU13 = row + 4;
  iU23 = row + 5;

  /* set A */

  /* U11 = T11 + zzL22 + yyL33 - 2yzL23 + 2zS21 - 2yS31 */
  FA(iU11, ATLS_T11) = 1.0;

  FA(iU11, ATLS_L11) = zz*dNL[1][0] + yy*dNL[2][0] - 2.0*yz*dNL[5][0];
  FA(iU11, ATLS_L22) = zz*dNL[1][1] + yy*dNL[2][1] - 2.0*yz*dNL[5][1];
  FA(iU11, ATLS_L33) = zz*dNL[1][2] + yy*dNL[2][2] - 2.0*yz*dNL[5][2];
  FA(iU11, ATLS_L12) = zz*dNL[1][3] + yy*dNL[2][3] - 2.0*yz*dNL[5][3];
  FA(iU11, ATLS_L13) = zz*dNL[1][4] + yy*dNL[2][4] - 2.0*yz*dNL[5][4];
  FA(iU11, ATLS_L23) = zz*dNL[1][5] + yy*dNL[2][5] - 2.0*yz*dNL[5][5];

  FA(iU11, ATLS_S31) = -2.0 * y;
  FA(iU11, ATLS_S21) =  2.0 * z;

  /* U22 = T22 + zzL11 + xxL33 - 2xzL13 -2zS12 + 2xS32 */
  FA(iU22, ATLS_T22) = 1.0;

  FA(iU22, ATLS_L11) = zz*dNL[0][0] + xx*dNL[2][0] - 2.0*xz*dNL[4][0];
  FA(iU22, ATLS_L22) = zz*dNL[0][1] + xx*dNL[2][1] - 2.0*xz*dNL[4][1];
  FA(iU22, ATLS_L33) = zz*dNL[0][2] + xx*dNL[2][2] - 2.0*xz*dNL[4][2];
  FA(iU22, ATLS_L12) = zz*dNL[0][3] + xx*dNL[2][3] - 2.0*xz*dNL[4][3];
  FA(iU22, ATLS_L13) = zz*dNL[0][4] + xx*dNL[2][4] - 2.0*xz*dNL[4][4];
  FA(iU22, ATLS_L23) = zz*dNL[0][5] + xx*dNL[2][5] - 2.0*xz*dNL[4][5];

  FA(iU22, ATLS_S12) = -2.0 * z;
  FA(iU22, ATLS_S32) =  2.0 * x;

  /* U33 = T33 + yyL11 + xxL22 - 2xyL12 - 2xS23 + 2yS13 */
  FA(iU33, ATLS_T33) = 1.0;

  FA(iU33, ATLS_L11) = yy*dNL[0][0] + xx*dNL[1][0] - 2.0*xy*dNL[3][0];
  FA(iU33, ATLS_L22) = yy*dNL[0][1] + xx*dNL[1][1] - 2.0*xy*dNL[3][1];
  FA(iU33, ATLS_L33) = yy*dNL[0][2] + xx*dNL[1][2] - 2.0*xy*dNL[3][2];
  FA(iU33, ATLS_L12) = yy*dNL[0][3] + xx*dNL[1][3] - 2.0*xy*dNL[3][3];
  FA(iU33, ATLS_L13) = yy*dNL[0][4] + xx*dNL[1][4] - 2.0*xy*dNL[3][4];
  FA(iU33, ATLS_L23) = yy*dNL[0][5] + xx*dNL[1][5] - 2.0*xy*dNL[3][5];

  FA(iU33, ATLS_S23) = -2.0 *  x;
  FA(iU33, ATLS_S13) =  2.0 *  y;

  /* U12 = T12 - xyL33 + xzL23 + yzL13 - zzL12 + zS2211 + xS31 - yS32   */
  FA(iU12, ATLS_T12)   = 1.0;

  FA(iU12, ATLS_L11)   = -xy*dNL[2][0] + xz*dNL[5][0] + yz*dNL[4][0] - zz*dNL[3][0];
  FA(iU12, ATLS_L22)   = -xy*dNL[2][1] + xz*dNL[5][1] + yz*dNL[4][1] - zz*dNL[3][1];
  FA(iU12, ATLS_L33)   = -xy*dNL[2][2] + xz*dNL[5][2] + yz*dNL[4][2] - zz*dNL[3][2];
  FA(iU12, ATLS_L12)   = -xy*dNL[2][3] + xz*dNL[5][3] + yz*dNL[4][3] - zz*dNL[3][3];
  FA(iU12, ATLS_L13)   = -xy*dNL[2][4] + xz*dNL[5][4] + yz*dNL[4][4] - zz*dNL[3][4];
  FA(iU12, ATLS_L23)   = -xy*dNL[2][5] + xz*dNL[5][5] + yz*dNL[4][5] - zz*dNL[3][5];

  FA(iU12, ATLS_S2211) =   z;
  FA(iU12, ATLS_S31)   =   x;
  FA(iU12, ATLS_S32)   =  -y;
  
  /* U13 = T13 - xzL22 + xyL23 - yyL13 + yzL12 + yS1133 + zS23 - xS21 */
  FA(iU13, ATLS_T13)   = 1.0;

  FA(iU13, ATLS_L11)   = -xz*dNL[1][0] + xy*dNL[5][0] - yy*dNL[4][0] + yz*dNL[3][0];
  FA(iU13, ATLS_L22)   = -xz*dNL[1][1] + xy*dNL[5][1] - yy*dNL[4][1] + yz*dNL[3][1];
  FA(iU13, ATLS_L33)   = -xz*dNL[1][2] + xy*dNL[5][2] - yy*dNL[4][2] + yz*dNL[3][2];
  FA(iU13, ATLS_L12)   = -xz*dNL[1][3] + xy*dNL[5][3] - yy*dNL[4][3] + yz*dNL[3][3];
  FA(iU13, ATLS_L13)   = -xz*dNL[1][4] + xy*dNL[5][4] - yy*dNL[4][4] + yz*dNL[3][4];
  FA(iU13, ATLS_L23)   = -xz*dNL[1][5] + xy*dNL[5][5] - yy*dNL[4][5] + yz*dNL[3][5];

  FA(iU13, ATLS_S1133) =   y;
  FA(iU13, ATLS_S23)   =   z;
  FA(iU13, ATLS_S21)   =  -x;
    
  /* U23 = T23 - yzL11 - xxL23 +xyL13 +xzL12 - xS2211 - xS1133 + yS12 - zS13  */
  FA(iU23, ATLS_T23)   = 1.0;

  FA(iU23, ATLS_L11)   = -yz*dNL[0][0] - xx*dNL[5][0] + xy*dNL[4][0] + xz*dNL[3][0];
  FA(iU23, ATLS_L22)   = -yz*dNL[0][1] - xx*dNL[5][1] + xy*dNL[4][1] + xz*dNL[3][1];
  FA(iU23, ATLS_L33)   = -yz*dNL[0][2] - xx*dNL[5][2] + xy*dNL[4][2] + xz*dNL[3][2];
  FA(iU23, ATLS_L12)   = -yz*dNL[0][3] - xx*dNL[5][3] + xy*dNL[4][3] + xz*dNL[3][3];
  FA(iU23, ATLS_L13)   = -yz*dNL[0][4] - xx*dNL[5][4] + xy*dNL[4][4] + xz*dNL[3][4];
  FA(iU23, ATLS_L23)   = -yz*dNL[0][5] - xx*dNL[5][5] + xy*dNL[4][5] + xz*dNL[3][5];

  FA(iU23, ATLS_S2211) =  -x;
  FA(iU23, ATLS_S1133) =  -x;
  FA(iU23, ATLS_S12)   =   y;
  FA(iU23, ATLS_S13)   =  -z;

#undef FA
}

// calculate Jacobian matrix
inline void
anisotropic_lmder1_fcn_jacobian(int m, int n,double *NL_ATLS, double *A, int lda) {
  double dNL[6][6];
  calc_pdtensor_derivs(NL_ATLS[NL_ATLS_LX], NL_ATLS[NL_ATLS_LY], NL_ATLS[NL_ATLS_LZ], 
		       NL_ATLS[NL_ATLS_LA], NL_ATLS[NL_ATLS_LB], NL_ATLS[NL_ATLS_LC], 
		       dNL);

  zero_dmatrix(A, m, n);

  int row = 0;
  std::vector<AnisotropicDataPoint>::const_iterator idp;
  for (idp = g_pASolver->adata_vector.begin(); idp != g_pASolver->adata_vector.end(); ++idp) {
    set_anisotropic_jacobian(A, m, n, row, idp->x, idp->y, idp->z, dNL);
    row += 6;
  }
}


// calculate function residuals
void
anisotropic_lmder1_fcn_R(int m, int n, double *NL_ATLS, double *R) {
  double ATLS[ATLS_NUM_PARAMS];
  calc_anisotropic_tls_parameters(NL_ATLS, ATLS);

  int row = 0;
  std::vector<AnisotropicDataPoint>::const_iterator idp;
  for (idp = g_pASolver->adata_vector.begin(); idp != g_pASolver->adata_vector.end(); ++idp) {
    double U[6];
    CalcAnisotropicTLSModelU(ATLS, idp->x, idp->y, idp->z, U);
    for (int j = 0; j < 6; ++j, ++row) {
      R[row] = U[j] - idp->U[j];
    }
  }
}

// callback for lmder1
// 
// calculates residual function for m equations and store 
// in R[i]
//
// calculates Jacobian and stores values in FORTRAN 2D
// array A(i,j) where i is the equation and j is the parameter
// fcn(     m,      n,         x,          fvec,      fjac,    ldfjac,    iflag)
static void 
anisotropic_lmder1_fcn(int *m, int *n, double *NL_ATLS, double *R, double *A, int *lda, int *iflag) {
  if (*iflag == 1) {
    /* calculate f(x) = fvec */
    anisotropic_lmder1_fcn_R(*m, *n, NL_ATLS, R);
    /* prnt_R(R, *m); */

  } else if (*iflag == 2) {
    /* calculate J(x) = fjac */
    anisotropic_lmder1_fcn_jacobian(*m, *n, NL_ATLS, A, *lda);
    /* prnt_A(A, *m, *n, NL_ATLS); */
  }
}

// calculate isotropic TLS model Jacobian matrix
inline void
isotropic_lmder1_fcn_jacobian(int m, int n,double *NL_ITLS, double *A, int lda) {
  double dNL[6][6];
  calc_pdtensor_derivs(NL_ITLS[NL_ITLS_LX], NL_ITLS[NL_ITLS_LY], NL_ITLS[NL_ITLS_LZ], 
		       NL_ITLS[NL_ITLS_LA], NL_ITLS[NL_ITLS_LB], NL_ITLS[NL_ITLS_LC], 
		       dNL);

  zero_dmatrix(A, m, n);

  int row = 0;
  std::vector<IsotropicDataPoint>::const_iterator idp;
  for (idp = g_pISolver->idata_vector.begin(); idp != g_pISolver->idata_vector.end(); ++idp, ++row) {
    set_isotropic_jacobian(A, m, n, row, idp->x, idp->y, idp->z, dNL);
  }
}

// calculate isotropic TLS model function residuals
inline void
isotropic_lmder1_fcn_R(int m, int n, double *NL_ITLS, double *R) {
  double ITLS[ITLS_NUM_PARAMS];
  calc_isotropic_tls_parameters(NL_ITLS, ITLS);

  int row = 0;
  std::vector<IsotropicDataPoint>::const_iterator idp;
  for (idp = g_pISolver->idata_vector.begin(); idp != g_pISolver->idata_vector.end(); ++idp, ++row) {
    double uiso;
    CalcIsotropicTLSModelUIso(ITLS, idp->x, idp->y, idp->z, &uiso);
    R[row] = uiso - idp->uiso;
  }
}

// callback for lmder1
//
// calculates residual function for m equations and store 
// in R[i]
//
// calculates Jacobian and stores values in FORTRAN 2D
// array A(i,j) where i is the equation and j is the parameter
// fcn(     m,      n,         x,          fvec,      fjac,    ldfjac,    iflag)
void 
isotropic_lmder1_fcn(int *m, int *n, double *NL_ITLS, double *R, double *A, int *lda, int *iflag) {
  if (*iflag == 1) {
    // calculate f(x) = fvec
    isotropic_lmder1_fcn_R(*m, *n, NL_ITLS, R);
  } else if (*iflag == 2) {
    // calculate J(x) = fjac
    isotropic_lmder1_fcn_jacobian(*m, *n, NL_ITLS, A, *lda);
  }
}

ConstrainedFitTLSModel::ConstrainedFitTLSModel()
  : num_rows(0), 
    num_cols(0), 
    fvec(0), 
    fjac(0), 
    ipvt(0), 
    lwa(0), 
    wa(0),
    max_num_atoms(0),
    iatom(0),
    tls_model(0) {
}

ConstrainedFitTLSModel::~ConstrainedFitTLSModel() {
  delete[] fvec;
  delete[] fjac;
  delete[] ipvt;
  delete[] wa;
}

void
ConstrainedFitTLSModel::set_size(int nrows, int ncols) {
  delete[] fvec;
  delete[] fjac;
  delete[] ipvt;
  delete[] wa;

  num_rows = nrows;
  num_cols = ncols;

  fvec = new double[num_rows];
  fjac = new double[num_rows * num_cols];
  ipvt = new int[num_cols];

  lwa = num_rows + 5 * num_cols;
  wa = new double[lwa];
}

void
ConstrainedFitIsotropicTLSModel::set_max_num_atoms(int num_atoms) {
  max_num_atoms = num_atoms;
  set_size(num_atoms, ITLS_NUM_PARAMS);
  idata_vector.reserve(num_atoms);
}

void
ConstrainedFitIsotropicTLSModel::reset_fit(TLSModel *tls_model, int num_atoms) {
  if (num_atoms > max_num_atoms) set_max_num_atoms(num_atoms);
  ConstrainedFitIsotropicTLSModel::tls_model = tls_model;
  iatom = 0;
  idata_vector.resize(num_atoms);
}

void
ConstrainedFitIsotropicTLSModel::set_data_point(double x, double y, double z, double uiso) {
  idata_vector[iatom].x = x - tls_model->origin_x;
  idata_vector[iatom].y = y - tls_model->origin_y;
  idata_vector[iatom].z = z - tls_model->origin_z;
  idata_vector[iatom].uiso = uiso;
  ++iatom;
}

void
ConstrainedFitIsotropicTLSModel::fit_params() {
  double mean_u_iso = 0.0;
  std::vector<IsotropicDataPoint>::const_iterator idp;
  for (idp = idata_vector.begin(); idp != idata_vector.end(); ++idp) {
    mean_u_iso += idp->uiso;
  }
  mean_u_iso = mean_u_iso / idata_vector.size();

  double NL_ITLS[ITLS_NUM_PARAMS];
  for (int i = 0; i < ITLS_NUM_PARAMS; ++i) NL_ITLS[i] = 0.0;
  NL_ITLS[NL_ITLS_T] = mean_u_iso;
  NL_ITLS[NL_ITLS_LX] = 5.0 * DEG2RAD;
  NL_ITLS[NL_ITLS_LY] = 5.0 * DEG2RAD;
  NL_ITLS[NL_ITLS_LZ] = 5.0 * DEG2RAD;

  g_pISolver = this;
  int info;
  int num_equations = idata_vector.size();
  int num_variables = ITLS_NUM_PARAMS;
  double tol = tolerance;
  lmder1_(isotropic_lmder1_fcn, &num_equations, &num_variables, NL_ITLS, fvec, fjac, &num_equations, &tol, &info, ipvt, wa, &lwa);
  g_pISolver = 0;

  calc_isotropic_tls_parameters(NL_ITLS, tls_model->get_params());
}

void
ConstrainedFitAnisotropicTLSModel::set_max_num_atoms(int num_atoms) {
  max_num_atoms = num_atoms;
  set_size(6 * num_atoms, ATLS_NUM_PARAMS);
  adata_vector.reserve(num_atoms);
}

void
ConstrainedFitAnisotropicTLSModel::reset_fit(TLSModel *tls_model, int num_atoms) {
  if (num_atoms > max_num_atoms) set_max_num_atoms(num_atoms);
  ConstrainedFitAnisotropicTLSModel::tls_model = tls_model;
  iatom = 0;
  adata_vector.resize(num_atoms);
}

void
ConstrainedFitAnisotropicTLSModel::set_data_point(double x, double y, double z, const double U[]) {
  adata_vector[iatom].x = x - tls_model->origin_x;
  adata_vector[iatom].y = y - tls_model->origin_y;
  adata_vector[iatom].z = z - tls_model->origin_z;
  for (int i = 0; i < 6; ++i) adata_vector[iatom].U[i] = U[i];
  ++iatom;
}

void
ConstrainedFitAnisotropicTLSModel::fit_params() {
  double mean_u_iso = 0.0;
  std::vector<AnisotropicDataPoint>::const_iterator idp;
  for (idp = adata_vector.begin(); idp != adata_vector.end(); ++idp) {
    mean_u_iso += (idp->U[0] + idp->U[1] + idp->U[2]) / 3.0;
  }
  mean_u_iso = mean_u_iso / adata_vector.size();

  double NL_ATLS[ATLS_NUM_PARAMS];
  for (int i = 0; i < ATLS_NUM_PARAMS; ++i) NL_ATLS[i] = 0.0;
  NL_ATLS[NL_ATLS_T11] = 1.0 * mean_u_iso;
  NL_ATLS[NL_ATLS_T22] = 1.0 * mean_u_iso;
  NL_ATLS[NL_ATLS_T33] = 1.0 * mean_u_iso;
  NL_ATLS[NL_ATLS_LX] = 5.0 * DEG2RAD;
  NL_ATLS[NL_ATLS_LY] = 5.0 * DEG2RAD;
  NL_ATLS[NL_ATLS_LZ] = 5.0 * DEG2RAD;

  g_pASolver = this;
  int info;
  int num_equations = 6 * adata_vector.size();
  int num_variables = ATLS_NUM_PARAMS;
  double tol = tolerance;
  lmder1_(anisotropic_lmder1_fcn, &num_equations, &num_variables, NL_ATLS, fvec, fjac, &num_equations, &tol, &info, ipvt, wa, &lwa);
  g_pASolver = 0;

  calc_anisotropic_tls_parameters(NL_ATLS, tls_model->get_params());
}

} // namespace TLSMD
