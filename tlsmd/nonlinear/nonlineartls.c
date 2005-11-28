/* nonlineartls.c 
 * Jay Painter <jpaint@u.washington.edu>
 * Sept 7, 2005
 *
 * Copyright 2002 by TLSMD Development Group (see AUTHORS file)
 * This code is part of the TLSMD distribution and governed by
 * its license.  Please see the LICENSE file that should have been
 * included as part of this package.
 *
 * Implementation of a Python module for the constrained non-linear 
 * fitting of TLS parameters given a list of atoms with 
 * crystallographically refined ADPs.  Uses MINPACK.
 */
#include "Python.h"
#include "structmember.h"

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <pthread.h>


/* LAPACK */
extern void dgesdd_(char *, int*, int*, double*, int*, double*, double*, int *, double*, int*, double*, int*, int*, int*);

/* prototype for MINPACK FORTRAN subroutine */
typedef void (*FCN)();
extern void lmder1_(FCN, int*, int*, double*, double*, double*, int*, double*, int*, int*, double*, int*);

/* macros */
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))

/* set to 1 for extra debugging informaiton */
/* #define _DEBUG 0 */


/* Python Interface: Exception class for this module */
static PyObject *NONLINEARTLS_ERROR = NULL;

/*
 * Mathmatical Constants
 */

#define PI     3.1415926535897931
#define PI2    (PI * PI)
#define PI3    (PI * PI * PI)

#define RAD2DEG  (180.0   / PI)
#define RAD2DEG2 (RAD2DEG*RAD2DEG)
#define DEG2RAD  (PI / 180.0)
#define DEG2RAD2 (DEG2RAD*DEG2RAD)

#define U2B ((8.0*PI2))
#define B2U ((1.0/U2B))

#define LSMALL (0.0001 * DEG2RAD2)


/*
 * Data Structures
 */

/* anisotropic U tensor parameter labels and indexes */
#define U11 0
#define U22 1
#define U33 2
#define U12 3
#define U13 4
#define U23 5
#define U_NUM_PARAMS 6
static char *U_PARAM_NAMES[] = {
  "u11", "u22", "u33", "u12", "u13", "u23"
};

/* isotropic TLS Model 
 * S1 == S21-S12; S2 == S13-S31; S3 == S32-S23
 */
#define ITLS_T     0
#define ITLS_L11   1
#define ITLS_L22   2
#define ITLS_L33   3
#define ITLS_L12   4
#define ITLS_L13   5
#define ITLS_L23   6
#define ITLS_S1    7
#define ITLS_S2    8
#define ITLS_S3    9
#define ITLS_NUM_PARAMS 10
static char *ITLS_PARAM_NAMES[] = {
  "it", 
  "il11", "il22", "il33", "il12", "il13", "il23", 
  "is1", "is2", "is3"
};

/* anisotropic TLS model parameter indexes and labels */
#define ATLS_T11   0
#define ATLS_T22   1
#define ATLS_T33   2
#define ATLS_T12   3
#define ATLS_T13   4
#define ATLS_T23   5
#define ATLS_L11   6
#define ATLS_L22   7
#define ATLS_L33   8
#define ATLS_L12   9
#define ATLS_L13   10
#define ATLS_L23   11
#define ATLS_S2211 12
#define ATLS_S1133 13
#define ATLS_S12   14
#define ATLS_S13   15
#define ATLS_S23   16
#define ATLS_S21   17
#define ATLS_S31   18
#define ATLS_S32   19
#define ATLS_NUM_PARAMS 20
static char *ATLS_PARAM_NAMES[] = {
  "t11", "t22", "t33", "t12", "t13", "t23",
  "l11", "l22", "l33", "l12", "l13", "l23",
  "s2211", "s1133", "s12", "s13", "s23", "s21", "s31", "s32"
};

/* isotropic non-linear TLS Model */
#define NL_ITLS_T     0
#define NL_ITLS_LX    1
#define NL_ITLS_LY    2
#define NL_ITLS_LZ    3
#define NL_ITLS_LA    4
#define NL_ITLS_LB    5
#define NL_ITLS_LC    6
#define NL_ITLS_S1    7
#define NL_ITLS_S2    8
#define NL_ITLS_S3    9
#define NL_ITLS_NUM_PARAMS 10
static char *NL_ITLS_PARAM_NAMES[] = {
  "nl_t",
  "nl_lx", "nl_ly", "nl_lz", "nl_la", "nl_lb", "nl_lc",
  "nl_s1", "nl_s2", "nl_s3"
};

/* anisotropic non-linear TLS model parameter indexes and labels */
#define NL_ATLS_T11   0
#define NL_ATLS_T22   1
#define NL_ATLS_T33   2
#define NL_ATLS_T12   3
#define NL_ATLS_T13   4
#define NL_ATLS_T23   5
#define NL_ATLS_LX    6
#define NL_ATLS_LY    7
#define NL_ATLS_LZ    8
#define NL_ATLS_LA    9
#define NL_ATLS_LB    10
#define NL_ATLS_LC    11
#define NL_ATLS_S2211 12
#define NL_ATLS_S1133 13
#define NL_ATLS_S12   14
#define NL_ATLS_S13   15
#define NL_ATLS_S23   16
#define NL_ATLS_S21   17
#define NL_ATLS_S31   18
#define NL_ATLS_S32   19
#define NL_ATLS_NUM_PARAMS 20
static char *NL_ATLS_PARAM_NAMES[] = {
  "nl_t11", "nl_t22", "nl_t33", "nl_t12", "nl_t13", "nl_t23",
  "nl_lx", "nl_ly", "nl_lz", "nl_la", "nl_lb", "nl_lc",
  "nl_s2211", "nl_s1133", "nl_s12", "nl_s13",  "nl_s23", "nl_s21", "nl_s31", "nl_s32"
};

/* Atom Structure */
#define NAME_LEN     8
#define FRAG_ID_LEN  8

struct Atom {
  char    name[NAME_LEN];
  char    frag_id[FRAG_ID_LEN];
  double  x, y, z;
  double  weight;
  double  sqrt_weight;
  double  xtls, ytls, ztls;
  double  u_iso;
  double  u_iso_tmp;
  double  U[6];
  double  Utmp[6];
};

/* TLS fit data structure used by both anisotropic 
 * and isotropic TLS models 
 */
struct TLSFit {
  int                istart;
  int                iend;
  struct Atom       *atoms;

  int                num_atoms;
  int                m;
  int                n;

  double             ox;
  double             oy;
  double             oz;

  double             NL_ATLS[NL_ATLS_NUM_PARAMS];
  double             NL_ITLS[NL_ITLS_NUM_PARAMS];

  double             ATLS[ATLS_NUM_PARAMS];
  double             ITLS[ITLS_NUM_PARAMS];

  double             alsqr;
  double             alsqr_mainchain;
  double             ilsqr;
  double             ilsqr_mainchain;
};

struct IHingeContext {
  struct Atom       *atoms;                 /* pointer to Atom array */

  int                istarta;               /* index of first atom in chain->atoms */
  int                ienda;                 /* index of last atom in chain->atoms */
  int                istartb;               /* index of first atom in chain->atoms */
  int                iendb;                 /* index of last atom in chain->atoms */

  double             msd_a;
  double             msd_b;
  double             msd_c;

  double             hdelta_ab;             /* hinge delta residual */
  double             hdelta_abo;            /* hinge delta residual */
  double             hdelta_c;
};

/* global pointer to current minimization problem */
struct TLSFit *g_pFit = NULL;


/* zero a M(m,n) double matrix */
inline void
zero_dmatrix(double *M, int m, int n)
{
  int i, sz;

  sz = m * n;
  for (i = 0; i < sz; i++) {
    M[i] = 0.0;
  }
}

/* return 1 if the atom is a mainchain atom */
inline int
atom_is_mainchain(struct Atom *atoms, int ia)
{
  if (strcmp(atoms[ia].name,"N")!=0 && strcmp(atoms[ia].name,"CA")!=0 && strcmp(atoms[ia].name,"C")!=0 && strcmp(atoms[ia].name,"0")!=0) {
    return 1;
  }

  return 0;
}

/* normalize the vector v
 */
inline void 
normalize(double v[3])
{
  double d;

  d = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  v[0] = v[0] / d;
  v[1] = v[1] / d;
  v[2] = v[2] / d;
}

/* compute the length of the vector v
 */
inline double
length(double v[3])
{
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

/* form the cross product of vectors u and v, and return the
 * result in w
 */
inline void 
cross(double u[3], double v[3], double w[3])
{
  w[0] = u[1]*v[2] - u[2]*v[1];
  w[1] = u[2]*v[0] - u[0]*v[2];
  w[2] = u[0]*v[1] - u[1]*v[0];
}

/* calculate the determinate of a symmetric 3x3 matrix
 */
inline double
det_symmetric_3(double U[6])
{
  return - U[4]*U[4]*U[1] + 2.0*U[3]*U[4]*U[5] - U[0]*U[5]*U[5] - U[3]*U[3]*U[2] + U[0]*U[1]*U[2];
}

/* matrix multiply symmetric 3x3 matrix U with V and return the result in M 
 * matrix format: u11,u22,u33,u12,u13,u23
 */
inline void
mult_symmetric_3(double U[6], double V[6], double M[6])
{
  int i;

  for (i = 0; i < 6; i++) {
    M[i] = 0.0;
  }

  M[0] = U[0]*V[0] + U[3]*V[3] + U[4]*V[4];
  M[1] = U[3]*V[3] + U[1]*V[1] + U[5]*V[5];
  M[2] = U[4]*V[4] + U[5]*V[5] + U[2]*V[2];
  M[3] = U[0]*V[3] + U[3]*V[1] + U[4]*V[5];
  M[4] = U[0]*V[4] + U[3]*V[5] + U[4]*V[2];
  M[5] = U[3]*V[4] + U[1]*V[5] + U[5]*V[2];
}


/* invert the symmetric 3x3 matrix in the argument U, and return the result in Ui
 * matrix format: u11,u22,u33,u12,u13,u23
 */
static int
invert_symmetric_3(double U[6], double Ui[6])
{
  double d;

  /* calculate determinate */
  d = - U[4]*U[4]*U[1] + 2.0*U[3]*U[4]*U[5] - U[0]*U[5]*U[5]
    - U[3]*U[3]*U[2] + U[0]*U[1]*U[2];

  if (d == 0.0) {
    return 0;
  }

  Ui[0] = (-U[5]*U[5] + U[1]*U[2]) / d;
  Ui[1] = (-U[4]*U[4] + U[0]*U[2]) / d;
  Ui[2] = (-U[3]*U[3] + U[0]*U[1]) / d;
  Ui[3] = ( U[4]*U[5] - U[3]*U[2]) / d;
  Ui[4] = (-U[4]*U[1] + U[3]*U[5]) / d;
  Ui[5] = ( U[3]*U[4] - U[0]*U[5]) / d;

  return 1;
}

/* calculates the value of the dP^2(U,V) function, which is the
 * square of the volumetric difference in the two trivariate 
 * Gaussian probability distribution functions defined by the U and
 * V symmetric variance-covariance tensors; this calculation requires
 * inverting U and V matrixes; if this fails, then this funtion returns
 * 0, otherwise it returns 1 with the dp2 value in the return parameter
 *
 */
inline int
calc_dp2(double U[6], double V[6], double *dp2)
{
  double detUi, detVi, detW, pU2, pV2, pUV;
  double Ui[6], Vi[6], W[6];

  if (!invert_symmetric_3(U, Ui)) {
    return 0;
  }
  if (!invert_symmetric_3(V, Vi)) {
    return 0;
  }

  W[0] = Ui[0] + Vi[0];
  W[1] = Ui[1] + Vi[1];
  W[2] = Ui[2] + Vi[2];
  W[3] = Ui[3] + Vi[3];
  W[4] = Ui[4] + Vi[4];
  W[5] = Ui[5] + Vi[5];

  detUi = det_symmetric_3(Ui);
  if (detUi <= 0.0) {
    return 0;
  }

  detVi = det_symmetric_3(Vi);
  if (detVi <= 0.0) {
    return 0;
  }

  detW = det_symmetric_3(W);
  if (detW <= 0.0) {
    return 0;
  }

  pU2 = sqrt(detUi / (64.0 * PI3));
  pV2 = sqrt(detVi / (64.0 * PI3));
  pUV = sqrt((detUi * detVi) / (8.0 * PI3 * detW));

  *dp2 = pU2 + pV2 - (2.0 * pUV);
  return 1;
}

/* calculates the centroid of the atoms indexed between istart and iend
 * then returns the centroid coordinates in x, y, z
 */
static void
calc_centroid(struct Atom *atoms, int istart, int iend, double *x, double *y, double *z)
{
  int ia;
  double n, cx, cy, cz;

  n  = 0.0;
  cx = 0.0;
  cy = 0.0;
  cz = 0.0;

  for (ia = istart; ia <= iend; ia++) {
    n  += 1.0;
    cx += atoms[ia].x;
    cy += atoms[ia].y;
    cz += atoms[ia].z;
  }

  if (n>0.0) {
    *x = cx / n;
    *y = cy / n;
    *z = cz / n;
  } else {
    *x = 0.0;
    *y = 0.0;
    *z = 0.0;
  }
}

/* calculate the 6x6 Jacobian of the positive definite tensor 
 * x,y,z are the root-mean-squre values of the principal tensor components
 * a,b,c are Alpha/Beta/Gamma convention Euler angles which orient the tensor
 */
inline void
calc_pdtensor_derivs(double x, double y, double z, double a, double b, double c, double dT[6][6])
{
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

/* calculate the 3x3 Gaussian covariance tensor given
 * x,y,z are the root-mean-squre values of the principal tensor components
 * a,b,c are Alpha/Beta/Gamma convention Euler angles which orient the tensor
 */
inline void
calc_pdtensor(double x, double y, double z, double a, double b, double c, double T[6])
{
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


/* calculates u_iso for a atom at position x,y,z relative to the origin of the 
 * isotropic TLS mode ITLS
 */
inline void
calc_isotropic_uiso(double *ITLS, double x, double y, double z, double *u_iso)
{
  double xx, yy, zz;

  xx = x*x;
  yy = y*y;
  zz = z*z;

  /* note: S1 == S21-S12; S2 == S13-S31; S3 == S32-S23 */

  *u_iso = ITLS[ITLS_T]                     +
           (       ITLS[ITLS_L11] * (zz+yy) + 
	           ITLS[ITLS_L22] * (xx+zz) +
	           ITLS[ITLS_L33] * (xx+yy) +
	    -2.0 * ITLS[ITLS_L12] * x*y     +
	    -2.0 * ITLS[ITLS_L13] * x*z     +
            -2.0 * ITLS[ITLS_L23] * y*z     +
	     2.0 * ITLS[ITLS_S1]  * z       +
	     2.0 * ITLS[ITLS_S2]  * y       +
             2.0 * ITLS[ITLS_S3]  * x )/3.0;
}

/* return the anisotropic TLS model predicted ADP in U for a atom 
 * located at coordinates x,y,z with respect to the ATLS origin
 */
inline void
calc_anisotropic_Utls(double *ATLS, double x, double y, double z, double U[6]) 
{
  double xx, yy, zz, xy, yz, xz;

  xx = x * x;
  yy = y * y;
  zz = z * z;
  xy = x * y;
  yz = y * z;
  xz = x * z;
    
  U[U11] =         
            ATLS[ATLS_T11]
    +       ATLS[ATLS_L22] * zz
    +       ATLS[ATLS_L33] * yy
    - 2.0 * ATLS[ATLS_L23] * yz
    - 2.0 * ATLS[ATLS_S31] * y
    + 2.0 * ATLS[ATLS_S21] * z;

  U[U22] =
            ATLS[ATLS_T22]
    +       ATLS[ATLS_L11] * zz
    +       ATLS[ATLS_L33] * xx
    - 2.0 * ATLS[ATLS_L13] * xz
    - 2.0 * ATLS[ATLS_S12] * z
    + 2.0 * ATLS[ATLS_S32] * x;

  U[U33] =
            ATLS[ATLS_T33]
    +       ATLS[ATLS_L11] * yy
    +       ATLS[ATLS_L22] * xx
    - 2.0 * ATLS[ATLS_L12] * xy
    - 2.0 * ATLS[ATLS_S23] * x
    + 2.0 * ATLS[ATLS_S13] * y;

  U[U12] =
            ATLS[ATLS_T12]
    -       ATLS[ATLS_L33]   * xy
    +       ATLS[ATLS_L23]   * xz
    +       ATLS[ATLS_L13]   * yz
    -       ATLS[ATLS_L12]   * zz
    +       ATLS[ATLS_S2211] * z
    +       ATLS[ATLS_S31]   * x
    -       ATLS[ATLS_S32]   * y;

  U[U13] =
            ATLS[ATLS_T13]
    -       ATLS[ATLS_L22]   * xz
    +       ATLS[ATLS_L23]   * xy
    -       ATLS[ATLS_L13]   * yy
    +       ATLS[ATLS_L12]   * yz
    +       ATLS[ATLS_S1133] * y
    +       ATLS[ATLS_S23]   * z
    -       ATLS[ATLS_S21]   * x;

  U[U23] =
            ATLS[ATLS_T23]
    -       ATLS[ATLS_L11]   * yz
    -       ATLS[ATLS_L23]   * xx
    +       ATLS[ATLS_L13]   * xy
    +       ATLS[ATLS_L12]   * xz
    -      (ATLS[ATLS_S2211] + ATLS[ATLS_S1133]) * x
    +       ATLS[ATLS_S12]   * y
    -       ATLS[ATLS_S13]   * z;
}

/* calculate the least squares residual of the isotropic TLS model */
inline double
calc_isotropic_lsqr(struct TLSFit *fit)
{
  int ia;
  double delta, lsqr, u_iso;
  struct Atom *atom;

  lsqr = 0.0;

  for (ia = g_pFit->istart; ia <= g_pFit->iend; ia++) {
    atom = &g_pFit->atoms[ia];

    calc_isotropic_uiso(fit->ITLS, atom->xtls, atom->ytls, atom->ztls, &u_iso);
    delta = u_iso - atom->u_iso;
    lsqr += delta * delta;
  }
  
  return lsqr;
}

/* calculate the least squares residual of the anisotropic TLS model */
inline double
calc_anisotropic_lsqr(struct TLSFit *fit)
{
  int i, ia;
  double delta, lsqr;
  double U[6];
  struct Atom *atom;

  lsqr = 0.0;

  for (ia = g_pFit->istart; ia <= g_pFit->iend; ia++) {
    atom = &g_pFit->atoms[ia];

    calc_anisotropic_Utls(fit->ATLS, atom->xtls, atom->ytls, atom->ztls, U);

    for (i = 0; i < 6; i++) {
      delta = U[i] - atom->U[i];
      lsqr += delta*delta;
    }
  }
  
  return lsqr;
}

/* calculate linear anisotropic TLS parameters from non-linear TLS parameters */
inline void
calc_isotropic_tls_parameters(double NL_ITLS[NL_ITLS_NUM_PARAMS], double ITLS[ITLS_NUM_PARAMS])
{
  int i;
  double lx, ly, lz;

  for (i = 0; i < ITLS_NUM_PARAMS; i++) {
    ITLS[i] = NL_ITLS[i];
  }

  ITLS[ITLS_T] = NL_ITLS[NL_ITLS_T] * NL_ITLS[NL_ITLS_T];

  lx = fabs(NL_ITLS[NL_ITLS_LX]);
  ly = fabs(NL_ITLS[NL_ITLS_LY]);
  lz = fabs(NL_ITLS[NL_ITLS_LZ]);

  if ((lx*lx)<LSMALL) {
    lx = 0.0;
  }
  if ((ly*ly)<LSMALL) {
    ly = 0.0;
  }
  if ((lz*lz)<LSMALL) {
    lz = 0.0;
  }

  calc_pdtensor(lx, ly, lz, NL_ITLS[NL_ITLS_LA], NL_ITLS[NL_ITLS_LB], NL_ITLS[NL_ITLS_LC], &ITLS[ITLS_L11]);
}

/* calculate linear anisotropic TLS parameters from non-linear TLS parameters */
inline void
calc_anisotropic_tls_parameters(double NL_ATLS[ATLS_NUM_PARAMS], double ATLS[ATLS_NUM_PARAMS])
{
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
set_isotropic_jacobian(double *A, int m, int n, int row, double x, double y, double z, double Trmsd, double dNL[6][6])
{
#define FA(__i, __j) A[__i + (m * __j)]

  double xx, yy, zz, xy, xz, yz;
   
  xx = x*x;
  yy = y*y;
  zz = z*z;
  xy = x*y;
  xz = x*z;
  yz = y*z;

  /* T iso */
  FA(row, NL_ITLS_T) = 2.0 * Trmsd;

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

/* Sets the six rows of matrix A starting at A[i,j] with the anistropic
 * TLS model coefficents for a atom located at t position x, y, z with
 * least-squares weight w.  Matrix A is filled to coumn j+12.
 *
 * The matrix A(m,n) is filled in FORTRAN-style, that is, assuming
 * a column-major memory layout.  That's because, athough the
 * matrix is contructed in C, the SVD subroutine from LAPACK is 
 * written in FORTRAN.
 */
inline void
set_anisotropic_jacobian(double *A, int m, int n, int row, double x, double y, double z, double dNL[6][6])
{
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

/* debugging nonlinear residual */
void
prnt_R(double *R, int m)
{
  int i;
  printf("=============================================================================================================\n");
  for (i = 0; i < m; i++) {
    printf("%f\n", R[i]);
  }
  printf("=============================================================================================================\n");
}

/* debugging nonlinear jacobian and comparison with numeric jacobian */
void 
prnt_A(double *A, int m, int n, double *NL_ATLS1)
{
#define FA(__i,__j) A[__i + (m * __j)]
#define FNA(__i,__j) NA[__i + (m * __j)]

  int i, j, ia, row, p;
  double h, NL_ATLS2[NL_ATLS_NUM_PARAMS], ATLS1[ATLS_NUM_PARAMS], ATLS2[ATLS_NUM_PARAMS], U1[6], U2[6], *NA;
  struct Atom *atom;

  h = 1E-8;

  calc_anisotropic_tls_parameters(NL_ATLS1, ATLS1);

  NA = (double *) malloc(sizeof(double) * m * n);

  row = 0;
  for (ia = g_pFit->istart; ia <= g_pFit->iend; ia++) {

    atom = &g_pFit->atoms[ia];
    calc_anisotropic_Utls(ATLS1, atom->xtls, atom->ytls, atom->ztls, U1);

    for (j = 0; j < ATLS_NUM_PARAMS; j++) {
	
      /* copy parameters */
      for (p = 0; p < ATLS_NUM_PARAMS; p++) {
	NL_ATLS2[p] = NL_ATLS1[p];
      }
      NL_ATLS2[j] += h;

      calc_anisotropic_tls_parameters(NL_ATLS2, ATLS2);
      calc_anisotropic_Utls(ATLS2, atom->xtls, atom->ytls, atom->ztls, U2);
	
      for (i = 0; i < 6; i++) {
	FNA(row+i, j) = (U2[i] - U1[i])/h;
      }
    }
    row += 6;
  }

  printf("=============================================================================================================\n");

  for (i = 0; i < m; i++) {

    printf("A");
    for (j = 0; j < n; j++) {
      printf("%8.4f ", FA(i, j));
    }
    printf("\n");

    printf("N");
    for (j = 0; j < n; j++) {
      printf("%8.4f ", FNA(i, j));
    }
    printf("\n");
  }

  printf("=============================================================================================================\n");

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      FA(i, j) = FNA(i, j);
    }
  }

  free(NA);

#undef FA
#undef FNA
}

/* calculate Jacobian matrix */
inline void
anisotropic_lmder1_fcn_jacobian(int m, int n,double *NL_ATLS, double *A, int lda)
{
  int row, ia, istart, iend;
  double dNL[6][6];
  struct Atom *atom;

  istart = g_pFit->istart;
  iend   = g_pFit->iend;

  calc_pdtensor_derivs(NL_ATLS[NL_ATLS_LX], NL_ATLS[NL_ATLS_LY], NL_ATLS[NL_ATLS_LZ], 
		       NL_ATLS[NL_ATLS_LA], NL_ATLS[NL_ATLS_LB], NL_ATLS[NL_ATLS_LC], 
		       dNL);

  zero_dmatrix(A, m, n);

  row = 0;
  for (ia = istart; ia <= iend; ia++) {
    atom = &g_pFit->atoms[ia];
    set_anisotropic_jacobian(A, m, n, row, atom->xtls, atom->ytls, atom->ztls, dNL);
    row += 6;
  }
}

/* calculate function residuals */
inline void
anisotropic_lmder1_fcn_R(int m, int n, double *NL_ATLS, double *R)
{
  int i, j, ia, istart, iend;
  double ATLS[ATLS_NUM_PARAMS], U[6];
  struct Atom *atom;

  istart = g_pFit->istart;
  iend   = g_pFit->iend;

  calc_anisotropic_tls_parameters(NL_ATLS, ATLS);

  i = 0;
  for (ia = istart; ia <= iend; ia++) {
    atom = &g_pFit->atoms[ia];
    calc_anisotropic_Utls(ATLS, atom->xtls, atom->ytls, atom->ztls, U);
    for (j = 0; j < 6; j++) {
      R[i] = atom->sqrt_weight * (U[j] - atom->U[j]);
      i++;
    }
  }
}

/* callback for lmder1
 *
 * calculates residual function for m equations and store 
 * in R[i]
 *
 * calculates Jacobian and stores values in FORTRAN 2D
 * array A(i,j) where i is the equation and j is the parameter
 *
 *     fcn(     m,      n,         x,          fvec,      fjac,    ldfjac,    iflag)
 */
static void 
anisotropic_lmder1_fcn(int *m, int *n, double *NL_ATLS, double *R, double *A, int *lda, int *iflag)
{
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

static int
anisotropic_nonlinear_fit(struct TLSFit *fit)
{
  int i, ia;
  double mean_u_iso;

  int info, lwa, *ipvt;
  double tol, *fvec, *fjac, *wa;

  /* number of atoms to fit */
  fit->num_atoms = fit->iend - fit->istart + 1;

  /* number of equations */
  fit->m = 6 * fit->num_atoms;

  /* number of parameters */
  fit->n = NL_ATLS_NUM_PARAMS; 

  /* allocate the memory LMDIF1 needs for working space */
  lwa   = 5*fit->n + fit->m;
  fvec  = (double *) malloc(sizeof(double) * fit->m);
  fjac  = (double *) malloc(sizeof(double) * fit->m * fit->n);

  ipvt  = (int *) malloc(sizeof(int) * fit->n);

  lwa   = 5*fit->n + fit->m;
  wa    = (double *) malloc(sizeof(double) * lwa);

  /* initalize the parameters with a initial estimate */
  for (i = 0; i < NL_ATLS_NUM_PARAMS; i++) {
    fit->NL_ATLS[i] = 0.0;
  }

  /* calculate the mean Uiso */
  mean_u_iso = 0.0;
  for (ia = fit->istart; ia <= fit->iend; ia++) {
    mean_u_iso += fit->atoms[ia].u_iso;
  }
  mean_u_iso = mean_u_iso / fit->num_atoms;

  fit->NL_ATLS[NL_ATLS_T11] = 1.0 * mean_u_iso;
  fit->NL_ATLS[NL_ATLS_T22] = 1.0 * mean_u_iso;
  fit->NL_ATLS[NL_ATLS_T33] = 1.0 * mean_u_iso;

  fit->NL_ATLS[NL_ATLS_LX] = 5.0 * DEG2RAD;
  fit->NL_ATLS[NL_ATLS_LY] = 5.0 * DEG2RAD;
  fit->NL_ATLS[NL_ATLS_LZ] = 5.0 * DEG2RAD;

  /* calculate centroid for use as the TLS origin, then 
   * set atom positions relative to the TLS origin 
   */
  calc_centroid(fit->atoms, fit->istart, fit->iend, &fit->ox, &fit->oy, &fit->oz);
  
  for (ia = fit->istart; ia <= fit->iend; ia++) {
    fit->atoms[ia].xtls = fit->atoms[ia].x - fit->ox;
    fit->atoms[ia].ytls = fit->atoms[ia].y - fit->oy;
    fit->atoms[ia].ztls = fit->atoms[ia].z - fit->oz;
  }

  /* perform minimization */
  g_pFit = fit;
  tol = 1E-3;

  /* lmder1(fcn,                  m,       n,       x,            fvec, fjac, ldfjac,  tol,  info,  ipvt, wa, lwa) */
  lmder1_(anisotropic_lmder1_fcn, &fit->m, &fit->n, fit->NL_ATLS, fvec, fjac, &fit->m, &tol, &info, ipvt, wa, &lwa);

  /* calculate linear TLS parameters from non linear TLS parameters */
  calc_anisotropic_tls_parameters(fit->NL_ATLS, fit->ATLS);

  /* calculate least squares residual */
  fit->alsqr = calc_anisotropic_lsqr(fit);

  /* free working memory */
  free(fvec);
  free(fjac);
  free(ipvt);
  free(wa);

  return info;
}

/* calculate isotropic TLS model Jacobian matrix */
inline void
isotropic_lmder1_fcn_jacobian(int m, int n,double *NL_ITLS, double *A, int lda)
{
  int row, ia, istart, iend;
  double dNL[6][6];
  struct Atom *atom;

  istart = g_pFit->istart;
  iend   = g_pFit->iend;

  calc_pdtensor_derivs(NL_ITLS[NL_ITLS_LX], NL_ITLS[NL_ITLS_LY], NL_ITLS[NL_ITLS_LZ], 
		       NL_ITLS[NL_ITLS_LA], NL_ITLS[NL_ITLS_LB], NL_ITLS[NL_ITLS_LC], 
		       dNL);

  zero_dmatrix(A, m, n);

  row = 0;
  for (ia = istart; ia <= iend; ia++) {
    atom = &g_pFit->atoms[ia];
    set_isotropic_jacobian(A, m, n, row, atom->xtls, atom->ytls, atom->ztls, NL_ITLS[NL_ITLS_T], dNL);
    row += 1;
  }
}

/* calculate isotropic TLS model function residuals */
inline void
isotropic_lmder1_fcn_R(int m, int n, double *NL_ITLS, double *R)
{
  int i, ia, istart, iend;
  double ITLS[ITLS_NUM_PARAMS], u_iso_tls;
  struct Atom *atom;

  istart = g_pFit->istart;
  iend   = g_pFit->iend;

  calc_isotropic_tls_parameters(NL_ITLS, ITLS);

  i = 0;
  for (ia = istart; ia <= iend; ia++) {
    atom = &g_pFit->atoms[ia];
    calc_isotropic_uiso(ITLS, atom->xtls, atom->ytls, atom->ztls, &u_iso_tls);

    R[i] = atom->sqrt_weight * (u_iso_tls - atom->u_iso);
    i++;
  }
}

/* callback for lmder1
 *
 * calculates residual function for m equations and store 
 * in R[i]
 *
 * calculates Jacobian and stores values in FORTRAN 2D
 * array A(i,j) where i is the equation and j is the parameter
 *
 *     fcn(     m,      n,         x,          fvec,      fjac,    ldfjac,    iflag)
 */
static void 
isotropic_lmder1_fcn(int *m, int *n, double *NL_ITLS, double *R, double *A, int *lda, int *iflag)
{
  if (*iflag == 1) {
    /* calculate f(x) = fvec */
    isotropic_lmder1_fcn_R(*m, *n, NL_ITLS, R);
  } else if (*iflag == 2) {
    /* calculate J(x) = fjac */
    isotropic_lmder1_fcn_jacobian(*m, *n, NL_ITLS, A, *lda);
  }
}

static int
isotropic_nonlinear_fit(struct TLSFit *fit)
{
  int i, ia;
  double mean_u_iso;

  int info, lwa, *ipvt;
  double tol, *fvec, *fjac, *wa;

  /* number of atoms to fit */
  fit->num_atoms = fit->iend - fit->istart + 1;

  /* number of equations */
  fit->m = fit->num_atoms;

  /* number of parameters */
  fit->n = NL_ITLS_NUM_PARAMS; 

  /* allocate the memory LMDIF1 needs for working space */
  lwa   = 5*fit->n + fit->m;
  fvec  = (double *) malloc(sizeof(double) * fit->m);
  fjac  = (double *) malloc(sizeof(double) * fit->m * fit->n);

  ipvt  = (int *) malloc(sizeof(int) * fit->n);

  lwa   = 5*fit->n + fit->m;
  wa    = (double *) malloc(sizeof(double) * lwa);

  /* initalize the parameters with a initial estimate */
  for (i = 0; i < NL_ITLS_NUM_PARAMS; i++) {
    fit->NL_ITLS[i] = 0.0;
  }

  /* calculate the mean Uiso */
  mean_u_iso = 0.0;
  for (ia = fit->istart; ia <= fit->iend; ia++) {
    mean_u_iso += fit->atoms[ia].u_iso;
  }
  mean_u_iso = mean_u_iso / fit->num_atoms;

  fit->NL_ITLS[NL_ITLS_T] = mean_u_iso;

  fit->NL_ITLS[NL_ITLS_LX] = 5.0 * DEG2RAD;
  fit->NL_ITLS[NL_ITLS_LY] = 5.0 * DEG2RAD;
  fit->NL_ITLS[NL_ITLS_LZ] = 5.0 * DEG2RAD;

  /* calculate centroid for use as the TLS origin, then 
   * set atom positions relative to the TLS origin 
   */
  calc_centroid(fit->atoms, fit->istart, fit->iend, &fit->ox, &fit->oy, &fit->oz);
  
  for (ia = fit->istart; ia <= fit->iend; ia++) {
    fit->atoms[ia].xtls = fit->atoms[ia].x - fit->ox;
    fit->atoms[ia].ytls = fit->atoms[ia].y - fit->oy;
    fit->atoms[ia].ztls = fit->atoms[ia].z - fit->oz;
  }

  /* perform minimization */
  g_pFit = fit;
  tol = 1E-3;

  /* lmder1(fcn,                m,       n,       x,            fvec, fjac, ldfjac,  tol,  info,  ipvt, wa, lwa) */
  lmder1_(isotropic_lmder1_fcn, &fit->m, &fit->n, fit->NL_ITLS, fvec, fjac, &fit->m, &tol, &info, ipvt, wa, &lwa);

  /* calculate linear TLS parameters from non linear TLS parameters */
  calc_isotropic_tls_parameters(fit->NL_ITLS, fit->ITLS);

  /* calculate least squares residual */
  fit->ilsqr = calc_isotropic_lsqr(fit);

  /* free working memory */
  free(fvec);
  free(fjac);
  free(ipvt);
  free(wa);

  return info;
}

static void
calc_isotropic_hinge_delta(struct IHingeContext *hinge)
{
  int i, num_atoms, ia;
  double u_iso_a, u_iso_b;
  double delta, dela, delb;
  double hdelta_ab, hdelta_abo;
  double hdelta_a, hdelta_b;
  double msd_a, msd_b;
  struct TLSFit fita, fitb;
  struct Atom *atoms;

  atoms = hinge->atoms;  

  /* fit segments A and B */
  fita.atoms  = hinge->atoms;
  fita.istart = hinge->istarta;
  fita.iend   = hinge->ienda;
  isotropic_nonlinear_fit(&fita);

  fitb.atoms  = hinge->atoms;
  fitb.istart = hinge->istartb;
  fitb.iend   = hinge->iendb;
  isotropic_nonlinear_fit(&fitb);

  /* calculate hinge delta */
  num_atoms = 0;
  hdelta_a = 0.0;
  for (ia = hinge->istarta; ia <= hinge->ienda; ia++) {
    num_atoms++;

    calc_isotropic_uiso(fita.ITLS, atoms[ia].x - fita.ox, atoms[ia].y - fita.oy, atoms[ia].z - fita.oz, &u_iso_a);
    delta = atoms[ia].u_iso - u_iso_a;
    hdelta_a += delta*delta;
  }
  msd_a = hdelta_a / (num_atoms);

  num_atoms = 0;
  hdelta_b = 0.0;
  for (ia = hinge->istartb; ia <= hinge->iendb; ia++) {
    num_atoms++;

    calc_isotropic_uiso(fitb.ITLS, atoms[ia].x - fitb.ox, atoms[ia].y - fitb.oy, atoms[ia].z - fitb.oz, &u_iso_b);
    delta = atoms[ia].u_iso - u_iso_b;
    hdelta_b += delta*delta;
  }
  msd_b = hdelta_b / (num_atoms);

  /* calculate inter-segment deltas */
  num_atoms  = 0;
  hdelta_ab  = 0.0;
  hdelta_abo = 0.0;

  /* segment a atoms */
  for (ia = hinge->istarta; ia <= hinge->ienda; ia++) {
    num_atoms++;

    calc_isotropic_uiso(fita.ITLS, atoms[ia].x - fita.ox, atoms[ia].y - fita.oy, atoms[ia].z - fita.oz, &u_iso_a);
    calc_isotropic_uiso(fitb.ITLS, atoms[ia].x - fitb.ox, atoms[ia].y - fitb.oy, atoms[ia].z - fitb.oz, &u_iso_b);

    dela = atoms[ia].u_iso - u_iso_a;
    delb = atoms[ia].u_iso - u_iso_b;

    delta = u_iso_b - u_iso_a;
    hdelta_ab += delta*delta;

    delta = delb*delb;
    hdelta_abo += delta;
  }

  /* segment b atoms */
  for (ia = hinge->istartb; ia <= hinge->iendb; ia++) {
    num_atoms++;

    calc_isotropic_uiso(fita.ITLS, atoms[ia].x - fita.ox, atoms[ia].y - fita.oy, atoms[ia].z - fita.oz, &u_iso_a);
    calc_isotropic_uiso(fitb.ITLS, atoms[ia].x - fitb.ox, atoms[ia].y - fitb.oy, atoms[ia].z - fitb.oz, &u_iso_b);

    dela = atoms[ia].u_iso - u_iso_a;
    delb = atoms[ia].u_iso - u_iso_b;

    delta = u_iso_a - u_iso_b;
    hdelta_ab += delta*delta;

    delta = dela*dela;
    hdelta_abo += delta;
  }
  
  hinge->msd_a = msd_a;
  hinge->msd_b = msd_b;
  hinge->msd_c = 0.0;

  hinge->hdelta_ab = hdelta_ab / (num_atoms);
  hinge->hdelta_abo = hdelta_abo / (num_atoms);
  hinge->hdelta_c = 0.0;

  hinge->msd_c = 0.0;
  for (i = 1; i < 7; i++) {
    delta = fita.ITLS[i] - fitb.ITLS[i];
    hinge->msd_c += delta*delta;
  }

}

/*
 *  PYTHON INTERFACE
 */

typedef struct {
  PyObject_HEAD
  PyObject       *xmlrpc_chain;
  struct Atom *atoms;
} NLTLSModel_Object;

static void
NLTLSModel_dealloc(NLTLSModel_Object* self)
{
  if (self->atoms) {
    free(self->atoms);
    self->atoms = NULL;
  }

  Py_XDECREF(self->xmlrpc_chain);
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
NLTLSModel_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  NLTLSModel_Object *self;
  
  self = (NLTLSModel_Object *)type->tp_alloc(type, 0);
  if (self == NULL) {
    return NULL;
  }

  self->xmlrpc_chain = NULL;
  self->atoms = NULL;

  return (PyObject *)self;
}

static PyObject *
NLTLSModel_set_xmlrpc_chain(PyObject *py_self, PyObject *args)
{
  NLTLSModel_Object *self;
  PyObject *xmlrpc_chain;
  PyObject *atm_desc;
  PyObject *tmp;

  int i, j, num_atoms;
  char *strx;

  self = (NLTLSModel_Object *) py_self;

  if (!PyArg_ParseTuple(args, "O", &xmlrpc_chain)) {
    goto error;
  }

  /* de-reference the old xmlrpc chain and refrence the new one */
  tmp = self->xmlrpc_chain;
  Py_INCREF(xmlrpc_chain);
  self->xmlrpc_chain = xmlrpc_chain;
  Py_XDECREF(tmp);

  /* free the atoms array */
  if (self->atoms) {
    free(self->atoms);
    self->atoms = NULL;
  }

  /* allocate and fill the new atoms array */
  num_atoms = PyList_Size(self->xmlrpc_chain);
  if (num_atoms > 0) {
    self->atoms = (struct Atom *)malloc(sizeof(struct Atom) * num_atoms);

    for (i = 0; i < num_atoms; i++) {
      atm_desc = PyList_GetItem(self->xmlrpc_chain, i);

      /* set name */
      tmp = PyDict_GetItemString(atm_desc, "name");
      if (tmp == NULL) {
	goto error;
      }
      strx = PyString_AsString(tmp);
      if (strx == NULL) {
	goto error;
      }
      strncpy(self->atoms[i].name, strx, NAME_LEN);
     
      /* set frag_id */
      tmp = PyDict_GetItemString(atm_desc, "frag_id");
      if (tmp == NULL) {
	goto error;
      }
      strx = PyString_AsString(tmp);
      if (strx == NULL) {
	goto error;
      }
      strncpy(self->atoms[i].frag_id, strx, FRAG_ID_LEN);

      /* set x, y, z coordinates */
      tmp = PyDict_GetItemString(atm_desc, "x");
      if (tmp == NULL) {
	goto error;
      }
      if (!PyFloat_Check(tmp)) {
	goto error;
      }
      self->atoms[i].x = PyFloat_AsDouble(tmp);

      tmp = PyDict_GetItemString(atm_desc, "y");
      if (tmp == NULL) {
	goto error;
      }
      if (!PyFloat_Check(tmp)) {
	goto error;
      }
      self->atoms[i].y = PyFloat_AsDouble(tmp);

      tmp = PyDict_GetItemString(atm_desc, "z");
      if (tmp == NULL) {
	goto error;
      }
      if (!PyFloat_Check(tmp)) {
	goto error;
      }
      self->atoms[i].z = PyFloat_AsDouble(tmp);

      tmp = PyDict_GetItemString(atm_desc, "u_iso");
      if (tmp == NULL) {
	goto error;
      }
      if (!PyFloat_Check(tmp)) {
	goto error;
      }
      self->atoms[i].u_iso = PyFloat_AsDouble(tmp);

      tmp = PyDict_GetItemString(atm_desc, "weight");
      if (tmp == NULL) {
	goto error;
      }
      if (!PyFloat_Check(tmp)) {
	goto error;
      }
      self->atoms[i].weight = PyFloat_AsDouble(tmp);
      self->atoms[i].sqrt_weight = sqrt(self->atoms[i].weight);

      /* get U tensor parameters */
      for (j = 0; j < U_NUM_PARAMS; j++) {
	tmp = PyDict_GetItemString(atm_desc, U_PARAM_NAMES[j]);
	if (tmp == NULL) {
	  goto error;
	}
	if (!PyFloat_Check(tmp)) {
	  goto error;
	}
	self->atoms[i].U[j] = PyFloat_AsDouble(tmp);
      }
    }
  }

  Py_INCREF(Py_None);
  return Py_None;

 error:
  return NULL;
}

static PyObject *
NLTLSModel_clear_xmlrpc_chain(PyObject *py_self, PyObject *args)
{
  NLTLSModel_Object *self;
  PyObject *tmp;

  self = (NLTLSModel_Object *) py_self;

  tmp = self->xmlrpc_chain;
  self->xmlrpc_chain = NULL;
  Py_XDECREF(tmp);

  /* free the atoms array */
  if (self->atoms) {
    free(self->atoms);
    self->atoms = NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
NLTLSModel_anisotropic_fit_segment(PyObject *py_self, PyObject *args)
{
  NLTLSModel_Object *self;
  struct TLSFit fit;
  PyObject *py_floatx, *py_intx, *rdict;
  int i, info;

  self = (NLTLSModel_Object *) py_self;

  if (!PyArg_ParseTuple(args, "ii", &fit.istart, &fit.iend)) {
    goto error;
  }
  if (self->atoms == NULL) {
    goto error;
  }

  fit.atoms = self->atoms;
  info = anisotropic_nonlinear_fit(&fit);

  /* construct return dictioary with results */
  rdict = PyDict_New();

  /* set minimization exit status */
  py_intx = PyInt_FromLong(info);
  PyDict_SetItemString(rdict, "info", py_intx);
  Py_DECREF(py_intx);
  
  /* TLS parmeters */
  for (i = 0; i < NL_ATLS_NUM_PARAMS; i++) {
    py_floatx = PyFloat_FromDouble(fit.NL_ATLS[i]);
    PyDict_SetItemString(rdict, NL_ATLS_PARAM_NAMES[i], py_floatx);
    Py_DECREF(py_floatx);
  }

  for (i = 0; i < ATLS_NUM_PARAMS; i++) {
    py_floatx = PyFloat_FromDouble(fit.ATLS[i]);
    PyDict_SetItemString(rdict, ATLS_PARAM_NAMES[i], py_floatx);
    Py_DECREF(py_floatx);
  }
   
  /* TLS origin */
  py_floatx = PyFloat_FromDouble(fit.ox);
  PyDict_SetItemString(rdict, "x", py_floatx);
  Py_DECREF(py_floatx);
  py_floatx = PyFloat_FromDouble(fit.oy);
  PyDict_SetItemString(rdict, "y", py_floatx);
  Py_DECREF(py_floatx);
  py_floatx = PyFloat_FromDouble(fit.oz);
  PyDict_SetItemString(rdict, "z", py_floatx);
  Py_DECREF(py_floatx);

  /* TLS least squares residual */
  py_floatx = PyFloat_FromDouble(fit.alsqr);
  PyDict_SetItemString(rdict, "alsqr", py_floatx);
  Py_DECREF(py_floatx);

  return rdict;

 error:
  return NULL;
}

static PyObject *
NLTLSModel_isotropic_fit_segment(PyObject *py_self, PyObject *args)
{
  NLTLSModel_Object *self;
  struct TLSFit fit;
  PyObject *py_floatx, *py_intx, *rdict;
  int i, info;

  self = (NLTLSModel_Object *) py_self;

  if (!PyArg_ParseTuple(args, "ii", &fit.istart, &fit.iend)) {
    goto error;
  }
  if (self->atoms == NULL) {
    goto error;
  }

  fit.atoms = self->atoms;
  info = isotropic_nonlinear_fit(&fit);

  /* construct return dictioary with results */
  rdict = PyDict_New();

  /* set minimization exit status */
  py_intx = PyInt_FromLong(info);
  PyDict_SetItemString(rdict, "info", py_intx);
  Py_DECREF(py_intx);
  
  /* TLS parameters */
  for (i = 0; i < ITLS_NUM_PARAMS; i++) {
    py_floatx = PyFloat_FromDouble(fit.ITLS[i]);
    PyDict_SetItemString(rdict, ITLS_PARAM_NAMES[i], py_floatx);
    Py_DECREF(py_floatx);
  }

  /* non-linear isotropic TLS parameters */
  for (i = 0; i < NL_ITLS_NUM_PARAMS; i++) {
    py_floatx = PyFloat_FromDouble(fit.NL_ITLS[i]);
    PyDict_SetItemString(rdict, NL_ITLS_PARAM_NAMES[i], py_floatx);
    Py_DECREF(py_floatx);
  }

  /* TLS origin */
  py_floatx = PyFloat_FromDouble(fit.ox);
  PyDict_SetItemString(rdict, "x", py_floatx);
  Py_DECREF(py_floatx);
  py_floatx = PyFloat_FromDouble(fit.oy);
  PyDict_SetItemString(rdict, "y", py_floatx);
  Py_DECREF(py_floatx);
  py_floatx = PyFloat_FromDouble(fit.oz);
  PyDict_SetItemString(rdict, "z", py_floatx);
  Py_DECREF(py_floatx);

  /* TLS least squares residual */
  py_floatx = PyFloat_FromDouble(fit.ilsqr);
  PyDict_SetItemString(rdict, "ilsqr", py_floatx);
  Py_DECREF(py_floatx);

  return rdict;

 error:
  return NULL;
}

static PyObject *
NLTLSModel_calc_isotropic_hinge_delta(PyObject *py_self, PyObject *args)
{
  NLTLSModel_Object *self;
  PyObject *py_floatx, *rdict;
  struct IHingeContext hinge;

  self = (NLTLSModel_Object *) py_self;

  /* fill in fields in the TLSFitContext structure */
  if (self->atoms==NULL) {
    goto error;
  }

  hinge.atoms = self->atoms;
  if (!PyArg_ParseTuple(args, "iiii", &hinge.istarta, &hinge.ienda, &hinge.istartb, &hinge.iendb)) {
    goto error;
  }

  calc_isotropic_hinge_delta(&hinge);

  /* construct return dictioary with results */
  rdict = PyDict_New();

  /* return msds of segments */
  py_floatx = PyFloat_FromDouble(hinge.msd_a);
  PyDict_SetItemString(rdict, "msd_a", py_floatx);
  Py_DECREF(py_floatx);

  py_floatx = PyFloat_FromDouble(hinge.msd_b);
  PyDict_SetItemString(rdict, "msd_b", py_floatx);
  Py_DECREF(py_floatx);

  py_floatx = PyFloat_FromDouble(hinge.msd_c);
  PyDict_SetItemString(rdict, "msd_c", py_floatx);
  Py_DECREF(py_floatx);

  /* return various hinge delta values */
  py_floatx = PyFloat_FromDouble(hinge.hdelta_ab);
  PyDict_SetItemString(rdict, "hdelta_ab", py_floatx);
  Py_DECREF(py_floatx);

  py_floatx = PyFloat_FromDouble(hinge.hdelta_abo);
  PyDict_SetItemString(rdict, "hdelta_abo", py_floatx);
  Py_DECREF(py_floatx);

  py_floatx = PyFloat_FromDouble(hinge.hdelta_c);
  PyDict_SetItemString(rdict, "hdelta_c", py_floatx);
  Py_DECREF(py_floatx);
  return rdict;

 error:
  return NULL;
}

static PyMethodDef NLTLSModel_methods[] = {
    {"set_xmlrpc_chain", 
     (PyCFunction) NLTLSModel_set_xmlrpc_chain, 
     METH_VARARGS,
     "Sets the Python list containing one dictionary for each atom." },

    {"clear_xmlrpc_chain", 
     (PyCFunction) NLTLSModel_clear_xmlrpc_chain, 
     METH_VARARGS,
     "Clears the Python list of atom descriptions." },

    {"anisotropic_fit_segment", 
     (PyCFunction) NLTLSModel_anisotropic_fit_segment, 
     METH_VARARGS,
     "Performs nonlinear fit of the anisotropic TLS model to the given atoms." },

    {"isotropic_fit_segment", 
     (PyCFunction) NLTLSModel_isotropic_fit_segment, 
     METH_VARARGS,
     "Performs nonlinear fit of the isotropic TLS model to the given atoms." },

    {"calc_isotropic_hinge_delta",
     (PyCFunction) NLTLSModel_calc_isotropic_hinge_delta, 
     METH_VARARGS,
     "Calculates the value of the hinge-delta function." },

    {NULL}  /* Sentinel */
};

static PyTypeObject NLTLSModel_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                          /*ob_size*/
    "NLTLSModel",                               /*tp_name*/
    sizeof(NLTLSModel_Object),                  /*tp_basicsize*/
    0,                                          /*tp_itemsize*/
    (destructor)NLTLSModel_dealloc,             /*tp_dealloc*/
    0,                                          /*tp_print*/
    0,                                          /*tp_getattr*/
    0,                                          /*tp_setattr*/
    0,                                          /*tp_compare*/
    0,                                          /*tp_repr*/
    0,                                          /*tp_as_number*/
    0,                                          /*tp_as_sequence*/
    0,                                          /*tp_as_mapping*/
    0,                                          /*tp_hash */
    0,                                          /*tp_call*/
    0,                                          /*tp_str*/
    0,                                          /*tp_getattro*/
    0,                                          /*tp_setattro*/
    0,                                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /*tp_flags*/
    "NLTLSModel objects",                       /* tp_doc */
    0,		                                /* tp_traverse */
    0,		                                /* tp_clear */
    0,		                                /* tp_richcompare */
    0,		                                /* tp_weaklistoffset */
    0,		                                /* tp_iter */
    0,		                                /* tp_iternext */
    NLTLSModel_methods,                         /* tp_methods */
    0,                                          /* tp_members */
    0,                                          /* tp_getset */
    0,                                          /* tp_base */
    0,                                          /* tp_dict */
    0,                                          /* tp_descr_get */
    0,                                          /* tp_descr_set */
    0,                                          /* tp_dictoffset */
    0,                                          /* tp_init */
    0,                                          /* tp_alloc */
    NLTLSModel_new,                             /* tp_new */
};


static PyMethodDef NONLINEARTLS_METHODS[] = {
  {NULL, NULL, 0, NULL}
};

DL_EXPORT(void)
initnonlineartls(void)
{
  PyObject *m;
  
  
  if (PyType_Ready(&NLTLSModel_Type) < 0)
    return;

  m = Py_InitModule("nonlineartls", NONLINEARTLS_METHODS);
  
  NONLINEARTLS_ERROR = PyErr_NewException("nonlineartls.error", NULL, NULL);
  Py_INCREF(NONLINEARTLS_ERROR);
  PyModule_AddObject(m, "error", NONLINEARTLS_ERROR);


  /* add the NLTLSModel class */
  Py_INCREF(&NLTLSModel_Type);
  PyModule_AddObject(m, "NLTLSModel", (PyObject *)&NLTLSModel_Type);
}
