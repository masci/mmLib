/* lineartls.c 
 * Jay Painter <jpaint@u.washington.edu>
 * Sept 7, 2005
 *
 * Copyright 2002 by TLSMD Development Group (see AUTHORS file)
 * This code is part of the TLSMD distribution and governed by
 * its license.  Please see the LICENSE file that should have been
 * included as part of this package.
 *
 * Implementation of a Python module for the unconstrained linear 
 * fitting of TLS parameters given a list of atoms with 
 * crystallographically refined ADPs.  Uses LAPACK.  Solved linear
 * equations by SVD.
 */
#include "Python.h"
#include "structmember.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <pthread.h>

/* LAPACK */
extern void
dgesdd_(char *, int*, int*, double*, int*, double*, double*, int *, double*, int*, double*, int*, int*, int*);

/* macros */
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))

#define PI    3.1415926535897931
#define SQRT2 1.4142135623730951
#define U2B   ((8.0*PI*PI))
#define B2U   ((1.0/U2B))


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

/* Isotropic TLS Model 
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

/* Anisotropic TLS model parameter indexes and labels */
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

/* Atom structure 
 * structure used to store the information of one atom 
 */
#define NAME_LEN     8
#define FRAG_ID_LEN  8
struct Atom {
  char    name[NAME_LEN];             /* atom name */
  char    frag_id[FRAG_ID_LEN];       /* fragment id (residue name) */
  int     ifrag;                      /* fragment index */
  double  x, y, z;
  double  xtls, ytls, ztls;
  double  u_iso;
  double  u_iso_tmp;
  double  U[6];
  double  Utmp[6];
  double  sqrt_weight;
  double  sqrt_weight_tmp;
};


/* structure used to store information on the chain the algorithm
 * is currently fitting; this structure contains the dynamically 
 * allocated memory used by the LAPACK SVD routine DGESDD and therefore
 * must be allocated and freed with new_chain/delete_chain
 */
struct Chain {
  struct Atom *atoms;
  int          num_atoms;
  
  double *A;
  double *b;
  double *S;
  double *U;
  double *VT;
  double *WORK;
  int     LWORK;
  int    *IWORK;
};

/* context structure for fitting one TLS segment using the  isotropic TLS model */
struct TLSFitContext {
  struct Chain      *chain;                 /* pointer to Chain structure */

  int                istart;                /* index of first atom in chain->atoms */
  int                iend;                  /* index of last atom in chain->atoms */

  double             ox;                    /* origin of TLS tensors */
  double             oy;
  double             oz;

  double             ITLS[ITLS_NUM_PARAMS]; /* isotropic TLS model params */
  double             ATLS[ATLS_NUM_PARAMS]; /* ansiotropic TLS model params */

  double             ilsqr;                 /* least-squares residual of isotropic TLS model */
  double             alsqr;                 /* least-squares residual of anisotropic TLS model */
};

struct IHingeContext {
  struct Chain      *chain;                 /* pointer to Chain structure */

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


/* incomplete gamma functions */

/* Returns the value ln[gamma(XX)] for XX>0.  Full accuracy is obtained for XX>1.
 * For 0<XX<1, the reflection formula can be used first.
 */
static void
gammaln(double xx, double *retval)
{
  int j;
  double x, ser, tmp;
  
  double cof[6] = {76.18009173, -86.50532033, 24.01409822, -1.231739516, 0.120858003E-2, -0.536382E-5};
  double stp  = 2.50662827465;
  double half = 0.5;
  double one  = 1.0;
  double fpf  = 5.0;

  x = xx - one;
  tmp = x + fpf;
  tmp = (x + half) * log(tmp) - tmp;
  ser = one;

  for (j = 0; j < 6; j++) {
    x += one;
    ser += cof[j] / x;
  }

  *retval = tmp + log(stp*ser);
}

static void
gser(double a, double x, double *gamser, double *gln)
{
  int i;
  double ap, sum, del;

  int itmax = 100;
  double eps = 3.0E-7;
 
  if (x <= 0.0) {
    if (x < 0.0) return;
    *gamser = 0.0;
    return;
  }

  ap = a;
  sum = 1.0 / a;
  del = sum;

  for (i = 0; i < itmax; i++) {
    ap += 1.0;
    del = del * x/ap;
    sum += del;

    if (fabs(del) < fabs(sum)*eps) break;
  }

  gammaln(a, gln);
  *gamser = sum * exp(-x + a*log(x) - *gln);
}

static void
gcf(double a, double x, double *gamser, double *gln)
{
  int n;
  double gold, a0, a1, b0, b1, fac, an, ana, g, anf;
  
  int itmax = 100;
  double eps = 3.0E-7;

  g    = 0.0;
  gold = 0.0;
  a0   = 1.0;
  a1   = x;
  b0   = 0.0;
  b1   = 1.0;
  fac  = 1.0;

  for (n = 1; n <= itmax; n++) {
    an = n;
    ana = an - a;
    a0 = (a1 + a0*ana) * fac;
    b0 = (b1 + b0*ana) * fac;
    anf = an * fac;
    a1 = x*a0 + anf*a1;
    b1 = x*b0 + anf*b1;

    if (a1 != 0.0) {
      fac = 1.0 / a1;
      g = b1*fac;
      if (fabs((g-gold)/g) < eps) break;
      gold = g;
    }
  }

  gammaln(a, gln);
  *gamser = exp(-x + a*log(x) - *gln) * g;
}

static void 
gammaq(double a, double x, double *retval)
{
  double gamser, gammcf, gln;
  
  if (x < 0.0 || a < 0.0) return;

  if (x < a+1.0) {
    gser(a, x, &gamser, &gln);
    *retval = 1.0 - gamser;
  } else {
    gcf(a, x, &gammcf, &gln);
    *retval = gammcf;
  }
}

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
  if (strcmp(atoms[ia].name,"N")==0)  return 1;
  if (strcmp(atoms[ia].name,"CA")==0) return 1;
  if (strcmp(atoms[ia].name,"C")==0)  return 1;
  if (strcmp(atoms[ia].name,"O")==0)  return 1;
  if (strcmp(atoms[ia].name,"CB")==0) return 1;

  return 0;
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

/* Sets the one row of matrix A starting at A[i,j] with the istropic
 * TLS model coefficents for a atom located at t position x, y, z with
 * least-squares weight w.  Matrix A is filled to coumn j+12.
 *
 * The matrix A(m,n) is filled in FORTRAN-style, that is, assuming
 * a column-major memory layout.  That's because, athough the
 * matrix is contructed in C, the SVD subroutine from LAPACK is 
 * written in FORTRAN.
 *
 */
inline void
set_ITLS_Ab(double *A, double *b, int m, int n, int row, double uiso, double x, double y, double z, double w)
{
#define FA(__i, __j) A[__i + (m * __j)]

  double xx, yy, zz, xy, xz, yz;
   
  xx = x*x;
  yy = y*y;
  zz = z*z;
  xy = x*y;
  xz = x*z;
  yz = y*z;

  /* set b */
  b[row] = w * uiso;

  /* T iso */
  FA(row, ITLS_T) = w * 1.0;

  /* l11, l22, l33, l12, l13, l23 */
  FA(row, ITLS_L11) = w * ((zz + yy) / 3.0);
  FA(row, ITLS_L22) = w * ((xx + zz) / 3.0);
  FA(row, ITLS_L33) = w * ((xx + yy) / 3.0);

  FA(row, ITLS_L12) = w * ((-2.0 * xy) / 3.0);
  FA(row, ITLS_L13) = w * ((-2.0 * xz) / 3.0);
  FA(row, ITLS_L23) = w * ((-2.0 * yz) / 3.0);

  FA(row, ITLS_S1)  = w * (( 2.0 * z) / 3.0);
  FA(row, ITLS_S2)  = w * (( 2.0 * y) / 3.0);
  FA(row, ITLS_S3)  = w * (( 2.0 * x) / 3.0);
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
set_ATLS_Ab(double *A, double *b, int m, int n, int row, double U[6], double x, double y, double z, double w)
{
#define FA(__i,__j) A[__i + (m * __j)]

  int rowU11, rowU22, rowU33, rowU12, rowU13, rowU23;
  double xx, yy, zz, xy, xz, yz;

  xx = x*x;
  yy = y*y;
  zz = z*z;
  xy = x*y;
  xz = x*z;
  yz = y*z;

  /* calculate row indexes */
  rowU11 = row;
  rowU22 = row + 1;
  rowU33 = row + 2;
  rowU12 = row + 3;
  rowU13 = row + 4;
  rowU23 = row + 5;

  /* set b  */
  b[rowU11] = w * U[0];
  b[rowU22] = w * U[1];
  b[rowU33] = w * U[2];
  b[rowU12] = w * U[3];
  b[rowU13] = w * U[4];
  b[rowU23] = w * U[5];

  /* set A */
  FA(rowU11, ATLS_T11) = w * 1.0;
  FA(rowU11, ATLS_L22) = w *        zz;
  FA(rowU11, ATLS_L33) = w *        yy;
  FA(rowU11, ATLS_L23) = w * -2.0 * yz;
  FA(rowU11, ATLS_S31) = w * -2.0 *  y;
  FA(rowU11, ATLS_S21) = w *  2.0 *  z;

  FA(rowU22, ATLS_T22) = w * 1.0;
  FA(rowU22, ATLS_L11) = w *        zz;
  FA(rowU22, ATLS_L33) = w *        xx;
  FA(rowU22, ATLS_L13) = w * -2.0 * xz;
  FA(rowU22, ATLS_S12) = w * -2.0 *  z;
  FA(rowU22, ATLS_S32) = w *  2.0 *  x;

  FA(rowU33, ATLS_T33) = w * 1.0;
  FA(rowU33, ATLS_L11) = w *        yy;
  FA(rowU33, ATLS_L22) = w *        xx;
  FA(rowU33, ATLS_L12) = w * -2.0 * xy;
  FA(rowU33, ATLS_S23) = w * -2.0 *  x;
  FA(rowU33, ATLS_S13) = w *  2.0 *  y;

  FA(rowU12, ATLS_T12)   = w * 1.0;
  FA(rowU12, ATLS_L33)   = w * -xy;
  FA(rowU12, ATLS_L23)   = w *  xz;
  FA(rowU12, ATLS_L13)   = w *  yz;
  FA(rowU12, ATLS_L12)   = w * -zz;
  FA(rowU12, ATLS_S2211) = w *   z;
  FA(rowU12, ATLS_S31)   = w *   x;
  FA(rowU12, ATLS_S32)   = w *  -y;
    
  FA(rowU13, ATLS_T13)   = w * 1.0;
  FA(rowU13, ATLS_L22)   = w * -xz;
  FA(rowU13, ATLS_L23)   = w *  xy;
  FA(rowU13, ATLS_L13)   = w * -yy;
  FA(rowU13, ATLS_L12)   = w *  yz;
  FA(rowU13, ATLS_S1133) = w *   y;
  FA(rowU13, ATLS_S23)   = w *   z;
  FA(rowU13, ATLS_S21)   = w *  -x;
    
  FA(rowU23, ATLS_T23)   = w * 1.0;
  FA(rowU23, ATLS_L11)   = w * -yz;
  FA(rowU23, ATLS_L23)   = w * -xx;
  FA(rowU23, ATLS_L13)   = w *  xy;
  FA(rowU23, ATLS_L12)   = w *  xz;
  FA(rowU23, ATLS_S2211) = w *  -x;
  FA(rowU23, ATLS_S1133) = w *  -x;
  FA(rowU23, ATLS_S12)   = w *   y;
  FA(rowU23, ATLS_S13)   = w *  -z;

#undef FA
}

/* solve for x given the singular value decomposition of a matrix
 * into U, S, and Vt 
 */
static void 
solve_SVD(int m, int n, double *x, double *b, double *U, double *S, double *VT, double *WORK)
{
#define FU(__i, __j)   U[__i + (m * __j)]
#define FVT(__i, __j)  VT[__i + (n * __j)]

  int i, j;
  double smax, scutoff, dtmp;

  /* invert the diagonal matrix S in place, and filter out any
   * unusually small singular values
   */

  for (smax = S[0], i = 1; i < n; i++) {
    smax = MAX(smax, S[i]);
  }
  scutoff = smax * 1E-10;
  for (i = 0; i < n; i++) {
    if (S[i] > scutoff) {
      S[i] = 1.0 / S[i];
    } else {
      S[i] = 0.0;
    }
  }

  /* matrix multiply Ut(n,m)*b(m) and store the result in x
   */
  for (i = 0; i < n; i++) {
    dtmp = 0.0;
    for (j = 0; j < m; j++) {
      dtmp += FU(j,i) * b[j];
    }
    WORK[i] = dtmp;
  }

  /* matrix multiply inverse-S by Ut*b */
  for (i = 0; i < n; i++) {
    WORK[i] *= S[i];
  }

  /* matrix multiple V*x */
  for (i = 0; i < n; i++) {
    dtmp = 0.0;
    for (j = 0; j < n; j++) {
      dtmp += FVT(j,i) * WORK[j];
    }
    x[i] = dtmp;
  }

#undef FU
#undef FVT
}


/* frees all memory from a allocated struct Chain 
 */
static void
delete_chain(struct Chain *chain)
{
  if (chain->atoms != NULL)     free(chain->atoms);
  if (chain->A != NULL)         free(chain->A);
  if (chain->b != NULL)         free(chain->b);
  if (chain->S != NULL)         free(chain->S);
  if (chain->U != NULL)         free(chain->U);
  if (chain->VT != NULL)        free(chain->VT);
  if (chain->WORK != NULL)      free(chain->WORK);
  if (chain->IWORK != NULL)     free(chain->IWORK);

  free(chain);
}

/* allocate a new struct Chain with working memory buffers large 
 * enough to accomodate the SVD solution of the TLS equations for
 * up to num_atoms using DGESDD from LAPACK
 */
static struct Chain *
new_chain(int num_atoms)
{
  char jobz;
  int num_rows, num_cols, info;
  double tmp_WORK;
  struct Chain *chain;

  /* six parameters per atom */
  num_rows = 6 * num_atoms;
  num_cols = ATLS_NUM_PARAMS;

  /* allocate chain struct and initialize */
  chain = malloc(sizeof(struct Chain));
  if (chain == NULL) {
    printf("new_chain: struct Chain\n");
    goto error;
  }

  chain->num_atoms = num_atoms;
  chain->A     = NULL;
  chain->b     = NULL;
  chain->S     = NULL;
  chain->U     = NULL;
  chain->VT    = NULL;
  chain->WORK  = NULL;
  chain->IWORK = NULL;

  /* allocate memory blocks */
  chain->atoms = malloc(sizeof(struct Atom) * num_atoms);
  if (chain->atoms == NULL) {
    printf("new_chain: chain->atoms\n");
    goto error;
  }

  chain->A = malloc(sizeof(double) * (num_rows * num_cols));
  if (chain->A == NULL) {
    printf("new_chain: chain->A\n");
    goto error;
  }

  chain->b = malloc(sizeof(double) * num_rows);
  if (chain->b == NULL) {
    printf("new_chain: chain->b\n");
    goto error;
  }

  chain->S = malloc(sizeof(double) * num_cols);
  if (chain->S == NULL) {
    printf("new_chain: chain->S\n");
    goto error;
  }

  /*  U(num_cols, num_rows) */
  chain->U = malloc(sizeof(double) * (num_cols * num_rows));
  if (chain->U == NULL) {
    printf("new_chain: chain->U\n");
    goto error;
  }
  
  /* V(num_cols, num_cols) */
  chain->VT = malloc(sizeof(double) * (num_cols * num_cols));
  if (chain->VT == NULL) {
    printf("new_chain: chain->VT\n");
    goto error;
  }

  /* iwork = Numeric.zeros((8*min(m,n),), 'i') */
  chain->IWORK = malloc(sizeof(int) * (8 * num_cols));
  if (chain->IWORK == NULL) {
    printf("new_chain: chain->IWORK\n");
    goto error;
  }

  /* calculate the ideal size of the WORK memory block */
  jobz = 'S';
  chain->LWORK = -1;
  info = 0;
  dgesdd_(&jobz,  
	  &num_rows,
	  &num_cols, 
	  chain->A, 
	  &num_rows,
	  chain->S, 
	  chain->U,
	  &num_rows, 
	  chain->VT,
	  &num_cols,
	  &tmp_WORK, 
	  &chain->LWORK,
	  chain->IWORK,
	  &info);
    
  chain->LWORK = tmp_WORK;

  chain->WORK = malloc(sizeof(double) * chain->LWORK);
  if (chain->WORK == NULL) {
    printf("new_chain: chain->WORK\n");
    goto error;
  }

  return chain;

  /* if malloc() fails */
 error:
  if (chain != NULL) delete_chain(chain);
  return NULL;
}

inline int
calc_istart_iend(struct Atom *atoms, char *frag_id1, char *frag_id2, int *istart, int *iend)
{
  return 1;
}

/* calculates the centroid of the atoms indexed between istart and iend
 * then returns the centroid coordinates in x, y, z
 */
static void
calc_centroid(struct Atom *atoms, int istart, int iend, double *x, double *y, double *z)
{
  int i, n;
  double cx, cy, cz;

  n  = 0;
  cx = 0.0;
  cy = 0.0;
  cz = 0.0;

  for (i = istart; i <= iend; i++) {
    n  += 1;
    cx += atoms[i].x;
    cy += atoms[i].y;
    cz += atoms[i].z;
  }

  if (n > 0) {
    *x = cx / n;
    *y = cy / n;
    *z = cz / n;
  } else {
    *x = 0.0;
    *y = 0.0;
    *z = 0.0;
  }
}

/**********************************************************************************************************************************************************/
/* ISOTROPIC MODELS */

static void
linear_isotropic_fit_parameters(struct TLSFitContext *fit)
{
  char jobz;
  int num_atoms, num_rows, num_cols, row, ia, istart, iend, info;
  double ox, oy, oz;
  double *A, *b;
  struct Atom *atoms;

  /* optimization */
  atoms = fit->chain->atoms;

  A  = fit->chain->A;
  b  = fit->chain->b;

  /* calculate the number of atoms to be fit */
  num_atoms = fit->iend - fit->istart + 1;
  istart    = fit->istart;
  iend      = fit->iend;

  /* calculate the centroid of the atoms which are to
   * be fit and use it as the origin of the TLS tensors
   */
  calc_centroid(fit->chain->atoms, fit->istart, fit->iend, &fit->ox, &fit->oy, &fit->oz);
  ox = fit->ox;
  oy = fit->oy;
  oz = fit->oz;

  num_rows = num_atoms;
  num_cols = ITLS_NUM_PARAMS;

  zero_dmatrix(A, num_rows, num_cols);

  for (ia = istart, row = 0; ia <= iend; ia++, row++) {
    set_ITLS_Ab(A, b, num_rows, num_cols, row, atoms[ia].u_iso, atoms[ia].x-ox, atoms[ia].y-oy, atoms[ia].z-oz, atoms[ia].sqrt_weight);
  }

  /* solve for isotropic TLS model */
  jobz = 'S';
  dgesdd_(&jobz,  
	  &num_rows,
	  &num_cols, 
	  A,
	  &num_rows,
	  fit->chain->S, 
	  fit->chain->U,
	  &num_rows,
	  fit->chain->VT,
	  &num_cols,
	  fit->chain->WORK, 
	  &fit->chain->LWORK,
	  fit->chain->IWORK,
	  &info);

  if (info != 0) {
    printf("DGESDD ERROR(isotropic): info = %d\n", info);
  }

  solve_SVD(num_rows, num_cols, fit->ITLS, b, fit->chain->U, fit->chain->S, fit->chain->VT, fit->chain->WORK);
}


static void
linear_isotropic_fit_parameters_filter_outliers(struct TLSFitContext *fit)
{
  int i, ia, ia_max_delta, orig_num_atoms, num_atoms, min_atoms;
  double delta, max_delta, msd, rmsd, u_iso;
  struct Atom *atoms;

  atoms = fit->chain->atoms;

  /* set all atoms to unit weights */
  orig_num_atoms = 0;
  for (ia = fit->istart; ia <= fit->iend; ia++) {
    orig_num_atoms++;
    atoms[ia].sqrt_weight = 1.0;
  }

  min_atoms = 0.90 * orig_num_atoms;
  min_atoms = MAX(10, min_atoms);

  for (i = 1; ; i++) {
    linear_isotropic_fit_parameters(fit);

    /* calculate MSD */
    num_atoms = 0;
    msd = 0.0;

    for (ia = fit->istart; ia <= fit->iend; ia++) {
      if (atoms[ia].sqrt_weight < 1.0) continue;

      num_atoms++;
      calc_isotropic_uiso(fit->ITLS, atoms[ia].x - fit->ox, atoms[ia].y - fit->oy, atoms[ia].z - fit->oz, &u_iso);
      delta = atoms[ia].u_iso - u_iso;
      msd += delta*delta;
    }

    /* can't have too few atoms! */
    if (num_atoms <= min_atoms) break;

    msd = msd / num_atoms;
    rmsd = sqrt(msd);

    /* reject outliers */
    ia_max_delta = 0;
    max_delta = 0.0;

    for (ia = fit->istart; ia <= fit->iend; ia++) {
      if (atoms[ia].sqrt_weight < 1.0) continue;

      calc_isotropic_uiso(fit->ITLS, atoms[ia].x - fit->ox, atoms[ia].y - fit->oy, atoms[ia].z - fit->oz, &u_iso);
      delta = fabs(atoms[ia].u_iso - u_iso);

      if (delta > max_delta) {
	max_delta = delta;
	ia_max_delta = ia;
      }
    }

    printf("iteration %d: rejecting atom %s:%s with delta=%f at rmsd=%f\n", i, atoms[ia_max_delta].frag_id, atoms[ia_max_delta].name, U2B*max_delta, U2B*rmsd);
    atoms[ia_max_delta].sqrt_weight = 0.0;
  }

  /* don't restore atom weights */
}

static void
linear_isotropic_fit_parameters_infit(struct TLSFitContext *fit)
{
  int i, ia, orig_num_atoms, num_atoms, num_atoms_pdelta, num_atoms_ndelta, max_ndelta_atoms;
  double delta, msd, rmsd, u_iso;
  struct Atom *atoms;

  atoms = fit->chain->atoms;

  /* save original atom weights */
  orig_num_atoms = 0;
  for (ia = fit->istart; ia <= fit->iend; ia++) {
    orig_num_atoms++;
    atoms[ia].sqrt_weight_tmp = atoms[ia].sqrt_weight;
    atoms[ia].sqrt_weight = 1.0;
  }

  /* allow only 1% of the atoms to be over-estimated by the TLS model */
  max_ndelta_atoms = orig_num_atoms * 0.10;
  max_ndelta_atoms = MAX(max_ndelta_atoms, 1);

  for (i = 1; i < 1000; i++) {
    linear_isotropic_fit_parameters(fit);

    /* calculate MSD */
    num_atoms = 0;
    num_atoms_pdelta = 0;
    num_atoms_ndelta = 0;
    msd = 0.0;

    for (ia = fit->istart; ia <= fit->iend; ia++) {
      num_atoms++;

      calc_isotropic_uiso(fit->ITLS, atoms[ia].x - fit->ox, atoms[ia].y - fit->oy, atoms[ia].z - fit->oz, &u_iso);
      delta = atoms[ia].u_iso - u_iso;
      msd += delta*delta;

      if (delta >= 0.0) {
	num_atoms_pdelta++;
      } else {
	num_atoms_ndelta++;
	atoms[ia].sqrt_weight += 1.0;
      }
    }

    msd = msd / num_atoms;
    rmsd = sqrt(msd);

    printf("iteration %d: rejecting +/- estimate %d/%d with at rmsd=%f\n", i, num_atoms_pdelta, num_atoms_ndelta, U2B*rmsd);

    if (num_atoms_ndelta <= max_ndelta_atoms) break;
  }

  /* reset atom weights */
  for (ia = fit->istart; ia <= fit->iend; ia++) {
    atoms[ia].sqrt_weight = atoms[ia].sqrt_weight_tmp;
  }
}

static void
linear_isotropic_fit_segment(struct TLSFitContext *fit)
{
  int ia, istart, iend, num_atoms, num_residues, ia_res_start;
  double ox, oy, oz;
  double u_tls, delta, chi2, sum_weight, tmp;
  struct Atom *atoms;

  /* optimization */
  atoms  = fit->chain->atoms;
  istart = fit->istart;
  iend   = fit->iend;

  /* initialize atoms in the segment */
  num_atoms = 0;
  sum_weight = 0.0;
  num_residues = 1;
  ia_res_start = istart;

  for (ia = istart; ia <= iend; ia++) {
    /* count the number of new residues */
    if (strcmp(atoms[ia_res_start].frag_id, atoms[ia].frag_id)!=0) {
      num_residues++;
      ia_res_start = ia;
    }
    num_atoms++;
    sum_weight += atoms[ia].sqrt_weight;
  }

  /* fit TLS parameters */
  linear_isotropic_fit_parameters(fit);

  ox = fit->ox;
  oy = fit->oy;
  oz = fit->oz;

  /* calculate residual */
  chi2 = 0.0;

  for (ia = istart; ia <= iend; ia++) {
    calc_isotropic_uiso(fit->ITLS, atoms[ia].x-ox, atoms[ia].y-oy, atoms[ia].z-oz, &u_tls);

    delta = u_tls - atoms[ia].u_iso;
    tmp = atoms[ia].sqrt_weight * delta;
    chi2 += tmp * tmp;
  }

  fit->ilsqr = num_residues * (chi2 / sum_weight);
  fit->ilsqr = num_residues * (chi2 / (num_atoms - ITLS_NUM_PARAMS));
}


/**********************************************************************************************************************************************************/
/* ANISOTROPIC MODELS */

static void
linear_anisotropic_fit_parameters(struct TLSFitContext *fit)
{
  char jobz;
  int num_atoms, num_rows, num_cols, row, ia, istart, iend, info;
  double ox, oy, oz;
  double *A, *b;
  struct Atom *atoms;

  /* optimization */
  atoms = fit->chain->atoms;
  A = fit->chain->A;
  b = fit->chain->b;

  /* calculate the number of atoms to be fit */
  num_atoms = fit->iend - fit->istart + 1;
  istart    = fit->istart;
  iend      = fit->iend;

  /* calculate the centroid of the atoms which are to
   * be fit and use it as the origin of the TLS tensors
   */
  calc_centroid(fit->chain->atoms, fit->istart, fit->iend, &fit->ox, &fit->oy, &fit->oz);
  ox = fit->ox;
  oy = fit->oy;
  oz = fit->oz;

  num_rows = num_atoms * 6;
  num_cols = ATLS_NUM_PARAMS;

  zero_dmatrix(A, num_rows, num_cols);

  for (ia = istart, row = 0; ia <= iend; ia++, row += 6) {
    set_ATLS_Ab(A, b, num_rows, num_cols, row, atoms[ia].U, atoms[ia].x-ox, atoms[ia].y-oy, atoms[ia].z-oz, atoms[ia].sqrt_weight);
  }

  jobz = 'S';
  dgesdd_(&jobz,  
	  &num_rows,
	  &num_cols, 
	  A, 
	  &num_rows,
	  fit->chain->S, 
	  fit->chain->U,
	  &num_rows, 
	  fit->chain->VT,
	  &num_cols,
	  fit->chain->WORK, 
	  &fit->chain->LWORK,
	  fit->chain->IWORK,
	  &info);

  if (info != 0) {
    printf("DGESDD ERROR(anisotropic): info = %d\n", info);
  }

  solve_SVD(num_rows, num_cols, fit->ATLS, b, fit->chain->U, fit->chain->S, fit->chain->VT, fit->chain->WORK);
}

static void
linear_anisotropic_fit_segment(struct TLSFitContext *fit)
{
  int num_atoms, num_residues, ia, istart, iend, ia_res_start;
  double ox, oy, oz;
  double Utls[6], delta, sum_weight, chi2, tmp; 
  struct Atom *atoms;

  /* optimization */
  atoms  = fit->chain->atoms;
  istart = fit->istart;
  iend   = fit->iend;

  num_atoms = 0;
  num_residues = 1;
  sum_weight = 0.0;
  ia_res_start = istart;

  for (ia = istart; ia <= iend; ia++) {
    /* count the number of new residues */
    if (strcmp(atoms[ia_res_start].frag_id, atoms[ia].frag_id)!=0) {
      num_residues++;
      ia_res_start = ia;
    }
    num_atoms++;
    sum_weight += atoms[ia].sqrt_weight;
  }

  /* fit TLS parameters */
  linear_anisotropic_fit_parameters(fit);

  ox = fit->ox;
  oy = fit->oy;
  oz = fit->oz;

  /* calculate residual */
  chi2 = 0.0;

  for (ia = istart; ia <= iend; ia++) {
    calc_anisotropic_Utls(fit->ATLS, atoms[ia].x-ox, atoms[ia].y-oy, atoms[ia].z-oz, Utls);

    /* calculate a isotropic residual */
    delta = ((Utls[0] + Utls[1] + Utls[2]) / 3.0) - ((atoms[ia].U[0] + atoms[ia].U[1] + atoms[ia].U[2]) / 3.0);
    tmp = atoms[ia].sqrt_weight * delta;
    chi2 += tmp * tmp;
  }

  fit->alsqr = num_residues * (chi2 / sum_weight);
}

/**********************************************************************************************************************************************************/
/* HINGE DETECTION */

static void
calc_isotropic_hinge_delta(struct IHingeContext *hinge)
{
  int num_atoms, ia;
  double u_iso_a, u_iso_b;
  double delta, dela, delb;
  double hdelta_ab, hdelta_abo;
  double hdelta_a, hdelta_b;
  double msd_a, msd_b;

  struct TLSFitContext fita, fitb;
  struct Atom *atoms;

  /* optimization */
  atoms = hinge->chain->atoms;

  fita.chain  = hinge->chain;
  fita.istart = hinge->istarta;
  fita.iend   = hinge->ienda;
  linear_isotropic_fit_segment(&fita);

  fitb.chain  = hinge->chain;
  fitb.istart = hinge->istartb;
  fitb.iend   = hinge->iendb;
  linear_isotropic_fit_segment(&fitb);

  /* Segment A */
  num_atoms = 0;
  hdelta_a = 0.0;

  for (ia = hinge->istarta; ia <= hinge->ienda; ia++) {
    num_atoms++;

    calc_isotropic_uiso(fita.ITLS, atoms[ia].x - fita.ox, atoms[ia].y - fita.oy, atoms[ia].z - fita.oz, &u_iso_a);
    delta = atoms[ia].u_iso - u_iso_a;
    hdelta_a += delta*delta;
  }
  msd_a = hdelta_a / (num_atoms);

  /* Segment B */
  num_atoms = 0;
  hdelta_b = 0.0;

  for (ia = hinge->istartb; ia <= hinge->iendb; ia++) {
    num_atoms++;

    calc_isotropic_uiso(fitb.ITLS, atoms[ia].x - fitb.ox, atoms[ia].y - fitb.oy, atoms[ia].z - fitb.oz, &u_iso_b);
    delta = atoms[ia].u_iso - u_iso_b;
    hdelta_b += delta*delta;
  }
  msd_b = hdelta_b / (num_atoms);


  /* Cross Prediction of A/B */
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

    hdelta_abo += delb*delb;

    delta = u_iso_b - u_iso_a;
    hdelta_ab += delta*delta;
  }

  /* segment b atoms */
  for (ia = hinge->istartb; ia <= hinge->iendb; ia++) {
    num_atoms++;

    calc_isotropic_uiso(fita.ITLS, atoms[ia].x - fita.ox, atoms[ia].y - fita.oy, atoms[ia].z - fita.oz, &u_iso_a);
    calc_isotropic_uiso(fitb.ITLS, atoms[ia].x - fitb.ox, atoms[ia].y - fitb.oy, atoms[ia].z - fitb.oz, &u_iso_b);

    dela = atoms[ia].u_iso - u_iso_a;
    delb = atoms[ia].u_iso - u_iso_b;

    hdelta_abo += dela*dela;

    delta = u_iso_a - u_iso_b;
    hdelta_ab += delta*delta;
  }
  
  hinge->msd_a = msd_a;
  hinge->msd_b = msd_b;
  hinge->msd_c = 0.0;

  hinge->hdelta_ab  = hdelta_ab / (num_atoms);
  hinge->hdelta_abo = hdelta_abo / (num_atoms);
  hinge->hdelta_c   = 0.0;
}

/**********************************************************************************************************************************************************/

/* 
 * PYTHON INTERFACE
 *
 * Two Python classes interface to the high-performance LSQ fitting
 * algorthims:
 *
 * LinearTLSModel: Linear Fit of TLS Models
 */

static PyObject *LINEARTLS_ERROR = NULL;


/* Python interface */
typedef struct {
  PyObject_HEAD
  PyObject       *xmlrpc_chain; /* list of dictionaies describing the atoms */
  struct Chain   *chain;        /* internal version of the chain */
} LinearTLSModel_Object;

static void
LinearTLSModel_dealloc(LinearTLSModel_Object* self)
{
  if (self->chain) {
    delete_chain(self->chain);
    self->chain = NULL;
  }

  Py_XDECREF(self->xmlrpc_chain);
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
LinearTLSModel_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  LinearTLSModel_Object *self;
  
  self = (LinearTLSModel_Object *)type->tp_alloc(type, 0);
  if (self == NULL) {
    return NULL;
  }

  self->xmlrpc_chain = NULL;
  self->chain = NULL;

  return (PyObject *)self;
}

static PyObject *
LinearTLSModel_set_xmlrpc_chain(PyObject *py_self, PyObject *args)
{
  LinearTLSModel_Object *self;
  PyObject *xmlrpc_chain;
  PyObject *atm_desc;
  PyObject *tmp;

  int i, j, num_atoms;
  char *strx;

  self = (LinearTLSModel_Object *) py_self;

  if (!PyArg_ParseTuple(args, "O", &xmlrpc_chain)) {
    goto error;
  }

  /* de-reference the old xmlrpc chain and refrence the new one */
  tmp = self->xmlrpc_chain;
  Py_INCREF(xmlrpc_chain);
  self->xmlrpc_chain = xmlrpc_chain;
  Py_XDECREF(tmp);

  /* free the atoms array */
  if (self->chain!=NULL) {
    delete_chain(self->chain);
    self->chain = NULL;
  }

  /* allocate and fill the new atoms array */
  num_atoms = PyList_Size(self->xmlrpc_chain);
  if (num_atoms > 0) {

    self->chain = new_chain(num_atoms);
    if (self->chain == NULL) {
      PyErr_SetString(LINEARTLS_ERROR, "unable to allocate struct Chain");
      goto error;
    }
    
    for (i = 0; i < num_atoms; i++) {
      atm_desc = PyList_GetItem(self->xmlrpc_chain, i);
      
      /* set name */
      tmp = PyDict_GetItemString(atm_desc, "name");
      if (tmp == NULL) {
	PyErr_SetString(LINEARTLS_ERROR, "name not in atm_desc");
	goto error;
      }
      strx = PyString_AsString(tmp);
      if (strx == NULL) {
	goto error;
      }
      strncpy(self->chain->atoms[i].name, strx, NAME_LEN);
	
      /* set frag_id */
      tmp = PyDict_GetItemString(atm_desc, "frag_id");
      if (tmp == NULL) {
	PyErr_SetString(LINEARTLS_ERROR, "frag_id not in atm_desc");
	goto error;
      }
      strx = PyString_AsString(tmp);
      if (strx == NULL) {
	goto error;
      }
      strncpy(self->chain->atoms[i].frag_id, strx, FRAG_ID_LEN);

      /* set x, y, z coordinates */
      tmp = PyDict_GetItemString(atm_desc, "ifrag");
      if (tmp == NULL) {
	PyErr_SetString(LINEARTLS_ERROR, "ifrag not in atm_desc");
	goto error;
      }
      if (!PyInt_Check(tmp)) {
	goto error;
      }
      self->chain->atoms[i].ifrag = PyInt_AsLong(tmp);

      /* set x, y, z coordinates */
      tmp = PyDict_GetItemString(atm_desc, "x");
      if (tmp == NULL) {
	PyErr_SetString(LINEARTLS_ERROR, "x not in atm_desc");
	goto error;
      }
      if (!PyFloat_Check(tmp)) {
	goto error;
      }
      self->chain->atoms[i].x = PyFloat_AsDouble(tmp);

      tmp = PyDict_GetItemString(atm_desc, "y");
      if (tmp == NULL) {
	PyErr_SetString(LINEARTLS_ERROR, "y not in atm_desc");
	goto error;
      }
      if (!PyFloat_Check(tmp)) {
	goto error;
      }
      self->chain->atoms[i].y = PyFloat_AsDouble(tmp);

      tmp = PyDict_GetItemString(atm_desc, "z");
      if (tmp == NULL) {
	PyErr_SetString(LINEARTLS_ERROR, "z not in atm_desc");
	goto error;
      }
      if (!PyFloat_Check(tmp)) {
	goto error;
      }
      self->chain->atoms[i].z = PyFloat_AsDouble(tmp);

      tmp = PyDict_GetItemString(atm_desc, "u_iso");
      if (tmp == NULL) {
	PyErr_SetString(LINEARTLS_ERROR, "u_iso not in atm_desc");
	goto error;
      }
      if (!PyFloat_Check(tmp)) {
	goto error;
      }
      self->chain->atoms[i].u_iso = PyFloat_AsDouble(tmp);

      /* get U tensor parameters */
      for (j = 0; j < U_NUM_PARAMS; j++) {
	tmp = PyDict_GetItemString(atm_desc, U_PARAM_NAMES[j]);
	if (tmp == NULL) {
	  PyErr_SetString(LINEARTLS_ERROR, "uXX not in atm_desc");
	  goto error;
	}
	if (!PyFloat_Check(tmp)) {
	  goto error;
	}
	self->chain->atoms[i].U[j] = PyFloat_AsDouble(tmp);
      }

      /* weight */
      tmp = PyDict_GetItemString(atm_desc, "weight");
      if (tmp == NULL) {
	PyErr_SetString(LINEARTLS_ERROR, "weight not in atm_desc");
	goto error;
      }
      if (!PyFloat_Check(tmp)) {
	goto error;
      }
      self->chain->atoms[i].sqrt_weight = sqrt(PyFloat_AsDouble(tmp));
    }
  }

  Py_INCREF(Py_None);
  return Py_None;

 error:
  return NULL;
}

static PyObject *
LinearTLSModel_isotropic_fit_segment(PyObject *py_self, PyObject *args)
{
  LinearTLSModel_Object *self;
  PyObject *py_floatx, *rdict;

  int i;
  struct TLSFitContext fit_context;

  self = (LinearTLSModel_Object *) py_self;

  /* fill in fields in the TLSFitContext structure */
  if (self->chain==NULL) {
    goto error;
  }
  fit_context.chain = self->chain;
  fit_context.ox = 0.0;
  fit_context.oy = 0.0;
  fit_context.oz = 0.0;

  for (i = 0; i <= ITLS_NUM_PARAMS; i++) {
    fit_context.ITLS[i] = 0.0;
  }
  for (i = 0; i <= ATLS_NUM_PARAMS; i++) {
    fit_context.ATLS[i] = 0.0;
  }
  fit_context.ilsqr = 0.0;

  if (!PyArg_ParseTuple(args, "ii", &fit_context.istart, &fit_context.iend)) {
    goto error;
  }

  /* fit the segment */
  linear_isotropic_fit_segment(&fit_context);

  /* construct return dictioary with results */
  rdict = PyDict_New();

  /* set minimization exit status */
  py_floatx = PyFloat_FromDouble(fit_context.ox);
  PyDict_SetItemString(rdict, "x", py_floatx);
  Py_DECREF(py_floatx);

  py_floatx = PyFloat_FromDouble(fit_context.oy);
  PyDict_SetItemString(rdict, "y", py_floatx);
  Py_DECREF(py_floatx);

  py_floatx = PyFloat_FromDouble(fit_context.oz);
  PyDict_SetItemString(rdict, "z", py_floatx);
  Py_DECREF(py_floatx);

  py_floatx = PyFloat_FromDouble(fit_context.ilsqr);
  PyDict_SetItemString(rdict, "ilsqr", py_floatx);
  Py_DECREF(py_floatx);
  
  for (i = 0; i < ITLS_NUM_PARAMS; i++) {
    py_floatx = PyFloat_FromDouble(fit_context.ITLS[i]);
    PyDict_SetItemString(rdict, ITLS_PARAM_NAMES[i], py_floatx);
    Py_DECREF(py_floatx);
  }
    
  return rdict;

 error:
  return NULL;
}

static PyObject *
LinearTLSModel_anisotropic_fit_segment(PyObject *py_self, PyObject *args)
{
  LinearTLSModel_Object *self;
  PyObject *py_floatx, *rdict;

  int i;
  struct TLSFitContext fit_context;

  self = (LinearTLSModel_Object *) py_self;

  /* fill in fields in the TLSFitContext structure */
  if (self->chain==NULL) {
    goto error;
  }

  fit_context.chain = self->chain;
  fit_context.ox = 0.0;
  fit_context.oy = 0.0;
  fit_context.oz = 0.0;

  for (i = 0; i <= ITLS_NUM_PARAMS; i++) {
    fit_context.ITLS[i] = 0.0;
  }
  for (i = 0; i <= ATLS_NUM_PARAMS; i++) {
    fit_context.ATLS[i] = 0.0;
  }
  fit_context.ilsqr = 0.0;

  if (!PyArg_ParseTuple(args, "ii", &fit_context.istart, &fit_context.iend)) {
    goto error;
  }

  /* fit the segment */
  linear_anisotropic_fit_segment(&fit_context);

  /* construct return dictioary with results */
  rdict = PyDict_New();

  /* set minimization exit status */
  py_floatx = PyFloat_FromDouble(fit_context.ox);
  PyDict_SetItemString(rdict, "x", py_floatx);
  Py_DECREF(py_floatx);

  py_floatx = PyFloat_FromDouble(fit_context.oy);
  PyDict_SetItemString(rdict, "y", py_floatx);
  Py_DECREF(py_floatx);

  py_floatx = PyFloat_FromDouble(fit_context.oz);
  PyDict_SetItemString(rdict, "z", py_floatx);
  Py_DECREF(py_floatx);

  py_floatx = PyFloat_FromDouble(fit_context.alsqr);
  PyDict_SetItemString(rdict, "alsqr", py_floatx);
  Py_DECREF(py_floatx);
  
  for (i = 0; i < ATLS_NUM_PARAMS; i++) {
    py_floatx = PyFloat_FromDouble(fit_context.ATLS[i]);
    PyDict_SetItemString(rdict, ATLS_PARAM_NAMES[i], py_floatx);
    Py_DECREF(py_floatx);
  }
    
  return rdict;

 error:
  return NULL;
}

static PyObject *
LinearTLSModel_calc_isotropic_hinge_delta(PyObject *py_self, PyObject *args)
{
  LinearTLSModel_Object *self;
  PyObject *py_floatx, *rdict;
  struct IHingeContext hinge;

  self = (LinearTLSModel_Object *) py_self;

  /* fill in fields in the TLSFitContext structure */
  if (self->chain==NULL) {
    goto error;
  }

  hinge.chain  = self->chain;
  hinge.hdelta_ab = 0.0;
  hinge.hdelta_abo = 0.0;
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

static PyObject *
LinearTLSModel_gammaq(PyObject *py_self, PyObject *args)
{
  LinearTLSModel_Object *self;
  PyObject *py_floatx;
  double a, x, gammaq_val;

  self = (LinearTLSModel_Object *) py_self;

  if (!PyArg_ParseTuple(args, "dd", &a, &x)) {
    goto error;
  }

  printf("a=%f x=%f\n",a, x);
  
  gammaq(a, x, &gammaq_val);

  py_floatx = PyFloat_FromDouble(gammaq_val);
  return py_floatx;

 error:
  return NULL;
}

static PyMethodDef LinearTLSModel_methods[] = {
    {"set_xmlrpc_chain", 
     (PyCFunction) LinearTLSModel_set_xmlrpc_chain, 
     METH_VARARGS,
     "Sets the Python list containing one dictionary for each atom." },

    {"isotropic_fit_segment",
     (PyCFunction) LinearTLSModel_isotropic_fit_segment, 
     METH_VARARGS,
     "Performs a linear fit of the isotropic TLS model to the given atoms." },

    {"anisotropic_fit_segment",
     (PyCFunction) LinearTLSModel_anisotropic_fit_segment, 
     METH_VARARGS,
     "Performs a linear fit of the anisotropic TLS model to the given atoms." },

    {"calc_isotropic_hinge_delta",
     (PyCFunction) LinearTLSModel_calc_isotropic_hinge_delta, 
     METH_VARARGS,
     "Calculates the value of the hinge-delta function." },

    {"gammaq",
     (PyCFunction) LinearTLSModel_gammaq, 
     METH_VARARGS,
     "Calculates the incomplete gamma function Q." },

    {NULL}  /* Sentinel */
};

static PyTypeObject LinearTLSModel_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "LinearTLSModel",          /*tp_name*/
    sizeof(LinearTLSModel_Object), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)LinearTLSModel_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "LinearTLSModel objects",  /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    LinearTLSModel_methods,    /* tp_methods */
    0,                         /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    0,                         /* tp_init */
    0,                         /* tp_alloc */
    LinearTLSModel_new,        /* tp_new */
};


static PyMethodDef LINEARTLS_METHODS[] = {
  {NULL, NULL, 0, NULL}
};

DL_EXPORT(void)
initlineartls(void)
{
  PyObject *m;
  
  
  if (PyType_Ready(&LinearTLSModel_Type) < 0)
    return;

  m = Py_InitModule("lineartls", LINEARTLS_METHODS);
  
  LINEARTLS_ERROR = PyErr_NewException("lineartls.error", NULL, NULL);
  Py_INCREF(LINEARTLS_ERROR);
  PyModule_AddObject(m, "error", LINEARTLS_ERROR);


  /* add the LinearTLSModel class */
  Py_INCREF(&LinearTLSModel_Type);
  PyModule_AddObject(m, "LinearTLSModel", (PyObject *)&LinearTLSModel_Type);
}
