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
#include <math.h>
#include <malloc.h>
#include <pthread.h>


/* LAPACK */
extern void
dgesdd_(char *, int*, int*, double*, int*, double*, double*, int *, double*, int*, double*, int*, int*, int*);


/* macros */
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))

/* set to 1 for extra debugging informaiton */
/* #define _DEBUG 0 */



/* Anisotropic ADP Parameters: U */

/* anisotropic U tensor parameter labels and indexes */
#define U11 0
#define U22 1
#define U33 2
#define U12 3
#define U13 4
#define U23 5

#define U_NUM_PARAMS 6

/* parameter name/labels used when they are passed in through 
 * Python dictionaries 
 */
static char *U_PARAM_NAMES[] = {
  "u11", "u22", "u33", "u12", "u13", "u23"
};


/* Isotropic TLS Model */
#define ITLS_T     0
#define ITLS_L11   1
#define ITLS_L22   2
#define ITLS_L33   3
#define ITLS_L12   4
#define ITLS_L13   5
#define ITLS_L23   6
#define ITLS_S12   7
#define ITLS_S21   8
#define ITLS_S13   9
#define ITLS_S31   10
#define ITLS_S23   11
#define ITLS_S32   12
#define ITLS_NUM_PARAMS 13
static char *ITLS_PARAM_NAMES[] = {
  "t",
  "l11", "l22", "l33", "l12", "l13", "l23",
  "s12", "s21", "s13", "s31", "s23", "s32"
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
  char    name[NAME_LEN];
  char    frag_id[FRAG_ID_LEN];
  double  x;
  double  y;
  double  z;
  double  u_iso;
  double  U[6];
  double  sqrt_weight;
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
  double *Aw;
  double *b;
  double *bw;
  double *bc;
  double *S;
  double *U;
  double *VT;
  double *WORK;
  int     LWORK;
  int    *IWORK;
};


/* context structure for fitting one TLS segment using the 
 * isotropic TLS model
 */
struct ITLSFitContext {
  struct Chain      *chain;                 /* pointer to Chain structure */

  int                istart;                /* index of first atom in 
					     * chain->atoms */
  int                iend;                  /* index of last atom in 
					     * chain->atoms */

  double             origin_x;              /* origin of TLS tensors */
  double             origin_y;
  double             origin_z;

  double             ITLS[ITLS_NUM_PARAMS]; /* isotropic TLS model params */
  double             ATLS[ATLS_NUM_PARAMS]; /* ansiotropic TLS model params */
  double             ilsqr;                 /* least-squares residual of
					     *  isotropic TLS model */
};


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

  FA(row, ITLS_S12) = w * ((-2.0 * z) / 3.0);
  FA(row, ITLS_S21) = w * (( 2.0 * z) / 3.0);

  FA(row, ITLS_S13) = w * (( 2.0 * y) / 3.0);
  FA(row, ITLS_S31) = w * ((-2.0 * y) / 3.0);

  FA(row, ITLS_S23) = w * ((-2.0 * x) / 3.0);
  FA(row, ITLS_S32) = w * (( 2.0 * x) / 3.0);
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
solve_SVD(int m,
	  int n,
	  double *x,
	  double *b,
	  double *U,
	  double *S,
	  double *VT,
	  double *WORK)
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
  scutoff = smax * 1E-12;
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


/* calculates bc = A*x and least-squares residual (b - bc)*(b-bc)t
 */
static void
calc_bc_lsqr(double *A,
	     int m,
	     int n,
	     double *x,
	     double *b,
	     double *bc,
	     double *lsqr)
{
#define FA(__i,__j)   A[__i + (m * __j)]

  int i, j;
  double dtmp, dltmp;

  dltmp = 0.0;

  /* now we have solved for x; calculate bc (b-calc) */
  for (i = 0; i < m; i++) {
    dtmp = 0.0;
    
    for (j = 0; j < n; j++) {
      dtmp += FA(i,j) * x[j];
    }
    
    bc[i] = dtmp;

    /* sum the lsqr residual (b - bc) */
    dtmp = b[i] - bc[i];
    dltmp += dtmp * dtmp;
  }

  /* return the least squares residual in lsqr */
  *lsqr = dltmp;

#undef FA
}


/* calculates bc = A*x and least-squares residual (b - bc)*(b-bc)t
 */
static void
calc_lsqr(double *A,
	  int m,
	  int n,
	  double *x,
	  double *b,
	  double *lsqr)
{
#define FA(__i,__j)   A[__i + (m * __j)]

  int i, j;
  double bcalc, dtmp, dltmp;

  dltmp = 0.0;

  /* now we have solved for x; calculate bc (b-calc) */
  for (i = 0; i < m; i++) {
    
    bcalc = 0.0;
    for (j = 0; j < n; j++) {
      bcalc += FA(i,j) * x[j];
    }

    /* sum the lsqr residual (b - bc) */
    dtmp = b[i] - bcalc;
    dltmp += dtmp * dtmp;
  }

  /* return the least squares residual in lsqr */
  *lsqr = dltmp;

#undef FA
}


/* frees all memory from a allocated struct Chain 
 */
static void
delete_chain(struct Chain *chain)
{
  if (chain->atoms != NULL)
    free(chain->atoms);
  if (chain->A != NULL)
    free(chain->A);
  if (chain->Aw != NULL)
    free(chain->Aw);
  if (chain->b != NULL)
    free(chain->b);
  if (chain->bw != NULL)
    free(chain->bw);
  if (chain->bc != NULL)
    free(chain->bc);
  if (chain->S != NULL)
    free(chain->S);
  if (chain->U != NULL)
    free(chain->U);
  if (chain->VT != NULL)
    free(chain->VT);
  if (chain->WORK != NULL)
    free(chain->WORK);
  if (chain->IWORK != NULL)
    free(chain->IWORK);
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
  chain->Aw    = NULL;
  chain->b     = NULL;
  chain->bw    = NULL;
  chain->bc    = NULL;
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

  chain->Aw = malloc(sizeof(double) * (num_rows * num_cols));
  if (chain->Aw == NULL) {
    printf("new_chain: chain->Aw\n");
    goto error;
  }

  chain->b = malloc(sizeof(double) * num_rows);
  if (chain->b == NULL) {
    printf("new_chain: chain->b\n");
    goto error;
  }

  chain->bw = malloc(sizeof(double) * num_rows);
  if (chain->bw == NULL) {
    printf("new_chain: chain->bw\n");
    goto error;
  }

  chain->bc = malloc(sizeof(double) * num_rows);
  if (chain->bc == NULL) {
    printf("new_chain: chain->bc\n");
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

#ifdef _DEBUG
  printf("new_chain: ideal WORK size %d\n", chain->LWORK);
#endif

  chain->WORK = malloc(sizeof(double) * chain->LWORK);
  if (chain->WORK == NULL) {
    printf("new_chain: chain->WORK\n");
    goto error;
  }

  return chain;

  /* if malloc() fails */
 error:
  if (chain != NULL) {
    delete_chain(chain);
  }
  
  return NULL;
}


/* calculates the centroid of the atoms indexed between istart and iend
 * then returns the centroid coordinates in x, y, z
 */
static void
calc_centroid(struct Atom *atoms, int istart, int iend, 
	      double *x, double *y, double *z)
{
  int i;
  double n, cx, cy, cz;

  n  = 0.0;
  cx = 0.0;
  cy = 0.0;
  cz = 0.0;

  for (i = istart; i<=iend; i++) {
    n  += 1.0;
    cx += atoms[i].x;
    cy += atoms[i].y;
    cz += atoms[i].z;
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

static void
fit_segment_ITLS(struct ITLSFitContext *itls_context)
{
  char jobz;
  int i, sz, num_atoms, num_rows, num_cols, row, ia, istart, iend, info;
  int natmp;

  double origin_x, origin_y, origin_z;

  double *A, *Aw, *b, *bw;
  struct Atom *atoms;

  /* calculate the number of atoms to be fit */
  num_atoms = itls_context->iend - itls_context->istart + 1;

  /* calculate the centroid of the atoms which are to
   * be fit and use it as the origin of the TLS tensors
   */
  calc_centroid(itls_context->chain->atoms,
		itls_context->istart,
		itls_context->iend,
		&itls_context->origin_x,
		&itls_context->origin_y,
		&itls_context->origin_z);

  /* optimization */
  atoms = itls_context->chain->atoms;

  istart = itls_context->istart;
  iend   = itls_context->iend;

  A  = itls_context->chain->A;
  Aw = itls_context->chain->Aw;
  b  = itls_context->chain->b;
  bw = itls_context->chain->bw;

  origin_x = itls_context->origin_x;
  origin_y = itls_context->origin_y;
  origin_z = itls_context->origin_z;


  /*
   * ISOTROPIC TLS MODEL
   */
  num_rows = num_atoms;
  num_cols = ITLS_NUM_PARAMS;

  zero_dmatrix(Aw, num_rows, num_cols);

  for (ia = istart, row = 0; ia <= iend; ia++, row++) {
    set_ITLS_Ab(Aw,
		bw,
		num_rows,
		num_cols, 
		row,
		atoms[ia].u_iso,
		atoms[ia].x - origin_x,
		atoms[ia].y - origin_y,
		atoms[ia].z - origin_z,
		atoms[ia].sqrt_weight);
  }

  /* make a copy of Aw into A because dgesdd_ destroys Aw */
  sz = num_rows * num_cols;
  for (i = 0; i < sz; i++) {
    A[i] = Aw[i];
  }

  /* solve for isotropic TLS model */
  jobz = 'S';
  dgesdd_(&jobz,  
	  &num_rows,
	  &num_cols, 
	  itls_context->chain->Aw,
	  &num_rows,
	  itls_context->chain->S, 
	  itls_context->chain->U,
	  &num_rows,
	  itls_context->chain->VT,
	  &num_cols,
	  itls_context->chain->WORK, 
	  &itls_context->chain->LWORK,
	  itls_context->chain->IWORK,
	  &info);

  if (info != 0) {
    printf("DGESDD ERROR(isotropic): info = %d\n", info);
  }

  solve_SVD(num_rows,
	    num_cols,
	    itls_context->ITLS,
	    itls_context->chain->bw,
	    itls_context->chain->U,
	    itls_context->chain->S,
	    itls_context->chain->VT,
	    itls_context->chain->WORK);

  /* now we have solved for x; now we can calculate the 
   * weighted residual; remember that A is actually a copy of Aw
   * residual */
  calc_lsqr(itls_context->chain->A, 
	    num_rows, 
	    num_cols,
	    itls_context->ITLS,
	    itls_context->chain->bw,
	    &itls_context->ilsqr);


  /*
   * ANISOTROPIC TLS MODEL
   */
  num_rows = num_atoms * 6;
  num_cols = ATLS_NUM_PARAMS;

  zero_dmatrix(Aw, num_rows, num_cols);

  natmp = 0;
  for (ia = istart, row = 0; ia <= iend; ia++, row += 6) {
    natmp += 1;

    set_ATLS_Ab(Aw,
		bw,
		num_rows,
		num_cols,
		row,
		atoms[ia].U,
		atoms[ia].x - origin_x,
		atoms[ia].y - origin_y,
		atoms[ia].z - origin_z,
		atoms[ia].sqrt_weight);
  }

  /* solve for isotropic TLS model */
  jobz = 'S';
  dgesdd_(&jobz,  
	  &num_rows,
	  &num_cols, 
	  itls_context->chain->Aw, 
	  &num_rows,
	  itls_context->chain->S, 
	  itls_context->chain->U,
	  &num_rows, 
	  itls_context->chain->VT,
	  &num_cols,
	  itls_context->chain->WORK, 
	  &itls_context->chain->LWORK,
	  itls_context->chain->IWORK,
	  &info);

  if (info != 0) {
    printf("DGESDD ERROR(anisotropic): info = %d\n", info);
  }

  solve_SVD(num_rows, 
	    num_cols,
	    itls_context->ATLS,
	    itls_context->chain->bw,
	    itls_context->chain->U,
	    itls_context->chain->S,
	    itls_context->chain->VT,
	    itls_context->chain->WORK);

#ifdef _DEBUG
  printf("fit_segment_ITLS(num_atoms=%d lsqr=%f)\n", 
	 natmp,
	 itls_context->ilsqr);
#endif
}



/* 
 * PYTHON INTERFACE
 *
 * Two Python classes interface to the high-performance LSQ fitting
 * algorthims:
 *
 * ITLSModel: Isotropic TLS Model
 * ATLSModel: Anisotropic TLS Model 
 */

static PyObject *LINEARTLS_ERROR = NULL;


/* Python interface */
typedef struct {
  PyObject_HEAD
  PyObject       *xmlrpc_chain; /* list of dictionaies describing the atoms */
  struct Chain   *chain;        /* internal version of the chain */
} ITLSModel_Object;

static void
ITLSModel_dealloc(ITLSModel_Object* self)
{
  if (self->chain) {
    delete_chain(self->chain);
    self->chain = NULL;
  }

  Py_XDECREF(self->xmlrpc_chain);
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
ITLSModel_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  ITLSModel_Object *self;
  
  self = (ITLSModel_Object *)type->tp_alloc(type, 0);
  if (self == NULL) {
    return NULL;
  }

  self->xmlrpc_chain = NULL;
  self->chain = NULL;

  return (PyObject *)self;
}

static PyObject *
ITLSModel_set_xmlrpc_chain(PyObject *py_self, PyObject *args)
{
  ITLSModel_Object *self;
  PyObject *xmlrpc_chain;
  PyObject *atm_desc;
  PyObject *tmp;

  int i, j, num_atoms;
  char *strx;

  self = (ITLSModel_Object *) py_self;

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
      tmp = PyDict_GetItemString(atm_desc, "sqrt_w");
      if (tmp == NULL) {
	PyErr_SetString(LINEARTLS_ERROR, "sqrt_w not in atm_desc");
	goto error;
      }
      if (!PyFloat_Check(tmp)) {
	goto error;
      }
      self->chain->atoms[i].sqrt_weight = PyFloat_AsDouble(tmp);

#ifdef _DEBUG
      printf("ATOMS[%d]: %s %s\n", 
	     i, 
	     self->chain->atoms[i].name, 
	     self->chain->atoms[i].frag_id);
#endif /* _DEBUG */

    }
  }

  Py_INCREF(Py_None);
  return Py_None;

 error:
  return NULL;
}

static PyObject *
ITLSModel_fit_segment(PyObject *py_self, PyObject *args)
{
  ITLSModel_Object *self;
  PyObject *py_floatx, *rdict;

  int i;
  struct ITLSFitContext fit_context;

  self = (ITLSModel_Object *) py_self;

  /* fill in fields in the ITLSFitContext structure */
  if (self->chain==NULL) {
    goto error;
  }
  fit_context.chain = self->chain;
  fit_context.origin_x = 0.0;
  fit_context.origin_y = 0.0;
  fit_context.origin_z = 0.0;
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
#ifdef _DEBUG
  printf("ITLSModel_fit_segment(istart=%d, iend=%d)\n",
	 fit_context.istart,
	 fit_context.iend);
#endif 
  fit_segment_ITLS(&fit_context);

  /* construct return dictioary with results */
  rdict = PyDict_New();

  /* set minimization exit status */
  py_floatx = PyFloat_FromDouble(fit_context.origin_x);
  PyDict_SetItemString(rdict, "x", py_floatx);
  Py_DECREF(py_floatx);

  py_floatx = PyFloat_FromDouble(fit_context.origin_y);
  PyDict_SetItemString(rdict, "y", py_floatx);
  Py_DECREF(py_floatx);

  py_floatx = PyFloat_FromDouble(fit_context.origin_z);
  PyDict_SetItemString(rdict, "z", py_floatx);
  Py_DECREF(py_floatx);

  py_floatx = PyFloat_FromDouble(fit_context.ilsqr);
  PyDict_SetItemString(rdict, "ilsqr", py_floatx);
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

static PyMethodDef ITLSModel_methods[] = {
    {"set_xmlrpc_chain", 
     (PyCFunction) ITLSModel_set_xmlrpc_chain, 
     METH_VARARGS,
     "Sets the Python list containing one dictionary for each atom." },

    {"fit_segment",
     (PyCFunction) ITLSModel_fit_segment, 
     METH_VARARGS,
     "Performs a TLS/ISO fit to the given atoms." },

    {NULL}  /* Sentinel */
};

static PyTypeObject ITLSModel_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "ITLSModel",          /*tp_name*/
    sizeof(ITLSModel_Object), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)ITLSModel_dealloc, /*tp_dealloc*/
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
    "ITLSModel objects",  /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    ITLSModel_methods,    /* tp_methods */
    0,                         /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    0,                         /* tp_init */
    0,                         /* tp_alloc */
    ITLSModel_new,        /* tp_new */
};


static PyMethodDef LINEARTLS_METHODS[] = {
  {NULL, NULL, 0, NULL}
};

DL_EXPORT(void)
initlineartls(void)
{
  PyObject *m;
  
  
  if (PyType_Ready(&ITLSModel_Type) < 0)
    return;

  m = Py_InitModule("lineartls", LINEARTLS_METHODS);
  
  LINEARTLS_ERROR = PyErr_NewException("lineartls.error", NULL, NULL);
  Py_INCREF(LINEARTLS_ERROR);
  PyModule_AddObject(m, "error", LINEARTLS_ERROR);


  /* add the ITLSModel class */
  Py_INCREF(&ITLSModel_Type);
  PyModule_AddObject(m, "ITLSModel", (PyObject *)&ITLSModel_Type);
}
