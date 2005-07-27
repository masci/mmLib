/* tlsmdmodule.c -- fast calculation of TLS segments of a Protein/DNA/RNA
 *                  polymer chain
 *
 */
#include "Python.h"
#include "structmember.h"

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <pthread.h>


/* macros */
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))

/* set to 1 for extra debugging informaiton */
#define _DEBUG 0


/* 
 * Anisotropic ADP Parameters: U
 */

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


/*
 * Isotropic TLS Model
 */
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
static inline void
set_ITLS_Ab(double *A, double *b, int m, int n, int row,
	    double uiso, double x, double y, double z, double w)
{
#define FA(_m, _n) A[_m + (m * _n)]

  int i, sz;
  double xx, yy, zz, xy, xz, yz;
   
  xx = x*x;
  yy = y*y;
  zz = z*z;
  xy = x*y;
  xz = x*z;
  yz = y*z;

  /* array size sz */
  sz = m * n;

  /* zero out array A */
  for (i = 0; i<=sz; i++) {
    A[i] = 0.0;
  }

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



/* 
 * Anisotropic TLS Model
 */

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

/* Sets the six rows of matrix A starting at A[i,j] with the anistropic
 * TLS model coefficents for a atom located at t position x, y, z with
 * least-squares weight w.  Matrix A is filled to coumn j+12.
 *
 * The matrix A(m,n) is filled in FORTRAN-style, that is, assuming
 * a column-major memory layout.  That's because, athough the
 * matrix is contructed in C, the SVD subroutine from LAPACK is 
 * written in FORTRAN.
 */
static inline void
set_ATLS_Ab(double *A, double *b, int m, int n, int row,
	    double U[6], double x, double y, double z, double w)
{
#define FA(_m, _n) A[_m + (m * _n)]

  int i, t, sz, rowU11, rowU22, rowU33, rowU12, rowU13, rowU23;
  double xx, yy, zz, xy, xz, yz;
  
  xx = x*x;
  yy = y*y;
  zz = z*z;
  xy = x*y;
  xz = x*z;
  yz = y*z;

  /* array size sz */
  sz = m * n;

  /* zero out array A */
  for (i = 0; i<=sz; i++) {
    A[i] = 0.0;
  }

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
}


/* calculate the anistropic TLS model prediction of anisiotropic
 * ADP U at position x, y, z
 */
static void
calc_ATLS_U(double ATLS[ATLS_NUM_PARAMS], 
	    double x, double y, double z, 
	    double U[6]) 
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


/* Common routines for the isotropic and ansiotropic TLS segment
 * fitting engine
 */
#define NAME_LEN     8
#define FRAG_ID_LEN  8

/* structure used to store the information of one atom */
struct Atom {
  char    name[NAME_LEN];
  char    frag_id[FRAG_ID_LEN];
  double  x;
  double  y;
  double  z;
  double  u_iso;
  double  U[6];
  double  weight;
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
  double *x;
  double *b;
  double *S;
  double *U;
  double *VT;
  double *WORK;
  int     LWORK;
  int    *IWORK;
};

/* frees all memory from a allocated struct Chain 
 */
void
delete_chain(struct Chain *chain)
{
  if (chain->atoms != NULL)
    free(chain->atoms);
  if (chain->A != NULL)
    free(chain->A);
  if (chain->x != NULL)
    free(chain->x);
  if (chain->b != NULL)
    free(chain->b);
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
struct Chain *
new_chain(int num_atoms)
{
  int num_rows, num_cols, nu, nvt;
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
  chain->A = NULL;
  chain->x = NULL;
  chain->b = NULL;
  chain->S = NULL;
  chain->U = NULL;
  chain->VT = NULL;
  chain->WORK = NULL;
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

  chain->x = malloc(sizeof(double) * num_cols);
  if (chain->x == NULL) {
    printf("new_chain: chain->x\n");
    goto error;
  }

  chain->b = malloc(sizeof(double) * num_rows);
  if (chain->b == NULL) {
    printf("new_chain: chain->b\n");
    goto error;
  }

  chain->S = malloc(sizeof(double) * MIN(num_rows, num_cols));
  if (chain->S == NULL) {
    printf("new_chain: chain->S\n");
    goto error;
  }

  /* nu = min(n,m)
   * u = Numeric.zeros((nu, m), t)
   */
  nu = MIN(num_rows, num_cols);
  chain->U = malloc(sizeof(double) * (nu * num_rows));
  if (chain->U == NULL) {
    printf("new_chain: chain->U\n");
    goto error;
  }
  
  /* nvt = min(n,m)
   * vt = Numeric.zeros((n, nvt), t)
   */
  nvt = MIN(num_rows, num_cols);
  chain->VT = malloc(sizeof(double) * (num_cols * nvt));
  if (chain->VT == NULL) {
    printf("new_chain: chain->VT\n");
    goto error;
  }

  /* iwork = Numeric.zeros((8*min(m,n),), 'i')
   */
  chain->IWORK = malloc(sizeof(int) * (8 * MIN(num_rows, num_cols)));
  if (chain->IWORK == NULL) {
    printf("new_chain: chain->IWORK\n");
    goto error;
  }

  /* calculate the size of the WORK memory block */

  return chain;

  /* if malloc() fails */
 error:
  if (chain != NULL) {
    delete_chain(chain);
  }
  
  return NULL;
}


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


/* calculates the centroid of the atoms indexed between istart and iend
 * then returns the centroid coordinates in x, y, z
 */
void
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

void
fit_segment_ITLS(struct ITLSFitContext *itls_context)
{
  int num_atoms, num_rows, num_cols, row, ia, istart, iend;
  double ox, oy, oz, lsqr;

  double *A, *b;
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

  A = itls_context->chain->A;
  b = itls_context->chain->b;

  ox = itls_context->origin_x;
  oy = itls_context->origin_y;
  oz = itls_context->origin_z;

#ifdef _DEBUG
  printf("fit_segment_ITLS(centrod=(%f, %f, %f))\n", ox, oy, oz);
#endif

  /* fill in coefficent matrix for the isotropic TLS model */
  num_rows = num_atoms;
  num_cols = ITLS_NUM_PARAMS;

  for (ia = istart, row = 0; ia <= iend; ia++, row++) {
    set_ITLS_Ab(A, 
		b,
		num_rows,
		num_cols, 
		row,
		atoms[ia].u_iso,
		atoms[ia].x - ox,
		atoms[ia].y - oy,
		atoms[ia].z - oz,
		atoms[ia].weight);
  }

  /* solve for isotropic TLS model */



  /* fill in coefficent matrix for the anisotropic TLS model */
  num_rows = num_atoms * 6;
  num_cols = ATLS_NUM_PARAMS;

  for (ia = istart, row = 0; ia <= iend; ia++, row += 6) {
    set_ATLS_Ab(A, 
		b,
		num_rows,
		num_cols,
		row,
		atoms[ia].U,
		atoms[ia].x - ox,
		atoms[ia].y - oy,
		atoms[ia].z - oz,
		atoms[ia].weight);
  }

  /* solve for isotropic TLS model */

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

static PyObject *TLSMDMODULE_ERROR = NULL;


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
  double xx;

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
      PyErr_SetString(TLSMDMODULE_ERROR, "unable to allocate struct Chain");
      goto error;
    }
    
    for (i = 0; i < num_atoms; i++) {
      atm_desc = PyList_GetItem(self->xmlrpc_chain, i);
      
      /* set name */
      tmp = PyDict_GetItemString(atm_desc, "name");
      if (tmp == NULL) {
	PyErr_SetString(TLSMDMODULE_ERROR, "name not in atm_desc");
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
	PyErr_SetString(TLSMDMODULE_ERROR, "frag_id not in atm_desc");
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
	PyErr_SetString(TLSMDMODULE_ERROR, "x not in atm_desc");
	goto error;
      }
      if (!PyFloat_Check(tmp)) {
	goto error;
      }
      self->chain->atoms[i].x = PyFloat_AsDouble(tmp);

      tmp = PyDict_GetItemString(atm_desc, "y");
      if (tmp == NULL) {
	PyErr_SetString(TLSMDMODULE_ERROR, "y not in atm_desc");
	goto error;
      }
      if (!PyFloat_Check(tmp)) {
	goto error;
      }
      self->chain->atoms[i].y = PyFloat_AsDouble(tmp);

      tmp = PyDict_GetItemString(atm_desc, "z");
      if (tmp == NULL) {
	PyErr_SetString(TLSMDMODULE_ERROR, "z not in atm_desc");
	goto error;
      }
      if (!PyFloat_Check(tmp)) {
	goto error;
      }
      self->chain->atoms[i].z = PyFloat_AsDouble(tmp);

      tmp = PyDict_GetItemString(atm_desc, "u_iso");
      if (tmp == NULL) {
	PyErr_SetString(TLSMDMODULE_ERROR, "u_iso not in atm_desc");
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
	  PyErr_SetString(TLSMDMODULE_ERROR, "uXX not in atm_desc");
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
	PyErr_SetString(TLSMDMODULE_ERROR, "sqrt_w not in atm_desc");
	goto error;
      }
      if (!PyFloat_Check(tmp)) {
	goto error;
      }
      self->chain->atoms[i].weight = PyFloat_AsDouble(tmp);

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
  PyObject *py_floatx, *py_intx, *rdict;

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


static PyMethodDef TLSMDMODULE_METHODS[] = {
  {NULL, NULL, 0, NULL}
};

DL_EXPORT(void)
inittlsmdmodule(void)
{
  PyObject *m;
  
  
  if (PyType_Ready(&ITLSModel_Type) < 0)
    return;

  m = Py_InitModule("tlsmdmodule", TLSMDMODULE_METHODS);
  
  TLSMDMODULE_ERROR = PyErr_NewException("tlsmdmodule.error", NULL, NULL);
  Py_INCREF(TLSMDMODULE_ERROR);
  PyModule_AddObject(m, "error", TLSMDMODULE_ERROR);


  /* add the ITLSModel class */
  Py_INCREF(&ITLSModel_Type);
  PyModule_AddObject(m, "ITLSModel", (PyObject *)&ITLSModel_Type);
}
