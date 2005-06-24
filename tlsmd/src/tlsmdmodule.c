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
/* #define _DEBUG 0 */


static PyObject *TLSMDMODULE_ERROR = NULL;


/*
 * Mathmatical Constants
 */

#define PI     3.1415926535897931
#define PI2    (PI * PI)
#define PI3    (PI * PI * PI)

#define RAD2DEG  (180.0   / PI)
#define RAD2DEG2 (180.0*180.0 / PI2)
#define DEG2RAD  (PI / 180.0)
#define DEG2RAD2 (PI2 / (180.0 * 180.0)


/* normalize the vector v
 */
static void 
normalize(double v[3])
{
  double d;

  d    = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  v[0] = v[0] / d;
  v[1] = v[1] / d;
  v[2] = v[2] / d;
}

/* compute the length of the vector v
 */
static double
length(double v[3])
{
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

/* form the cross product of vectors u and v, and return the
 * result in w
 */
static void 
cross(double u[3], double v[3], double w[3])
{
  w[0] = u[1]*v[2] - u[2]*v[1];
  w[1] = u[2]*v[0] - u[0]*v[2];
  w[2] = u[0]*v[1] - u[1]*v[0];
}

/* calculate the determinate of a symmetric 3x3 matrix
 */
static double
det_symmetric_3(double U[6])
{
  return - U[4]*U[4]*U[1] + 2.0*U[3]*U[4]*U[5] - U[0]*U[5]*U[5]
         - U[3]*U[3]*U[2] + U[0]*U[1]*U[2];
}

/* invert the symmetric 3x3 matrix in the argument U, and return the
 * result in Ui 
 * matrix format: u11,u22,u33,u12,u13,u23
 */
static void
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

static int
invert_symmetric_3(double U[6], double Ui[6])
{
  double d;
  double I[6];

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
static int
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


/* TLS_ISO_Solver: Isotropic TLS Model solver object 
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


void
calc_tls_dT11()
{
}

void
calc_u_tls(double *TLS, double x, double y, double z, double U[6]) 
{
  double xx, yy, zz, xy, yz, xz;

  xx = x * x;
  yy = y * y;
  zz = z * z;
  xy = x * y;
  yz = y * z;
  xz = x * z;
    
  U[U11] =         
            TLS[ATLS_T11]
    +       TLS[ATLS_L22] * zz
    +       TLS[ATLS_L33] * yy
    - 2.0 * TLS[ATLS_L23] * yz
    - 2.0 * TLS[ATLS_S31] * y
    + 2.0 * TLS[ATLS_S21] * z;

  U[U22] =
            TLS[ATLS_T22]
    +       TLS[ATLS_L11] * zz
    +       TLS[ATLS_L33] * xx
    - 2.0 * TLS[ATLS_L13] * xz
    - 2.0 * TLS[ATLS_S12] * z
    + 2.0 * TLS[ATLS_S32] * x;

  U[U33] =
            TLS[ATLS_T33]
    +       TLS[ATLS_L11] * yy
    +       TLS[ATLS_L22] * xx
    - 2.0 * TLS[ATLS_L12] * xy
    - 2.0 * TLS[ATLS_S23] * x
    + 2.0 * TLS[ATLS_S13] * y;

  U[U12] =
            TLS[ATLS_T12]
    -       TLS[ATLS_L33]   * xy
    +       TLS[ATLS_L23]   * xz
    +       TLS[ATLS_L13]   * yz
    -       TLS[ATLS_L12]   * zz
    +       TLS[ATLS_S2211] * z
    +       TLS[ATLS_S31]   * x
    -       TLS[ATLS_S32]   * y;

  U[U13] =
            TLS[ATLS_T13]
    -       TLS[ATLS_L22]   * xz
    +       TLS[ATLS_L23]   * xy
    -       TLS[ATLS_L13]   * yy
    +       TLS[ATLS_L12]   * yz
    +       TLS[ATLS_S1133] * y
    +       TLS[ATLS_S23]   * z
    -       TLS[ATLS_S21]   * x;

  U[U23] =
            TLS[ATLS_T23]
    -       TLS[ATLS_L11]   * yz
    -       TLS[ATLS_L23]   * xx
    +       TLS[ATLS_L13]   * xy
    +       TLS[ATLS_L12]   * xz
    -      (TLS[ATLS_S2211] + TLS[ATLS_S1133]) * x
    +       TLS[ATLS_S12]   * y
    -       TLS[ATLS_S13]   * z;
}

/* anisotropic TLS model */
#define NAME_LEN     8
#define FRAG_ID_LEN  8

struct TLSAtom {
  char    name[NAME_LEN];
  char    frag_id[FRAG_ID_LEN];
  double  x;
  double  y;
  double  z;
  double  u_iso;
  double  U[6];
};

struct TISO_SegmentFitData {
  int                istart;
  int                iend;
  double             parameters[ATLS_NUM_PARAMS];
  struct TLSAtom    *atoms;
};


/* XXX: global pointer to current minimization problem */
struct TISO_SegmentFitData *g_pFit = NULL;

/* prototype for MINPACK FORTRAN subroutine */
typedef void (*FCN)();

extern void lmdif1_(
    FCN, int*, int*, double*, double*, double*, int*, int*, double*, int*);

extern void lmder1_(
    FCN, int*, int*, double*, double*, double*, int*, int*, double*, int*);


double calc_lsqr(int m, double *fvec)
{
  int i;
  double lsqr;

  lsqr = 0.0;

  for (i = 0; i < m; i++) {
    lsqr += fvec[i] * fvec[i];
  }

  return lsqr;
}

void lmdir1_fcn(int *m, int *n,double *x, double *fvec,
		double **fjac, int *ldfjac, int *iflag)
{
  if (*iflag == 1) {
    /* calculate f(x) = fvec */
  } else {
    /* calculate J(x) = fjac */
  }
}


void
lmdif1_fcn(int *m, int *n, double *x, double *fvec, int *iflag)
{
  int i, j;
  double dp2;
  double U[6];
  struct TISO_SegmentFitData *pFit;

  /* XXX: not thread safe, get context */
  pFit = g_pFit;

  for (i = pFit->istart, j = 0; i <= pFit->iend; i++, j++) {

    calc_u_tls(x, pFit->atoms[i].x, pFit->atoms[i].y, pFit->atoms[i].z, U);

    if (!calc_dp2(U, pFit->atoms[i].U, &dp2)) {
      *iflag = -1;
      return;
    }

    fvec[j] = dp2;
  }
}


/* Python interface */
typedef struct {
  PyObject_HEAD
  PyObject       *xmlrpc_chain;
  struct TLSAtom *atoms;
} TLS_ISO_Solver_Object;

static void
TLS_ISO_Solver_dealloc(TLS_ISO_Solver_Object* self)
{
  if (self->atoms) {
    free(self->atoms);
    self->atoms = NULL;
  }

  Py_XDECREF(self->xmlrpc_chain);
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
TLS_ISO_Solver_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  TLS_ISO_Solver_Object *self;
  
  self = (TLS_ISO_Solver_Object *)type->tp_alloc(type, 0);
  if (self == NULL) {
    return NULL;
  }

  self->xmlrpc_chain = NULL;
  self->atoms = NULL;

  return (PyObject *)self;
}

static PyObject *
TLS_ISO_Solver_set_xmlrpc_chain(PyObject *py_self, PyObject *args)
{
  TLS_ISO_Solver_Object *self;
  PyObject *xmlrpc_chain;
  PyObject *atm_desc;
  PyObject *tmp;

  int i, j, num_atoms;
  char *strx;
  double xx;

  self = (TLS_ISO_Solver_Object *) py_self;

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
    self->atoms = (struct TLSAtom *)malloc(sizeof(struct TLSAtom) * num_atoms);

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

#ifdef _DEBUG
      printf("ATOMS[%d]: %s %s\n", 
	     i, 
	     self->atoms[i].name, 
	     self->atoms[i].frag_id);
#endif /* _DEBUG */

    }
  }

  Py_INCREF(Py_None);
  return Py_None;

 error:
  return NULL;
}

static PyObject *
TLS_ISO_Solver_clear_xmlrpc_chain(PyObject *py_self, PyObject *args)
{
  TLS_ISO_Solver_Object *self;
  PyObject *tmp;

  self = (TLS_ISO_Solver_Object *) py_self;

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
TLS_ISO_Solver_fit_segment(PyObject *py_self, PyObject *args)
{
  TLS_ISO_Solver_Object *self;
  char *frag_id1;
  char *frag_id2;
  struct TISO_SegmentFitData fit;

  int i, n, m, *iwa, lwa, info;
  double tol, *fvec, *wa, *params;

  double mean_u_iso;

  PyObject *py_floatx, *py_intx, *rdict;


  self = (TLS_ISO_Solver_Object *) py_self;

  if (!PyArg_ParseTuple(args, "ii", &fit.istart, &fit.iend)) {
    goto error;
  }
  if (self->atoms == NULL) {
    goto error;
  }

  fit.atoms = self->atoms;

  /* number of parameters */
  n = ATLS_NUM_PARAMS; 
  /* number of equations */
  m = fit.iend - fit.istart + 1;

  if (n > m) {
    goto error;
  }

  /* allocate the memory LMDIF1 needs for working space */
  lwa    = (m * n) + 5*n + m;

  params = (double *) malloc(sizeof(double) * n);
  fvec   = (double *) malloc(sizeof(double) * m);
  iwa    = (int *)    malloc(sizeof(int)    * n);
  wa     = (double *) malloc(sizeof(double) * lwa);

  /* initalize the parameters with a initial estimate */
  for (i = 0; i < ATLS_NUM_PARAMS; i++) {
    params[i] = 0.0;
  }

  mean_u_iso = 0.0;
  for (i = fit.istart; i <= fit.iend; i++) {
    mean_u_iso += fit.atoms[i].u_iso;
  }
  mean_u_iso = mean_u_iso / m;

  params[ATLS_T11] = mean_u_iso;
  params[ATLS_T22] = mean_u_iso;
  params[ATLS_T33] = mean_u_iso;

  /* permorm minimization */
  g_pFit = &fit;
  tol = 1e-6;
  lmdif1_(lmdif1_fcn, &m, &n, params, fvec, &tol, &info, iwa, wa, &lwa);

  /* construct return dictioary with results */
  rdict = PyDict_New();

  /* set minimization exit status */
  py_intx = PyInt_FromLong(info);
  PyDict_SetItemString(rdict, "info", py_intx);
  Py_DECREF(py_intx);
  
  for (i = 0; i < ATLS_NUM_PARAMS; i++) {
    py_floatx = PyFloat_FromDouble(params[i]);
    PyDict_SetItemString(rdict, ATLS_PARAM_NAMES[i], py_floatx);
    Py_DECREF(py_floatx);
  }
    
  /* free working memory */
  free(params);
  free(fvec);
  free(iwa);
  free(wa);

  return rdict;

 error:
  return NULL;
}

static PyMethodDef TLS_ISO_Solver_methods[] = {
    {"set_xmlrpc_chain", 
     (PyCFunction) TLS_ISO_Solver_set_xmlrpc_chain, 
     METH_VARARGS,
     "Sets the Python list containing one dictionary for each atom." },

    {"clear_xmlrpc_chain", 
     (PyCFunction) TLS_ISO_Solver_clear_xmlrpc_chain, 
     METH_VARARGS,
     "Clears the Python list of atom descriptions." },

    {"fit_segment", 
     (PyCFunction) TLS_ISO_Solver_fit_segment, 
     METH_VARARGS,
     "Performs a TLS/ISO fit to the given atoms." },

    {NULL}  /* Sentinel */
};

static PyTypeObject TLS_ISO_Solver_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "TLS_ISO_Solver",          /*tp_name*/
    sizeof(TLS_ISO_Solver_Object), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)TLS_ISO_Solver_dealloc, /*tp_dealloc*/
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
    "TLS_ISO_Solver objects",  /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    TLS_ISO_Solver_methods,    /* tp_methods */
    0,                         /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    0,                         /* tp_init */
    0,                         /* tp_alloc */
    TLS_ISO_Solver_new,        /* tp_new */
};


static PyMethodDef TLSMDMODULE_METHODS[] = {
  {NULL, NULL, 0, NULL}
};

DL_EXPORT(void)
inittlsmdmodule(void)
{
  PyObject *m;
  
  
  if (PyType_Ready(&TLS_ISO_Solver_Type) < 0)
    return;

  m = Py_InitModule("tlsmdmodule", TLSMDMODULE_METHODS);
  
  TLSMDMODULE_ERROR = PyErr_NewException("tlsmdmodule.error", NULL, NULL);
  Py_INCREF(TLSMDMODULE_ERROR);
  PyModule_AddObject(m, "error", TLSMDMODULE_ERROR);


  /* add the TLS_ISO_Solver class */
  Py_INCREF(&TLS_ISO_Solver_Type);
  PyModule_AddObject(m, "TLS_ISO_Solver", (PyObject *)&TLS_ISO_Solver_Type);
}
