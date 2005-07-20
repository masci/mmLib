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
#define ITLS_TISO  0
#define ITLS_L11   1
#define ITLS_L22   2
#define ITLS_L33   3
#define ITLS_L12   4
#define ITLS_L13   5
#define ITLS_L23   6
#define ITLS_S12   7
#define ITLS_S13   8
#define ITLS_S23   9
#define ITLS_S21   10
#define ITLS_S31   11
#define ITLS_S32   12




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


/* calculate the anistropic TLS model prediction of anisiotropic
 * ADP U at position x, y, z
 */
static void
calc_u_atls(double ATLS[ATLS_NUM_PARAMS], 
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



/*
 * Common routines for the isotropic and ansiotropic TLS segment
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
};

/* structure used to store information on the chain the algorithm
 * is currently fitting; this structure contains the dynamically 
 * allocated memory used by the LAPACK SVD routine DGESDD and therefore
 * must be allocated and freed with new_chain/delete_chain
 */
struct Chain {
  struct Atom *atoms;
  
  double *A;
  double *x;
  double *b;
  double *S;
  double *U;
  double *VT;
  double *WORK;
  int     LWORK;
};

struct Chain *
new_chain(int num_atoms)
{
  struct Chain *chain;

  chain = malloc(sizeof(Chain));
  if (chain==NULL)
    return NULL;
  bzero(chain, sizeof(Chain));

  chain->atoms = malloc(sizeof(Atom) * num_atoms);
  if (chain->atoms==NULL)
    goto error;
    
  chain->A = malloc(sizeof(double) * 20 * num_atoms);
  if (chain->A==NULL)
    goto error;

  chain->x = malloc(sizeof(double) * 20);
  if (chain->x==NULL)
    goto error;

  chain->b = malloc(sizeof(double) * num_atoms);
  if (chain->b==NULL)
    goto error;
  
  chain->S = malloc(sizeof(double) * 20);
  if (chain->S==NULL)
    goto error;

  chain->U = malloc(sizeof(double) * num_atoms * num_atoms);
  if (chain->U==NULL)
    goto error;

  chain->VT = malloc(sizeof(double) * 20 * 20);
  if (chain->VT==NULL)
    goto error;

  /* calculate the size of the WORK memory block */


  return chain;

  /* if malloc() fails */
 error:
  if (chain->atoms!=NULL)
    free(chain->atoms);
  if (chain->A!=NULL)
    free(chain->A);
  if (chain->x!=NULL)
    free(chain->x);
  if (chain->b!=NULL)
    free(chain->b);
  if (chain->S!=NULL)
    free(chain->S);
  if (chain->U!=NULL)
    free(chain->U);
  if (chain->VT!=NULL)
    free(chain->VT);
  if (chain->WORK!=NULL)
    free(chain->WORK);

  free(chain);
}

void
delete_chain(struct Chain *chain)
{
  if (chain->atoms!=NULL)
    free(chain->atoms);
  if (chain->A!=NULL)
    free(chain->A);
  if (chain->x!=NULL)
    free(chain->x);
  if (chain->b!=NULL)
    free(chain->b);
  if (chain->S!=NULL)
    free(chain->S);
  if (chain->U!=NULL)
    free(chain->U);
  if (chain->VT!=NULL)
    free(chain->VT);
  if (chain->WORK!=NULL)
    free(chain->WORK);

  free(chain);
}



/* context structure for fitting one TLS segment using the 
 * isotropic TLS model
 */
struct ITLSFitContext {
  int                istart;                /* index of first atom in atoms */
  int                iend;                  /* index of last atom in atoms */
  struct Atom       *atoms;                 /* pointer to atoms array */

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
itls_fit_segment(struct ITLSFitContext *itls_context)
{
  int num_atoms;
  double i, ia, x, y, z;

  /* calculate the number of atoms to be fit */
  num_atoms = itls_context->iend - itls_context->istart + 1;

  /* calculate the centroid of the atoms which are to
   * be fit and use it as the origin of the TLS tensors
   */
  calc_centroid(itls_context->atoms,
		itls_context->istart,
		itls_context->iend,
		&itls_context->origin_x,
		&itls_context->origin_y,
		&itls_context->origin_z);


  /* fill in coefficent matrix for both the isotropic TLS model
   * and the anisotropic TLS model
   */



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
  if (self->atoms) {
    free(self->atoms);
    self->atoms = NULL;
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
ITLSModel_clear_xmlrpc_chain(PyObject *py_self, PyObject *args)
{
  ITLSModel_Object *self;
  PyObject *tmp;

  self = (ITLSModel_Object *) py_self;

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
ITLSModel_fit_segment(PyObject *py_self, PyObject *args)
{
  ITLSModel_Object *self;
  char *frag_id1;
  char *frag_id2;
  struct TISO_SegmentFitData fit;

  int i, n, m, *iwa, lwa, info;
  double tol, *fvec, *wa, *params;

  double mean_u_iso;

  PyObject *py_floatx, *py_intx, *rdict;


  self = (ITLSModel_Object *) py_self;

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

static PyMethodDef ITLSModel_methods[] = {
    {"set_xmlrpc_chain", 
     (PyCFunction) ITLSModel_set_xmlrpc_chain, 
     METH_VARARGS,
     "Sets the Python list containing one dictionary for each atom." },

    {"clear_xmlrpc_chain", 
     (PyCFunction) ITLSModel_clear_xmlrpc_chain, 
     METH_VARARGS,
     "Clears the Python list of atom descriptions." },

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
