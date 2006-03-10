// tlsmdmodule.cpp
// Jay Painter <jpaint@u.washington.edu>
// March 10, 2006
// 
// Copyright 2006 by TLSMD Development Group (see AUTHORS file)
// This code is part of the TLSMD distribution and governed by
// its license.  Please see the LICENSE file that should have been
// included as part of this package.
//
// Implementation of a Python module for the unconstrained linear 
// fitting of TLS parameters given a list of atoms with 
// crystallographically refined ADPs.  Uses LAPACK.
#include "Python.h"
#include "structmember.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#include "tls_model_engine.h"

// names for return dictionaries
char *U_PARAM_NAMES[] = {
  "u11", "u22", "u33", "u12", "u13", "u23"
};

char *ITLS_PARAM_NAMES[] = {
  "it",
  "il11", "il22", "il33", "il12", "il13", "il23",
  "is1", "is2", "is3"
};

char *ATLS_PARAM_NAMES[] = {
  "t11", "t22", "t33", "t12", "t13", "t23",
  "l11", "l22", "l33", "l12", "l13", "l23",
  "s2211", "s1133", "s12", "s13", "s23", "s21", "s31", "s32"
};

// PYTHON INTERFACE
// LinearTLSModel: Linear Fit of Isotropic and Anisotropic TLS Models

static PyObject *LINEARTLS_ERROR = NULL;


typedef struct {
  PyObject_HEAD
  TLSModelEngine *tls_model_engine;
} LinearTLSModel_Object;


static void
LinearTLSModel_dealloc(LinearTLSModel_Object* self) {
  if (self->tls_model_engine) {
    delete self->tls_model_engine;
    self->tls_model_engine = 0;
  }
  self->ob_type->tp_free((PyObject*)self);
}


static PyObject *
LinearTLSModel_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  LinearTLSModel_Object *self;
  self = (LinearTLSModel_Object *)type->tp_alloc(type, 0);
  if (self == NULL) {
    return NULL;
  }
  self->tls_model_engine = new TLSModelEngine();
  return (PyObject *)self;
}


static PyObject *
LinearTLSModel_set_xmlrpc_chain(PyObject *py_self, PyObject *args)
{
  LinearTLSModel_Object *self;
  self = (LinearTLSModel_Object *) py_self;

  PyObject *xmlrpc_chain;
  int num_atoms;

  if (!PyArg_ParseTuple(args, "O", &xmlrpc_chain)) {
    return NULL;
  }

  /* allocate and fill the new atoms array */
  num_atoms = PyList_Size(xmlrpc_chain);
  self->tls_model_engine->set_num_atoms(num_atoms);
  
  if (num_atoms > 0) {
    Atom *atom = self->tls_model_engine->chain.atoms;
    for (int i = 0; i < num_atoms; ++i, ++atom) {
      char *strx;
      PyObject *atm_desc, *tmp;

      atm_desc = PyList_GetItem(xmlrpc_chain, i);
      
      /* set name */
      tmp = PyDict_GetItemString(atm_desc, "name");
      if (tmp == NULL) {
	PyErr_SetString(LINEARTLS_ERROR, "name not in atm_desc");
	return NULL;
      }
      strx = PyString_AsString(tmp);
      if (strx == NULL) {
	return NULL;
      }
      strncpy(atom->name, strx, NAME_LEN);
	
      /* set frag_id */
      tmp = PyDict_GetItemString(atm_desc, "frag_id");
      if (tmp == NULL) {
	PyErr_SetString(LINEARTLS_ERROR, "frag_id not in atm_desc");
	return NULL;
      }
      strx = PyString_AsString(tmp);
      if (strx == NULL) {
	return NULL;
      }
      strncpy(atom->frag_id, strx, FRAG_ID_LEN);

      /* set x, y, z coordinates */
      tmp = PyDict_GetItemString(atm_desc, "ifrag");
      if (tmp == NULL) {
	PyErr_SetString(LINEARTLS_ERROR, "ifrag not in atm_desc");
	return NULL;
      }
      if (!PyInt_Check(tmp)) {
	return NULL;
      }
      atom->ifrag = PyInt_AsLong(tmp);

      /* set x, y, z coordinates */
      tmp = PyDict_GetItemString(atm_desc, "x");
      if (tmp == NULL) {
	PyErr_SetString(LINEARTLS_ERROR, "x not in atm_desc");
	return NULL;
      }
      if (!PyFloat_Check(tmp)) {
	return NULL;
      }
      atom->x = PyFloat_AsDouble(tmp);

      tmp = PyDict_GetItemString(atm_desc, "y");
      if (tmp == NULL) {
	PyErr_SetString(LINEARTLS_ERROR, "y not in atm_desc");
	return NULL;
      }
      if (!PyFloat_Check(tmp)) {
	return NULL;
      }
      atom->y = PyFloat_AsDouble(tmp);

      tmp = PyDict_GetItemString(atm_desc, "z");
      if (tmp == NULL) {
	PyErr_SetString(LINEARTLS_ERROR, "z not in atm_desc");
	return NULL;
      }
      if (!PyFloat_Check(tmp)) {
	return NULL;
      }
      atom->z = PyFloat_AsDouble(tmp);

      tmp = PyDict_GetItemString(atm_desc, "u_iso");
      if (tmp == NULL) {
	PyErr_SetString(LINEARTLS_ERROR, "u_iso not in atm_desc");
	return NULL;
      }
      if (!PyFloat_Check(tmp)) {
	return NULL;
      }
      atom->u_iso = PyFloat_AsDouble(tmp);

      /* get U tensor parameters */
      for (int j = 0; j < U_NUM_PARAMS; j++) {
	tmp = PyDict_GetItemString(atm_desc, U_PARAM_NAMES[j]);
	if (tmp == NULL) {
	  PyErr_SetString(LINEARTLS_ERROR, "uXX not in atm_desc");
	  return NULL;
	}
	if (!PyFloat_Check(tmp)) {
	  return NULL;
	}
	atom->U[j] = PyFloat_AsDouble(tmp);
      }

      /* weight */
      tmp = PyDict_GetItemString(atm_desc, "weight");
      if (tmp == NULL) {
	PyErr_SetString(LINEARTLS_ERROR, "weight not in atm_desc");
	return NULL;
      }
      if (!PyFloat_Check(tmp)) {
	return NULL;
      }
      atom->sqrt_weight = sqrt(PyFloat_AsDouble(tmp));
    }
  }

  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject *
LinearTLSModel_isotropic_fit_segment(PyObject *py_self, PyObject *args)
{
  LinearTLSModel_Object *self;
  self = (LinearTLSModel_Object *) py_self;

  int istart, iend;
  if (!PyArg_ParseTuple(args, "ii", &istart, &iend)) {
    return NULL;
  }

  double residual;
  self->tls_model_engine->isotropic_fit_segment(istart, iend, &residual);

  /* construct return dictioary with results */
  PyObject *rdict = PyDict_New();

  /* set minimization exit status */
  PyObject *py_floatx;

  py_floatx = PyFloat_FromDouble(self->tls_model_engine->itls.origin_x);
  PyDict_SetItemString(rdict, "x", py_floatx);
  Py_DECREF(py_floatx);

  py_floatx = PyFloat_FromDouble(self->tls_model_engine->itls.origin_y);
  PyDict_SetItemString(rdict, "y", py_floatx);
  Py_DECREF(py_floatx);

  py_floatx = PyFloat_FromDouble(self->tls_model_engine->itls.origin_z);
  PyDict_SetItemString(rdict, "z", py_floatx);
  Py_DECREF(py_floatx);

  py_floatx = PyFloat_FromDouble(residual);
  PyDict_SetItemString(rdict, "ilsqr", py_floatx);
  Py_DECREF(py_floatx);

  int num_params = ITLS_NUM_PARAMS;
  char **param_name = ITLS_PARAM_NAMES;
  double *param = self->tls_model_engine->itls.ITLS;
  for (int i = 0; i < num_params; ++i, ++param, ++param_name) {
    py_floatx = PyFloat_FromDouble(*param);
    PyDict_SetItemString(rdict, *param_name, py_floatx);
    Py_DECREF(py_floatx);
  }
    
  return rdict;
}


static PyObject *
LinearTLSModel_anisotropic_fit_segment(PyObject *py_self, PyObject *args)
{
  LinearTLSModel_Object *self;
  self = (LinearTLSModel_Object *) py_self;

  int istart, iend;
  if (!PyArg_ParseTuple(args, "ii", &istart, &iend)) {
    return NULL;
  }

  double residual;
  self->tls_model_engine->anisotropic_fit_segment(istart, iend, &residual);

  /* construct return dictioary with results */
  PyObject *rdict = PyDict_New();

  /* set minimization exit status */
  PyObject *py_floatx;

  py_floatx = PyFloat_FromDouble(self->tls_model_engine->atls.origin_x);
  PyDict_SetItemString(rdict, "x", py_floatx);
  Py_DECREF(py_floatx);

  py_floatx = PyFloat_FromDouble(self->tls_model_engine->atls.origin_y);
  PyDict_SetItemString(rdict, "y", py_floatx);
  Py_DECREF(py_floatx);

  py_floatx = PyFloat_FromDouble(self->tls_model_engine->atls.origin_z);
  PyDict_SetItemString(rdict, "z", py_floatx);
  Py_DECREF(py_floatx);

  py_floatx = PyFloat_FromDouble(residual);
  PyDict_SetItemString(rdict, "alsqr", py_floatx);
  Py_DECREF(py_floatx);

  int num_params = ATLS_NUM_PARAMS;
  char **param_name = ATLS_PARAM_NAMES;
  double *param = self->tls_model_engine->atls.ATLS;
  for (int i = 0; i < num_params; ++i, ++param, ++param_name) {
    py_floatx = PyFloat_FromDouble(*param);
    PyDict_SetItemString(rdict, *param_name, py_floatx);
    Py_DECREF(py_floatx);
  }
    
  return rdict;
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


extern "C" DL_EXPORT(void)
initlineartls(void)
{
  if (PyType_Ready(&LinearTLSModel_Type) < 0)
    return;

  PyObject *m;
  m = Py_InitModule("lineartls", LINEARTLS_METHODS);
  
  LINEARTLS_ERROR = PyErr_NewException("lineartls.error", NULL, NULL);
  Py_INCREF(LINEARTLS_ERROR);
  PyModule_AddObject(m, "error", LINEARTLS_ERROR);


  /* add the LinearTLSModel class */
  Py_INCREF(&LinearTLSModel_Type);
  PyModule_AddObject(m, "LinearTLSModel", (PyObject *)&LinearTLSModel_Type);
}
