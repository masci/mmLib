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

#include "tls_model.h"
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

char *NL_ITLS_PARAM_NAMES[] = {
  "nl_t",
  "nl_lx", "nl_ly", "nl_lz", "nl_la", "nl_lb", "nl_lc",
  "nl_s1", "nl_s2", "nl_s3"
};

char *NL_ATLS_PARAM_NAMES[] = {
  "nl_t11", "nl_t22", "nl_t33", "nl_t12", "nl_t13", "nl_t23",
  "nl_lx", "nl_ly", "nl_lz", "nl_la", "nl_lb", "nl_lc",
  "nl_s2211", "nl_s1133", "nl_s12", "nl_s13",  "nl_s23", "nl_s21", "nl_s31", "nl_s32"
};

static void
AddTLSModelToPyDict(const TLSMD::TLSModel &tls_model, PyObject *rdict) {
  PyObject* py_floatx = PyFloat_FromDouble(tls_model.origin_x);
  PyDict_SetItemString(rdict, "x", py_floatx);
  Py_DECREF(py_floatx);

  py_floatx = PyFloat_FromDouble(tls_model.origin_y);
  PyDict_SetItemString(rdict, "y", py_floatx);
  Py_DECREF(py_floatx);

  py_floatx = PyFloat_FromDouble(tls_model.origin_z);
  PyDict_SetItemString(rdict, "z", py_floatx);
  Py_DECREF(py_floatx);
}

static void
AddFitTLSModelResultToPyDict(const TLSMD::FitTLSModelResult &tls_result,  PyObject *rdict) {
  PyObject* py_intx = PyInt_FromLong(tls_result.get_num_atoms());
  PyDict_SetItemString(rdict, "num_atoms", py_intx);
  Py_DECREF(py_intx);

  py_intx = PyInt_FromLong(tls_result.get_num_residues());
  PyDict_SetItemString(rdict, "num_residues", py_intx);
  Py_DECREF(py_intx);

  PyObject* py_floatx = PyFloat_FromDouble(tls_result.get_residual());
  PyDict_SetItemString(rdict, "residual", py_floatx);
  Py_DECREF(py_floatx);
}

static void
AddIsotropicTLSModelToPyDict(const TLSMD::IsotropicTLSModel &itls_model, PyObject *rdict) {
  AddTLSModelToPyDict(itls_model, rdict);

  int num_params = itls_model.num_params();
  char** param_name = ITLS_PARAM_NAMES;
  const double* param = itls_model.get_params();
  for (int i = 0; i < num_params; ++i, ++param, ++param_name) {
    PyObject* py_floatx = PyFloat_FromDouble(*param);
    PyDict_SetItemString(rdict, *param_name, py_floatx);
    Py_DECREF(py_floatx);
  }
}

static void
AddAnisotropicTLSModelToPyDict(const TLSMD::AnisotropicTLSModel &atls_model, PyObject *rdict) {
  AddTLSModelToPyDict(atls_model, rdict);

  int num_params = atls_model.num_params();
  char** param_name = ATLS_PARAM_NAMES;
  const double* param = atls_model.get_params();
  for (int i = 0; i < num_params; ++i, ++param, ++param_name) {
    PyObject* py_floatx = PyFloat_FromDouble(*param);
    PyDict_SetItemString(rdict, *param_name, py_floatx);
    Py_DECREF(py_floatx);
  }
}

static PyObject*
IsotropicFitTLSModelResultToPyDict(const TLSMD::IsotropicFitTLSModelResult& itls_result) {
  PyObject* rdict = PyDict_New();
  AddIsotropicTLSModelToPyDict(itls_result.itls_model, rdict);
  AddFitTLSModelResultToPyDict(itls_result, rdict);
  return rdict;
}

static PyObject*
AnisotropicFitTLSModelResultToPyDict(const TLSMD::AnisotropicFitTLSModelResult& atls_result) {
  PyObject* rdict = PyDict_New();
  AddAnisotropicTLSModelToPyDict(atls_result.atls_model, rdict);
  AddFitTLSModelResultToPyDict(atls_result, rdict);
  return rdict;
}

//
// TLSModelAnalyzer Class
//
static PyObject *TLSMDMODULE_ERROR = NULL;

typedef struct {
  PyObject_HEAD
  TLSMD::TLSModelEngine *tls_model_engine;
} TLSModelAnalyzer_Object;

static void
TLSModelAnalyzer_dealloc(TLSModelAnalyzer_Object* self) {
  if (self->tls_model_engine) {
    delete self->tls_model_engine;
    self->tls_model_engine = 0;
  }
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
TLSModelAnalyzer_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  TLSModelAnalyzer_Object *self;
  self = (TLSModelAnalyzer_Object *)type->tp_alloc(type, 0);
  if (self == NULL) {
    return NULL;
  }
  self->tls_model_engine = new TLSMD::TLSModelEngine();
  return (PyObject *)self;
}

static PyObject *
TLSModelAnalyzer_set_xmlrpc_chain(PyObject *py_self, PyObject *args)
{
  TLSModelAnalyzer_Object *self;
  self = (TLSModelAnalyzer_Object *) py_self;

  PyObject *xmlrpc_chain;
  if (!PyArg_ParseTuple(args, "O", &xmlrpc_chain)) {
    return NULL;
  }

  /* allocate and fill the new atoms array */
  int num_atoms = PyList_Size(xmlrpc_chain);
  self->tls_model_engine->set_num_atoms(num_atoms);
  
  if (num_atoms > 0) {
    int ia = 0;
    std::vector<TLSMD::Atom> &atoms = self->tls_model_engine->chain.atoms;
    std::vector<TLSMD::Atom>::iterator atom;
    for (atom = atoms.begin(); atom != atoms.end(); ++atom, ++ia) {
      char *strx;
      PyObject *atm_desc, *tmp;

      atm_desc = PyList_GetItem(xmlrpc_chain, ia);
      
      /* set name */
      tmp = PyDict_GetItemString(atm_desc, "name");
      if (tmp == NULL) {
	PyErr_SetString(TLSMDMODULE_ERROR, "name not in atm_desc");
	return NULL;
      }
      strx = PyString_AsString(tmp);
      if (strx == NULL) {
	return NULL;
      }
      atom->name.assign(strx);
	
      /* set frag_id */
      tmp = PyDict_GetItemString(atm_desc, "frag_id");
      if (tmp == NULL) {
	PyErr_SetString(TLSMDMODULE_ERROR, "frag_id not in atm_desc");
	return NULL;
      }
      strx = PyString_AsString(tmp);
      if (strx == NULL) {
	return NULL;
      }
      atom->frag_id.assign(strx);

      /* set x, y, z coordinates */
      tmp = PyDict_GetItemString(atm_desc, "ifrag");
      if (tmp == NULL) {
	PyErr_SetString(TLSMDMODULE_ERROR, "ifrag not in atm_desc");
	return NULL;
      }
      if (!PyInt_Check(tmp)) {
	return NULL;
      }
      atom->ifrag = PyInt_AsLong(tmp);

      /* set x, y, z coordinates */
      tmp = PyDict_GetItemString(atm_desc, "x");
      if (tmp == NULL) {
	PyErr_SetString(TLSMDMODULE_ERROR, "x not in atm_desc");
	return NULL;
      }
      if (!PyFloat_Check(tmp)) {
	return NULL;
      }
      atom->x = PyFloat_AsDouble(tmp);

      tmp = PyDict_GetItemString(atm_desc, "y");
      if (tmp == NULL) {
	PyErr_SetString(TLSMDMODULE_ERROR, "y not in atm_desc");
	return NULL;
      }
      if (!PyFloat_Check(tmp)) {
	return NULL;
      }
      atom->y = PyFloat_AsDouble(tmp);

      tmp = PyDict_GetItemString(atm_desc, "z");
      if (tmp == NULL) {
	PyErr_SetString(TLSMDMODULE_ERROR, "z not in atm_desc");
	return NULL;
      }
      if (!PyFloat_Check(tmp)) {
	return NULL;
      }
      atom->z = PyFloat_AsDouble(tmp);

      tmp = PyDict_GetItemString(atm_desc, "u_iso");
      if (tmp == NULL) {
	PyErr_SetString(TLSMDMODULE_ERROR, "u_iso not in atm_desc");
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
	  PyErr_SetString(TLSMDMODULE_ERROR, "uXX not in atm_desc");
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
	PyErr_SetString(TLSMDMODULE_ERROR, "weight not in atm_desc");
	return NULL;
      }
      if (!PyFloat_Check(tmp)) {
	return NULL;
      }
      atom->weight = PyFloat_AsDouble(tmp);
      atom->sqrt_weight = sqrt(PyFloat_AsDouble(tmp));
    }
  }

  self->tls_model_engine->chain.map_frag_ids();

  Py_INCREF(Py_None);
  return Py_None;
}

static bool
PythonSegmentListToSegmentSet(PyObject *segment_list, TLSMD::Chain::SegmentSet* segment_set) {
  int num_segments = PyList_Size(segment_list);

  for (int is = 0; is < num_segments; ++is) {
    PyObject* range = PyList_GetItem(segment_list, is);
    
    char *cfrag_id1, *cfrag_id2;
    if (!PyArg_ParseTuple(range, "ss", &cfrag_id1, &cfrag_id2)) {
      return false;
    }
    std::string frag_id1(cfrag_id1);
    std::string frag_id2(cfrag_id2);
    
    try {
	segment_set->add_segment(frag_id1, frag_id2);
    } catch(TLSMD::Chain::FragmentIDMap::FragmentIDNotFound fnf) {
      std::string msg;
      msg = "fragment id not found: " + fnf.frag_id();
      PyErr_SetString(TLSMDMODULE_ERROR, msg.c_str());
      return false;
    }
  }

  return true;
}

static PyObject*
TLSModelAnalyzer_isotropic_fit(PyObject *py_self, PyObject *args) {
  TLSModelAnalyzer_Object *self;
  self = (TLSModelAnalyzer_Object *) py_self;

  PyObject *segment_list;
  if (!PyArg_ParseTuple(args, "O", &segment_list)) return NULL;

  TLSMD::Chain::SegmentSet segment_set(&self->tls_model_engine->chain);
  if (!PythonSegmentListToSegmentSet(segment_list, &segment_set)) return NULL;

  TLSMD::IsotropicFitTLSModelResult itls_result;
  self->tls_model_engine->isotropic_fit(segment_set, itls_result);
  return IsotropicFitTLSModelResultToPyDict(itls_result);
}

static PyObject*
TLSModelAnalyzer_isotropic_fit_segment(PyObject *py_self, PyObject *args) {
  TLSModelAnalyzer_Object *self;
  self = (TLSModelAnalyzer_Object *) py_self;

  char *cfrag_id1, *cfrag_id2;
  if (!PyArg_ParseTuple(args, "ss", &cfrag_id1, &cfrag_id2)) {
    return NULL;
  }
  std::string frag_id1(cfrag_id1);
  std::string frag_id2(cfrag_id2);

  TLSMD::IsotropicFitTLSModelResult itls_result;
  try {
    self->tls_model_engine->isotropic_fit_segment(frag_id1, frag_id2, itls_result);
  } catch(TLSMD::Chain::FragmentIDMap::FragmentIDNotFound fnf) {
    std::string msg;
    msg = "fragment id not found: " + fnf.frag_id();
    PyErr_SetString(TLSMDMODULE_ERROR, msg.c_str());
    return false;
  }
  return IsotropicFitTLSModelResultToPyDict(itls_result);
}

static PyObject*
TLSModelAnalyzer_anisotropic_fit(PyObject *py_self, PyObject *args) {
  TLSModelAnalyzer_Object *self;
  self = (TLSModelAnalyzer_Object *) py_self;

  PyObject *segment_list;
  if (!PyArg_ParseTuple(args, "O", &segment_list)) return NULL;

  TLSMD::Chain::SegmentSet segment_set(&self->tls_model_engine->chain);
  if (!PythonSegmentListToSegmentSet(segment_list, &segment_set)) return NULL;

  TLSMD::AnisotropicFitTLSModelResult atls_result;
  self->tls_model_engine->anisotropic_fit(segment_set, atls_result);
  return AnisotropicFitTLSModelResultToPyDict(atls_result);
}

static PyObject*
TLSModelAnalyzer_anisotropic_fit_segment(PyObject *py_self, PyObject *args) {
  TLSModelAnalyzer_Object *self;
  self = (TLSModelAnalyzer_Object *) py_self;

  char *cfrag_id1, *cfrag_id2;
  if (!PyArg_ParseTuple(args, "ss", &cfrag_id1, &cfrag_id2)) {
    return NULL;
  }
  std::string frag_id1(cfrag_id1);
  std::string frag_id2(cfrag_id2);

  TLSMD::AnisotropicFitTLSModelResult atls_result;
  try {
    self->tls_model_engine->anisotropic_fit_segment(frag_id1, frag_id2, atls_result);
  } catch(TLSMD::Chain::FragmentIDMap::FragmentIDNotFound fnf) {
    std::string msg;
    msg = "fragment id not found: " + fnf.frag_id();
    PyErr_SetString(TLSMDMODULE_ERROR, msg.c_str());
    return false;
  }

  return AnisotropicFitTLSModelResultToPyDict(atls_result);
}

static PyObject*
TLSModelAnalyzer_constrained_isotropic_fit(PyObject *py_self, PyObject *args) {
  TLSModelAnalyzer_Object *self;
  self = (TLSModelAnalyzer_Object *) py_self;

  PyObject *segment_list;
  if (!PyArg_ParseTuple(args, "O", &segment_list)) return NULL;

  TLSMD::Chain::SegmentSet segment_set(&self->tls_model_engine->chain);
  if (!PythonSegmentListToSegmentSet(segment_list, &segment_set)) return NULL;

  TLSMD::IsotropicFitTLSModelResult itls_result;
  self->tls_model_engine->constrained_isotropic_fit(segment_set, itls_result);
  return IsotropicFitTLSModelResultToPyDict(itls_result);
}

static PyObject*
TLSModelAnalyzer_constrained_isotropic_fit_segment(PyObject *py_self, PyObject *args) {
  TLSModelAnalyzer_Object *self;
  self = (TLSModelAnalyzer_Object *) py_self;

  char *cfrag_id1, *cfrag_id2;
  if (!PyArg_ParseTuple(args, "ss", &cfrag_id1, &cfrag_id2)) {
    return NULL;
  }
  std::string frag_id1(cfrag_id1);
  std::string frag_id2(cfrag_id2);

  TLSMD::IsotropicFitTLSModelResult itls_result;
  try {
    self->tls_model_engine->constrained_isotropic_fit_segment(frag_id1, frag_id2, itls_result);
  } catch(TLSMD::Chain::FragmentIDMap::FragmentIDNotFound fnf) {
    std::string msg;
    msg = "fragment id not found: " + fnf.frag_id();
    PyErr_SetString(TLSMDMODULE_ERROR, msg.c_str());
    return false;
  }
  return IsotropicFitTLSModelResultToPyDict(itls_result);
}

static PyObject*
TLSModelAnalyzer_constrained_anisotropic_fit(PyObject *py_self, PyObject *args)
{
  TLSModelAnalyzer_Object *self;
  self = (TLSModelAnalyzer_Object *) py_self;

  PyObject *segment_list;
  if (!PyArg_ParseTuple(args, "O", &segment_list)) return NULL;

  TLSMD::Chain::SegmentSet segment_set(&self->tls_model_engine->chain);
  if (!PythonSegmentListToSegmentSet(segment_list, &segment_set)) return NULL;

  TLSMD::AnisotropicFitTLSModelResult atls_result;
  self->tls_model_engine->constrained_anisotropic_fit(segment_set, atls_result);
  return AnisotropicFitTLSModelResultToPyDict(atls_result);
}

static PyObject*
TLSModelAnalyzer_constrained_anisotropic_fit_segment(PyObject *py_self, PyObject *args) {
  TLSModelAnalyzer_Object *self;
  self = (TLSModelAnalyzer_Object *) py_self;

  char *cfrag_id1, *cfrag_id2;
  if (!PyArg_ParseTuple(args, "ss", &cfrag_id1, &cfrag_id2)) {
    return NULL;
  }
  std::string frag_id1(cfrag_id1);
  std::string frag_id2(cfrag_id2);

  TLSMD::AnisotropicFitTLSModelResult atls_result;
  try {
    self->tls_model_engine->constrained_anisotropic_fit_segment(frag_id1, frag_id2, atls_result);
  } catch(TLSMD::Chain::FragmentIDMap::FragmentIDNotFound fnf) {
    std::string msg;
    msg = "fragment id not found: " + fnf.frag_id();
    PyErr_SetString(TLSMDMODULE_ERROR, msg.c_str());
    return false;
  }
  return AnisotropicFitTLSModelResultToPyDict(atls_result);
}

static PyMethodDef TLSModelAnalyzer_methods[] = {
    {"set_xmlrpc_chain", 
     (PyCFunction) TLSModelAnalyzer_set_xmlrpc_chain, 
     METH_VARARGS,
     "Sets the Python list containing one dictionary for each atom." },

    {"isotropic_fit",
     (PyCFunction) TLSModelAnalyzer_isotropic_fit, 
     METH_VARARGS,
     "Performs a linear fit of the isotropic TLS model to a list of 2-tuple fragment id ranges." },

    {"isotropic_fit_segment",
     (PyCFunction) TLSModelAnalyzer_isotropic_fit_segment, 
     METH_VARARGS,
     "Performs a linear fit of the isotropic TLS model to the given atoms." },

    {"anisotropic_fit",
     (PyCFunction) TLSModelAnalyzer_anisotropic_fit, 
     METH_VARARGS,
     "Performs a linear fit of the anisotropic TLS model to a list of 2-tuple fragment id ranges." },

    {"anisotropic_fit_segment",
     (PyCFunction) TLSModelAnalyzer_anisotropic_fit_segment, 
     METH_VARARGS,
     "Performs a linear fit of the anisotropic TLS model to the given atoms." },

    {"constrained_isotropic_fit",
     (PyCFunction) TLSModelAnalyzer_constrained_isotropic_fit, 
     METH_VARARGS,
     "Performs a constrained fit of the isotropic TLS model to a list of 2-tuple fragment id ranges." },

    {"constrained_isotropic_fit_segment",
     (PyCFunction) TLSModelAnalyzer_constrained_isotropic_fit_segment, 
     METH_VARARGS,
     "Performs a constrained fit of the isotropic TLS model to the given atoms." },

    {"constrained_anisotropic_fit",
     (PyCFunction) TLSModelAnalyzer_constrained_anisotropic_fit, 
     METH_VARARGS,
     "Performs a constrained fit of the anisotropic TLS model to a list of 2-tuple fragment id ranges." },

    {"constrained_anisotropic_fit_segment",
     (PyCFunction) TLSModelAnalyzer_constrained_anisotropic_fit_segment, 
     METH_VARARGS,
     "Performs a constrained fit of the anisotropic TLS model to the given atoms." },

    {NULL}  /* Sentinel */
};

static PyTypeObject TLSModelAnalyzer_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "TLSModelAnalyzer",          /*tp_name*/
    sizeof(TLSModelAnalyzer_Object), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)TLSModelAnalyzer_dealloc, /*tp_dealloc*/
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
    "TLSModelAnalyzer objects",  /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    TLSModelAnalyzer_methods,    /* tp_methods */
    0,                         /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    0,                         /* tp_init */
    0,                         /* tp_alloc */
    TLSModelAnalyzer_new,        /* tp_new */
};

static PyMethodDef TLSMDMODULE_METHODS[] = {
  {NULL, NULL, 0, NULL}
};

extern "C" DL_EXPORT(void)
inittlsmdmodule(void)
{
  if (PyType_Ready(&TLSModelAnalyzer_Type) < 0)
    return;

  PyObject *m;
  m = Py_InitModule("tlsmdmodule", TLSMDMODULE_METHODS);
  
  TLSMDMODULE_ERROR = PyErr_NewException("tlsmdmodule.error", NULL, NULL);
  Py_INCREF(TLSMDMODULE_ERROR);
  PyModule_AddObject(m, "error", TLSMDMODULE_ERROR);


  /* add the TLSModelAnalyzer class */
  Py_INCREF(&TLSModelAnalyzer_Type);
  PyModule_AddObject(m, "TLSModelAnalyzer", (PyObject *)&TLSModelAnalyzer_Type);
}
