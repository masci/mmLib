/* pdbmodule.c - PDB parser/accelorator for mmLib
 *
 */
#include "Python.h"
#include <stdio.h>
#include <math.h>
#include <GL/glut.h>

static PyObject *GLAccelError = NULL;

void 
normalize(float v[3])
{
  float d;

  d    = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  v[0] = v[0] / d;
  v[1] = v[1] / d;
  v[2] = v[2] / d;
}

void 
cross(float u[3], float v[3], float cp[3])
{
  cp[0] = u[1]*v[2] - u[2]*v[1];
  cp[1] = u[2]*v[0] - u[0]*v[2];
  cp[2] = u[0]*v[1] - u[1]*v[0];
}


/* functions for the rendering of atomic thermal peanuts
 */

void
peanut_func(float U[6], float v[3], float w[3]) 
{
  float d;

  d =     U[0]*v[0]*v[0] +     U[1]*v[1]*v[1] + U[2]*v[2]*v[2]
    + 2.0*U[3]*v[0]*v[1] + 2.0*U[4]*v[0]*v[2] + 2.0*U[5]*v[1]*v[2];

  /* abort before we take the sqare root of a negitive
   * number */
  if (d < 0.0) {
    w[0] = 0.0;
    w[1] = 0.0;
    w[2] = 0.0;

    return;
  }

  d = sqrt(d);

  w[0] = d * v[0];
  w[1] = d * v[1];
  w[2] = d * v[2];
}

void
peanut_normal(float U[6], float v[3], float n[3])
{
  float d;

  d =     U[0]*v[0]*v[0] +     U[1]*v[1]*v[1] + U[2]*v[2]*v[2]
    + 2.0*U[3]*v[0]*v[1] + 2.0*U[4]*v[0]*v[2] + 2.0*U[5]*v[1]*v[2];

  /* abort before we take the sqare root of a negitive
   * number */
  if (d < 0.0) {
    n[0] = 0.0;
    n[1] = 0.0;
    n[2] = 0.0;
    
    return;
  }

  d = 1.0 / (2.0 * sqrt(d));

  n[0] = d * (2.0*U[0]*v[0] +     U[3]*v[1] +    U[4]*v[2]);
  n[1] = d * (    U[3]*v[0] + 2.0*U[1]*v[1] +    U[5]*v[2]);
  n[2] = d * (    U[4]*v[0] +     U[5]*v[1] +2.0*U[2]*v[2]);

  normalize(n);
}

void
peanut_triangle(float U[6], float v1[3], float v2[3], float v3[3])
{
  float v[3], n[3];

  peanut_func(U, v1, v);
  peanut_normal(U, v, n);
  glNormal3f(n[0], n[1], n[2]);
  glVertex3f(v[0], v[1], v[2]);

  peanut_func(U, v2, v);
  peanut_normal(U, v, n);
  glNormal3f(n[0], n[1], n[2]);
  glVertex3f(v[0], v[1], v[2]);

  peanut_func(U, v3, v);
  peanut_normal(U, v, n);
  glNormal3f(n[0], n[1], n[2]);
  glVertex3f(v[0], v[1], v[2]);
}

void
peanut_tesselate(float U[6], float v1[3], float v2[3], float v3[3], int depth)
{
  float v12[3], v23[3], v31[3];

  v12[0] = ((v2[0] - v1[0]) / 2.0) + v1[0];
  v12[1] = ((v2[1] - v1[1]) / 2.0) + v1[1];
  v12[2] = ((v2[2] - v1[2]) / 2.0) + v1[2];

  v23[0] = ((v3[0] - v2[0]) / 2.0) + v2[0];
  v23[1] = ((v3[1] - v2[1]) / 2.0) + v2[1];
  v23[2] = ((v3[2] - v2[2]) / 2.0) + v2[2];

  v31[0] = ((v1[0] - v3[0]) / 2.0) + v3[0];
  v31[1] = ((v1[1] - v3[1]) / 2.0) + v3[1];
  v31[2] = ((v1[2] - v3[2]) / 2.0) + v3[2];

  normalize(v12);
  normalize(v23);
  normalize(v31);

  depth -= 1;

  if (depth > 0) {
    peanut_tesselate(U, v1,  v12, v31, depth);
    peanut_tesselate(U, v2,  v23, v12, depth);
    peanut_tesselate(U, v3,  v31, v23, depth);
    peanut_tesselate(U, v12, v23, v31, depth);
  } else {
    peanut_triangle(U, v1,  v12, v31);
    peanut_triangle(U, v2,  v23, v12);
    peanut_triangle(U, v3,  v31, v23);
    peanut_triangle(U, v12, v23, v31);
  }
}

static PyObject *
glaccel_Upeanut(PyObject *self, PyObject *args)
{
  int          depth;
  float        x, y, z;
  float        U[6];
  float        v1[3], v2[3], v3[3];

  if (!PyArg_ParseTuple(args, "fffffffffi", 
			&x, &y, &z,
			&U[0], &U[1], &U[2], &U[3], &U[4], &U[5],
			&depth)) {
    return NULL;
  }

  
  glPushMatrix();
  glTranslatef(x, y, z);
  glEnable(GL_LIGHTING);

  glBegin(GL_TRIANGLES);

  /* x, y, z */
  v1[0] = 1.0;
  v1[1] = 0.0;
  v1[2] = 0.0;
  
  v2[0] = 0.0;
  v2[1] = 1.0;
  v2[2] = 0.0;
  
  v3[0] = 0.0;
  v3[1] = 0.0;
  v3[2] = 1.0;
  peanut_tesselate(U, v1, v2, v3, depth);

  /* x, -z, y */
  v1[0] = 1.0;
  v1[1] = 0.0;
  v1[2] = 0.0;
  
  v2[0] = 0.0;
  v2[1] = 0.0;
  v2[2] =-1.0;
  
  v3[0] = 0.0;
  v3[1] = 1.0;
  v3[2] = 0.0;
  peanut_tesselate(U, v1, v2, v3, depth);

  /* x, z, -y */
  v1[0] = 1.0;
  v1[1] = 0.0;
  v1[2] = 0.0;
  
  v2[0] = 0.0;
  v2[1] = 0.0;
  v2[2] = 1.0;
  
  v3[0] = 0.0;
  v3[1] =-1.0;
  v3[2] = 0.0;
  peanut_tesselate(U, v1, v2, v3, depth);

  /* x, -y, -z */
  v1[0] = 1.0;
  v1[1] = 0.0;
  v1[2] = 0.0;
  
  v2[0] = 0.0;
  v2[1] =-1.0;
  v2[2] = 0.0;
  
  v3[0] = 0.0;
  v3[1] = 0.0;
  v3[2] =-1.0;
  peanut_tesselate(U, v1, v2, v3, depth);

  /* -x, z, y */
  v1[0] =-1.0;
  v1[1] = 0.0;
  v1[2] = 0.0;
  
  v2[0] = 0.0;
  v2[1] = 0.0;
  v2[2] = 1.0;
  
  v3[0] = 0.0;
  v3[1] = 1.0;
  v3[2] = 0.0;
  peanut_tesselate(U, v1, v2, v3, depth);

  /* -x, y, -z */
  v1[0] =-1.0;
  v1[1] = 0.0;
  v1[2] = 0.0;
  
  v2[0] = 0.0;
  v2[1] = 1.0;
  v2[2] = 0.0;
  
  v3[0] = 0.0;
  v3[1] = 0.0;
  v3[2] =-1.0;
  peanut_tesselate(U, v1, v2, v3, depth);

  /* -x, -y, z */
  v1[0] =-1.0;
  v1[1] = 0.0;
  v1[2] = 0.0;
  
  v2[0] = 0.0;
  v2[1] =-1.0;
  v2[2] = 0.0;
  
  v3[0] = 0.0;
  v3[1] = 0.0;
  v3[2] = 1.0;
  peanut_tesselate(U, v1, v2, v3, depth);

  /* -x, -z, -y */
  v1[0] =-1.0;
  v1[1] = 0.0;
  v1[2] = 0.0;
  
  v2[0] = 0.0;
  v2[1] = 0.0;
  v2[2] =-1.0;
  
  v3[0] = 0.0;
  v3[1] =-1.0;
  v3[2] = 0.0;
  peanut_tesselate(U, v1, v2, v3, depth);

  glEnd();
  glPopMatrix();

  Py_INCREF(Py_None);
  return Py_None;
}



/* functions for the rendering of atomic thermal ellipsoids
 */

void
ellipse_func(float U[6], float C, float v[3], float w[3]) 
{
  float d;

  d =     U[0]*v[0]*v[0] +     U[1]*v[1]*v[1] + U[2]*v[2]*v[2]
    + 2.0*U[3]*v[0]*v[1] + 2.0*U[4]*v[0]*v[2] + 2.0*U[5]*v[1]*v[2];

  /* abort before we take the sqare root of a negitive
   * number */
  if (d < 0.0) {
    w[0] = 0.0;
    w[1] = 0.0;
    w[2] = 0.0;

    return;
  }

  d = C * sqrt(1.0 / d);

  w[0] = d * v[0];
  w[1] = d * v[1];
  w[2] = d * v[2];
}

void
ellipse_normal(float U[6], float C, float v[3], float n[3])
{
  float d;

  d =     U[0]*v[0]*v[0] +     U[1]*v[1]*v[1] + U[2]*v[2]*v[2]
    + 2.0*U[3]*v[0]*v[1] + 2.0*U[4]*v[0]*v[2] + 2.0*U[5]*v[1]*v[2];

  d = 1.0 / d;

  /* abort before we take the sqare root of a negitive
   * number */
  if (d < 0.0) {
    n[0] = 0.0;
    n[1] = 0.0;
    n[2] = 0.0;
    
    return;
  }

  d = d * sqrt(d);

  n[0] = 0.5 * C * d * (2.0*U[0]*v[0] +     U[3]*v[1] +    U[4]*v[2]);
  n[1] = 0.5 * C * d * (    U[3]*v[0] + 2.0*U[1]*v[1] +    U[5]*v[2]);
  n[2] = 0.5 * C * d * (    U[4]*v[0] +     U[5]*v[1] +2.0*U[2]*v[2]);

  normalize(n);
}

void
ellipse_triangle(float U[6], float C, float v1[3], float v2[3], float v3[3])
{
  float v[3], n[3];

  ellipse_func(U, C, v1, v);
  ellipse_normal(U, C, v, n);
  glNormal3f(n[0], n[1], n[2]);
  glVertex3f(v[0], v[1], v[2]);

  ellipse_func(U, C, v2, v);
  ellipse_normal(U, C, v, n);
  glNormal3f(n[0], n[1], n[2]);
  glVertex3f(v[0], v[1], v[2]);

  ellipse_func(U, C, v3, v);
  ellipse_normal(U, C, v, n);
  glNormal3f(n[0], n[1], n[2]);
  glVertex3f(v[0], v[1], v[2]);
}

void
ellipse_tesselate(float U[6], float C, 
		  float v1[3], float v2[3], float v3[3], int depth)
{
  float v12[3], v23[3], v31[3];

  v12[0] = ((v2[0] - v1[0]) / 2.0) + v1[0];
  v12[1] = ((v2[1] - v1[1]) / 2.0) + v1[1];
  v12[2] = ((v2[2] - v1[2]) / 2.0) + v1[2];

  v23[0] = ((v3[0] - v2[0]) / 2.0) + v2[0];
  v23[1] = ((v3[1] - v2[1]) / 2.0) + v2[1];
  v23[2] = ((v3[2] - v2[2]) / 2.0) + v2[2];

  v31[0] = ((v1[0] - v3[0]) / 2.0) + v3[0];
  v31[1] = ((v1[1] - v3[1]) / 2.0) + v3[1];
  v31[2] = ((v1[2] - v3[2]) / 2.0) + v3[2];

  normalize(v12);
  normalize(v23);
  normalize(v31);

  depth -= 1;

  if (depth > 0) {
    ellipse_tesselate(U, C, v1,  v12, v31, depth);
    ellipse_tesselate(U, C, v2,  v23, v12, depth);
    ellipse_tesselate(U, C, v3,  v31, v23, depth);
    ellipse_tesselate(U, C, v12, v23, v31, depth);
  } else {
    ellipse_triangle(U, C, v1,  v12, v31);
    ellipse_triangle(U, C, v2,  v23, v12);
    ellipse_triangle(U, C, v3,  v31, v23);
    ellipse_triangle(U, C, v12, v23, v31);
  }
}

static PyObject *
glaccel_Uellipse(PyObject *self, PyObject *args)
{
  int          depth;
  float        x, y, z;
  float        U[6];
  float        C;
  float        v1[3], v2[3], v3[3];

  if (!PyArg_ParseTuple(args, "ffffffffffi", 
			&x, &y, &z,
			&U[0], &U[1], &U[2], &U[3], &U[4], &U[5],
			&C, &depth)) {
    return NULL;
  }

  
  glPushMatrix();
  glTranslatef(x, y, z);
  glEnable(GL_LIGHTING);

  glBegin(GL_TRIANGLES);

  /* x, y, z */
  v1[0] = 1.0;
  v1[1] = 0.0;
  v1[2] = 0.0;
  
  v2[0] = 0.0;
  v2[1] = 1.0;
  v2[2] = 0.0;
  
  v3[0] = 0.0;
  v3[1] = 0.0;
  v3[2] = 1.0;
  ellipse_tesselate(U, C, v1, v2, v3, depth);

  /* x, -z, y */
  v1[0] = 1.0;
  v1[1] = 0.0;
  v1[2] = 0.0;
  
  v2[0] = 0.0;
  v2[1] = 0.0;
  v2[2] =-1.0;
  
  v3[0] = 0.0;
  v3[1] = 1.0;
  v3[2] = 0.0;
  ellipse_tesselate(U, C, v1, v2, v3, depth);

  /* x, z, -y */
  v1[0] = 1.0;
  v1[1] = 0.0;
  v1[2] = 0.0;
  
  v2[0] = 0.0;
  v2[1] = 0.0;
  v2[2] = 1.0;
  
  v3[0] = 0.0;
  v3[1] =-1.0;
  v3[2] = 0.0;
  ellipse_tesselate(U, C, v1, v2, v3, depth);

  /* x, -y, -z */
  v1[0] = 1.0;
  v1[1] = 0.0;
  v1[2] = 0.0;
  
  v2[0] = 0.0;
  v2[1] =-1.0;
  v2[2] = 0.0;
  
  v3[0] = 0.0;
  v3[1] = 0.0;
  v3[2] =-1.0;
  ellipse_tesselate(U, C, v1, v2, v3, depth);

  /* -x, z, y */
  v1[0] =-1.0;
  v1[1] = 0.0;
  v1[2] = 0.0;
  
  v2[0] = 0.0;
  v2[1] = 0.0;
  v2[2] = 1.0;
  
  v3[0] = 0.0;
  v3[1] = 1.0;
  v3[2] = 0.0;
  ellipse_tesselate(U, C, v1, v2, v3, depth);

  /* -x, y, -z */
  v1[0] =-1.0;
  v1[1] = 0.0;
  v1[2] = 0.0;
  
  v2[0] = 0.0;
  v2[1] = 1.0;
  v2[2] = 0.0;
  
  v3[0] = 0.0;
  v3[1] = 0.0;
  v3[2] =-1.0;
  ellipse_tesselate(U, C, v1, v2, v3, depth);

  /* -x, -y, z */
  v1[0] =-1.0;
  v1[1] = 0.0;
  v1[2] = 0.0;
  
  v2[0] = 0.0;
  v2[1] =-1.0;
  v2[2] = 0.0;
  
  v3[0] = 0.0;
  v3[1] = 0.0;
  v3[2] = 1.0;
  ellipse_tesselate(U, C, v1, v2, v3, depth);

  /* -x, -z, -y */
  v1[0] =-1.0;
  v1[1] = 0.0;
  v1[2] = 0.0;
  
  v2[0] = 0.0;
  v2[1] = 0.0;
  v2[2] =-1.0;
  
  v3[0] = 0.0;
  v3[1] =-1.0;
  v3[2] = 0.0;
  ellipse_tesselate(U, C, v1, v2, v3, depth);

  glEnd();
  glPopMatrix();

  Py_INCREF(Py_None);
  return Py_None;
}


/* starndard Python module registration code
 */

static PyMethodDef GLAccelMethods[] = {

  {"Upeanut",
   glaccel_Upeanut,
   METH_VARARGS,
   "Renders a thermal peanut."},
  
  {"Uellipse",
   glaccel_Uellipse,
   METH_VARARGS,
   "Renders a thermal ellipsoid."},
  
  {NULL, NULL, 0, NULL}
};


DL_EXPORT(void)
initglaccel(void)
{
  PyObject *m;
  
  m = Py_InitModule("glaccel", GLAccelMethods);
  
  GLAccelError = PyErr_NewException("glaccel.error", NULL, NULL);
  Py_INCREF(GLAccelError);
  PyModule_AddObject(m, "error", GLAccelError);
}
