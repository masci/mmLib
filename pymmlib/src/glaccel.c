/* glaccel.c - OpenGL C accelerator for GLViewer.py
 *
 */
#include "Python.h"
#include <stdio.h>
#include <math.h>

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

static PyObject *GLAccelError = NULL;

/* normalize the vector v
 */
inline void 
normalize(float v[3])
{
  float d;

  d    = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  v[0] = v[0] / d;
  v[1] = v[1] / d;
  v[2] = v[2] / d;
}

/* compute the length of the vector v
 */
inline float
length(float v[3])
{
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

/* form the cross product of vectors u and v, and return the
 * result in w
 */
inline void 
cross(float u[3], float v[3], float w[3])
{
  w[0] = u[1]*v[2] - u[2]*v[1];
  w[1] = u[2]*v[0] - u[0]*v[2];
  w[2] = u[0]*v[1] - u[1]*v[0];
}

/* invert the symmetric 3x3 matrix in the argument U, and return the
 * result in Ui 
 * matrix format: u11,u22,u33,u12,u13,u23
 */
inline int
invert_symmetric_3(float U[6], float Ui[6])
{
  float d;

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

/* create a rotation matrix which will align the Z axis of the OpenGL 
 * coordinate system with the vector u.  R is a column-major 4x4 matrix
 * suited for use as a argument for glMultMatrixf().
 */
inline void
gl_rmatrixz(float u[3], float R[16])
{
  float d;
  float n[3];

  n[0] = u[0];
  n[1] = u[1];
  n[2] = u[2];

  normalize(n);
  d = sqrt(n[0]*n[0] + n[1]*n[1]);

  if (d == 0.0) {
    R[0] = n[2];          /* R11 */
    R[1] = 0.0;           /* R21 */
    R[2] = 0.0;           /* R31 */
    
    R[4] = 0.0;           /* R12 */
    R[5] = 1.0;           /* R22 */
    R[6] = 0.0;           /* R32 */
    
    R[8] = 0.0;          /* R13 */
    R[9] = 0.0;           /* R23 */
    R[10] = n[2];         /* R33 */
  } else {
    R[0] = n[2]*n[0]/d;   /* R11 */
    R[1] = n[2]*n[1]/d;   /* R21 */
    R[2] = -d;            /* R31 */
    
    R[4] = -n[1]/d;       /* R12 */
    R[5] = n[0]/d;        /* R22 */
    R[6] = 0.0;           /* R32 */
    
    R[8] = n[0];          /* R13 */
    R[9] = n[1];          /* R23 */
    R[10] = n[2];         /* R33 */
  }

  R[3]  = 0.0;
  R[7]  = 0.0;
  R[11] = 0.0;
  R[15] = 1.0;
}

/* functions for the rendering of a sphere
 */ 
static PyObject *
glaccel_sphere(PyObject *self, PyObject *args)
{
  int          quality;
  float        x, y, z, radius;

  if (!PyArg_ParseTuple(args, "ffffi", &x, &y, &z, &radius, &quality)) {
    return NULL;
  }
  
  glEnable(GL_LIGHTING);
  glPushMatrix();
  glTranslatef(x, y, z);
  glutSolidSphere(radius, quality, quality);
  glPopMatrix();

  Py_INCREF(Py_None);
  return Py_None;
}

/* functions for the rendering of a tube
 */ 
static PyObject *
glaccel_tube(PyObject *self, PyObject *args)
{
  float        v1[3], v2[3];
  float        radius;
  float        tube_v[3], tube_len;
  float        R[16];
  GLUquadric  *tube_quad;

  if (!PyArg_ParseTuple(args, "fffffff", &v1[0], &v1[1], &v1[2], &v2[0], &v2[1], &v2[2], &radius)) {
    return NULL;
  }
  
  tube_v[0]   = v2[0] - v1[0];
  tube_v[1]   = v2[1] - v1[1];
  tube_v[2]   = v2[2] - v1[2];
  tube_len    = length(tube_v);

  glEnable(GL_LIGHTING);
  tube_quad = gluNewQuadric();
  /*  gluQuadricNormals(tube_quad, GLU_SMOOTH); */

  glPushMatrix();

  /* rotation matrix to align the current OpenGL coordinate
   * system such that the vector axis is along the z axis
   */
  gl_rmatrixz(tube_v, R);
  
  R[12] = v1[0];
  R[13] = v1[1];
  R[14] = v1[2];
  
  glMultMatrixf(R);

  gluCylinder(tube_quad, radius, radius, tube_len , 10, 1);

  glPopMatrix();

  gluDeleteQuadric(tube_quad);

  Py_INCREF(Py_None);
  return Py_None;
}


/* functions for the rendering of a rod
 */ 
static PyObject *
glaccel_rod(PyObject *self, PyObject *args)
{
  float        v1[3], v2[3];
  float        radius;
  float        rod_v[3], rod_len;
  float        R[16];
  GLUquadric  *rod_quad;

  if (!PyArg_ParseTuple(args, "fffffff", &v1[0], &v1[1], &v1[2], &v2[0], &v2[1], &v2[2], &radius)) {
    return NULL;
  }
  
  rod_v[0]   = v2[0] - v1[0];
  rod_v[1]   = v2[1] - v1[1];
  rod_v[2]   = v2[2] - v1[2];
  rod_len    = length(rod_v);

  glEnable(GL_LIGHTING);
  rod_quad = gluNewQuadric();

  glPushMatrix();

  /* rotation matrix to align the current OpenGL coordinate
   * system such that the vector axis is along the z axis
   */
  gl_rmatrixz(rod_v, R);
  
  R[12] = v1[0];
  R[13] = v1[1];
  R[14] = v1[2];
  
  glMultMatrixf(R);

  gluDisk(rod_quad, 0.0, radius, 10, 5);
  gluCylinder(rod_quad, radius, radius, rod_len , 10, 1);
  glTranslatef(0.0, 0.0, rod_len);
  gluDisk(rod_quad, 0.0, radius, 10, 5);

  glPopMatrix();

  gluDeleteQuadric(rod_quad);

  Py_INCREF(Py_None);
  return Py_None;
}


/* functions for the rendering of atomic thermal peanuts
 */
static void
peanut_func(float U[6], float v[3], float w[3]) 
{
  float d;

  d = U[0]*v[0]*v[0] + U[1]*v[1]*v[1] + U[2]*v[2]*v[2] + 2.0*U[3]*v[0]*v[1] + 2.0*U[4]*v[0]*v[2] + 2.0*U[5]*v[1]*v[2];

  /* abort before we take the sqare root of a negative
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

static void
peanut_normal(float U[6], float v[3], float n[3])
{
  float d;

  d = U[0]*v[0]*v[0] + U[1]*v[1]*v[1] + U[2]*v[2]*v[2] + 2.0*U[3]*v[0]*v[1] + 2.0*U[4]*v[0]*v[2] + 2.0*U[5]*v[1]*v[2];

  /* abort before we take the sqare root of a negative
   * number */
  if (d < 0.0) {
    n[0] = v[0];
    n[1] = v[1];
    n[2] = v[2];
    normalize(n);

    return;
  }

  d = 1.0 / (2.0 * sqrt(d));

  n[0] = d * (2.0*U[0]*v[0] +     U[3]*v[1] +    U[4]*v[2]);
  n[1] = d * (    U[3]*v[0] + 2.0*U[1]*v[1] +    U[5]*v[2]);
  n[2] = d * (    U[4]*v[0] +     U[5]*v[1] +2.0*U[2]*v[2]);

  normalize(n);
}

static void
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

static void
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

  if (!PyArg_ParseTuple(args, "fffffffffi", &x, &y, &z,	&U[0], &U[1], &U[2], &U[3], &U[4], &U[5], &depth)) {
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



/* functions for the rendering of atomic thermal ellipsoids */

inline void
ellipse_func(float U[6], float C, float v[3], float w[3]) 
{
  float d;

  d = U[0]*v[0]*v[0] + U[1]*v[1]*v[1] + U[2]*v[2]*v[2] + 2.0*U[3]*v[0]*v[1] + 2.0*U[4]*v[0]*v[2] + 2.0*U[5]*v[1]*v[2];

  /* abort before we take the sqare root of a negative number */
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

inline void
ellipse_normal(float U[6], float C, float v[3], float n[3])
{
  float d;

  d = U[0]*v[0]*v[0] + U[1]*v[1]*v[1] + U[2]*v[2]*v[2] + 2.0*U[3]*v[0]*v[1] + 2.0*U[4]*v[0]*v[2] + 2.0*U[5]*v[1]*v[2];
  d = 1.0 / d;

  /* abort before we take the sqare root of a negative  number */
  if (d < 0.0) {
    n[0] = v[0];
    n[1] = v[1];
    n[2] = v[2];
    normalize(n);

    return;
  }

  d = d * sqrt(d);

  n[0] = 0.5 * C * d * (2.0*U[0]*v[0] +     U[3]*v[1] +    U[4]*v[2]);
  n[1] = 0.5 * C * d * (    U[3]*v[0] + 2.0*U[1]*v[1] +    U[5]*v[2]);
  n[2] = 0.5 * C * d * (    U[4]*v[0] +     U[5]*v[1] +2.0*U[2]*v[2]);

  normalize(n);
}

inline void
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

static void
ellipse_tesselate(float U[6], float C, float v1[3], float v2[3], float v3[3], int depth)
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
  float        U[6], Ui[6];
  float        C;
  float        v1[3], v2[3], v3[3];

  if (!PyArg_ParseTuple(args, "ffffffffffi", &x, &y, &z, &U[0], &U[1], &U[2], &U[3], &U[4], &U[5], &C, &depth)) {
    return NULL;
  }

  if (!invert_symmetric_3(U, Ui)) {
    Py_INCREF(Py_None);
    return Py_None;
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
  ellipse_tesselate(Ui, C, v1, v2, v3, depth);

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
  ellipse_tesselate(Ui, C, v1, v2, v3, depth);

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
  ellipse_tesselate(Ui, C, v1, v2, v3, depth);

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
  ellipse_tesselate(Ui, C, v1, v2, v3, depth);

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
  ellipse_tesselate(Ui, C, v1, v2, v3, depth);

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
  ellipse_tesselate(Ui, C, v1, v2, v3, depth);

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
  ellipse_tesselate(Ui, C, v1, v2, v3, depth);

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
  ellipse_tesselate(Ui, C, v1, v2, v3, depth);

  glEnd();
  glPopMatrix();

  Py_INCREF(Py_None);
  return Py_None;
}


/* starndard Python module registration code
 */

static PyMethodDef GLAccelMethods[] = {

  {"sphere",
   glaccel_sphere,
   METH_VARARGS,
   "Renders a sphere."},

  {"tube",
   glaccel_tube,
   METH_VARARGS,
   "Renders a tube."},

  {"rod",
   glaccel_rod,
   METH_VARARGS,
   "Renders a rod."},

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
