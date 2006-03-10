// Copyright 2006 by TLSMD Development Group (see AUTHORS file)
// This code is part of the TLSMD distribution and governed by
// its license.  Please see the LICENSE file that should have been
// included as part of this package.
#include "tls_model.h"

void
IsotropicTLSModel::calc_uiso(double x, double y, double z, double *u_iso) {
  x -= origin_x;
  y -= origin_y;
  z -= origin_z;

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
void 
IsotropicTLSModel::reset_fit(int num_atoms) {
  row = 0;
  Axb.set_matrix_size(num_atoms, ITLS_NUM_PARAMS);
  Axb.zero();
}

void
IsotropicTLSModel::set_atom_params(double x, double y, double z, double uiso, double w) {
  x -= origin_x;
  y -= origin_y;
  z -= origin_z;

  int m;
  double *A, *b;
  A = Axb.A;
  b = Axb.b;
  m = Axb.num_rows;
#define FA(__i, __j) A[__i + (m * __j)]

  double xx, yy, zz, xy, xz, yz;
  xx = x*x;
  yy = y*y;
  zz = z*z;
  xy = x*y;
  xz = x*z;
  yz = y*z;

  // set b
  b[row] = w * uiso;

  // T iso
  FA(row, ITLS_T) = w * 1.0;

  // l11, l22, l33, l12, l13, l23
  FA(row, ITLS_L11) = w * ((zz + yy) / 3.0);
  FA(row, ITLS_L22) = w * ((xx + zz) / 3.0);
  FA(row, ITLS_L33) = w * ((xx + yy) / 3.0);
  
  FA(row, ITLS_L12) = w * ((-2.0 * xy) / 3.0);
  FA(row, ITLS_L13) = w * ((-2.0 * xz) / 3.0);
  FA(row, ITLS_L23) = w * ((-2.0 * yz) / 3.0);

  FA(row, ITLS_S1)  = w * (( 2.0 * z) / 3.0);
  FA(row, ITLS_S2)  = w * (( 2.0 * y) / 3.0);
  FA(row, ITLS_S3)  = w * (( 2.0 * x) / 3.0);

  next_atom();
}


// return the anisotropic TLS model predicted ADP in U for a atom 
// located at coordinates x,y,z with respect to the ATLS origin
void
AnisotropicTLSModel::calc_U(double x, double y, double z, double U[6]) {
  x -= origin_x;
  y -= origin_y;
  z -= origin_z;

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


// Sets the six rows of matrix A starting at A[i,j] with the anistropic
// TLS model coefficents for a atom located at t position x, y, z with
// least-squares weight w.  Matrix A is filled to coumn j+12.
//
// The matrix A(m,n) is filled in FORTRAN-style, that is, assuming
// a column-major memory layout.  That's because, athough the
// matrix is contructed in C, the SVD subroutine from LAPACK is 
// written in FORTRAN.
void 
AnisotropicTLSModel::reset_fit(int num_atoms) {
  row = 0;
  Axb.set_matrix_size(num_atoms * 6, ATLS_NUM_PARAMS);
  Axb.zero();
}

void
AnisotropicTLSModel::set_atom_params(double x, double y, double z, double U[6], double w) {
  x -= origin_x;
  y -= origin_y;
  z -= origin_z;

  int m;
  double *A, *b;
  A = Axb.A;
  b = Axb.b;
  m = Axb.num_rows;
#define FA(__i,__j) A[__i + (m * __j)]

  double xx, yy, zz, xy, xz, yz;
  xx = x*x;
  yy = y*y;
  zz = z*z;
  xy = x*y;
  xz = x*z;
  yz = y*z;

  // calculate row indexes
  int rowU11, rowU22, rowU33, rowU12, rowU13, rowU23;
  rowU11 = row;
  rowU22 = row + 1;
  rowU33 = row + 2;
  rowU12 = row + 3;
  rowU13 = row + 4;
  rowU23 = row + 5;

  // set b
  b[rowU11] = w * U[0];
  b[rowU22] = w * U[1];
  b[rowU33] = w * U[2];
  b[rowU12] = w * U[3];
  b[rowU13] = w * U[4];
  b[rowU23] = w * U[5];

  // set A
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

  next_atom();
}
