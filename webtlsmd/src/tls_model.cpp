// Copyright 2006 by TLSMD Development Group (see AUTHORS file)
// This code is part of the TLSMD distribution and governed by
// its license.  Please see the LICENSE file that should have been
// included as part of this package.
#include "tls_model.h"

namespace TLSMD {

void
CalcIsotropicTLSModelUIso(const double ITLS[], double x, double y, double z, double *uiso) {
  double xx, yy, zz;
  xx = x*x;
  yy = y*y;
  zz = z*z;

  // note: S1 == S21-S12; S2 == S13-S31; S3 == S32-S23
  *uiso = ITLS[ITLS_T]                      +
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

// return the anisotropic TLS model predicted ADP in U for a atom 
// located at coordinates x,y,z with respect to the ATLS origin
void
CalcAnisotropicTLSModelU(const double ATLS[], double x, double y, double z, double U[]) {
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

TLSModel::TLSModel()
  : origin_x(0.0), origin_y(0.0), origin_z(0.0) {
}

void
IsotropicTLSModel::calc_uiso(double x, double y, double z, double *uiso) const {
  CalcIsotropicTLSModelUIso(ITLS, x - origin_x, y - origin_y, z - origin_z, uiso);
}

void
AnisotropicTLSModel::calc_U(double x, double y, double z, double U[6]) const {
  CalcAnisotropicTLSModelU(ATLS, x - origin_x, y - origin_y, z - origin_z, U);
}

FitTLSModel::FitTLSModel() 
  : max_num_atoms(0), tls_model(0), row(0) {
}

void
FitTLSModel::fit_params() {
  Axb.svd();
  Axb.solve_for_x(tls_model->get_params());
}

void
FitIsotropicTLSModel::set_max_num_atoms(int num_atoms) {
  max_num_atoms = num_atoms;
  Axb.set_max_matrix_size(num_atoms, ITLS_NUM_PARAMS);
}

void 
FitIsotropicTLSModel::reset_fit(TLSModel *tls_model, int num_atoms) {
  if (num_atoms > max_num_atoms) set_max_num_atoms(num_atoms);
  FitTLSModel::tls_model = tls_model;
  row = 0;
  Axb.set_matrix_size(num_atoms, ITLS_NUM_PARAMS);
  Axb.zero();
}

// Sets the one row of matrix A starting at A[i,j] with the istropic
// TLS model coefficents for a atom located at t position x, y, z with
// least-squares weight w.  Matrix A is filled to coumn j+12.
// 
// The matrix A(m,n) is filled in FORTRAN-style, that is, assuming
// a column-major memory layout.  That's because, athough the
// matrix is contructed in C, the SVD subroutine from LAPACK is 
// written in FORTRAN.
void
FitIsotropicTLSModel::set_data_point(double x, double y, double z, double uiso, double sqrt_weight) {
  x -= tls_model->origin_x;
  y -= tls_model->origin_y;
  z -= tls_model->origin_z;

  double *A = Axb.A;
  double *b = Axb.b;
  int m = Axb.num_rows;
#define FA(__i, __j) A[__i + (m * __j)]

  double xx, yy, zz, xy, xz, yz;
  xx = x*x;
  yy = y*y;
  zz = z*z;
  xy = x*y;
  xz = x*z;
  yz = y*z;

  // set b
  b[row] = sqrt_weight * uiso;

  // T iso
  FA(row, ITLS_T) = sqrt_weight * 1.0;

  // l11, l22, l33, l12, l13, l23
  FA(row, ITLS_L11) = sqrt_weight * ((zz + yy) / 3.0);
  FA(row, ITLS_L22) = sqrt_weight * ((xx + zz) / 3.0);
  FA(row, ITLS_L33) = sqrt_weight * ((xx + yy) / 3.0);
  
  FA(row, ITLS_L12) = sqrt_weight * ((-2.0 * xy) / 3.0);
  FA(row, ITLS_L13) = sqrt_weight * ((-2.0 * xz) / 3.0);
  FA(row, ITLS_L23) = sqrt_weight * ((-2.0 * yz) / 3.0);

  FA(row, ITLS_S1)  = sqrt_weight * (( 2.0 * z) / 3.0);
  FA(row, ITLS_S2)  = sqrt_weight * (( 2.0 * y) / 3.0);
  FA(row, ITLS_S3)  = sqrt_weight * (( 2.0 * x) / 3.0);

  ++row;
}

void
FitAnisotropicTLSModel::set_max_num_atoms(int num_atoms) {
  max_num_atoms = num_atoms;
  Axb.set_max_matrix_size(num_atoms * 6, ATLS_NUM_PARAMS);
}

void 
FitAnisotropicTLSModel::reset_fit(TLSModel *tls_model, int num_atoms) {
  if (num_atoms > max_num_atoms) set_max_num_atoms(num_atoms);
  FitTLSModel::tls_model = tls_model;
  row = 0;
  Axb.set_matrix_size(6 * num_atoms, ATLS_NUM_PARAMS);
  Axb.zero();
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
FitAnisotropicTLSModel::set_data_point(double x, double y, double z, const double U[6], double sqrt_weight) {
  x -= tls_model->origin_x;
  y -= tls_model->origin_y;
  z -= tls_model->origin_z;

  double *A = Axb.A;
  double *b = Axb.b;
  int m = Axb.num_rows;
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
  b[rowU11] = sqrt_weight * U[0];
  b[rowU22] = sqrt_weight * U[1];
  b[rowU33] = sqrt_weight * U[2];
  b[rowU12] = sqrt_weight * U[3];
  b[rowU13] = sqrt_weight * U[4];
  b[rowU23] = sqrt_weight * U[5];

  // set A
  FA(rowU11, ATLS_T11) = sqrt_weight * 1.0;
  FA(rowU11, ATLS_L22) = sqrt_weight *        zz;
  FA(rowU11, ATLS_L33) = sqrt_weight *        yy;
  FA(rowU11, ATLS_L23) = sqrt_weight * -2.0 * yz;
  FA(rowU11, ATLS_S31) = sqrt_weight * -2.0 *  y;
  FA(rowU11, ATLS_S21) = sqrt_weight *  2.0 *  z;

  FA(rowU22, ATLS_T22) = sqrt_weight * 1.0;
  FA(rowU22, ATLS_L11) = sqrt_weight *        zz;
  FA(rowU22, ATLS_L33) = sqrt_weight *        xx;
  FA(rowU22, ATLS_L13) = sqrt_weight * -2.0 * xz;
  FA(rowU22, ATLS_S12) = sqrt_weight * -2.0 *  z;
  FA(rowU22, ATLS_S32) = sqrt_weight *  2.0 *  x;

  FA(rowU33, ATLS_T33) = sqrt_weight * 1.0;
  FA(rowU33, ATLS_L11) = sqrt_weight *        yy;
  FA(rowU33, ATLS_L22) = sqrt_weight *        xx;
  FA(rowU33, ATLS_L12) = sqrt_weight * -2.0 * xy;
  FA(rowU33, ATLS_S23) = sqrt_weight * -2.0 *  x;
  FA(rowU33, ATLS_S13) = sqrt_weight *  2.0 *  y;

  FA(rowU12, ATLS_T12)   = sqrt_weight * 1.0;
  FA(rowU12, ATLS_L33)   = sqrt_weight * -xy;
  FA(rowU12, ATLS_L23)   = sqrt_weight *  xz;
  FA(rowU12, ATLS_L13)   = sqrt_weight *  yz;
  FA(rowU12, ATLS_L12)   = sqrt_weight * -zz;
  FA(rowU12, ATLS_S2211) = sqrt_weight *   z;
  FA(rowU12, ATLS_S31)   = sqrt_weight *   x;
  FA(rowU12, ATLS_S32)   = sqrt_weight *  -y;
    
  FA(rowU13, ATLS_T13)   = sqrt_weight * 1.0;
  FA(rowU13, ATLS_L22)   = sqrt_weight * -xz;
  FA(rowU13, ATLS_L23)   = sqrt_weight *  xy;
  FA(rowU13, ATLS_L13)   = sqrt_weight * -yy;
  FA(rowU13, ATLS_L12)   = sqrt_weight *  yz;
  FA(rowU13, ATLS_S1133) = sqrt_weight *   y;
  FA(rowU13, ATLS_S23)   = sqrt_weight *   z;
  FA(rowU13, ATLS_S21)   = sqrt_weight *  -x;
    
  FA(rowU23, ATLS_T23)   = sqrt_weight * 1.0;
  FA(rowU23, ATLS_L11)   = sqrt_weight * -yz;
  FA(rowU23, ATLS_L23)   = sqrt_weight * -xx;
  FA(rowU23, ATLS_L13)   = sqrt_weight *  xy;
  FA(rowU23, ATLS_L12)   = sqrt_weight *  xz;
  FA(rowU23, ATLS_S2211) = sqrt_weight *  -x;
  FA(rowU23, ATLS_S1133) = sqrt_weight *  -x;
  FA(rowU23, ATLS_S12)   = sqrt_weight *   y;
  FA(rowU23, ATLS_S13)   = sqrt_weight *  -z;
#undef FA

  row += 6;
}

} // namespace TLSMD

