// Copyright 2006 by TLSMD Development Group (see AUTHORS file)
// This code is part of the TLSMD distribution and governed by
// its license.  Please see the LICENSE file that should have been
// included as part of this package.
#include "dgesdd.h"

#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))

namespace TLSMD {

/* LAPACK */
extern "C" void
dgesdd_(char *, int*, int*, double*, int*, double*, double*, int *, double*, int*, double*, int*, int*, int*);

DGESDD::DGESDD() : A(0), b(0), S(0), U(0), VT(0), WORK(0), LWORK(0), IWORK(0) {
}

DGESDD::~DGESDD() {
  delete[] A;
  delete[] b;
  delete[] S;
  delete[] U;
  delete[] VT;
  delete[] WORK;
  delete[] IWORK;
}

void
DGESDD::zero() {
  int i, sz;
  
  if (A == 0) return;
  sz = num_rows * num_cols;
  for (i = 0; i < sz; ++i) A[i] = 0.0;
}

void
DGESDD::set_max_matrix_size(int nrows, int ncols) {
  char jobz;
  int info;
  double tmp_WORK;

  // set number of rows and columns
  max_rows = nrows;
  max_cols = ncols;

  A = new double[max_rows * max_cols];
  b = new double[max_rows];
  S = new double[max_cols];
  U = new double[max_rows * max_cols];
  VT = new double[max_cols * max_rows];
  IWORK = new int[8 * max_cols];

  // calculate the ideal size of the WORK memory block
  jobz = 'S';
  LWORK = -1;
  info = 0;
  dgesdd_(&jobz, &max_rows, &max_cols,  A, &max_rows, S, U, &max_rows, VT, &max_cols, &tmp_WORK, &LWORK, IWORK, &info);
  LWORK = (int) tmp_WORK;

  WORK = new double[LWORK];
}

void 
DGESDD::set_matrix_size(int nrows, int ncols) {
  num_rows = nrows;
  num_cols = ncols;
}

bool
DGESDD::svd() {
  char jobz;
  int info;

  jobz = 'S';
  dgesdd_(&jobz, &num_rows, &num_cols, A, &num_rows, S, U, &num_rows, VT, &num_cols, WORK, &LWORK, IWORK, &info);
  if (info != 0) {
    return false;
  }
  return true;
}

// solve for x given the singular value decomposition of a matrix
// into U, S, and Vt 
void 
DGESDD::solve_for_x(double* x) {
#define FU(__i, __j)   U[__i + (m * __j)]
#define FVT(__i, __j)  VT[__i + (n * __j)]

  int m, n, i, j;
  double smax, scutoff, dtmp;

  m = num_rows;
  n = num_cols;

  // invert the diagonal matrix S in place, and filter out any
  // unusually small singular values
  for (smax = S[0], i = 1; i < n; ++i) {
    smax = MAX(smax, S[i]);
  }
  scutoff = smax * 1E-20;
  for (i = 0; i < n; ++i) {
    if (S[i] > scutoff) {
      S[i] = 1.0 / S[i];
    } else {
      S[i] = 0.0;
    }
  }

  // matrix multiply Ut(n,m)*b(m) and store the result in x
  for (i = 0; i < n; i++) {
    dtmp = 0.0;
    for (j = 0; j < m; ++j) {
      dtmp += FU(j,i) * b[j];
    }
    WORK[i] = dtmp;
  }

  // matrix multiply inverse-S by Ut*b
  for (i = 0; i < n; ++i) {
    WORK[i] *= S[i];
  }

  // matrix multiple V*x
  for (i = 0; i < n; ++i) {
    dtmp = 0.0;
    for (j = 0; j < n; ++j) {
      dtmp += FVT(j,i) * WORK[j];
    }
    x[i] = dtmp;
  }

#undef FU
#undef FVT
}

} // namespace TLSMD
