// Copyright 2006 by TLSMD Development Group (see AUTHORS file)
// This code is part of the TLSMD distribution and governed by
// its license.  Please see the LICENSE file that should have been
// included as part of this package.
#ifndef __DGESDD_H__
#define __DGESDD_H__

// wrapper for LAPACK DGESDD subroutine

namespace TLSMD {

class DGESDD {
 public:

  DGESDD();
  ~DGESDD();

  void zero();
  void set_max_matrix_size(int nrows, int ncols);
  void set_matrix_size(int nrows, int ncols);
  bool svd();
  void solve_for_x(double *x);

  int max_rows, max_cols;
  int num_rows, num_cols;

  double *A;
  double *b;
  double *S;
  double *U;
  double *VT;

 private:
  double *WORK;
  int     LWORK;
  int    *IWORK;
};

} // namespace TLSMD

#endif // __DGESDD_H__
