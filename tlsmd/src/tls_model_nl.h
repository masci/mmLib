// Copyright 2006 by TLSMD Development Group (see AUTHORS file)
// This code is part of the TLSMD distribution and governed by
// its license.  Please see the LICENSE file that should have been
// included as part of this package.
#ifndef __TLS_MODEL_NL_H__
#define __TLS_MODEL_NL_H__

#include <vector>
#include "structure.h"
#include "tls_model.h"

// isotropic non-linear TLS Model
#define NL_ITLS_T     0
#define NL_ITLS_LX    1
#define NL_ITLS_LY    2
#define NL_ITLS_LZ    3
#define NL_ITLS_LA    4
#define NL_ITLS_LB    5
#define NL_ITLS_LC    6
#define NL_ITLS_S1    7
#define NL_ITLS_S2    8
#define NL_ITLS_S3    9
#define NL_ITLS_NUM_PARAMS 10

// anisotropic non-linear TLS model parameter indexes and labels
#define NL_ATLS_T11   0
#define NL_ATLS_T22   1
#define NL_ATLS_T33   2
#define NL_ATLS_T12   3
#define NL_ATLS_T13   4
#define NL_ATLS_T23   5
#define NL_ATLS_LX    6
#define NL_ATLS_LY    7
#define NL_ATLS_LZ    8
#define NL_ATLS_LA    9
#define NL_ATLS_LB    10
#define NL_ATLS_LC    11
#define NL_ATLS_S2211 12
#define NL_ATLS_S1133 13
#define NL_ATLS_S12   14
#define NL_ATLS_S13   15
#define NL_ATLS_S23   16
#define NL_ATLS_S21   17
#define NL_ATLS_S31   18
#define NL_ATLS_S32   19
#define NL_ATLS_NUM_PARAMS 20

namespace TLSMD {

struct IsotropicDataPoint {
  double x, y, z;
  double uiso;
};

struct AnisotropicDataPoint {
  double x, y, z;
  double U[6];
};

class ConstrainedFitTLSModel : public IFitTLSModel {
 public:
  ConstrainedFitTLSModel();
  virtual ~ConstrainedFitTLSModel();

 protected:
  void set_size(int nrows, int ncols);

  static const double tolerance = 1E-4;

  int num_rows;
  int num_cols;

  double *fvec;
  double *fjac;
  int *ipvt;
  int lwa;
  double *wa;

  int max_num_atoms;
  int iatom;
  TLSModel *tls_model;
};

class ConstrainedFitIsotropicTLSModel : public ConstrainedFitTLSModel {
 public:
  virtual void set_max_num_atoms(int num_atoms);
  virtual void reset_fit(TLSModel *tls_model, int num_atoms);
  virtual void set_atom_data_point(const Atom& atom) {
    set_data_point(atom.x, atom.y, atom.z, atom.u_iso);
  }
  virtual void fit_params();
  void set_data_point(double x, double y, double z, double uiso);

  std::vector<IsotropicDataPoint> idata_vector;
};

class ConstrainedFitAnisotropicTLSModel : public ConstrainedFitTLSModel {
 public:
  virtual void set_max_num_atoms(int num_atoms);
  virtual void reset_fit(TLSModel *tls_model, int num_atoms);
  virtual void set_atom_data_point(const Atom& atom) {
    set_data_point(atom.x, atom.y, atom.z, atom.U);
  }
  virtual void fit_params();
  void set_data_point(double x, double y, double z, const double U[6]);

  std::vector<AnisotropicDataPoint> adata_vector;
};

} // namespace TLSMD

#endif // __TLS_MODEL_NL_H__
