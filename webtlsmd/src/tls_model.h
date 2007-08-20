// Copyright 2006 by TLSMD Development Group (see AUTHORS file)
// This code is part of the TLSMD distribution and governed by
// its license.  Please see the LICENSE file that should have been
// included as part of this package.
#ifndef __TLS_MODEL_H__
#define __TLS_MODEL_H__

#include "dgesdd.h"
#include "structure.h"

#define PI    3.1415926535897931
#define SQRT2 1.4142135623730951
#define U2B   ((8.0*PI*PI))
#define B2U   ((1.0/U2B))

// anisotropic U tensor parameter labels and indexes
#define U11 0
#define U22 1
#define U33 2
#define U12 3
#define U13 4
#define U23 5
#define U_NUM_PARAMS 6

// Isotropic TLS Model 
// S1 == S21-S12; S2 == S13-S31; S3 == S32-S23
#define ITLS_T     0
#define ITLS_L11   1
#define ITLS_L22   2
#define ITLS_L33   3
#define ITLS_L12   4
#define ITLS_L13   5
#define ITLS_L23   6
#define ITLS_S1    7
#define ITLS_S2    8
#define ITLS_S3    9
#define ITLS_NUM_PARAMS 10

// Anisotropic TLS model parameter indexes and labels
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

namespace TLSMD {

class TLSModel {
 public:
  TLSModel();
  virtual ~TLSModel() {}

  virtual int num_params() const = 0;
  virtual double* get_params() = 0;
  virtual const double* get_params() const = 0;

  void set_origin(const double x, const double y, const double z) {
    origin_x = x; origin_y = y, origin_z = z;
  }

  double origin_x;
  double origin_y;
  double origin_z;
};

class IsotropicTLSModel : public TLSModel {
 public:
  virtual int num_params() const { return ITLS_NUM_PARAMS; }
  virtual double* get_params() { return ITLS; }
  virtual const double* get_params() const { return ITLS; }
  void calc_uiso(double x, double y, double z, double *uiso) const;
  double ITLS[ITLS_NUM_PARAMS];
};

class AnisotropicTLSModel : public TLSModel {
 public:
  virtual int num_params() const { return ATLS_NUM_PARAMS; }
  virtual double* get_params() { return ATLS; }
  virtual const double* get_params() const { return ATLS; }
  void calc_U(double x, double y, double z, double U[]) const;
  double ATLS[ATLS_NUM_PARAMS];
};

class IFitTLSModel {
 public:
  virtual ~IFitTLSModel() {}
  virtual void set_max_num_atoms(int num_atoms)  = 0;
  virtual void reset_fit(TLSModel *tls_model, int num_atoms) = 0;
  virtual void set_atom_data_point(const Atom &atom) = 0;
  virtual void fit_params() = 0;
};

class FitTLSModel : public IFitTLSModel {
 public:
  FitTLSModel();
  virtual void fit_params();

 protected:
  int max_num_atoms;
  TLSModel *tls_model;
  int row;
  DGESDD Axb;
};

class FitIsotropicTLSModel : public FitTLSModel {
 public:
  virtual void set_max_num_atoms(int num_atoms);
  virtual void reset_fit(TLSModel *tls_model, int num_atoms);
  virtual void set_atom_data_point(const Atom &atom) {
    set_data_point(atom.x, atom.y, atom.z, atom.u_iso, atom.sqrt_weight);
  }
  void set_data_point(double x, double y, double z, double uiso, double sqrt_weight);
};

class FitAnisotropicTLSModel : public FitTLSModel {
 public:
  virtual void set_max_num_atoms(int num_atoms);
  virtual void reset_fit(TLSModel *tls_model, int num_atoms);
  virtual void set_atom_data_point(const Atom &atom) {
    set_data_point(atom.x, atom.y, atom.z, atom.U, atom.sqrt_weight);
  }
  void set_data_point(double x, double y, double z, const double U[6], double sqrt_weight);
};

void CalcIsotropicTLSModelUIso(const double ITLS[], double x, double y, double z, double *uiso);
void CalcAnisotropicTLSModelU(const double ATLS[], double x, double y, double z, double U[]);

} // namespace TLSMD

#endif // __TLS_MODEL_H__
