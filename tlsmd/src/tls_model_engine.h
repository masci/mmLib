// Copyright 2006 by TLSMD Development Group (see AUTHORS file)
// This code is part of the TLSMD distribution and governed by
// its license.  Please see the LICENSE file that should have been
// included as part of this package.
#ifndef __TLS_MODEL_ENGINE__
#define __TLS_MODEL_ENGINE__

#include <string>

#include "structure.h"
#include "tls_model.h"
#include "tls_model_nl.h"

namespace TLSMD {

class FitTLSModelResult {
 public:
  FitTLSModelResult() : num_atoms_(0), num_residues_(0), residual_(0.0) {}
  virtual ~FitTLSModelResult() {}

  int get_num_atoms() const { return num_atoms_; }
  void set_num_atoms(int num_atoms) { num_atoms_ = num_atoms; }
 
  int get_num_residues() const { return num_residues_; }
  void set_num_residues(int num_residues) { num_residues_ = num_residues; }
  
  double get_residual() const { return residual_; }
  void set_residual(double residual) { residual_ = residual; }

  virtual TLSModel& get_tls_model() = 0;

 private:
  int num_atoms_;
  int num_residues_;
  double residual_;
};

class IsotropicFitTLSModelResult : public FitTLSModelResult {
 public:
  virtual TLSModel& get_tls_model() { return itls_model; }
  IsotropicTLSModel itls_model;
};

class AnisotropicFitTLSModelResult : public FitTLSModelResult {
 public:
  virtual TLSModel& get_tls_model() { return atls_model; }
  AnisotropicTLSModel atls_model;
};

class TLSModelEngine {
public:
  void set_num_atoms(int num_atoms);

  void isotropic_fit_segment(const std::string& frag_id1,
			     const std::string& frag_id2,
			     IsotropicFitTLSModelResult& itls_result);
			     
  void anisotropic_fit_segment(const std::string& frag_id1,
			       const std::string& frag_id2,
			       AnisotropicFitTLSModelResult& atls_result);

  void constrained_isotropic_fit_segment(const std::string& frag_id1,
					 const std::string& frag_id2,
					 IsotropicFitTLSModelResult& itls_result);

  void constrained_anisotropic_fit_segment(const std::string& frag_id1,
					   const std::string& frag_id2,
					   AnisotropicFitTLSModelResult& atls_result);

  Chain chain;

 private:
  FitIsotropicTLSModel fit_itls;
  FitAnisotropicTLSModel fit_atls;
  ConstrainedFitIsotropicTLSModel cfit_itls;
  ConstrainedFitAnisotropicTLSModel cfit_atls;
};

} // namespace TLSMD

#endif // __TLS_MODEL_ENGINE__
