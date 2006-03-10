// Copyright 2006 by TLSMD Development Group (see AUTHORS file)
// This code is part of the TLSMD distribution and governed by
// its license.  Please see the LICENSE file that should have been
// included as part of this package.
#ifndef __TLS_MODEL_ENGINE__
#define __TLS_MODEL_ENGINE__

#include "structure.h"
#include "tls_model.h"

class TLSModelEngine {
public:
  void set_num_atoms(int num_atoms);

  void isotropic_fit_group(int group_id);
  double isotropic_group_residual(int group_id);
  void isotropic_fit_segment(int istart, int iend, double *residual);

  void anisotropic_fit_group(int group_id);
  double anisotropic_group_residual(int group_id);
  void anisotropic_fit_segment(int istart, int iend, double *residual);

  Chain chain;
  IsotropicTLSModel itls;
  AnisotropicTLSModel atls;

 private:
  void fit_group(int group_id, TLSModel &tls, double parameters[]);
};

#endif // __TLS_MODEL_ENGINE__
