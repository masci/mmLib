// Copyright 2006 by TLSMD Development Group (see AUTHORS file)
// This code is part of the TLSMD distribution and governed by
// its license.  Please see the LICENSE file that should have been
// included as part of this package.
#ifndef __TLS_MODEL_ENGINE__
#define __TLS_MODEL_ENGINE__

#include "structure.h"
#include "tls_model.h"

namespace TLSMD {

class TLSModelEngine {
public:
  void set_num_atoms(int num_atoms);

  void isotropic_fit_segment(int istart, int iend, IsotropicTLSModel &itls_model, double *residual);
  void anisotropic_fit_segment(int istart, int iend, AnisotropicTLSModel &atls_model, double *residual);

  Chain chain;

 private:
  FitIsotropicTLSModel fit_itls;
  FitAnisotropicTLSModel fit_atls;

};

} // namespace TLSMD

#endif // __TLS_MODEL_ENGINE__
