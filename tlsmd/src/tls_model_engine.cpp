// Copyright 2006 by TLSMD Development Group (see AUTHORS file)
// This code is part of the TLSMD distribution and governed by
// its license.  Please see the LICENSE file that should have been
// included as part of this package.
#include <stdio.h>
#include "tls_model_engine.h"

namespace TLSMD {

void
FitTLSModel(Chain &chain, int group_id, IFitTLSModel &tls_fit, TLSModel &tls_model) {
  // calculate the number of atoms to be fit
  int num_atoms = chain.calc_group_num_atoms(group_id);

  // calculate the centroid of the atoms which are to
  // be fit and use it as the origin of the TLS tensors
  double x, y, z;
  chain.calc_group_centroid(group_id, &x, &y, &z);
  tls_model.set_origin(x, y, z);

  // set the datapoints for the fit
  tls_fit.reset_fit(&tls_model, num_atoms);
 
  Atom *atom = chain.atoms;
  for (int ia = 0; ia < chain.num_atoms; ++ia, ++atom) {
    if (!atom->in_group(group_id)) continue;
    tls_fit.set_atom_data_point(atom);
  }

  tls_fit.fit_params();
}

double
IsotropicTLSResidual(Chain &chain, int group_id, IsotropicTLSModel &itls_model) {
  double chi2 = 0.0;
  double sum_weight = 0.0;

  Atom *atom = chain.atoms;
  for (int ia = 0; ia < chain.num_atoms; ++ia, ++atom) {
    if (!atom->in_group(group_id)) continue;
    
    double uiso_tls;
    itls_model.calc_uiso(atom->x, atom->y, atom->z, &uiso_tls);
    
    double delta = uiso_tls - atom->u_iso;
    chi2 += atom->weight * (delta * delta);
    sum_weight += atom->weight;
  }

  int num_residues = chain.calc_group_num_residues(group_id);
  return num_residues * (chi2 / sum_weight);
}

double
AnisotropicTLSResidual(Chain &chain, int group_id, AnisotropicTLSModel &atls_model) {
  double chi2 = 0.0;
  double sum_weight = 0.0;

  Atom *atom = chain.atoms;
  for (int ia = 0; ia < chain.num_atoms; ++ia, ++atom) {
    if (!atom->in_group(group_id)) continue;
    
    double Utls[6];
    atls_model.calc_U(atom->x, atom->y, atom->z, Utls);
 
    // calculated residule is of trace only
    double delta = ((Utls[0] + Utls[1] + Utls[2]) - (atom->U[0] + atom->U[1] + atom->U[2])) / 3.0;
    chi2 += atom->weight * (delta * delta);
    sum_weight += atom->weight;
  }

  int num_residues = chain.calc_group_num_residues(group_id);
  return num_residues * (chi2 / sum_weight);
}

void 
TLSModelEngine::set_num_atoms(int num_atoms) {
  chain.set_num_atoms(num_atoms);
  fit_itls.set_max_num_atoms(num_atoms);
  fit_atls.set_max_num_atoms(num_atoms);
  cfit_itls.set_max_num_atoms(num_atoms);
  cfit_atls.set_max_num_atoms(num_atoms);
}

void
TLSModelEngine::isotropic_fit_segment(int istart, int iend, IsotropicTLSModel &itls_model, double *residual) {
  int group_id = 1;
  chain.set_group_range(group_id, istart, iend);
  FitTLSModel(chain, group_id, fit_itls, itls_model);
  *residual = IsotropicTLSResidual(chain, group_id, itls_model);
}

void
TLSModelEngine::anisotropic_fit_segment(int istart, int iend, AnisotropicTLSModel &atls_model, double *residual)
{
  int group_id = 1;
  chain.set_group_range(group_id, istart, iend);
  FitTLSModel(chain, group_id, fit_atls, atls_model);
  *residual = AnisotropicTLSResidual(chain, group_id, atls_model);
}

void
TLSModelEngine::constrained_isotropic_fit_segment(int istart, int iend, IsotropicTLSModel &itls_model, double *residual) {
  int group_id = 1;
  chain.set_group_range(group_id, istart, iend);
  FitTLSModel(chain, group_id, cfit_itls, itls_model);
  *residual = IsotropicTLSResidual(chain, group_id, itls_model);
}

void
TLSModelEngine::constrained_anisotropic_fit_segment(int istart, int iend, AnisotropicTLSModel &atls_model, double *residual)
{
  int group_id = 1;
  chain.set_group_range(group_id, istart, iend);
  FitTLSModel(chain, group_id, cfit_atls, atls_model);
  *residual = AnisotropicTLSResidual(chain, group_id, atls_model);
}

} // namespace TLSMD
