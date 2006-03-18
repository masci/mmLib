// Copyright 2006 by TLSMD Development Group (see AUTHORS file)
// This code is part of the TLSMD distribution and governed by
// its license.  Please see the LICENSE file that should have been
// included as part of this package.
#include <stdio.h>
#include "tls_model_engine.h"

namespace TLSMD {

void
FitTLSModel(const Chain &chain, int group_id, IFitTLSModel &tls_fit, TLSModel &tls_model) {
  // calculate the number of atoms to be fit
  int num_atoms = chain.calc_group_num_atoms(group_id);

  // calculate the centroid of the atoms which are to
  // be fit and use it as the origin of the TLS tensors
  double x, y, z;
  chain.calc_group_centroid(group_id, &x, &y, &z);
  tls_model.set_origin(x, y, z);

  // set the datapoints for the fit
  tls_fit.reset_fit(&tls_model, num_atoms);
 
  std::vector<Atom>::const_iterator atom;
  for (atom = chain.atoms.begin(); atom != chain.atoms.end(); ++atom) {
    if (!atom->in_group(group_id)) continue;
    tls_fit.set_atom_data_point(*atom);
  }

  tls_fit.fit_params();
}

double
IsotropicTLSResult(const Chain& chain, int group_id, IsotropicFitTLSModelResult& itls_result) {
  int num_atoms = 0;
  double chi2 = 0.0;
  double sum_weight = 0.0;

  std::vector<Atom>::const_iterator atom;
  for (atom = chain.atoms.begin(); atom != chain.atoms.end(); ++atom) {
    if (!atom->in_group(group_id)) continue;
    ++num_atoms;

    double uiso_tls;
    itls_result.itls_model.calc_uiso(atom->x, atom->y, atom->z, &uiso_tls);
    
    double delta = uiso_tls - atom->u_iso;
    chi2 += atom->weight * (delta * delta);
    sum_weight += atom->weight;
  }

  itls_result.set_num_atoms(num_atoms);

  int num_residues = chain.calc_group_num_residues(group_id);
  itls_result.set_num_residues(num_residues);

  double residual = num_residues * (chi2 / sum_weight);
  itls_result.set_residual(residual);

  return residual;
}

void
AnisotropicTLSResult(const Chain& chain, int group_id, AnisotropicFitTLSModelResult& atls_result) {
  int num_atoms = 0;
  double chi2 = 0.0;
  double sum_weight = 0.0;

  std::vector<Atom>::const_iterator atom;
  for (atom = chain.atoms.begin(); atom != chain.atoms.end(); ++atom) {
    if (!atom->in_group(group_id)) continue;
    ++num_atoms;

    double Utls[6];
    atls_result.atls_model.calc_U(atom->x, atom->y, atom->z, Utls);
 
    // calculated residule is of trace only
    double delta = ((Utls[0] + Utls[1] + Utls[2]) - (atom->U[0] + atom->U[1] + atom->U[2])) / 3.0;
    chi2 += atom->weight * (delta * delta);
    sum_weight += atom->weight;
  }

  atls_result.set_num_atoms(num_atoms);

  int num_residues = chain.calc_group_num_residues(group_id);
  atls_result.set_num_residues(num_residues);

  double residual = num_residues * (chi2 / sum_weight);
  atls_result.set_residual(residual);
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
TLSModelEngine::isotropic_fit_segment(const std::string& frag_id1,
				      const std::string& frag_id2, 
				      IsotropicFitTLSModelResult& itls_result) {
  int group_id = 1;
  chain.set_group_range(group_id, frag_id1, frag_id2);
  FitTLSModel(chain, group_id, fit_itls, itls_result.get_tls_model());
  IsotropicTLSResult(chain, group_id, itls_result);
}

void
TLSModelEngine::anisotropic_fit_segment(const std::string& frag_id1,
					const std::string& frag_id2,
					AnisotropicFitTLSModelResult& atls_result) {
  int group_id = 1;
  chain.set_group_range(group_id, frag_id1, frag_id2);
  FitTLSModel(chain, group_id, fit_atls, atls_result.get_tls_model());
  AnisotropicTLSResult(chain, group_id, atls_result);
}

void
TLSModelEngine::constrained_isotropic_fit_segment(const std::string& frag_id1,
						  const std::string& frag_id2,
						  IsotropicFitTLSModelResult& itls_result) {
  int group_id = 1;
  chain.set_group_range(group_id, frag_id1, frag_id2);
  FitTLSModel(chain, group_id, cfit_itls, itls_result.get_tls_model());
  IsotropicTLSResult(chain, group_id, itls_result);
}

void
TLSModelEngine::constrained_anisotropic_fit_segment(const std::string& frag_id1,
						    const std::string& frag_id2,
						    AnisotropicFitTLSModelResult& atls_result) {
  int group_id = 1;
  chain.set_group_range(group_id, frag_id1, frag_id2);
  FitTLSModel(chain, group_id, cfit_atls, atls_result.get_tls_model());
  AnisotropicTLSResult(chain, group_id, atls_result);
}

} // namespace TLSMD
