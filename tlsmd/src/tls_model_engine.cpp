// Copyright 2006 by TLSMD Development Group (see AUTHORS file)
// This code is part of the TLSMD distribution and governed by
// its license.  Please see the LICENSE file that should have been
// included as part of this package.
#include <stdio.h>
#include "tls_model_engine.h"


void 
TLSModelEngine::set_num_atoms(int num_atoms) {
  chain.set_num_atoms(num_atoms);
  itls.Axb.set_max_matrix_size(num_atoms, ITLS_NUM_PARAMS);
  atls.Axb.set_max_matrix_size(num_atoms * 6, ATLS_NUM_PARAMS);
}

void
TLSModelEngine::fit_group(int group_id, TLSModel &tls, double parameters[]) {
  // calculate the number of atoms to be fit
  int num_atoms = chain.calc_group_num_atoms(group_id);
  tls.reset_fit(num_atoms);
  
  // calculate the centroid of the atoms which are to
  // be fit and use it as the origin of the TLS tensors
  chain.calc_group_centroid(group_id, &tls.origin_x, &tls.origin_y, &tls.origin_z);

  Atom *atom = chain.atoms;
  for (int ia = 0; ia < chain.num_atoms; ++ia, ++atom) {
    if (!atom->in_group(group_id)) continue;
    tls.set_atom_params(atom);
  }

  tls.Axb.svd();
  tls.Axb.solve_for_x(parameters);
}

void
TLSModelEngine::isotropic_fit_group(int group_id) {
  fit_group(group_id, itls, itls.ITLS);
}

double
TLSModelEngine::isotropic_group_residual(int group_id) {
  double chi2 = 0.0;
  double sum_weight = 0.0;

  Atom *atom = chain.atoms;
  for (int ia = 0; ia < chain.num_atoms; ++ia, ++atom) {
    if (!atom->in_group(group_id)) continue;
    
    double u_iso_tls;
    itls.calc_uiso(atom->x, atom->y, atom->z, &u_iso_tls);
 
    double tmp = atom->sqrt_weight * (u_iso_tls - atom->u_iso);
    chi2 += tmp * tmp;
    sum_weight += atom->sqrt_weight;
  }

  int num_residues = chain.calc_group_num_residues(group_id);
  return num_residues * (chi2 / sum_weight);
}

void
TLSModelEngine::isotropic_fit_segment(int istart, int iend, double *residual) {
  int group_id = 1;
  chain.set_group_range(group_id, istart, iend);
  isotropic_fit_group(group_id);
  *residual = isotropic_group_residual(group_id);
}

void
TLSModelEngine::anisotropic_fit_group(int group_id)
{
  fit_group(group_id, atls, atls.ATLS);
}

double
TLSModelEngine::anisotropic_group_residual(int group_id) {
  double chi2 = 0.0;
  double sum_weight = 0.0;

  Atom *atom = chain.atoms;
  for (int ia = 0; ia < chain.num_atoms; ++ia, ++atom) {
    if (!atom->in_group(group_id)) continue;
    
    double *U = atom->U;
    double Utls[6];
    atls.calc_U(atom->x, atom->y, atom->z, Utls);
 
    // calculated residule is of trace only
    double delta = ((Utls[0] + Utls[1] + Utls[2]) / 3.0) - ((U[0] + U[1] + U[2]) / 3.0);
    double tmp = atom->sqrt_weight * delta;
    chi2 += tmp * tmp;
  }

  int num_residues = chain.calc_group_num_residues(group_id);
  return num_residues * (chi2 / sum_weight);
}

void
TLSModelEngine::anisotropic_fit_segment(int istart, int iend, double *residual)
{
  int group_id = 1;
  chain.set_group_range(group_id, istart, iend);

  isotropic_fit_group(group_id);
  *residual = isotropic_group_residual(group_id);
}

