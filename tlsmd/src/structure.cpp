// Copyright 2006 by TLSMD Development Group (see AUTHORS file)
// This code is part of the TLSMD distribution and governed by
// its license.  Please see the LICENSE file that should have been
// included as part of this package.
#include <string.h>
#include <stdio.h>

#include "structure.h"

namespace TLSMD {

Atom::Atom() 
  : ifrag(0), 
    x(0.0), 
    y(0.0), 
    z(0.0), 
    u_iso(0.0), 
    weight(0.0),
    sqrt_weight(0.0), 
    group_id(0) {
  for (int i = 0; i < 6; ++i) U[i] = 0.0;
}

bool 
Atom::is_mainchain() const {
  if (name == "N"|| name == "CA" || name == "C" || name == "O" || name == "CB")  return true;
  return false;
}

Chain::Chain() {
}

Chain::~Chain() {
}

void
Chain::set_num_atoms(int na) {
  atoms.resize(na);
}

void
Chain::set_group_range(int group_id, int istart, int iend) {
  int ia = 0;
  std::vector<Atom>::iterator atom;
  for (atom = atoms.begin(); atom != atoms.end(); ++atom, ++ia) {
    if (ia >= istart && ia <= iend) {
      atom->set_group(group_id);
    } else {
      atom->set_group(0);
    }
  }
  calc_group_num_atoms(group_id);
}

int
Chain::calc_group_num_atoms(int group_id) const {
  int num_atoms = 0;
  std::vector<Atom>::const_iterator atom;
  for (atom = atoms.begin(); atom != atoms.end(); ++atom) {
    if (atom->in_group(group_id)) ++num_atoms;
  }
  return num_atoms;
}

int
Chain::calc_group_num_residues(int group_id) const {
  int num_residues = 0;
  std::string const *cur_frag_id = 0;
  
  std::vector<Atom>::const_iterator atom;
  for (atom = atoms.begin(); atom != atoms.end(); ++atom) {
    if (!atom->in_group(group_id)) continue;
    if (cur_frag_id == 0 || cur_frag_id->compare(atom->frag_id) != 0) {
      ++num_residues;
      cur_frag_id = &atom->frag_id;
    }
  }
  return num_residues;
}

void
Chain::calc_group_centroid(int group_id, double *x, double *y, double *z) const {
  int natoms = 0;
  double cx = 0.0;
  double cy = 0.0;
  double cz = 0.0;

  std::vector<Atom>::const_iterator atom;
  for (atom = atoms.begin(); atom != atoms.end(); ++atom) {
    if (!atom->in_group(group_id)) continue;
    ++natoms;
    cx += atom->x;
    cy += atom->y;
    cz += atom->z;
  }

  if (natoms > 0) {
    *x = cx / natoms;
    *y = cy / natoms;
    *z = cz / natoms;
  } else {
    *x = 0.0;
    *y = 0.0;
    *z = 0.0;
  }
}

double
Chain::calc_group_mean_uiso(int group_id) const {
  int natoms = 0;
  double sum_uiso = 0.0;

  std::vector<Atom>::const_iterator atom;
  for (atom = atoms.begin(); atom != atoms.end(); ++atom) {
    if (atom->in_group(group_id)) {
      ++natoms;
      sum_uiso += atom->u_iso;
    }
  }
  if (natoms == 0) return 0.0;
  return sum_uiso / natoms;
}

} // TLSMD
