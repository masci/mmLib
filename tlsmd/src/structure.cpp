// Copyright 2006 by TLSMD Development Group (see AUTHORS file)
// This code is part of the TLSMD distribution and governed by
// its license.  Please see the LICENSE file that should have been
// included as part of this package.
#include <string.h>
#include <stdio.h>

#include "structure.h"

bool 
Atom::is_mainchain() {
  if (strcmp(name, "N")==0)  return true;
  if (strcmp(name, "CA")==0) return true;
  if (strcmp(name, "C")==0)  return true;
  if (strcmp(name, "O")==0)  return true;
  if (strcmp(name, "CB")==0) return true;
  return false;
}

Chain::Chain() : atoms(0), num_atoms(0) {
}

Chain::~Chain() {
  delete[] atoms;
}

void
Chain::set_num_atoms(int na) {
  if (atoms != 0) {
    delete[] atoms;
    atoms = 0;
  }
  num_atoms = na;
  atoms = new Atom[num_atoms];
}

void
Chain::set_group_range(int group_id, int istart, int iend) {
  Atom *atom = atoms;
  for (int ia = 0; ia < num_atoms; ++ia, ++atom) {
    if (ia >= istart && ia <= iend) {
      atom->group_id = group_id;
    } else {
      atom->group_id = 0;
    }
  }
  calc_group_num_atoms(group_id);
}

int
Chain::calc_group_num_atoms(int group_id) {
  int natoms = 0;
  Atom *atom = atoms;
  for (int ia = 0; ia < num_atoms; ++ia, ++atom) {
    if (atom->in_group(group_id)) ++natoms;
  }
  //printf("calc_group_num_atoms %d\n", natoms);
  return natoms;
}

int
Chain::calc_group_num_residues(int group_id) {
  int num_residues = 0;
  char *cur_frag_id = 0;

 Atom *atom = atoms;
 for (int ia = 0; ia < num_atoms; ++ia, ++atom) {
   if (!atom->in_group(group_id)) continue;
   if (cur_frag_id == 0 || strcmp(atom->frag_id, cur_frag_id) != 0) {
     ++num_residues;
     cur_frag_id = atom->frag_id;
   }
 }
 //printf("calc_group_num_residues %d\n", num_residues);
 return num_residues;
}

void
Chain::calc_group_centroid(int group_id, double *x, double *y, double *z) {
  int n = 0;
  double cx = 0.0;
  double cy = 0.0;
  double cz = 0.0;

  Atom *atom = atoms;
  for (int ia = 0; ia < num_atoms; ++ia, ++atom) {
    if (!atom->in_group(group_id)) continue;
    ++n;
    cx += atom->x;
    cy += atom->y;
    cz += atom->z;
  }

  if (n > 0) {
    *x = cx / n;
    *y = cy / n;
    *z = cz / n;
  } else {
    *x = 0.0;
    *y = 0.0;
    *z = 0.0;
  }
}
