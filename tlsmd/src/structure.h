// Copyright 2006 by TLSMD Development Group (see AUTHORS file)
// This code is part of the TLSMD distribution and governed by
// its license.  Please see the LICENSE file that should have been
// included as part of this package.
#ifndef __STRUCTURE_H__
#define __STRUCTURE_H__

#define NAME_LEN     8
#define FRAG_ID_LEN  8

class Atom {
public:
  bool is_mainchain();
  bool in_group(int gid) { return gid == group_id; }

  char    name[NAME_LEN];             /* atom name */
  char    frag_id[FRAG_ID_LEN];       /* fragment id (residue name) */
  int     ifrag;                      /* fragment index */
  double  x, y, z;
  double  xtls, ytls, ztls;
  double  u_iso;
  double  u_iso_tmp;
  double  U[6];
  double  Utmp[6];
  double  sqrt_weight;
  double  sqrt_weight_tmp;
  int     group_id;
};

class Chain {
public:
  Chain();
  ~Chain();

  void set_num_atoms(int na);
  void set_group_range(int group_id, int istart, int iend);
  int calc_group_num_atoms(int group_id);
  int calc_group_num_residues(int group_id);
  void calc_group_centroid(int group_id, double *x, double *y, double *z);
  
  Atom *atoms;
  int num_atoms;
};

#endif // __STRUCTURE_H__
