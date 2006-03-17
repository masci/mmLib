// Copyright 2006 by TLSMD Development Group (see AUTHORS file)
// This code is part of the TLSMD distribution and governed by
// its license.  Please see the LICENSE file that should have been
// included as part of this package.
#ifndef __STRUCTURE_H__
#define __STRUCTURE_H__

#include <string>
#include <vector>

namespace TLSMD {

class Atom {
public:
  Atom();

  bool is_mainchain() const;
  bool in_group(int gid) const { return gid == group_id; }
  void set_group(int gid) { group_id = gid; }

  std::string name;
  std::string frag_id;
  int ifrag;
  double x, y, z;
  double u_iso;
  double U[6];
  double weight;
  double sqrt_weight;

 private:
  int group_id;
};

class Chain {
public:
  Chain();
  ~Chain();

  void set_num_atoms(int na);
  void set_group_range(int group_id, int istart, int iend);
  int calc_group_num_atoms(int group_id) const;
  int calc_group_num_residues(int group_id) const;
  void calc_group_centroid(int group_id, double *x, double *y, double *z) const;
  double calc_group_mean_uiso(int group_id) const;
  
  std::vector<Atom> atoms;
};

} // namespace TLSMD

#endif // __STRUCTURE_H__
