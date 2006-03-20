// Copyright 2006 by TLSMD Development Group (see AUTHORS file)
// This code is part of the TLSMD distribution and governed by
// its license.  Please see the LICENSE file that should have been
// included as part of this package.
#ifndef __STRUCTURE_H__
#define __STRUCTURE_H__

#include <string>
#include <vector>

namespace TLSMD {

class GroupID {
 public:
  GroupID() : group_id_(0) {}
  GroupID(int gid) : group_id_(gid) {}

 private:
  int group_id_;
};


class ChainSegment {
 public:
  ChainSegment(char* frag_id1_, char* frag_id2_) 
    : frag_id1(frag_id1_), frag_id2(frag_id2_) {}
  ChainSegment(const std::string& frag_id1_, const std::string& frag_id2_) 
    : frag_id1(frag_id1_), frag_id2(frag_id2_) {}
  ChainSegment(const ChainSegment& other) 
    : frag_id1(other.frag_id1), frag_id2(other.frag_id2) {}

  std::string frag_id1;
  std::string frag_id2;
};

class ChainSegmentSet {
 public:
  void add_segment(char* frag_id1, char* frag_id2) { 
    chain_segments.push_back(ChainSegment(frag_id1, frag_id2)); 
  }
  void add_segment(const std::string& frag_id1, const std::string& frag_id2) {
    chain_segments.push_back(ChainSegment(frag_id1, frag_id2));
  }
  void add_segment(const ChainSegment& chain_segment) {
    chain_segments.push_back(chain_segment);
  }

  std::vector<ChainSegment> chain_segments;
};



class Atom {
public:
  Atom();
  Atom(const Atom& other) 
    : name(other.name), frag_id(other.frag_id), x(other.x), y(other.y), z(other.z), 
    u_iso(other.u_iso), weight(other.weight), sqrt_weight(other.sqrt_weight) {
    for (int i = 0; i < 6; ++i) U[i] = other.U[i];
  }

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
  void set_group_range(int group_id, const std::string& frag_id1, const std::string& frag_id2);
  int calc_group_num_atoms(int group_id) const;
  int calc_group_num_residues(int group_id) const;
  void calc_group_centroid(int group_id, double *x, double *y, double *z) const;
  double calc_group_mean_uiso(int group_id) const;
  
  std::vector<Atom> atoms;

 private:
  
};

} // namespace TLSMD

#endif // __STRUCTURE_H__
