// Copyright 2006 by TLSMD Development Group (see AUTHORS file)
// This code is part of the TLSMD distribution and governed by
// its license.  Please see the LICENSE file that should have been
// included as part of this package.
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <sstream>

#include "structure.h"

namespace TLSMD {

inline bool
frag_id_le(const std::string& frag_id1, const std::string& frag_id2) {
  std::istringstream ifrag_id1(frag_id1);
  std::istringstream ifrag_id2(frag_id2);
  int seq_num1, seq_num2;
  
  ifrag_id1 >> seq_num1;
  ifrag_id2 >> seq_num2;

  if (seq_num1 < seq_num2) return true;
  if (seq_num1 > seq_num2) return false;

  char icode1;
  bool icode1_exists = ifrag_id1.get(icode1);

  char icode2;
  bool icode2_exists = ifrag_id2.get(icode2);

  if (icode1_exists && icode2_exists) {
    return icode1 <= icode2;
  } else if (icode1_exists && !icode2_exists) {
    return false;
  }
  
  return true;
}

inline bool
frag_id_lt(const std::string& frag_id1, const std::string& frag_id2) {
  std::istringstream ifrag_id1(frag_id1);
  std::istringstream ifrag_id2(frag_id2);
  int seq_num1, seq_num2;
  
  ifrag_id1 >> seq_num1;
  ifrag_id2 >> seq_num2;

  if (seq_num1 < seq_num2) return true;
  if (seq_num2 > seq_num1) return false;

  char icode1;
  bool icode1_exists = ifrag_id1.get(icode1);

  char icode2;
  bool icode2_exists = ifrag_id2.get(icode2);

  if (icode1_exists && icode2_exists) {
    return icode1 < icode2;
  } else if (!icode1_exists && icode2_exists) {
    return true;
  }
  
  return false;
}


Atom::Atom() 
  : name(""),
    frag_id(""),
    ifrag(0), 
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


Chain::FragmentIDMap::FragmentIDMap(const Chain& chain) : seq_num_table_() {
  const Atom& atom1 = chain.atoms.front();
  const Atom& atom2 = chain.atoms.back();
  reset(atom1.frag_id, atom2.frag_id);
}

Chain::FragmentIDMap::FragmentIDMap(const std::string& frag_id1, const std::string& frag_id2) : seq_num_table_() {
  reset(frag_id1, frag_id2);
}

std::vector<Atom>::iterator& 
Chain::FragmentIDMap::operator[](const std::string& frag_id) {
  int seq_num_idx, icode_idx;
  bool has_icode = splitidx(frag_id, seq_num_idx, icode_idx);
  int sz = seq_num_table_.size();
  assert(seq_num_idx >= 0 && seq_num_idx <= sz);
  SeqNumCell_& cell = seq_num_table_[seq_num_idx];
  if (has_icode) {
    sz = cell.icode_table_.size();
    if (icode_idx >= sz) {
      std::cout << "ADDING ICODE " << seq_num_idx + seq_num1_ << ":" << icode_idx;
      cell.icode_table_.resize(icode_idx + 1);
    }
    return cell.icode_table_[icode_idx];
  }
  return cell.value_;
}

const std::vector<Atom>::iterator&
Chain::FragmentIDMap::operator[](const std::string& frag_id) const {
  int seq_num_idx, icode_idx;
  bool has_icode = splitidx(frag_id, seq_num_idx, icode_idx);
  int sz = seq_num_table_.size();
  assert(seq_num_idx >= 0 && seq_num_idx <= sz);
  const SeqNumCell_& cell = seq_num_table_[seq_num_idx];
  if (has_icode) {
    sz = cell.icode_table_.size();
    assert(icode_idx < sz);
    return cell.icode_table_[icode_idx];
  }
  return cell.value_;
}

bool
Chain::FragmentIDMap::has_key(const std::string& frag_id) const {
  std::vector<Atom>::iterator def;
  return (*this)[frag_id] != def;
}

void 
Chain::FragmentIDMap::reset(const std::string& frag_id1, const std::string& frag_id2) {
  int seq_num1, seq_num2;
  char icode1, icode2;
  split(frag_id1, seq_num1, icode1);
  split(frag_id2, seq_num2, icode2);
  seq_num1_ = seq_num1;
  seq_num_table_.resize(seq_num2 - seq_num1 + 1);
}

bool
Chain::FragmentIDMap::split(const std::string& frag_id, int& seq_num, char& icode) const {
  std::istringstream ifrag_id(frag_id);
  ifrag_id >> seq_num;
  return ifrag_id.get(icode);
}

bool
Chain::FragmentIDMap::splitidx(const std::string& frag_id, int& seq_num_idx, int& icode_idx) const {
  int seq_num;
  char icode;
  bool has_icode = split(frag_id, seq_num, icode);
  seq_num_idx = seq_num - seq_num1_;
  if (has_icode) {
    assert(icode > 64 && icode < 91);
    icode_idx = icode - 64;
  }
  return has_icode;
}

Chain::SegmentSet::SegmentSet(Chain* chain) : chain_(chain), segments_() {}

void
Chain::SegmentSet::add_segment(const std::string& frag_id1, const std::string& frag_id2) {
  assert(chain_->has_frag_id(frag_id1) && chain_->has_frag_id(frag_id2));
  add_segment(chain_->frag_id_begin(frag_id1), chain_->frag_id_end(frag_id2));
}

void
Chain::SegmentSet::add_segment_slow(const std::string& frag_id1, const std::string& frag_id2) {
  std::vector<Atom>::iterator atom;
  for (atom = chain_->atoms.begin(); atom != chain_->atoms.end() && frag_id_lt(atom->frag_id, frag_id1); ++atom);
  std::vector<Atom>::iterator first = atom;
  for (; atom != chain_->atoms.end() && frag_id_le(atom->frag_id, frag_id2); ++atom);
  std::vector<Atom>::iterator last = atom;
  add_segment(first, last);
}

Chain::Chain() : atoms(), frag_id_begin_map_(0), frag_id_end_map_(0) {}

Chain::~Chain() {
  delete frag_id_begin_map_;
  delete frag_id_end_map_;
}

void
Chain::set_num_atoms(int na) {
  atoms.resize(na);
}

void
Chain::map_frag_ids() {
  delete frag_id_begin_map_;
  delete frag_id_end_map_;

  frag_id_begin_map_ = new FragmentIDMap(*this);
  frag_id_end_map_ = new FragmentIDMap(*this);

  std::vector<Atom>::iterator atom = atoms.begin();
  const std::string* last_frag_id = &atom->frag_id;
  (*frag_id_begin_map_)[atom->frag_id] = atom;

  for (; atom != atoms.end(); ++atom) {
    if (last_frag_id->compare(atom->frag_id) != 0) {
      (*frag_id_begin_map_)[atom->frag_id] = atom;
      (*frag_id_end_map_)[*last_frag_id] = atom;
      last_frag_id = &atom->frag_id;
    }
  }
  (*frag_id_end_map_)[*last_frag_id] = atom;
}

int
CalcNumAtoms(Chain::SegmentSet::AtomIterator atom, Chain::SegmentSet::AtomIterator end) {
  int num_atoms = 0;
  for (; atom != end; ++atom) ++num_atoms;
  return num_atoms;
}

int
CalcNumResidues(Chain::SegmentSet::AtomIterator atom, Chain::SegmentSet::AtomIterator end) {
  int num_residues = 0;
  std::string const *cur_frag_id = 0;
  
  for (; atom != end; ++atom) {
    if (cur_frag_id == 0 || cur_frag_id->compare(atom->frag_id) != 0) {
      ++num_residues;
      cur_frag_id = &atom->frag_id;
    }
  }
  return num_residues;
}

void
CalcCentroid(Chain::SegmentSet::AtomIterator begin, 
	     Chain::SegmentSet::AtomIterator end,
	     double *x, double *y, double *z) {

  int natoms = 0;
  double cx = 0.0;
  double cy = 0.0;
  double cz = 0.0;

  Chain::SegmentSet::AtomIterator atom;
  for (atom = begin; atom != end; ++atom) {
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

} // TLSMD
