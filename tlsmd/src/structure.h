// Copyright 2006 by TLSMD Development Group (see AUTHORS file)
// This code is part of the TLSMD distribution and governed by
// its license.  Please see the LICENSE file that should have been
// included as part of this package.
#ifndef __STRUCTURE_H__
#define __STRUCTURE_H__

#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <iostream>
#include <sstream>

namespace TLSMD {

  class Atom {
  public:
    Atom();
    Atom(const Atom& other) 
      : name(other.name), frag_id(other.frag_id), x(other.x), y(other.y), z(other.z), 
      u_iso(other.u_iso), weight(other.weight), sqrt_weight(other.sqrt_weight) {
      for (int i = 0; i < 6; ++i) U[i] = other.U[i];
    }
    Atom& operator=(const Atom& other) {
      name = other.name;
      frag_id = other.frag_id;
      x = other.x;
      y = other.y;
      z = other.z;
      u_iso = other.u_iso;
      for (int i = 0; i < 6; ++i) U[i] = other.U[i];
      weight = other.weight;
      sqrt_weight = other.sqrt_weight;
      return(*this);
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

    class FragmentIDMap {
      // custom hash table for fragment ids
    public:
      FragmentIDMap(const Chain& chain);
      FragmentIDMap(const std::string& frag_id1, const std::string& frag_id2);
      std::vector<Atom>::iterator& operator[](const std::string& frag_id);
      const std::vector<Atom>::iterator& operator[](const std::string& frag_id) const;
      bool has_key(const std::string& frag_id) const;
    private:
      void reset(const std::string& frag_id1, const std::string& frag_id2);
      bool split(const std::string& frag_id, int& seq_num, char& icode) const;
      bool splitidx(const std::string& frag_id, int& seq_num_idx, int& icode_idx) const;
      struct SeqNumCell_ {
	std::vector<Atom>::iterator value_;
	std::vector<std::vector<Atom>::iterator> icode_table_;
      };
      int seq_num1_;
      std::vector<SeqNumCell_> seq_num_table_;
    }; // FragmentIDMap
    
    class Segment {
      // begin/end iterators defining a segment in std::vector<Atom>
    public:
      Segment(std::vector<Atom>::iterator first, std::vector<Atom>::iterator last) 
	: begin_(first), end_(last) {}
      Segment& operator=(const Segment& other) {
	  begin_ = other.begin_;
	  end_ = other.end_;
	  return(*this);
	}
      std::vector<Atom>::iterator begin() { return begin_; }
      std::vector<Atom>::iterator end() { return end_; }
    private:
      std::vector<Atom>::iterator begin_;
      std::vector<Atom>::iterator end_;
    }; // Segment
    
    class SegmentSet {
    public:

      class AtomIterator {
      public:
	AtomIterator() {}
	// begin constructor
	AtomIterator(std::vector<Segment>::iterator first, std::vector<Segment>::iterator last) 
	  : segment_iter_(first), segment_iter_end_(last) { atom_iter_ = segment_iter_->begin(); }
	// end constructor 
	  AtomIterator(std::vector<Segment>::iterator slast, std::vector<Atom>::iterator alast) 
	    : segment_iter_(slast), segment_iter_end_(slast), atom_iter_(alast) {}
	// copy constructor
	AtomIterator(const AtomIterator& other) 
	  : segment_iter_(other.segment_iter_), segment_iter_end_(other.segment_iter_end_), atom_iter_(other.atom_iter_) {}
	AtomIterator& operator=(const AtomIterator& other) {
	  segment_iter_ = other.segment_iter_;
	  segment_iter_end_ = other.segment_iter_end_;
	  atom_iter_ = other.atom_iter_;
	  return(*this);
	}
	bool operator!=(const AtomIterator& other) {
	  return(atom_iter_ != other.atom_iter_);
	}
	AtomIterator& operator++() {
	  ++atom_iter_;
	  if (atom_iter_ != segment_iter_->end()) {
	    return(*this);
	  }
	  ++segment_iter_;
	  if (segment_iter_ != segment_iter_end_) {
	    atom_iter_ = segment_iter_->begin();
	  }
	  return(*this);
	}
	AtomIterator operator++(int) {
	  AtomIterator tmp(*this);
	  ++(*this);
	  return(tmp);
	}
	Atom& operator*() {
	  return(*atom_iter_);
	}
	Atom* operator->() {
	  return(&(*atom_iter_));
	}
      private:
	std::vector<Segment>::iterator segment_iter_;
	std::vector<Segment>::iterator segment_iter_end_;
	std::vector<Atom>::iterator atom_iter_;
      }; // AtomIterator
      
    public:
      SegmentSet(Chain* chain);
      void add_segment(std::vector<Atom>::iterator first, std::vector<Atom>::iterator last) {
	segments_.push_back(Segment(first, last));
      }
      void add_segment(const std::string& frag_id1, const std::string& frag_id2);
      AtomIterator begin() { 
	return AtomIterator(segments_.begin(), segments_.end()); 
      }
      AtomIterator end() {
	Segment& segment = segments_.back();
	return AtomIterator(segments_.end(), segment.end());
      }
    private:
      SegmentSet(const SegmentSet&);
      Chain* chain_;
      std::vector<Segment> segments_;
    }; // SegmentSet 

  public:
    Chain();
    ~Chain();

    void set_num_atoms(int na);
    void map_frag_ids();
    bool has_frag_id(const std::string& frag_id) const;
    const std::vector<Atom>::iterator& frag_id_begin(const std::string& frag_id) const;
    const std::vector<Atom>::iterator& frag_id_end(const std::string& frag_id) const;

    std::vector<Atom> atoms;

  private:
    FragmentIDMap* frag_id_begin_map_;
    FragmentIDMap* frag_id_end_map_;
  };
  
  int CalcNumAtoms(Chain::SegmentSet::AtomIterator atom, Chain::SegmentSet::AtomIterator end);
  int CalcNumResidues(Chain::SegmentSet::AtomIterator atom, Chain::SegmentSet::AtomIterator end);
  void CalcCentroid(Chain::SegmentSet::AtomIterator atom, Chain::SegmentSet::AtomIterator end,
		    double *x, double *y, double *z);

} // namespace TLSMD

#endif // __STRUCTURE_H__
