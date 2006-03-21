// Copyright 2006 by TLSMD Development Group (see AUTHORS file)
// This code is part of the TLSMD distribution and governed by
// its license.  Please see the LICENSE file that should have been
// included as part of this package.
#ifndef __STRUCTURE_H__
#define __STRUCTURE_H__

#include <iostream>
#include <string>
#include <vector>

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
    
    class Segment {
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
    };
    
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
      };
      
    public:
      SegmentSet(Chain* chain) : chain_(chain), segments_() {}

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
    }; 

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
  
  int CalcNumAtoms(Chain::SegmentSet::AtomIterator atom, Chain::SegmentSet::AtomIterator end);
  int CalcNumResidues(Chain::SegmentSet::AtomIterator atom, Chain::SegmentSet::AtomIterator end);
  void CalcCentroid(Chain::SegmentSet::AtomIterator atom, Chain::SegmentSet::AtomIterator end,
		    double *x, double *y, double *z);

} // namespace TLSMD

#endif // __STRUCTURE_H__
