#!/usr/bin/env python
## This program performs some basic tests in the mmLib.Structure object model

import sys
import time
import copy
import test_util
from mmLib.FileLoader import LoadStructure, SaveStructure, decode_format
from mmLib.Structure import *


class Stats(dict):
    def __init__(self):
        self["chain_count"]    = 0
        self["fragment_count"] = 0
        self["atom_count"]     = 0
        self["bond_count"]     = 0

    def print_stats(self):
        print "Number of Chains-----:",self["chain_count"]
        print "Number of Fragments--:",self["fragment_count"]
        print "Number of Atoms------:",self["atom_count"]
        print "Number of Bonds------:",self["bond_count"]


def bond_test(bond, atom, stats):
    """Tests the mmLib.Structure.Bond object
    """
    stats["bond_count"] += 1
    
    assert isinstance(bond, Bond)
    assert bond.atom1 != bond.atom2
    assert bond.atom1 == atom or bond.atom2 == atom
    assert bond.atom1.model_id == bond.atom2.model_id

    assert bond.atom1.alt_loc == bond.atom2.alt_loc or \
           (bond.atom1.alt_loc == "" and bond.atom2.alt_loc != "") or \
           (bond.atom1.alt_loc != "" and bond.atom2.alt_loc == "")

    assert bond.get_partner(bond.atom1) == bond.atom2
    assert bond.get_partner(bond.atom2) == bond.atom1


def atom_test(atom, stats):
    """Tests the mmLib.Structure.Atom object.
    """
    stats["atom_count"] += 1
    stats["testing"]     = atom

    len(atom)
    
    alt_loc = atom.fragment.chain.structure.default_alt_loc

    visited_atm_list = []
    for atm in atom.iter_alt_loc():
        assert isinstance(atm, Atom)
        assert atm in atom
        assert atm not in visited_atm_list
        visited_atm_list.append(atm)

        assert atm[atom.alt_loc] == atom

        assert atm.get_fragment() == atom.get_fragment()
        assert atm.get_chain() == atom.get_chain()
        assert atm.get_structure() == atom.get_structure()
        assert atm.name == atom.name
        assert atm.res_name == atom.res_name
        assert atm.fragment_id == atom.fragment_id
        assert atm.chain_id == atom.chain_id

    atom.calc_anisotropy()
    atom.calc_anisotropy3()
    #atom.calc_CCuij(atom)
    #atom.calc_Suij(atom)



def fragment_test(frag, stats):
    """Tests the mmLib.Structure.Fragment/AminoAcidResidue/NucleicAcidResidue
    objects.
    """
    stats["fragment_count"] += 1
    stats["testing"] = frag

    len(frag)

    alt_loc = frag.chain.structure.default_alt_loc

    visited_atm_list = []
    for atm in frag.iter_atoms():
        assert isinstance(atm, Atom)
        assert atm in frag
        assert atm.name in frag
        
        assert atm not in visited_atm_list
        visited_atm_list.append(atm)

        assert frag[atm.name] == atm
        assert frag.get_atom(atm.name) == atm
        assert atm.get_fragment() == frag
        assert atm.get_chain() == frag.get_chain()
        assert atm.get_structure() == frag.get_structure()
        assert atm.alt_loc == alt_loc or atm.alt_loc == ""

    ## test iter_bonds
    num_bonds = 0
    for bond in frag.iter_bonds():
        assert isinstance(bond, Bond)
        assert bond.atom1 in frag or bond.atom2 in frag
        num_bonds += 1
        
    f = frag.get_offset_fragment(-1)
    f = frag.get_offset_fragment(1)

    ## test fragment API for obvious errors
    if isinstance(frag, Residue):
        r = frag.get_offset_residue(-1)
        r = frag.get_offset_residue(1)

    if isinstance(frag, AminoAcidResidue):
        x = frag.calc_mainchain_bond_length()
        x = frag.calc_mainchain_bond_angle()
        x = frag.calc_torsion_psi()
        x = frag.calc_torsion_phi()
        x = frag.calc_torsion_omega()
        x = frag.is_cis()
        x = frag.calc_pucker_torsion()
        x = frag.calc_torsion_chi1()
        x = frag.calc_torsion_chi2()
        x = frag.calc_torsion_chi3()
        x = frag.calc_torsion_chi4()
        x = frag.calc_torsion_chi()

    stats["testing"] = None



def chain_test(chain, stats):
    """Tests the mmLib.Structure.Chain object.
    """
    stats["chain_count"] += 1
    stats["testing"]      = chain

    len(chain)

    for frag in chain.iter_fragments():
        assert isinstance(frag, Fragment)
        assert frag in chain
        assert frag.fragment_id in chain
        assert chain[frag.fragment_id] == frag
        assert chain[chain.index(frag)] == frag
        assert chain.get_fragment(frag.fragment_id) == frag
        assert frag.get_chain() == chain
        assert frag.get_structure() == chain.get_structure()

    for res in chain.iter_amino_acids():
        assert isinstance(res, AminoAcidResidue)
        assert res.get_chain() == chain
        assert res.get_structure() == chain.get_structure()

    for res in chain.iter_nucleic_acids():
        assert isinstance(res, NucleicAcidResidue)
        assert res.get_chain() == chain
        assert res.get_structure() == chain.get_structure()

    for atm in chain.iter_atoms():
        assert isinstance(atm, Atom)
        assert atm.get_chain() == chain
        assert atm.get_structure() == chain.get_structure()

    for bond in chain.iter_bonds():
        assert isinstance(bond, Bond)

    stats["testing"] = None



def struct_test(struct, stats):
    """Tests the mmLib.Structure.Structure object.
    """
    stats["testing"] = struct

    len(struct)

    ## get a lit of all alt_loc ids
    alt_loc_list = struct.alt_loc_list()
    print "alt_loc_list: ",alt_loc_list

    ## make sure the default alt_loc was used when constructing the
    ## structure
    for atm in struct.iter_atoms():
        assert atm.alt_loc == "" or atm.alt_loc == struct.default_alt_loc

    old_model = struct.model
    for model in struct.iter_models():
        struct.set_model(model)
        assert model in struct

        for chain in struct.iter_chains():
            assert isinstance(chain, Chain)
            assert chain in struct
            assert chain.chain_id in struct
            assert struct[chain.chain_id] == chain
            assert struct[struct.index(chain)] == chain
            assert struct.get_chain(chain.chain_id) == chain
            assert chain.get_structure() == struct 

        for frag in struct.iter_fragments():
            assert isinstance(frag, Fragment)
            assert frag.get_structure() == struct

        for res in struct.iter_amino_acids():
            assert isinstance(res, AminoAcidResidue)
            assert res.get_structure() == struct

        for res in struct.iter_nucleic_acids():
            assert isinstance(res, NucleicAcidResidue)
            assert res.get_structure() == struct

        for atm in struct.iter_atoms():
            assert isinstance(atm, Atom)
            assert atm.get_structure() == struct

        for atm in struct.iter_all_atoms():
            assert isinstance(atm, Atom)
            assert atm.get_structure() == struct

        for bond in struct.iter_bonds():
            assert isinstance(bond, Bond)

        old_alt_loc  = struct.default_alt_loc
        for alt_loc in alt_loc_list:
            struct.set_alt_loc(alt_loc)
            for atm in struct.iter_atoms():
                assert atm.alt_loc == "" or atm.alt_loc == alt_loc
        struct.set_alt_loc(old_alt_loc)

    struct.set_model(old_model)
    
    stats["testing"] = None



def run_structure_tests(struct, stats):
    """Run basic API tests on the mmLib.Structure object and print statistics.
    """
    struct_test(struct, stats)

    for chain in struct.iter_chains():
        chain_test(chain, stats)

    for frag in struct.iter_fragments():
        fragment_test(frag, stats)

    for atm in struct.iter_all_atoms():
        atom_test(atm, stats)

    bond_dict = {}
    for atm in struct.iter_all_atoms():
        for bond in atm.iter_bonds():
            if not bond_dict.has_key(bond):
                bond_dict[bond] = True
                bond_test(bond, atm, stats)


def copy_verify(struct, stats, struct_cp, stats_cp):
    assert stats["chain_count"] == stats_cp["chain_count"]
    assert stats["fragment_count"] == stats_cp["fragment_count"]
    assert stats["atom_count"] == stats_cp["atom_count"]
    assert stats["bond_count"] == stats_cp["bond_count"]
    

def file_verify(path, struct, stats):
    """Use some independent parsers to verify some simple stats between the
    structure and the file description.
    """
    if decode_format(path) == "PDB":
        fil_stats = test_util.pdb_stats(path)

    elif decode_format(path) == "CIF":
        fil_stats = test_util.cif_stats(path)

    print "File Atom Count------:",fil_stats["atoms"]
    print "Struct Atom Count----:",stats["atom_count"]

    assert fil_stats["atoms"] == stats["atom_count"]



def save_verify(struct, stats):
    """Save structure in all supported formats, then reload it and
    compare structures.
    """
    ## pdb
    print "[temp.pdb]"
    SaveStructure(fil="temp.pdb", struct=struct, format="PDB")
    fil_stats = test_util.pdb_stats("temp.pdb")

    struct_stats = Stats()
    pdb_struct = LoadStructure(fil="temp.pdb")
    run_structure_tests(pdb_struct, struct_stats)

    struct_stats.print_stats()
    
    ## mmCIF
    print "[temp.cif]"
    SaveStructure(fil="temp.cif", struct=struct, format="CIF")
    fil_stats = test_util.cif_stats("temp.cif")

    struct_stats = Stats()
    cif_struct = LoadStructure(fil="temp.cif")
    run_structure_tests(cif_struct, struct_stats)

    struct_stats.print_stats()
    
    


def main(walk_path, start_path):
    print "Running Python Macromolecular Library Test Program"
    print "--------------------------------------------------"
    print "This program will throw a AssertionError (or worse!) if it"
    print "runs into any problems."
    print

    for path in test_util.walk_pdb_cif(walk_path, start_path):
        print "[%s]" % (path)
        time1 = time.time()
        
        struct = LoadStructure(fil = path)        
        stats = Stats()

        ## test the mmLib.Structure object API and
        ## with massive sanity checking
        print "[loaded struct]"
        try:
            run_structure_tests(struct, stats)
        except AssertionError:
            print "*** AssertionError while testing: %s ***" % (
                str(stats["testing"]))
            raise
        except:
            print "*** Error while testing: %s ***" % (
                str(stats["testing"]))
            raise

        stats.print_stats()


        ## copy the structure and re-run those tests
        print "[copy struct]"
        struct_cp = copy.deepcopy(struct)
        stats_cp  = Stats()

        try:
            run_structure_tests(struct_cp, stats_cp)
        except AssertionError:
            print "*** AssertionError while testing: %s ***" % (
                str(stats_cp["testing"]))
            raise
        except:
            print "*** Error while testing: %s ***" % (
                str(stats_cp["testing"]))
            raise

        stats_cp.print_stats()
        copy_verify(struct, stats, struct_cp, stats_cp)
        

        ## verify the number of atoms in the mmLib.Structure object
        ## matches the number of atoms in the source file 
        print "[file verify]"
        file_verify(path, struct, stats)

        ## test file saving
        print "[save verify]"
        save_verify(struct, stats)

        time2 = time.time()
        print "Tests Time (sec)-----:",int(time2-time1)

def usage():
    print "usage: mmlib_test.py <PDB/mmCIF file or directory of files>"
    sys.exit(1)

if __name__ == "__main__":
    try:
        path = sys.argv[1]
    except IndexError:
        usage()

    try:
        start_path = sys.argv[2]
    except IndexError:
        start_path = None
    
    main(path, start_path)

