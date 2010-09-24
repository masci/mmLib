#!/usr/bin/env python
## Copyright 2002-2010 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""This program performs some basic tests in the mmLib.Structure object model.
"""

import sys
import time
import copy
import weakref
import gc

try:
    import numpy
except ImportError:
    import mmLib.NumericCompat as numpy

import test_util
from mmLib.FileIO import get_file_extension
from mmLib import Structure, FileIO


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

    assert isinstance(bond, Structure.Bond)
    assert bond.atom1 != bond.atom2
    assert bond.atom1 == atom or bond.atom2 == atom
    assert bond.atom1.model_id == bond.atom2.model_id

    assert bond.atom1.alt_loc == bond.atom2.alt_loc or \
           (bond.atom1.alt_loc == "" and bond.atom2.alt_loc != "") or \
           (bond.atom1.alt_loc != "" and bond.atom2.alt_loc == "")

    assert bond.get_partner(bond.atom1) == bond.atom2
    assert bond.get_partner(bond.atom2) == bond.atom1

    assert bond.get_atom1() == bond.atom1
    assert bond.get_atom2() == bond.atom2

    assert bond.get_fragment1() == bond.atom1.fragment
    assert bond.get_fragment2() == bond.atom2.fragment

    assert bond.get_chain1() == bond.atom1.fragment.chain
    assert bond.get_chain2() == bond.atom2.fragment.chain

    assert bond.get_model1() == bond.atom1.fragment.chain.model
    assert bond.get_model2() == bond.atom2.fragment.chain.model

    assert bond.get_structure1() == bond.atom1.fragment.chain.model.structure
    assert bond.get_structure2() == bond.atom2.fragment.chain.model.structure


def atom_test(atom, stats):
    """Tests the mmLib.Structure.Atom object.
    """
    stats["atom_count"] += 1
    stats["testing"]     = atom

    len(atom)

    alt_loc = atom.get_structure().get_default_alt_loc()

    atom.get_fragment()
    atom.get_chain()
    atom.get_model()
    atom.get_structure()

    visited_atm_list = []
    for atm in atom.iter_alt_loc():
        assert isinstance(atm, Structure.Atom)
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

    try:
        atom.calc_anisotropy()
    except ZeroDivisionError:
        pass

    try:
        atom.calc_anisotropy3()
    except ZeroDivisionError:
        pass

def fragment_test(frag, stats):
    """Tests the mmLib.Structure.Fragment/AminoAcidResidue/NucleicAcidResidue
    objects.
    """
    stats["fragment_count"] += 1
    stats["testing"] = frag

    len(frag)

    frag.get_chain()
    frag.get_model()
    frag.get_structure()

    alt_loc = frag.get_structure().default_alt_loc

    ## test iter_atoms
    visited_atm_list = []
    for atm in frag.iter_atoms():
        assert isinstance(atm, Structure.Atom)
        assert atm in frag
        assert atm.name in frag

        assert atm not in visited_atm_list
        visited_atm_list.append(atm)

        assert frag[atm.name] == atm
        assert frag.get_atom(atm.name) == atm
        assert atm.get_fragment() == frag
        assert atm.get_chain() == frag.get_chain()
        assert atm.get_model() == frag.get_model()
        assert atm.get_structure() == frag.get_structure()

        assert atm.alt_loc == alt_loc or atm.alt_loc == ""

    ## test iter_all_atoms
    for atm in frag.iter_all_atoms():
        pass

    ## test iter_bonds
    num_bonds = 0
    for bond in frag.iter_bonds():
        assert isinstance(bond, Structure.Bond)
        assert bond.atom1 in frag or bond.atom2 in frag
        num_bonds += 1

    f = frag.get_offset_fragment(-1)
    f = frag.get_offset_fragment(1)

    ## test fragment API for obvious errors
    if isinstance(frag, Structure.Residue):
        r = frag.get_offset_residue(-1)
        r = frag.get_offset_residue(1)

    if isinstance(frag, Structure.AminoAcidResidue):
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

    chain.get_model()
    chain.get_structure()

    for frag in chain.iter_fragments():
        assert isinstance(frag, Structure.Fragment)
        assert frag in chain
        assert frag.fragment_id in chain
        assert chain[frag.fragment_id] == frag
        assert chain[chain.index(frag)] == frag
        assert chain.get_fragment(frag.fragment_id) == frag

        assert frag.get_chain() == chain
        assert frag.get_model() == chain.get_model()
        assert frag.get_structure() == chain.get_structure()

    chain.has_amino_acids()
    chain.count_amino_acids()
    for frag in chain.iter_amino_acids():
        assert isinstance(frag, Structure.AminoAcidResidue)

    chain.has_nucleic_acids()
    chain.count_nucleic_acids()
    for frag in chain.iter_nucleic_acids():
        assert isinstance(frag, Structure.NucleicAcidResidue)

    chain.has_standard_residues()
    chain.count_standard_residues()
    for frag in chain.iter_standard_residues():
        assert isinstance(frag, Structure.AminoAcidResidue) or \
               isinstance(frag, Structure.NucleicAcidResidue)

    chain.has_non_standard_residues()
    chain.count_non_standard_residues()
    for frag in chain.iter_non_standard_residues():
        assert not isinstance(frag, Structure.AminoAcidResidue) and \
               not isinstance(frag, Structure.NucleicAcidResidue)

    chain.has_waters()
    chain.count_waters()
    for frag in chain.iter_waters():
        assert frag.is_water()

    for atm in chain.iter_atoms():
        assert isinstance(atm, Structure.Atom)
        assert atm.get_chain() == chain
        assert atm.get_model() == chain.get_model()
        assert atm.get_structure() == chain.get_structure()

    for atm in chain.iter_all_atoms():
        assert isinstance(atm, Structure.Atom)
        assert atm.get_chain() == chain
        assert atm.get_model() == chain.get_model()
        assert atm.get_structure() == chain.get_structure()

    for bond in chain.iter_bonds():
        assert isinstance(bond, Structure.Bond)

    stats["testing"] = None


def model_test(model, stats):
    stats["testing"] = model

    model.get_structure()

    for chain in model.iter_chains():
        assert isinstance(chain, Structure.Chain)
        assert chain.get_model() == model
        assert chain.get_structure() == model.get_structure()

    for frag in model.iter_fragments():
        assert isinstance(frag, Structure.Fragment)
        assert frag.get_model() == model
        assert frag.get_structure() == model.get_structure()

    for atm in model.iter_atoms():
        assert isinstance(atm, Structure.Atom)
        assert atm.get_model() == model
        assert atm.get_structure() == model.get_structure()

    for atm in model.iter_all_atoms():
        assert isinstance(atm, Structure.Atom)
        assert atm.get_model() == model
        assert atm.get_structure() == model.get_structure()

    model.has_amino_acids()
    model.count_amino_acids()
    for frag in model.iter_amino_acids():
        assert isinstance(frag, Structure.AminoAcidResidue)

    model.has_nucleic_acids()
    model.count_nucleic_acids()
    for frag in model.iter_nucleic_acids():
        assert isinstance(frag, Structure.NucleicAcidResidue)

    model.has_standard_residues()
    model.count_standard_residues()
    for frag in model.iter_standard_residues():
        assert isinstance(frag, Structure.AminoAcidResidue) or \
               isinstance(frag, Structure.NucleicAcidResidue)

    model.has_non_standard_residues()
    model.count_non_standard_residues()
    for frag in model.iter_non_standard_residues():
        assert not isinstance(frag, Structure.AminoAcidResidue) and \
               not isinstance(frag, Structure.NucleicAcidResidue)

    model.has_waters()
    model.count_waters()
    for frag in model.iter_waters():
        assert frag.is_water()

    model.count_chains()
    model.count_fragments()
    model.count_atoms()
    model.count_all_atoms()
    model.count_bonds()

    ## test AlphaHelix
    for helix in model.iter_alpha_helicies():
        assert isinstance(helix, Structure.AlphaHelix)
        assert helix.get_model() == model
        assert helix.get_structure() == model.get_structure()

        for frag in helix.iter_fragments():
            assert isinstance(frag, Structure.Fragment)

        for atm in helix.iter_atoms():
            assert isinstance(atm, Structure.Atom)

        for atm in helix.iter_all_atoms():
            assert isinstance(atm, Structure.Atom)

    ## test BetaSheet
    for sheet in model.iter_beta_sheets():
        assert isinstance(sheet, Structure.BetaSheet)
        assert sheet.get_model() == model
        assert sheet.get_structure() == model.get_structure()

        for frag in sheet.iter_fragments():
            assert isinstance(frag, Structure.Fragment)

        for atm in sheet.iter_atoms():
            assert isinstance(atm, Structure.Atom)

        for atm in sheet.iter_all_atoms():
            assert isinstance(atm, Structure.Atom)

        for strand in sheet.iter_strands():
            assert isinstance(strand, Structure.Strand)
            assert strand.get_model() == model
            assert strand.get_structure() == model.get_structure()

            for frag in strand.iter_fragments():
                assert isinstance(frag, Structure.Fragment)

            for atm in strand.iter_atoms():
                assert isinstance(atm, Structure.Atom)

            for atm in strand.iter_all_atoms():
                assert isinstance(atm, Structure.Atom)

    ## test Site
    for site in model.iter_sites():
        assert isinstance(site, Structure.Site)
        assert site.get_model() == model
        assert site.get_structure() == model.get_structure()

        for frag in site.iter_fragments():
            assert isinstance(frag, Structure.Fragment)

        for atm in site.iter_atoms():
            assert isinstance(atm, Structure.Atom)

        for atm in site.iter_all_atoms():
            assert isinstance(atm, Structure.Atom)


def struct_test(struct, stats):
    """Tests the mmLib.Structure.Structure object.
    """
    stats["testing"] = struct

    ## get a lit of all alt_loc ids
    alt_loc_list = struct.alt_loc_list()
    print "alt_loc_list: ",alt_loc_list

    ## iterate over all atoms
    for atm in struct.iter_all_atoms():
        assert isinstance(atm, Structure.Atom)
        assert atm.get_structure() == struct

    ## make sure the default alt_loc was used when constructing the
    ## structure
    for atm in struct.iter_atoms():
        assert atm.alt_loc == "" or atm.alt_loc == struct.default_alt_loc

    for frag in struct.iter_all_fragments():
        assert isinstance(frag,Structure. Fragment)        

    struct.has_amino_acids()
    struct.count_amino_acids()
    for frag in struct.iter_amino_acids():
        assert isinstance(frag, Structure.AminoAcidResidue)

    struct.has_nucleic_acids()
    struct.count_nucleic_acids()
    for frag in struct.iter_nucleic_acids():
        assert isinstance(frag, Structure.NucleicAcidResidue)

    struct.has_standard_residues()
    struct.count_standard_residues()
    for frag in struct.iter_standard_residues():
        assert isinstance(frag, Structure.AminoAcidResidue) or \
               isinstance(frag, Structure.NucleicAcidResidue)

    struct.has_non_standard_residues()
    struct.count_non_standard_residues()
    for frag in struct.iter_non_standard_residues():
        assert not isinstance(frag, Structure.AminoAcidResidue) and \
               not isinstance(frag, Structure.NucleicAcidResidue)

    struct.has_waters()
    struct.count_waters()
    for frag in struct.iter_waters():
        assert frag.is_water()

    old_model = struct.get_default_model()

    for model in struct.iter_models():
        struct.set_default_model(model)

        for chain in struct.iter_chains():
            assert isinstance(chain, Structure.Chain)
            assert chain in struct
            assert chain.chain_id in struct
            assert struct[chain.chain_id] == chain
            assert struct[struct.index(chain)] == chain
            assert struct.get_chain(chain.chain_id) == chain
            assert chain.get_structure() == struct 

        for frag in struct.iter_fragments():
            assert isinstance(frag, Structure.Fragment)
            assert frag.get_structure() == struct

        for res in struct.iter_amino_acids():
            assert isinstance(res, Structure.AminoAcidResidue)
            assert res.get_structure() == struct

        for res in struct.iter_nucleic_acids():
            assert isinstance(res, Structure.NucleicAcidResidue)
            assert res.get_structure() == struct

        for atm in struct.iter_atoms():
            assert isinstance(atm, Structure.Atom)
            assert atm.get_structure() == struct

        for atm in struct.iter_all_atoms():
            assert isinstance(atm, Structure.Atom)
            assert atm.get_structure() == struct

        for bond in struct.iter_bonds():
            assert isinstance(bond, Structure.Bond)

        ## check default_alt_loc setting
        old_alt_loc  = struct.default_alt_loc

        for alt_loc in alt_loc_list:
            struct.set_default_alt_loc(alt_loc)
            for atm in struct.iter_atoms():
                assert atm.alt_loc == "" or atm.alt_loc == alt_loc

        struct.set_default_alt_loc(old_alt_loc)

        struct.count_models()
        struct.count_chains()
        struct.count_fragments()
        struct.count_amino_acids()
        struct.count_nucleic_acids()
        struct.count_atoms()
        struct.count_all_atoms()
        struct.count_bonds()

    struct.set_model(old_model)

    stats["testing"] = None


def run_structure_tests(struct, stats):
    """Run basic API tests on the mmLib.Structure object and print statistics.
    """
    struct_test(struct, stats)

    for model in struct.iter_models():
        model_test(model, stats)

    for chain in struct.iter_all_chains():
        chain_test(chain, stats)

    for frag in struct.iter_all_fragments():
        fragment_test(frag, stats)

    for atm in struct.iter_all_atoms():
        atom_test(atm, stats)

    bond_dict = {}
    for atm in struct.iter_all_atoms():
        for bond in atm.iter_bonds():
            if not bond_dict.has_key(bond):
                bond_dict[bond] = True
                bond_test(bond, atm, stats)

    ## crazy test: remove all atoms, then add them back
    print "ATOM REMOVEAL TEST"
    print "STRUCT NUM ATOMS %d" % (struct.count_all_atoms())
    alist = []
    for atm in struct.iter_all_atoms():
        alist.append(atm)
    for atm in alist:
        struct.remove_atom(atm)
    print "AFTER REMOVAL NUM ATOMS %d" % (struct.count_all_atoms())
    for atm in alist:
        struct.add_atom(atm)
    print "STRUCT NUM ATOMS(AFTERADD) %d" % (struct.count_all_atoms())


def cmp_struct(struct1, struct2):
    """Compare two Structure objects.
    """
    assert struct1.count_models()    == struct2.count_models()
    assert struct1.count_chains()    == struct2.count_chains()
    assert struct1.count_fragments() == struct2.count_fragments()
    assert struct1.count_atoms()     == struct2.count_atoms()
    assert struct1.count_all_atoms() == struct2.count_all_atoms()

    def cmp_atoms(atm1, atm2):
        try:
            assert atm1.element == atm2.element
        except AssertionError:
            print "%s != %s" % (atm1.element, atm2.element)
            raise

        assert atm1.name        == atm2.name
        assert atm1.fragment_id == atm2.fragment_id
        assert atm1.res_name    == atm2.res_name
        assert atm1.chain_id    == atm2.chain_id
        assert atm1.model_id    == atm2.model_id
        assert atm1.alt_loc     == atm2.alt_loc

        assert (atm1.position == None and atm2.position == None) or \
               numpy.allclose(atm1.position, atm2.position)

        assert (atm1.sig_position == None and atm2.sig_position == None) or \
               numpy.allclose(atm1.sig_position, atm2.sig_position)

        assert atm1.occupancy == atm2.occupancy or \
               numpy.allclose(atm1.occupancy, atm2.occupancy)

        assert atm1.sig_occupancy == atm2.sig_occupancy or \
               numpy.allclose(atm1.sig_occupancy, atm2.sig_occupancy)

        assert atm1.temp_factor == atm2.temp_factor or \
               numpy.allclose(atm1.temp_factor, atm2.temp_factor)

        assert atm1.sig_temp_factor == atm2.sig_temp_factor or \
               numpy.allclose(atm1.sig_temp_factor, atm2.sig_temp_factor)

        assert atm1.charge == atm2.charge or \
               numpy.allclose(atm1.charge, atm2.charge)

        U1 = atm1.get_U()
        U2 = atm2.get_U()
        assert (U1 == None and U2 == None) or numpy.allclose(U1, U2)

    for atm1 in struct1.iter_all_atoms():
        model2 = struct2.get_model(atm1.model_id)
        chain2 = model2.get_chain(atm1.chain_id)
        frag2  = chain2.get_fragment(atm1.fragment_id)
        atm2   = frag2.get_atom(atm1.name, atm1.alt_loc)

        try:
            cmp_atoms(atm1, atm2)
        except (AssertionError, TypeError):
            print
            print "ERROR: cmp_atom(%s, %s)" % (atm1, atm2)
            print
            print "atm1.position = %s" % (atm1.position)
            print "atm2.position = %s" % (atm2.position)
            print
            print "atm1.temp_factor = %s" % (atm1.temp_factor)
            print "atm2.temp_factor = %s" % (atm2.temp_factor)
            raise

        for bond in atm1.iter_bonds():
            atm1p = bond.get_partner(atm1)


def file_verify(path, struct, stats):
    """Use some independent parsers to verify some simple stats between the
    structure and the file description.
    """
    print "[file verify %s]" % path

    if get_file_extension(path) == "PDB":
        fil_stats = test_util.pdb_stats(path)

    elif get_file_extension(path) == "CIF":
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
    FileIO.SaveStructure(fil = "temp.pdb", struct = struct, format = "PDB")
    pdb_struct = FileIO.LoadStructure(
        fil           = "temp.pdb",
        library_bonds = True);
    cmp_struct(struct, pdb_struct)

    ## mmCIF
    print "[temp.cif]"
    FileIO.SaveStructure(fil = "temp.cif", struct = struct, format = "CIF")

    cif_struct = FileIO.LoadStructure(
        fil           = "temp.cif",
        library_bonds = True)
    cmp_struct(struct, cif_struct)


WEAKREF_LIST = []
WEAKREF_PATH = {}
def weakref_callback(ref):
    print "STRUCTURE DESTROYED: ",WEAKREF_PATH[ref]
    del WEAKREF_PATH[ref]
    WEAKREF_LIST.remove(ref)


def main(walk_path, start_path):
    print "Running Python Macromolecular Library Test Program"
    print "--------------------------------------------------"
    print "This program will throw a AssertionError (or worse!) if it"
    print "runs into any problems."
    print

    for path in test_util.walk_pdb_cif(walk_path, start_path):
        print "[%s]" % (path)
        time1 = time.time()

        stats  = Stats()
        struct = FileIO.LoadStructure(
            fil              = path,
            build_properties = ("library_bonds",))   

        ## track memory leaks
        sref = weakref.ref(struct, weakref_callback)
        WEAKREF_PATH[sref] = path
        WEAKREF_LIST.append(sref)

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
        cmp_struct(struct, struct_cp)

        ## verify the number of atoms in the mmLib.Structure object
        ## matches the number of atoms in the source file 
        file_verify(path, struct, stats)

        ## test file saving
        from mmLib.CIF import DataBlock
        if (hasattr(struct, "cif_data")
                and isinstance(struct.cif_data, DataBlock)):
            print "[save verify skipped]"
        else:
            print "[save verify]"
            save_verify(struct, stats)

        time2 = time.time()
        print "Tests Time (sec)-----:",int(time2-time1)

        ## dereference!
        struct    = None
        struct_cp = None
        stats     = None
        stats_cp  = None

        ## force garbage collection
        gc.collect()
        print


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
