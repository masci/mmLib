#!/usr/bin/env python
## This program performs some basic tests in the mmLib.Structure object model

import sys
from mmLib.FileLoader import LoadStructure
from mmLib.Structure import *

def atom_test(atom):
    model = atom.fragment.chain.structure.default_model
    alt_loc = atom.fragment.chain.structure.default_alt_loc

    print "Atom Test",atom

    visited_atm_list = []
    for atm in atom.iter_alt_loc():
        assert isinstance(atm, Atom)
        assert atm in atom
        assert atm not in visited_atm_list
        visited_atm_list.append(atm)
        assert atom[atm.alt_loc, atm.model] == atm
        assert atom[atm.alt_loc].model == model
        assert atom[atm.model].alt_loc == alt_loc or \
               atom[atm.model].alt_loc == ""
        assert atm.get_fragment() == atom.get_fragment()
        assert atm.get_chain() == atom.get_chain()
        assert atm.get_structure() == atom.get_structure()
        assert atm.name == atom.name
        assert atm.res_name == atom.res_name
        assert atm.fragment_id == atom.fragment_id
        assert atm.chain_id == atom.chain_id

    for bond in atom.iter_bonds():
        assert isinstance(bond, Bond)
        assert bond.atom1 == atom or bond.atom2 == atom
        assert bond.atom1.model == bond.atom2.model
        assert bond.atom1.alt_loc == bond.atom2.alt_loc or \
               (bond.atom1.alt_loc == "" and bond.atom2.alt_loc != "") or \
               (bond.atom1.alt_loc != "" and bond.atom2.alt_loc == "")

    a = atom.calc_anisotropy()

    print "Atom Test Complete"

def fragment_test(frag):
    print "Fragment Test",frag

    model = frag.chain.structure.default_model
    alt_loc = frag.chain.structure.default_alt_loc

    visited_atm_list = []
    for atm in frag.iter_atoms():
        assert isinstance(atm, Atom)
        assert atm not in visited_atm_list
        visited_atm_list.append(atm)
        assert atm in frag
        try:
            assert frag[atm.name] == atm
        except AssertionError:
            print "[ERROR::fragment_test] atm=%s frag[atm.name=" % (
                atm, frag[atm.name])
        assert frag.get_atom(atm.name) == atm
        assert atm.get_fragment() == frag
        assert atm.get_chain() == frag.get_chain()
        assert atm.get_structure() == frag.get_structure()
        assert atm.alt_loc == alt_loc or atm.alt_loc == ""
        assert atm.model == model

    ## test iter_bonds
    num_bonds = 0
    for bond in frag.iter_bonds():
        assert isinstance(bond, Bond)
        assert bond.atom1 in frag or bond.atom2 in frag
        num_bonds += 1

    ## simple test to check that the bonds are being created properly
    ## amino acids should have at least one bond per atom, usually more
   # if isinstance(frag, AminoAcidResidue):
   #     assert num_bonds >= len(frag)-1
        
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

    print "Fragment Test Complete"

def chain_test(chain):
    print "Chain Test",chain

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

    print "Chain Test Complete"

def struct_test(struct):
    print "Structure Test"
    print struct
    
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

    for bond in struct.iter_bonds():
        assert isinstance(bond, Bond)

    print "Structure Test Complete"

def main(path):
    print "Testing File: ",path

    struct = LoadStructure(path)
    struct_test(struct)

    for chain in struct.iter_chains():
        chain_test(chain)

    for frag in struct.iter_fragments():
        fragment_test(frag)

    for atm in struct.iter_atoms():
        atom_test(atm)

if __name__ == "__main__":
    import os

    try:
        path = sys.argv[1]
    except IndexError:
        print "usage: ansotropy.py <PDB/mmCIF file>"

    if os.path.isfile(path):
        main(path)
    elif os.path.isdir(path):
        for name in os.listdir(path):
            name = os.path.join(path, name)

            if not os.path.isfile(name):
                continue

            try:
                main(name)
            except:
                print "ERROR: ",name
                raise
