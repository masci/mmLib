#!/usr/bin/env python
# this example tests the API of the structure hierarchy

import sys
from mmLib.FileLoader import LoadStructure
from mmLib.Structure import *

def atom_test(atom):
    print "Atom Test"
    print atom

    for atm in atom.iter_alt_loc():
        assert isinstance(atm, Atom)
        assert atom[atm.alt_loc] == atm
        assert atm.get_alt_loc(atm.alt_loc) == atm
        assert atm.get_fragment() == atom.get_fragment()
        assert atm.get_chain() == atom.get_chain()
        assert atm.get_structure() == atom.get_structure()

    for bond in atom.iter_bonds():
        assert isinstance(bond, Bond)

    print "Atom Test Complete"

def fragment_test(frag):
    print "Fragment Test"
    print frag

    for atm in frag.iter_atoms():
        assert isinstance(atm, Atom)
        assert frag[atm.name] == atm
        assert frag.get_atom(atm.name) == atm
        assert atm.get_fragment() == frag
        assert atm.get_chain() == frag.get_chain()
        assert atm.get_structure() == frag.get_structure()

    for bond in frag.iter_bonds():
        assert isinstance(bond, Bond)

    print "Fragment Test Complete"

def chain_test(chain):
    print "Chain Test"
    print chain

    for frag in chain.iter_fragments():
        assert isinstance(frag, Fragment)
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
    struct = LoadStructure(path)
    struct_test(struct)

    for chain in struct.iter_chains():
        chain_test(chain)

    for frag in struct.iter_fragments():
        fragment_test(frag)

    for atm in struct.iter_atoms():
        atom_test(atm)

if __name__ == "__main__":
    try:
        path = sys.argv[1]
    except IndexError:
        print "usage: ansotropy.py <PDB/mmCIF file>"

    main(path)
