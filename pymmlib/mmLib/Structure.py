## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
from   __future__           import generators
import fpformat
import weakref

from   mmTypes              import *
from   AtomMath             import *
from   Library              import Library
from   StructureAlgorithms  import *



### <LISTS>
class ChainList(OrderedTupleList):
    def add(self, chain):
        assert isinstance(chain, Chain)
        self.list.append((chain.chain_id or "ZZZ", chain))
        self.list.sort()

    def get(self, chain_id):
        for item in self.list:
            if item[-1].chain_id == chain_id: return item[-1]
        return None


class PolymerList(OrderedTupleList):
    def add(self, poly):
        assert isinstance(poly, Polymer)
        self.list.append((poly.polymer_id, poly))
        self.list.sort()
            
    def get(self, polymer_id):
        try:
            return self.list[polymer_id][-1]
        except IndexError:
            print "[ERROR] BAD polymer_id=",polymer_id
        return None


class FragmentList(OrderedTupleList):
    def __init__(self):
        OrderedTupleList.__init__(self)
        self.cache = {}

    def add(self, frag):
        assert isinstance(frag, Fragment)
        fid = FragmentID(frag.fragment_id)
        self.list.append((fid.res_seq, fid.icode, frag))
        self.cache[frag.fragment_id] = frag
        self.list.sort()

    def get(self, fragment_id):
        try:
            return self.cache[fragment_id]
        except KeyError:
            print "[ERROR] BAD frag_id=",fragment_id
        return None


class AtomList(OrderedList):
    def __iter__(self, alt_loc):
        for atm in self.list:
            if   not atm.alt_loc:        yield atm
            elif atm.alt_loc == alt_loc: yield atm

    def add(self, atom):
        assert isinstance(atom, Atom)

        self.list.append(atom)

    def get(self, name, alt_loc):
        for atm in self.list:
            if atm.name == name and \
               (not atm.alt_loc or atm.alt_loc == alt_loc):
                return atm
        return None


class AltLocList(WeakrefList):
    def add(self, atom):
        assert isinstance(atom, Atom)

        for i in range(len(self)):
            if atom.alt_loc < self[i].alt_loc:
                self.insert(i, atom)
                return

        self.append(atom)


class BondList(OrderedList):
    def add(self, bond):
        assert isinstance(bond, Bond)

        self.list.append(bond)
### </LISTS>



### <STRUCTURE>
class Structure(object):
    """The Structure object is the parent container object for the entire
    macromolecular data structure."""

    def __init__(self,
                 library         = None,
                 default_alt_loc = ""):

        ## element & monomer library
        self.library         = library

        ## default alt_loc
        self.default_alt_loc = default_alt_loc

        ## chain list
        self.chain_list = ChainList()

    def __str__(self):
        return "Structure"

    def __getitem__(self, chain_id):
        """Same as getChain, but raises KeyError if the requested chain_id
        is not found."""
        chain = self.getChain(chain_id)
        if not chain: raise KeyError
        return chain

    def __iter__(self):
        return self.iterChains()

    def iterChains(self):
        """Iterates over all Chain objects in alphabetical order according
        to their chain_id."""
        return iter(self.chain_list)

    def getChain(self, chain_id):
        """Returns the Chain object matching the chain_id charactor."""
        return self.chain_list.get(chain_id)

    def iterPolymers(self):
        """Iterate over all Segment objects.  The iteration is preformed in
        order according to the Fragment object ordering the Segment objects
        span."""
        for chain in self.iterChains():
            for poly in chain.iterPolymers():
                yield poly

    def iterPolypeptides(self):
        for chain in self.iterChains():
            for poly in chain.iterPolypeptides():
                yield poly
 
    def iterDNA(self):
        for chain in self.iterChains():
            for dna in chain.iterDNA():
                yield dna

    def getFragment(self, chain_id, fragment_id):
        """Returns the PDB fragment uniquely identified by its chain_id,
        res_seq, and icode."""
        chain = self.getChain(chain_id)
        if chain:
            return chain.getFragment(fragment_id)
        return None
            
    def iterFragments(self):
        """Iterates over all Fragment objects.  The iteration is performed
        in order according the the parent Chain's chain_id, and the
        Fragment's positioin within the chain."""
        for chain in self.iterChains():
            for frag in chain.iterFragments():
                yield frag

    def iterAminoAcids(self):
        """Same as iterFragments() but only iterates over Fragments of the
        subclass AminoAcidResidue."""
        for chain in self.iterChains():
            for aa in chain.iterAminoAcids():
                yield aa

    def iterNucleicAcids(self):
        """Same as iterFragments() but only iterates over Fragments of the
        subclas NucleicAcidResidue."""
        for chain in self.iterChains():
            for aa in chain.iterNucleicAcids():
                yield aa

    def getAtom(self, chain_id, fragment_id, name):
        """Returns the Atom object matching the given set of arguments."""
        chain = self.getChain(chain_id)
        if chain:
            return chain.getAtom(fragment_id, name)
        return None

    def iterAtoms(self):
        """Iterates over all Atom objects.  The iteration is preformed in
        order according to the Chain and Fragment ordering rules the Atom
        object is a part of."""
        for chain in self.iterChains():
            for atm in chain.iterAtoms():
                yield atm

    def iterBonds(self):
        """Iterates over all Bond objects.  The iteration is preformed by
        iterating over all Atom objects in the same order as iterAtoms(),
        then iterating over each Atom's Bond objects."""
        visited = []

        for atm in self.iterAtoms():
            for bond in atm.iterBonds():
                if bond not in visited:
                    yield bond
                    visited.insert(0, bond)


class StructureMember(object):
    """All objects which are contained within a Structure object should be
    subclassed from StructureMember."""
    
    def __init__(self, structure):
        assert isinstance(structure, Structure)
        
        self.getStructure = weakref.ref(structure)
### </STRUCTURE>



### <CHAIN>
class Chain(StructureMember):
    """Chain objects conatain a ordered list of Fragment objects."""

    def __init__(self,
                 structure,
                 chain_id    = ""):
        
        StructureMember.__init__(self, structure)

        self.chain_id      = chain_id
        self.polymer_list  = PolymerList()
        self.fragment_list = FragmentList()

    def __str__(self):
        return "Chain(%s, %s...%s)" % (self.chain_id,
                                       self.fragment_list[0],
                                       self.fragment_list[-1])

    def __iter__(self):
        """Same as iterFragments()."""
        return self.iterFragments()
        
    def __getitem__(self, fragment_id):
        """This is a alternative to calling getFragment, but a KeyError
        exception is raised if the Fragment is not found."""
        frag = self.getFragment(fragment_id)
        if not frag: raise KeyError
        return frag

    def getPolymer(self, polymer_id):
        """Returns the Polymer object matching the given polymer_id.  Returns
        None if no match is found."""
        return self.polymer_list.get(polymer_id)

    def iterPolymers(self):
        """Iterate over all Polymer objects with the Chain."""
        return iter(self.polymer_list)

    def iterPolypeptides(self):
        """Iterates over all Polypeptide object within the Chain."""
        for poly in self.iterPolymers():
            if isinstance(poly, Polypeptide):
                yield poly

    def iterDNA(self):
        """Iterates over all DNA objects within the Chain."""
        for poly in self.iterPolymers():
            if isinstance(poly, DNA):
                yield poly

    def getFragment(self, fragment_id):
        """Returns the PDB fragment uniquely identified by its chain_id,
        res_seq, and icode."""
        return self.fragment_list.get(fragment_id)

    def iterFragments(self):
        """Iterates over all Fragment objects.  The iteration is performed
        in order according to the Fragment's position within the Chain
        object."""
        return iter(self.fragment_list)

    def iterAminoAcids(self):
        """Same as iterFragments(), but only iterates over AminoAcidResidue
        objects."""
        for frag in self.iterFragments():
            if isinstance(frag, AminoAcidResidue):
                yield frag

    def iterNucleicAcids(self):
        """Same as iterFragments(), but only iterates over NucleicAcidResidue
        objects."""
        for frag in self.iterFragments():
            if isinstance(frag, NucleicAcidResidue):
                yield frag

    def iterAtoms(self):
        """Iterates over all Atom objects within the Chain."""
        for frag in self.iterFragments():
            for atm in frag.iterAtoms():
                yield atm
                
    def iterBonds(self):
        """Iterates over all Bond objects attached to Atom objects within the
        Chain"""
        visited = []
        for atm in self.iterAtoms():
            for bond in atm.iterBonds():
                if bond not in visited:
                    yield bond
                    visited.insert(0, bond)
### </CHAIN>



### <POLYMER>
class Polymer(StructureMember):
    """Abstract class for grouping a series of bonded Residue objects
    into a Polymer."""

    def __init__(self,
                 structure,
                 chain_id     = "",
                 polymer_id   = 0):
        
        StructureMember.__init__(self, structure)
        
        self.chain_id      = chain_id
        self.polymer_id    = polymer_id
        self.fragment_list = FragmentList()

    def __str__(self):
        try:
            frag1 = self.fragment_list[0]
            frag2 = self.fragment_list[-1]
        except IndexError:
            frag1 = None
            frag2 = None
        
        return "%s(%d, %s...%s)" % (self.__class__.__name__,
                                    self.polymer_id,
                                    frag1, frag2)
    
    def __iter__(self):
        """Same as iterResidues()."""
        return self.iterResidues()
        
    def __getitem__(self, fragment_id):
        """This is a alternative to calling getFragment, but a KeyError
        exception is raised if the fragment is not found."""
        frag = self.getFragment(fragment_id)
        if not frag: raise KeyError
        return frag

    def getFragment(self, fragment_id):
        """Returns the Residue object matching the argument res_seq and icode
        within the Polymer.  Returns None if no match is found."""
        return self.fragment_list.get(fragment_id)

    def iterFragments(self):
        """Iterates over all Residue objects within the Polymer."""
        return iter(self.fragment_list)

    def getResidue(self, fragment_id):
        """Returns the Residue object matching the argument res_seq and icode
        within the Polymer.  Returns None if no match is found."""
        return self.fragment_list.get(fragment_id)

    def iterResidues(self):
        """Iterates over all Residue objects within the Polymer."""
        return iter(self.fragment_list)

    def iterAtoms(self):
        """Iterates over all Atom objects with the Polymer."""
        for frag in self.iterFragments():
            for atm in frag.iterAtoms():
                yield atm

    def iterBonds(self):
        """Iterates over all Bond objects attached to Atom objects within the
        Polymer."""
        visited = []
        for atm in self.iterAtoms():
            for bond in atm.iterBonds():
                if bond not in visited:
                    yield bond
                    visited.insert(0, bond)

class Polypeptide(Polymer):
    """Subclass of the Polymer object grouping a chain of
    AminoAcidResidue objects into a polypeptide."""


class DNA(Polymer):
    """Subclass of the Polymer object grouping a chain of
    NucleicAcidResidue object into DNA."""
### </POLYMER>



### <FRAGMENT>
class Fragment(StructureMember):
    """Fragment objects are a basic unit for organizing small groups of
    Atoms.  Amino acid residues are fragments, as well as nucleic
    acids and other small molecules.  In terms of a PDB file, they are
    all the atoms from a unique residue in a chain.  Fragments have the
    following attributes:

    Fragment.res_name     - the fragment/residue name
    Fragment.res_seq      - the sequence id of the fragment/residue
    Fragment.chain_id     - the ID of the chain containing this fragment
    """

    def __init__(self,
                 structure,
                 res_name    = "",
                 fragment_id = "",
                 chain_id    = ""):
        
        StructureMember.__init__(self, structure)

        self.res_name    = res_name
        self.fragment_id = fragment_id
        self.chain_id    = chain_id

        self.atom_list = AtomList()

    def __str__(self):
        return "%s(%s,%s,%s)" % (self.__class__.__name__,
                                 self.res_name,
                                 self.fragment_id,
                                 self.chain_id)

    def __iter__(self):
        """Same as iterAtoms() method."""
        return self.iterAtoms()

    def __getitem__(self, name):
        """This is a alternative to calling the getAtom method, except
        a KeyError exception is raised if the Atom object is not found."""
        atm = self.getAtom(name)
        if not atm: raise KeyError
        return atm

    def getAtom(self, name):
        """Returns the matching Atom object contained in the Fragment.
        Returns None if a match is not found."""
        return self.atom_list.get(name, self.getStructure().default_alt_loc)
    
    def iterAtoms(self):
        """Iterates over all Atom objects contained in the Fragment.
        There is no defined order for the iteration."""
        return self.atom_list.__iter__(self.getStructure().default_alt_loc)

    def iterBonds(self):
        """Iterates over all Bond objects.  The iteration is preformed by
        iterating over all Atom objects in the same order as iterAtoms(),
        then iterating over each Atom's Bond objects."""
        visited = []

        for atm in self.iterAtoms():
            for bond in atm.iterBonds():
                if bond not in visited:
                    yield bond
                    visited.insert(0, bond)

    def getOffsetFragment(self, offset):
        chain = self.getStructure().getChain(self.chain_id)
        if not chain: return None

        i = chain.fragment_list.index(self) + offset
        if i < 0 or i >= len(chain.fragment_list): return None
        return chain.fragment_list[i]


class Residue(Fragment):
    """A subclass of Fragment representing one residue in a polymer chain."""

    polymer_class          = Polymer

    def __init__(self,
                 structure,
                 res_name    = "",
                 fragment_id = "",
                 chain_id    = "",
                 polymer_id  = 0):

        Fragment.__init__(self,
                          structure,
                          res_name    = res_name,
                          fragment_id = fragment_id,
                          chain_id    = chain_id)

        self.polymer_id = polymer_id

    def __str__(self):
        return "%s(%s,%s,%s,%d)" % (self.__class__.__name__,
                                    self.res_name,
                                    self.fragment_id,
                                    self.chain_id,
                                    self.polymer_id)

    def getOffsetResidue(self, offset):
        chain = self.getStructure().getChain(self.chain_id)
        if not chain: return None

        poly = chain.getPolymer(self.polymer_id)
        if not poly: return None

        i = poly.fragment_list.index(self) + offset
        if i < 0 or i >= len(poly.fragment_list): return None
        return poly.fragment_list[i]



class AminoAcidResidue(Residue, AminoAcidResidueAlgorithms):
    """A subclass of Residue representing one amino acid residue in a
    polypeptide chain."""

    polymer_class          = Polypeptide


class NucleicAcidResidue(Residue):
    """A subclass of Residue representing one nuclic acid in a strand of
    DNA or RNA."""

    polymer_class          = DNA
### </FRAGMENT>



### <ATOM>
class Bond(BondAlgorithms):
    """Indicates two atoms are bonded together."""
    def __init__(self, atom1, atom2):
        assert isinstance(atom1, Atom)
        assert isinstance(atom2, Atom)

        assert atom1 != atom2

        self.getAtom1 = weakref.ref(atom1)
        self.getAtom2 = weakref.ref(atom2)

        atom1.bond_list.add(self)
        atom2.bond_list.add(self)

    def __str__(self):
        return "Bond(%s...%s)" % (self.getAtom1(), self.getAtom2())

    def getPartner(self, atm):
        if   atm == self.getAtom1(): return self.getAtom2()
        elif atm == self.getAtom2(): return self.getAtom1()
        return None


class Atom(StructureMember, AtomAlgorithms):
    """Class representing a single atom.  Atoms have the following default
    attributes.  If a attribue has the value None, then the attribue was
    never set.  If the attribue has a default, then it is required.

    Atom[alt_loc]    - Atom objects in alternate locations can be accessed
                       by using Python's dictionary syntax with the alt_loc
                       charactor
    iter(Atom)       - iterates over all alt_loc versions of the Atom
    Atom.name        - label of the atom
    Atom.alt_loc     - alternate location indicater for the atom
    Atom.res_name    - the name of the resiude/fragment this atom is part of
    Atom.res_seq     - the residue/fragment sequence number
    Atom.icode       - the insertion code for the residue/fragment
    Atom.chain_id    - the chain ID of the chain containing this atom
    Atom.element     - symbol for the element
    Atom.position    - a Vector (ScientificPython)
    Atom.occupancy
    Atom.temp_factor
    Atom.U           - a 6-tuple of the anisotropic values
    Atom.charge      - charge on the atom
    """
    def __init__(self,
                 structure,
                 name        = "",
                 alt_loc     = "",
                 res_name    = "",
                 fragment_id = "",
                 chain_id    = ""):

        assert isinstance(structure, Structure)

        StructureMember.__init__(self, structure)

        self.name         = name
        self.alt_loc      = alt_loc
        self.res_name     = res_name
        self.fragment_id  = fragment_id
        self.chain_id     = chain_id

        self.element      = None
        self.position     = None
        self.occupancy    = None
        self.temp_factor  = None
        self.U            = None
        self.charge       = None
        
        self.alt_loc_list = None
        self.bond_list    = BondList()

    def __str__(self):
        return "Atom(%s,%s,%s,%s,%s)" % (self.name,
                                         self.alt_loc,
                                         self.res_name,
                                         self.fragment_id,
                                         self.chain_id)

    def __getitem__(self, alt_loc):
        """This is a alternative to calling getAltLoc, but a KeyError
        exception is raised if the alt_loc Atom is not found."""
        atm = self.getAltLoc(alt_loc)
        if not atm: raise KeyError
        return atm

    def __iter__(self):
        """Same as iterAltLoc()."""
        return self.iterAltLoc()

    def getAltLoc(self, alt_loc):
        """Returns the Atom object matching the alt_loc argument."""
        if self.alt_loc_list:
            alt_loc = alt_loc or self.getStructure().default_alt_loc
            for atm in self.alt_loc_list:
                if atm.alt_loc == alt_loc: return atm
            return None
        else:
            return self

    def iterAltLoc(self):
        """Iterate over all alt_loc versions of this atom in the
        alphabetical order of the alt_loc labels.  If there are no
        alt_loc versions, do not iterate."""
        if self.alt_loc_list:
            for atm in self.alt_loc_list:
                yield atm
        else:
            yield self
            
    def getBond(self, atm):
        """Returns the Bond connecting self with the argument atom."""
        for bond in self.bond_list:
            if atm == bond.getAtom1() or atm == bond.getAtom2():
                return bond
        return None

    def iterBonds(self):
        """Iterates over all the Bond edges connected to self."""
        for bond in self.bond_list:
            yield bond

    def iterBondedAtoms(self):
        """Iterates over all the Atoms bonded to self."""
        for bond in self.iterBonds():
            yield bond.getPartner(self)
### </ATOM>
