## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Classes for representing biological macromolecules.
"""
from   __future__ import generators
import fpformat
import weakref
from mmTypes import *
from AtomMath import *
from Library import Library
from ExpData import ExpData
from UnitCell import UnitCell
from SpaceGroup import SpaceGroup, SymmetryOperator


class Structure(object):
    """The Structure object is the parent container object for the entire
    macromolecular data structure.

    Attributes:
    library          (mmLib.Library object) the monomer library used by
                              the structure

    exp_data         Experimental Data class containing non-structural
                     data from the structure source file
                     
    unit_cell        UnitCell class containing the unit cell descripton

    space_group      SpaceGroup class containing space group description

    default_alt_loc  (string) all objects in the structure hierarcy
                              default to iterating, retrieving, and
                              calculating the structure using the alt_loc
                              conformation set in default_alt_loc
    sites            (list)   list of Site objects containing all defined
                              sites in the structure
    alpha_helices    (list)   list of AlphaHelix objects containing all
                              alpha helices in the structure
    beta_sheets      (list)   list of BetaSheet objects containing all
                              beta sheets in the structure
    turns            (list)   list of Turn objects
    """
    def __init__(self, library = None):
        if library:
            self.library = library
        else:
            self.library = Library()

        self.exp_data        = ExpData()
        self.unit_cell       = UnitCell()
        self.space_group     = SpaceGroup()

        self.default_alt_loc = "A"
        self.sites           = []
        self.alpha_helices   = []
        self.beta_sheets     = []
        self.turns           = []
        self.__chain_list    = []

    def __str__(self):
        tstr =  "Structure::%s\n" % (self.exp_data.get("id", ""))
        tstr += "%s\n" % (self.exp_data.get("title", "no title")[:40])
        try:
            tstr += "res(A): %1.2f\n" % (self.exp_data["res_high"])
            tstr += "R: %1.3f\n" % (self.exp_data["R_fact"])
            tstr += "Rfree: %1.3f\n" % (self.exp_data["free_R_fact"])
        except KeyError:
            pass
        return tstr
 
    def __len__(self):
        """Returns the number of stored Chain objects.
        """
        return len(self.__chain_list)
    
    def __getitem__(self, x):
        """Same as get_chain, but raises KeyError if the requested chain_id
        is not found.
        """
        if type(x) == StringType:
            for chain in self.__chain_list:
                if chain.chain_id == x:
                    return chain
            raise KeyError, x

        elif type(x) == IntType:
            return self.__chain_list[x]

        raise TypeError, x

    def __delitem__(self, x):
        self.__chain_list.remove(self[x])

    def __iter__(self):
        """Iterates the Chain objects in the Structure.
        """
        return iter(self.__chain_list)

    def __contains__(self, x):
        if isinstance(x, Chain):
            return x in self.__chain_list
        elif type(x) == StringType:
            return self[x] in self.__chain_list
        raise TypeError, x

    def index(self, chain):
        """Returns the numeric index of the Chain object.
        """
        assert isinstance(chain, Chain)
        return self.__chain_list.index(chain)

    def remove(self, chain):
        """Removes the 
        """
        assert isinstance(chain, Chain)
        del chain.get_structure
        self.__chain_list.remove(chain)

    def sort(self):
        self.__chain_list.sort()

    def add_chain(self, chain, delay_sort = True):
        assert isinstance(chain, Chain)
        chain.get_structure = weakref.ref(self)
        self.__chain_list.append(chain)

        if not delay_sort:
            self.__chain_list.sort()

    def get_chain(self, chain_id):
        """Returns the Chain object matching the chain_id charactor.
        """
        try:
            return self[chain_id]
        except KeyError:
            return NOne

    def iter_chains(self):
        """Iterates over all Chain objects in alphabetical order according
        to their chain_id.
        """
        return iter(self)

    def iter_fragments(self):
        """Iterates over all Fragment objects.  The iteration is performed
        in order according the the parent Chain's chain_id, and the
        Fragment's positioin within the chain.
        """
        for chain in self.iter_chains():
            for frag in chain.iter_fragments():
                yield frag

    def iter_amino_acids(self):
        """Same as iter_fragments() but only iterates over Fragments of the
        subclass AminoAcidResidue.
        """
        for chain in self.iter_chains():
            for aa in chain.iter_amino_acids():
                yield aa

    def iter_nucleic_acids(self):
        """Same as iter_fragments() but only iterates over Fragments of the
        subclas NucleicAcidResidue.
        """
        for chain in self.iter_chains():
            for aa in chain.iter_nucleic_acids():
                yield aa

    def iter_atoms(self):
        """Iterates over all Atom objects.  The iteration is preformed in
        order according to the Chain and Fragment ordering rules the Atom
        object is a part of.
        """
        for chain in self.iter_chains():
            for atm in chain.iter_atoms():
                yield atm

    def iter_bonds(self):
        """Iterates over all Bond objects.  The iteration is preformed by
        iterating over all Atom objects in the same order as iter_atoms(),
        then iterating over each Atom's Bond objects.
        """
        visited = []

        for atm in self.iter_atoms():
            for bond in atm.iter_bonds():
                if bond not in visited:
                    yield bond
                    visited.insert(0, bond)


class Chain(object):
    """Chain objects conatain a ordered list of Fragment objects.
    """
    def __init__(self, chain_id = ""):
        self.chain_id = chain_id

        ## the sequence list contains a list 3-letter residue names
        self.sequence = []

        ## fragments are contained in the list and also cached in
        ## a dictionary for fast random-access lookup
        self.__fragment_list  = []

    def __str__(self):
        return "Chain(%s, %s...%s)" % (self.chain_id,
                                       self.__fragment_list[0],
                                       self.__fragment_list[-1])

    def __lt__(self, other):
        assert isinstance(other, Chain)
        return self.chain_id < other.chain_id
        
    def __le__(self, other):
        assert isinstance(other, Chain)
        return self.chain_id <= other.chain_id
        
    def __gt__(self, other):
        assert isinstance(other, Chain)
        return self.chain_id > other.chain_id

    def __ge__(self, other):
        assert isinstance(other, Chain)
        return self.chain_id >= other.chain_id

    def __len__(self):
        """Return the number of fragments in the Chain."""
        return len(self.__fragment_list)

    def __getitem__(self, x):
        """Retrieve a Fragment within the Chain.  This can take a integer
        index of the Fragment's position within the chain, the fragment_id
        string of the Fragment to retrieve, or a slice of the Chain to
        return a new Chain object containing the sliced subset of Fragments.
        """
        if type(x) == IntType:
            return self.__fragment_list[x]
        
        elif type(x) == StringType:
            for frag in self:
                if frag.fragment_id == x:
                    return frag
            raise KeyError, x

        elif type(x) == SliceType:
            ## slices are a difficult concept to get right here, the
            ## "Python way" is to return a object of the same type
            ## so I return a new Chain object with only the fragments
            ## in the slice range, but I don't copy the fragments nor do
            ## I set the fragments's get_chain() function to return the
            ## new Chain object
            chain               = Chain(self.chain_id)
            chain.get_structure = self.get_structure

            frag_start = self[x.start]
            frag_stop  = self[x.stop]

            for frag in self:
                if frag >= frag_start and frag <= frag_stop:
                    chain.__fragment_list.append(frag)
                    if frag.fragment_id in self.sequence:
                        chain.sequence.append(frag.fragment_id)

            return chain

        raise TypeError, x

    def __delitem__(self, x):
        """Delete Fraagment from the chain.  This can take a reference to the
        Fragment object to delete, the fragment_id of the Fragment to delete,
        or the integer index of the Fragment within the Chain.
        """
        frag = self[x]
        self.__fragment_list.remove(frag)

    def __iter__(self):
        """Iterate all Fragments contained in the Chain.
        """
        return iter(self.__fragment_list)

    def __contains__(self, x):
        if isinstance(x, Fragment):
            return x in self.__fragment_list
        elif type(x) == StringType:
            return self[x] in self.__fragment_list
        raise TypeError, x

    def index(self, frag):
        """Return the 0-based index of the framgent in the chain list.
        """
        return self.__fragment_list.index(frag)

    def remove(self, frag):
        """Remove the Fragment from the chain.
        """
        del frag.get_chain
        self.__fragment_list.remove(frag)

    def sort(self):
        """Sort the Fragments in the chain into proper order.
        """
        self.__fragment_list.sort()
        
    def add_fragment(self, frag, delay_sort = False):
        """Adds a Fragment instance to the chain.  If delay_sort is True,
        then the fragment is not inserted in the proper position within the
        chain.
        """
        assert isinstance(frag, Fragment)
        assert frag.chain_id == self.chain_id

        frag.get_chain = weakref.ref(self)
        self.__fragment_list.append(frag)

        if not delay_sort:
            self.__fragment_list.sort()

    def get_fragment(self, fragment_id):
        """Returns the PDB fragment uniquely identified by its fragment_id.
        """
        try:
            return self[fragment_id]
        except KeyError:
            return None

    def iter_fragments(self):
        """Iterates over all Fragment objects.  The iteration is performed
        in order according to the Fragment's position within the Chain
        object.
        """
        return iter(self)

    def iter_amino_acids(self):
        """Same as iter_fragments(), but only iterates over AminoAcidResidue
        objects.
        """
        for frag in self.iter_fragments():
            if isinstance(frag, AminoAcidResidue):
                yield frag

    def iter_nucleic_acids(self):
        """Same as iter_fragments(), but only iterates over NucleicAcidResidue
        objects.
        """
        for frag in self.iter_fragments():
            if isinstance(frag, NucleicAcidResidue):
                yield frag

    def has_standard_residues(self):
        """Returns True if the chain contains standard residues as defined
        by the PDB.  Standard residues are amino and nucleic acid resiudes.
        """
        for frag in self.iter_fragments():
            if frag.is_standard_residue():
                return True
        return False

    def iter_standard_residues(self):
        """Iterates over standard residues in the chain, as defined by the
        PDB.  Standard residues are amino and nucleic acid residues.
        """
        for frag in self.iter_fragments():
            if frag.is_standard_residue():
                yield frag

    def iter_non_standard_residues(self):
        """Iterates over non-standard residues in the chain, as defined
        by the PDB.  Non-standard residues are anything which is not a
        amino or nucleic acid.
        """
        for frag in self.iter_fragments():
            if not frag.is_standard_residue():
                yield frag

    def iter_atoms(self):
        """Iterates over all Atom objects within the Chain.
        """
        for frag in self.iter_fragments():
            for atm in frag.iter_atoms():
                yield atm
                
    def iter_bonds(self):
        """Iterates over all Bond objects attached to Atom objects within the
        Chain.
        """
        visited = []
        for atm in self.iter_atoms():
            for bond in atm.iter_bonds():
                if bond not in visited:
                    yield bond
                    visited.insert(0, bond)

    def set_chain_id(self, chain_id):
        """Sets a new ID for the Chain object, updating the chain_id
        for all objects in the Structure hierarchy.
        """
        ## check for conflicting chain_id in the structure
        try:             self.get_structure()[chain_id]
        except KeyError: pass
        else:            raise ValueError, chain_id

        ## set the new chain_id in all the additional groups

        ## set the new chain_id for the chain object (self)
        self.chain_id = chain_id

        ## set the chain_id in all the fragment and atom children
        for frag in self:
            frag.chain_id = chain_id

            for atm in frag:
                atm.chain_id = chain_id

        ## resort the parent structure
        self.get_structure().sort()

    def calc_sequence(self):
        """Attempts to calculate the residue sequence contained in the
        Chain object.  This is a simple algorithm: find the longest running
        sequence of the same bio-residue, and that's the sequence.  Returns
        a list of 3-letter residues codes of the calculated sequence.
        """
        residue_class = None
        sequence_list = []
        for frag in self.iter_standard_residues():
            if residue_class:
                if not isinstance(frag, residue_class):
                    break
            else:
                residue_class = frag.__class__

            sequence_list.append(frag.res_name)

        return sequence_list
            

class Fragment(object):
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
                 res_name    = "",
                 fragment_id = "",
                 chain_id    = ""):
        
        self.res_name    = res_name
        self.fragment_id = fragment_id
        self.chain_id    = chain_id

        self.__atom_list = []

    def __str__(self):
        return "%s(%s,%s,%s)" % (self.__class__.__name__,
                                 self.res_name,
                                 self.fragment_id,
                                 self.chain_id)

    def __lt__(self, other):
        assert isinstance(other, Fragment)
        return FragmentID(self.fragment_id) < FragmentID(other.fragment_id)
        
    def __le__(self, other):
        assert isinstance(other, Fragment)
        return FragmentID(self.fragment_id) <= FragmentID(other.fragment_id)

    def __gt__(self, other):
        assert isinstance(other, Fragment)
        return FragmentID(self.fragment_id) > FragmentID(other.fragment_id)

    def __ge__(self, other):
        assert isinstance(other, Fragment)
        return FragmentID(self.fragment_id) >= FragmentID(other.fragment_id)

    def __len__(self):
        return len(self.__atom_list)

    def __getitem__(self, x):
        if   type(x) == IntType:
            return self.__atom_list[x]

        elif type(x) == StringType:
            alt_loc = self.get_chain().get_structure().default_alt_loc
            
            for atom in self.__atom_list:
                if atom.name == x:
                    if atom.alt_loc == alt_loc:
                        return atom
                    elif atom.alt_loc == None or atom.alt_loc == "":
                        return atom
                    
            raise KeyError, x

        raise TypeError, x

    def __delitem__(self, x):
        self.__atom_list.remove(self[x])

    def __iter__(self):
        alt_loc = self.get_chain().get_structure().default_alt_loc
        for atom in self.__atom_list:
            if atom.alt_loc == alt_loc:
                yield atom
            elif atom.alt_loc == None or atom.alt_loc == "":
                yield atom

    def __contains__(self, x):
        if isinstance(x, Atom):
            return x in self.__atom_list
        elif type(x) == StringType:
            return self[x] in self.__atom_list
        raise TypeError, x

    def index(self, atom):
        assert isinstance(atom, Atom)
        return self.__atom_list.index(atom)

    def remove(self, atom):
        assert isinstance(atom, Atom)
        del atom.get_fragment
        atom.alt_loc_list.remove(atom)
        self.__atom_list.remove(atom)

    def add_atom(self, atom):
        assert isinstance(atom, Atom)
        assert atom.chain_id    == self.chain_id
        assert atom.fragment_id == self.fragment_id
        assert atom not in self.__atom_list
        
        atom.get_fragment = weakref.ref(self)

        for a in self.__atom_list:
            if atom.name == a.name:
                a.add_alt_loc(atom)
                break

        self.__atom_list.append(atom)

    def get_atom(self, name):
        """Returns the matching Atom object contained in the Fragment.
        Returns None if a match is not found.
        """
        try:
            return self[name]
        except KeyError:
            return None
    
    def iter_atoms(self):
        """Iterates over all Atom objects contained in the Fragment.
        There is no defined order for the iteration.
        """
        return iter(self)

    def iter_bonds(self):
        """Iterates over all Bond objects.  The iteration is preformed by
        iterating over all Atom objects in the same order as iter_atoms(),
        then iterating over each Atom's Bond objects."""
        visited = []

        for atm in self.iter_atoms():
            for bond in atm.iter_bonds():
                if bond not in visited:
                    yield bond
                    visited.insert(0, bond)

    def get_offset_fragment(self, offset):
        """Returns the fragment in the same chain at integer offset from
        self.  Returns None if no fragment is found.
        """
        assert type(offset) == IntType
        
        chain = self.get_chain()
        i     = chain.index(self) + offset
        try:               return chain[i]
        except IndexError: return None

    def get_structure(self):
        """Returns the parent structure.
        """
        return self.get_chain().get_structure()

    def set_fragment_id(self, fragment_id):
        """Sets a new ID for the Fragment object, updating the fragment_id
        for all objects in the Structure hierarchy.
        """
        ## check for conflicting chain_id in the structure
        try:             self.get_chain()[fragment_id]
        except KeyError: pass
        else:            raise ValueError, fragment_id

        ## set the new fragment_id in all the additional groups

        ## set the new chain_id for the chain object (self)
        self.fragment_id = fragment_id

        ## set the chain_id in all the fragment and atom children
        for atm in self:
            atm.fragment_id = fragment_id

        ## resort the parent chain
        self.get_chain().sort()

    def create_bonds(self):
        """Contructs bonds within a fragment.  Bond definitions are retrieved
        from the monomer library.
        """
        mon = self.get_structure().library.get_monomer(self.res_name)
        if not mon:
            return

        for (name1, name2) in mon.bond_list:
            try:
                atm1 = self[name1]
                atm2 = self[name2]
            except KeyError:
                continue
            else:
                atm1.create_bond(atm2, bond_alt_loc = True)

    def is_standard_residue(self):
        """Returns True if the Fragment/Residue object is one of the
        PDB defined standard residues.  PDB standard residues are amino
        and nucleic acid residues.
        """
        return isinstance(self, AminoAcidResidue) or \
               isinstance(self, NucleicAcidResidue)


    def is_water(self):
        """Returns True if the Fragment is a water molecule, returns False
        otherwise.
        """
        if self.get_structure().library.is_water(self.res_name):
            return True
        return False


class Residue(Fragment):
    """A subclass of Fragment representing one residue in a polymer chain.
    """
    def __init__(self,
                 res_name    = "",
                 fragment_id = "",
                 chain_id    = ""):

        Fragment.__init__(self,
                          res_name    = res_name,
                          fragment_id = fragment_id,
                          chain_id    = chain_id)

    def __str__(self):
        return "%s(%s,%s,%s)" % (self.__class__.__name__,
                                 self.res_name,
                                 self.fragment_id,
                                 self.chain_id)

    def get_offset_residue(self, offset):
        """Returns the residue along the chain at the given integer offset
        from self.  Returns None if there is no residue at that offset, or
        if the fragment found is not the same type of residue as self.
        """
        assert type(offset) == IntType
        
        frag = Fragment.get_offset_fragment(self, offset)
        if type(self) == type(frag):
            return frag
        return None

    def create_bonds(self):
        """Contructs bonds within a fragment.  Bond definitions are retrieved
        from the monomer library.  This version also constructs the bonds
        between adjectent residues.
        """
        Fragment.create_bonds(self)

        next_res = self.get_offset_residue(1)
        if not next_res:
            return

        mon1 = self.get_structure().library.get_monomer(self.res_name)
        mon2 = self.get_structure().library.get_monomer(next_res.res_name)
        if not (mon1 and mon2):
            return

        for (name1, name2) in mon1.get_polymer_bond_list(self, next_res):
            try:
                atm1 = self[name1]
                atm2 = next_res[name2]
            except KeyError:
                continue
            else:
                atm1.create_bond(atm2, bond_alt_loc = True)


class AminoAcidResidue(Residue):
    """A subclass of Residue representing one amino acid residue in a
    polypeptide chain.
    """
    def calc_mainchain_bond_length(self):
        """Calculates the main chain bond lengths: (N-CA, CA-C, C-O, CA-CB,
        CA-(next)N).  The result is returned as a 5-tuple in that order.  Bond
        lengths involving missing atoms are returned as None in the tuple.
        """
        aN  = self.get_atom('N')
        aCA = self.get_atom('CA')
        aC  = self.get_atom('C')
        aO  = self.get_atom('O')
        aCB = self.get_atom('CB')

        try:
            naN = self.get_offset_residue(1).get_atom('N')
        except AttributeError:
            naN = None
     
        N_CA  = calculateDistance(aN, aCA)
        CA_C  = calculateDistance(aCA, aC)
        C_O   = calculateDistance(aC, aO)
        C_nN  = calculateDistance(aC, naN)
        CA_CB = calculateDistance(aCA, aCB)
        return (N_CA, CA_C, C_O, CA_CB, C_nN)

    def calc_mainchain_bond_angle(self):
        """Calculates main chain bond angles (N-CA-C, N-CA-CB, CB-CA-C,
        CA-C-O, CA-C-(next)N, C-(next residue)N-(next residue)CA) and
        returnst the result as a 6-tuple in that order.  Angles involving
        missing atoms are returned as None in the tuple.
        """
        aN       = self.get_atom('N')
        aCA      = self.get_atom('CA')
        aC       = self.get_atom('C')
        aO       = self.get_atom('O')
        aCB      = self.get_atom('CB')

        naN      = None
        naCA     = None
        next_res = self.get_offset_residue(1)
        if next_res:
            naN  = next_res.get_atom('N')
            naCA = next_res.get_atom('CA')

        N_CA_C   = calculateAngle(aN, aCA, aC)
        CA_C_O   = calculateAngle(aCA, aC, aO)
        N_CA_CB  = calculateAngle(aN, aCA, aCB)
        CB_CA_C  = calculateAngle(aCB, aCA, aC)
        CA_C_nN  = calculateAngle(aCA, aC, naN)
        C_nN_nCA = calculateAngle(aC, naN, naCA)

        return (N_CA_C, N_CA_CB, CB_CA_C, CA_C_O, CA_C_nN, C_nN_nCA) 

    def calc_torsion_psi(self):
        """Calculates the Psi torsion angle of the amino acid.  Raises a
        CTerminal exception if called on a C-terminal residue which does
        not have a Psi torsion angle.
        """
        next_res = self.get_offset_residue(1)
        if not next_res:
            return None

        aN  = self.get_atom('N')
        aCA = self.get_atom('CA')
        aC  = self.get_atom('C')
        naN = next_res.get_atom('N')
        return calculateTorsionAngle(aN, aCA, aC, naN)

    def calc_torsion_phi(self):
        """Calculates the Phi torsion angle of the amino acid.  Raises a
        NTerminal exception if called on a N-terminal residue which does
        not have a Phi torsion angle.
        """
        prev_res = self.get_offset_residue(-1)
        if not prev_res:
            return None

        paC = prev_res.get_atom('C')
        aN  = self.get_atom('N')
        aCA = self.get_atom('CA')
        aC  = self.get_atom('C')
        return calculateTorsionAngle(paC, aN, aCA, aC)

    def calc_torsion_omega(self):
        """Calculates the Omega torsion angle of the amino acid. Raises a
        CTerminal exception if called on a C-terminal residue which does
        not have a Omega torsion angle.
        """
        next_res = self.get_offset_residue(1)
        if not next_res:
            return None

        aCA  = self.get_atom('CA')
        aC   = self.get_atom('C')
        naN  = next_res.get_atom('N')
        naCA = next_res.get_atom('CA')
        return calculateTorsionAngle(aCA, aC, naN, naCA)

    def is_cis(self):
        """Returns true if this is a CIS amino acid, otherwise returns false.
        It uses calc_torsion_omega.
        """
        omega = self.calc_torsion_omega()
        return abs(omega) > (math.pi / 2.0)

    def calc_pucker_torsion(self, conf_id = None):
        """Calculates the Pucker torsion of a ring system.  Returns None
        for Amino Acids which do not have Pucker torsion angles.
        """
        mon = self.get_structure().library.get_monomer(self.res_name)
        if not mon or not mon.pucker_definition:
            return None

        a1 = self.get_atom(mon.pucker_definition[0])
        a2 = self.get_atom(mon.pucker_definition[1])
        a3 = self.get_atom(mon.pucker_definition[2])
        a4 = self.get_atom(mon.pucker_definition[3])
        return calculateTorsionAngle(a1, a2, a3, a4)

    def calc_torsion_chi1(self):
        mon = self.get_structure().library.get_monomer(self.res_name)
        if not mon.chi1_definition:
            return None
        
        a1 = self.get_atom(mon.chi1_definition[0])
        a2 = self.get_atom(mon.chi1_definition[1])
        a3 = self.get_atom(mon.chi1_definition[2])
        a4 = self.get_atom(mon.chi1_definition[3])
        return calculateTorsionAngle(a1, a2, a3, a4)

    def calc_torsion_chi2(self):
        mon = self.get_structure().library.get_monomer(self.res_name)
        if not mon.chi2_definition:
            return None
        
        a1 = self.get_atom(mon.chi2_definition[0])
        a2 = self.get_atom(mon.chi2_definition[1])
        a3 = self.get_atom(mon.chi2_definition[2])
        a4 = self.get_atom(mon.chi2_definition[3])
        return calculateTorsionAngle(a1, a2, a3, a4)

    def calc_torsion_chi3(self):
        mon = self.get_structure().library.get_monomer(self.res_name)
        if not mon.chi3_definition:
            return None
        
        a1 = self.get_atom(mon.chi3_definition[0])
        a2 = self.get_atom(mon.chi3_definition[1])
        a3 = self.get_atom(mon.chi3_definition[2])
        a4 = self.get_atom(mon.chi3_definition[3])
        return calculateTorsionAngle(a1, a2, a3, a4)

    def calc_torsion_chi4(self):
        mon = self.get_structure().library.get_monomer(self.res_name)
        if not mon.chi4_definition:
            return None
        
        a1 = self.get_atom(mon.chi4_definition[0])
        a2 = self.get_atom(mon.chi4_definition[1])
        a3 = self.get_atom(mon.chi4_definition[2])
        a4 = self.get_atom(mon.chi4_definition[3])
        return calculateTorsionAngle(a1, a2, a3, a4)

    def calc_torsion_chi(self):
        """Calculates CHI side-chain torsion angles according to the
        amino acid specific definitions in the AminoAcids library.
        Returns the 4-tuple (CHI1, CHI2, CHI3, CHI4).  Angles involving
        missing atoms, or angles which do not exist for the amino acid
        are returned as None in the tuple.
        """
        chi1 = self.calc_torsion_chi1()
        chi2 = self.calc_torsion_chi2()
        chi3 = self.calc_torsion_chi3()
        chi4 = self.calc_torsion_chi4()
        return (chi1, chi2, chi3, chi4)


class NucleicAcidResidue(Residue):
    """A subclass of Residue representing one nuclic acid in a strand of
    DNA or RNA.
    """
    pass


class Atom(object):
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
    Atom.occupancy   - [1.0 - 0.0] float 
    Atom.temp_factor - float represting B-style temp factor
    Atom.U           - a 6-tuple of the anisotropic values
    Atom.charge      - charge on the atom
    """
    def __init__(self,
                 name        = "",
                 alt_loc     = "",
                 res_name    = "",
                 fragment_id = "",
                 chain_id    = ""):

        self.name = name
        self.alt_loc = alt_loc
        self.res_name = res_name
        self.fragment_id = fragment_id
        self.chain_id = chain_id
        self.element = None
        self.position = None
        self.sig_position = None
        self.occupancy = None
        self.sig_occupancy = None
        self.temp_factor = None
        self.sig_temp_factor = None
        self.U = None
        self.sig_U = None
        self.charge = None

        self.__alt_loc_list = WeakrefList()
        self.__alt_loc_list.append(self)
        
        self.bond_list    = []

    def __str__(self):
        return "Atom(%s,%s,%s,%s,%s)" % (self.name,
                                         self.alt_loc,
                                         self.res_name,
                                         self.fragment_id,
                                         self.chain_id)

    def __lt__(self, other):
        assert isinstance(other, Atom)
        return self.alt_loc < other.alt_loc

    def __le__(self, other):
        assert isinstance(other, Atom)
        return self.alt_loc <= other.alt_loc
        
    def __gt__(self, other):
        assert isinstance(other, Atom)
        return self.alt_loc > other.alt_loc

    def __ge__(self, other):
        assert isinstance(other, Atom)
        return self.alt_loc >= other.alt_loc

    def __len__(self):
        return len(self.__alt_loc_list)

    def __getitem__(self, x):
        """This is a alternative to calling get_alt_loc, but a KeyError
        exception is raised if the alt_loc Atom is not found.
        """
        if type(x) == StringType:
            for atom in self.__alt_loc_list:
                if atom.alt_loc == x:
                    return atom
            raise KeyError, x

        raise TypeError, x

    def __delitem__(self, x):
        self.__alt_loc_list.remove(self[x])

    def __iter__(self):
        return iter(self.__alt_loc_list)

    def __contains__(self, x):
        if isinstance(x, Atom):
            return x in self.__alt_loc_list
        elif type(x) == StringType:
            return self[x] in self.__alt_loc_list
        raise TypeError, x

    def add_alt_loc(self, atom):
        assert isinstance(atom, Atom)
        assert atom not in self.__alt_loc_list

        self.__alt_loc_list.append(atom)
        atom.__alt_loc_list = self.__alt_loc_list
        self.__alt_loc_list.sort()

    def get_alt_loc(self, alt_loc):
        """Returns the Atom object matching the alt_loc argument.
        """
        assert type(alt_loc) == StringType

        try:
            return self[alt_loc]
        except KeyError:
            return None

    def iter_alt_loc(self):
        """Iterate over all alt_loc versions of this atom in the
        alphabetical order of the alt_loc labels.  If there are no
        alt_loc versions, do not iterate.
        """
        return iter(self)

    def create_bond(self, atom, bond_alt_loc = False):
        """Creates a bond between two Atom objects.  If bond_alt_loc
        is True, then the bond is also formed between the alternate
        locations of the same atom.
        """
        assert isinstance(atom, Atom)

        def make_bond(a1, a2):
            if a1.get_bond(a2):
                return
            bond = Bond(a1, a2)
            a1.bond_list.append(bond)
            a2.bond_list.append(bond)

        if bond_alt_loc:
            ## this handles constructing the bonds for alternate
            ## conformations correctly
            alist1 = [a.alt_loc for a in self if a.alt_loc]
            alist2 = [a.alt_loc for a in atom if a.alt_loc]
            
            if not (alist1 or alist2):
                make_bond(self, atom)

            elif alist1 and alist2:
                for alt_loc in alist1:
                    make_bond(self[alt_loc], atom[alt_loc])

            elif alist1 and not alist2:
                for alt_loc in alist1:
                    make_bond(self[alt_loc], atom)

            else:
                for alt_loc in alist2:
                    make_bond(self, atom[alt_loc])  

        else:
            make_bond(atom)

    def get_bond(self, atom):
        """Returns the Bond connecting self with the argument atom.
        """
        assert isinstance(atom, Atom)

        for bond in self.bond_list:
            if atom == bond.get_atom1() or atom == bond.get_atom2():
                return bond
        return None

    def iter_bonds(self):
        """Iterates over all the Bond edges connected to self.
        """
        for bond in self.bond_list:
            yield bond

    def iter_bonded_atoms(self):
        """Iterates over all the Atoms bonded to self.
        """
        for bond in self.iter_bonds():
            yield bond.get_partner(self)

    def get_chain(self):
        """Return the parent Chain object.
        """
        return self.get_fragment().get_chain()

    def get_structure(self):
        """Return the parent Structure object.
        """
        return self.get_chain().get_structure()

    def set_U(self, u11, u22, u33, u12, u13, u23):
        """Sets the symmetric U tensor from the six unique values.
        """
        self.U = array([[u11, u12, u13],
                        [u12, u22, u23],
                        [u13, u23, u33]])

    def set_sig_U(self, u11, u22, u33, u12, u13, u23):
        """Sets the symmetric sig_U tensor from the six unique values.
        """
        self.sig_U = array([[u11, u12, u13],
                            [u12, u22, u23],
                            [u13, u23, u33]])

    def calc_anisotropy(self):
        """Calculates the anisotropy of that atom.  Anisotropy is defined
        as the ratio of the minimum/maximum eigenvalues of the 3x3
        symmetric tensor defined by U.
        """
        ## no Anisotropic values, we have a spherical atom
        if not self.U: return 1.0
        evals = eigenvalues(self.U)
        ansotropy = min(evals) / max(evals)
        return ansotropy
        
    def iter_atoms_by_distance(self, max_distance = None):
        """Iterates all atoms in the Structure object from the closest to the
        farthest up to the cutoff distance max_distance if given.  Yields
        the 2-tuple (dist, atm).
        """
        list = []

        if max_distance:
            for atm in self.get_structure().iter_atoms():
                d = calculateDistance(self, atm)
                if d <= max_distance:
                    list.append((calculateDistance(self, atm), atm))
        else:
            for atm in self.get_structure().iter_atoms():
                list.append((calculateDistance(self, atm), atm))

        list.sort()
        return iter(list)


class Bond(object):
    """Indicates two atoms are bonded together.
    """
    def __init__(self, atom1, atom2):
        assert isinstance(atom1, Atom)
        assert isinstance(atom2, Atom)
        assert atom1 != atom2

        self.get_atom1 = weakref.ref(atom1)
        self.get_atom2 = weakref.ref(atom2)

    def __str__(self):
        return "Bond(%s...%s)" % (self.get_atom1(), self.get_atom2())

    def get_partner(self, atm):
        if   atm == self.get_atom1(): return self.get_atom2()
        elif atm == self.get_atom2(): return self.get_atom1()
        return None

    def get_fragment(self):
        return self.get_atom1().get_fragment()

    def get_chain(self):
        return self.get_atom1().get_chain()

    def get_structure(self):
        return self.get_atom1().get_structure()


class FragmentList(list):
    """Provides the functionallity of a Python list class for containing
    Fragment instances.
    """
    def __setitem__(self, x, item):
        assert isinstance(item, Fragment)
        list.__setitem__(self, x, item)

    def append(self, item):
        assert isinstance(item, Fragment)
        list.append(self, item)

    def insert(self, x, item):
        assert isinstance(item, Fragment)
        list.insert(self, x, item)


class Site(FragmentList):
    """List of Fragments within a structure involved in a SITE description.
    """
    def __init__(self, name):
        FragmentList.__init__(self)
        self.name = name


class AlphaHelix(FragmentList):
    """List of Fragments within a structure which are part of a alpha
    helix.
    """
    pass


class BetaSheet(FragmentList):
    """List of Fragments within a structure which are part of a beta
    sheet.
    """
    pass


class AtomList(list):
    """Provides the functionallity of a Python list class for containing
    Atom instances.  It also provides class methods for performing some
    useful calculations on the list of atoms.
    """
    def __setitem__(self, x, item):
        assert isinstance(item, Atom)
        list.__setitem__(self, x, item)

    def append(self, item):
        assert isinstance(item, Atom)
        list.append(self, item)

    def insert(self, x, item):
        assert isinstance(item, Atom)
        list.insert(self, x, item)
    
    def calc_centroid(self):
        """Calculates the centroid of all contained Atom instances and
        returns a Vector to the centroid.
        """
        centroid = Vector(0.0, 0.0, 0.0)
        for atm in self:
            centroid += atm.position
        return centroid / len(self)
        
    def calc_adv_temp_factor(self):
        """Calculates the adverage temperature factor of all contained
        Atom instances and returns the adverage temperature factor.
        """
        adv_tf = 0.0
        for atm in self:
            adv_tf += atm.temp_factor
        return adv_tf / len(self)

    def calc_adv_U(self):
        """Calculates the adverage U matrix of all contained Atom
        instances and returns the 3x3 symmetric U matrix of that
        adverage.
        """
        adv_U = array([[0.0,0.0,0.0],
                       [0.0,0.0,0.0],
                       [0.0,0.0,0.0]])

        for atm in self:
            ## use the atom's U matrix if it exists, otherwise use the
            ## temperature factor
            if atm.U:
                adv_U += atm.U
            else:
                adv_U[0,0] += atm.temp_factor
                adv_U[1,1] += atm.temp_factor
                adv_U[2,2] += atm.temp_factor
        return adv_U / len(self)

    
