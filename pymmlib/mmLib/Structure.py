## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Classes for representing biological macromolecules.
"""
from __future__ import generators
import fpformat
import time
from mmTypes import *
from AtomMath import *
from Library import Library
from UnitCell import UnitCell
from mmCIFDB import *


class Structure(object):
    """The Structure object is the parent container object for the entire
    macromolecular data structure.  It contains a list of the Chain objects
    in the structure hierarchy, and contains these additional data
    objects:

    library(mmLib.Library) The monomer library used by the structure.
                              
    cifdb(mmLib.mmCIFDB) A mmCIF database with additional structure data.

    unit_cell(mmLib.UnitCell) Unit cell/Spacegroup for the structure.

    default_model(int) The default NMR-style model used when iterating
    or retreiving the structure

    default_alt_loc(string) The default alternate location identifier used
    when iterating or retreiving Atom objects in the structure.
    """
    def __init__(self, library = None):
        self.library         = library or Library()
        self.cifdb           = mmCIFDB("XXX")
        self.unit_cell       = UnitCell()
        self.default_model   = 1
        self.default_alt_loc = "A"
        self.chain_list      = []

    def __str__(self):
        tstr = "Structure(id=%s, meth=)" % (self.cifdb.get_entry_id())
        return tstr
 
    def __len__(self):
        """Returns the number of stored Chain objects.
        """
        return len(self.chain_list)
    
    def __getitem__(self, x):
        """Same as get_chain, but raises KeyError if the requested chain_id
        is not found.
        """
        if type(x) == StringType:
            for chain in self.chain_list:
                if chain.chain_id == x:
                    return chain
            raise KeyError, x

        elif type(x) == IntType:
            return self.chain_list[x]

        raise TypeError, x

    def __delitem__(self, x):
        self.chain_list.remove(self[x])

    def __iter__(self):
        """Iterates the Chain objects in the Structure.
        """
        return iter(self.chain_list)

    def __contains__(self, x):
        if isinstance(x, Chain):
            return x in self.chain_list
        elif type(x) == StringType:
            return self[x] in self.chain_list
        raise TypeError, x

    def index(self, chain):
        """Returns the numeric index of the Chain object.
        """
        assert isinstance(chain, Chain)
        return self.chain_list.index(chain)

    def remove(self, chain):
        """Removes the chain from the structure.
        """
        assert isinstance(chain, Chain)
        del chain.structure
        self.chain_list.remove(chain)

    def sort(self):
        self.chain_list.sort()

    def add_chain(self, chain, delay_sort = True):
        """Adds a Chain object to the Structure, and set the chain's
        structure attribute to self.
        """
        assert isinstance(chain, Chain)

        for chn in self:
            if chn.chain_id == chain.chain_id:
                print "[ERROR] add_chain overwrite",chn,chain
                return False

        chain.structure = self
        self.chain_list.append(chain)
        if not delay_sort:
            self.chain_list.sort()
        return True

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
        visited = {}
        for atm in self.iter_atoms():
            for bond in atm.iter_bonds():
                if bond not in visited:
                    yield bond
                    visited[bond] = True

    def iter_alpha_helicies(self):
        """Iterates over all alpha helicies in the Structure.
        """
        try:
            struct_conf = self.cifdb["struct_conf"]
        except KeyError:
            return

        for row in struct_conf.iter_rows(("conf_type_id", "HELX_P")):
            try:
                yield AlphaHelix(self, row["id"])
            except KeyError:
                continue

    def iter_beta_sheets(self):
        """Iterate over all beta sheets in the Structure.
        """
        try:
            struct_sheet = self.cifdb["struct_sheet"]
        except KeyError:
            return

        for row in struct_sheet:
            try:
                yield BetaSheet(self, row["id"])
            except KeyError:
                continue

    def iter_turns(self):
        pass

    def iter_sites(self):
        """Iterate over all active/important sites defined in the Structure.
        """
        try:
            struct_site_gen = self.cifdb["struct_site_gen"]
        except KeyError:
            return

        site_ids = []
        for row in struct_site_gen:
            try:
                site_id = row["site_id"]
            except KeyError:
                continue

            if site_id not in site_ids:
                site_ids.insert(0, site_id)
                yield Site(self, site_id)

    
class Chain(object):
    """Chain objects conatain a ordered list of Fragment objects.
    """
    def __init__(self, chain_id = ""):
        assert type(chain_id) == StringType
        self.chain_id = chain_id

        ## the sequence list contains a list 3-letter residue names
        self.sequence = []

        ## fragments are contained in the list and also cached in
        ## a dictionary for fast random-access lookup
        self.fragment_list  = []

    def __str__(self):
        return "Chain(%s, %s...%s)" % (self.chain_id,
                                       self.fragment_list[0],
                                       self.fragment_list[-1])

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
        return len(self.fragment_list)

    def __getitem__(self, x):
        """Retrieve a Fragment within the Chain.  This can take a integer
        index of the Fragment's position within the chain, the fragment_id
        string of the Fragment to retrieve, or a slice of the Chain to
        return a new Chain object containing the sliced subset of Fragments.
        """
        if type(x) == IntType:
            return self.fragment_list[x]
        
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
            chain           = Chain(self.chain_id)
            chain.structure = self.structure

            frag_start = self[x.start]
            frag_stop  = self[x.stop]

            for frag in self:
                if frag >= frag_start and frag <= frag_stop:
                    chain.fragment_list.append(frag)
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
        self.fragment_list.remove(frag)

    def __iter__(self):
        """Iterate all Fragments contained in the Chain.
        """
        return iter(self.fragment_list)

    def __contains__(self, x):
        if isinstance(x, Fragment):
            return x in self.fragment_list
        elif type(x) == StringType:
            return self[x] in self.fragment_list
        raise TypeError, x

    def index(self, frag):
        """Return the 0-based index of the framgent in the chain list.
        """
        return self.fragment_list.index(frag)

    def remove(self, frag):
        """Remove the Fragment from the chain.
        """
        del frag.chain
        self.fragment_list.remove(frag)

    def sort(self):
        """Sort the Fragments in the chain into proper order.
        """
        self.fragment_list.sort()
        
    def add_fragment(self, frag, delay_sort = False):
        """Adds a Fragment instance to the chain.  If delay_sort is True,
        then the fragment is not inserted in the proper position within the
        chain.
        """
        assert isinstance(frag, Fragment)
        assert frag.chain_id == self.chain_id

        for frg in self:
            if frg.fragment_id == frag.fragment_id:
                print "[ERROR] add_fragment overwrite",frg,frag
                return False

        frag.chain = self
        self.fragment_list.append(frag)
        if not delay_sort:
            self.fragment_list.sort()
        return True

    def get_structure(self):
        """Returns the parent structure.  Also available by self.structure.
        """
        return self.structure

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
        visited = {}
        for atm in self.iter_atoms():
            for bond in atm.iter_bonds():
                if bond not in visited:
                    yield bond
                    visited[bond] = True

    def set_chain_id(self, chain_id):
        """Sets a new ID for the Chain object, updating the chain_id
        for all objects in the Structure hierarchy.
        """
        ## check for conflicting chain_id in the structure
        try:             self.structure[chain_id]
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
        self.structure.sort()

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

        assert type(res_name) == StringType
        assert type(fragment_id) == StringType
        assert type(chain_id) == StringType
        
        self.res_name    = res_name
        self.fragment_id = fragment_id
        self.chain_id    = chain_id

        self.atom_list = []

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
        """Returns the number of atoms contained in the fragment in the
        default model and default alt_loc.
        """
        i = 0
        for atm in self:
            i += 1
        return i

    def __getitem__(self, x):
        """Lookup a atom contained in a fragment by its name, or by its index
        within the fragment's private atom_list.  If the atom is not found,
        a exception is raised.  The type of exception depends on the argument
        type.  If the argument was a integer, then a IndexError is raised.
        If the argument was a string, then a KeyError is raised.
        """
        if type(x) == IntType:
            return self.atom_list[x]

        elif type(x) == StringType:
            for atm in self.atom_list:
                if atm.name == x:
                    ret_atm = atm.get_default()
                    if ret_atm:
                        return ret_atm
                    raise KeyError, x
            raise KeyError, x

        raise TypeError, str(x)

    def __delitem__(self, x):
        """Removes a atom from the fragment.  The argument can be the name
        of the atom, or the private integer index of the atom within the
        fragment's private atom_list.
        """
        self.atom_list.remove(self[x])

    def __iter__(self):
        """Iterates the atoms within the fragment.  If the fragment contains
        atoms in alternate conformations, only the atoms with the structure's
        default_alt_loc are iterated.
        """
        ## this is a difficult function to write and perserve what is left
        ## of iteration performance (not much now that I'm tring to handle
        ## all cases of illegal altLoc combinations)
        ## assume atom_list is properly sorted -- this should be true
        model = self.chain.structure.default_model
        alt_loc = self.chain.structure.default_alt_loc

        yield_me_if_no_match = None

        for atm in self.atom_list:
            if atm.model != model:
                continue
            if yield_me_if_no_match:
                if yield_me_if_no_match.name == atm.name:
                    if atm.alt_loc == alt_loc:
                        yield_me_if_no_match = None
                        yield atm
                else:
                    yield yield_me_if_no_match
                    yield_me_if_no_match = None
                    if atm.alt_loc == alt_loc:
                        yield atm
                    elif not atm.alt_loc:
                        yield_me_if_no_match = atm
            else:
                if atm.alt_loc == alt_loc:
                    yield atm
                elif not atm.alt_loc:
                    yield_me_if_no_match = atm

        if yield_me_if_no_match:
            yield yield_me_if_no_match

    def __contains__(self, x):
        """Returns True if the atom is contained in the fragment.  The argument
        can be a atom instance of a the name of the atom.
        """
        if isinstance(x, Atom):
            return x in self.atom_list
        return self[x] in self.atom_list

    def index(self, atom):
        """Returns the index of the atom within the fragment's private
        atom_list.
        """
        assert isinstance(atom, Atom)
        return self.atom_list.index(atom)

    def remove(self, atom):
        """Removes the atom from the fragment, and deletes the atom's
        reference to the framgnet.
        """
        assert isinstance(atom, Atom)
        del atom.fragment
        atom.alt_loc_list.remove(atom)
        self.atom_list.remove(atom)

    def add_atom(self, atom):
        """Adds a atom to the fragment, and sets the atom's atom.fragment
        attribute to the fragment.
        """
        assert isinstance(atom, Atom)
        assert atom.chain_id    == self.chain_id
        assert atom.fragment_id == self.fragment_id
        assert atom not in self.atom_list

        for atm in self.atom_list:
            if atm.name == atom.name and \
               atm.model == atom.model and \
               atm.alt_loc == atom.alt_loc:
                debug("[ERROR] add_atom overwrite " + str(atm) + str(atom))
                return False

        atom.fragment = self

        ## if this atom is a alternate conformation of a atom already in
        ## the fragment, add it to that atom's private __alt_loc_list
        for atm in self.atom_list:
            if atm.name == atom.name:
                atm.add_alt_loc(atom)
                break

        ## now add the atom to the fragment
        self.atom_list.append(atom)
        return True

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
        visited = {}
        for atm in self.iter_atoms():
            for bond in atm.iter_bonds():
                if bond not in visited:
                    yield bond
                    visited[bond] = True

    def get_offset_fragment(self, offset):
        """Returns the fragment in the same chain at integer offset from
        self.  Returns None if no fragment is found.
        """
        assert type(offset) == IntType

        i = self.chain.index(self) + offset
        if i < 0:
            return None
        try:
            return self.chain[i]
        except IndexError:
            return None

    def get_structure(self):
        """Returns the parent structure.  This is also available by the
        attribute self.chain.structure.
        """
        return self.chain.structure

    def get_chain(self):
        """Returns the parent chain, this is also available by the attribute
        self.chain.
        """
        return self.chain

    def set_fragment_id(self, fragment_id):
        """Sets a new ID for the Fragment object, updating the fragment_id
        for all objects in the Structure hierarchy.
        """
        ## check for conflicting chain_id in the structure
        try:             self.chain[fragment_id]
        except KeyError: pass
        else:            raise ValueError, fragment_id

        ## set the new fragment_id in all the additional groups

        ## set the new chain_id for the chain object (self)
        self.fragment_id = fragment_id

        ## set the chain_id in all the fragment and atom children
        for atm in self:
            atm.fragment_id = fragment_id

        ## resort the parent chain
        self.chain.sort()

    def create_bonds(self):
        """Contructs bonds within a fragment.  Bond definitions are retrieved
        from the monomer library.
        """
        mon = self.chain.structure.library.get_monomer(self.res_name)
        if not mon:
            return

        for (name1, name2) in mon.bond_list:
            try:
                atm1 = self[name1]
                atm2 = self[name2]
            except KeyError:
                continue
            else:
                atm1.create_bonds(atm2, standard_res_bond = True)

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
        if self.chain.structure.library.is_water(self.res_name):
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

        library = self.chain.structure.library
        mon1 = library.get_monomer(self.res_name)
        mon2 = library.get_monomer(next_res.res_name)
        if not (mon1 and mon2):
            return

        for (name1, name2) in mon1.get_polymer_bond_list(self, next_res):
            try:
                atm1 = self[name1]
                atm2 = next_res[name2]
            except KeyError:
                continue
            else:
                atm1.create_bonds(atm2, standard_res_bond = True)


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
        """Returns True if this is a CIS amino acid, otherwise returns False.
        It uses calc_torsion_omega, and if there are missing atoms this method
        will return None.
        """
        prev_res = self.get_offset_residue(-1)
        if prev_res == None:
            return None

        prev_omega = prev_res.calc_torsion_omega()
        if prev_omega == None:
            return None

        if abs(prev_omega) <= (math.pi/2.0):
            return True

        return False

    def calc_pucker_torsion(self, conf_id = None):
        """Calculates the Pucker torsion of a ring system.  Returns None
        for Amino Acids which do not have Pucker torsion angles.
        """
        mon = self.chain.structure.library.get_monomer(self.res_name)
        if not mon or not mon.pucker_definition:
            return None

        a1 = self.get_atom(mon.pucker_definition[0])
        a2 = self.get_atom(mon.pucker_definition[1])
        a3 = self.get_atom(mon.pucker_definition[2])
        a4 = self.get_atom(mon.pucker_definition[3])
        return calculateTorsionAngle(a1, a2, a3, a4)

    def calc_torsion_chi1(self):
        mon = self.chain.structure.library.get_monomer(self.res_name)
        if not mon.chi1_definition:
            return None
        
        a1 = self.get_atom(mon.chi1_definition[0])
        a2 = self.get_atom(mon.chi1_definition[1])
        a3 = self.get_atom(mon.chi1_definition[2])
        a4 = self.get_atom(mon.chi1_definition[3])
        return calculateTorsionAngle(a1, a2, a3, a4)

    def calc_torsion_chi2(self):
        mon = self.chain.structure.library.get_monomer(self.res_name)
        if not mon.chi2_definition:
            return None
        
        a1 = self.get_atom(mon.chi2_definition[0])
        a2 = self.get_atom(mon.chi2_definition[1])
        a3 = self.get_atom(mon.chi2_definition[2])
        a4 = self.get_atom(mon.chi2_definition[3])
        return calculateTorsionAngle(a1, a2, a3, a4)

    def calc_torsion_chi3(self):
        mon = self.chain.structure.library.get_monomer(self.res_name)
        if not mon.chi3_definition:
            return None
        
        a1 = self.get_atom(mon.chi3_definition[0])
        a2 = self.get_atom(mon.chi3_definition[1])
        a3 = self.get_atom(mon.chi3_definition[2])
        a4 = self.get_atom(mon.chi3_definition[3])
        return calculateTorsionAngle(a1, a2, a3, a4)

    def calc_torsion_chi4(self):
        mon = self.chain.structure.library.get_monomer(self.res_name)
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
                 model       = 0,
                 alt_loc     = "",
                 res_name    = "",
                 fragment_id = "",
                 chain_id    = ""):

        assert type(name) == StringType
        assert type(model) == IntType
        assert type(alt_loc) == StringType
        assert type(res_name) == StringType
        assert type(fragment_id) == StringType
        assert type(chain_id) == StringType

        self.name = name
        self.model = model
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

        self.alt_list = [self]
        self.bond_list = []

    def __str__(self):
        return 'Atom(n=%s,mdl=%d,alt=%s,rn=%s,fid=%s,cid=%s)' % (
            self.name, self.model, self.alt_loc, self.res_name,
            self.fragment_id, self.chain_id)
    
    def __lt__(self, other):
        assert isinstance(other, Atom)
        if self.model == other.model:
            return self.alt_loc < other.alt_loc
        else:
            return self.model < other.model

    def __le__(self, other):
        assert isinstance(other, Atom)
        if self.model == other.model:
            return self.alt_loc <= other.alt_loc
        else:
            return self.model <= other.model
        
    def __gt__(self, other):
        assert isinstance(other, Atom)
        if self.model == other.model:
            return self.alt_loc > other.alt_loc
        else:
            return self.model > other.model

    def __ge__(self, other):
        assert isinstance(other, Atom)
        if self.model == other.model:
            return self.alt_loc >= other.alt_loc
        else:
            return self.model >= other.model

    def __len__(self):
        """Returns the number of alternate conformations of this atom.
        """
        return len(self.alt_list)

    def __getitem__(self, x):
        """This is a alternative to calling get_alt_loc, but a KeyError
        exception is raised if the alt_loc Atom is not found.  Posiible
        arguments are:
        """
        default_model = self.fragment.chain.structure.default_model
        default_alt_loc = self.fragment.chain.structure.default_alt_loc
        
        if type(x) == IntType:
            default_atm = None
            for atm in self.alt_list:
                if atm.model != x:
                    continue
                if atm.alt_loc == default_alt_loc:
                    return atm
                elif not atm.alt_loc:
                    default_atm = atm
            if default_atm:
                return default_atm
            raise KeyError, x

        elif type(x) == StringType:
            for atm in self.alt_list:
                if atm.model == default_model and atm.alt_loc == x:
                    return atm
            raise KeyError, x

        elif type(x) == TupleType and \
             len(x) == 2 and \
             type(x[0]) == StringType and \
             type(x[1]) == IntType:
            for atm in self.alt_list:
                if atm.alt_loc == x[0] and atm.model == x[1]:
                    return atm
            raise KeyError, x

        raise TypeError, x

    def __delitem__(self, x):
        """Deletes the alternate conformation of the atom matching the
        argument.  The argument is the same as the argument for __getitem__.
        """
        self.remove(self[x])

    def __iter__(self):
        """Iterate all alternate locations of this atom within the default
        model.
        """
        default_model = self.fragment.chain.structure.default_model
        for atm in self.alt_list:
            if atm.model == default_model:
                yield atm

    def __contains__(self, x):
        """Returns True if the argument matches a alternate conformation of
        the atom.  The argument can be a alt_loc label, or a Atom object.
        """
        if isinstance(x, Atom):
            return x in self.alt_list
        return self[x] in self.alt_list

    def remove(self, atom):
        """Remove the atom from the list of alternate conformations of this
        atom.
        """
        atom.alt_list.remove(atm)
        atom.alt_list = []
        atom.fragment.remove(atom)
        del atom.fragment

    def add_alt_loc(self, atom):
        """Add atom as a alternate conformation of the current atom.
        """
        assert isinstance(atom, Atom)
        assert atom.name == self.name
        assert atom.fragment_id == self.fragment_id
        assert atom.chain_id == self.chain_id
        assert atom not in self.alt_list

        for atm in self.alt_list:
            if atm.model == atom.model and \
               atm.alt_loc == atom.alt_loc:
                print "[ERROR] add_alt_loc overwrite",atm,atom
                return False

        self.alt_list.append(atom)
        atom.alt_list = self.alt_list
        self.alt_list.sort()
        return True

    def get_default(self):
        """Returns the atom matching the Structure's default_model
        and defualt_alt_loc.  Returns None if there is no default atom
        witht he current default settings.
        """
        default_model = self.fragment.chain.structure.default_model
        default_alt_loc = self.fragment.chain.structure.default_alt_loc
        default_atm = None
        for atm in self.alt_list:
            if atm.model != default_model:
                continue
            if atm.alt_loc == default_alt_loc:
                return atm
            elif not atm.alt_loc:
                default_atm = atm
        if default_atm:
            return default_atm
        return None

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
        alphabetical order of the alt_loc labels.
        """
        return iter(self)

    def get_model(self, model):
        """Returns the atom in the argument model number.  Uses the Structure
        default_alt_loc.  Return None if the atom is not found.
        """
        assert type(model) == IntType
        try:
            return self[model]
        except KeyError:
            return None

    def iter_model(self):
        """Iterates over all models of this atom matching the structure
        wide default_alt_loc.
        """
        default_alt_loc = self.fragment.chain.structure.default_alt_loc
        for atm in self.alt_list:
            if not atm.alt_loc or atm.alt_loc == default_alt_loc:
                yield atm

    def create_bond(self,
                    atom,
                    bond_type = None,
                    atm1_symop = None,
                    atm2_symop = None,
                    standard_res_bond = False):
        """Creates a bond between this atom and the argumentatom.  The
        arugment bond_type is a string, atm1_symop and atm2_symop are
        symmetry operations to be applied to self and the argument atom
        before distance calculations, and standard_res_bond is a flag
        used to indicate this bond is a standard bond.
        """
        assert isinstance(atom, Atom)
        assert self.model == atom.model
        assert ((self.alt_loc == atom.alt_loc) or
                (self.alt_loc == "" and atom.alt_loc != "") or
                (self.alt_loc != "" and atom.alt_loc == ""))
            
        bond = Bond(self, atom,
                    bond_type = bond_type,
                    atm1_symop = atm1_symop,
                    atm2_symop = atm2_symop,
                    standard_res_bond = standard_res_bond)

        self.bond_list.append(bond)
        atom.bond_list.append(bond)

    def create_bonds(self,
                     atom,
                     bond_type = None,
                     atm1_symop = None,
                     atm2_symop = None,
                     standard_res_bond = False):
        """Like create_bonds, but it bonds all alternate locations and models
        of this atom.
        """
        assert isinstance(atom, Atom)

        for atm1 in self.alt_list:
            for atm2 in atom.alt_list:
                if atm1.model != atm2.model:
                    continue

                if not atm1.alt_loc or not atm2.alt_loc or \
                   atm1.alt_loc == atm2.alt_loc:
                    atm1.create_bond(atm2,
                                     bond_type=bond_type,
                                     atm1_symop=atm1_symop,
                                     atm2_symop=atm2_symop,
                                     standard_res_bond=standard_res_bond)

    def get_bond(self, atom):
        """Returns the Bond connecting self with the argument atom.
        """
        assert isinstance(atom, Atom)

        for bond in self.bond_list:
            if atom == bond.atom1 or atom == bond.atom2:
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

    def get_fragment(self):
        """Return the parent Fragment object.  This is also available by the
        attribute self.fragment.
        """
        return self.fragment

    def get_chain(self):
        """Return the parent Chain object.  This is also available by the
        attribute self.fragment.chain.
        """
        return self.fragment.chain

    def get_structure(self):
        """Return the parent Structure object.  This is also abailable by the
        attribute self.fragment.chain.structure
        """
        return self.fragment.chain.structure

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
    def __init__(self, atom1, atom2,
                 bond_type = None,
                 atm1_symop = None,
                 atm2_symop = None,
                 standard_res_bond = False):

        assert isinstance(atom1, Atom)
        assert isinstance(atom2, Atom)
        assert atom1 != atom2

        self.atom1 = atom1
        self.atom2 = atom2
        self.bond_type = bond_type
        self.atm1_symop = atm1_symop
        self.atm2_symop = atm2_symop
        self.standard_res_bond = standard_res_bond

    def __str__(self):
        return "Bond(%s...%s)" % (self.atom1, self.atom2)

    def get_partner(self, atm):
        if   atm == self.atom1: return self.atom2
        elif atm == self.atom2: return self.atom1
        return None

    def get_atom1(self):
        """Returns atom #1 of the pair of bonded atoms.  This is also
        accessable by bond.atom1.
        """
        return self.atom1

    def get_atom2(self):
        """Returns atom #2 of the pair of bonded atoms.
        """
        return self.atom2


class AlphaHelix:
    def __init__(self, structure, helix_id):
        self.structure = structure
        self.helix_id = helix_id

    def __str__(self):
        try:
            (start_frag, end_frag) = self.get_start_end_fragments()
        except KeyError:
            return "AlphaHelix(id=%s,start=*,end=*)" % (self.helix_id)
        return "AlphaHelix(id=%s,start=%s,end=%s)" % (
            self.helix_id, start_frag, end_frag)

    def get_start_end_fragments(self):
        """Return the tuple of the fragment at the start of the
        alpha helix, and the fragment at the end.  Raises KeyError
        if either fragment cannot be found.
        """
        struct_conf = self.structure.cifdb["struct_conf"]
        row = struct_conf.get_row(("id", self.helix_id))
        
        chain_id1 = row["beg_label_asym_id"]
        frag_id1 = row["beg_label_seq_id"]
        chain_id2 = row["end_label_asym_id"]
        frag_id2 = row["end_label_seq_id"]
        
        frag1 = self.structure[chain_id1][frag_id1]
        frag2 = self.structure[chain_id2][frag_id2]

        return (frag1, frag2)

    def iter_fragments(self):
        """Iterate the Fragment objects contained in the AlphaHelix.
        """
        try:
            (start_frag, end_frag) = self.get_start_end_fragments()
        except KeyError:
            return
        frag_iter = self.structure.iter_fragments()
        for frag in frag_iter:
            if frag == start_frag:
                yield frag
                break
        for frag in frag_iter:
            yield frag
            if frag == end_frag:
                break

    def iter_atoms(self):
        """Iterate all Atoms in the AlphaHelix.
        """
        for frag in self.iter_fragments():
            for atm in frag.iter_atoms():
                yield atm


class BetaSheet:
    """List of Fragments within a structure which are part of a beta
    sheet.
    """
    def __init__(self, structure, sheet_id):
        self.structure = structure
        self.sheet_id = sheet_id

    def __str__(self):
        return "BetaSheet(id=%s)" % (self.sheet_id)


class Site:
    """List of Fragments within a structure involved in a SITE description.
    """
    def __init__(self, structure, site_id):
        self.structure = structure
        self.site_id = site_id

    def __str__(self):
        return "Site(id=%s)" % (self.site_id)

    def iter_fragments(self):
        struct_site_gen = self.structure.cifdb["struct_site_gen"]
        for row in struct_site_gen.iter_rows(("site_id", self.site_id)):
            try:
                chain_id = row.mget("_asym_id", "label_asym_id")
                frag_id = row.mget("auth_seq_id","label_seq_id")
                yield self.structure[chain_id][frag_id]
            except KeyError:
                continue


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

    
