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


class StructureError(Exception):
    """Base class of errors raised by Structure objects.
    """
    pass

class ModelOverwrite(StructureError):
    """Raised by Structure.add_model() when a Model added to a Structure
    has the same model_id of a Model already in the Structure.
    """
    pass

class ChainOverwrite(StructureError):
    """Raised by Structure.add_chain() or by Model.add_chain() when a
    Chain added to a Structure has the same chain_id of a Chain already
    in the Structure.
    """
    pass

class FragmentOverwrite(StructureError):
    """Raised by Chain.add_fragment() when a Fragment added to a Chain
    has the same fragment_id as a Fragment already in the Chain.
    """
    pass

class AtomOverwrite(StructureError):
    """Raised by Structure.add_atom() or Fragment.add_atom() when a Atom
    added to a Structure or Fragment has the same chain_id, fragment_id,
    name, and alt_loc as a Atom already in the Structure or Fragment.
    """
    def __init__(self, text):
        self.text = text
    def __str__(self):
        return self.text


class Structure(object):
    """The Structure object is the parent container object for the entire
    macromolecular data structure.  It contains a list of the Chain objects
    in the structure hierarchy, and contains these additional data
    objects:

    library(mmLib.Library) The monomer library used by the structure.
                              
    cifdb(mmLib.mmCIFDB) A mmCIF database with additional structure data.

    unit_cell(mmLib.UnitCell) Unit cell/Spacegroup for the structure.

    default_alt_loc(string) The default alternate location identifier used
    when iterating or retreiving Atom objects in the structure.
    """
    def __init__(self, **args):
        self.library         = args.get("library")   or Library()
        self.cifdb           = args.get("cifdb")     or mmCIFDB("XXXX")
        self.unit_cell       = args.get("unit_cell") or UnitCell()

        self.default_alt_loc = "A"
        self.model           = None
        self.model_list      = []
        self.model_dict      = {}

    def __str__(self):
        return "Struct(%s)" % (self.cifdb.get_entry_id())
 
    def __deepcopy__(self, memo):
        structure = Structure(
            library   = self.library,
            cifdb     = copy.deepcopy(self.cifdb, memo),
            unit_cell = copy.deepcopy(self.unit_cell, memo))

        for model in self.iter_models():
            structure.add_model(copy.deepcopy(model, memo))

        return structure

    def __len__(self):
        """Returns the number of stored Chain objects.
        """
        try:
            return len(self.model)
        except TypeError:
            return 0
    
    def __getitem__(self, chain_idx):
        """Same as get_chain, but raises KeyError if the requested chain_id
        is not found.
        """
        try:
            return self.model[chain_idx]
        except TypeError:
            raise KeyError, chain_idx

    def __delitem__(self, chain_idx):
        """Removes the Chain from the default Model by its chain_id or
        index in the Model.
        """
        try:
            self.model.remove(self[chain_idx])
        except AttributeError:
            raise KeyError, chain_idx
            
    def __iter__(self):
        """Iterates the Chain objects in the Structure.
        """
        try:
            return iter(self.model)
        except TypeError:
            return iter([])

    def __contains__(self, model_chain_idx):
        """Returns True if item is a Model in the Structure, or a
        Chain or chain_id in the default Model.
        """
        if isinstance(model_chain_idx, Model):
            return self.model_list.__contains__(model_chain_idx)
        elif isinstance(model_chain_idx, Chain) or \
             type(model_chain_idx) == StringType:
            try:
                return self.model.__contains__(model_chain_idx)
            except AttributeError:
                raise KeyError, model_chain_idx
        raise TypeError, model_chain_idx

    def index(self, model_chain):
        """If item is a Model, returns the index of the Model in the
        Structure, or if the item is a Chain, returns the index of the
        Chain in the default Model.
        """
        if isinstance(model_chain, Model):
            return self.model_list.index(model_chain)
        elif isinstance(model_chain, Chain):
            try:
                return self.model.index(model_chain)
            except AttributeError:
                raise ValueError, model_chain
        raise TypeError, model_chain

    def remove(self, model_chain):
        """Removes a Model or a default Model's Chain from the Structure.
        """
        if isinstance(model_chain, Model):
            self.model_list.remove(model_chain)
            del self.model_dict[model_chain.model_id]
            del model_chain.structure
            for chain in model_chain.iter_chains():
                del chain.structure
        elif isinstance(model_chain, Chain):
            try:
                self.model.remove(model_chain)
            except AttributeError:
                raise ValueError, model_chain
        raise TypeError, model_chain

    def sort(self):
        """Sorts all Models and Chains in the Structure occording to standard
        model_id, chain_id, and fragment_id sorting rules.
        """
        self.model_list.sort()
        for model in self.model_list:
            model.sort()

    def add_model(self, model, delay_sort = True):
        """Adds a Model to a Structure.  Raises the ModelOverwrite exception
        if the model_id of the Model matches the model_id of a Model
        already in the Structure.  If there are no Models in the Structure,
        the Model is used as the default Model.
        """
        assert isinstance(model, Model)

        if self.model_dict.has_key(model.model_id):
            raise ModelOverwrite()

        if self.model == None:
            self.model = model

        self.model_list.append(model)
        self.model_dict[model.model_id] = model

        model.structure = self
        for chain in model.iter_chains():
            chain.structure = self

    def add_chain(self, chain, delay_sort = True):
        """Adds a Chain object to the Structure.
        """
        assert isinstance(chain, Chain)

        try:
            model = self.model_dict[chain.model_id]
        except KeyError:
            model = Model(model_id = chain.model_id)
            self.add_model(model)
            model.add_chain(chain)

        if delay_sort == False:
            self.model.chain_list.sort()

    def add_atom(self, atom):
        """
        """
        assert isinstance(atom, Atom)

        ## add new model if necesary
        try:
            model = self.model_dict[atom.model_id]
        except KeyError:
            model = Model(model_id = atom.model_id)
            self.add_model(model, delay_sort = True)

        ## add new chain if necessary
        try:
            chain = model.chain_dict[atom.chain_id]
        except KeyError:
            chain = Chain(model_id = atom.model_id, chain_id = atom.chain_id)
            model.add_chain(chain, delay_sort = True)

        ## add new fragment if necessary 
        try:
            frag = chain[atom.fragment_id]
        except KeyError:
            if self.library.is_amino_acid(atom.res_name):
                frag = AminoAcidResidue(
                    res_name = atom.res_name,
                    fragment_id = atom.fragment_id,
                    chain_id = atom.chain_id)
            elif self.library.is_nucleic_acid(atom.res_name):
                frag = NucleicAcidResidue(
                    res_name = atom.res_name,
                    fragment_id = atom.fragment_id,
                    chain_id = atom.chain_id)
            else:
                frag = Fragment(
                    res_name = atom.res_name,
                    fragment_id = atom.fragment_id,
                    chain_id = atom.chain_id)

            chain.add_fragment(frag, delay_sort = True)
        else:
            if frag.res_name != atom.res_name:
                raise FragmentOverwrite()

        frag.add_atom(atom)

    def get_structure(self):
        """Returns self.
        """
        return self

    def get_chain(self, chain_id):
        """Returns the Chain object matching the chain_id charactor.
        """
        try:
            return self[chain_id]
        except KeyError:
            return None

    def set_model(self, model_id):
        """Sets the default Model for the Structure to model_id.  Returns
        False if a Model with the proper model_id does
        not exist in the Structure.
        """
        try:
            self.model = self.model_dict[model_id]
        except KeyError:
            return False
        return False

    def set_alt_loc(self, alt_loc):
        """Sets the default alt_loc for the Stucture.
        """
        self.default_alt_loc = alt_loc
        for frag in self.iter_fragments():
            frag.set_alt_loc(alt_loc)

    def iter_models(self):
        """Iterates over all Model objects.
        """
        return iter(self.model_list)

    def iter_chains(self):
        """Iterates over all Chain objects in the current Model, in
        alphabetical order according to their chain_id.
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

    def iter_atoms(self, **args):
        """Iterates over all Atom objects in the current Model, using the
        default alt_loc.  The iteration is preformed in order according to
        the Chain and Fragment ordering rules the Atom object is a part of.
        """
        for chain in self.iter_chains():
            for atm in chain.iter_atoms(**args):
                yield atm

    def iter_all_atoms(self):
        """Iterates over all Atom objects in the Structure.  The iteration
        is performed according to common PDB ordering rules, over all Models
        and all Altlocs.
        """
        for model in self.iter_models():
            for frag in model.iter_fragments():
                for atm in frag.iter_all_atoms():
                    yield atm

    def iter_bonds(self):
        """Iterates over all Bond objects.  The iteration is preformed by
        iterating over all Atom objects in the same order as iter_atoms(),
        then iterating over each Atom's Bond objects.
        """
        visited = {}
        for atm in self.iter_atoms():
            for bond in atm.iter_bonds():
                if visited.has_key(bond):
                    continue
                yield bond
                visited[bond] = True

    def alt_loc_list(self):
        """Return the unique list of Atom alternate location IDs found in
        the Structure.
        """
        al_list = []
        for atm in self.iter_all_atoms():
            if atm.alt_loc != "" and atm.alt_loc not in al_list:
                al_list.append(atm.alt_loc)
        return al_list

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

    def add_bonds_from_library(self):
        """Builds bonds for all Fragments in the Structure from bond
        tables for monomers retrieved from the Library implementation
        of the Structure.
        """
        for model in self.iter_models():
            for frag in self.iter_fragments():
                frag.create_bonds()


class Model(object):
    """Multiple models support.
    """
    def __init__(self, model_id = 1):
        assert type(model_id) == IntType

        self.model_id   = model_id
        self.chain_dict = {}
        self.chain_list = []

    def __str__(self):
        return "Model(model_id=%d)" % (self.model_id)

    def __deepcopy__(self, memo):
        model = Model(model_id = self.model_id)
        for chain in self:
            model.add_chain(copy.deepcopy(chain, memo))
        return model

    def __lt__(self, other):
        assert isinstance(other, Model)
        return int(self.model_id) < int(other.model_id)
        
    def __le__(self, other):
        assert isinstance(other, Model)
        return int(self.model_id) <= int(other.model_id)
        
    def __gt__(self, other):
        assert isinstance(other, Model)
        return int(self.model_id) > int(other.model_id)

    def __ge__(self, other):
        assert isinstance(other, Model)
        return int(self.model_id) >= int(other.model_id)

    def __len__(self):
        """Returns the number of stored Chain objects.
        """
        return len(self.chain_list)
    
    def __getitem__(self, chain_idx):
        """Same as get_chain, but raises KeyError if the requested chain_id
        is not found.
        """
        if type(chain_idx) == StringType:
            return self.chain_dict[chain_idx]
        elif type(chain_idx) == IntType:
            return self.chain_list[chain_idx]
        raise TypeError, chain_idx

    def __delitem__(self, chain_idx):
        """Removes a Chain from the Model,, given the chain_id or index of
        the Chain in the Model.
        """
        self.remove(self[chain_idx])

    def __iter__(self):
        """Iterates the Chain objects in the Model.
        """
        return iter(self.chain_list)

    def __contains__(self, chain_idx):
        """Returns True if the argument Chain or chain_id is in the Model.
        """
        if isinstance(chain_idx, Chain):
            return self.chain_list.__contains__(chain_idx)
        elif type(chain_idx) == StringType:
            return self.chain_dict.__contains__(chain_idx)
        raise TypeError, x

    def index(self, chain):
        """Returns the numeric index of the Chain object in the Model.
        """
        assert isinstance(chain, Chain)
        return self.chain_list.index(chain)

    def remove(self, chain):
        """Removes the Chain from the Model.
        """
        assert isinstance(chain, Chain)

        self.chain_list.remove(chain)
        del self.chain_dict[chain.chain_id]
        del chain.model
        try:
            del chain.structure
        except AttributeError:
            pass

    def sort(self):
        """Sorts all Chains in the Model by their chain_id.
        """
        self.chain_list.sort()
        for chain in self.chain_list:
            chain.sort()

    def add_chain(self, chain, delay_sort = False):
        """Adds a Chain to the Model.
        """
        assert isinstance(chain, Chain)

        if self.chain_dict.has_key(chain.chain_id):
            raise ChainOverwrite()

        self.chain_list.append(chain)
        self.chain_dict[chain.chain_id] = chain
        chain.model = self
        try:
            chain.structure = self.structure
        except AttributeError:
            pass

        if delay_sort == False:
            self.chain_list.sort()
            
    def get_chain(self, chain_id):
        """Returns the Chain object matching the chain_id charactor.
        """
        try:
            return self[chain_id]
        except KeyError:
            return None

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

    def iter_atoms(self, **args):
        """Iterates over all Atom objects.  The iteration is preformed in
        order according to the Chain and Fragment ordering rules the Atom
        object is a part of.
        """
        for chain in self.iter_chains():
            for atm in chain.iter_atoms(**args):
                yield atm

    def iter_bonds(self):
        """Iterates over all Bond objects.  The iteration is preformed by
        iterating over all Atom objects in the same order as iter_atoms(),
        then iterating over each Atom's Bond objects.
        """
        visited = {}
        for atm in self.iter_atoms():
            for bond in atm.iter_bonds():
                if visited.has_key(bond):
                    continue
                yield bond
                visited[bond] = True

    
class Chain(object):
    """Chain objects conatain a ordered list of Fragment objects.
    """
    def __init__(self, model_id = 1, chain_id = ""):
        assert type(model_id) == IntType
        assert type(chain_id) == StringType

        self.model_id = model_id
        self.chain_id = chain_id

        ## the sequence list contains a list 3-letter residue names
        self.sequence = None

        ## fragments are contained in the list and also cached in
        ## a dictionary for fast random-access lookup
        self.fragment_list  = []
        self.fragment_dict  = {}

    def __str__(self):
        try:
            return "Chain(%d:%s, %s...%s)" % (
                self.model_id, self.chain_id,
                self.fragment_list[0],
                self.fragment_list[-1])
        except IndexError:
             return "Chain(%d:%s)" % (self.model_id, self.chain_id)

    def __deepcopy__(self, memo):
        chain = Chain(model_id = self.model_id,
                      chain_id = self.chain_id)
        for fragment in self:
            chain.add_fragment(copy.deepcopy(fragment, memo))
        return chain
    
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

    def __getitem__(self, fragment_idx):
        """Retrieve a Fragment within the Chain.  This can take a integer
        index of the Fragment's position within the chain, the fragment_id
        string of the Fragment to retrieve, or a slice of the Chain to
        return a new Chain object containing the sliced subset of Fragments.
        """
        if type(fragment_idx) == IntType:
            return self.fragment_list[fragment_idx]
        elif type(fragment_idx) == StringType:
            return self.fragment_dict[fragment_idx]
        raise TypeError, fragment_idx

    def __delitem__(self, fragment_idx):
        """Delete Fraagment from the chain.  This can take a reference to the
        Fragment object to delete, the fragment_id of the Fragment to delete,
        or the integer index of the Fragment within the Chain.
        """
        self.remove(self[fragment_idx])

    def __iter__(self):
        """Iterate all Fragments contained in the Chain.
        """
        return iter(self.fragment_list)

    def __contains__(self, fragment_idx):
        if isinstance(fragment_idx, Fragment):
            return self.fragment_list.__contains__(fragment_idx)
        elif type(fragment_idx) == StringType:
            return self.fragment_dict.__contains__(fragment_idx)
        raise TypeError, fragment_idx

    def index(self, fragment):
        """Return the 0-based index of the framgent in the chain list.
        """
        return self.fragment_list.index(fragment)

    def remove(self, fragment):
        """Remove the Fragment from the chain.
        """
        self.fragment_list.remove(frag)
        del self.fragment_dict[fragment.fragment_id]
        del fragment.chain

    def sort(self):
        """Sort the Fragments in the chain into proper order.
        """
        self.fragment_list.sort()
        
    def add_fragment(self, fragment, delay_sort = False):
        """Adds a Fragment instance to the chain.  If delay_sort is True,
        then the fragment is not inserted in the proper position within the
        chain.
        """
        assert isinstance(fragment, Fragment)
        assert fragment.chain_id == self.chain_id

        if self.fragment_dict.has_key(fragment.fragment_id):
            raise FragmentOverwrite()

        self.fragment_list.append(fragment)
        self.fragment_dict[fragment.fragment_id] = fragment
        fragment.chain = self
        
        if delay_sort == False:
            self.fragment_list.sort()

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

    def iter_atoms(self, **args):
        """Iterates over all Atom objects within the Chain.
        """
        for frag in self.iter_fragments():
            for atm in frag.iter_atoms(**args):
                yield atm
                
    def iter_bonds(self):
        """Iterates over all Bond objects attached to Atom objects within the
        Chain.
        """
        visited = {}
        for atm in self.iter_atoms():
            for bond in atm.iter_bonds():
                if visited.has_key(bond):
                    continue
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
        structure = self.get_structure()
        if structure != None:
            sequence = Sequence(library = structure.library)
        else:
            sequence = Sequence()

        residue_class = None
        for frag in self.iter_standard_residues():
            if residue_class:
                if not isinstance(frag, residue_class):
                    break
            else:
                residue_class = frag.__class__
            sequence.append(frag.res_name)

        return sequence


class Sequence(list):
    """A polymer sequence 
    """
    def __init__(self, **args):
        self.library  = args.get("library")       or Library()
        sequence_list = args.get("sequence_list") or []
        list.__init__(self, sequence_list)

    def __str__(self):
        return self.sequence_one_letter_code()

    def sequence_one_letter_code(self):
        """Return the one letter code representation of the sequence as
        a string.
        """
        one_letter_code = ""
        
        for res_name in self:
            mon = self.library.get_monomer(res_name)

            if mon == None or not mon.is_standard_residue():
                break

            if mon.one_letter_code == "":
                one_letter_code += "(%s)" % (res_name)
            else:
                one_letter_code += mon.one_letter_code

        return one_letter_code


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

        assert type(res_name)    == StringType
        assert type(fragment_id) == StringType
        assert type(chain_id)    == StringType

        self.res_name = res_name
        self.fragment_id = fragment_id
        self.chain_id = chain_id

        self.default_alt_loc = "A"

        self.atom_order_list = []
        self.alt_loc_dict    = {}

        self.atom_list       = []
        self.atom_dict       = {}

    def __str__(self):
        return "Frag(%s,%s,%s)" % (
            self.res_name,
            self.fragment_id,
            self.chain_id)
    
    def __deepcopy__(self, memo):
        fragment = Fragment(
            res_name    = self.res_name,
            fragment_id = self.fragment_id,
            chain_id    = self.chain_id)

        for atom in self.iter_all_atoms():
            fragment.add_atom(copy.deepcopy(atom, memo))

        return fragment
    
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
        return len(self.atom_list)
    
    def __getitem__(self, name_idx):
        """Lookup a atom contained in a fragment by its name, or by its index
        within the fragment's private atom_list.  If the atom is not found,
        a exception is raised.  The type of exception depends on the argument
        type.  If the argument was a integer, then a IndexError is raised.
        If the argument was a string, then a KeyError is raised.
        """
        if type(name_idx) == StringType:
            return self.atom_dict[name_idx]
        elif type(name_idx) == IntType:
            return self.atom_list[name_idx]
        raise TypeError, name_idx

    def __delitem__(self, name_idx):
        """Removes a Atom from the Fragment.
        """
        pass

    def __iter__(self):
        """Iterates the atoms within the fragment.  If the fragment contains
        atoms in alternate conformations, only the atoms with the structure's
        default_alt_loc are iterated.
        """
        return iter(self.atom_list)

    def __contains__(self, atom_idx):
        """Return True if the Atom object is contained in the fragment.
        """
        if isinstance(atom_idx, Atom):
            return self.atom_list.__contains__(atom_idx)
        elif type(atom_idx) == StringType:
            return self.atom_dict.__contains__(atom_idx)
        raise TypeError, atom_idx

    def index(self, atom):
        """Returns the sequential index of the atom.
        """
        return self.atom_list.index(atom)

    def add_atom(self, atom):
        """Adds a atom to the fragment, and sets the atom's atom.fragment
        attribute to the fragment.
        """
        assert isinstance(atom, Atom)
        assert atom.chain_id == self.chain_id
        assert atom.fragment_id == self.fragment_id
        assert atom.res_name == self.res_name

        name    = atom.name
        alt_loc = atom.alt_loc

        if alt_loc == "":

            try:
                altloc = self.alt_loc_dict[name]

            except KeyError:
                ## case 1:
                ##     add atom without alt_loc partners to the fragment
                ## procedure:
                ##     check if a atom with the same name is already in the
                ##     fragment, and raise a AtomOverwrite exception if
                ##     it is, otherwise, add the atom to the fragment

                try:
                    atomA = self.atom_dict[name]

                except KeyError:
                    ## case 1:
                    ##     add atom to the fragment
                    self.atom_order_list.append(atom)
                    self.atom_list.append(atom)
                    self.atom_dict[name] = atom
                    
                else:
                    ## case 1.5:
                    ##     multiple atoms with the same name, without
                    ##     alt_loc labels, but they are really alt_loc
                    ##     partners
                    assert atomA != atom

                    if atomA.occupancy < 1.0 and atom.occupancy < 1.0:
                        iA = self.atom_order_list.index(atomA)

                        self.alt_loc_dict[name] = altloc = Altloc()
                        self.atom_order_list[iA] = altloc

                        altloc.add_atom(atomA)
                        altloc.add_atom(atom)
                        self.set_alt_loc(self.default_alt_loc)

                    else:
                        raise AtomOverwrite(
                            "overwrite %s with %s" % (atomA, atom))

            else:
                ## case 2:
                ##    adding atom without alt_loc, but partner atoms
                ##    are already in the fragment with alt_loc
                ## procedure:
                ##    set the atom.alt_loc to the next reasonable alt_loc
                ##    and add it to the fragment
                altloc.add_atom(atom)
                self.set_alt_loc(self.default_alt_loc)

        else: ## alt_loc != ""

            try:
                altloc = self.alt_loc_dict[name]

            except KeyError:
                ## case 2:
                ##     add a atom with alt_loc partners to the
                ##     fragment for the first time
                ## procedure:
                ##    *check for atoms without alt_locs already in the
                ##     fragment, and 
                ##    *create new Altloc, and place it in the
                ##     alt_loc_dict under the atom name
                ##    *add the atom to the atom_order_list to preserve
                ##     sequential order of added atoms
                ##    *place atom in the atom_list and atom_dict 

                try:
                    atomA = self.atom_dict[name]

                except KeyError:
                    ## case 2:
                    ##     add a atom with alt_loc partners to the
                    ##     fragment for the first time
                    self.alt_loc_dict[name] = altloc = Altloc()
                    altloc.add_atom(atom)

                    self.atom_order_list.append(altloc)

                    self.atom_list.append(atom)
                    self.atom_dict[name] = atom

                else:
                    ## case 3:
                    ##     add atom with alt_loc, but there is already a
                    ##     atom in the fragment with a null alt_loc which
                    ##     needs to be given a valid alt_loc and placed
                    ##     in the Altloc container before adding the new
                    ##     atom
                    iA = self.atom_order_list.index(atomA)

                    self.alt_loc_dict[name] = altloc = Altloc()
                    self.atom_order_list[iA] = altloc

                    altloc.add_atom(atomA)
                    altloc.add_atom(atom)

            else:
                ## case 4:
                ##     add a atom with alt_loc partners to the
                ##     fragment when there are already alt_loc
                ##     partner atoms in the fragment
                altloc.add_atom(atom)

            self.set_alt_loc(self.default_alt_loc)
            
        atom.fragment = self

    def set_alt_loc(self, alt_loc):
        """Sets the default alt_loc of the Fragment.
        """
        self.default_alt_loc = alt_loc

        ishift = 0
        for i in range(len(self.atom_order_list)):
            atm = self.atom_order_list[i]

            if isinstance(atm, Atom):
                ## case 1: atom has no alt_locs
                try:
                    self.atom_list[i-ishift] = atm
                except IndexError:
                    self.atom_list.append(atm)
                self.atom_dict[atm.name] = atm

            else:
                try:
                    atmx = atm[alt_loc]
                except KeyError:
                    ## case 2: atom has alt_loc partners, but not one
                    ##         for this given alt_loc
                    try:
                        del self.atom_list[i-ishift]
                    except IndexError:
                        pass
                    for atmx in atm.values():
                        try:
                            del self.atom_dict[atmx.name]
                        except KeyError:
                            pass
                        break
                    ishift += 1
                else:
                    ## case 3: atom has alt_loc partners, and one for
                    ##         this alt_loc too
                    try:
                        self.atom_list[i-ishift] = atmx
                    except IndexError:
                        self.atom_list.append(atmx)
                    self.atom_dict[atmx.name] = atmx

    def get_atom(self, name):
        """Returns the matching Atom object contained in the Fragment.
        Returns None if a match is not found.
        """
        try:
            return self[name]
        except KeyError:
            return None
    
    def iter_atoms(self, **args):
        """Iterates over all Atom objects contained in the Fragment matching
        the current model and default alt_loc.
        """
        return iter(self)

    def iter_all_atoms(self):
        """Iterates of all Atoms in the Fragment includeing Altlocs.
        """
        for atm in self.atom_order_list:
            if isinstance(atm, Atom):
                yield atm
            else:
                for atmx in atm:
                    yield atmx

    def iter_bonds(self):
        """Iterates over all Bond objects.  The iteration is preformed by
        iterating over all Atom objects in the same order as iter_atoms(),
        then iterating over each Atom's Bond objects."""
        visited = {}
        for atm in self.iter_atoms():
            for bond in atm.iter_bonds():
                if not visited.has_key(bond):
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
        try:
            self.chain[fragment_id]
        except KeyError:
            pass
        else:
            raise ValueError, fragment_id

        ## set the new fragment_id in all the additional groups

        ## set the new chain_id for the chain object (se
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
        if mon == None:
            return

        for bond in mon.bond_list:
            try:
                atm1 = self[bond["atom1"]]
                atm2 = self[bond["atom2"]]
            except KeyError:
                continue
            else:
                atm1.create_bonds(atom = atm2, standard_res_bond = True)

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
        return self.chain.structure.library.is_water(self.res_name)


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
        return "Res(%s,%s,%s)" % (self.res_name,
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
        if next_res == None:
            return

        library = self.chain.structure.library
        mon1 = library.get_monomer(self.res_name)
        mon2 = library.get_monomer(next_res.res_name)
        if mon1 == None or mon2 == None:
            return

        for (name1, name2) in mon1.get_polymer_bond_list(self, next_res):
            try:
                atm1 = self[name1]
                atm2 = next_res[name2]
            except KeyError:
                continue
            else:
                atm1.create_bonds(atom = atm2, standard_res_bond = True)


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
     
        N_CA  = calc_distance(aN, aCA)
        CA_C  = calc_distance(aCA, aC)
        C_O   = calc_distance(aC, aO)
        C_nN  = calc_distance(aC, naN)
        CA_CB = calc_distance(aCA, aCB)
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

        N_CA_C   = calc_angle(aN, aCA, aC)
        CA_C_O   = calc_angle(aCA, aC, aO)
        N_CA_CB  = calc_angle(aN, aCA, aCB)
        CB_CA_C  = calc_angle(aCB, aCA, aC)
        CA_C_nN  = calc_angle(aCA, aC, naN)
        C_nN_nCA = calc_angle(aC, naN, naCA)

        return (N_CA_C, N_CA_CB, CB_CA_C, CA_C_O, CA_C_nN, C_nN_nCA) 

    def calc_torsion_psi(self):
        """Calculates the Psi torsion angle of the amino acid.  Raises a
        CTerminal exception if called on a C-terminal residue which does
        not have a Psi torsion angle.
        """
        next_res = self.get_offset_residue(1)
        if next_res == None:
            return None

        aN  = self.get_atom('N')
        aCA = self.get_atom('CA')
        aC  = self.get_atom('C')
        naN = next_res.get_atom('N')
        return calc_torsion_angle(aN, aCA, aC, naN)

    def calc_torsion_phi(self):
        """Calculates the Phi torsion angle of the amino acid.  Raises a
        NTerminal exception if called on a N-terminal residue which does
        not have a Phi torsion angle.
        """
        prev_res = self.get_offset_residue(-1)
        if prev_res == None:
            return None

        paC = prev_res.get_atom('C')
        aN  = self.get_atom('N')
        aCA = self.get_atom('CA')
        aC  = self.get_atom('C')
        return calc_torsion_angle(paC, aN, aCA, aC)

    def calc_torsion_omega(self):
        """Calculates the Omega torsion angle of the amino acid. Raises a
        CTerminal exception if called on a C-terminal residue which does
        not have a Omega torsion angle.
        """
        next_res = self.get_offset_residue(1)
        if next_res == None:
            return None

        aCA  = self.get_atom('CA')
        aC   = self.get_atom('C')
        naN  = next_res.get_atom('N')
        naCA = next_res.get_atom('CA')
        return calc_torsion_angle(aCA, aC, naN, naCA)

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

    def calc_torsion(self, torsion_angle_name):
        """Calculates the given torsion angle for the monomer.  The torsion
        angles are defined by name in monomers.cif.
        """
        mon = self.chain.structure.library.get_monomer(self.res_name)
        try:
            (atom1_name,
             atom2_name,
             atom3_name,
             atom4_name) = mon.torsion_angle_dict[torsion_angle_name]
        except KeyError:
            return None

        atom1 = self.get_atom(atom1_name)
        atom2 = self.get_atom(atom2_name)
        atom3 = self.get_atom(atom3_name)
        atom4 = self.get_atom(atom4_name)
        
        return calc_torsion_angle(atom1, atom2, atom3, atom4)

    def calc_torsion_chi1(self):
        return self.calc_torsion("chi1")

    def calc_torsion_chi2(self):
        return self.calc_torsion("chi2")

    def calc_torsion_chi3(self):
        return self.calc_torsion("chi3")

    def calc_torsion_chi4(self):
        return self.calc_torsion("chi4")

    def calc_torsion_chi(self):
        """Calculates CHI side-chain torsion angles according to the
        amino acid specific definitions in the AminoAcids library.
        Returns the 4-tuple (CHI1, CHI2, CHI3, CHI4).  Angles involving
        missing atoms, or angles which do not exist for the amino acid
        are returned as None in the tuple.
        """
        return (self.calc_torsion("chi1"),
                self.calc_torsion("chi2"),
                self.calc_torsion("chi3"),
                self.calc_torsion("chi4"))
        
    def calc_pucker_torsion(self):
        """Calculates the Pucker torsion of a ring system.  Returns None
        for Amino Acids which do not have Pucker torsion angles.
        """
        return self.calc_torsion("pucker")

    
class NucleicAcidResidue(Residue):
    """A subclass of Residue representing one nuclic acid in a strand of
    DNA or RNA.
    """
    pass


class Altloc(dict):
    """
    """
    def __deepcopy__(self, memo):
        altloc = AltLoc()
        for atom in self.values():
            altloc.add_atom(copy.deepcopy(atom, memo))
        return altloc
    
    def __iter__(self):
        """Iterates over all Altloc representations of this Atom.
        """
        alt_locs = self.keys()
        alt_locs.sort()
        for alt_loc in alt_locs:
            yield self[alt_loc]
    
    def add_atom(self, atom):
        """Adds a atom to the Altloc.
        """
        if self.has_key(atom.alt_loc) or atom.alt_loc == "":
            atom.alt_loc = self.calc_next_alt_loc_id(atom)

        self[atom.alt_loc] = atom
        atom.altloc = self
    
    def calc_next_alt_loc_id(self, atom):
        """Returns the next vacant alt_loc letter to be used for a key.
        This is part of a half-ass algorithm to deal with disordered
        input Atoms with improper alt_loc tagging.
        """
        if len(self) == 0:
            return "A"
        for alt_loc in string.uppercase:
            if not self.has_key(alt_loc):
                return alt_loc
        raise AtomOverwrite(
            "exhausted availible alt_loc labels for "+str(atom))
        

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
    Atom.position    - a array[3] (Numeric Python)
    Atom.occupancy   - [1.0 - 0.0] float 
    Atom.temp_factor - float represting B-style temp factor
    Atom.U           - a 6-tuple of the anisotropic values
    Atom.charge      - charge on the atom
    """
    def __init__(
        self,
        name            = "",
        alt_loc         = "",
        res_name        = "",
        fragment_id     = "",
        chain_id        = "",
        model_id        = 1,
        element         = "",
        position        = None,
        x               = None,
        y               = None,
        z               = None,
        sig_position    = None,
        sig_x           = None,
        sig_y           = None,
        sig_z           = None,
        temp_factor     = None,
        sig_temp_factor = None,
        occupancy       = None,
        sig_occupancy   = None,
        charge          = None,

        U   = None,
        u11 = None, u22 = None, u33 = None,
        u12 = None, u13 = None, u23 = None,

        sig_U   = None,
        sig_u11 = None, sig_u22 = None, sig_u33 = None,
        sig_u12 = None, sig_u13 = None, sig_u23 = None,

        **args):

        assert type(name)        == StringType
        assert type(model_id)    == IntType
        assert type(alt_loc)     == StringType
        assert type(res_name)    == StringType
        assert type(fragment_id) == StringType
        assert type(chain_id)    == StringType

        self.name            = name
        self.alt_loc         = alt_loc
        self.res_name        = res_name
        self.fragment_id     = fragment_id
        self.chain_id        = chain_id
        self.model_id        = model_id
        self.element         = element
        self.temp_factor     = temp_factor
        self.sig_temp_factor = sig_temp_factor
        self.occupancy       = occupancy
        self.sig_occupancy   = sig_occupancy
        self.charge          = charge
        
        if type(position) != NoneType:
            self.position = position
        elif x != None or y != None or z != None:
            self.position = array([x, y, z])
        else:
            self.position = None
        
        if type(sig_position) != NoneType:
            self.sig_position = sig_position
        elif sig_x != None or sig_y != None or sig_z != None:
            self.sig_position = array([sig_x, sig_y, sig_z])
        else:
            self.sig_position = None

        if type(U) != NoneType:
            self.U = U
        elif u11 != None:
            self.U = array(
                [ [u11, u12, u13],
                  [u12, u22, u23],
                  [u13, u23, u33] ])
        else:
            self.U = None

        if type(sig_U) != NoneType:
            self.sig_U = sig_U
        elif sig_u11 != None:
            self.sig_U = array(
                [ [sig_u11, sig_u12, sig_u13],
                  [sig_u12, sig_u22, sig_u23],
                  [sig_u13, sig_u23, sig_u33] ])
        else:
            self.sig_U = None

        self.bond_list = []

    def __str__(self):
        return "Atom(%4s%2s%4s%2s%4s%2d)" % (
            self.name, self.alt_loc, self.res_name,
            self.chain_id, self.fragment_id, self.model_id)

    def __deepcopy__(self, memo):
        atom_cpy = Atom(
            name            = self.name,
            alt_loc         = self.alt_loc,
            res_name        = self.res_name,
            fragment_id     = self.fragment_id,
            chain_id        = self.chain_id,
            model_id        = self.model_id,
            element         = self.element,
            position        = self.position,
            sig_position    = self.sig_position,
            temp_factor     = self.temp_factor,
            sig_temp_factor = self.sig_temp_factor,
            occupancy       = self.occupancy,
            sig_occupancy   = self.sig_occupancy,
            charge          = self.charge,
            U               = copy.deepcopy(self.U, memo),
            sig_U           = copy.deepcopy(self.sig_U, memo))
        
        for bond in self.bond_list:
            bond_cpy = copy.deepcopy(bond, memo)
            atom_cpy.bond_list.append(bond_cpy)
            
            if bond_cpy.atom1 == None:
                bond_cpy.atom1 = atom_cpy
            elif bond_cpy.atom2 == None:
                bond_cpy.atom2 = atom_cpy

        return atom_cpy

    def __len__(self):
        """Returns the number of alternate conformations of this atom.
        """
        try:
            return len(self.altloc)
        except AttributeError:
            return 0

    def __getitem__(self, alt_loc):
        """This is a alternative to calling get_alt_loc, but a KeyError
        exception is raised if the alt_loc Atom is not found.  Posiible
        arguments are:
        """
        if type(alt_loc) == StringType:
            try:
                return self.altloc[alt_loc]
            except AttributeError:
                if self.alt_loc == alt_loc:
                    return self
                raise KeyError
        raise TypeError, alt_loc

    def __delitem__(self, x):
        """Deletes the alternate conformation of the atom matching the
        argument.  The argument is the same as the argument for __getitem__.
        """
        self.remove(self[x])

    def __iter__(self):
        """Iterates over all Altloc representations of this Atom.
        """
        try:
            alt_locs = self.altloc.keys()
        except AttributeError:
            yield self
        else:
            alt_locs.sort()
            for alt_loc in alt_locs:
                yield self.altloc[alt_loc]

    def __contains__(self, atom_alt_loc):
        """Returns True if the argument matches a alternate conformation of
        the Atom.  The argument can be a alt_loc label, or a Atom object.
        """
        if isinstance(atom_alt_loc, Atom):
            try:
                return self.altloc[atom_alt_loc.alt_loc] == atom_alt_loc
            except AttributeError:
                return atom_alt_loc == self
        elif type(atom_alt_loc) == StringType:
            try:
                return self.altloc.__contains__(atom_alt_loc)
            except AttributeError:
                return atom_alt_loc == self.alt_loc

    def remove(self, atom):
        """Removes the argument Atom from the Altloc. 
        """
        try:
            self.fragment.remove(atom)
        except AttributeError:
            try:
                del self.altloc[atom.alt_loc]
            except AttributeError:
                pass
            else:
                del self.altloc

    def get_alt_loc(self, alt_loc):
        """Returns the Atom object matching the alt_loc argument.
        """
        try:
            return self[alt_loc]
        except KeyError:
            return None

    def iter_alt_loc(self):
        """Iterate over all alt_loc versions of this atom in the
        alphabetical order of the alt_loc labels, within the current model.
        """
        return iter(self)

    def get_model(self, model):
        """Returns the atom in the argument model number.  Uses the Structure
        default_alt_loc.  Return None if the atom is not found.
        """
        pass
    
    def iter_model(self):
        """Iterates over all models of this atom matching the structure
        wide default_alt_loc.
        """
        pass

    def create_bond(self,
                    atom              = None,
                    bond_type         = None,
                    atom1_symop       = None,
                    atom2_symop       = None,
                    standard_res_bond = False):

        """Creates a bond between this atom and the argumentatom.  The
        arugment bond_type is a string, atom1_symop and atom2_symop are
        symmetry operations to be applied to self and the argument atom
        before distance calculations, and standard_res_bond is a flag
        used to indicate this bond is a standard bond.
        """
        assert isinstance(atom, Atom)
        assert ((self.alt_loc == atom.alt_loc) or
                (self.alt_loc == "" and atom.alt_loc != "") or
                (self.alt_loc != "" and atom.alt_loc == ""))

        bond = Bond(atom1             = self,
                    atom2             = atom,
                    bond_type         = bond_type,
                    atom1_symop       = atom1_symop,
                    atom2_symop       = atom2_symop,
                    standard_res_bond = standard_res_bond)

        self.bond_list.append(bond)
        atom.bond_list.append(bond)

    def create_bonds(self,
                     atom              = None,
                     bond_type         = None,
                     atom1_symop       = None,
                     atom2_symop       = None,
                     standard_res_bond = False):
        """Like create_bonds, but it bonds all alternate locations of this
        atom.
        """
        assert isinstance(atom, Atom)

        try:
            self_altloc = self.altloc
        except AttributeError:
            try:
                atom_altloc = atom.altloc
            except AttributeError:
                ## case 1: self has no alt_loc, atom no alt_loc
                self.create_bond(
                    atom = atom,
                    bond_type = bond_type,
                    atom1_symop = atom1_symop,
                    atom2_symop = atom2_symop,
                    standard_res_bond = standard_res_bond)
            else:
                ## case 2: self.has no alt_loc, atom has alt_loc
                for atmx in atom_altloc.values():
                    self.create_bond(
                        atom = atmx,
                        bond_type = bond_type,
                        atom1_symop = atom1_symop,
                        atom2_symop = atom2_symop,
                        standard_res_bond = standard_res_bond)
        else:
            try:
                atom_altloc = atom.altloc
            except AttributeError:
                ## case 3: self has alt_loc, atom has no alt_loc
                for (alt_loc, atmx) in self_altloc.items():
                    atmx.create_bond(
                        atom = atom,
                        bond_type = bond_type,
                        atom1_symop = atom1_symop,
                        atom2_symop = atom2_symop,
                        standard_res_bond = standard_res_bond)
            else:
                ## case 4: self has alt_loc, atom has alt_loc
                for (alt_loc, atmx) in self_altloc.items():
                    try:
                        atmx.create_bond(
                            atom = atom_altloc[alt_loc],
                            bond_type = bond_type,
                            atom1_symop = atom1_symop,
                            atom2_symop = atom2_symop,
                            standard_res_bond = standard_res_bond)
                    except KeyError:
                        continue

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

    def calc_Uiso(self):
        """Calculates the Uiso tensor from the Atom's temperature factor.
        """
        if self.temp_factor == None:
            return None
        return identity(3) * (self.temp_factor / (24.0 * math.pi * math.pi))

    def get_U(self):
        """Returns the Atoms's U tensor if it exists, otherwise returns
        the isotropic U tensor calculated by self.calc_Uiso
        """
        if self.U != None:
            return self.U
        return self.calc_Uiso()

    def calc_anisotropy(self):
        """Calculates the anisotropy of that atom.  Anisotropy is defined
        as the ratio of the minimum/maximum eigenvalues of the 3x3
        symmetric tensor defined by U.
        """
        ## no Anisotropic values, we have a spherical atom
        if self.U == None:
            return 1.0

        evals = eigenvalues(self.U)
        ansotropy = min(evals) / max(evals)
        return ansotropy

    def calc_anisotropy3(self):
        """Calculates the eigenvalues of the U matrix and returns the
        3-tuple of the eigenvalue ratios: (e1/e2, e1/e3, e2/e3)
        """
        ## no Anisotropic values, we have a spherical atom
        if self.U == None:
            return (1.0, 1.0, 1.0)

        e1, e2, e3 = eigenvalues(self.U)
        elist = [e1, e2, e3]
        elist.sort()
        e1, e2, e3 = elist
        
        return (min(e1, e2) / max(e1, e2),
                min(e1, e3) / max(e1, e3),
                min(e2, e3) / max(e2, e3))

    def calc_CCuij(self, atom_U):
        """Calculates the correlation coefficent between this Atom and the
        argument Atom.  The argument atom may also be a U matrix
        """
        U = self.get_U()
        if isinstance(atom_U, Atom):
            V = atom_U.get_U()
        else:
            V = atom_U
        return calc_CCuij(U, V)

    def calc_Suij(self, atom_U):
        """Compares self Atom with argument Atom or U matrix, and returns
        a value greator than 1.0 when the two Atoms are more alike than a
        isotropic Atom.
        """
        U = self.get_U()
        if isinstance(atom_U, Atom):
            V = atom_U.get_U()
        else:
            V = atom_U
        return calc_Suij(U, V)
        
    def calc_DP2uij(self, atom_U):
        """
        """
        U = self.get_U()
        if isinstance(atom_U, Atom):
            V = atom_U.get_U()
        else:
            V = atom_U
        return calc_DP2uij(U, V)
        
    def iter_atoms_by_distance(self, max_distance = None):
        """Iterates all atoms in the Structure object from the closest to the
        farthest up to the cutoff distance max_distance if given.  Yields
        the 2-tuple (dist, atm).
        """
        listx = []

        if max_distance:
            for atm in self.get_structure().iter_atoms():
                d = calc_distance(self, atm)
                if d <= max_distance:
                    listx.append((calc_distance(self, atm), atm))
        else:
            for atm in self.get_structure().iter_atoms():
                listx.append((calc_distance(self, atm), atm))

        listx.sort()
        return iter(listx)


class Bond(object):
    """Indicates two atoms are bonded together.
    """
    def __init__(
        self,
        atom1             = None,
        atom2             = None,
        bond_type         = None,
        atom1_symop       = None,
        atom2_symop       = None,
        standard_res_bond = False):
        
        self.atom1             = atom1
        self.atom2             = atom2
        self.bond_type         = bond_type
        self.atom1_symop       = atom1_symop
        self.atom2_symop       = atom2_symop
        self.standard_res_bond = standard_res_bond

    def __str__(self):
        return "Bond(%s %s)" % (self.atom1, self.atom2)

    def __deepcopy__(self, memo):
        return Bond(
            bond_type         = self.bond_type,
            atom1_symop       = self.atom1_symop,
            atom2_symop       = self.atom2_symop,
            standard_res_bond = self.standard_res_bond)
    
    def get_partner(self, atm):
        """Returns the other atom involved in the bond.
        """
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


class AlphaHelix(object):
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


class BetaSheet(object):
    """List of Fragments within a structure which are part of a beta
    sheet.
    """
    def __init__(self, structure, sheet_id):
        self.structure = structure
        self.sheet_id = sheet_id

    def __str__(self):
        return "BetaSheet(id=%s)" % (self.sheet_id)


class Site(object):
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
    def __setitem__(self, i, fragment):
        assert isinstance(fragment, Fragment)
        list.__setitem__(self, i, fragment)
    
    def append(self, fragment):
        assert isinstance(fragment, Fragment)
        list.append(self, fragment)

    def insert(self, i, fragment):
        assert isinstance(fragment, Fragment)
        list.insert(self, i, fragment)


class AtomList(list):
    """Provides the functionallity of a Python list class for containing
    Atom instances.  It also provides class methods for performing some
    useful calculations on the list of atoms.
    """
    def __setitem__(self, i, atom):
        assert isinstance(atom, Atom)
        list.__setitem__(self, i, atom)

    def append(self, atom):
        assert isinstance(atom, Atom)
        list.append(self, atom)

    def insert(self, i, atom):
        assert isinstance(atom, Atom)
        list.insert(self, i, atom)
    
    def calc_centroid(self):
        """Calculates the centroid of all contained Atom instances and
        returns a Vector to the centroid.
        """
        num      = 0
        centroid = zeros(3, Float)

        for atm in self:
            if type(atm.position) != NoneType:
                centroid += atm.position
                num += 1

        return centroid / num
        
    def calc_adv_temp_factor(self):
        """Calculates the adverage temperature factor of all contained
        Atom instances and returns the adverage temperature factor.
        """
        num_tf = 0
        adv_tf = 0.0

        for atm in self:
            if atm.temp_factor != None:
                adv_tf += atm.temp_factor
                num_tf += 1

        return adv_tf / num_tf

    def calc_adv_U(self):
        """Calculates the adverage U matrix of all contained Atom
        instances and returns the 3x3 symmetric U matrix of that
        adverage.
        """
        num_U = 0
        adv_U = zeros((3,3), Float)

        for atm in self:
            ## use the atom's U matrix if it exists, otherwise use the
            ## temperature factor

            if atm.U != None:
                adv_U += atm.U
                num_U += 1

        return adv_U / num_U

    def calc_adv_anisotropy(self):
        """Calculates the adverage anisotropy for all Atoms in the AtomList.
        """
        num_atoms = 0
        adv_aniso = 0.0

        for atm in self:
            try:
                adv_aniso += atm.calc_anisotropy()
            except ZeroDivisionError:
                pass
            else:
                num_atoms += 1

        return adv_aniso / num_atoms
        
    def calc_adv_anisotropy3(self):
        """Calculates the adverage anisotropy 3-tuple for all Atoms
        in the AtomList.
        """
        num_atoms  = 0
        adv_aniso1 = 0.0
        adv_aniso2 = 0.0
        adv_aniso3 = 0.0

        for atm in self:
            try:
                a1, a2, a3 = atm.calc_anisotropy3()
            except ZeroDivisionError:
                pass
            else:
                adv_aniso1 += a1
                adv_aniso2 += a2
                adv_aniso3 += a3
                num_atoms  += 1

        return (adv_aniso1 / num_atoms,
                adv_aniso2 / num_atoms,
                adv_aniso3 / num_atoms)

    def calc_adv_temp_factor(self):
        """Calculates the adverage temperature factor of all contained
        Atom instances and returns the adverage temperature factor.
        """
        num_tf = 0
        adv_tf = 0.0

        for atm in self:
            if atm.temp_factor != None:
                adv_tf += atm.temp_factor
                num_tf += 1

        return adv_tf / num_tf

    
### <testing>
if __name__ == "__main__":
    struct = Structure()


    for mx in range(1, 4):
        mid = str(mx)
        for cid in ["A", "B", "C", "D"]:
            for fx in range(1, 4): 
                fid = str(fx)

                for name in ["N", "CA", "C", "O"]:

                    for alt_loc in ["A", "B", "C"]:

                        if alt_loc == "C" and name == "CA": continue

                        atm = Atom(
                            name = name,
                            alt_loc = alt_loc,
                            model = mid,
                            chain_id = cid,
                            res_name = "GLY",
                            fragment_id = fid)

                        struct.add_atom(atm)

    for cx in struct.iter_chains():
        print "iter_chains: ",cx
    for fx in struct.iter_fragments():
        print "iter_fragments: ",fx
        for ax in fx.iter_all_atoms():
            print "iter_all_atoms: ",ax

    for ax in struct.iter_atoms():
        print "iter_atoms: ",ax

    struct.set_alt_loc("C")
    for ax in struct.iter_atoms():
        print "iter_atoms: ",ax


    ## make some exceptions happen
    atm = Atom(
        name = "CA",
        alt_loc = "B",
        model = "1",
        chain_id = "A",
        res_name = "GLY",
        fragment_id = "1")
    
    struct.add_atom(atm)


### </testing>
