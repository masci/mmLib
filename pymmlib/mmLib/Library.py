## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

"""Monomer and element library data classes.  The Library classes are used
for the identification and construction of biopolymers and ligands."""


class Element(object):
    """Class for holding the properties of a atomic element.
    """
    def __init__(self,
                 name                    = "",
                 symbol                  = "",
                 group                   = "",
                 period                  = "",
                 atomic_number           = 0,
                 atomic_weight           = 0.0,
                 atomic_radius           = 0.0,
                 covalent_radius         = 0.0,
                 van_der_waals_radius    = 2.0,
                 electronegativity       = 0.0):

        self.name                 = name
        self.symbol               = symbol
        self.group                = group
        self.period               = period
        self.atomic_number        = atomic_number
        self.atomic_weight        = atomic_weight
        self.atomic_radius        = atomic_radius
        self.covalent_radius      = covalent_radius
        self.van_der_waals_radius = van_der_waals_radius
        self.electronegativity    = electronegativity

    def __str__(self):
        return "Element=%s" % (self.name)


class Monomer:
    """Base class for all monomer library entries.
    """
    def __init__(self,
                 name            = "",
                 full_name       = "",
                 one_letter_name = "",
                 atom_list       = [],
                 bond_list       = []):

        self.name              = name
        self.full_name         = full_name
        self.one_letter_name   = one_letter_name
        self.atom_list         = atom_list
        self.bond_list         = bond_list

    def is_amino_acid(self):
        """Returns True if the Monomer is a amino acid, otherwise
        returns False.
        """
        return isinstance(self, AminoAcid)

    def is_nucleic_acid(self):
        """Returns True if the Monomer is a nucleic acid, otherwise
        returns False.
        """
        return isinstance(self, NucleicAcid)

    def is_water(self):
        """Returns True if the Monomer is a water molecule,
        otherwise returns False.
        """
        return isinstance(self, Water)

    def get_polymer_bond_list(self, mon1, mon2):
        """Returns a list of 2-tuples.  Each 2-tuple (mon1_name, mon2_name)
        represents one bond between the atom named mon1_name in mon1 and
        the atom named mon2_name in mon2.
        """
        return []


class AminoAcid(Monomer):
    def __init__(self,
                 name                = "",
                 atom_list           = [],
                 bond_list           = [],
                 full_name           = "",
                 one_letter_name     = "",
                 chi1_definition     = None,
                 chi2_definition     = None,
                 chi3_definition     = None,
                 chi4_definition     = None,
                 pucker_definition   = None):

        Monomer.__init__(self,
                         name            = name,
                         full_name       = full_name,
                         one_letter_name = one_letter_name,
                         atom_list       = atom_list,
                         bond_list       = bond_list)

        self.chi1_definition   = chi1_definition
        self.chi2_definition   = chi2_definition
        self.chi3_definition   = chi3_definition
        self.chi4_definition   = chi4_definition
        self.pucker_definition = pucker_definition

    def __str__(self):
	return "AminoAcid(%s)" % (self.name)

    def get_polymer_bond_list(self, mon1, mon2):
        return [("C", "N")]
        

class NucleicAcid(Monomer):
    """Empty definition class for building a Nucleic Acid library.
    """
    def __str__(self):
        return "NucleicAcid=%s" % (self.name)


class Water(Monomer):
    """Monomer class for water molecules.
    """
    pass


class Library:
    """Interface to chemical and monomer libraries.
    """
    def __init__(self):
        from Data.Elements     import ElementMap
        from Data.AminoAcids   import AminoAcidMap
        from Data.NucleicAcids import NucleicAcidMap

        self.element_map      = ElementMap
        self.amino_acid_map   = AminoAcidMap
        self.nucleic_acid_map = NucleicAcidMap

    def __getitem__(self, key):
        try:
            return self.amino_acid_map[key]
        except KeyError:
            pass

        try:
            return self.nucleic_acid_map[key]
        except KeyError:
            pass

        raise KeyError

    def get_element(self, element):
        """Return the corresponding Element description instance for
        the given element symbol.  Returns None if no description is
        found.
        """
        try:
            return self.element_map[element]
        except KeyError:
            pass
        return None

    def get_monomer(self, monomer_id):
        """Returns the corresponding Monomer descripton instance for the
        given monomer_id.  Returns None if no description is found.
        """
        try:
            return self[monomer_id]
        except KeyError:
            pass
        return None
        
    def is_amino_acid(self, res_name):
        return res_name in self.amino_acid_map

    def is_nucleic_acid(self, res_name):
        return res_name in self.nucleic_acid_map

    def is_water(self, res_name):
        return res_name in ["HOH", "WAT"]

