## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.


class Monomer:
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

    def getPolymerBondList(self, mon1, mon2):
        """Returns a list of 2-tuples.  Each 2-tuple (mon1_name, mon2_name)
        represents one bond between the atom named mon1_name in mon1 and
        the atom named mon2_name in mon2."""
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

    def getPolymerBondList(self, mon1, mon2):
        return [("C", "N")]
        

class NucleicAcid(Monomer):
    pass
