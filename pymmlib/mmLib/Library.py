## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

from Elements     import ElementMap
from AminoAcids   import AminoAcidMap
from NucleicAcids import NucleicAcidMap


class Library:
    def __init__(self):
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
        
    def isAminoAcid(self, res_name):
        return res_name in self.amino_acid_map

    def isNucleicAcid(self, res_name):
        return res_name in self.nucleic_acid_map

    def isWater(self, res_name):
        return res_name in ["HOH", "WAT"]

