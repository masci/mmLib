## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

## NOTES:
## Standard Van der Walls radii are from J.Phys.Chem., 68, 441, 1964.

"""Monomer and element library data classes.  The Library classes are used
for the identification and construction of biopolymers and ligands.
"""


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
        self.element_map = ElementMap
        self.amino_acid_map = AminoAcidMap
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


### begin elements
    
H  = Element(
    name                 = "Hydrogen",
    symbol               = "H",
    atomic_number        = 1,
    atomic_weight        = 1.007940,
    van_der_waals_radius = 1.20)

He = Element(
    name                 = "Helium",
    symbol               = "He",
    atomic_number        = 2,
    atomic_weight        = 4.002602,
    van_der_waals_radius = 1.40)

Li = Element(
    name                 = "Lithium",
    symbol               = "Li",
    atomic_number        = 3,
    atomic_weight        = 6.941000,
    van_der_waals_radius = 1.82)

Be = Element(
    name                 = "Beryllium",
    symbol               = "Be",
    atomic_number        = 4,
    atomic_weight        = 9.012182)

B  = Element(
    name                 = "Boron",
    symbol               = "B",
    atomic_number        = 5,
    atomic_weight        = 10.811000)

C  = Element(
    name                 = "Carbon",
    symbol               = "C",
    atomic_number        = 6,
    atomic_weight        = 12.010700,
    van_der_waals_radius = 1.70)

N  = Element(
    name                 = "Nitrogen",
    symbol               = "N",
    atomic_number        = 7,
    atomic_weight        = 14.006700,
    van_der_waals_radius = 1.55)

O  = Element(
    name                 = "Oxygen",
    symbol               = "O",
    atomic_number        = 8,
    atomic_weight        = 15.999400,
    van_der_waals_radius = 1.52)

F  = Element(
    name                 = "Fluorine",
    symbol               = "F",
    atomic_number        = 9,
    atomic_weight        = 18.998403,
    van_der_waals_radius = 1.47)

Ne = Element(
    name                 = "Neon",
    symbol               = "Ne",
    atomic_number        = 10,
    atomic_weight        = 20.179700,
    van_der_waals_radius = 1.54)

Na = Element(
    name                 = "Sodium",
    symbol               = "Na",
    atomic_number        = 11,
    atomic_weight        = 22.989770,
    van_der_waals_radius = 2.27)

Mg = Element(
    name                 = "Magnesium",
    symbol               = "Mg",
    atomic_number        = 12,
    atomic_weight        = 24.305000,
    van_der_waals_radius = 1.73)

Al = Element(
    name                 = "Aluminium",
    symbol               = "Al",
    atomic_number        = 13,
    atomic_weight        = 26.981538)

Si = Element(
    name                 = "Silicon",
    symbol               = "Si",
    atomic_number        = 14,
    atomic_weight        = 28.085500,
    van_der_waals_radius = 2.10)

P  = Element(
    name                 = "Phosphorus",
    symbol               = "P",
    atomic_number        = 15,
    atomic_weight        = 30.973761,
    van_der_waals_radius = 1.80)

S  = Element(
    name                 = "Sulfur",
    symbol               = "S",
    atomic_number        = 16,
    atomic_weight        = 32.065000,
    van_der_waals_radius = 1.80)

Cl = Element(
    name                 = "Chlorine",
    symbol               = "Cl",
    atomic_number        = 17,
    atomic_weight        = 35.453000,
    van_der_waals_radius = 1.75)

Ar = Element(
    name                 = "Argon",
    symbol               = "Ar",
    atomic_number        = 18,
    atomic_weight        = 39.948000,
    van_der_waals_radius = 1.88)

K  = Element(
    name                 = "Potassium",
    symbol               = "K",
    atomic_number        = 19,
    atomic_weight        = 39.098300,
    van_der_waals_radius = 2.75)

Ca = Element(
    name                 = "Calcium",
    symbol               = "Ca",
    atomic_number        = 20,
    atomic_weight        = 40.078000)

Sc = Element(
    name                 = "Scandium",
    symbol               = "Sc",
    atomic_number        = 21,
    atomic_weight        = 44.955910)

Ti = Element(
    name                 = "Titanium",
    symbol               = "Ti",
    atomic_number        = 22,
    atomic_weight        = 47.867000)

V  = Element(
    name                 = "Vanadium",
    symbol               = "V",
    atomic_number        = 23,
    atomic_weight        = 50.941500)

Cr = Element(
    name                 = "Chromium",
    symbol               = "Cr",
    atomic_number        = 24,
    atomic_weight        = 51.996100)

Mn = Element(
    name                 = "Manganese",
    symbol               = "Mn",
    atomic_number        = 25,
    atomic_weight        = 54.938049)

Fe = Element(
    name                 = "Iron",
    symbol               = "Fe",
    atomic_number        = 26,
    atomic_weight        = 55.845000)

Co = Element(
    name                 = "Cobalt",
    symbol               = "Co",
    atomic_number        = 27,
    atomic_weight        = 58.933200)

Ni = Element(
    name                 = "Nickel",
    symbol               = "Ni",
    atomic_number        = 28,
    atomic_weight        = 58.693400,
    van_der_waals_radius = 1.63)

Cu = Element(
    name                 = "Copper",
    symbol               = "Cu",
    atomic_number        = 29,
    atomic_weight        = 63.546000,
    van_der_waals_radius = 1.40)

Zn = Element(
    name                 = "Zinc",
    symbol               = "Zn",
    atomic_number        = 30,
    atomic_weight        = 65.409000,
    van_der_waals_radius = 1.39)

Ga = Element(
    name                 = "Gallium",
    symbol               = "Ga",
    atomic_number        = 31,
    atomic_weight        = 69.723000,
    van_der_waals_radius = 1.87)

Ge = Element(
    name                 = "Germanium",
    symbol               = "Ge",
    atomic_number        = 32,
    atomic_weight        = 72.640000)

Ge = Element(
    name                 = "Germanium",
    symbol               = "Ge",
    atomic_number        = 32,
    atomic_weight        = 72.640000)

As = Element(
    name                 = "Arsenic",
    symbol               = "As",
    atomic_number        = 33,
    atomic_weight        = 74.921600,
    van_der_waals_radius = 1.85)

Se = Element(
    name                 = "Selenium",
    symbol               = "Se",
    atomic_number        = 34,
    atomic_weight        = 78.960000,
    van_der_waals_radius = 1.90)

Br = Element(
    name                 = "Bromine",
    symbol               = "Br",
    atomic_number        = 35,
    atomic_weight        = 79.904000,
    van_der_waals_radius = 1.85)

Kr = Element(
    name                 = "Krypton",
    symbol               = "Kr",
    atomic_number        = 36,
    atomic_weight        = 83.798000,
    van_der_waals_radius = 2.02)

Rb = Element(
    name                 = "Rubidium",
    symbol               = "Rb",
    atomic_number        = 37,
    atomic_weight        = 85.467800)

Rb = Element(
    name                 = "Rubidium",
    symbol               = "Rb",
    atomic_number        = 37,
    atomic_weight        = 85.467800)

Sr = Element(
    name                 = "Strontium",
    symbol               = "Sr",
    atomic_number        = 38,
    atomic_weight        = 87.620000)

Sr = Element(
    name                 = "Strontium",
    symbol               = "Sr",
    atomic_number        = 38,
    atomic_weight        = 87.620000)

Y  = Element(
    name                 = "Yttrium",
    symbol               = "Y",
    atomic_number        = 39,
    atomic_weight        = 88.905850)

Zr = Element(
    name                 = "Zirconium",
    symbol               = "Zr",
    atomic_number        = 40,
    atomic_weight        = 91.224000)

Zr = Element(
    name                 = "Zirconium",
    symbol               = "Zr",
    atomic_number        = 40,
    atomic_weight        = 91.224000)

Nb = Element(
    name                 = "Niobium",
    symbol               = "Nb",
    atomic_number        = 41,
    atomic_weight        = 92.906380)

Nb = Element(
    name                 = "Niobium",
    symbol               = "Nb",
    atomic_number        = 41,
    atomic_weight        = 92.906380)

Mo = Element(
    name                 = "Molybdenum",
    symbol               = "Mo",
    atomic_number        = 42,
    atomic_weight        = 95.940000)

Mo = Element(
    name                 = "Molybdenum",
    symbol               = "Mo",
    atomic_number        = 42,
    atomic_weight        = 95.940000)

Tc = Element(
    name                 = "Technetium",
    symbol               = "Tc",
    atomic_number        = 43,
    atomic_weight        = 98.000000)

Tc = Element(
    name                 = "Technetium",
    symbol               = "Tc",
    atomic_number        = 43,
    atomic_weight        = 98.000000)

Ru = Element(
    name                 = "Ruthenium",
    symbol               = "Ru",
    atomic_number        = 44,
    atomic_weight        = 101.070000)

Ru = Element(
    name                 = "Ruthenium",
    symbol               = "Ru",
    atomic_number        = 44,
    atomic_weight        = 101.070000)

Rh = Element(
    name                 = "Rhodium",
    symbol               = "Rh",
    atomic_number        = 45,
    atomic_weight        = 102.905500)

Pd = Element(
    name                 = "Palladium",
    symbol               = "Pd",
    atomic_number        = 46,
    atomic_weight        = 106.420000,
    van_der_waals_radius = 1.63)

Ag = Element(
    name                 = "Silver",
    symbol               = "Ag",
    atomic_number        = 47,
    atomic_weight        = 107.868200,
    van_der_waals_radius = 1.72)

Cd = Element(
    name                 = "Cadmium",
    symbol               = "Cd",
    atomic_number        = 48,
    atomic_weight        = 112.411000,
    van_der_waals_radius = 1.58)

In = Element(
    name                 = "Indium",
    symbol               = "In",
    atomic_number        = 49,
    atomic_weight        = 114.818000)

In = Element(
    name                 = "Indium",
    symbol               = "In",
    atomic_number        = 49,
    atomic_weight        = 114.818000,
    van_der_waals_radius = 1.93)

Sn = Element(
    name                 = "Tin",
    symbol               = "Sn",
    atomic_number        = 50,
    atomic_weight        = 118.710000,
    van_der_waals_radius = 2.17)

Sb = Element(
    name                 = "Antimony",
    symbol               = "Sb",
    atomic_number        = 51,
    atomic_weight        = 121.760000)

Te = Element(
    name                 = "Tellurium",
    symbol               = "Te",
    atomic_number        = 52,
    atomic_weight        = 127.600000,
    van_der_waals_radius = 2.06)

I  = Element(
    name                 = "Iodine",
    symbol               = "I",
    atomic_number        = 53,
    atomic_weight        = 126.904470,
    van_der_waals_radius = 1.98)

Xe = Element(
    name                 = "Xenon",
    symbol               = "Xe",
    atomic_number        = 54,
    atomic_weight        = 131.293000,
    van_der_waals_radius = 2.16)

Cs = Element(
    name                 = "Caesium",
    symbol               = "Cs",
    atomic_number        = 55,
    atomic_weight        = 132.905450)

Ba = Element(
    name                 = "Barium",
    symbol               = "Ba",
    atomic_number        = 56,
    atomic_weight        = 137.327000)

La = Element(
    name                 = "Lanthanum",
    symbol               = "La",
    atomic_number        = 57,
    atomic_weight        = 138.905500)

Ce = Element(
    name                 = "Cerium",
    symbol               = "Ce",
    atomic_number        = 58,
    atomic_weight        = 140.116000)

Pr = Element(
    name                 = "Praseodymium",
    symbol               = "Pr",
    atomic_number        = 59,
    atomic_weight        = 140.907650)

Nd = Element(
    name                 = "Neodymium",
    symbol               = "Nd",
    atomic_number        = 60,
    atomic_weight        = 144.240000)

Pm = Element(
    name                 = "Promethium",
    symbol               = "Pm",
    atomic_number        = 61,
    atomic_weight        = 145.000000)

Sm = Element(
    name                 = "Samarium",
    symbol               = "Sm",
    atomic_number        = 62,
    atomic_weight        = 150.360000)

Eu = Element(
    name                 = "Europium",
    symbol               = "Eu",
    atomic_number        = 63,
    atomic_weight        = 151.964000)

Gd = Element(
    name                 = "Gadolinium",
    symbol               = "Gd",
    atomic_number        = 64,
    atomic_weight        = 157.250000)

Tb = Element(
    name                 = "Terbium",
    symbol               = "Tb",
    atomic_number        = 65,
    atomic_weight        = 158.925340)

Dy = Element(
    name                 = "Dysprosium",
    symbol               = "Dy",
    atomic_number        = 66,
    atomic_weight        = 162.500000)

Ho = Element(
    name                 = "Holmium",
    symbol               = "Ho",
    atomic_number        = 67,
    atomic_weight        = 164.930320)

Er = Element(
    name                 = "Erbium",
    symbol               = "Er",
    atomic_number        = 68,
    atomic_weight        = 167.259000)

Tm = Element(
    name                 = "Thulium",
    symbol               = "Tm",
    atomic_number        = 69,
    atomic_weight        = 168.934210)

Yb = Element(
    name                 = "Ytterbium",
    symbol               = "Yb",
    atomic_number        = 70,
    atomic_weight        = 173.040000)

Lu = Element(
    name                 = "Lutetium",
    symbol               = "Lu",
    atomic_number        = 71,
    atomic_weight        = 174.967000)

Hf = Element(
    name                 = "Hafnium",
    symbol               = "Hf",
    atomic_number        = 72,
    atomic_weight        = 178.490000)

Ta = Element(
    name                 = "Tantalum",
    symbol               = "Ta",
    atomic_number        = 73,
    atomic_weight        = 180.947900)

W  = Element(
    name                 = "Tungsten",
    symbol               = "W",
    atomic_number        = 74,
    atomic_weight        = 183.840000)

Re = Element(
    name                 = "Rhenium",
    symbol               = "Re",
    atomic_number        = 75,
    atomic_weight        = 186.207000)

Os = Element(
    name                 = "Osmium",
    symbol               = "Os",
    atomic_number        = 76,
    atomic_weight        = 190.230000)

Ir = Element(
    name                 = "Iridium",
    symbol               = "Ir",
    atomic_number        = 77,
    atomic_weight        = 192.217000)

Pt = Element(
    name                 = "Platinum",
    symbol               = "Pt",
    atomic_number        = 78,
    atomic_weight        = 195.078000,
    van_der_waals_radius = 1.72)

Au = Element(
    name                 = "Gold",
    symbol               = "Au",
    atomic_number        = 79,
    atomic_weight        = 196.966550,
    van_der_waals_radius = 1.66)

Hg = Element(
    name                 = "Mercury",
    symbol               = "Hg",
    atomic_number        = 80,
    atomic_weight        = 200.590000,
    van_der_waals_radius = 1.55)

Tl = Element(
    name                 = "Thallium",
    symbol               = "Tl",
    atomic_number        = 81,
    atomic_weight        = 204.383300,
    van_der_waals_radius = 1.96)

Pb = Element(
    name                 = "Lead",
    symbol               = "Pb",
    atomic_number        = 82,
    atomic_weight        = 207.200000,
    van_der_waals_radius = 2.02)

Bi = Element(
    name                 = "Bismuth",
    symbol               = "Bi",
    atomic_number        = 83,
    atomic_weight        = 208.980380)

Po = Element(
    name                 = "Polonium",
    symbol               = "Po",
    atomic_number        = 84,
    atomic_weight        = 209.000000)

At = Element(
    name                 = "Astatine",
    symbol               = "At",
    atomic_number        = 85,
    atomic_weight        = 210.000000)

Rn = Element(
    name                 = "Radon",
    symbol               = "Rn",
    atomic_number        = 86,
    atomic_weight        = 222.000000)

Fr = Element(
    name                 = "Francium",
    symbol               = "Fr",
    atomic_number        = 87,
    atomic_weight        = 223.000000)

Ra = Element(
    name                 = "Radium",
    symbol               = "Ra",
    atomic_number        = 88,
    atomic_weight        = 226.000000)

Ac = Element(
    name                 = "Actinium",
    symbol               = "Ac",
    atomic_number        = 89,
    atomic_weight        = 227.000000)

Th = Element(
    name                 = "Thorium",
    symbol               = "Th",
    atomic_number        = 90,
    atomic_weight        = 232.038100)

Pa = Element(
    name                 = "Protactinium",
    symbol               = "Pa",
    atomic_number        = 91,
    atomic_weight        = 231.035880)

U  = Element(
    name                 = "Uranium",
    symbol               = "U",
    atomic_number        = 92,
    atomic_weight        = 238.028910,
    van_der_waals_radius = 1.86)


## this map includes upper-case versions of the element strings
ElementMap = {
    "H"   : H,
    "He"  : He,
    "HE"  : He,
    "Li"  : Li,
    "LI"  : Li,
    "Be"  : Be,
    "BE"  : Be,
    "B"   : B,
    "C"   : C,
    "N"   : N,
    "O"   : O,
    "F"   : F,
    "Ne"  : Ne,
    "NE"  : Ne,
    "Na"  : Na,
    "NA"  : Na,
    "Mg"  : Mg,
    "MG"  : Mg,
    "Al"  : Al,
    "AL"  : Al,
    "Si"  : Si,
    "SI"  : Si,
    "P"   : P,
    "S"   : S,
    "Cl"  : Cl,
    "CL"  : Cl,
    "Ar"  : Ar,
    "AR"  : Ar,
    "K"   : K,
    "Ca"  : Ca,
    "CA"  : Ca,
    "Sc"  : Sc,
    "SC"  : Sc,
    "Ti"  : Ti,
    "TI"  : Ti,
    "V"   : V,
    "Cr"  : Cr,
    "CR"  : Cr,
    "Mn"  : Mn,
    "MN"  : Mn,
    "Fe"  : Fe,
    "FE"  : Fe,
    "Co"  : Co,
    "CO"  : Co,
    "Ni"  : Ni,
    "NI"  : Ni,
    "Cu"  : Cu,
    "CU"  : Cu,
    "Zn"  : Zn,
    "ZN"  : Zn,
    "Ga"  : Ga,
    "GA"  : Ga,
    "Ge"  : Ge,
    "GE"  : Ge,
    "As"  : As,
    "AS"  : As,
    "Se"  : Se,
    "SE"  : Se,
    "Br"  : Br,
    "BR"  : Br,
    "Kr"  : Kr,
    "KR"  : Kr,
    "Rb"  : Rb,
    "RB"  : Rb,
    "Sr"  : Sr,
    "SR"  : Sr,
    "Y"   : Y,
    "Zr"  : Zr,
    "ZR"  : Zr,
    "Nb"  : Nb,
    "NB"  : Nb,
    "Mo"  : Mo,
    "MO"  : Mo,
    "Tc"  : Tc,
    "TC"  : Tc,
    "Ru"  : Ru,
    "RU"  : Ru,
    "Rh"  : Rh,
    "RH"  : Rh,
    "Pd"  : Pd,
    "PD"  : Pd,
    "Ag"  : Ag,
    "AG"  : Ag,
    "Cd"  : Cd,
    "CD"  : Cd,
    "In"  : In,
    "IN"  : In,
    "Sn"  : Sn,
    "SN"  : Sn,
    "Sb"  : Sb,
    "SB"  : Sb,
    "Te"  : Te,
    "TE"  : Te,
    "I"   : I,
    "Xe"  : Xe,
    "XE"  : Xe,
    "Cs"  : Cs,
    "CS"  : Cs,
    "Ba"  : Ba,
    "BA"  : Ba,
    "La"  : La,
    "LA"  : La,
    "Ce"  : Ce,
    "CE"  : Ce,
    "Pr"  : Pr,
    "PR"  : Pr,
    "Nd"  : Nd,
    "ND"  : Nd,
    "Pm"  : Pm,
    "PM"  : Pm,
    "Sm"  : Sm,
    "SM"  : Sm,
    "Eu"  : Eu,
    "EU"  : Eu,
    "Gd"  : Gd,
    "GD"  : Gd,
    "Tb"  : Tb,
    "TB"  : Tb,
    "Dy"  : Dy,
    "DY"  : Dy,
    "Ho"  : Ho,
    "HO"  : Ho,
    "Er"  : Er,
    "ER"  : Er,
    "Tm"  : Tm,
    "TM"  : Tm,
    "Yb"  : Yb,
    "YB"  : Yb,
    "Lu"  : Lu,
    "LU"  : Lu,
    "Hf"  : Hf,
    "HF"  : Hf,
    "Ta"  : Ta,
    "TA"  : Ta,
    "W"   : W,
    "Re"  : Re,
    "RE"  : Re,
    "Os"  : Os,
    "OS"  : Os,
    "Ir"  : Ir,
    "IR"  : Ir,
    "Pt"  : Pt,
    "PT"  : Pt,
    "Au"  : Au,
    "AU"  : Au,
    "Hg"  : Hg,
    "HG"  : Hg,
    "Tl"  : Tl,
    "TL"  : Tl,
    "Pb"  : Pb,
    "PB"  : Pb,
    "Bi"  : Bi,
    "BI"  : Bi,
    "Po"  : Po,
    "PO"  : Po,
    "At"  : At,
    "AT"  : At,
    "Rn"  : Rn,
    "RN"  : Rn,
    "Fr"  : Fr,
    "FR"  : Fr,
    "Ra"  : Ra,
    "RA"  : Ra,
    "Ac"  : Ac,
    "AC"  : Ac,
    "Th"  : Th,
    "TH"  : Th,
    "Pa"  : Pa,
    "PA"  : Pa,
    "U"   : U }


## begin amino acids

## <ALANINE>
ALA = AminoAcid(
    name            = "ALA",
    full_name       = "Alanine",
    one_letter_name = "A",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("HB3", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "HB3"),
                       ("CA", "C"),
                       ("C", "O")])
## </ALANINE>


### <VALINE>
VAL = AminoAcid(
    name            = "VAL",
    full_name       = "Valine",
    one_letter_name = "V",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB", "H"),
                       ("CG1", "C"),
                       ("HG11", "H"),
                       ("HG12", "H"),
                       ("HG13", "H"),
                       ("CG2", "C"),
                       ("HG21", "H"),
                       ("HG22", "H"),
                       ("HG23", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB"),
                       ("CB", "CG1"),
                       ("CG1", "HG11"),
                       ("CG1", "HG12"),
                       ("CG1", "HG13"),
                       ("CB", "CG2"),
                       ("CG2", "HG21"),
                       ("CG2", "HG22"),
                       ("CG2", "HG23"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG1"))
### </VALINE>


### <LEUSINE>
LEU = AminoAcid(
    name            = "LEU",
    full_name       = "Leusine",
    one_letter_name = "L",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("HG", "H"),
                       ("CD1", "C"),
                       ("HD11", "H"),
                       ("HD12", "H"),
                       ("HD13", "H"),
                       ("CD2", "C"),
                       ("HD21", "H"),
                       ("HD22", "H"),
                       ("HD23", "H"),
                       ("C", "C"),
                       ("O", "O")],

    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "HG"),
                       ("CG", "CD1"),
                       ("CD1", "HD11"),
                       ("CD1", "HD12"),
                       ("CD1", "HD13"),
                       ("CG", "CD2"),
                       ("CD2", "HD21"),
                       ("CD2", "HD22"),
                       ("CD2", "HD23"),
                       ("CA", "C"),
                       ("C", "O")],

    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "CD1"))
### </LEUSINE>


### <ISOLEUSINE>
ILE = AminoAcid(
    name            = "ILE",
    full_name       = "Isoleusine",
    one_letter_name = "I",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB", "H"),
                       ("CG1", "C"),
                       ("HG11", "H"),
                       ("HG12", "H"),
                       ("CD1", "C"),
                       ("HD11", "H"),
                       ("HD12", "H"),
                       ("HD13", "H"),
                       ("CG2", "C"),
                       ("HG21", "H"),
                       ("HG22", "H"),
                       ("HG23", "H"),
                       ("C", "C"),
                       ("O", "O")],

    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB"),
                       ("CB", "CG1"),
                       ("CG1", "HG11"),
                       ("CG1", "HG12"),
                       ("CG1", "CD1"),
                       ("CD1", "HD11"),
                       ("CD1", "HD12"),
                       ("CD1", "HD13"),
                       ("CB", "CG2"),
                       ("CG2", "HG21"),
                       ("CG2", "HG22"),
                       ("CG2", "HG23"),
                       ("CA", "C"),
                       ("C", "O")],

    chi1_definition = ("N", "CA", "CB", "CG1"),
    chi2_definition = ("CA", "CB", "CG1", "CD1"))
### <ISOLEUSINE>


### <GLYSINE>
GLY = AminoAcid(
    name            = "GLY",
    full_name       = "Glysine",
    one_letter_name = "G",
    
    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA1", "H"),
                       ("HA2", "H"),
                       ("C", "C"),
                       ("O", "O")],

    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA1"),
                       ("CA", "HA2"),
                       ("CA", "C"),
                       ("C", "O")])
### </GLYSINE>


### <PROLINE>
PRO = AminoAcid(
    name            = "PRO",
    full_name       = "Proline",
    one_letter_name = "P",

    atom_list       = [("N", "N"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("HG1", "H"),
                       ("HG2", "H"),
                       ("CD", "C"),
                       ("HD1", "H"),
                       ("HD2", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "HG1"),
                       ("CG", "HG2"),
                       ("CG", "CD"),
                       ("CD", "HD1"),
                       ("CD", "HD2"),
                       ("CD", "N"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition   = ("N", "CA", "CB", "CG"),
    pucker_definition = ("C", "CA", "CB", "CG"))
### </PROLINE>


### <CYSTINE>
CYS = AminoAcid(
    name            = "CYS",
    full_name       = "Cystine",
    one_letter_name = "C",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("SG", "S"),
                       ("C", "C"),
                       ("O", "O")],

    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "SG"),
                       ("CA", "C"),
                       ("C", "O")],

    chi1_definition = ("N", "CA", "CB", "SG"))
### </CYSTINE>


### <METHIONINE>
MET = AminoAcid(
    name            = "MET",
    full_name       = "Methionine",
    one_letter_name = "M",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("HG1", "H"),
                       ("HG2", "H"),
                       ("SD", "S"),
                       ("CE", "C"),
                       ("HE1", "H"),
                       ("HE2", "H"),
                       ("HE3", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "HG1"),
                       ("CG", "HG2"),
                       ("CG", "SD"),
                       ("SD", "CE"),
                       ("CE", "HE1"),
                       ("CE", "HE2"),
                       ("CE", "HE3"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "SD"),
    chi3_definition = ("CB", "CG", "SD", "CE"))
### </METHIONINE>


### <HISTIDINE>
HIS = AminoAcid(
    name            = "HIS",
    full_name       = "Histidine",
    one_letter_name = "H",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("ND1", "N"),
                       ("HD1", "H"),
                       ("CE1", "C"),
                       ("HE1", "H"),
                       ("NE2", "N"),
                       ("HE2", "H"),
                       ("CD2", "C"),
                       ("HD2", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "CD2"),
                       ("CG", "ND1"),
                       ("ND1", "HD1"),
                       ("ND1", "CE1"),
                       ("CE1", "HE1"),
                       ("CE1", "NE2"),
                       ("NE2", "HE2"),
                       ("NE2", "CD2"),
                       ("CD2", "HD2"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "ND1"))
### </HISTIDINE>


### <PHENYLALALINE>
PHE = AminoAcid(
    name            = "PHE",
    full_name       = "Phenylalaline",
    one_letter_name = "F",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("CD1", "C"),
                       ("HD1", "H"),
                       ("CE1", "C"),
                       ("HE1", "H"),
                       ("CZ", "C"),
                       ("HZ", "H"),
                       ("CE2", "C"),
                       ("HE2", "H"),
                       ("CD2", "C"),
                       ("HD2", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "CD2"),
                       ("CG", "CD1"),
                       ("CD1", "HD1"),
                       ("CD1", "CE1"),
                       ("CE1", "HE1"),
                       ("CE1", "CZ"),
                       ("CZ", "HZ"),
                       ("CZ", "CE2"),
                       ("CE2", "HE2"),
                       ("CE2", "CD2"),
                       ("CD2", "HD2"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "CD1"))
### </PHENYLALALINE>


### <TYROSINE>
TYR = AminoAcid(
    name            = "TYR",
    full_name       = "Tyrosine",
    one_letter_name = "Y",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("CD1", "C"),
                       ("HD1", "H"),
                       ("CE1", "C"),
                       ("HE1", "H"),
                       ("CZ", "C"),
                       ("OH", "O"),
                       ("HH", "H"),
                       ("CE2", "C"),
                       ("HE2", "H"),
                       ("CD2", "C"),
                       ("HD2", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "CD2"),
                       ("CG", "CD1"),
                       ("CD1", "HD1"),
                       ("CD1", "CE1"),
                       ("CE1", "HE1"),
                       ("CE1", "CZ"),
                       ("CZ", "OH"),
                       ("OH", "HH"),
                       ("CZ", "CE2"),
                       ("CE2", "HE2"),
                       ("CE2", "CD2"),
                       ("CD2", "HD2"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "CD1"))
### </TYROSINE>


### <TRPTOPHAN>
TRP = AminoAcid(
    name            = "TRP",
    full_name       = "Trptophan",
    one_letter_name = "W",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("CD1", "C"),
                       ("HD1", "H"),
                       ("NE1", "N"),
                       ("HE1", "H"),
                       ("CE2", "C"),
                       ("CD2", "C"),
                       ("CE3", "C"),
                       ("HE3", "H"),
                       ("CZ3", "C"),
                       ("HZ3", "H"),
                       ("CH2", "C"),
                       ("HH2", "H"),
                       ("CZ2", "C"),
                       ("HZ2", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "CD2"),
                       ("CG", "CD1"),
                       ("CD1", "HD1"),
                       ("CD1", "NE1"),
                       ("NE1", "HE1"),
                       ("NE1", "CE2"),
                       ("CE2", "CZ2"),
                       ("CE2", "CD2"),
                       ("CD2", "CE3"),
                       ("CE3", "HE3"),
                       ("CE3", "CZ3"),
                       ("CZ3", "HZ3"),
                       ("CZ3", "CH2"),
                       ("CH2", "HH2"),
                       ("CH2", "CZ2"),
                       ("CZ2", "HZ2"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "CD1"))
### </TRPTOPHAN>


### <ASPARAGINE>
ASN = AminoAcid(
    name            = "ASN",
    full_name       = "Asparagine",
    one_letter_name = "N",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("OD1", "O"),
                       ("ND2", "N"),
                       ("HD21", "H"),
                       ("HD22", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "OD1"),
                       ("CG", "ND2"),
                       ("ND2", "HD21"),
                       ("ND2", "HD22"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "OD1"))
### </ASPARAGINE>


### <GLUTAMINE>
GLN = AminoAcid(
    name            = "GLN",
    full_name       = "Glutamine",
    one_letter_name = "Q",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("HG1", "H"),
                       ("HG2", "H"),
                       ("CD", "C"),
                       ("OE1", "O"),
                       ("NE2", "N"),
                       ("HE21", "H"),
                       ("HE22", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "HG1"),
                       ("CG", "HG2"),
                       ("CG", "CD"),
                       ("CD", "OE1"),
                       ("CD", "NE2"),
                       ("NE2", "HE21"),
                       ("NE2", "HE22"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "CD"),
    chi3_definition = ("CB", "CG", "CD", "OE1"))
### </GLUTAMINE>


### <SERINE>
SER = AminoAcid(
    name            = "SER",
    full_name       = "Serine",
    one_letter_name = "S",


    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("OG", "O"),
                       ("HG", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "OG"),
                       ("OG", "HG"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "OG"))
### <SERINE>


### <THREONINE>
THR = AminoAcid(
    name            = "THR",
    full_name       = "Threonine",
    one_letter_name = "T",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB", "H"),
                       ("OG1", "O"),
                       ("HG1", "H"),
                       ("CG2", "C"),
                       ("HG21", "H"),
                       ("HG22", "H"),
                       ("HG23", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB"),
                       ("CB", "OG1"),
                       ("OG1", "HG1"),
                       ("CB", "CG2"),
                       ("CG2", "HG21"),
                       ("CG2", "HG22"),
                       ("CG2", "HG23"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "OG1"))
### <THREONINE>


### <LYSINE>
LYS = AminoAcid(
    name            = "LYS",
    full_name       = "lysine",
    one_letter_name = "K",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("HG1", "H"),
                       ("HG2", "H"),
                       ("CD", "C"),
                       ("HD1", "H"),
                       ("HD2", "H"),
                       ("CE", "C"),
                       ("HE1", "H"),
                       ("HE2", "H"),
                       ("NZ", "N"),
                       ("HZ1", "H"),
                       ("HZ2", "H"),
                       ("HZ3", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "HG1"),
                       ("CG", "HG2"),
                       ("CG", "CD"),
                       ("CD", "HD1"),
                       ("CD", "HD2"),
                       ("CD", "CE"),
                       ("CE", "HE1"),
                       ("CE", "HE2"),
                       ("CE", "NZ"),
                       ("NZ", "HZ1"),
                       ("NZ", "HZ2"),
                       ("NZ", "HZ3"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "CD"),
    chi3_definition = ("CB", "CG", "CD", "CE"),
    chi4_definition = ("CG", "CD", "CE", "NZ"))
### </LYSINE>


### <ARGININE>
ARG = AminoAcid(
    name            = "ARG",
    full_name       = "Arginine",
    one_letter_name = "R",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("HG1", "H"),
                       ("HG2", "H"),
                       ("CD", "C"),
                       ("HD1", "H"),
                       ("HD2", "H"),
                       ("NE", "N"),
                       ("HE", "H"),
                       ("CZ", "C"),
                       ("NH1", "N"),
                       ("HH11", "H"),
                       ("HH12", "H"),
                       ("NH2", "N"),
                       ("HH21", "H"),
                       ("HH22", "H"),
                       ("C", "C"),
                       ("O", "O")],

    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "HG1"),
                       ("CG", "HG2"),
                       ("CG", "CD"),
                       ("CD", "HD1"),
                       ("CD", "HD2"),
                       ("CD", "NE"),
                       ("NE", "HE"),
                       ("NE", "CZ"),
                       ("CZ", "NH1"),
                       ("NH1", "HH11"),
                       ("NH1", "HH12"),
                       ("CZ", "NH2"),
                       ("NH2", "HH21"),
                       ("NH2", "HH22"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "CD"),
    chi3_definition = ("CB", "CG", "CD", "NE"),
    chi4_definition = ("CG", "CD", "NE", "CZ"))
### </ARGININE>


### <ASPATATE>
ASP = AminoAcid(
    name            = "ASP",
    full_name       = "Aspatate",
    one_letter_name = "D",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("OD1", "O"),
                       ("OD2", "O"),
                       ("C", "C"),
                       ("O", "O")],

    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "OD1"),
                       ("CG", "OD2"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "OD1"))
### </ASPATATE>


### <GLUTAMATE>
GLU = AminoAcid(
    name            = "GLU",
    full_name       = "Glutamate",
    one_letter_name = "E",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("HG1", "H"),
                       ("HG2", "H"),
                       ("CD", "C"),
                       ("OE1", "O"),
                       ("OE2", "O"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "HG1"),
                       ("CG", "HG2"),
                       ("CG", "CD"),
                       ("CD", "OE1"),
                       ("CD", "OE2"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "CD"),
    chi3_definition = ("CA", "CB", "CG", "OE1"))
### </GLUTAMATE>


## for accessing the AminoAcid classes
AminoAcidMap = {
    "ALA" : ALA,
    "VAL" : VAL,
    "LEU" : LEU,
    "ILE" : ILE,
    "GLY" : GLY,
    "PRO" : PRO,
    "CYS" : CYS,
    "MET" : MET,
    "HIS" : HIS,
    "PHE" : PHE,
    "TYR" : TYR,
    "TRP" : TRP,
    "ASN" : ASN,
    "GLN" : GLN,
    "SER" : SER,
    "THR" : THR,
    "LYS" : LYS,
    "ARG" : ARG,
    "ASP" : ASP,
    "GLU" : GLU,
    }


## begin nucleic acids
NucleicAcidMap = {}




### <TESTING>
if __name__ == "__main__":
    for ename in ElementNames:
        e = ElementMap[ename]
        print '%s = Element(' % (e.symbol.ljust(2))
        print '    name                 = "%s",' % (e.name)
        print '    symbol               = "%s",' % (e.symbol)
        print '    atomic_number        = %d,' % (e.atomic_number)
        print '    atomic_weight        = %3.6f)' % (e.atomic_weight)
        print
### </TESTING>
