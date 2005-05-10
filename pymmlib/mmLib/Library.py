## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Monomer and element library data classes.  The Library classes are used
for the identification and construction of biopolymers and ligands.
"""
import os
import sys

from mmTypes import *
from mmCIF   import mmCIFFile


###############################################################################
## Library Data Locations
##
(MMLIB_PATH, JUNK) = os.path.split(__file__)

ELEMENT_DATA_PATH       = os.path.join(MMLIB_PATH, "Data", "elements.cif")
MMLIB_MONOMER_DATA_PATH = os.path.join(MMLIB_PATH, "Data", "monomers.cif")
RCSB_MONOMER_DATA_PATH  = os.path.join(MMLIB_PATH, "Data", "Monomers") 
    
###############################################################################
## Caches
##
ELEMENT_CIF_FILE       = None
ELEMENT_CACHE          = {}
MONOMER_RES_NAME_CACHE = {}

###############################################################################
## Constants
##
ELEMENT_SYMBOL_DICT = {
    "H" : True, "h" : True,
    "He": True, "he": True, "HE": True,
    "Li": True, "li": True, "LI": True,
    "Be": True, "be": True, "BE": True,
    "B" : True, "b" : True,
    "C" : True, "c" : True,
    "N" : True, "n" : True,
    "O" : True, "o" : True,
    "F" : True, "f" : True,
    "Ne": True, "ne": True, "NE": True,
    "Na": True, "na": True, "NA": True,
    "Mg": True, "mg": True, "MG": True,
    "Al": True, "al": True, "AL": True,
    "Si": True, "si": True, "SI": True,
    "P" : True, "p" : True,
    "S" : True, "s" : True,
    "Cl": True, "cl": True, "CL": True,
    "Ar": True, "ar": True, "AR": True,
    "K" : True, "k" : True,
    "Ca": True, "ca": True, "CA": True,
    "Sc": True, "sc": True, "SC": True,
    "Ti": True, "ti": True, "TI": True,
    "V" : True, "v" : True,
    "Cr": True, "cr": True, "CR": True,
    "Mn": True, "mn": True, "MN": True,
    "Fe": True, "fe": True, "FE": True,
    "Co": True, "co": True, "CO": True,
    "Ni": True, "ni": True, "NI": True,
    "Cu": True, "cu": True, "CU": True,
    "Zn": True, "zn": True, "ZN": True,
    "Ga": True, "ga": True, "GA": True,
    "Ge": True, "ge": True, "GE": True,
    "As": True, "as": True, "AS": True,
    "Se": True, "se": True, "SE": True,
    "Br": True, "br": True, "BR": True,
    "Kr": True, "kr": True, "KR": True,
    "Rb": True, "rb": True, "RB": True,
    "Sr": True, "sr": True, "SR": True,
    "Y" : True, "y" : True,
    "Zr": True, "zr": True, "ZR": True,
    "Nb": True, "nb": True, "NB": True,
    "Mo": True, "mo": True, "MO": True,
    "Tc": True, "tc": True, "TC": True,
    "Ru": True, "ru": True, "RU": True,
    "Rh": True, "rh": True, "RH": True,
    "Pd": True, "pd": True, "PD": True,
    "Ag": True, "ag": True, "AG": True,
    "Cd": True, "cd": True, "CD": True,
    "In": True, "in": True, "IN": True,
    "Sn": True, "sn": True, "SN": True,
    "Sb": True, "sb": True, "SB": True,
    "Te": True, "te": True, "TE": True,
    "I" : True, "i" : True,
    "Xe": True, "xe": True, "XE": True,
    "Cs": True, "cs": True, "CS": True,
    "Ba": True, "ba": True, "BA": True,
    "La": True, "la": True, "LA": True,
    "Ce": True, "ce": True, "CE": True,
    "Pr": True, "pr": True, "PR": True,
    "Nd": True, "nd": True, "ND": True,
    "Pm": True, "pm": True, "PM": True,
    "Sm": True, "sm": True, "SM": True,
    "Eu": True, "eu": True, "EU": True,
    "Gd": True, "gd": True, "GD": True,
    "Tb": True, "tb": True, "TB": True,
    "Dy": True, "dy": True, "DY": True,
    "Ho": True, "ho": True, "HO": True,
    "Er": True, "er": True, "ER": True,
    "Tm": True, "tm": True, "TM": True,
    "Yb": True, "yb": True, "YB": True,
    "Lu": True, "lu": True, "LU": True,
    "Hf": True, "hf": True, "HF": True,
    "Ta": True, "ta": True, "TA": True,
    "W" : True, "w" : True,
    "Re": True, "re": True, "RE": True,
    "Os": True, "os": True, "OS": True,
    "Ir": True, "ir": True, "IR": True,
    "Pt": True, "pt": True, "PT": True,
    "Au": True, "au": True, "AU": True,
    "Hg": True, "hg": True, "HG": True,
    "Tl": True, "tl": True, "TL": True,
    "Pb": True, "pb": True, "PB": True,
    "Bi": True, "bi": True, "BI": True,
    "Po": True, "po": True, "PO": True,
    "At": True, "at": True, "AT": True,
    "Rn": True, "rn": True, "RN": True,
    "Fr": True, "fr": True, "FR": True,
    "Ra": True, "ra": True, "RA": True,
    "Ac": True, "ac": True, "AC": True,
    "Th": True, "th": True, "TH": True,
    "Pa": True, "pa": True, "PA": True,
    "U" : True, "u" : True }

AMINO_ACID3_LIST = [
    "GLY", "ALA", "VAL", "LEU", "ILE", "PRO", "PHE", "TYR", "TRP",
    "MET", "CYS", "SER", "THR", "ASP", "GLU", "HIS", "LYS", "ARG",
    "ASN", "GLN"
    ]

AMINO_ACID31_DICT = {
    "GLY":"G", "ALA":"A", "VAL":"V", "LEU":"L", "ILE":"I", "PRO":"P",
    "PHE":"F", "TYR":"Y", "TRP":"W", "MET":"M", "CYS":"C", "SER":"S",
    "THR":"T", "ASP":"D", "GLU":"E", "HIS":"H", "LYS":"K", "ARG":"R",
    "ASN":"N", "GLN":"Q"
    }

AMINO_ACID13_DICT = {
    'A': 'ALA', 'C': 'CYS', 'E': 'GLU', 'D': 'ASP', 'G': 'GLY',
    'F': 'PHE', 'I': 'ILE', 'H': 'HIS', 'K': 'LYS', 'M': 'MET',
    'L': 'LEU', 'N': 'ASN', 'Q': 'GLN', 'P': 'PRO', 'S': 'SER',
    'R': 'ARG', 'T': 'THR', 'W': 'TRP', 'V': 'VAL', 'Y': 'TYR'}

###############################################################################
## Library Description Objects
##

class ElementDesc(object):
    """Element description class returned by library_get_element_desc().
    """
    def __init__(self):
        self.cif_data             = None
        self.name                 = None
        self.symbol               = None
        self.group                = None
        self.period               = None
        self.atomic_number        = None
        self.atomic_weight        = None
        self.atomic_radius        = None
        self.covalent_radius      = None
        self.van_der_waals_radius = None
        self.covalent_radius      = None
        self.electronegativity    = None
        self.color_rgbf           = None


class MonomerDesc(object):
    """Monomer description class returned by library_get_monomer_desc().
    """
    def __init__(self):
        self.res_name           = None
        self.full_name          = None
        self.one_letter_code    = None
        self.type               = None
        self.pdbx_type          = None
        self.formula            = None
        self.rcsb_class_1       = None
        self.chem_type          = None
        self.atom_list          = []
        self.atom_dict          = {}
        self.bond_list          = []
        self.torsion_angle_dict = {}

        self.amino_acid         = False
        self.nucleic_acid       = False
        self.water              = False
        
    def is_amino_acid(self):
        """Returns True if the Monomer is a amino acid, otherwise
        returns False.
        """
        return self.amino_acid
            
    def is_nucleic_acid(self):
        """Returns True if the Monomer is a nucleic acid, otherwise
        returns False.
        """
        return self.nucleic_acid

    def is_standard_residue(self):
        """
        """
        return self.amino_acid or self.nucleic_acid

    def is_non_standard_residue(self):
        """
        """
        return not self.amino_acid and not self.nucleic_acid

    def is_water(self):
        """Returns True if the Monomer is a water molecule,
        otherwise returns False.
        """
        return self.water


###############################################################################
## Library API
##

def library_construct_element_desc(symbol):
    """Constructs the ElementDesc object for the given element symbol.
    """
    ## only load elements.cif once
    global ELEMENT_CIF_FILE
    if ELEMENT_CIF_FILE==None:
        ELEMENT_CIF_FILE = mmCIFFile()
        ELEMENT_CIF_FILE.load_file(ELEMENT_DATA_PATH)

    cif_data = ELEMENT_CIF_FILE.get_data(symbol)
    if cif_data==None:
        warning("element description not found for %s" % (symbol))
        return None

    ## create element description
    element_desc = ElementDesc()

    element_desc.cif_data = cif_data

    element = cif_data.get_table("element")
    element_desc.name            = element["name"]
    element_desc.symbol          = element["symbol"]
    element_desc.number          = int(element["number"])
    element_desc.atomic_weight   = float(element["atomic_weight"])
    element_desc.vdw_radius      = float(element["van_der_walls_radius"])
    element_desc.covalent_radius = float(element.get("covalent_radius", 0.0))
    
    rgb8 = element["color_rgb"]
    element_desc.color_rgbf = (int(rgb8[1:3], 16) / 255.0,
                               int(rgb8[3:5], 16) / 255.0,
                               int(rgb8[5:7], 16) / 255.0)

    return element_desc


def library_get_element_desc(symbol):
    """Loads/caches/returns a instance of the ElementDesc class for the given
    element symbol.  The source of the element data is the
    mmLib/Data/elements.cif file.
    """
    assert type(symbol)==StringType

    try:
        return ELEMENT_CACHE[symbol]
    except KeyError:
        pass

    element_desc = library_construct_element_desc(symbol)
    if element_desc==None:
        warning("element description not found for %s" % (symbol))
        return None
    
    ELEMENT_CACHE[symbol] = element_desc
    return element_desc


def library_construct_monomer_desc(res_name):
    """Constructs the MonomerDesc object for the given residue name.
    """
    ## this hack is necessary beause most PDB files do not have
    ## a indicator in the 3-letter-code for RNA or DNA
    dna_map = {
        "C+": "C", "Cr": "C", "+C": "C",
        "G+": "G", "Gr": "G", "+G": "G",
        "A+": "A", "Ar": "A", "+A": "A",
        "T+": "T", "Tr": "T", "+T": "T",
        "U+": "U", "Ur": "U", "+U": "U",
        }

    try:
        lookup_name = dna_map[res_name]
    except KeyError:
        lookup_name = res_name.upper()

    ## form path to locate the monomer library file
    try:
        r0 = lookup_name[0]
    except IndexError:
        return None

    fil_name = "%s.cif" % (lookup_name.upper())
    path     = os.path.join(RCSB_MONOMER_DATA_PATH, r0, fil_name)

    if not os.path.isfile(path):
        warning("monomer description not found for %s" % (res_name))
        return None

    ## generate monomer description    
    mon_desc = MonomerDesc()

    ## data from RCSB library
    rcsb_cif_file = mmCIFFile()
    rcsb_cif_file.load_file(open(path, "r"))
    rcsb_cif_data = rcsb_cif_file[0]

    chem_comp = rcsb_cif_data.get_table("chem_comp")[0]

    mon_desc.res_name     = chem_comp.get("res_name")
    mon_desc.full_name    = chem_comp.get("name")
    mon_desc.type         = chem_comp.get("type")
    mon_desc.pdbx_type    = chem_comp.get("pdbx_type")
    mon_desc.formula      = chem_comp.get("formula")
    mon_desc.rcsb_class_1 = chem_comp.get("rcsb_class_1")

    chem_comp_atom = rcsb_cif_data.get_table("chem_comp_atom")
    if chem_comp_atom!=None:
        for cif_row in chem_comp_atom:
            name   = cif_row["atom_id"]
            symbol = cif_row["type_symbol"]
            
            mon_desc.atom_list.append({"name": name, "symbol": symbol})
            mon_desc.atom_dict[name] = symbol

    chem_comp_bond = rcsb_cif_data.get_table("chem_comp_bond")
    if chem_comp_bond!=None:
        for cif_row in chem_comp_bond:
            atom1 = cif_row["atom_id_1"]
            atom2 = cif_row["atom_id_2"]
            mon_desc.bond_list.append({"atom1": atom1, "atom2": atom2}) 

    ## data from mmLib supplimental library in mmLib/Monomers/monomers.cif
    mmlib_cif_file = mmCIFFile()
    mmlib_cif_file.load_file(open(MMLIB_MONOMER_DATA_PATH, "r"))

    mmlib_cif_data = mmlib_cif_file.get_data(res_name)
    if mmlib_cif_data!=None:
        ## get additional chemical information on amino acids
        chem_comp = mmlib_cif_data.get_table("chem_comp")
        if chem_comp!=None:
            mon_desc.one_letter_code = chem_comp["one_letter_code"]
            mon_desc.chem_type = chem_comp["chem_type"]

        ## get torsion angle definitions
        torsion_angles = mmlib_cif_data.get_table("torsion_angles")
        if torsion_angles!=None:
            for cif_row in torsion_angles:
                mon_desc.torsion_angle_dict[cif_row["name"]] = (
                    cif_row["atom1"],
                    cif_row["atom2"],
                    cif_row["atom3"],
                    cif_row["atom4"])              


    ## set some derived flags on the monomer description
    mon_type = mon_desc.type.upper()
    
    if mon_type=="L-PEPTIDE LINKING":
        mon_desc.amino_acid = True

    elif mon_type=="DNA LINKING" or mon_type=="RNA LINKING":
        mon_desc.nucleic_acid = True

    elif mon_type=="HOH" or mon_type=="WAT":
        mon_desc.water = True

    return mon_desc

    
def library_get_monomer_desc(res_name):
    """Loads/caches/returns the monomer description objec MonomerDesc
    for the given monomer residue name.
    """
    assert type(res_name)==StringType

    try:
        return MONOMER_RES_NAME_CACHE[res_name]
    except KeyError:
        pass

    mon_desc = library_construct_monomer_desc(res_name)
    if mon_desc==None:
        return None

    MONOMER_RES_NAME_CACHE[res_name] = mon_desc
    return mon_desc


def library_is_amino_acid(res_name):
    """Returns True if the res_name is a amino acid.
    """
    assert type(res_name)==StringType

    mdesc = library_get_monomer_desc(res_name)
    if mdesc==None:
        return False

    return mdesc.is_amino_acid()


def library_is_nucleic_acid(res_name):
    """Returns True if the res_name is a nucleic acid.
    """
    assert type(res_name)==StringType

    mdesc = library_get_monomer_desc(res_name)
    if mdesc==None:
        return False

    return mdesc.is_nucleic_acid()


def library_is_water(res_name):
    """Return True if the res_name is water.
    """
    assert type(res_name)==StringType
    if res_name=="HOH" or res_name=="WAT":
        return True
    return False


def library_guess_element_from_name(name0, res_name):
    """Try everything I can possibly think of to extract the element
    symbol from the atom name.  If availible, use the monmer dictionary
    to help narrow down the search.
    """
    ## strip any space from the name, and return now if there
    ## is nothing left to work with
    name = name0.strip()
    if name=="":
        return None

    ## try the easy way out -- look up the atom in the monomer dictionary
    mdesc = library_get_monomer_desc(res_name)
    if mdesc!=None:
        if mdesc.atom_dict.has_key(name):
            symbol = mdesc.atom_dict[name]
            if symbol!=None:
                return symbol

        if mdesc.is_amino_acid() and name=="OXT":
            return "O"

        if mdesc.is_amino_acid():
            warning("invalid amino acid atom name %s" % (name))

    ## ok, that didn't work...

    ## set the space_flag to true if the name starts with a space
    ## which can indicate the name of the atom is only 1 charactor
    ## long
    if name0.startswith(" "):
        space_flag = True
    else:
        space_flag = False

    ## remove all non-alpha chars from the name
    alpha_name = ""
    for c in name:
        if c.isalpha()==True:
            alpha_name += c

    ## look up two possible element symbols in the library:
    ## e1 is the possible one-charactor symbol
    ## e2 is the possible two-charactor symbol
    if len(alpha_name)==0:
        return None

    e1_symbol = alpha_name[0]
    e1_valid  = ELEMENT_SYMBOL_DICT.has_key(e1_symbol)

    if len(alpha_name)>1:
        e2_symbol = alpha_name[:2]
        e2_valid  = ELEMENT_SYMBOL_DICT.has_key(e2_symbol)
    else:
        e2_symbol = None
        e2_valid  = False

    ## e1 or e2 must return somthing for us to proceed, otherwise,
    ## there's just no possible element symbol contained in the atom
    ## name
    if e1_valid==False and e2_valid==False:
        return None

    elif e1_valid==True and e2_valid==False:
        return e1_symbol

    elif e1_valid==False and e2_valid==True:
        return e2_symbol

    ## if we get here, then e1 and e2 are both valid elements

    ## we're out of choices, go by the space_flag: if there is a space
    ## before the atom name, then use the 1-char element symbol;
    ## if there is no space, then use the 2-char element symbol
    if space_flag==True:
        return e1_symbol

    return e2_symbol


## <TESTING>
if __name__ == "__main__":
    h = library_get_element_desc("H")

    for cif_data in ELEMENT_CIF_FILE:
        if len(cif_data.name)==1:
            print '    "%s" : True, "%s" : True,' % (
                cif_data.name, cif_data.name.lower())
        else:
            print '    "%s": True, "%s": True, "%s": True,' % (
                cif_data.name, cif_data.name.lower(), cif_data.name.upper())
## </TESTING>
