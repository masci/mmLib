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

ELEMENT_CACHE          = {}
MONOMER_RES_NAME_CACHE = {}

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

    cif_file = mmCIFFile()
    cif_file.load_file(ELEMENT_DATA_PATH)
        
    cif_data = cif_file.get_data(symbol)
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
    
    ELEMENT_CACHE[symbol] = element_desc
    return element_desc


def library_get_monomer_desc(res_name):
    """Loads/caches/returns the monomer description objec MonomerDesc
    for the given monomer residue name.
    """
    assert type(res_name)==StringType

    try:
        return MONOMER_RES_NAME_CACHE[res_name]
    except KeyError:
        pass

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


## <TESTING>
if __name__ == "__main__":
    pass
## </TESTING>
