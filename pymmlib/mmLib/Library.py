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
from mmCIF import mmCIFFile

    
###############################################################################
## mmLib Library Interface
##

class ElementInterface(object):
    """Public inteface for element properties.
    """
    def __init__(self, **args):
        self.name                 = args.get("name", "")
        self.symbol               = args.get("symbol", "")
        self.group                = args.get("group", "")
        self.period               = args.get("period", "")
        self.atomic_number        = args.get("atomic_number", 0)
        self.atomic_weight        = args.get("atomic_weight", 0.0)
        self.atomic_radius        = args.get("atomic_radius", 0.0)
        self.covalent_radius      = args.get("covalent_radius", 0.0)
        self.van_der_waals_radius = args.get("van_der_waals_radius", 0.0)
        self.covalent_radius      = args.get("covalent_radius", 0.0)
        self.electronegativity    = args.get("electronegativity", 0.0)
        self.color                = args.get("color", (1.0, 1.0, 1.0))

    def __str__(self):
        return "Element=%s" % (self.name)


class MonomerInterface(object):
    """Public interface for monomer properties.
    """
    def __init__(self, **args):
        self.res_name          = args.get("res_name", "")
        self.full_name         = args.get("full_name", "")
        self.one_letter_code   = args.get("one_letter_code", "")
        self.type              = args.get("type", "")
        self.pdbx_type         = args.get("pdbx_type", "")
        self.formula           = args.get("formula", "")
        self.rcsb_class_1      = args.get("rcsb_class_1", "")
        self.chem_type         = args.get("chem_type", "")
        
        self.atom_list         = args.get("atom_list", [])
        self.atom_dict         = args.get("atom_dict", {})
        self.bond_list         = args.get("bond_list", [])
        self.torsion_angle_dict= args.get("torsion_angle_dict", {})
        
    def is_amino_acid(self):
        """Returns True if the Monomer is a amino acid, otherwise
        returns False.
        """
        pass
            
    def is_nucleic_acid(self):
        """Returns True if the Monomer is a nucleic acid, otherwise
        returns False.
        """
        pass

    def is_standard_residue(self):
        """
        """
        return self.is_amino_acid() or self.is_nucleic_acid()
    
    def is_water(self):
        """Returns True if the Monomer is a water molecule,
        otherwise returns False.
        """
        pass

    def get_polymer_bond_list(self, mon1, mon2):
        """Returns a list of 2-tuples.  Each 2-tuple (mon1_name, mon2_name)
        represents one bond between the atom named mon1_name in mon1 and
        the atom named mon2_name in mon2.
        """
        pass


class LibraryInterface(object):
    """The public interface for implementing new Library classes
    """
    def get_element(self, element):
        """Return the corresponding Element description instance for
        the given element symbol.  Returns None if no description is
        found.
        """
        pass
    
    def get_monomer(self, res_name):
        """Returns the Monomer class instance for the monomer defined by
        the 3 letter code res_name.  Returns None if no description is found.
        """
        pass
        
    def is_amino_acid(self, res_name):
        """Returns True if the monomer defined by the 3 letter code
        res_name is a L-amino acid, otherwise returns False.
        """
        pass

    def is_nucleic_acid(self, res_name):
        """Returns True if the monomer defined by res_name is a
        nucleic acid, otherwise returns False.
        """
        pass

    def is_standard_residue(self, res_name):
        """
        """
        return self.is_amino_acid(res_name) or self.is_nucleic_acid(res_name)

    def is_water(self, res_name):
        """Returns True if the 3 letter code res_name is one of the
        synonyms for water, otherwise returns False.
        """
        pass


###############################################################################
## mmLib Library Implementation
##

class Element(ElementInterface):
    """Element description class.  A element instance is loaded from data
    in the mmLib/Data/elements.xml file.
    """
    def __init__(self, **args):
        ElementInterface.__init__(self, **args)

        self.cif_data = args["cif_data"]
        element = self.cif_data["element"]
        self.name = element["name"]
        self.number = int(element["number"])
        self.atomic_weight = float(element["atomic_weight"])
        self.van_der_waals_radius = float(element["van_der_walls_radius"])
        self.covalent_radius = float(element.get("covalent_radius", 0.0))

        rgb = element["color_rgb"]
        self.color = (int(rgb[1:3], 16) / 255.0,
                      int(rgb[3:5], 16) / 255.0,
                      int(rgb[5:7], 16) / 255.0)


class Monomer(MonomerInterface):
    """Monomer definition from RCSB component dictionary.
    """
    def __init__(self, **args):
        MonomerInterface.__init__(self, **args)

        cif_file = mmCIFFile()
        cif_file.load_file(open(args["path"], "r"))

        cif_data = cif_file[0]

        chem_comp = cif_data["chem_comp"][0]

        self.res_name     = chem_comp.get("res_name", self.res_name)
        self.full_name    = chem_comp.get("name", self.full_name)
        self.type         = chem_comp.get("type", self.type)
        self.pdbx_type    = chem_comp.get("pdbx_type", self.pdbx_type)
        self.formula      = chem_comp.get("formula", self.formula)
        self.rcsb_class_1 = chem_comp.get("rcsb_class_1", self.rcsb_class_1)

        try:
            chem_comp_atom = cif_data["chem_comp_atom"]
        except KeyError:
            pass
        else:
            for cif_row in chem_comp_atom:
                name   = cif_row["atom_id"]
                symbol = cif_row["type_symbol"]
            
                self.atom_list.append({"name": name, "symbol": symbol})
                self.atom_dict[name] = symbol

        try:
            chem_comp_bond = cif_data["chem_comp_bond"]
        except KeyError:
            pass
        else:
            for cif_row in chem_comp_bond:
                atom1 = cif_row["atom_id_1"]
                atom2 = cif_row["atom_id_2"]
                self.bond_list.append({"atom1": atom1, "atom2": atom2}) 
        
        ## set items contained in mmlib's mon1_lib, which is a
        ## mmcif file which supplements the RCSB's component
        ## dictionary
        mon1_lib = args["mon1_lib"]

        try:
            self.mon_cif_data = mon1_lib[self.res_name]
        except KeyError:
            pass
        else:

            ## get additional chemical information on amino acids
            try:
                chem_comp2 = self.mon_cif_data["chem_comp"]
            except KeyError:
                pass
            else:
                self.one_letter_code = chem_comp2["one_letter_code"]
                self.chem_type = chem_comp2["chem_type"]

            ## get torsion angle definitions
            try:
                torsion_angles = self.mon_cif_data["torsion_angles"]
            except KeyError:
                pass
            else:
                for cif_row in torsion_angles:
                    self.torsion_angle_dict[cif_row["name"]] = (
                        cif_row["atom1"],
                        cif_row["atom2"],
                        cif_row["atom3"],
                        cif_row["atom4"])              


    def is_amino_acid(self):
        """Returns True if the Monomer is a amino acid, otherwise
        returns False.
        """
        return self.type == "L-PEPTIDE LINKING"
            
    def is_nucleic_acid(self):
        """Returns True if the Monomer is a nucleic acid, otherwise
        returns False.
        """
        return self.type == "DNA LINKING" or self.type == "RNA LINKING"

    def is_water(self):
        """Returns True if the Monomer is a water molecule,
        otherwise returns False.
        """
        return self.name == "HOH"

    def get_polymer_bond_list(self, mon1, mon2):
        """Returns a list of 2-tuples.  Each 2-tuple (mon1_name, mon2_name)
        represents one bond between the atom named mon1_name in mon1 and
        the atom named mon2_name in mon2.
        """
        if self.is_amino_acid()==True:
            return [("C", "N")]
        elif self.is_nucleic_acid()==True:
            return [("O3*", "P")]
        return []


class Library(LibraryInterface):
    """Interface to the mmLib element and monomer library.  The library
    data is in XML and mmCIF files in the Data directory of the mmLib
    distribution.
    """
    def __init__(self):
        LibraryInterface.__init__(self)
        
        self.element_dict = {}
        self.monomer_dict = {}

        ## set library paths -- this involves Python trickery
        (path, x) = os.path.split(__file__)
        
        ## set elements.xml path
        self.elm_lib_path  = os.path.join(path, "Data", "elements.cif")
        self.mon1_lib_path = os.path.join(path, "Data", "monomers.cif")
        self.mon_lib_path  = os.path.join(path, "Data", "Monomers") 

        ## doms to be loaded as needed
        self.elm_cif_file = mmCIFFile()
        self.elm_cif_file.load_file(self.elm_lib_path)
        
        ## open mmlib's monomer suppliment library
        self.mon1_lib = mmCIFFile()
        self.mon1_lib.load_file(self.mon1_lib_path)

    def get_element(self, symbol):
        """Loads and caches a instance of the Element class for the given
        element symbol.  The source of the element data is the elements.xml
        file.
        """
        assert type(symbol)==StringType

        try:
            return self.element_dict[symbol]
        except KeyError:
            pass

        try:
            cif_data = self.elm_cif_file[symbol]
        except KeyError:
            element = None
        else:
            element = Element(cif_data=cif_data, symbol=symbol)
            
        self.element_dict[symbol] = element
        return element

    def get_monomer(self, res_name):
        """Loads and caches the monomer description.
        """
        assert type(res_name) == StringType
        
        try:
            return self.monomer_dict[res_name]
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
        path     = os.path.join(self.mon_lib_path, r0, fil_name)
        
        if os.path.isfile(path):
            mon = Monomer(res_name = res_name,
                          mon1_lib = self.mon1_lib,
                          path     = path)
        else:
            mon = None

        self.monomer_dict[res_name] = mon
        return mon
        
    def is_amino_acid(self, res_name):
        assert type(res_name) == StringType
        try:
            return self.get_monomer(res_name).is_amino_acid()
        except AttributeError:
            return False

    def is_nucleic_acid(self, res_name):
        assert type(res_name) == StringType
        try:
            return self.get_monomer(res_name).is_nucleic_acid()
        except AttributeError:
            return False

    def is_water(self, res_name):
        assert type(res_name) == StringType
        return res_name in ["HOH", "WAT"]


## <TESTING>
if __name__ == "__main__":
    pass
## </TESTING>
