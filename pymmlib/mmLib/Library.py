## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Monomer and element library data classes.  The Library classes are used
for the identification and construction of biopolymers and ligands.
"""

import os
import sys
import xml.dom.minidom
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
        self.electronegativity    = args.get("electronegativity", 0.0)
        self.color                = args.get("color", (1.0, 1.0, 1.0))

    def __str__(self):
        return "Element=%s" % (self.name)


class MonomerInterface(object):
    """Public interface for monomer properties.
    """
    def __init__(self, **args):
        self.name              = args.get("name", "")
        self.full_name         = args.get("full_name", "")
        self.one_letter_name   = args.get("one_letter_name", "")
        self.rcsb_type         = args.get("rcsb_type", "")
        self.ccp4_group        = args.get("ccp4_group", "")
        self.pdbx_type         = args.get("pdbx_type", "")
        self.formula           = args.get("formula", "")
        self.rcsb_class        = args.get("rcsb_class", "")
        
        self.atom_list         = args.get("atom_list", [])
        self.atom_dict         = args.get("atom_dict", {})
        self.bond_list         = args.get("bond_list", [])
        self.chi1_definition   = args.get("chi1_definition", None)
        self.chi2_definition   = args.get("chi2_definition", None)
        self.chi3_definition   = args.get("chi3_definition", None)
        self.chi4_definition   = args.get("chi4_definition", None)
        self.pucker_definition = args.get("pucker_definition", None)

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

        self.doc = args["xml_doc"]
        self.elm = None

        ## find xml document element for this element and
        ## store it in self.elm
        self.symbol = self.symbol.capitalize()

        for elx in self.doc.documentElement.getElementsByTagName("Element"):
            if elx.getAttribute("symbol") == self.symbol:
                self.elm = elx

        ## unable to find element description
        if self.elm == None:
            #print "Unable to find element: %s" % (self.symbol)
            return

        self.name   = self.elm.getAttribute("name")
        self.number = int(self.elm.getAttribute("number"))
        
        for elx in self.elm.getElementsByTagName("AtomicWeight"):
            self.atomic_weight = elx.getAttribute("value")

        for elx in self.elm.getElementsByTagName("VanderWallRadius"):
            self.van_der_waals_radius = elx.getAttribute("value")

        for elx in self.elm.getElementsByTagName("Color"):
            rgb        = elx.getAttribute("rgb")
            self.color = (int(rgb[1:3], 16) / 255.0,
                          int(rgb[3:5], 16) / 255.0,
                          int(rgb[5:7], 16) / 255.0)


class Monomer(MonomerInterface):
    """Base class for all monomer library entries.  Load monomer
    description from monomer library files located in Data/Monomers
    """
    def __init__(self, **args):
        MonomerInterface.__init__(self, **args)

        ## load the XML library file
        doc       = xml.dom.minidom.parse(open(args["path"], "r"))
        elm       = doc.documentElement
        self.name = str(elm.getAttribute("name"))

        for elx in elm.getElementsByTagName("FullName"):
            self.full_name = str(elx.getAttribute("value"))

        for elx in elm.getElementsByTagName("CCP4Group"):
            self.ccp4_group = str(elx.getAttribute("value"))
            
        for elx in elm.getElementsByTagName("RCSBType"):
            self.rcsb_type = str(elx.getAttribute("value"))

        for elx in elm.getElementsByTagName("PDBXType"):
            self.pdbx_type = str(elx.getAttribute("value"))

        for elx in elm.getElementsByTagName("Formula"):
            self.formula = str(elx.getAttribute("value"))

        for elx in elm.getElementsByTagName("RCSBClass"):
            self.rcsb_class = str(elx.getAttribute("value"))
            
        for elx in elm.getElementsByTagName("Atom"):
            name   = str(elx.getAttribute("atom"))
            symbol = str(elx.getAttribute("symbol"))
            
            self.atom_list.append({"name": name, "symbol": symbol})
            self.atom_dict[name] = symbol

        for elx in elm.getElementsByTagName("Bond"):
            self.bond_list.append(
                {"atom1": str(elx.getAttribute("atom1")),
                 "atom2": str(elx.getAttribute("atom2"))})

        ## flag setting logic
        self.amino_acid_flag = (
            self.rcsb_type.upper().count("L-PEPTIDE") or
            self.ccp4_group.upper().count("L-PEPTIDE LINKING") or
            self.rcsb_class.upper().count("STANDARD ALPHA AMINO ACIDS")) > 0

        self.nucleic_acid_flag = (
            self.rcsb_type.upper().count("DNA") or
            self.rcsb_type.upper().count("RNA") or
            self.ccp4_group.upper().count("DNA") or
            self.ccp4_group.upper().count("RNA") or
            self.rcsb_class.upper().count("NUCLEIC")) > 0

        self.water_flag = self.name == "HOH" or self.name == "WAT"
        
    def is_amino_acid(self):
        """Returns True if the Monomer is a amino acid, otherwise
        returns False.
        """
        return self.amino_acid_flag
            
    def is_nucleic_acid(self):
        """Returns True if the Monomer is a nucleic acid, otherwise
        returns False.
        """
        return self.nucleic_acid_flag

    def is_water(self):
        """Returns True if the Monomer is a water molecule,
        otherwise returns False.
        """
        return self.water_flag
    
    def get_polymer_bond_list(self, mon1, mon2):
        """Returns a list of 2-tuples.  Each 2-tuple (mon1_name, mon2_name)
        represents one bond between the atom named mon1_name in mon1 and
        the atom named mon2_name in mon2.
        """
        if self.amino_acid_flag:
            return [("C", "N")]
        elif self.nucleic_acid_flag:
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
        self.elm_lib_path = os.path.join(path, "Data", "elements.xml")
        self.mon_lib_path = os.path.join(path, "Data", "Monomers.xml") 

        ## doms to be loaded as needed
        self.elm_doc = xml.dom.minidom.parse(open(self.elm_lib_path, "r"))

    def get_element(self, symbol):
        """Loads and caches a instance of the Element class for the given
        element symbol.  The source of the element data is the elements.xml
        file.
        """
        assert type(symbol) == StringType

        try:
            return self.element_dict[symbol]
        except KeyError:
            pass
        
        element = Element(xml_doc = self.elm_doc, symbol = symbol)
        if element.elm == None:
            element = None
        
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

        ## form path to locate the monomer library file
        try:
            r0 = res_name[0].lower()
        except IndexError:
            return None

        ## this hack is necessary beause most PDB files do not have
        ## a indicator in the 3-letter-code for RNA or DNA
        if res_name in ["C","G","A","T","U"]:
            fil_name = "%sr.xml" % (res_name)
        else:
            fil_name = "%s.xml" % (res_name)
            
        path = os.path.join(self.mon_lib_path, r0, fil_name)
        if os.path.isfile(path):
            mon = Monomer(name = res_name, path = path)
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
