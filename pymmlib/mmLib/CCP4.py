## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys

import mmCIF
from   Structure import *


CCP4Error = "CCP4 Error"


try:
    CCP4_ROOT = os.environ["CCP4"]
except KeyError:
    print "Set your CCP4 environment variable."
    sys.exit(1)



class CCP4MonomerFile:
    """The CCP4 monomer description files are wrappted with some HTML
    that needs to be ignored.  Here we construct a minimal file-like
    object to strip out the HTML so it can be passed to the mmCIF parser
    and parsed like a normal mmCIF file."""

    def __init__(self, path, mode):
        self.fil = open(path, mode)


    def readline(self):
        while 1:
            ln = self.fil.readline()

            try:
                if ln[0] != "<": return ln
            except IndexError:
                return ln

            ln2 = ln.rstrip()

            try:
                if ln2[-1] != ">": return ln
            except IndexError:
                return ln



class CCP4MonomerLibrary:

    def __init__(self, ccp4_root = CCP4_ROOT):
        self.ccp4_root = ccp4_root
        self.cif_file_cache = {}

    
    def loadMonomerDescription(self, mon_name):
        """Locates and reads the mmCIF monomer description found in the CCP4
        monomer dictionary."""
        path = os.path.join(self.ccp4_root, "lib", "data", "monomers",
                            mon_name[0].lower(), "%s.cif" % (mon_name))

        ## return None if a description for this monomer name is not found
        if not os.path.isfile(path):
            raise CCP4Error, "monomer path=%s not found" % (path)

        cif_file = mmCIF.mmCIFFile()
        cif_file.loadFile(CCP4MonomerFile(path, "r"))
        return cif_file


    def getMonomerDescription(self, mon_name):
        try:
            return self.cif_file_cache[mon_name]
        except KeyError:
            pass
        cf = self.loadMonomerDescription(mon_name)
        self.cif_file_cache[mon_name] = cf
        return cf


    def buildMonomer(self, mon_name):
        cif_file = self.getMonomerDescription(mon_name)
        cif_data = cif_file.getData("comp_%s" % (mon_name))

        ## create the Amino Acid
        aa_res = AminoAcidResidue()
        aa_res.setName(mon_name)

        ## add the atoms
        for chem_comp_atom in cif_data.chem_comp_atom.getRowList():
            atm = Atom()
            aa_res.appendChild(atm)
            atm.setElement(chem_comp_atom.type_symbol)
            atm.setAtomLabel(chem_comp_atom.atom_id)
            atm.setCharge(chem_comp_atom.partial_charge)

        for chem_comp_bond in cif_data.chem_comp_bond.getRowList():
            atm1 = aa_res.getAtom(chem_comp_bond.atom_id_1)
            atm2 = aa_res.getAtom(chem_comp_bond.atom_id_2)
            bond = atm1.createBond(atm2)
            bond.setDistance(chem_comp_bond.value_dist,
                             chem_comp_bond.value_dist_esd)

        return aa_res


    def getBondAngle(self, mon_name, label1, label2, label3):
        """Returns the bond angle and estimated standard deviation from the
        CCP4 monomer library."""
        cif_file = self.getMonomerDescription(mon_name)
        cif_data = cif_file.getData("comp_%s" % (mon_name))

        
    def getBondDistance(self, mon_name, label1, label2):
        """Returns the bond distance and estimated standard deviation from the
        CCP4 monomer library."""
        cif_file = self.getMonomerDescription(mon_name)
        cif_data = cif_file.getData("comp_%s" % (mon_name))
        
        for chem_comp_bond in cif_data.chem_comp_bond.getRowList():
            if chem_comp_bond.atom_id_1 != label1 and \
               chem_comp_bond.atom_id_1 != label2:
                continue

            if chem_comp_bond.atom_id_2 != label1 and \
               chem_comp_bond.atom_id_2 != label2:
                continue

            return (chem_comp_bond.value_dist, chem_comp_bond.value_dist_esd)

        return None


    def setResidueInfo(self, res):
        """Sets all bond distance information possible for the residue
        from the CCP4 monomer library."""
        mon_name = res.getName()
        
        cif_file = self.getMonomerDescription(mon_name)
        cif_data = cif_file.getData("comp_%s" % (mon_name))
        
        for chem_comp_bond in cif_data.chem_comp_bond.getRowList():
            try:
                atm1 = res.getAtomByLabel(chem_comp_bond.atom_id_1)
                atm2 = res.getAtomByLabel(chem_comp_bond.atom_id_2)
            except MissingAtom:
                continue
            
            bond = atm1.getBond(atm2)
            if bond == None:
                bond = atm1.createBond(atm2)
            
            bond.setDistance(chem_comp_bond.value_dist,
                             chem_comp_bond.value_dist_esd)



##
## <testing>
##

if __name__ == "__main__":
    from AminoAcids import *
    
    mlib = CCP4MonomerLibrary(CCP4_ROOT)
    for name in AminoAcidNames:
        try:
            print name, mlib.bondDistance(name, "CA", "CB")
        except CCP4Error:
            pass

##
## </testing>
##
