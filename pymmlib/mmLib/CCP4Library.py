## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

"""Subclass of mmLib.Library providing a supplemental monomer library by
dynamically reading CCP4's monomer library."""

import os
import sys
import copy

import mmCIF
from   Library    import Monomer, Library



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


class CCP4Library(Library):
    def __init__(self, ccp4_root):
        Library.__init__(self)

        self.mon_lib_path    = os.path.join(ccp4_root,
                                            "lib",
                                            "data",
                                            "monomers")
        self.monomer_cache   = {}

    def __getitem__(self, key):
        try:
            return self.monomer_cache[key]
        except KeyError:
            pass

        try:
            mon = Library.__getitem__(self, key)
        except KeyError:
            mon = None
            
        ccp4mon = self.loadMonomer(key, copy.deepcopy(mon))
        if not ccp4mon:
            raise KeyError

        self.monomer_cache[key] = ccp4mon
        return ccp4mon
    
    def loadMonomer(self, name, mon):
        """Locates and reads the mmCIF monomer description found in the CCP4
        monomer dictionary."""
        path = os.path.join(self.mon_lib_path,
                            name[0].lower(),
                            "%s.cif" % (name))

        if not os.path.isfile(path):
            print "CCP4Library: file not found %s" % (path)
            return None

        if not mon:
            mon = Monomer(name = name)

        cfile = mmCIF.mmCIFFile()
        cfile.loadFile(CCP4MonomerFile(path, "r"))

        ldata = cfile.getData("comp_list")
        cdata = cfile.getData("comp_%s" % (name))

        full_name = ""
        if hasattr(ldata, "chem_comp"):
            for cc in ldata.chem_comp:
                if cc.id == name:
                    mon.full_name = cc.name
                    break

        if hasattr(cdata, "chem_comp_atom"):
            mon.atom_list = []
            for cca in cdata.chem_comp_atom:
                mon.atom_list.append((cca.atom_id, cca.type_symbol))

        if hasattr(cdata, "chem_comp_bond"):
            mon.bond_list = []
            for ccb in cdata.chem_comp_bond:
                mon.bond_list.append((ccb.atom_id_1, ccb.atom_id_2))

        return mon


##
## <testing>
##
if __name__ == "__main__":
    lib = CCP4Library("/home/jpaint/ccp4/ccp4-4.2.2")
    print lib["GLY"].getPolymerBondList(None, None)
##
## </testing>
##
