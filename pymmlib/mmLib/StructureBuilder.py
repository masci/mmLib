## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import PDB
import mmCIF

from Library      import Library
from Structure    import *


class StructureBuilder:
    def __init__(self,
                 fil,
                 library          = None,
                 build_properties = ()):


        ## if no custom library was defined, use the default
        if not library:
            library = Library()

        ## contstruct the Structure graph we are building
        self.structure = Structure(
            default_alt_loc = "A",
            library         = library)
        
        ## what items are going to be built into the Structure graph
        ## follow up with adding structural components which depend on
        ## other components
        self.build_properties = build_properties

        ## these caches are used to speed the contruction of the Structure
        self.atom_cache       = {}
        self.fragment_cache   = {}
        self.chain_cache      = []
        self.alt_loc_cache    = []
        
        ## invokes the file parser, builds all the atoms
        self.parseFormat(fil)

        ## build mandatory parts of the structure
        self.buildFragments()
        self.buildChains()
        
        ## build optional parts of the structure
        if "polymers"  in self.build_properties: self.buildPolymers()
        if "bonds"     in self.build_properties: self.buildBonds()
        #if "molecules" in self.build_properties: self.buildMolecules()

        ## free caches
        del self.atom_cache
        del self.fragment_cache
        del self.chain_cache
        del self.alt_loc_cache

    def addFragment(self, resName, resSeq, iCode, chainID):
        frag_id = (resName, resSeq, iCode, chainID)

        try:
            return self.fragment_cache[frag_id]
        except KeyError:
            pass

        frag = None

        if self.structure.library.isAminoAcid(resName):
            frag = AminoAcidResidue(self.structure,
                                    res_name   = resName,
                                    res_seq    = resSeq,
                                    icode      = iCode,
                                    chain_id   = chainID)

        elif self.structure.library.isNucleicAcid(resName):
            frag = NucleicAcidResidue(self.structure,
                                      res_name   = resName,
                                      res_seq    = resSeq,
                                      icode      = iCode,
                                      chain_id   = chainID)

        else:
            frag = Fragment(self.structure,
                            res_name   = resName,
                            res_seq    = resSeq,
                            icode      = iCode,
                            chain_id   = chainID)

        self.fragment_cache[frag_id] = frag

        return frag

    def addAtom(self, name, altLoc, resName, resSeq, iCode, chainID):
        atm_id = (name, altLoc, resName, resSeq, iCode, chainID)

        try:
            return self.atom_cache[atm_id]
        except KeyError:
            pass

        atm = Atom(self.structure,
                   name        = name,
                   alt_loc     = altLoc,
                   res_name    = resName,
                   res_seq     = resSeq,
                   icode       = iCode,
                   chain_id    = chainID)

        ## update caches
        self.atom_cache[atm_id] = atm

        if altLoc not in self.alt_loc_cache:
            self.alt_loc_cache.append(altLoc)
            self.alt_loc_cache.sort()
            
        if chainID not in self.chain_cache:
            self.chain_cache.append(chainID)
            self.chain_cache.sort()

        ## link alternate confirmations
        if altLoc:
            alt_loc_list = None

            for (atm_id2, atm2) in self.atom_cache.items():
                if atm_id[0]  == atm_id2[0] and \
                   atm_id[1]  != atm_id2[1] and \
                   atm_id[2:] == atm_id2[2:]:
                    alt_loc_list = atm2.alt_loc_list

            if alt_loc_list:
                atm.alt_loc_list = alt_loc_list
                atm.alt_loc_list.add(atm)
            else:
                atm.alt_loc_list = AltLocList()
                atm.alt_loc_list.add(atm)
        
        return atm

    def buildFragments(self):
        ## first sort all the atoms into a guarenteed order which will
        ## make things fast
        asort = []
        for (atm_id, atm) in self.atom_cache.items():
            (name, altLoc, resName, resSeq, iCode, chainID) = atm_id
            asort.append((chainID, resSeq, iCode, resName, name, altLoc, atm))
        asort.sort()
            
        ## create the fragments and connect them to their atoms
        for (chainID, resSeq, iCode, resName, name, altLoc, atm) in asort:
            frag = self.addFragment(resName, resSeq, iCode, chainID)
            frag.atom_list.add(atm)

    def buildChains(self):
        ## create the Chain objects and load them into a map of:
        ##   chain_map[chain_id] -> Chain object
        chain_map = {}
        for chain_id in self.chain_cache:
            chain_map[chain_id] = chain = Chain(self.structure, chain_id)
            self.structure.chain_list.add(chain)

        ## iterate through all fragments and connect them to their
        ## Chain object according to chain_id
        for frag in self.fragment_cache.values():
            chain = chain_map[frag.chain_id]
            chain.fragment_list.add(frag)
            
    def buildPolymers(self):
        ## iterate through the chains, and search each chain for segments
        ## of continous polymers; add these segments to the segment_list
        ## as the 2-tuple (chain, [frag1, frag2, frag3, ...])
        segment_list = []
        for chain in self.structure.iterChains():

            residue_class = None
            residue_list  = []

            for frag in chain.iterFragments():
                if not residue_class:
                    if isinstance(frag, Residue):
                        residue_class = frag.__class__
                        residue_list.append(frag)

                elif isinstance(frag, residue_class):
                    residue_list.append(frag)

                elif residue_list and isinstance(frag, Residue):
                    segment_list.append((chain, residue_list))
                    residue_class = frag.__class__
                    residue_list  = [frag]

                elif residue_list:
                    segment_list.append((chain, residue_list))
                    residue_class = None
                    residue_list  = []

            if residue_list:
                segment_list.append((chain, residue_list))

        ## now determine segment type and construct:
        ##  1) the correct Polymer object to contain it
        ##  2) load the Polymer and add it to the Chain.polymer_list
        for i in range(len(segment_list)):
            (chain, residue_list) = segment_list[i]

            res0 = residue_list[0]
            poly = res0.polymer_class(self.structure, chain_id=chain.chain_id)

            chain.polymer_list.add(poly)
            poly.polymer_id = chain.polymer_list.index(poly)
            
            for res in residue_list:
                res.polymer_id = poly.polymer_id
                poly.fragment_list.add(res)
                
            
    def buildBonds(self):
        def bond_atoms(atm1, atm2):
            def make_bond(a1, a2):
                bond = a1.getBond(a2)
                if bond: return bond
                return Bond(a1, a2)

            ## this handles constructing the bonds for alternate
            ## conformations correctly
            alist1 = [aatm.alt_loc for aatm in atm1 if aatm.alt_loc]
            alist2 = [aatm.alt_loc for aatm in atm2 if aatm.alt_loc]

            if not (alist1 or alist2):
                make_bond(atm1, atm2)

            elif alist1 and alist2:
                for alt_loc in alist1:
                    make_bond(atm1[alt_loc], atm2[alt_loc])

            elif alist1 and not alist2:
                for alt_loc in alist1:
                    make_bond(atm1[alt_loc], atm2)

            else:
                for alt_loc in alist2:
                    make_bond(atm1, atm2[alt_loc])  

        ## BUILD BONDS INSIDE FRAGMENTS
        for frag in self.structure.iterFragments():

            ## lookup the definition of the monomer in the library
            try:
                mon = self.structure.library[frag.res_name]
            except KeyError:
                continue

            for (name1, name2) in mon.bond_list:
                try:
                    atm1 = frag[name1]
                    atm2 = frag[name2]
                except KeyError:
                    continue

                bond_atoms(atm1, atm2)

        ## BUILD BONDS BETWEEN RESIDUES IN POLYMERS
        for poly in self.structure.iterPolymers():
            for res1 in poly.iterResidues():

                res2 = res1.getResidue(1)
                if not res2: continue

                ## lookup the definition of the monomer in the library
                try:
                    mon1 = self.structure.library[res1.res_name]
                    mon2 = self.structure.library[res2.res_name]
                except KeyError:
                    continue

                for (name1, name2) in mon.getPolymerBondList(res1, res2):
                    try:
                        atm1 = res1[name1]
                        atm2 = res2[name2]
                    except KeyError:
                        continue

                    bond_atoms(atm1, atm2)
                
    def buildMolecules(self):
        pass

    def loadAtom(self, atm_map):
        name    = atm_map["name"]
        altLoc  = atm_map["alt_loc"]
        resName = atm_map["res_name"]
        resSeq  = atm_map["res_seq"]
        iCode   = atm_map["icode"]
        chainID = atm_map["chain_id"]

        atm_id = (name, altLoc, resName, resSeq, iCode, chainID)

        ## <XXX> check the data we're loading and make sure a few
        ##       rules are followed; this is necessary, but a bit
        ##       hackish to have this here

        ## don't allow the same atom to be loaded twice
        if self.atom_cache.has_key(atm_id):
            print "[DUPLICATE ATOM ERROR]",
            print "name=%s altLoc=%s resName=%s resSeq=%s iCode=%s "\
                  "chainID=%s" % (name, altLoc, resName, resSeq,
                                  iCode, chainID)
            return

        ## don't allow atoms with blank altLoc if one has a altLoc
        if altLoc: 
            tmp_id = (name, "", resName, resSeq, iCode, chainID)
            if self.atom_cache.has_key(tmp_id):
                print "[ALTLOC LABELING ERROR]",atm_id

        ##
        ## </XXX>

        atm = self.addAtom(name, altLoc, resName, resSeq, iCode, chainID)

        ## set position
        atm.position    = Vector(atm_map["x"], atm_map["y"], atm_map["z"])

        ## set occupancy
        atm.occupancy   = atm_map["occupancy"]

        ## set temperature factor
        atm.temp_factor = atm_map["temp_factor"]

        ## try to set the U values if they exist
        try:
            atm.U = (atm_map["U[1][1]"],
                     atm_map["U[2][2]"],
                     atm_map["U[3][3]"],
                     atm_map["U[1][2]"],
                     atm_map["U[1][3]"],
                     atm_map["U[2][3]"])
        except KeyError:
            pass

        ## try to set the charge, if it exists
        try:
            atm.charge  = atm_map["charge"]
        except KeyError:
            pass


class CopyStructureBuilder(StructureBuilder):
    """Builds a new Structure object by copying from a current Structure
    object.  This builder can take any member of the Structure object which
    has a iterAtoms() method."""
    
    def parseFormat(self, fil):
        for atm in fil.iterAtoms():
            atm_map = {}

            atm_map["name"]        = atm.name
            atm_map["alt_loc"]     = atm.alt_loc
            atm_map["res_name"]    = atm.res_name
            atm_map["res_seq"]     = atm.res_seq
            atm_map["icode"]       = atm.icode
            atm_map["chain_id"]    = atm.chain_id
            
            atm_map["x"]           = atm.position[0]
            atm_map["y"]           = atm.position[1]
            atm_map["z"]           = atm.position[2]
            atm_map["occupancy"]   = atm.occupancy
            atm_map["temp_factor"] = atm.temp_factor
            
            try:
                atm_map["element"] = atm.element
            except AttributeError:
                pass
            
            try:
                atm_map["charge"]  = atm.charge
            except AttributeError:
                pass
            
            self.loadAtom(atm_map)
            

class PDBStructureBuilder(StructureBuilder):
    """Builds a new Structure object by loading a PDB file."""

    def parseFormat(self, fil):
        pdb_file = PDB.PDBFile()
        pdb_file.loadFile(fil)

        record_list = pdb_file.pdb_list

        for i in range(len(record_list)):
            rec = record_list[i]
            
            if isinstance(rec, PDB.ATOM):
                atm_map                = {}
                atm_map["name"]        = getattr(rec, "name")       or ""
                atm_map["element"]     = getattr(rec, "element")    or ""
                atm_map["alt_loc"]     = getattr(rec, "altLoc")     or ""
                atm_map["res_name"]    = getattr(rec, "resName")    or ""
                atm_map["res_seq"]     = getattr(rec, "resSeq")     or 0
                atm_map["icode"]       = getattr(rec, "iCode")      or ""
                atm_map["chain_id"]    = getattr(rec, "chainID")    or ""
                atm_map["x"]           = getattr(rec, "x")          or 0.0
                atm_map["y"]           = getattr(rec, "y")          or 0.0
                atm_map["z"]           = getattr(rec, "z")          or 0.0
                atm_map["occupancy"]   = getattr(rec, "occupancy")  or 0.0
                atm_map["temp_factor"] = getattr(rec, "tempFactor") or 0.0
                atm_map["charge"]      = getattr(rec, "charge")     or 0.0
               
##                 atm_map["U[1][1]"] = float(cifattr(aniso, "U[1][1]"))
##                 atm_map["U[2][2]"] = float(cifattr(aniso, "U[2][2]"))
##                 atm_map["U[3][3]"] = float(cifattr(aniso, "U[3][3]"))
##                 atm_map["U[1][2]"] = float(cifattr(aniso, "U[1][2]"))
##                 atm_map["U[1][3]"] = float(cifattr(aniso, "U[1][3]"))
##                 atm_map["U[2][3]"] = float(cifattr(aniso, "U[2][3]"))

                self.loadAtom(atm_map)

class mmCIFStructureBuilder(StructureBuilder):
    """Builds a new Structure object by loading a mmCIF file."""
    
    def parseFormat(self, fil):
        cif_file = mmCIF.mmCIFFile()
        cif_file.loadFile(fil)
        
        for cif_data in cif_file.getDataList():
            self.loadCIFData(cif_data)

    def loadCIFData(self, cif_data):

        def cifattr(cif_row, attr, default = ""):
            val = getattr(cif_row, attr, default)
            if val == "?" or val == ".":
                val = ""
            return val

        for atom_site in cif_data.atom_site.getRowList():
            atm_map = {}

            atm_map["name"]      = cifattr(atom_site, "label_atom_id")
            atm_map["alt_loc"]   = cifattr(atom_site, "label_alt_id")
            
            atm_map["res_name"]  = \
                cifattr(atom_site, "auth_comp_id") or \
                cifattr(atom_site, "label_comp_id")

            res_seq              = \
                cifattr(atom_site, "auth_seq_id") or \
                cifattr(atom_site, "label_seq_id")

            try:
               atm_map["res_seq"]  = int(res_seq)
            except ValueError:
                atm_map["icode"]   = res_seq[-1]
                atm_map["res_seq"] = int(res_seq[:-1])
            else:
                atm_map["icode"]   = ""

            atm_map["chain_id"]    = \
                cifattr(atom_site, "auth_asym_id") or \
                cifattr(atom_site, "label_asym_id")

            atm_map["element"]   = cifattr(atom_site, "type_symbol")

            atm_map["x"]         = \
                float(cifattr(atom_site, "Cartn_x"))

            atm_map["y"]         = \
                float(cifattr(atom_site, "Cartn_y"))

            atm_map["z"]         = \
                float(cifattr(atom_site, "Cartn_z"))

            atm_map["occupancy"] = \
                float(cifattr(atom_site, "occupancy"))

            atm_map["temp_factor"] = \
                float(cifattr(atom_site, "B_iso_or_equiv"))

            if hasattr(cif_data, "atom_site_anisotrop"):
                try:
                    (aniso, ) = cif_data.atom_site_anisotrop.selectRowList(
                        ("id", atom_site.id))

                except ValueError:
                    pass

                else:
                    atm_map["U[1][1]"] = float(cifattr(aniso, "U[1][1]"))
                    atm_map["U[2][2]"] = float(cifattr(aniso, "U[2][2]"))
                    atm_map["U[3][3]"] = float(cifattr(aniso, "U[3][3]"))
                    atm_map["U[1][2]"] = float(cifattr(aniso, "U[1][2]"))
                    atm_map["U[1][3]"] = float(cifattr(aniso, "U[1][3]"))
                    atm_map["U[2][3]"] = float(cifattr(aniso, "U[2][3]"))
                    

            self.loadAtom(atm_map)



### <TESTING>
if __name__ == "__main__":
    import sys
    struct = PDBStructureBuilder(sys.argv[1],
                                   build_properties=("polymers","bonds")
                                   ).structure

    for res in struct.iterAminoAcids():
        print "---"
        print res.getResidue(-1)
        print res
        print res.getResidue(1)
        print "---"

    s2 = CopyStructureBuilder(res,
                              build_properties=("polymers","bonds")
                              ).structure
    struct = None

    for atm in s2.iterAtoms():
        print atm

    print "# exit"
### </TESTING>

        
