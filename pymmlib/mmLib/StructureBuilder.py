## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

"""Classes for building a mmLib.Structure representation of biological
macromolecules."""

import PDB
import mmCIF

from Library      import Library
from Structure    import *


class StructureBuilder:
    """Builder class for the mmLib.Structure object hierarchy.
    StructureBuilder must be subclassed with a working parseFormat()
    method to implement a working builder."""
    
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
        
        ## invokes the file parser, builds all the atoms
        self.parseFormat(fil)

        ## build mandatory parts of the structure
        self.buildStructure()
        
        ## build optional parts of the structure
        if "sequence"  in self.build_properties: self.buildSequence()
        if "bonds"     in self.build_properties: self.buildBonds()

        ## free caches
        del self.atom_cache

    def addAtom(self, name, altLoc, resName, fragmentID, chainID):
        atm_id = (name, altLoc, fragmentID, chainID)

        try:
            return self.atom_cache[atm_id]
        except KeyError:
            pass

        atm = Atom(name        = name,
                   alt_loc     = altLoc,
                   res_name    = resName,
                   fragment_id = fragmentID,
                   chain_id    = chainID)

        self.atom_cache[atm_id] = atm

        return atm

    def buildStructure(self):
        """Construct required parts of the mmLib.Structure hierarchy."""

        def fragment_factory(atm):
            if self.structure.library.isAminoAcid(atm.res_name):
                frag_class = AminoAcidResidue
                
            elif self.structure.library.isNucleicAcid(atm.res_name):
                frag_class = NucleicAcidResidue

            else:
                frag_class = Fragment

            return frag_class(res_name    = atm.res_name,
                              fragment_id = atm.fragment_id,
                              chain_id    = atm.chain_id)

        for atm in self.atom_cache.values():
            ## retrieve/create Chain
            try:
                chain = self.structure[atm.chain_id]
            except KeyError:
                chain = Chain(atm.chain_id)
                self.structure.addChain(chain, delay_sort = True)

            ## retrieve/create Fragment
            try:
                frag = chain[atm.fragment_id]
            except KeyError:
                frag = fragment_factory(atm)
                chain.addFragment(frag, delay_sort = True)

            frag.addAtom(atm)

        ## sort structural objects into their correct order
        self.structure.sort()
        for chain in self.structure.iterChains():
            chain.sort()
        
    def buildSequence(self):
        """The residue sequence in a chain can be calculated by the
        algorithm in mmLib.Structure.Chain.calcSequence(), or by the
        sequence loaded from the input file.  This function use the
        file data if it exists, otherwise it will try to calculate
        the sequence."""
        for chain in self.structure.iterChains():
            chain.calcSequence()

    def buildBonds(self):
        """Bonds are constructed using a combination of monomer library
        bond definitions, and data loaded from parseFormat()."""
        for frag in self.structure.iterFragments():
            frag.createBonds()

    def loadInfo(self, info_map):
        """Called by the implementation of parseFormat to load descriptive
        information about the structure."""
        
        try: self.structure.id = info_map["id"]
        except KeyError: pass

        try: self.structure.date = info_map["date"]
        except KeyError: pass

        try: self.structure.keywords = info_map["keywords"]
        except KeyError: pass
       
        try: self.structure.pdbx_keywords = info_map["pdbx_keywords"]
        except KeyError: pass 
      
        try: self.structure.title = info_map["title"]
        except KeyError: pass  

        try: self.structure.R_fact = info_map["R_fact"]
        except KeyError: pass  

        try: self.structure.free_R_fact = info_map["free_R_fact"]
        except KeyError: pass  

        try: self.structure.res_high = info_map["res_high"]
        except KeyError: pass  

        try: self.structure.res_low = info_map["res_low"]
        except KeyError: pass  

    def loadAtom(self, atm_map):
        """Called by the implementation of parseFormat to load all the
        data for a single atom.  The data is contained in the atm_map
        argument, and is not well documented at this point.  Look at
        this function and you'll figure it out."""
        name       = atm_map["name"]
        altLoc     = atm_map["alt_loc"]
        resName    = atm_map["res_name"]
        fragmentID = atm_map["fragment_id"]
        chainID    = atm_map["chain_id"]

        atm_id = (name, altLoc, fragmentID, chainID)

        ## <XXX> check the data we're loading and make sure a few
        ##       rules are followed; this is necessary, but a bit
        ##       hackish to have this here

        ## don't allow the same atom to be loaded twice
        if self.atom_cache.has_key(atm_id):
            print "[DUPLICATE ATOM ERROR]", atm_id
            return

        ## don't allow atoms with blank altLoc if one has a altLoc
        if altLoc:
            tmp_id = (name, "", fragmentID, chainID)
            if self.atom_cache.has_key(tmp_id):
                print "[ALTLOC LABELING ERROR]", atm_id

        ##
        ## </XXX>

        atm = self.addAtom(name, altLoc, resName, fragmentID, chainID)

        ## additional properties
        atm.element     = atm_map["element"]
        atm.position    = Vector(atm_map["x"], atm_map["y"], atm_map["z"])
        atm.occupancy   = atm_map["occupancy"]
        atm.temp_factor = atm_map["temp_factor"]

        try:
            atm.U = (atm_map["U[1][1]"],
                     atm_map["U[2][2]"],
                     atm_map["U[3][3]"],
                     atm_map["U[1][2]"],
                     atm_map["U[1][3]"],
                     atm_map["U[2][3]"])
        except KeyError:
            pass

        try:
            atm.charge  = atm_map["charge"]
        except KeyError:
            pass

    def loadSite(self, site_id, site_list):
        """Called by the implementation of parseFormat to load information
        about one site in the structure.  Sites are groups of residues
        which are of special interest.  This usually means active sites
        of enzymes and such."""
        site = Site(site_id)
        for (chain_id, frag_id) in site_list:
            site.addFragment(chain_id, frag_id)
        self.structure.sites.append(site)


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
            atm_map["fragment_id"] = atm.fragment_id
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

    def __process_ATOM(self, atm_map, rec):
        name    = getattr(rec, "name")       or ""
        element = getattr(rec, "element")    or ""

        ## get the element symbol from the first letter in the
        ## atom name
        if not self.structure.library.element_map.has_key(element):
            for c in name:
                if c in string.letters:
                    element = c
                    break

        atm_map["name"]        = name
        atm_map["element"]     = element
        atm_map["alt_loc"]     = getattr(rec, "altLoc")     or ""
        atm_map["res_name"]    = getattr(rec, "resName")    or ""
        atm_map["fragment_id"] = str(getattr(rec, "resSeq") or "") + \
                                 (getattr(rec, "iCode")      or "")
        atm_map["chain_id"]    = getattr(rec, "chainID")    or ""
        atm_map["x"]           = getattr(rec, "x")          or 0.0
        atm_map["y"]           = getattr(rec, "y")          or 0.0
        atm_map["z"]           = getattr(rec, "z")          or 0.0
        atm_map["occupancy"]   = getattr(rec, "occupancy")  or 0.0
        atm_map["temp_factor"] = getattr(rec, "tempFactor") or 0.0
        atm_map["charge"]      = getattr(rec, "charge")     or 0.0

    def __process_ANISOU(self, atm_map, rec):
        atm_map["U[1][1]"] = getattr(rec, "u[0][0]") / 10000.0
        atm_map["U[2][2]"] = getattr(rec, "u[1][1]") / 10000.0
        atm_map["U[3][3]"] = getattr(rec, "u[2][2]") / 10000.0
        atm_map["U[1][2]"] = getattr(rec, "u[0][1]") / 10000.0
        atm_map["U[1][3]"] = getattr(rec, "u[0][2]") / 10000.0
        atm_map["U[2][3]"] = getattr(rec, "u[1][2]") / 10000.0

    def __process_SITE(self, site_map, rec):
        site_id = getattr(rec, "siteID")

        if not site_map.has_key(site_id):
            site_map[site_id] = []

        try:
            chain_id = rec.chainID1
            frag_id  = str(rec.seq1) + getattr(rec, "icode1", "")
        except AttributeError:
            return
        else:
            site_map[site_id].append((chain_id, frag_id))

        try:
            chain_id = rec.chainID2
            frag_id  = str(rec.seq2) + getattr(rec, "icode2", "")
        except AttributeError:
            return
        else:
            site_map[site_id].append((chain_id, frag_id))

        try:
            chain_id = rec.chainID3
            frag_id  = str(rec.seq3) + getattr(rec, "icode3", "")
        except AttributeError:
            return
        else:
            site_map[site_id].append((chain_id, frag_id))

        try:
            chain_id = rec.chainID4
            frag_id  = str(rec.seq4) + getattr(rec, "icode4", "")
        except AttributeError:
            return
        else:
            site_map[site_id].append((chain_id, frag_id))

    
    def parseFormat(self, fil):
        pdb_file = PDB.PDBFile()
        pdb_file.loadFile(fil)

        site_map      = {}

        atm_map       = {}
        record_list   = pdb_file.pdb_list

        for i in range(len(record_list)):
            rec = record_list[i]
            
            if isinstance(rec, PDB.ATOM):
                ## load the last atom before moving to this one
                if atm_map:
                    self.loadAtom(atm_map)
                    atm_map = {}

                try:
                    self.__process_ATOM(atm_map, rec)
                except:
                    print "ERROR WITH ATOM RECORD:"
                    print rec

            elif isinstance(rec, PDB.ANISOU):
                try:
                    self.__process_ANISOU(atm_map, rec)
                except:
                    print "ERROR WITH ANISOU RECORD:"
                    print rec

            elif isinstance(rec, PDB.SITE):
                self.__process_SITE(site_map, rec)
                
        ## if we were in the process of building a atom,
        ## then load the final atm_map 
        if atm_map:
            self.loadAtom(atm_map)

        ## load the SITE records collected during parsing
        for (site_id, site_list) in site_map.items():
            self.loadSite(site_id, site_list)


class mmCIFStructureBuilder(StructureBuilder):
    """Builds a new Structure object by loading a mmCIF file."""
    
    def parseFormat(self, fil):
        cif_file = mmCIF.mmCIFFile()
        cif_file.loadFile(fil)
        
        for cif_data in cif_file:
            self.__load_CIF_Data(cif_data)

    def __load_CIF_Data(self, cif_data):

        def cifattr(cif_row, attr, default = ""):
            val = cif_row.get(attr, default)
            if val == "?" or val == ".":
                val = ""
            return val

        ## read structure header/title info
        info_map = {}

        try: info_map["id"] = cif_data["entry"][0]["id"]
        except KeyError: print "missing entry.id"

        try: info_map["date"] = \
             cif_data["database_pdb_rev"][0]["date_original"]
        except KeyError: print "missing database_pdb_rev.date_original"

        try: info_map["keywords"] = cif_data["struct_keywords"][0]["text"]
        except KeyError: print "missing struct_keywords.text"

        try: info_map["pdbx_keywords"] = \
             cif_data["struct_keywords"][0]["pdbx_keywords"]
        except KeyError: print "missing struct_keywords.pdbx_keywords"

        try: info_map["title"] = cif_data["struct"][0]["title"]
        except KeyError: print "missing struct.title"

        try: info_map["R_fact"] = \
            float(cif_data["refine"][0]["ls_R_factor_R_work"])
        except KeyError:   print "missing refine.ls_R_factor_R_work"
        except ValueError: print "missing refine.ls_R_factor_R_work"
        
        try: info_map["free_R_fact"] = \
             float(cif_data["refine"][0]["ls_R_factor_R_free"])
        except KeyError:   print "missing refine.ls_R_factor_R_free"
        except ValueError: print "missing refine.ls_R_factor_R_free"

        try: info_map["res_high"] = \
             float(cif_data["refine"][0]["ls_d_res_high"])
        except KeyError:   print "missing refine.ls_d_res_high"
        except ValueError: print "missing refine.ls_d_res_high"

        try: info_map["res_low"] = \
             float(cif_data["refine"][0]["ls_d_res_low"])
        except KeyError:   print "missing refine.ls_d_res_low"
        except ValueError: print "missing refine.ls_d_res_low"
        
        self.loadInfo(info_map)

        ## read atom coordinate data
        for atom_site in cif_data["atom_site"]:
            atm_map = {}

            atm_map["name"]      = cifattr(atom_site, "label_atom_id")
            atm_map["alt_loc"]   = cifattr(atom_site, "label_alt_id")
            
            atm_map["res_name"]  = \
                cifattr(atom_site, "auth_comp_id") or \
                cifattr(atom_site, "label_comp_id")

            atm_map["fragment_id"] = \
                cifattr(atom_site, "auth_seq_id") or \
                cifattr(atom_site, "label_seq_id")

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

            if cif_data.has_key("atom_site_anisotrop"):
                ctable = cif_data["atom_site_anisotrop"]
                try:
                    (aniso, ) = ctable.selectRowList(("id", atom_site["id"]))
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


        ## load SITE data
        if cif_data.has_key("struct_site_gen"):
            site_map = {}
        
            for struct_site_gen in cif_data["struct_site_gen"]:
                ## extract data for chain_id and seq_id
                site_id  = struct_site_gen["site_id"]

                chain_id = struct_site_gen["auth_asym_id"]
                if chain_id in ["?", "."]:
                    chain_id = struct_site_gen["label_asym_id"]

                frag_id  = struct_site_gen["auth_seq_id"]
                if frag_id in ["?", "."]:
                    frag_id = struct_site_gen["label_seq_id"]

                ## skip bad rows
                if chain_id in ["?", "."] or frag_id in ["?", "."]:
                    continue
                    
                if not site_map.has_key(site_id):
                    site_map[site_id] = []

                site_map[site_id].append((chain_id, frag_id))

            for (site_id, site_list) in site_map.items():
                self.loadSite(site_id, site_list)
            

### <TESTING>
if __name__ == "__main__":
    import sys
    struct = PDBStructureBuilder(sys.argv[1],
                                 build_properties=()
                                 ).structure

    for atm in struct.iterAtoms():
        if atm.alt_loc:
            print atm

    print "# exit"
### </TESTING>

        
