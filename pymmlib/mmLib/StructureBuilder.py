## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Classes for building a mmLib.Structure representation of biological
macromolecules.
"""
from mmTypes import *
import PDB
import mmCIF
from Library import Library
from Structure import *


class StructureBuilder:
    """Builder class for the mmLib.Structure object hierarchy.
    StructureBuilder must be subclassed with a working parse_format()
    method to implement a working builder.
    """
    def __init__(self, fil, library = None, build_properties = ()):
        ## contstruct the Structure graph we are building
        self.structure = Structure(library = library)
        
        ## what items are going to be built into the Structure graph
        ## follow up with adding structural components which depend on
        ## other components
        self.build_properties = build_properties

        ## if anything goes wrong, setting self.halt=1 will stop the madness
        self.halt = 0

        ## build the structure by executing this fixed sequence of methods
        self.read_start(fil)
        if not self.halt: self.read_start_finalize()
        if not self.halt: self.read_atoms()
        if not self.halt: self.read_atoms_finalize()
        if not self.halt: self.read_metadata()
        if not self.halt: self.read_metadata_finalize()
        if not self.halt: self.read_end()
        if not self.halt: self.read_end_finalize()
        ## self.structure is now built and ready for use

    def read_start(self, fil):
        """This methods needs to be reimplemented in a functional subclass.
        This function is called with the file object (or any other object
        passed in to build a Structure from) to begin the reading process.
        This is usually used to open the source file.
        """
        pass

    def read_start_finalize(self):
        """Called after the read_start method.  Does nothing currently,
        but may be used in the future.
        """
        pass

    def read_atoms(self):
        """This method needs to be reimplemented in a fuctional subclass.
        The subclassed read_atoms method should call load_atom once for
        every atom in the sturcture, and should not call any other
        load_* methods.
        """
        pass

    def load_atom(self, atm_map):
        """Called repeatedly by the implementation of read_atoms to
        load all the data for a single atom.  The data is contained
        in the atm_map argument, and is not well documented at this
        point.  Look at this function and you'll figure it out.
        """
        ## XXX -- I presently do not support more than one NMR
        ##        style MODEL; this is first on the list for the
        ##        next version
        if atm_map.has_key("model_num") and atm_map["model_num"] > 1:
            debug("NMR-style multi-models not supported yet")
            return
        ## /XXX

        ## <really, really, required>
        name = atm_map["name"]
        fragment_id = atm_map["fragment_id"]
        chain_id = atm_map["chain_id"]
        ### </really, really, required>
        
        res_name = atm_map.get("res_name", "")
        alt_loc = atm_map.get("alt_loc", "")

        ## form the cache ID for the atom and
        ## allocate the atom cache if it does not exist
        atm_id = (name, alt_loc, fragment_id, chain_id)
        if not hasattr(self, "atom_cache"):
            self.atom_cache = {}

        ## don't allow the same atom to be loaded twice
        if self.atom_cache.has_key(atm_id):
            debug("duplicate atom "+str(atm_id))
            return

        ## don't allow atoms with blank altLoc if one has a altLoc
        if alt_loc and self.atom_cache.has_key((name,"",fragment_id,chain_id)):
            debug("altloc labeling error "+name+fragment_id+chain_id)
            return

        atm = self.atom_cache[atm_id] = \
              Atom(name = name, alt_loc = alt_loc, res_name = res_name,
                   fragment_id = fragment_id, chain_id = chain_id)

        ## extract element symbol from atom name if it isn't given
        try:
            atm.element = atm_map["element"]
        except KeyError:
            for c in atm.name:
                if c in string.letters:
                    atm.element = c
                    break

        try: atm.position = Vector(atm_map["x"],atm_map["y"],atm_map["z"])
        except KeyError: pass

        try: atm.occupancy = atm_map["occupancy"]
        except KeyError: pass

        try: atm.temp_factor = atm_map["temp_factor"]
        except KeyError: pass

        try: atm.charge  = atm_map["charge"]
        except KeyError: pass

        try:
            atm.sig_position = Vector(
                atm_map["sig_x"],atm_map["sig_y"],atm_map["sig_z"])
        except KeyError:
            pass

        try:
            atm.set_U(atm_map["U[1][1]"], atm_map["U[2][2]"],
                      atm_map["U[3][3]"], atm_map["U[1][2]"],
                      atm_map["U[1][3]"], atm_map["U[2][3]"])
        except KeyError:
            pass

        try:
            atm.set_sig_U(atm_map["sig_U[1][1]"], atm_map["sig_U[2][2]"],
                          atm_map["sig_U[3][3]"], atm_map["sig_U[1][2]"],
                          atm_map["sig_U[1][3]"], atm_map["sig_U[2][3]"])
        except KeyError:
            pass

    def read_atoms_finalize(self):
        """After loading all atom records, use the list of atom records to
        build the structure.
        """
        def fragment_factory(atm):
            if self.structure.library.is_amino_acid(atm.res_name):
                frag_class = AminoAcidResidue
                
            elif self.structure.library.is_nucleic_acid(atm.res_name):
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
                self.structure.add_chain(chain, delay_sort = True)

            ## retrieve/create Fragment
            try:
                frag = chain[atm.fragment_id]
            except KeyError:
                frag = fragment_factory(atm)
                chain.add_fragment(frag, delay_sort = True)

            frag.add_atom(atm)

        ## sort structural objects into their correct order
        self.structure.sort()
        for chain in self.structure.iter_chains():
            chain.sort()

        ## iterate through all the atoms, and choose a default alt_loc
        alt_loc_list = []
        for atm1 in self.structure.iter_atoms():
            for atm2 in atm1:
                if atm2.alt_loc and atm2.alt_loc not in alt_loc_list:
                    alt_loc_list.append(atm2.alt_loc)
        if alt_loc_list:
            alt_loc_list.sort()
            self.structure.default_alt_loc = alt_loc_list[0]
        
        ## we're done with the atom cache, delete it
        if hasattr(self, "atom_cache"):
            del self.atom_cache
                    
    def read_metadata(self):
        """This method needs to be reimplemented in a fuctional subclass.
        The subclassed read_metadata method should call the various
        load_* methods to set non-atom coordinate data for the Structure.
        """
        pass

    def load_info(self, info_map):
        """Called by the implementation of parse_format to load descriptive
        information about the structure.
        """
        def set_ed(smap, skey, dkey):
            try:
                self.structure.exp_data[dkey] = smap[skey]
            except KeyError:
                pass

        ## ID of structure
        set_ed(info_map, "id", "id")

        ## source file object: pdb_file or cif_file
        set_ed(info_map, "pdb_file", "pdb_file")
        set_ed(info_map, "cif_file", "cif_file")

        ## try to determine what kind of experiment produced the data
        set_ed(info_map, "exp_method", "exp_method")

        ## structure data: Generic
        set_ed(info_map, "title", "title")
        set_ed(info_map, "date", "date")
        set_ed(info_map, "author_list", "author_list")
        set_ed(info_map, "method", "method")

        ## structure data: X-Ray Generic
        set_ed(info_map, "R_fact", "R_fact")
        set_ed(info_map, "free_R_fact", "free_R_fact")
        set_ed(info_map, "res_high", "res_high")
        set_ed(info_map, "res_low", "res_low")
        set_ed(info_map, "", "")
        set_ed(info_map, "", "")
        
        ## structure date: May Be File-Format Specific (Questionable)
        set_ed(info_map, "keywords", "keywords")
        set_ed(info_map, "pdbx_keywords", "pdbx_keywords")

    def load_unit_cell(self, ucell_map):
        """Called by the implementation of load_metadata to load the
        unit cell pararameters for the structure.
        """
        self.structure.unit_cell = UnitCell(
            ucell_map["a"],     ucell_map["b"],    ucell_map["c"],
            ucell_map["alpha"], ucell_map["beta"], ucell_map["gamma"])

        try: self.structure.space_group = SpaceGroup(ucell_map["space_group"])
        except KeyError: pass

    def load_site(self, site_map):
        """Called by the implementation of load_metadata to load information
        about one site in the structure.  Sites are groups of residues
        which are of special interest.  This usually means active sites
        of enzymes and such.  The structure of the site_map is this: the
        keys are the ID of the site (integer), and the values are a list
        of maps with values under the keys "id"(redundent), "chain_id",
        and "fragment_id".
        """
        for (site_id, sentry_list) in site_map.items():
            site = Site(name = site_id)
            self.structure.sites.append(site)

            for sentry in sentry_list:
                chain_id = sentry["chain_id"]
                fragment_id = sentry["fragment_id"]
                try:
                    frag = self.structure[chain_id][fragment_id]
                except KeyError:
                    debug("load_site: chain_id=%s fragment_id=%s not found"%(
                        chain_id, fragment_id))
                else:
                    site.append(frag)

    def load_bond(self, bond_map):
        """Call by the implementation of load_metadata to load bond
        information on the structure.
        """
        pass

    def read_metadata_finalize(self):
        """Called after the the metadata loading is complete.
        """
        ## calculate sequences for all chains
        for chain in self.structure.iter_chains():
            chain.calc_sequence()

        ## build bonds within structure
        for frag in self.structure.iter_fragments():
            frag.create_bonds()

    def read_end(self):
        """This method needs to be reimplemented in a fuctional subclass.
        The subclassed read_end method can be used for any clean up from
        the file loading process you need, or may be left unimplemented.
        """
        pass

    def read_end_finalize(self):
        """Called for final cleanup after structure source readinging is
        done.  Currently, this method does nothing but may be used in
        future versions.
        """
        pass

    def setmap(self, conv_func, smap, skey, dmap, dkey, default = None):
        """The setmap methods are meant to help with the creation of the
        maps passed into the load_* methods of this class.  The smap[skey]
        is the data source, and dmap[dkey] is the data destination.  If
        a default value is given, then the dmap[dkey] is set to that value
        if the smap[skey] does not exist.  The conv_func is used to convert
        the smap[skey] data before setting dmap[dkey].  The skey can also
        be a list of keys to try instead of a single key.  Returns True
        if dmap[dkey] is set, otherwise this method returns False.
        """
        src_key = None

        if type(skey) == StringType:
            if smap.has_key(skey):
                src_key = skey

        elif type(skey) == ListType:
            for key in skey:
                if smap.has_key(key):
                    src_key = key
                    break

        if src_key != None:
            data = smap[src_key]
            try:
                data = conv_func(data)
            except ValueError:
                debug("setmap() bad value for conv_func: %s" % (data))
                return False
            dmap[dkey] = data
            return True

        elif default != None:
            dmap[dkey] = default
            return True

        return False

    def setmaps(self, smap, skey, dmap, dkey, default = None):
        """Calls setmap with a conv_func of str().
        """
        return self.setmap(str, smap, skey, dmap, dkey, default)

    def setmapi(self, smap, skey, dmap, dkey, default = None):
        """Calls setmap with a conv_func of int().
        """
        return self.setmap(int, smap, skey, dmap, dkey, default)

    def setmapf(self, smap, skey, dmap, dkey, default = None):
        """Calls setmap with a conv_func of float().
        """
        return self.setmap(float, smap, skey, dmap, dkey, default)


class CopyStructureBuilder(StructureBuilder):
    """Builds a new Structure object by copying from a current Structure
    object.  This builder can take any member of the Structure object which
    has a iter_atoms() method.
    """
    def read_start(self, fil):
        self.copy_struct = fil

    def read_atoms(self):
        for atm in self.copy_struct.iter_atoms():
            atm_map = {}
            atm_map["name"] = atm.name
            atm_map["alt_loc"] = atm.alt_loc
            atm_map["res_name"] = atm.res_name
            atm_map["fragment_id"] = atm.fragment_id
            atm_map["chain_id"] = atm.chain_id
            atm_map["x"] = atm.position[0]
            atm_map["y"] = atm.position[1]
            atm_map["z"] = atm.position[2]
            atm_map["occupancy"] = atm.occupancy
            atm_map["temp_factor"] = atm.temp_factor
            try: atm_map["element"] = atm.element
            except AttributeError: pass
            try: atm_map["charge"]  = atm.charge
            except AttributeError: pass
            self.load_atom(atm_map)
            
    def read_metadata(self):
        pass


class PDBStructureBuilder(StructureBuilder):
    """Builds a new Structure object by loading a PDB file.
    """
    def read_start(self, fil):
        self.pdb_file = PDB.PDBFile()
        self.pdb_file.load_file(fil)
        self.structure.exp_data["pdb_file"] = self.pdb_file

    def read_atoms(self):
        model_num = None
        atm_map = {}        

        for rec in self.pdb_file:
            if isinstance(rec, PDB.ATOM):
                if atm_map:
                    self.load_atom(atm_map)
                    atm_map = {}
                try:
                    self.process_ATOM(atm_map, rec)
                except:
                    debug("bad PDB.ATOM record: "+str(rec))
                else:
                    if model_num != None:
                        atm_map["model_num"] = model_num

            elif isinstance(rec, PDB.SIGATM):
                self.process_SIGATM(atm_map, rec)

            elif isinstance(rec, PDB.ANISOU):
                try:
                    self.process_ANISOU(atm_map, rec)
                except:
                    debug("bad PDB.ANISOU record: "+str(rec))

            elif isinstance(rec, PDB.SIGUIJ):
                self.process_SIGUIJ(atm_map, rec)

            elif isinstance(rec, PDB.MODEL):
                model_num = rec.get("serial")

            elif isinstance(rec, PDB.ENDMDL):
                model_num = None

        ## load last atom read
        if atm_map:
            self.load_atom(atm_map)

    def read_metadata(self):
        info_map = {}
        site_map = {}
        ucell_map = {}

        ## gather metadata
        for rec in self.pdb_file:
            if isinstance(rec, PDB.SITE):
                self.process_SITE(site_map, rec)

            elif isinstance(rec, PDB.CRYST1):
                self.process_CRYST1(ucell_map, rec)

            elif isinstance(rec, PDB.HEADER):
                self.process_HEADER(info_map, rec)

            elif isinstance(rec, PDB.TITLE):
                self.process_TITLE(info_map, rec)

            elif isinstance(rec, PDB.COMPND):
                self.process_COMPND(info_map, rec)

            elif isinstance(rec, PDB.AUTHOR):
                self.process_AUTHOR(info_map, rec)

        ## load info
        info_map["pdb_file"] = self.pdb_file
        self.load_info(info_map)

        ## load the SITE information if found
        if site_map:
            self.load_site(site_map)

        ## load the unit cell parameters if found
        if ucell_map:
            self.load_unit_cell(ucell_map)

    def process_ATOM(self, atm_map, rec):
        self.setmaps(rec, "name", atm_map, "name")
        self.setmaps(rec, "element", atm_map, "element")
        self.setmaps(rec, "altLoc", atm_map, "alt_loc")
        self.setmaps(rec, "resName", atm_map, "res_name")
        self.setmaps(rec, "chainID", atm_map, "chain_id")

        fragment_id = self.get_fragment_id(rec)
        if fragment_id != None:
            atm_map["fragment_id"] = fragment_id
            
        self.setmapf(rec, "x", atm_map, "x")
        self.setmapf(rec, "y", atm_map, "y")
        self.setmapf(rec, "z", atm_map, "z")
        self.setmapf(rec, "occupancy", atm_map, "occupancy")
        self.setmapf(rec, "tempFactor", atm_map, "temp_factor")
        self.setmapf(rec, "charge", atm_map, "charge")

    def process_SIGATM(self, atm_map, rec):
        self.setmapf(rec, "sigX", atm_map, "sig_x")
        self.setmapf(rec, "sigY", atm_map, "sig_y")
        self.setmapf(rec, "sigZ", atm_map, "sig_z")
        self.setmapf(rec, "sigOccupancy", atm_map, "sig_occupancy")
        self.setmapf(rec, "sigTempFactor", atm_map, "sig_temp_factor")

    def process_ANISOU(self, atm_map, rec):
        atm_map["U[1][1]"] = rec.get("u[0][0]", 0.0) / 10000.0
        atm_map["U[2][2]"] = rec.get("u[1][1]", 0.0) / 10000.0
        atm_map["U[3][3]"] = rec.get("u[2][2]", 0.0) / 10000.0
        atm_map["U[1][2]"] = rec.get("u[0][1]", 0.0) / 10000.0
        atm_map["U[1][3]"] = rec.get("u[0][2]", 0.0) / 10000.0
        atm_map["U[2][3]"] = rec.get("u[1][2]", 0.0) / 10000.0

    def process_SIGUIJ(self, atm_map, rec):
        atm_map["sig_U[1][1]"] = rec.get("sig[1][1]", 0.0) / 10000.0
        atm_map["sig_U[2][2]"] = rec.get("sig[2][2]", 0.0) / 10000.0
        atm_map["sig_U[3][3]"] = rec.get("sig[3][3]", 0.0) / 10000.0
        atm_map["sig_U[1][2]"] = rec.get("sig[1][2]", 0.0) / 10000.0
        atm_map["sig_U[1][3]"] = rec.get("sig[1][3]", 0.0) / 10000.0
        atm_map["sig_U[2][3]"] = rec.get("sig[2][3]", 0.0) / 10000.0

    def process_CRYST1(self, ucell_map, rec):
        ucell_map["a"] = rec["a"]
        ucell_map["b"] = rec["b"]
        ucell_map["c"] = rec["c"]
        ucell_map["alpha"] = rec["alpha"]
        ucell_map["beta"] = rec["beta"]
        ucell_map["gamma"] = rec["gamma"]
        ucell_map["space_group"] = rec["sgroup"]
        ucell_map["z"] = rec.get("z", "")

    def process_SITE(self, site_map, rec):
        for i in range(1, 5):
            chain_key = "chainID%d" % (i)
            frag_key = "seq%d" % (i)
            icode_key = "icode%d" % (i)

            if not rec.has_key(chain_key) or \
               not rec.has_key(frag_key) or \
               not rec.has_key(icode_key):
                break

            sentry = {}
            self.setmaps(rec, "siteID", sentry, "id")
            self.setmaps(rec, chain_key, sentry, "chain_id")
            sentry["fragment_id"]=self.get_fragment_id(rec,frag_key,icode_key)

            try:
                site_map[sentry["id"]].append(sentry)
            except KeyError:
                site_map[sentry["id"]] = [sentry]

    def process_SSBOND(self, bond_map, rec):
        bond_map["type"] = "disulfide"
        bond_map["chain_id1"] = rec["chainID1"]
        bond_map["fragment_id1"] = self.get_fragment_id(rec,"seqNum1","iCode1")
        bond_map["chain_id2"] = rec["chainID2"]
        bond_map["fragment_id2"] = self.get_fragment_id(rec,"seqNum2","iCode2")

        self.setmaps(rec, "sym1", bond_map, "symop1")
        self.setmaps(rec, "sym2", bond_map, "symop2")

    def process_LINK(self, bond_map, rec):
        pass

    def process_HEADER(self, info_map, rec):
        self.setmaps(rec, "idCode", info_map, "id")
        self.setmaps(rec, "depDate", info_map, "date")
        self.setmaps(rec, "classification", info_map, "pdbx_keywords")

    def process_TITLE(self, info_map, rec):
        self.setmaps(rec, "title", info_map, "title")

    def process_COMPND(self, info_map, rec):
        self.setmaps(rec, "title", info_map, "title")

    def process_AUTHOR(self, info_map, rec):
        if rec.has_key("authorList"):
            for auth in rec["authorList"].split(","):
                if not auth:
                    continue
                try:
                    info_map["author_list"].append(auth.strip())
                except KeyError:
                    info_map["author_list"] = [auth.strip()]

    def get_fragment_id(self, rec, res_seq = "resSeq", icode = "iCode"):
        fragment_id = None
        
        if rec.has_key(res_seq):
            fragment_id = str(rec[res_seq])
            if rec.has_key(icode):
                fragment_id += rec[icode]
                
        return fragment_id


class mmCIFStructureBuilder(StructureBuilder):
    """Builds a new Structure object by loading a mmCIF file.
    """
    def setmap(self, conv_func, smap, skey, dmap, dkey, default = None):
        """mmCIF files have those annoying question marks which we want
        to treat as no data.
        """
        src_key = None
        
        if type(skey) == StringType:
            if smap.has_key(skey) and smap[skey] != "?":
                src_key = skey

        elif type(skey) == ListType:
            for key in skey:
                if smap.has_key(key) and smap[key] != "?":
                    src_key = key
                    break

        if src_key != None:
            data = smap[src_key]
            try:
                data = conv_func(data)
            except ValueError:
                debug("setmap() bad value for conv_func: %s" % (data))
                return False
            dmap[dkey] = data
            return True

        elif default != None:
            dmap[dkey] = default
            return True

        return False

    def read_start(self, fil):
        self.cif_file = mmCIF.mmCIFFile()
        self.cif_file.load_file(fil)

        ## only the first data block in the mmCIF file will be read
        self.cif_data = self.cif_file[0]

    def read_atoms(self):
        for atom_site in self.cif_data["atom_site"]:
            atm_map = {}
            
            self.setmaps(atom_site, ["auth_atom_id","label_atom_id"],
                         atm_map, "name")

            self.setmaps(atom_site, ["auth_alt_id","label_alt_id"],
                         atm_map, "alt_loc")

            self.setmaps(atom_site, ["auth_comp_id","label_comp_id"],
                         atm_map, "res_name")

            self.setmaps(atom_site, ["auth_seq_id","label_seq_id"],
                         atm_map, "fragment_id")

            self.setmaps(atom_site, ["auth_asym_id","label_asym_id"],
                         atm_map, "chain_id")

            self.setmaps(atom_site, "type_symbol", atm_map, "element")

            self.setmapf(atom_site, "Cartn_x", atm_map, "x")
            self.setmapf(atom_site, "Cartn_y", atm_map, "y")
            self.setmapf(atom_site, "Cartn_z", atm_map, "z")
            self.setmapf(atom_site, "occupancy", atm_map, "occupancy")
            self.setmapf(atom_site, "B_iso_or_equiv", atm_map, "temp_factor")

            self.setmapf(atom_site, "Cartn_x_esd", atm_map, "sig_x")
            self.setmapf(atom_site, "Cartn_y_esd", atm_map, "sig_y")
            self.setmapf(atom_site, "Cartn_z_esd", atm_map, "sig_z")
            self.setmapf(atom_site, "occupancy_esd", atm_map, "sig_occupancy")
            self.setmapf(atom_site, "B_iso_or_equiv_esd",
                         atm_map, "sig_temp_factor")
            
            self.setmapi(atom_site, "pdbx_PDB_model_num", atm_map, "model_num")

            if self.cif_data.has_key("atom_site_anisotrop"):
                ctable = self.cif_data["atom_site_anisotrop"]
                try:
                    (aniso, ) = ctable.select_row_list(("id", atom_site["id"]))
                except ValueError:
                    pass
                else:
                    self.setmapf(aniso, "U[1][1]", atm_map, "U[1][1]")
                    self.setmapf(aniso, "U[2][2]", atm_map, "U[2][2]")
                    self.setmapf(aniso, "U[3][3]", atm_map, "U[3][3]")
                    self.setmapf(aniso, "U[1][2]", atm_map, "U[1][2]")
                    self.setmapf(aniso, "U[1][3]", atm_map, "U[1][3]")
                    self.setmapf(aniso, "U[2][3]", atm_map, "U[2][3]")

                    self.setmapf(aniso, "U[1][1]_esd", atm_map, "sig_U[1][1]")
                    self.setmapf(aniso, "U[2][2]_esd", atm_map, "sig_U[2][2]")
                    self.setmapf(aniso, "U[3][3]_esd", atm_map, "sig_U[3][3]")
                    self.setmapf(aniso, "U[1][2]_esd", atm_map, "sig_U[1][2]")
                    self.setmapf(aniso, "U[1][3]_esd", atm_map, "sig_U[1][3]")
                    self.setmapf(aniso, "U[2][3]_esd", atm_map, "sig_U[2][3]")

            self.load_atom(atm_map)

    def read_metadata(self):
        def set_im(conv_func, ctable, col, dmap, dkey):
            try:
                dmap[dkey] = conv_func(self.cif_data[ctable][0][col])
            except KeyError:
                return False
            except IndexError:
                return False
            except ValueError:
                return False
            return True
        def set_ims(ctable, col, dmap, dkey):
            return set_im(str, ctable, col, dmap, dkey)
        def set_imi(ctable, col, dmap, dkey):
            return set_im(int, ctable, col, dmap, dkey)
        def set_imf(ctable, col, dmap, dkey):
            return set_im(float, ctable, col, dmap, dkey)
        
        ## entry.id is required, really
        info_map = {}
        if not set_ims("entry", "id", info_map, "id"):
            debug("missing entry.id, read_metatadata cannot continue")
            return

        ## reference to source mmCIFFile object
        info_map["cif_file"] = self.cif_file

        ## generic data
        set_ims("database_pdb_rev", "date_original", info_map, "date")
        set_ims("struct_keywords", "text", info_map, "keywords")
        set_ims("struct_keywords", "pdbx_keywords", info_map, "pdbx_keywords")
        set_ims("struct", "title", info_map, "title")
        
        ## X-ray experimental data
        set_imf("refine", "ls_R_factor_R_work", info_map, "R_fact")
        set_imf("refine", "ls_R_factor_R_free", info_map, "free_R_fact")
        set_imf("refine", "ls_d_res_high", info_map, "res_high")
        set_imf("refine", "ls_d_res_low", info_map, "res_low")

        self.load_info(info_map)

        ## SITE
        if self.cif_data.has_key("struct_site_gen"):
            site_map = {}
        
            for struct_site_gen in self.cif_data["struct_site_gen"]:
                sentry = {}
                self.setmaps(struct_site_gen, "site_id", sentry, "id")
                self.setmaps(struct_site_gen, ["auth_asym_id","label_asym_id"],
                             sentry, "chain_id")
                self.setmaps(struct_site_gen, ["auth_seq_id","label_seq_id"],
                             sentry, "fragment_id")

                try:
                    site_map[sentry["id"]].append(sentry)
                except KeyError:
                    site_map[sentry["id"]] = [sentry]
                    
            self.load_site(site_map)

        ## UNIT CELL
        ucell_map = {}
        if self.cif_data.has_key("cell"):
            try:
                (cell,) = self.cif_data["cell"].select_row_list(
                    ("entry_id", info_map["id"]))
            except ValueError:
                pass
            else:
                self.setmapf(cell, "length_a", ucell_map, "a")
                self.setmapf(cell, "length_b", ucell_map, "b")
                self.setmapf(cell, "length_c", ucell_map, "c")
                self.setmapf(cell, "angle_alpha", ucell_map, "alpha")
                self.setmapf(cell, "angle_beta", ucell_map, "beta")
                self.setmapf(cell, "angle_gamma", ucell_map, "gamma")
                self.setmapi(cell, "Z_PDB", ucell_map, "z")

        if self.cif_data.has_key("symmetry"):
            try:
                (symm,) = self.cif_data["symmetry"].select_row_list(
                    ("entry_id", info_map["id"]))
            except ValueError:
                pass
            else:
                self.setmaps(symm, "space_group_name_H-M",
                             ucell_map, "space_group")
        
        if ucell_map:
            self.load_unit_cell(ucell_map)
            

### <TESTING>
if __name__ == "__main__":
    import sys
    struct = PDBStructureBuilder(sys.argv[1],
                                 build_properties=()
                                 ).structure

    for atm in struct.iter_atoms():
        if atm.alt_loc:
            print atm

    print "# exit"
### </TESTING>

        
