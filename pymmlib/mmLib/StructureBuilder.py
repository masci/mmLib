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
    def __init__(self,
                 fil = None,
                 struct = None,
                 library = None,
                 build_properties = ()):

        ## contstruct the Structure graph we are building
        self.struct = struct or Structure(library = library)

        ## what items are going to be built into the Structure graph
        ## follow up with adding structural components which depend on
        ## other components
        self.build_properties = build_properties

        ## caches used while building
        self.cache_chain = None
        self.cache_frag = None

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
        ## self.struct is now built and ready for use

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
        self.name_service_list = []


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
        ## required items
        name = atm_map.get("name", "")
        model = atm_map.get("model_num", 1)
        fragment_id = atm_map.get("fragment_id", "")
        chain_id = atm_map.get("chain_id", "")
        res_name = atm_map.get("res_name", "")
        alt_loc = atm_map.get("alt_loc", "")

        atm = Atom(name = name,
                   model = model,
                   alt_loc = alt_loc,
                   res_name = res_name,
                   fragment_id = fragment_id,
                   chain_id = chain_id)

        try: atm.element = atm_map["element"]
        except KeyError: pass

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

        ## survey the atom and structure and determine if the atom requres
        ## being passed to the naming service

        ## absence of requred fields
        if not atm.chain_id or not atm.fragment_id or not atm.name:
            self.name_service_list.append(atm)
            return atm

        try:
            frag = self.struct[atm.chain_id][atm.fragment_id]
        except KeyError:
            pass
        else:
            if frag.res_name != atm.res_name:
                self.name_service_list.append(atm)
                return atm
                
            for atmx in frag.atom_list:
                if atm.name == atmx.name and \
                   atm.alt_loc == atmx.alt_loc and \
                   atm.model == atmx.model:
                    self.name_service_list.append(atm)
                    return atm

        self.place_atom(atm)
        return atm

    def place_atom(self, atm):
        """Places the atom into the structure, adding the new Chain and
        Fragment if necessary.
        """
        ## pack atom into its fragment, create necessary parents
        ## add chain
        if not self.cache_chain or self.cache_chain.chain_id != atm.chain_id:
            try:
                self.cache_chain = self.struct[atm.chain_id]
            except KeyError:
                self.cache_chain = Chain(atm.chain_id)
                self.struct.add_chain(self.cache_chain, delay_sort = True)

        ## add fragment
        if not self.cache_frag or \
               self.cache_frag.fragment_id != atm.fragment_id or \
               self.cache_frag.chain_id != atm.chain_id:
            try:
                self.cache_frag = self.cache_chain[atm.fragment_id]
            except KeyError:
                if self.struct.library.is_amino_acid(atm.res_name):
                    self.cache_frag = AminoAcidResidue(
                        res_name = atm.res_name,
                        fragment_id = atm.fragment_id,
                        chain_id = atm.chain_id)
                elif self.struct.library.is_nucleic_acid(atm.res_name):
                    self.cache_frag = NucleicAcidResidue(
                        res_name = atm.res_name,
                        fragment_id = atm.fragment_id,
                        chain_id = atm.chain_id)
                else:
                    self.cache_frag = Fragment(
                        res_name = atm.res_name,
                        fragment_id = atm.fragment_id,
                        chain_id = atm.chain_id)
                    
                self.cache_chain.add_fragment(
                    self.cache_frag, delay_sort = True)

        ## sanity check: if this fails, then the name service failed
        assert atm.chain_id == self.cache_chain.chain_id
        assert atm.fragment_id == self.cache_frag.fragment_id
        assert atm.res_name == self.cache_frag.res_name

        ## add atom
        self.cache_frag.add_atom(atm)

    def name_service(self):
        """Runs the name service on all atoms needing to be named.
        """
        if not self.name_service_list:
            return

        ## checks for a matching atom in the list
        def chk_matching_atom(atm, atm_list):
            for atmx in atm_list:
                if atm.name == atmx.name and \
                   atm.model == atmx.model and \
                   atm.alt_loc == atmx.alt_loc:
                    return True
            return False

        ## returns the next available chain_id in self.struct
        ## XXX: it's possible to run out of chain IDs!
        def next_chain_id():
            for chain_id in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz":
                try:
                    self.struct[chain_id]
                except KeyError:
                    return chain_id
                
        chain = []
        chain_list = [chain]

        ## split atoms into chains
        for atm in self.name_service_list:
            try:
                catm = chain[0]
            except IndexError:
                chain.append(atm)
                continue

            if atm.chain_id == catm.chain_id and \
               atm.res_name == catm.res_name and \
               (atm.res_name == "HOH" or not chk_matching_atom(atm, chain)):

                chain.append(atm)

            else:
                chain = [atm]
                chain_list.append(chain)
            
        ## split chains into fragments
        for chain in chain_list:
            chain_id = next_chain_id()

            ## special handling for waters: 1 chain, each water a fragment
            if chain[0].res_name == "HOH":
                for i in range(len(chain)):
                    atm = chain[i]
                    atm.chain_id = chain_id
                    atm.fragment_id = str(i+1)
                    self.place_atom(atm)

            else:
                for atm in chain:
                    atm.chain_id = chain_id
                    atm.fragment_id = "1"
                    self.place_atom(atm)
                
    def read_atoms_finalize(self):
        """After loading all atom records, use the list of atom records to
        build the structure.
        """
        ## name atoms which didn't fit into the Structure hierarch with
        ## their names from the file
        self.name_service()

        ## sort structural objects into their correct order
        self.struct.sort()
        for chain in self.struct.iter_chains():
            chain.sort()

        ## iterate through all the atoms, and choose a default alt_loc
        alt_loc_list = []
        for atm1 in self.struct.iter_atoms():
            for atm2 in atm1:
                if atm2.alt_loc and atm2.alt_loc not in alt_loc_list:
                    alt_loc_list.append(atm2.alt_loc)
        if alt_loc_list:
            alt_loc_list.sort()
            self.struct.default_alt_loc = alt_loc_list[0]
                    
    def read_metadata(self):
        """This method needs to be reimplemented in a fuctional subclass.
        The subclassed read_metadata method should call the various
        load_* methods to set non-atom coordinate data for the Structure.
        """
        pass

    def load_unit_cell(self, ucell_map):
        """Called by the implementation of load_metadata to load the
        unit cell pararameters for the structure.
        """
        for key in ("a", "b", "c", "alpha", "beta", "gamma"):
            if not ucell_map.has_key(key):
                debug("ucell_map missing: %s" % (key))
                return

        if ucell_map.has_key("space_group"):
            self.struct.unit_cell = UnitCell(
                a = ucell_map["a"],
                b = ucell_map["b"],
                c = ucell_map["c"],
                alpha = ucell_map["alpha"],
                beta = ucell_map["beta"],
                gamma = ucell_map["gamma"],
                space_group = ucell_map["space_group"])
        else:
            self.struct.unit_cell = UnitCell(
                a = ucell_map["a"],
                b = ucell_map["b"],
                c = ucell_map["c"],
                alpha = ucell_map["alpha"],
                beta = ucell_map["beta"],
                gamma = ucell_map["gamma"])

    def load_sites(self, site_map):
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
            self.struct.sites.append(site)

            for sentry in sentry_list:
                chain_id = sentry["chain_id"]
                fragment_id = sentry["fragment_id"]
                try:
                    frag = self.struct[chain_id][fragment_id]
                except KeyError:
                    debug("load_site: chain_id=%s fragment_id=%s not found"%(
                        chain_id, fragment_id))
                else:
                    site.append(frag)

    def load_bonds(self, bond_map):
        """Call by the implementation of load_metadata to load bond
        information on the structure.  The keys of the bond map are a 2-tuple
        of the bonded Atom instances, and the value is a dictionary
        containing information on the type of bond, which may also
        be a symmetry operator.

        [bond_map]
        keys: (atm1, atm2)
        values: bond_data_map(s)

        [bond_data_map]
        bond_type -> text description of bond type: covalent, salt bridge,
                     hydrogen, cispeptide
                     
        atm1_symop -> symmetry operation (if any) to be applied to atm1
        atm2_symop -> same as above, for atom 2

        The symmetry operations themselves are a 3x4 array of floating point
        values composed of the 3x3 rotation matrix and the 3x1 translation.
        """
        for ((atm1, atm2), bd_map) in bond_map.items():
            atm1.create_bonds(atm2,
                              bond_type = bd_map.get("bond_type"),
                              atm1_symop = bd_map.get("atm1_symop"),
                              atm2_symop = bd_map.get("atm2_symop"),
                              standard_res_bond = False)

    def read_metadata_finalize(self):
        """Called after the the metadata loading is complete.
        """
        pass
    
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
        ## calculate sequences for all chains
        if "calc_sequence" in self.build_properties:
            for chain in self.struct.iter_chains():
                chain.calc_sequence()

        ## build bonds within structure
        if "bonds" in self.build_properties:
            for frag in self.struct.iter_fragments():
                frag.create_bonds()

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

    def read_atoms(self):
        ## map PDB atom serial numbers to the structure atom classes
        self.atom_serial_map = {}

        ## state vars
        atm_map = {}
        model_num = None

        ## small function to load the atom and add it to the serial map
        def load_atom(atm_map):
            atm = self.load_atom(atm_map)
            try:
                self.atom_serial_map[atm_map["serial"]] = atm
            except KeyError:
                pass

        ## loop over all records
        for rec in self.pdb_file:

            if isinstance(rec, PDB.ATOM):
                if atm_map:
                    load_atom(atm_map)

                atm_map = {}
                self.process_ATOM(atm_map, rec)
                if model_num != None:
                    atm_map["model_num"] = model_num

            elif isinstance(rec, PDB.SIGATM):
                self.process_SIGATM(atm_map, rec)

            elif isinstance(rec, PDB.ANISOU):
                self.process_ANISOU(atm_map, rec)

            elif isinstance(rec, PDB.SIGUIJ):
                self.process_SIGUIJ(atm_map, rec)

            elif isinstance(rec, PDB.MODEL):
                model_num = rec.get("serial")

            elif isinstance(rec, PDB.ENDMDL):
                model_num = None

        ## load last atom read
        if atm_map:
            load_atom(atm_map)

    def read_metadata(self):
        ucell_map = {}
        bond_map = {}

        ## gather metadata
        for rec in self.pdb_file:
            if isinstance(rec, PDB.SITE):
                self.process_SITE(rec)

            elif isinstance(rec, PDB.HELIX):
                self.process_HELIX(rec)

            elif isinstance(rec, PDB.CRYST1):
                self.process_CRYST1(ucell_map, rec)

            elif isinstance(rec, PDB.HEADER):
                self.process_HEADER(rec)

            elif isinstance(rec, PDB.TITLE):
                self.process_TITLE(rec)

            elif isinstance(rec, PDB.COMPND):
                self.process_COMPND(rec)

            elif isinstance(rec, PDB.AUTHOR):
                self.process_AUTHOR(rec)

            elif isinstance(rec, PDB.SSBOND):
                self.process_SSBOND(bond_map, rec)
                
            elif isinstance(rec, PDB.LINK):
                self.process_LINK(bond_map, rec)
                
            elif isinstance(rec, PDB.CONECT):
                self.process_CONECT(bond_map, rec)

        self.load_bonds(bond_map)
        self.load_unit_cell(ucell_map)

    def process_ATOM(self, atm_map, rec):
        self.setmapi(rec, "serial", atm_map, "serial")
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

        ## attempt to guess the element type from the name of atoms if
        ## the element field is not defined in the ATOM record
        try:
            name = atm_map["name"]
        except KeyError:
            return

        if name[0] == " ":
            atm_map["name"] = name[1:]
            for c in name:
                if c in string.letters:
                    element = c
                    break
                
        elif name[0] in string.letters:
            element = name[:2]

        ## if a element field exists, check to see if it matches
        ## the derived element name
        try:
            if atm_map["element"][0] != element[0]:
                atm_map["element"] = element
        except KeyError:
            atm_map["element"] = element

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

    def process_HELIX(self, rec):
        struct_conf = self.struct.cifdb.confirm_table("struct_conf")
        row = mmCIFRow()
        struct_conf.append(row)

        row["conf_type_id"] = "HELX_P"
        self.setmapi(rec, "serNum", row, "id")
        self.setmaps(rec, "helixID", row, "pdbx_PDB_helix_id")
        self.setmaps(rec, "initResName", row, "beg_auth_comp_id")
        self.setmaps(rec, "initChainID", row, "beg_auth_asym_id")
        seq_id = self.get_fragment_id(rec, "initSeqNum", "initICode")
        if seq_id:
            row["beg_auth_seq_id"] = seq_id
        self.setmaps(rec, "initICode", row, "pdbx_beg_PDB_ins_code")

        self.setmaps(rec, "endResName", row, "end_auth_comp_id")
        self.setmaps(rec, "endChainID", row, "end_auth_asym_id")
        seq_id = self.get_fragment_id(rec, "endSeqNum", "endICode")
        if seq_id:
            row["end_auth_seq_id"] = seq_id
        self.setmaps(rec, "endICode", row, "end_beg_PDB_ins_code")

        self.setmaps(rec, "helixClass", row, "pdbx_PDB_helix_class")
        self.setmaps(rec, "comment", row, "details")
        self.setmaps(rec, "length", row, "pdbx_PDB_helix_length")

    def process_SITE(self, rec):
        struct_site_gen = self.struct.cifdb.confirm_table("struct_site_gen")

        for i in range(1, 5):
            chain_key = "chainID%d" % (i)
            res_name = "resName%d" % (i)
            frag_key = "seq%d" % (i)
            icode_key = "icode%d" % (i)

            if not (rec.has_key(chain_key) or rec.has_key(frag_key)):
                break

            row = mmCIFRow()
            struct_site_gen.append(row)

            self.setmaps(rec, "siteID", row, "site_id")
            self.setmaps(rec, res_name, row, "auth_comp_id")
            self.setmaps(rec, chain_key, row, "auth_asym_id")
            seq_id = self.get_fragment_id(rec, frag_key, icode_key)
            if seq_id:
                row["auth_seq_id"] = seq_id
            self.setmaps(rec, icode_key, row, "pdbx_auth_ins_code")

    def process_SSBOND(self, bond_map, rec):
        try:
            chain_id1 = rec["chainID1"]
            frag_id1 = self.get_fragment_id(rec,"seqNum1","iCode1")

            chain_id2 = rec["chainID2"]
            frag_id2 = self.get_fragment_id(rec,"seqNum2","iCode2")
        
            frag1 = self.struct[chain_id1][frag_id1]
            frag2 = self.struct[chain_id2][frag_id2]

            atm1 = frag1["SG"]
            atm2 = frag2["SG"]
        except KeyError, x:
            debug("pdb invalid SSBOND record: %s" % str(x))
            return

        if id(atm1) < id(atm2):
            bnd = (atm1, atm2)
        else:
            bnd = (atm2, atm1)

        if not bond_map.has_key(bnd):
            bond_map[bnd] = {}

        bond_map[bnd]["bond_type"] = "disulf"

        if rec.has_key("sym1"):
            bond_map[bnd]["symop1"] = rec["sym1"]
        if rec.has_key("sym2"):
            bond_map[bnd]["symop2"] = rec["sym2"]

    def process_LINK(self, bond_map, rec):
        try:
            chain_id1 = rec["chainID1"]
            frag_id1 = self.get_fragment_id(rec,"resSeq1","iCode1")
            name1 = rec["name1"]
            alt_loc1 = rec["altLoc1"]

            chain_id2 = rec["chainID2"]
            frag_id2 = self.get_fragment_id(rec,"resSeq2","iCode2")
            name2 = rec["name2"]
            alt_loc2 = rec["altLoc2"]

            atm1 = self.struct[chain_id1][frag_id1][name1]
            atm2 = self.struct[chain_id2][frag_id2][name2]

            if alt_loc1:
                atm1 = atm1[alt_loc1]
            if alt_loc2:
                atm2 = atm2[alt_loc2]
        except KeyError, x:
            debug("pdb invalid LINK record: %s", str(x))
            return

        if id(atm1) < id(atm2):
            bnd = (atm1, atm2)
        else:
            bnd = (atm2, atm1)

        if not bond_map.has_key(bnd):
            bond_map[bnd] = {}

        bond_map[bnd]["bond_type"] = "covale"

        if alt_loc1:
            bond_map[bnd]["alt_loc1"] = alt_loc1
        if alt_loc2:
            bond_map[bnd]["alt_loc2"] = alt_loc2

        if rec.has_key("sym1"):
            bond_map[bnd]["symop1"] = rec["sym1"]
        if rec.has_key("sym2"):
            bond_map[bnd]["symop2"] = rec["sym2"]

    def process_HYDBND(self, bond_map, rec):
        try:
            name = rec["name1"]
            alt_loc = rec.get("altLoc1", "")
            chain_id = rec["chainID1"]
            frag_id = self.get_fragment_id(rec, "resSeq1", "iCode1")
            atm1 = self.struct[chain_id][frag_id][name][alt_loc]
        except KeyError:
            return
        try:
            name = rec["nameH"]
            alt_loc = rec.get("altLocH", "")
            chain_id = rec["chainH"]
            frag_id = self.get_fragment_id(rec, "resSeqH", "iCodeH")
            atmh = self.struct[chain_id][frag_id][name][alt_loc]
        except KeyError:
            atmh = None
        try:
            name = rec["name2"]
            alt_loc = rec.get("altLoc2", "")
            chain_id = rec["chainID2"]
            frag_id = self.get_fragment_id(rec, "resSeq2", "iCode2")
            atm2 = self.struct[chain_id][frag_id][name][alt_loc]
        except KeyError:
            return
        
        if id(atm1) < id(atm2):
            bnd = (atm1, atm2)
        else:
            bnd = (atm2, atm1)

        if not bond_map.has_key(bnd):
            bond_map[bnd] = {}

        bond_map[bnd]["bond_type"] = "hydrog"

        if rec.has_key("sym1"):
            bond_map[bnd]["symop1"] = rec["sym1"]
        if rec.has_key("sym2"):
            bond_map[bnd]["symop2"] = rec["sym2"]

    def process_S(self, bond_map, rec):
        try:
            chain_id1 = rec["chainID1"]
            frag_id1 = self.get_fragment_id(rec,"resSeq1","iCode1")
            name1 = rec["name1"]
            alt_loc1 = rec["altLoc1"]

            chain_id2 = rec["chainID2"]
            frag_id2 = self.get_fragment_id(rec,"resSeq2","iCode2")
            name2 = rec["name2"]
            alt_loc2 = rec["altLoc2"]

            atm1 = self.struct[chain_id1][frag_id1][name1]
            atm2 = self.struct[chain_id2][frag_id2][name2]

            if alt_loc1:
                atm1 = atm1[alt_loc1]
            if alt_loc2:
                atm2 = atm2[alt_loc2]
        except KeyError, x:
            debug("pdb invalid LINK record: %s", str(x))
            return

        if id(atm1) < id(atm2):
            bnd = (atm1, atm2)
        else:
            bnd = (atm2, atm1)

        if not bond_map.has_key(bnd):
            bond_map[bnd] = {}

        bond_map[bnd]["bond_type"] = "saltbr"

        if alt_loc1:
            bond_map[bnd]["alt_loc1"] = alt_loc1
        if alt_loc2:
            bond_map[bnd]["alt_loc2"] = alt_loc2

        if rec.has_key("sym1"):
            bond_map[bnd]["symop1"] = rec["sym1"]
        if rec.has_key("sym2"):
            bond_map[bnd]["symop2"] = rec["sym2"]
        
    def process_CONECT(self, bond_map, rec):
        try:
            atm1 = self.atom_serial_map[rec["serial"]]
        except KeyError:
            debug("pdb CONNECT record atom serial number not found.")
            return

        def helper_func(field_list, bond_type):
            for field in field_list:
                try:
                    atm2 = self.atom_serial_map[rec[field]]
                except KeyError:
                    continue

                if id(atm1) < id(atm2):
                    bnd = (atm1, atm2)
                else:
                    bnd = (atm2, atm1)

                if not bond_map.has_key(bnd):
                    bond_map[bnd] = {}

                if not bond_map[bnd].has_key("bond_type"):
                    bond_map[bnd]["bond_type"] = bond_type

        helper_func(
            ["serialBond1","serialBond2",
             "serialBond3","serialBond4"], "covale")
        helper_func(
            ["serialHydBond1","serialHydBond2",
             "serialHydBond3","serialHydBond4"], "hydrog")
        helper_func(
            ["serialSaltBond1","serialSaltBond2"], "saltbr")

    def process_HEADER(self, rec):
        self.struct.cifdb.set_single(
            "struct_keywords", "pdbx_keywords", rec.get("classification"))
        self.struct.cifdb.set_single(
            "database_pdb_rev", "date_original", rec.get("depDate"))
        self.struct.cifdb.set_single(
            "entry", "id", rec.get("idCode"))

    def process_TITLE(self, rec):
        try:
            title = self.struct.cifdb["struct"]["title"] + rec.get("title")
        except KeyError:
            title = rec.get("title")
        self.struct.cifdb.set_single("struct", "title", title)

    def process_COMPND(self, rec):
        ## multi-record
        return

    def process_AUTHOR(self, rec):
        audit_author = self.struct.cifdb.confirm_table("audit_author")
        for author in rec.get("authorList", "").split(","):
            row = mmCIFRow()
            audit_author.append(row)
            row["name"] = author.strip()

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
        ## optionally use the "auth" mmCIF labels 
        if "auth" in self.build_properties:
            self.atom_id = "auth_atom_id"
            self.alt_id = "auth_alt_id"
            self.comp_id = "auth_comp_id"
            self.seq_id = "auth_seq_id"
            self.asym_id = "auth_asym_id"
            self.ptnr1_atom_id = "ptnr1_auth_atom_id"
            self.ptnr1_comp_id = "ptnr1_auth_comp_id"
            self.ptnr1_asym_id = "ptnr1_auth_asym_id"
            self.ptnr1_seq_id = "ptnr1_auth_seq_id"
            self.ptnr2_atom_id = "ptnr2_auth_atom_id"
            self.ptnr2_comp_id = "ptnr2_auth_comp_id"
            self.ptnr2_asym_id = "ptnr2_auth_asym_id"
            self.ptnr2_seq_id = "ptnr2_auth_seq_id"
        else:
            self.atom_id = "label_atom_id"
            self.alt_id = "label_alt_id"
            self.comp_id = "label_comp_id"
            self.seq_id = "label_seq_id"
            self.asym_id = "label_asym_id"
            self.ptnr1_atom_id = "ptnr1_label_atom_id"
            self.ptnr1_comp_id = "ptnr1_label_comp_id"
            self.ptnr1_asym_id = "ptnr1_label_asym_id"
            self.ptnr1_seq_id = "ptnr1_label_seq_id"
            self.ptnr2_atom_id = "ptnr2_label_atom_id"
            self.ptnr2_comp_id = "ptnr2_label_comp_id"
            self.ptnr2_asym_id = "ptnr2_label_asym_id"
            self.ptnr2_seq_id = "ptnr2_label_seq_id"

        self.atom_site_id_map = {}
            
        self.cif_file = mmCIF.mmCIFFile()
        self.cif_file.load_file(fil)

        ## for a mmCIF file for a structure, assume the first data item
        ## contains the structure
        self.cif_data = self.cif_file[0]

    def read_atoms(self):
        for atom_site in self.cif_data.get("atom_site", []):
            self.read_atom_site(atom_site)

    def read_atom_site(self, atom_site):
        atm_map = {}

        self.setmaps(atom_site, self.atom_id, atm_map, "name")
        self.setmaps(atom_site, self.alt_id, atm_map, "alt_loc")
        self.setmaps(atom_site, self.comp_id, atm_map, "res_name")
        self.setmaps(atom_site, self.seq_id, atm_map, "fragment_id")
        self.setmaps(atom_site, self.asym_id, atm_map, "chain_id")

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
                     atm_map,   "sig_temp_factor")

        self.setmapi(atom_site, "pdbx_PDB_model_num",
                     atm_map,   "model_num")

        if self.cif_data.has_key("atom_site_anisotrop"):
            ctable = self.cif_data["atom_site_anisotrop"]

            aniso = ctable.get_row(("id", atom_site["id"]))
            if aniso:
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

        atm = self.load_atom(atm_map)
        try:
            self.atom_site_id_map[atom_site["id"]] = atm
        except KeyError:
            pass

    def read_metadata(self):
        ## set the mmCIF database to use the strucutre's data lock
        self.struct.cifdb.add_tables(self.cif_data)
        self.read_unit_cell()
        self.read_struct_conn()

    def read_unit_cell(self):
        """Load unit cell and symmetry tables.
        """
        ucell_map = {}
        try:
            cell = self.cif_data["cell"].get_row(
                ("entry_id", self.cif_data["entry"]["id"]))
        except KeyError:
            pass
        else:
            self.setmapf(cell, "length_a", ucell_map, "a")
            self.setmapf(cell, "length_b", ucell_map, "b")
            self.setmapf(cell, "length_c", ucell_map, "c")
            self.setmapf(cell, "angle_alpha", ucell_map, "alpha")
            self.setmapf(cell, "angle_beta", ucell_map, "beta")
            self.setmapf(cell, "angle_gamma", ucell_map, "gamma")
            self.setmapi(cell, "Z_PDB", ucell_map, "z")

        try:
            symm = self.cif_data["symmetry"].get_row(
                ("entry_id",self.cif_data["entry"]["id"]))
        except KeyError:
            pass
        else:
            self.setmaps(symm, "space_group_name_H-M",
                         ucell_map, "space_group")
        
        self.load_unit_cell(ucell_map)

    def read_struct_conn(self):
        """Read bond information form the struct_conn and struct_conn_type
        sections.
        """
        bond_map = {}

        for row in self.cif_data.get("struct_conn", []):

            asym_id1 = row.get(self.ptnr1_asym_id)
            seq_id1 = row.get(self.ptnr1_seq_id)
            comp_id1 = row.get(self.ptnr1_comp_id)
            atom_id1 = row.get(self.ptnr1_atom_id)
            symm1 = row.get("ptnr1_symmetry")

            asym_id2 = row.get(self.ptnr2_asym_id)
            seq_id2 = row.get(self.ptnr2_seq_id)
            comp_id2 = row.get(self.ptnr2_comp_id)
            atom_id2 = row.get(self.ptnr2_atom_id)
            symm2 = row.get("ptnr2_symmetry")

            conn_type = row.get("conn_type_id")

            ## check for these special mmCIF tokens
            if conn_type == "disulf":
                atom_id1 = atom_id2 = "SG"

            as1 = self.cif_data["atom_site"].get_row(
                (self.asym_id, asym_id1),
                (self.seq_id, seq_id1),
                (self.comp_id, comp_id1),
                (self.atom_id, atom_id1))

            as2 = self.cif_data["atom_site"].get_row(
                (self.asym_id, asym_id2),
                (self.seq_id, seq_id2),
                (self.comp_id, comp_id2),
                (self.atom_id, atom_id2))

            if not as1 or not as2:
                debug("read_struct_conn: ignoring: " + conn_type)
                continue

            try:
                atm1 = self.atom_site_id_map[as1["id"]]
                atm2 = self.atom_site_id_map[as2["id"]]
            except KeyError:
                debug("read_struct_conn: atm not found: " + conn_type)
                continue

            if id(atm1) < id(atm2):
                bnd = (atm1, atm2)
            else:
                bnd = (atm2, atm1)

            try:
                bond_map[bnd]["bond_type"] = conn_type
            except KeyError:
                bond_map[bnd] = {"bond_type": conn_type}

            if symm1:
                bond_map[bnd]["symop1"] = symm1
            if symm2:
                bond_map[bnd]["symop2"] = symm2

            print bond_map[bnd]

        ## load the bonds
        self.load_bonds(bond_map)


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

        
