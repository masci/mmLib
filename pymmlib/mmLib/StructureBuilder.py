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


class StructureBuilder(object):
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
        ## create atom object
        atm = Atom(name        = atm_map.get("name", ""),
                   model       = atm_map.get("model_num", 1),
                   alt_loc     = atm_map.get("alt_loc", ""),
                   res_name    = atm_map.get("res_name", ""),
                   fragment_id = atm_map.get("fragment_id", ""),
                   chain_id    = atm_map.get("chain_id", ""))

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
        if atm.chain_id    == "" or \
           atm.fragment_id == "" or \
           atm.name        == "":
            #print "NS1:",atm
            self.name_service_list.append(atm)
            return atm

        try:
            self.place_atom(atm)
        except AtomOverwrite:
            #print "NS2:",atm
            self.name_service_list.append(atm)

        return atm

    def place_atom(self, atm):
        """Places the atom into the structure, adding the new Chain and
        Fragment if necessary.
        """
        assert isinstance(atm, Atom)
        
        ## pack atom into its fragment, create necessary parents
        ## add chain
        if self.cache_chain          == None or \
           self.cache_chain.chain_id != atm.chain_id:

            try:
                self.cache_chain = self.struct[atm.chain_id]
            except KeyError:
                self.cache_chain = Chain(atm.chain_id)
                self.struct.add_chain(self.cache_chain, delay_sort = True)

        ## add fragment
        if self.cache_frag             == None or \
           self.cache_frag.fragment_id != atm.fragment_id or \
           self.cache_frag.chain_id    != atm.chain_id:

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

        self.cache_frag.add_atom(atm)

    def name_service(self):
        """Runs the name service on all atoms needing to be named.
        """
        if len(self.name_service_list) == 0:
            return

        ## returns the next available chain_id in self.struct
        ## XXX: it's possible to run out of chain IDs!
        def next_chain_id(suggest_chain_id):
            if suggest_chain_id != "":
                try:
                    self.struct[suggest_chain_id]
                except KeyError:
                    return suggest_chain_id
            
            for chain_id in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                try:
                    self.struct[chain_id]
                except KeyError:
                    return chain_id

            print "NO MORE CHAIN NAMES"
            sys.exit(1)

        ## cr = (chain_id, res_name)
        ##
        ## cr_dict[cr_key] = model_dict
        ##
        ## model_dict[model] = frag_list
        ##
        ## frag_list = [ frag1, frag2, frag3, ...]
        ##
        ## frag = [atm1, atm2, atm3, ...]

        cr_dict      = {}
        cr_key_list  = []
        
        frag_id   = None
        frag      = None
        name_dict = None

        ## split atoms into fragments
        for atm in self.name_service_list:
            atm_id      = (atm.name, atm.alt_loc)
            atm_frag_id = (atm.model,
                           atm.chain_id,
                           atm.fragment_id,
                           atm.res_name)

            if atm_frag_id == frag_id and not name_dict.has_key(atm_id):
                frag.append(atm)
                name_dict[atm_id] = True

            else:
                cr_key = (atm.chain_id, atm.res_name)
                
                ### debug
                if frag:
                    debug("name_service: fragment detected in cr=%s" % (
                        str(cr_key)))
                    for a in frag:
                        debug("  " + str(a))
                ### /debug
                
                try:
                    model_dict = cr_dict[cr_key]
                except KeyError:
                    model_dict = cr_dict[cr_key] = {}
                    cr_key_list.append(cr_key)

                try:
                    frag_list = model_dict[atm.model]
                except KeyError:
                    frag_list = model_dict[atm.model] = []

                name_dict = {atm_id: True}

                frag_id  = atm_frag_id
                frag     = [atm]
                frag_list.append(frag)

        ## free self.name_service_list and other vars to save some memory
        del self.name_service_list

        for cr_key in cr_key_list:
            ### debug
            debug("name_service: chain_id / res_name keys")
            debug("  cr_key: chain_id='%s' res_name='%s'" % (
                cr_key[0],cr_key[1]))
            ### /debug

            ## get the next chain ID, use the cfr group's
            ## loaded chain_id if possible
            new_chain_id = next_chain_id(cr_key[0])

            ## get model dictionary
            model_dict = cr_dict[cr_key]

            ## inspect the model dictionary to determine the number
            ## of fragments in each model -- they should be the same
            ## and have a 1:1 cooraspondance; if not, match up the
            ## fragments as much as possible
            max_frags = -1
            for (model, frag_list) in model_dict.items():
                frag_list_len = len(frag_list)

                if max_frags == -1:
                    max_frags = frag_list_len
                    continue

                if max_frags != frag_list_len:
                    strx = "name_service: model fragments not identical"
                    debug(strx)
                    warning(strx)
                    max_frags = max(max_frags, frag_list_len)

            ## now iterate through the fragment lists in parallel and assign
            ## the new chain_id and fragment_id
            for i in range(max_frags):
                new_fragment_id = str(i+1)
                
                for frag_list in model_dict.values():
                    try:
                        frag = frag_list[i]
                    except KeyError:
                        continue

                    ## assign new chain_id and fragment_id, than place the
                    ## atom in the structure
                    for atm in frag:
                        atm.chain_id = new_chain_id
                        atm.fragment_id = new_fragment_id
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
        if not "no_bonds" in self.build_properties:
            for frag in self.struct.iter_fragments():
                frag.create_bonds()


class PDBStructureBuilder(StructureBuilder):
    """Builds a new Structure object by loading a PDB file.
    """
    def pdb_error(self, rec_name, text):
        warning("[ERROR PDB::%s] %s" % (rec_name, text))

    def get_fragment_id(self, rec, res_seq = "resSeq", icode = "iCode"):
        fragment_id = None
        if rec.has_key(res_seq):
            fragment_id = str(rec[res_seq])
            if rec.has_key(icode):
                fragment_id += rec[icode]
        return fragment_id
    
    def read_start(self, fil):
        self.pdb_file = PDB.PDBFile()
        self.pdb_file.load_file(fil)

    def load_atom(self, atm_map):
        """Override load_atom to maintain a serial_num->atm map.
        """
        if self.model_num != None:
            atm_map["model_num"] = self.model_num

        atm = StructureBuilder.load_atom(self, atm_map)

        try:
            serial = atm_map["serial"]
        except KeyError:
            pass
        else:
            self.atom_serial_map[serial] = atm

    def read_atoms(self):
        ## map PDB atom serial numbers to the structure atom classes
        self.atom_serial_map = {}
        ## current atom map
        self.atm_map = {}
        ## current model number
        self.model_num = None

        def filter_func(rec):
            if isinstance(rec, PDB.ATOM) or \
               isinstance(rec, PDB.ANISOU) or \
               isinstance(rec, PDB.SIGUIJ) or \
               isinstance(rec, PDB.TER) or \
               isinstance(rec, PDB.MODEL) or \
               isinstance(rec, PDB.ENDMDL):
                return True
            return False

        ## process the coordinate records
        self.pdb_file.record_processor(self, filter_func)

        ## load last atom read
        if self.atm_map:
            self.load_atom(self.atm_map)
        del self.atm_map
        
    def read_metadata(self):
        ## store extracted bond information
        self.bond_map = {}

        def filter_func(rec):
            if isinstance(rec, PDB.ATOM) or \
               isinstance(rec, PDB.ANISOU) or \
               isinstance(rec, PDB.SIGUIJ) or \
               isinstance(rec, PDB.TER) or \
               isinstance(rec, PDB.MODEL) or \
               isinstance(rec, PDB.ENDMDL):
                return False
            return True

        ## process the non-coordinate records
        self.pdb_file.record_processor(self, filter_func)

        ## load chemical bond informaton
        self.load_bonds(self.bond_map)
        del self.bond_map

    def process_ATOM(self, rec):
        ## load current atom since this record indicates a new atom
        if self.atm_map:
            self.load_atom(self.atm_map)
            self.atm_map = {}

        ## name and element
        try:
            name = rec["name"]
        except KeyError:
            self.atm_map["name"] = ""
        else:
            self.atm_map["name"] = name.strip()
            
            element = self.name2element(name)
            if element != None:
                self.atm_map["element"] = element

        ## additional atom information
        setmapi(rec, "serial", self.atm_map, "serial")
        setmaps(rec, "altLoc", self.atm_map, "alt_loc")
        setmaps(rec, "resName", self.atm_map, "res_name")
        setmaps(rec, "chainID", self.atm_map, "chain_id")

        fragment_id = self.get_fragment_id(rec)
        if fragment_id != None:
            self.atm_map["fragment_id"] = fragment_id
            
        setmapf(rec, "x", self.atm_map, "x")
        setmapf(rec, "y", self.atm_map, "y")
        setmapf(rec, "z", self.atm_map, "z")
        setmapf(rec, "occupancy", self.atm_map, "occupancy")
        setmapf(rec, "tempFactor", self.atm_map, "temp_factor")
        setmaps(rec, "charge", self.atm_map, "charge")

    def name2element(self, name0):
        ## set the space_flag to true if the name starts with a space
        ## which can indicate the name of the atom is only 1 charactor
        ## long
        if name0.startswith(" "):
            space_flag = True
        else:
            space_flag = False

        ## strip any space from the name, and return now if there
        ## is nothing left to work with
        name = name0.strip()
        if name == "":
            self.pdb_error("ATOM", "record with blank name field")
            return None

        ## set the digit0_flag to true if the first charactor of the
        ## name is a digit -- commonly used in hydrogen naming
        if name[0] in string.digits:
            digit0_flag = True
        else:
            digit0_flag = False

        ## get the first 1 or 2 consecutive alphabetic charactors
        x = ""
        for c in name:
            if c in string.ascii_letters:
                x += c
                if len(x) == 2:
                    break
            elif x != "":
                break

        ## if no charactors were extracted from the name, then, well,
        ## that's bad
        if x == "":
            self.pdb_error(
                "ATOM",
                "unable to extract element from atom name: %s" % (name0))
            return None

        ## if the space_flag is false, it is possible the element is
        ## two digits, if there is a space, it is likely a one charactor
        ## element name
        if space_flag == True:
            element = x[0]
        else:
            element = x
            
        ## check to see if the element is valid
        if self.struct.library.element_map.has_key(element):
            return element
        elif self.struct.library.element_map.has_key(element[:1]):
            return element[:1]
        else:
            self.pdb_error(
                "ATOM",
                "unable to extract valid element from atom name: %s" % (name0))

        return None

    def process_HETATM(self, rec):
        self.process_ATOM(rec)

    def process_SIGATM(self, rec):
        setmapf(rec, "sigX", self.atm_map, "sig_x")
        setmapf(rec, "sigY", self.atm_map, "sig_y")
        setmapf(rec, "sigZ", self.atm_map, "sig_z")
        setmapf(rec, "sigOccupancy", self.atm_map, "sig_occupancy")
        setmapf(rec, "sigTempFactor", self.atm_map, "sig_temp_factor")

    def process_ANISOU(self, rec):
        self.atm_map["U[1][1]"] = rec.get("u[0][0]", 0.0) / 10000.0
        self.atm_map["U[2][2]"] = rec.get("u[1][1]", 0.0) / 10000.0
        self.atm_map["U[3][3]"] = rec.get("u[2][2]", 0.0) / 10000.0
        self.atm_map["U[1][2]"] = rec.get("u[0][1]", 0.0) / 10000.0
        self.atm_map["U[1][3]"] = rec.get("u[0][2]", 0.0) / 10000.0
        self.atm_map["U[2][3]"] = rec.get("u[1][2]", 0.0) / 10000.0

    def process_SIGUIJ(self, rec):
        self.atm_map["sig_U[1][1]"] = rec.get("sig[1][1]", 0.0) / 10000.0
        self.atm_map["sig_U[2][2]"] = rec.get("sig[2][2]", 0.0) / 10000.0
        self.atm_map["sig_U[3][3]"] = rec.get("sig[3][3]", 0.0) / 10000.0
        self.atm_map["sig_U[1][2]"] = rec.get("sig[1][2]", 0.0) / 10000.0
        self.atm_map["sig_U[1][3]"] = rec.get("sig[1][3]", 0.0) / 10000.0
        self.atm_map["sig_U[2][3]"] = rec.get("sig[2][3]", 0.0) / 10000.0

    def process_MODEL(self, rec):
        self.model_num = rec.get("serial")

    def process_ENDMDL(self, rec):
        self.model_num = None

    def process_HEADER(self, rec):
        self.struct.cifdb.set_single(
            "struct_keywords", "pdbx_keywords", rec.get("classification"))
        self.struct.cifdb.set_single(
            "database_pdb_rev", "date_original", rec.get("depDate"))
        self.struct.cifdb.set_single(
            "entry", "id", rec.get("idCode"))

    def preprocess_TITLE(self, title):
        self.struct.cifdb.set_single("struct", "title", title)

    def preprocess_COMPND(self, compnd_list):
        entity = self.struct.cifdb.confirm_table("entity")
        entity_keywords = self.struct.cifdb.confirm_table("entity_keywords")

        for compnd in compnd_list:
            erow = mmCIFRow()
            ekrow = mmCIFRow()

            setmaps(compnd, "MOLECULE", erow, "pdbx_description")
            if erow:
                entity.append(erow)

            setmaps(compnd, "FRAGMENT", ekrow, "pdbx_fragment")
            setmaps(compnd, "EC", ekrow, "pdbx_ec")
            setmaps(compnd, "MUTATION", ekrow, "pdbx_mutation")
            if ekrow:
                entity_keywords.append(ekrow)

    def preprocess_SOURCE(self, source_list):
        entity_src_nat = self.struct.cifdb.confirm_table("entity_src_nat")
        entity_src_gen = self.struct.cifdb.confirm_table("entity_src_gen")

        for source in source_list:
            nrow = mmCIFRow()
            grow = mmCIFRow()

            setmaps(source, "FRAGMENT",
                    grow, "pdbx_gene_src_fragment")
            setmaps(source, "ORGANISM_SCIENTIFIC",
                    grow, "pdbx_gene_src_scientific_name")
            setmaps(source, "ORGANISM_COMMON",
                    grow, "pdbx_gene_src_common_name")            
            setmaps(source, "GENUS",
                    grow, "pdbx_gene_src_genus")
            setmaps(source, "GENUS",
                    grow, "pdbx_gene_src_genus")
            setmaps(source, "SPECIES",
                    grow, "pdbx_gene_src_species")
            setmaps(source, "STRAIN",
                    grow, "pdbx_gene_src_strain")
            setmaps(source, "VARIANT",
                    grow, "pdbx_gene_src_variant")
            setmaps(source, "CELL_LINE",
                    grow, "pdbx_gene_src_cell_line")
            setmaps(source, "ATCC",
                    grow, "pdbx_gene_src_atcc")
            setmaps(source, "ORGAN",
                    grow, "pdbx_gene_src_organ")
            setmaps(source, "TISSUE",
                    grow, "pdbx_gene_src_tissue")
            setmaps(source, "CELL",
                    grow, "pdbx_gene_src_cell")
            setmaps(source, "ORGANELLE",
                    grow, "pdbx_gene_src_organelle")
            setmaps(source, "SECRETION",
                    nrow, "pdbx_secretion")
            setmaps(source, "CELLULAR_LOCATION",
                    grow, "pdbx_gene_src_cellular_location")
            setmaps(source, "PLASMID",
                    nrow, "pdbx_plasmid_name")
            setmaps(source, "GENE",
                    grow, "pdbx_gene_src_gene")
            setmaps(source, "EXPRESSION_SYSTEM",
                    grow, "pdbx_host_org_scientific_name")
            setmaps(source, "EXPRESSION_SYSTEM_COMMON",
                    grow, "pdbx_host_org_common_name")
            setmaps(source, "EXPRESSION_SYSTEM_GENUS",
                    grow, "pdbx_host_org_genus")
            setmaps(source, "EXPRESSION_SYSTEM_SPECIES",
                    grow, "pdbx_host_org_species")
            setmaps(source, "EXPRESSION_SYSTEM_STRAIN",
                    grow, "pdbx_host_org_strain")
            setmaps(source, "EXPRESSION_SYSTEM_VARIANT",
                    grow, "pdbx_host_org_variant")
            setmaps(source, "EXPRESSION_SYSTEM_CELL_LINE",
                    grow, "pdbx_host_org_cell_line")
            setmaps(source, "EXPRESSION_SYSTEM_ATCC_NUMBER",
                    grow, "pdbx_host_org_atcc")
            setmaps(source, "EXPRESSION_SYSTEM_ORGAN",
                    grow, "pdbx_host_org_organ")
            setmaps(source, "EXPRESSION_SYSTEM_TISSUE",
                    grow, "pdbx_host_org_tissue")
            setmaps(source, "EXPRESSION_SYSTEM_CELL",
                    grow, "pdbx_host_org_cell")
            setmaps(source, "EXPRESSION_SYSTEM_ORGANELLE",
                    grow, "pdbx_host_org_organelle")
            setmaps(source, "EXPRESSION_SYSTEM_CELLULAR_LOCATION",
                    grow, "pdbx_host_org_cellular_location")
            setmaps(source, "EXPRESSION_SYSTEM_VECTOR_TYPE",
                    grow, "pdbx_host_org_vector_type")
            setmaps(source, "EXPRESSION_SYSTEM_VECTOR",
                    grow, "pdbx_host_org_vector")
            setmaps(source, "EXPRESSION_SYSTEM_PLASMID",
                    grow, "plasmid")
            setmaps(source, "EXPRESSION_SYSTEM_GENE",
                    grow, "pdbx_host_org_gene")
            setmaps(source, "OTHER_DETAILS",
                    grow, "pdbx_description")

            if nrow:
                entity_src_nat.append(nrow)
            if grow:
                entity_src_gen.append(grow)

    def preprocess_KEYWDS(self, keywds_list):
        struct_keywords = self.struct.cifdb.confirm_table("struct_keywords")
        for keywds in keywds_list:
            struct_keywords.append(mmCIFRow({"text": keywds}))

    def preprocess_AUTHOR(self, author_list):
        audit_author = self.struct.cifdb.confirm_table("audit_author")
        for author in author_list:
            audit_author.append(mmCIFRow({"name": author}))

    def preprocess_EXPDTA(self, expdta_list):
        exptl = self.struct.cifdb.confirm_table("exptl")
        for (technique, details) in expdta_list:
            row = mmCIFRow({"method": technique})
            if details:
                row["details"] = details
            exptl.append(row)

    def process_CRYST1(self, rec):
        ucell_map = {}

        setmapf(rec, "a", ucell_map, "a")
        setmapf(rec, "b", ucell_map, "b")
        setmapf(rec, "c", ucell_map, "c")
        setmapf(rec, "alpha", ucell_map, "alpha")
        setmapf(rec, "beta", ucell_map, "beta")
        setmapf(rec, "gamma", ucell_map, "gamma")

        setmaps(rec, "sgroup", ucell_map, "space_group")
        setmapi(rec, "z", ucell_map, "z")

        self.load_unit_cell(ucell_map)

    def process_HELIX(self, rec):
        struct_conf = self.struct.cifdb.confirm_table("struct_conf")
        row = mmCIFRow()
        struct_conf.append(row)

        row["conf_type_id"] = "HELX_P"
        setmapi(rec, "serNum", row, "id")
        setmaps(rec, "helixID", row, "pdbx_PDB_helix_id")
        setmaps(rec, "initResName", row, "beg_auth_comp_id")
        setmaps(rec, "initChainID", row, "beg_auth_asym_id")
        seq_id = self.get_fragment_id(rec, "initSeqNum", "initICode")
        if seq_id:
            row["beg_auth_seq_id"] = seq_id
        setmaps(rec, "initICode", row, "pdbx_beg_PDB_ins_code")

        setmaps(rec, "endResName", row, "end_auth_comp_id")
        setmaps(rec, "endChainID", row, "end_auth_asym_id")
        seq_id = self.get_fragment_id(rec, "endSeqNum", "endICode")
        if seq_id:
            row["end_auth_seq_id"] = seq_id
        setmaps(rec, "endICode", row, "end_beg_PDB_ins_code")

        setmaps(rec, "helixClass", row, "pdbx_PDB_helix_class")
        setmaps(rec, "comment", row, "details")
        setmaps(rec, "length", row, "pdbx_PDB_helix_length")

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

            setmaps(rec, "siteID", row, "site_id")
            setmaps(rec, res_name, row, "auth_comp_id")
            setmaps(rec, chain_key, row, "auth_asym_id")
            seq_id = self.get_fragment_id(rec, frag_key, icode_key)
            if seq_id:
                row["auth_seq_id"] = seq_id
            setmaps(rec, icode_key, row, "pdbx_auth_ins_code")

    def bond_processor(self, **args):
        """Complicated method.  Required arguments are:
            rec = PDB record
            atm1/2 = Atom object, if you want to override the lookup
            chain_id_field1/2: PDB field name for the chain ID
            res_seq1/2_field: PDB field for the residue sequence num
            icode1/2_field: PDB field for the residue insertion code
            name1/2_field: PDB field for the atom name
            atl_loc1/2: PDB filed name for the atom alt_loc
            symop1/2_field: PDB field name for the atom symmetry operation

            chain_id1/2: override the chain ID
            frag_id1/2: override the fragmetn ID
            name1/2: override the atom name
            alt_loc1/2: override the atom alt_loc
        """
        rec = args["rec"]

        def get_atom(chain_id, frag_id, name, alt_loc):
            try:
                atm = self.struct[chain_id][frag_id][name]
            except KeyError:
                return None
            except TypeError:
                return None
            
            if alt_loc:
                try:
                    atm = atm[alt_loc]
                except KeyError:
                    pass

            return atm

        ## get atm1
        try:
            atm1 = args["atm1"]
        except KeyError:
            chain_id1 = args.get("chain_id1") or \
                        rec.get(args["chain_id1_field"])

            frag_id1 = args.get("frag_id1") or \
                       self.get_fragment_id(rec,
                           args["res_seq1_field"],args["icode1_field"])

            name1 = args.get("name1") or \
                    rec.get("name1_field")

            alt_loc1 = args.get("alt_loc1") or \
                       rec.get(args["alt_loc1_field"])

            atm1 = get_atom(chain_id1, frag_id1, name1, alt_loc1)

        ## get atm2
        try:
            atm2 = args["atm2"]
        except KeyError:
            chain_id2 = args.get("chain_id2") or \
                        rec.get(args["chain_id2_field"])

            frag_id2 = args.get("frag_id2") or \
                       self.get_fragment_id(rec,
                           args["res_seq2_field"],args["icode2_field"])

            name2 = args.get("name2") or \
                    rec.get("name2_field")

            alt_loc2 = args.get("alt_loc2") or \
                       rec.get(args["alt_loc2_field"])

            atm2 = get_atom(chain_id2, frag_id2, name2, alt_loc2)

        ## unable to retrieve the atoms?
        if not (atm1 and atm2):
            return None

        ## the bond map is keyed from the 2-tuple of the atoms involved in
        ## the bond; they are sorted by their object ID just to have
        ## a definate order
        if id(atm1) < id(atm2):
            bkey = (atm1, atm2)
        else:
            bkey = (atm2, atm1)

        try:
            bond = self.bond_map[bkey]
        except KeyError:
            bond = self.bond_map[bkey] = {}

        ## set bond type
        bond["bond_type"] = args["bond_type"]

        ## symmetry operations
        symop1 = args.get("symop1") or \
                 rec.get(args["symop1_field"])
        symop2 = args.get("symop2") or \
                 rec.get(args["symop2_field"])

        if symop1:
            bond["symop1"] = symop1
        if symop2:
            bond["symop2"] = symop2

        return bkey

    def process_SSBOND(self, rec):
        x = self.bond_processor(
            rec = rec,
            
            bond_type = "disulf",
            
            chain_id1_field = "chainID1",
            res_seq1_field = "seqNum1",
            icode1_field = "iCode1",
            name1 = "SG",
            alt_loc1_field = None,
            symop1_field = "sym1",
            
            chain_id2_field = "chainID2",
            res_seq2_field = "seqNum2",
            icode2_field = "iCode2",
            name2 = "SG",
            alt_loc2_field = None,
            symop2_field = "sym2")

        if not x:
            self.pdb_error("SSBOND", "Atom not found")

    def process_LINK(self, rec):
        x = self.bond_processor(
            rec = rec,
            
            bond_type = "covale",
            
            chain_id1_field = "chainID1",
            res_seq1_field = "resSeq1",
            icode1_field = "iCode1",
            name1_field = "name1",
            alt_loc1_field = "altLoc1",
            symop1_field = "sym1",
            
            chain_id2_field = "chainID2",
            res_seq2_field = "resSeq2",
            icode2_field = "iCode2",
            name2_field = "name2",
            alt_loc2_field = "altLoc2",
            symop2_field = "sym2")

        if not x:
            self.pdb_error("LINK", "Atom not found")

    def process_HYDBND(self, rec):
        ## retrieve the hydrogen atom
        try:
            name = rec["nameH"]
            alt_loc = rec.get("altLocH", "")
            chain_id = rec["chainH"]
            frag_id = self.get_fragment_id(rec, "resSeqH", "iCodeH")
            atmh = self.struct[chain_id][frag_id][name][alt_loc]
        except KeyError:
            atmh = None
        
        x = self.bond_processor(
            rec = rec,
            
            bond_type = "hydrog",
            
            chain_id1_field = "chainID1",
            res_seq1_field = "resSeq1",
            icode1_field = "iCode1",
            name1_field = "name1",
            alt_loc1_field = "altLoc1",
            symop1_field = "sym1",
            
            chain_id2_field = "chainID2",
            res_seq2_field = "resSeq2",
            icode2_field = "iCode2",
            name2_field = "name2",
            alt_loc2_field = "altLoc2",
            symop2_field = "sym2")
        
        if not x:
            self.pdb_error("HYDBND", "Atom not found")

    def process_SLTBRG(self, rec):
        x = self.bond_processor(
            rec = rec,
            
            bond_type = "saltbr",
            
            chain_id1_field = "chainID1",
            res_seq1_field = "resSeq1",
            icode1_field = "iCode1",
            name1_field = "name1",
            alt_loc1_field = "altLoc1",
            symop1_field = "sym1",
            
            chain_id2_field = "chainID2",
            res_seq2_field = "resSeq2",
            icode2_field = "iCode2",
            name2_field = "name2",
            alt_loc2_field = "altLoc2",
            symop2_field = "sym2")

        if not x:
            self.pdb_error("SLTBRG", "Atom not found")
        
    def process_CONECT(self, rec):
        try:
            serial = rec["serial"]
        except KeyError:
            self.pdb_error("CONECT", "missing serial field")
            return

        try:
            atm1 = self.atom_serial_map[serial]
        except KeyError:
            self.pdb_error("CONECT", "incorrect serial number")
            return

        def helper_func(field_list, bond_type):
            for field in field_list:
                try:
                    serial2 = rec[field]
                except KeyError:
                    continue

                try:
                    atm2 = self.atom_serial_map[serial2]
                except KeyError:
                    self.pdb_error("CONECT", "incorrect serial number")
                    continue

                x = self.bond_processor(
                    rec = rec,
                    bond_type = bond_type,
                    atm1 = atm1,
                    atm2 = atm2,
                    symop1_field = None,
                    symop2_field = None)

        helper_func(
            ["serialBond1","serialBond2",
             "serialBond3","serialBond4"], "covale")
        helper_func(
            ["serialHydBond1","serialHydBond2",
             "serialHydBond3","serialHydBond4"], "hydrog")
        helper_func(
            ["serialSaltBond1","serialSaltBond2"], "saltbr")


class mmCIFStructureBuilder(StructureBuilder):
    """Builds a new Structure object by loading a mmCIF file.
    """
    def setmaps(self, smap, skey, dmap, dkey, default = None):
        """For string converisons, treat question marks as blank.
        """
        ## get data from smap
        x = None
        try:
            x = smap[skey]
        except KeyError:
            pass

        ## set default under conditions
        if x == None or x == "?":
            if default == None:
                return False
            else:
                dmap[dkey] = default
                return True
            
        ## set dmap
        try:
            dmap[dkey] = str(x)
        except ValueError:
            return False

        return True

    def read_start(self, fil):
        ## optionally use the "auth" mmCIF labels 
        if not "label" in self.build_properties:
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

        ## parse the mmCIF file
        self.cif_file = mmCIF.mmCIFFile()
        self.cif_file.load_file(fil)

        ## for a mmCIF file for a structure, assume the first data item
        ## contains the structure
        self.cif_data = self.cif_file[0]

        ## maintain a map of atom_site.id -> atm
        self.atom_site_id_map = {}
        

    def read_atoms(self):
        try:
            atom_site_table = self.cif_data["atom_site"]
        except KeyError:
            warning("read_atoms: atom_site table not found")
            return

        try:
            aniso_table = self.cif_data["atom_site_anisotrop"]
        except KeyError:
            aniso_table = None
        else:
            aniso_dict  = aniso_table.row_index_dict("id")
        
        for atom_site in atom_site_table:
            try:
                atom_site_id = atom_site["id"]
            except KeyError:
                warning("unable to find id for atom_site row")
                continue

            atm_map = {}

            self.setmaps(atom_site, self.atom_id, atm_map, "name")
            self.setmaps(atom_site, self.alt_id, atm_map, "alt_loc")
            self.setmaps(atom_site, self.comp_id, atm_map, "res_name")
            self.setmaps(atom_site, self.seq_id, atm_map, "fragment_id")
            self.setmaps(atom_site, self.asym_id, atm_map, "chain_id")

            self.setmaps(atom_site, "type_symbol", atm_map, "element")
            setmapf(atom_site, "Cartn_x", atm_map, "x")
            setmapf(atom_site, "Cartn_y", atm_map, "y")
            setmapf(atom_site, "Cartn_z", atm_map, "z")
            setmapf(atom_site, "occupancy", atm_map, "occupancy")
            setmapf(atom_site, "B_iso_or_equiv", atm_map, "temp_factor")
            setmapf(atom_site, "Cartn_x_esd", atm_map, "sig_x")
            setmapf(atom_site, "Cartn_y_esd", atm_map, "sig_y")
            setmapf(atom_site, "Cartn_z_esd", atm_map, "sig_z")
            setmapf(atom_site, "occupancy_esd", atm_map, "sig_occupancy")

            setmapf(atom_site, "B_iso_or_equiv_esd",
                         atm_map,   "sig_temp_factor")

            setmapi(atom_site, "pdbx_PDB_model_num",
                         atm_map,   "model_num")

            if aniso_table != None:
                try:
                    aniso = aniso_dict[atom_site_id]
                except KeyError:
                    warning("unable to find aniso row for atom")
                else:
                    setmapf(aniso, "U[1][1]", atm_map, "U[1][1]")
                    setmapf(aniso, "U[2][2]", atm_map, "U[2][2]")
                    setmapf(aniso, "U[3][3]", atm_map, "U[3][3]")
                    setmapf(aniso, "U[1][2]", atm_map, "U[1][2]")
                    setmapf(aniso, "U[1][3]", atm_map, "U[1][3]")
                    setmapf(aniso, "U[2][3]", atm_map, "U[2][3]")

                    setmapf(aniso, "U[1][1]_esd", atm_map, "sig_U[1][1]")
                    setmapf(aniso, "U[2][2]_esd", atm_map, "sig_U[2][2]")
                    setmapf(aniso, "U[3][3]_esd", atm_map, "sig_U[3][3]")
                    setmapf(aniso, "U[1][2]_esd", atm_map, "sig_U[1][2]")
                    setmapf(aniso, "U[1][3]_esd", atm_map, "sig_U[1][3]")
                    setmapf(aniso, "U[2][3]_esd", atm_map, "sig_U[2][3]")

            atm = self.load_atom(atm_map)
            self.atom_site_id_map[atom_site_id] = atm

    def read_metadata(self):
        ## copy selected mmCIF tables to the structure's mmCIF database
        skip_tables = ["atom_site",
                       "atom_site_anisotrop",
                       "atom_sites_alt"]
        
        for table in self.cif_data:
            if table.name not in skip_tables:
                self.struct.cifdb.add_table(table)

        ## read unit cell table
        self.read_unit_cell()
        ## read bond information
        self.read_struct_conn()

    def read_unit_cell(self):
        """Load unit cell and symmetry tables.
        """
        ucell_map = {}
        
        try:
            entry_id = self.cif_data["entry"]["id"]
        except KeyError:
            warning("read_unit_cell: entry id not found")
            return

        try:
            cell_table = self.cif_data["cell"]
        except KeyError:
            warning("read_unit_cell: cell table not found")
        else:
            cell = cell_table.get_row(("entry_id", entry_id))
            if cell != None:
                setmapf(cell, "length_a", ucell_map, "a")
                setmapf(cell, "length_b", ucell_map, "b")
                setmapf(cell, "length_c", ucell_map, "c")
                setmapf(cell, "angle_alpha", ucell_map, "alpha")
                setmapf(cell, "angle_beta", ucell_map, "beta")
                setmapf(cell, "angle_gamma", ucell_map, "gamma")
                setmapi(cell, "Z_PDB", ucell_map, "z")

        try:
            symmetry_table = self.cif_data["symmetry"]
        except KeyError:
            warning("read_unit_cell: symmetry table not found")
        else:
            symm = symmetry_table.get_row(("entry_id", entry_id))
            if symm != None:
                self.setmaps(symm, "space_group_name_H-M",
                             ucell_map, "space_group")
        
        self.load_unit_cell(ucell_map)

    def read_struct_conn(self):
        """Read bond information form the struct_conn and struct_conn_type
        sections.
        """
        ## only read these types of bonds for now
        bond_type_list = [
            "disulf", "covale",
            ]

        try:
            atom_site = self.cif_data["atom_site"]
        except KeyError:
            warning("read_struct_conn: atom_site table not found")
            return

        try:
            struct_conn_table = self.cif_data["struct_conn"]
        except KeyError:
            warning("read_struct_conn: struct_conn table not found")
            return

        bond_map = {}

        for row in struct_conn_table:
            conn_type = row.get("conn_type_id")
            if conn_type not in bond_type_list:
                continue
            
            asym_id1 = row.get(self.ptnr1_asym_id)
            seq_id1  = row.get(self.ptnr1_seq_id)
            comp_id1 = row.get(self.ptnr1_comp_id)
            atom_id1 = row.get(self.ptnr1_atom_id)
            symm1    = row.get("ptnr1_symmetry")

            asym_id2 = row.get(self.ptnr2_asym_id)
            seq_id2  = row.get(self.ptnr2_seq_id)
            comp_id2 = row.get(self.ptnr2_comp_id)
            atom_id2 = row.get(self.ptnr2_atom_id)
            symm2    = row.get("ptnr2_symmetry")

            ## check for these special mmCIF tokens
            if conn_type == "disulf":
                atom_id1 = atom_id2 = "SG"

            as1 = atom_site.get_row(
                (self.asym_id, asym_id1),
                (self.seq_id,  seq_id1),
                (self.comp_id, comp_id1),
                (self.atom_id, atom_id1))

            as2 = atom_site.get_row(
                (self.asym_id, asym_id2),
                (self.seq_id,  seq_id2),
                (self.comp_id, comp_id2),
                (self.atom_id, atom_id2))

            if not as1 or not as2:
                warning("read_struct_conn: atom not found id: " + \
                        row.get("id","[No ID]"))
                
                warning("atm1: asym=%s seq=%s comp=%s atom=%s symm=%s" % (
                    asym_id1, seq_id1, comp_id1, atom_id1, symm1))
                
                warning("atm2: asym=%s seq=%s comp=%s atom=%s symm=%s" % (
                    asym_id2, seq_id2, comp_id2, atom_id2, symm2))

                continue

            try:
                atm1 = self.atom_site_id_map[as1["id"]]
                atm2 = self.atom_site_id_map[as2["id"]]
            except KeyError:
                warning("read_struct_conn: atom_site_id_map incorrect id: " + \
                        row.get("id", "[No ID]"))

                warning("atm1: asym=%s seq=%s comp=%s atom=%s symm=%s" % (
                    asym_id1, seq_id1, comp_id1, atom_id1, symm1))
                
                warning("atm2: asym=%s seq=%s comp=%s atom=%s symm=%s" % (
                    asym_id2, seq_id2, comp_id2, atom_id2, symm2))

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

            print "Bond[%s,%s] = %s" % (atm1, atm2, conn_type)


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

        
