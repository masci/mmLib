## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Convert a Structure object into its mmCIFFile description.
"""
from __future__ import generators
import copy


from mmTypes          import *
from mmCIF            import *
from Structure        import *
from StructureBuilder import *


mmCIFStandardColumnsMap = {
    "entry": ["id"],

    "entity": ["id",
               "type",
               "details"],

    "cell": ["entry_id",
             "length_a",
             "length_b",
             "length_c",
             "angle_alpha",
             "angle_beta",
             "angle_gamma",
             "PDB_Z"],

    "symmetry": ["entry_id",
                 "space_group_name_H-M",
                 "cell_setting",
                 "Int_Tables_number"],

    "entity_poly_seq": ["entity_id",
                        "num",
                        "mon_id"],

    "entity_poly": ["entity_id",
                    "type",
                    "nstd_linkage",
                    "nstd_monomer",
                    "pdbx_seq_one_letter_code"],

    "audit_author": ["name"],

    "atom_site": ["group_PDB",
                  "id",
                  "type_symbol",
                  "label_entity_id",
                  "label_asym_id",
                  "label_seq_id",
                  "label_comp_id",
                  "label_alt_id",
                  "label_atom_id",
                  "Cartn_x",
                  "Cartn_y",
                  "Cartn_z", 
                  "occupancy",
                  "B_iso_or_equiv",
                  "Cartn_x_esd",
                  "Cartn_y_esd",
                  "Cartn_z_esd",
                  "occupancy_esd",
                  "B_iso_or_equiv_esd",
                  "auth_asym_id",
                  "auth_seq_id",
                  "auth_comp_id",
                  "auth_alt_id",
                  "auth_atom_id",
                  "pdbx_PDB_model_num"],

    "atom_site_anisotrop": ["id",
                            "type_symbol",
                            "label_entity_id",
                            "U[1][1]",
                            "U[1][2]",
                            "U[1][3]",
                            "U[2][2]",
                            "U[2][3]",
                            "U[3][3]",
                            "U[1][1]_esd",
                            "U[1][2]_esd",
                            "U[1][3]_esd",
                            "U[2][2]_esd",
                            "U[2][3]_esd",
                            "U[3][3]_esd",
                            "pdbx_auth_seq_id",
                            "pdbx_auth_comp_id",
                            "pdbx_auth_asym_id",
                            "pdbx_auth_atom_id"]
    }







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

    def read_start(self, fil, update_cb = None):
        ## parse the mmCIF file
        self.cif_file = mmCIFFile()
        self.cif_file.load_file(fil, update_cb)

        ## for a mmCIF file for a structure, assume the first data item
        ## contains the structure
        self.cif_data = self.cif_file[0]
        self.set_atom_site_data_columns()

        ## maintain a map of atom_site.id -> atm
        self.atom_site_id_map = {}

    def set_atom_site_auth(self):
        """Read atom_site.auth_ labels for atom definitions.
        """
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
        
    def set_atom_site_label(self):
        """Read atom_site.label_ items for atom definitions.
        """
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

    def set_atom_site_data_columns(self):
        """Choose to use atom_site.auth_ labels, or atom_site.label_
        """
        try:
            atom_site_table = self.cif_data["atom_site"]
        except KeyError:
            return

        ## count the number of columns which exist for the auth_ style
        ## columns and label_ style columns
        auth_cols  = ["auth_atom_id", "auth_comp_id", "auth_seq_id",
                      "auth_asym_id"]

        label_cols = ["label_atom_id", "label_comp_id", "label_seq_id",
                      "label_asym_id"]
        
        auth_count = 0
        for col in auth_cols:
            if col in atom_site_table.columns:
                auth_count += 1
        
        label_count = 0
        for col in label_cols:
            if col in atom_site_table.columns:
                label_count += 1
        
        if auth_count>=label_count:
            self.set_atom_site_auth()
        else:
            self.set_atom_site_label()

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
            self.setmaps(atom_site, self.alt_id,  atm_map, "alt_loc")
            self.setmaps(atom_site, self.comp_id, atm_map, "res_name")
            self.setmaps(atom_site, self.seq_id,  atm_map, "fragment_id")
            self.setmaps(atom_site, self.asym_id, atm_map, "chain_id")

            self.setmaps(atom_site, "type_symbol", atm_map, "element")
            setmapf(atom_site, "cartn_x", atm_map, "x")
            setmapf(atom_site, "cartn_y", atm_map, "y")
            setmapf(atom_site, "cartn_z", atm_map, "z")
            setmapf(atom_site, "occupancy", atm_map, "occupancy")
            setmapf(atom_site, "b_iso_or_equiv", atm_map, "temp_factor")
            setmapf(atom_site, "cartn_x_esd", atm_map, "sig_x")
            setmapf(atom_site, "cartn_y_esd", atm_map, "sig_y")
            setmapf(atom_site, "cartn_z_esd", atm_map, "sig_z")
            setmapf(atom_site, "occupancy_esd", atm_map, "sig_occupancy")

            setmapf(atom_site, "b_iso_or_equiv_esd",
                    atm_map,   "sig_temp_factor")

            setmapi(atom_site, "pdbx_pdb_model_num",
                    atm_map,   "model_id")

            if aniso_table != None:
                try:
                    aniso = aniso_dict[atom_site_id]
                except KeyError:
                    warning("unable to find aniso row for atom")
                else:
                    setmapf(aniso, "u[1][1]", atm_map, "u11")
                    setmapf(aniso, "u[2][2]", atm_map, "u22")
                    setmapf(aniso, "u[3][3]", atm_map, "u33")
                    setmapf(aniso, "u[1][2]", atm_map, "u12")
                    setmapf(aniso, "u[1][3]", atm_map, "u13")
                    setmapf(aniso, "u[2][3]", atm_map, "u23")

                    setmapf(aniso, "u[1][1]_esd", atm_map, "sig_u12")
                    setmapf(aniso, "u[2][2]_esd", atm_map, "sig_u22")
                    setmapf(aniso, "u[3][3]_esd", atm_map, "sig_u33")
                    setmapf(aniso, "u[1][2]_esd", atm_map, "sig_u12")
                    setmapf(aniso, "u[1][3]_esd", atm_map, "sig_u13")
                    setmapf(aniso, "u[2][3]_esd", atm_map, "sig_u23")

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
                setmapi(cell, "z_pdb", ucell_map, "z")

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

        ## load the bonds
        self.load_bonds(bond_map)

    def read_struct_conf(self):
        """Reads the struct_conf table getting information on alpha
        helicies and turns in the structure.
        """
        try:
            struct_conf = self.cif_data["struct_conf"]
        except KeyError:
            return

        ## iterate over struct_conf and create the helix_list
        helix_list = []
        
        for row in struct_conf:
            ## check for required fields
            try:
                row["id"]
                row["conf_type_id"]
            except KeyError:
                continue

            ## if this is a alpha helix
            if row["conf_type_id"].startswith("HELIX"):
                helix = {"helix_id":    row["id"],
                         "helix_class": row["conf_type_id"]}

                           


class mmCIFFileBuilder(object):
    """Builds a mmCIF file from a Structure object.
    """
    def __init__(self, struct, cif_file):
        self.struct = struct

        try:
            self.entry_id = self.struct.cifdb["entry"]["id"]
        except KeyError:
            self.entry_id = "XXXX"

        self.cif_data = mmCIFData(self.entry_id)
        cif_file.append(self.cif_data)

        ## maps Fragment -> entity_id number
        self.entity_id_map = {}

        ## tables which are not generated from the structure hierarchy
        ## can be copied directly from the structure's cif database
        for table in self.struct.cifdb:
            if not mmCIFStandardColumnsMap.has_key(table.name):
                new_table = copy.deepcopy(table)
                new_table.autoset_columns()
                self.cif_data.append(new_table)

        ## these tables need to be formed from the atom structure
        self.add__entry()
        self.add__entity()
        self.add__cell()
        self.add__symmetry()
        self.add__atom_site()

    def get_table(self, name):
        """Returns the self.cif_data[name] mmCIFTable, or it creates
        it and adds it to self.cif_data if it does not exist.
        """
        try:
            return self.cif_data[name]
        except KeyError:
            pass

        columns = copy.deepcopy(mmCIFStandardColumnsMap[name])
        
        table = mmCIFTable(name, columns)
        self.cif_data.append(table)
        return table

    def add__entry(self):
        """Add the _entry table.
        """
        entry = self.get_table("entry")
        row = mmCIFRow()
        entry.append(row)
        row["id"] = self.entry_id

    def add__entity(self):
        """Adds the entity table.  The entity names are faked here, since
        it's really not clear to me how the names are chosen by the PDB.
        """

        def iter_all_chains():
            for model in self.struct.iter_models():
                for chain in model.iter_chains():
                    yield chain

        ## maps fragment -> entity_id
        entity = self.get_table("entity")

        ## ADD BIO-POLYMERS
        ## list of polymer entities (entity_id, sequence1)
        poly_entity_list = []
        for chain in iter_all_chains():

            ## if the chain is a bio-polymer, it is one entity; come up
            ## with a name from its sequence and add it to the
            ## entity map
            if not chain.has_standard_residues():
                continue

            ## calculate sequence and compare the sequence to chains
            ## already added so we can re-use the entity ID
            entity_id = None

            if chain.sequence == None:
                sequence = chain.calc_sequence()
            else:
                sequence = chain.sequence

            for (eid, seq1) in poly_entity_list:
                if seq1 == sequence:
                    entity_id = eid
                    break

            ## new entity!
            if entity_id == None:

                ## figure out what type of biopolymer this is
                if self.struct.library.is_amino_acid(sequence[0]):
                    details   = "%d residue polypeptide" % (len(sequence))
                    poly_type = "polypeptide(L)"

                elif self.struct.library.is_nucleic_acid(sequence[0]):
                    details   = "%d residue DNA/RNA" % (len(sequence))
                    poly_type = "polydeoxyribonucleotide"
                    
                else:
                    details   = "unknown polymer"
                    poly_type = "other"
                
                row = mmCIFRow()
                entity.append(row)

                entity_id = entity.index(row) + 1
                poly_entity_list.append((entity_id, sequence))
                
                row["id"]      = entity_id
                row["type"]    = "polymer"
                row["details"] = details

##                 ## add the new sequence to the entity_poly_seq table
##                 entity_poly_seq = self.get_table("entity_poly_seq")
##                 for i in range(len(sequence)):
##                     row = mmCIFRow()
##                     entity_poly_seq.append(row)

##                     row["entity_id"] = entity_id
##                     row["num"]       = i + 1
##                     row["mon_id"]    = sequence[i]

                ## add the new sequence to the entity_poly table
                entity_poly = self.get_table("entity_poly")
                row = mmCIFRow()
                entity_poly.append(row)
                row["entity_id"]                = entity_id                
                row["type"]                     = poly_type
                row["pdbx_seq_one_letter_code"] = \
                    sequence.sequence_one_letter_code()


            ## loop over all residues and map the Residues to their entity_id
            for res in chain.iter_standard_residues():
                self.entity_id_map[res] = entity_id


        ## ADD HET ATOMS (Water, metal)
        er_map = {}

        for chain in iter_all_chains():
            for frag in chain.iter_non_standard_residues():

                ## already assigned a entity_id for this fragment_id
                if er_map.has_key(frag.res_name):
                    self.entity_id_map[frag] = er_map[frag.res_name]

                ## we need to assign a entity_id for this fragment_id
                ## and add a row for it in the entity table
                else:
                    row = mmCIFRow()
                    entity.append(row)

                    entity_id = row["id"] = entity.index(row) + 1

                    if frag.is_water():
                        row["type"] = "water"
                        row["details"] = ""
                    else:
                        row["type"] = "non-polymer"
                        row["details"] = frag.res_name

                    er_map[frag.res_name] = entity_id
                    self.entity_id_map[frag] = entity_id

    def add__cell(self):
        """Adds the _cell table.
        """
        try:
            unit_cell = self.struct.unit_cell
        except AttributeError:
            return
        
        cell = self.get_table("cell")
        row = mmCIFRow()
        cell.append(row)

        row["entry_id"]    = self.cif_data["entry"]["id"]
        row["length_a"]    = unit_cell.a
        row["length_b"]    = unit_cell.b
        row["length_c"]    = unit_cell.c
        row["angle_alpha"] = unit_cell.calc_alpha_deg()
        row["angle_beta"]  = unit_cell.calc_beta_deg()
        row["angle_gamma"] = unit_cell.calc_gamma_deg()

    def add__symmetry(self):
        """Adds the _symmetry table.
        """
        try: 
            space_group = self.struct.unit_cell.space_group
        except AttributeError:
            return

        cell = self.get_table("symmetry")
        row = mmCIFRow()
        cell.append(row)

        row["entry_id"]             = self.cif_data["entry"]["id"]
        row["space_group_name_H-M"] = space_group.pdb_name
        row["Int_Tables_number"]    = space_group.number

    def add__atom_site(self):
        """Adds the _atom_site table.
        """
        atom_site = self.get_table("atom_site")        
        atom_id = 0

        for atm in self.struct.iter_all_atoms():
            asrow = mmCIFRow()
            atom_site.append(asrow)

            atom_id += 1
            asrow["id"] = atom_id
            self.set_atom_site_row(asrow, atm)

    def set_atom_site_row(self, asrow, atm):
        if atm.get_fragment().is_standard_residue():
            asrow["group_PDB"] = "ATOM"
        else:
            asrow["group_PDB"] = "HETATM"

        asrow["label_entity_id"]    = self.entity_id_map[atm.get_fragment()]
        asrow["label_atom_id"]      = atm.name
        asrow["label_alt_id"]       = atm.alt_loc
        asrow["label_comp_id"]      = atm.res_name
        asrow["label_seq_id"]       = atm.fragment_id
        asrow["label_asym_id"]      = atm.chain_id

        asrow["auth_atom_id"]       = atm.name
        asrow["auth_alt_id"]        = atm.alt_loc
        asrow["auth_comp_id"]       = atm.res_name
        asrow["auth_seq_id"]        = atm.fragment_id
        asrow["auth_asym_id"]       = atm.chain_id

        asrow["type_symbol"]        = atm.element
        asrow["Cartn_x"]            = atm.position[0]
        asrow["Cartn_y"]            = atm.position[1]
        asrow["Cartn_z"]            = atm.position[2]
        asrow["occupancy"]          = atm.occupancy
        asrow["B_iso_or_equiv"]     = atm.temp_factor
        asrow["pdbx_PDB_model_num"] = atm.model_id
        
        if atm.sig_position:
            asrow["Cartn_x_esd"] = atm.sig_position[0]
            asrow["Cartn_y_esd"] = atm.sig_position[1]
            asrow["Cartn_z_esd"] = atm.sig_position[2]
            asrow["occupancy_esd"] = atm.sig_occupancy
            asrow["B_iso_or_equiv_esd"] = atm.sig_temp_factor

        if atm.U:
            aniso = self.get_table("atom_site_anisotrop")
            anrow = mmCIFRow()
            aniso.append(anrow)

            anrow["id"] = asrow["id"]
            anrow["type_symbol"] = asrow["type_symbol"]
            anrow["label_entity_id"] = asrow["label_entity_id"]
            anrow["pdbx_auth_seq_id"] = asrow["auth_seq_id"]
            anrow["pdbx_auth_comp_id"] = asrow["auth_comp_id"]
            anrow["pdbx_auth_asym_id"] = asrow["auth_asym_id"]
            anrow["pdbx_auth_atom_id"] = asrow["auth_atom_id"]
            anrow["pdbx_auth_alt_id"] = asrow["auth_alt_id"]
            anrow["U[1][1]"] = atm.U[0,0]
            anrow["U[2][2]"] = atm.U[1,1]
            anrow["U[3][3]"] = atm.U[2,2]
            anrow["U[1][2]"] = atm.U[0,1]
            anrow["U[1][3]"] = atm.U[0,2]
            anrow["U[2][3]"] = atm.U[1,2]

            if atm.sig_U:
                anrow["U[1][1]_esd"] = atm.sig_U[0,0]
                anrow["U[2][2]_esd"] = atm.sig_U[1,1]
                anrow["U[3][3]_esd"] = atm.sig_U[2,2]
                anrow["U[1][2]_esd"] = atm.sig_U[0,1]
                anrow["U[1][3]_esd"] = atm.sig_U[0,2]
                anrow["U[2][3]_esd"] = atm.sig_U[1,2]

