## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Builds a mmCIFFile object from a Structure object.
"""
from __future__ import generators
import copy
from mmTypes import *
from mmCIF import *


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
                    details = "unknown polymer"
                    type    = "other"
                
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

