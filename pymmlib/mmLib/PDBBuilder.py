## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Convert a Structure object to a PDBFile object.
"""
from __future__ import generators
import string
import types
import fpformat
from PDB import *


class PDBFileBuilder(object):
    """Builds a PDBFile object from a Structure object.
    """
    def __init__(self, struct, pdb_file):
        self.struct = struct
        self.pdb_file = pdb_file

        self.atom_count = 0
        self.atom_serial_num = 0
        self.atom_serial_map = {}

        self.add_title_section()
        self.add_primary_structure_section()
        self.add_heterogen_section()
        self.add_secondary_structure_section()
        self.add_connectivity_annotation_section()
        self.add_miscellaneous_fatures_section()
        self.add_crystallographic_oordinate_transformation_section()
        self.add_coordinate_section()
        self.add_connectivity_section()
        self.bookkeeping_section()

    def next_serial_number(self):
        self.atom_serial_num += 1
        return self.atom_serial_num

    def new_atom_serial(self, atm):
        """Gets the next available atom serial number for the given atom
        instance, and stores a map from atm->atom_serial_num for use
        when creating PDB records which require serial number identification
        of the atoms.
        """
        assert isinstance(atm, Atom)
        
        try:
            return self.atom_serial_map[atm]
        except KeyError:
            pass
        atom_serial_num = self.next_serial_number()
        self.atom_serial_map[atm] = atom_serial_num
        return atom_serial_num

    def set_from_cifdb(self, rec, field, ctbl, ccol):
        try:
            rec[field] = self.struct.cifdb[ctbl][ccol]
        except KeyError:
            pass

    def add_title_section(self):
        """ HEADER, TITLE, EXPDTA, AUTHOR
        """
        header = HEADER()
        self.pdb_file.append(header)
        self.set_from_cifdb(header, "idCode",
                            "entry", "id")
        self.set_from_cifdb(header, "depDate",
                            "database_pdb_rev", "date_original")
        self.set_from_cifdb(header, "classification",
                            "struct_keywords", "pdbx_keywords")
        
        title = TITLE()
        self.pdb_file.append(title)
        self.set_from_cifdb(title, "title", "struct", "title")

    def add_primary_structure_section(self):
        """DBREF,SEQADV,SEQRES,MODRES
        """
        pass
    
    def add_heterogen_section(self):
        """HET,HETNAM,HETSYN,FORMUL
        """
        pass

    def add_secondary_structure_section(self):
        """HELIX,SHEET,TURN
        """
        pass

    def add_connectivity_annotation_section(self):
        """SSBOND,LINK,SLTBRG,CISPEP
        """
        pass

    def add_miscellaneous_fatures_section(self):
        """SITE
        """
        pass

    def add_crystallographic_oordinate_transformation_section(self):
        """CRYST1,ORIGXn,SCALEn,MTRIXn,TVECT
        """
        cryst1 = CRYST1()
        self.pdb_file.append(cryst1)

        unit_cell = self.struct.unit_cell

        cryst1["a"] = self.struct.unit_cell.a
        cryst1["b"] = self.struct.unit_cell.b
        cryst1["c"] = self.struct.unit_cell.c
        cryst1["alpha"] = self.struct.unit_cell.calc_alpha_deg()
        cryst1["beta"] = self.struct.unit_cell.calc_beta_deg()
        cryst1["gamma"] = self.struct.unit_cell.calc_gamma_deg()
        cryst1["sgroup"] = self.struct.unit_cell.space_group.pdb_name

    def add_coordinate_section(self):
        """ MODEL,ATOM,SIGATM,ANISOU,SIGUIJ,TER,HETATM,ENDMDL 
        """
        if len(self.struct.model_list) > 1:
            ## case 1: multiple models
            orig_model = self.struct.model
            
            for model in self.struct.iter_models():
                self.struct.model = model

                model_rec = MODEL()
                self.pdb_file.append(model_rec)
                model_rec["serial"] = model.model_id

                self.add_atom_records()

                endmdl = ENDMDL()
                self.pdb_file.append(endmdl)

            self.struct.model = orig_model

        else:
            ## case 2: single model
            self.add_atom_records()

    def add_connectivity_section(self):
        """CONECT
        """
        pass
    
    def bookkeeping_section(self):
        """MASTER,END
        """
        ## END
        end = END()
        self.pdb_file.append(end)


    def add_atom_records(self):
        """With a default model set, output all the ATOM and associated
        records for the model.
        """
        ## atom records for standard groups
        for chain in self.struct.iter_chains():
            res = None
            
            for res in chain.iter_standard_residues():
                for atm in res.iter_all_atoms():
                    self.add_ATOM("ATOM", atm)

            ## chain termination record
            if res:
                ter_rec = TER()
                self.pdb_file.append(ter_rec)
                fid = FragmentID(res.fragment_id)
                ter_rec["serial"]  = self.next_serial_number()
                ter_rec["resName"] = res.res_name
                ter_rec["chainID"] = res.chain_id
                ter_rec["resSeq"]  = fid.res_seq
                ter_rec["iCode"]   = fid.icode

        ## hetatm records for non-standard groups
        for chain in self.struct.iter_chains():
            for frag in chain.iter_non_standard_residues():
                for atm in frag.iter_all_atoms():
                    self.add_ATOM("HETATM", atm)

    def add_ATOM(self, rec_type, atm):
        """Adds ATOM/SIGATM/ANISOU/SIGUIJ/TER/HETATM 
        """
        self.atom_count += 1

        if rec_type == "ATOM":
            atom_rec = ATOM()
        elif rec_type == "HETATM":
            atom_rec = HETATM()

        self.pdb_file.append(atom_rec)

        serial = self.new_atom_serial(atm)
        fid = FragmentID(atm.fragment_id)

        atom_rec["serial"]      = serial
        atom_rec["chainID"]     = atm.chain_id
        atom_rec["resName"]     = atm.res_name
        atom_rec["resSeq"]      = fid.res_seq
        atom_rec["iCode"]       = fid.icode
        atom_rec["name"]        = atm.name
        atom_rec["element"]     = atm.element
        atom_rec["altLoc"]      = atm.alt_loc
        atom_rec["x"]           = atm.position[0]
        atom_rec["y"]           = atm.position[1]
        atom_rec["z"]           = atm.position[2]
        atom_rec["occupancy"]   = atm.occupancy
        atom_rec["tempFactor"]  = atm.temp_factor
        atom_rec["charge"]      = atm.charge

        def atom_common(arec1, arec2):
            arec2["serial"]  = arec1["serial"]
            arec2["chainID"] = arec1["chainID"]
            arec2["resName"] = arec1["resName"]
            arec2["resSeq"]  = arec1["resSeq"]
            arec2["iCode"]   = arec1["iCode"]
            arec2["name"]    = arec1["name"]
            arec2["altLoc"]  = arec1["altLoc"]
            arec2["element"] = arec1["element"]
            arec2["charge"]  = arec1["charge"]

        if atm.sig_position != None:
            sigatm_rec = SIGATM()
            self.pdb_file.append(sigatm_rec)
            atom_common(atom_rec, sigatm_rec)
            sigatm_rec["sigX"] = atm.sig_position[0]
            sigatm_rec["sigY"] = atm.sig_position[1]
            sigatm_rec["sigZ"] = atm.sig_position[2]
            sigatm_rec["sigOccupancy"] = atm.sig_temp_factor
            sigatm_rec["sigTempFactor"] = atm.sig_occupancy

        if atm.U != None:
            anisou_rec = ANISOU()
            self.pdb_file.append(anisou_rec)
            atom_common(atom_rec, anisou_rec)
            anisou_rec["u[0][0]"] = int(round(atm.U[0,0] * 10000.0))
            anisou_rec["u[1][1]"] = int(round(atm.U[1,1] * 10000.0))
            anisou_rec["u[2][2]"] = int(round(atm.U[2,2] * 10000.0))
            anisou_rec["u[0][1]"] = int(round(atm.U[0,1] * 10000.0))
            anisou_rec["u[0][2]"] = int(round(atm.U[0,2] * 10000.0))
            anisou_rec["u[1][2]"] = int(round(atm.U[1,2] * 10000.0))

        if atm.sig_U != None:
            siguij_rec = SIGUIJ()
            self.pdb_file.append(siguij_rec)
            atom_common(atom_rec, siguij_rec)
            siguij_rec["u[0][0]"] = int(round(atm.U[0,0] * 10000.0))
            siguij_rec["u[1][1]"] = int(round(atm.U[1,1] * 10000.0))
            siguij_rec["u[2][2]"] = int(round(atm.U[2,2] * 10000.0))
            siguij_rec["u[0][1]"] = int(round(atm.U[0,1] * 10000.0))
            siguij_rec["u[0][2]"] = int(round(atm.U[0,2] * 10000.0))
            siguij_rec["u[1][2]"] = int(round(atm.U[1,2] * 10000.0))

