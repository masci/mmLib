## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import string
import types

import mmCIF
import PDB

from Scientific.Geometry import Vector, Tensor
from Structure           import *


FileLoaderError = "FileLoaderError"


## various codes used in PDB/mmCIF files
InsertionCodes    = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
ChainIDs          = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
ConformationsIDs  = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


## columns used for writing the atom_site table in mmCIF
AtomSiteColumns = [
    "group_PDB",
    "id",
    "type_symbol",
    "label_atom_id",
    "label_alt_id",
    "label_comp_id",
    "label_asym_id",
    "label_entity_id", 
    "label_seq_id",
    "pdbx_PDB_ins_code",
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
    "auth_seq_id",
    "auth_comp_id",
    "auth_asym_id",
    "auth_atom_id",
    "pdbx_PDB_model_num"]


AtomSiteAnisotropColumns = [
    "id",
    "U[1][1]",
    "U[2][2]",
    "U[3][3]",
    "U[1][2]",
    "U[1][3]",
    "U[2][3]"]


def SortChainIDList(list):
    """Given a list of chain IDs this function sorts them into proper
    order.  The only difference from the standard list sort is to put
    a blank chain at the end of the list."""
    list.sort()
    if list and not list[0]: list = list[1:]
    list.append("")
    return list



class mmCIFToStructureLoader:
    """Constructs a mmLib.Structure data structure from a mmCIF file.  It
    uses the mmLib.mmCIF parser for parsing the mmCIF file."""

    def load(self, fil):
        self.structure = Structure()
        cif_file = mmCIF.mmCIFFile()
        cif_file.loadFile(fil)

        for cif_data in cif_file.getDataList():
            ## skip mmCIF files without atom_site tables
            if not hasattr(cif_data, "atom_site"):
                continue
            
            self.cif_data = cif_data
            self.load_cif_data()

        ## XXX: hack, construct bonds
        for pep_chain in self.structure.polypeptideChainIterator():
            pep_chain.constructPeptideBonds()
        
        for aa_res in self.structure.aminoAcidResidueIterator():
            aa_res.constructBonds()

        return self.structure

    
    def fragment_key(self, atom_site):
        """Construct a tuple as a unique set of values to lookup the fragment
        in the PDB records list."""

        ## these fields I consider mandatory
        asym_id = atom_site.label_asym_id

        ## process the sequence id/icode
        seq_id  = atom_site.auth_seq_id
        icode   = ""
        if seq_id and seq_id[-1] in InsertionCodes:
            icode = seq_id[-1]
            seq_id = seq_id[:-1]

        try:
            seq_id = int(seq_id)
        except ValueError:
            raise FileLoaderError, "atom_site.auth_seq_id invalid type"
        
        comp_id = ""

        ## these fields can be abscent
        try:
            comp_id = atom_site.auth_comp_id
        except AttributeError:
            try:
                resName = atom_site.label_comp_id
            except AttributeError:
                pass

        return (asym_id, seq_id, icode, comp_id)


    def load_cif_data(self):
        ## presort and create a some useful maps
        ##
        ## self.cif_data -> the current mmCIF data section we're loading
        ##
        ## self.atom_map
        ##   key:   tuple of values uniquely identifing a atom in the CIF file
        ##   value: a list of all  records giving additional information
        ##          on the atom
        ##
        ## self.fragment_map
        ##   key:   tuple of values uniquely identifing a fragment, usually
        ##          a Amino Acid residue
        ##   value: a list of all ATOM_SITE rows involved in the fragment
        ##
        ## self.chain_map
        ##   key:   the 1-char chain id
        ##   value: a list of all the fragment lists involved in the chain
        ##
        self.fragment_map = {}
        self.chain_map = {}

        for atom_site in self.cif_data.atom_site.getRowList():
            fkey = self.fragment_key(atom_site)                
            try:
                self.fragment_map[fkey].append(atom_site)
            except KeyError:
                self.fragment_map[fkey] = [atom_site]

            asym_id = fkey[0]
            try:
                self.chain_map[asym_id].append(self.fragment_map[fkey])
            except KeyError:
                self.chain_map[asym_id] = [self.fragment_map[fkey]]

        ##
        ## NOTE: asym_id == chain_id from now on
        ##

        ## create a list of the chain_id keys sorted alphabeticly
        chain_id_list = SortChainIDList(self.chain_map.keys())

        ## fragment_order is a list of tuples (resSeq, iCode)
        ## in-order for this chain
        ##
        ## python's brilliant sort routine will sort this just the
        ## way we want it: by chainID, resSeq, then iCode, then resName
        fragment_order = self.fragment_map.keys()
        fragment_order.sort()

        for chain_id in chain_id_list:

            cur_solvent  = None
            cur_chain    = None
            cur_molecule = None

            for fkey in fragment_order:
                (asym_id, seq_id, icode, comp_id) = fkey

                if asym_id != chain_id:
                    continue

                frag_type = self.figure_fragment_type(fkey)

                if   frag_type == "amino acid residue":

                    if cur_molecule == None:
                        cur_molecule = Molecule()
                        self.structure.appendChild(cur_molecule)
                        cur_molecule.setType("polymer")
                        cur_molecule.setPolymerType("polypeptide")

                    if cur_chain == None:
                        cur_chain = PolypeptideChain()
                        cur_molecule.appendChild(cur_chain)
                        cur_chain.setID(asym_id)
                    
                    aa_res = self.construct_amino_acid_residue(fkey)
                    cur_chain.appendChild(aa_res)


                elif frag_type == "nucleic acid residue":

                    if cur_molecule == None:
                        cur_molecule = Molecule()
                        self.structure.appendChild(cur_molecule)
                        cur_molecule.setType("polymer")
                        cur_molecule.setPolymerType("nucleicacid")

                    if cur_chain == None:
                        cur_chain = NucleicAcidChain()
                        cur_molecule.appendChild(cur_chain)
                        cur_chain.setID(asym_id)
                    
                    nres = self.construct_nucleic_acid_residue(fkey)
                    cur_chain.appendChild(nres)


                elif frag_type == "water":
                    wtr = self.construct_water(fkey)
                    if cur_solvent == None:
                        cur_solvent = Solvent()
                        self.structure.appendChild(cur_solvent)
                    cur_solvent.appendChild(wtr)
                    

                else:
                    unk = self.construct_unknown(fkey)
                    self.structure.appendChild(unk)


    def figure_fragment_type(self, fkey):
        """Returns the type of fragment: amino acid residue,
        nucleic acid residue, water, unknown.  It also preforms
        a sanity check to make sure all atoms listed have the same
        fragment name (resName)."""
        resname_list = []
        for atom_site in self.fragment_map[fkey]:
            try:
                comp_id = atom_site.auth_comp_id
            except AttributeError:
                try:
                    resName = atom_site.label_comp_id
                except AttributeError:
                    comp_id = ""

            if comp_id not in resname_list:
                resname_list.append(comp_id)
            
        if len(resname_list) > 1:
            raise FileLoaderError, "Bad mmCIF file."

        elif len(resname_list) == 0:
            return "unknown"

        resName = resname_list[0]
        if resName in AminoAcidNames:
            return "amino acid residue"

        elif resName in NucleicAcidNames:
            return "nucleic acid residue"

        elif resName == "HOH":
            return "water"

        return "unknown"


    def construct_unknown(self, fkey):
        (asym_id, seq_id, icode, comp_id) = fkey
        molecule = Molecule()
        molecule.setName(comp_id)
        molecule.setType("non-polymer")
        self.add_atoms(fkey, molecule)
        return molecule

                
    def construct_water(self, fkey):
        molecule = Molecule()
        molecule.setType("water")
        self.add_atoms(fkey, molecule)
        return molecule


    def construct_nucleic_acid_residue(self, fkey):
        (asym_id, seq_id, icode, comp_id) = fkey
        nres = NucleicAcidResidue()
        nres.setName(comp_id)
        nres.setSequenceID(str(seq_id) + icode)
        self.add_atoms(fkey, nres)
        return nres


    def construct_amino_acid_residue(self, fkey):
        (asym_id, seq_id, icode, comp_id) = fkey
        aa_res = AminoAcidResidue()
        aa_res.setName(comp_id)
        aa_res.setSequenceID(str(seq_id) + icode)
        self.add_atoms(fkey, aa_res)
        return aa_res


    def add_atoms(self, fkey, parent):
        conformation_map = {}

        for atom_site in self.fragment_map[fkey]:
            atom = self.construct_atom(atom_site)

            if hasattr(atom_site, "label_alt_id") and atom_site.label_alt_id:
                try:
                    alt_id = atom_site.label_alt_id
                except AttributeError:
                    alt_id = ""
                else:
                    if alt_id and alt_id not in InsertionCodes:
                        alt_id = ""

                if alt_id:
                    try:
                        conf = conformation_map[alt_id]
                    except KeyError:
                        conf = Conformation()
                        conformation_map[alt_id] = conf
                        parent.appendChild(conf)
                        conf.setID(alt_id)

                    conf.appendChild(atom)
                else:
                    parent.appendChild(atom)


    def construct_atom(self, atom_site):
        atom = Atom()
        
        atom.setElement(atom_site.type_symbol)
        atom.setAtomLabel(atom_site.label_atom_id)
        x = float(atom_site.Cartn_x)
        y = float(atom_site.Cartn_y)
        z = float(atom_site.Cartn_z)
        atom.setPosition(Vector(x, y, z))
        atom.setOccupancy(atom_site.occupancy)
        atom.setTemperatureFactor(atom_site.B_iso_or_equiv)

        ## set anisotropic values for the atom
        if hasattr(self.cif_data, "atom_site_anisotrop"):
            try:
                (aniso, ) = self.cif_data.atom_site_anisotrop.selectRowList(
                    ("id", atom_site.id))

            except ValueError:
                pass

            else:
                u00 = getattr(aniso, "U[1][1]")
                u11 = getattr(aniso, "U[2][2]")
                u22 = getattr(aniso, "U[3][3]")
                u01 = getattr(aniso, "U[1][2]")
                u02 = getattr(aniso, "U[1][3]")
                u12 = getattr(aniso, "U[2][3]")

                u00 = float(u00)
                u11 = float(u11)
                u22 = float(u22)
                u01 = float(u01)
                u02 = float(u02)
                u12 = float(u12)
                
                atom.setU( (u00, u11, u22, u01, u02, u12) )

        return atom



class PDBToStructureLoader:
    """Constructs a mmLib.Structure data structure from a PDB file.  It
    uses the mmLib.PDB parser for parsing the PDB file."""

    def load(self, fil):
        self.structure = Structure()
        self.pdb_file = PDB.PDBFile()
        self.pdb_file.loadFile(fil)

        ## presort and create a some useful maps
        ##
        ## self.pdb_record_list -> all non-atom oriented records
        ##
        ## self.atom_map
        ##   key:   tuple of values uniquely identifing a atom in the PDB file
        ##   value: a list of all PDB records giving additional information
        ##          on the atom
        ##
        ## self.fragment_map
        ##   key:   tuple of values uniquely identifing a fragment, usually
        ##          a Amino Acid residue
        ##   value: a list of all ATOM records involved in the fragment
        ##
        ## self.chain_map
        ##   key:   the 1-char chain id
        ##   value: a list of all the fragment lists involved in the chain
        ##
        self.pdb_record_list = []
        
        self.atom_map = {}
        self.fragment_map = {}
        self.chain_map = {}

        for rec in self.pdb_file.getRecordList():

            if isinstance(rec, PDB.ATOM):
                akey = self.atom_key(rec)
                try:
                    self.atom_map[akey].append(rec) 
                except KeyError:
                    self.atom_map[akey] = [rec]
                
                fkey = self.fragment_key(rec)                
                try:
                    self.fragment_map[fkey].append(rec)
                except KeyError:
                    self.fragment_map[fkey] = [rec]

                chainID = fkey[0]
                try:
                    self.chain_map[chainID].append(self.fragment_map[fkey])
                except KeyError:
                    self.chain_map[chainID] = [self.fragment_map[fkey]]


            elif isinstance(rec, PDB.ANISOU) or \
                 isinstance(rec, PDB.SIGATM) or \
                 isinstance(rec, PDB.SIGUIJ):
                
                akey = self.atom_key(rec)
                try:
                    self.atom_map[akey].append(rec) 
                except KeyError:
                    self.atom_map[akey] = [rec]

            else:
                self.pdb_record_list.append(rec)
                
        ## create a list of the chainID keys sorted alphabeticly
        chain_id_list = SortChainIDList(self.chain_map.keys())

        ## fragment_order is a list of tuples (resSeq, iCode)
        ## in-order for this chain
        ##
        ## python's brilliant sort routine will sort this just the
        ## way we want it: by chainID, resSeq, then iCode, then resName
        fragment_order = self.fragment_map.keys()
        fragment_order.sort()

        for chain_id in chain_id_list:

            cur_solvent  = None
            cur_chain    = None
            cur_molecule = None

            for fkey in fragment_order:
                (chainID, resSeq, iCode, resName) = fkey

                if chainID != chain_id:
                    continue

                frag_type = self.figure_fragment_type(fkey)


                if   frag_type == "amino acid residue":

                    if cur_molecule == None:
                        cur_molecule = Molecule()
                        self.structure.appendChild(cur_molecule)
                        cur_molecule.setType("polymer")
                        cur_molecule.setPolymerType("polypeptide")

                    if cur_chain == None:
                        cur_chain = PolypeptideChain()
                        cur_molecule.appendChild(cur_chain)
                        cur_chain.setID(chainID)
                    
                    aa_res = self.construct_amino_acid_residue(fkey)
                    cur_chain.appendChild(aa_res)


                elif frag_type == "nucleic acid residue":

                    if cur_molecule == None:
                        cur_molecule = Molecule()
                        self.structure.appendChild(cur_molecule)
                        cur_molecule.setType("polymer")
                        cur_molecule.setPolymerType("nucleicacid")

                    if cur_chain == None:
                        cur_chain = NucleicAcidChain()
                        cur_molecule.appendChild(cur_chain)
                        cur_chain.setID(chainID)
                    
                    nres = self.construct_nucleic_acid_residue(fkey)
                    cur_chain.appendChild(nres)


                elif frag_type == "water":
                    wtr = self.construct_water(fkey)
                    if cur_solvent == None:
                        cur_solvent = Solvent()
                        self.structure.appendChild(cur_solvent)
                    cur_solvent.appendChild(wtr)


                else:
                    unk = self.construct_unknown(fkey)
                    self.structure.appendChild(unk)

        ## XXX: hack, construct bonds
        for pep_chain in self.structure.polypeptideChainIterator():
            pep_chain.constructPeptideBonds()
        
        for aa_res in self.structure.aminoAcidResidueIterator():
            aa_res.constructBonds()

        return self.structure
    

    def fragment_key(self, atom_rec):
        """Construct a tuple as a unique set of values to lookup the fragment
        in the PDB records list."""
        ## these fields I consider mandatory
        name   = atom_rec.name
        resSeq = atom_rec.resSeq

        ## these fields can be abscent
        try:
            chainID = atom_rec.chainID
        except AttributeError:
            chainID = ""
        
        try:
            resName = atom_rec.resName
        except AttributeError:
            resName = ""

        try:
            iCode = atom_rec.iCode
        except AttributeError:
            iCode = ""

        return (chainID, resSeq, iCode, resName)


    def atom_key(self, rec):
        """Construct a tuple as a unique set of values to identify records
        describing a atom."""
        (chainID, resSeq, resName, iCode)  = self.fragment_key(rec)
        name = rec.name
        try:
            altLoc = rec.altLoc
        except AttributeError:
            altLoc = ""
        return (chainID, resSeq, iCode, resName, name, altLoc)


    def figure_fragment_type(self, fkey):
        """Returns the type of fragment: amino acid residue,
        nucleic acid residue, water, unknown.  It also preforms
        a sanity check to make sure all atoms listed have the same
        fragment name (resName)."""
        resname_list = []
        for atom_rec in self.fragment_map[fkey]:
            if not hasattr(atom_rec, "resName"):
                continue
            if atom_rec.resName not in resname_list:
                resname_list.append(atom_rec.resName)
            
        if len(resname_list) > 1:
            raise FileLoaderError, "Bad PDB file."

        elif len(resname_list) == 0:
            return "unknown"

        resName = resname_list[0]

        if resName in AminoAcidNames:
            return "amino acid residue"

        elif resName == "HOH":
            return "water"

        elif resName in NucleicAcidNames:
            return "nucleic acid residue"

        return "unknown"


    def construct_unknown(self, fkey):
        (chainID, resSeq, iCode, resName) = fkey
        molecule = Molecule()
        molecule.setName(resName)
        molecule.setType("non-polymer")
        self.add_atoms(fkey, molecule)
        return molecule

                
    def construct_water(self, fkey):
        molecule = Molecule()
        molecule.setType("water")
        self.add_atoms(fkey, molecule)
        return molecule


    def construct_nucleic_acid_residue(self, fkey):
        (chainID, resSeq, iCode, resName) = fkey
        nres = NucleicAcidResidue()
        nres.setName(resName)
        nres.setSequenceNumber(resSeq)
        nres.setSequenceID(str(resSeq) + iCode)
        self.add_atoms(fkey, nres)
        return nres


    def construct_amino_acid_residue(self, fkey):
        (chainID, resSeq, iCode, resName) = fkey
        aa_res = AminoAcidResidue()
        aa_res.setName(resName)
        aa_res.setSequenceNumber(resSeq)
        aa_res.setSequenceID(str(resSeq) + iCode)
        self.add_atoms(fkey, aa_res)
        return aa_res


    def add_atoms(self, fkey, parent):
        conformation_map = {}

        for atom_rec in self.fragment_map[fkey]:
            atom = self.construct_atom(atom_rec)

            if hasattr(atom_rec, "altLoc") and atom_rec.altLoc:
                try:
                    conf = conformation_map[atom_rec.altLoc]
                except KeyError:
                    conf = conformation_map[atom_rec.altLoc] = Conformation()
                    parent.appendChild(conf)
                    conf.setID(atom_rec.altLoc)

                conf.appendChild(atom)
            else:
                parent.appendChild(atom)


    def construct_atom(self, atom_rec):
        atom = Atom()
        akey = self.atom_key(atom_rec)

        for rec in self.atom_map[akey]:
            if isinstance(rec, PDB.ATOM):
                ## some PDB files do not have the proper element column, so
                ## the name of the element can be taken from the first letter
                ## of the atom name(label)
                try:
                    atom.setElement(rec.element)
                except AttributeError:
                    atom.setElement(rec.name[0])

                atom.setAtomLabel(rec.name)
                atom.setPosition(Vector(rec.x, rec.y, rec.z))
                atom.setOccupancy(rec.occupancy)
                atom.setTemperatureFactor(rec.tempFactor)

                ## some PDB files do not have charge columns
                try:
                    atom.setCharge(rec.charge)
                except AttributeError:
                    pass

            elif isinstance(rec, PDB.ANISOU):
                try:
                    u00 = getattr(rec, "u[0][0]")
                    u11 = getattr(rec, "u[1][1]")
                    u22 = getattr(rec, "u[2][2]")
                    u01 = getattr(rec, "u[0][1]")
                    u02 = getattr(rec, "u[0][2]")
                    u12 = getattr(rec, "u[1][2]")
                except AttributeError:
                    continue

                ## ANISOU records are scaled by 10^4 and presented as
                ## integers; we need to convert them back
                u00 = u00 / 10000.0
                u11 = u11 / 10000.0
                u22 = u22 / 10000.0
                u01 = u01 / 10000.0
                u02 = u02 / 10000.0
                u12 = u12 / 10000.0

                atom.setU( (u00, u11, u22, u01, u02, u12) )

            else:
                print "PDB::construct_atom ignoring record=%s" % (rec._name)

        return atom



class StructureTommCIFSaver:
    """Converts Structure to mmCIFFile class and saves it to the given path."""

    def save(self, fil, structure):
        self.cif_file = mmCIF.mmCIFFile()
        self.structure = structure

        ## waters shall always be in chain S
        self.water_chain = "S"
        self.water_res_seq = 1

        ## list of chain ID's available for use
        ## "S" should not be in the list -- it's the default solvent chain
        self.available_chain_ids = ChainIDs

        ## first step removes the chainID's in use
        for chain in self.structure.chainIterator():
            if chain.getID() in self.available_chain_ids:
                i = self.available_chain_ids.index(chain.getID())
                self.available_chain_ids = self.available_chain_ids[i+1:]

        ## create the mmCIF data block and mmCIF tables we need
        self.cif_data = mmCIF.mmCIFData("mol")
        self.cif_file.addData(self.cif_data)
                
        tbl = mmCIF.mmCIFTable("atom_site", AtomSiteColumns)
        self.cif_data.addTable(tbl)
        
        tbl = mmCIF.mmCIFTable("atom_site_anisotrop",AtomSiteAnisotropColumns)
        self.cif_data.addTable(tbl)

        ## write out the structure
        for mol in self.structure.moleculeIterator():
            if mol.getType() == "polymer":
                self.molecule_polymer(mol)

            elif mol.getType() == "non-polymer":
                self.molecule_non_polymer(mol)

            elif mol.getType() == "water":
                self.molecule_water(mol)
                
        self.cif_file.saveFile(fil)
    

    def molecule_polymer(self, mol):
        for chain in mol.chainIterator():
            label_seq_id = 1
            
            for aa_res in chain.aminoAcidResidueIterator():

                for atm in aa_res.atomIterator():
                    label_comp_id = aa_res.getName()
                    label_asym_id = chain.getID()
                    auth_seq_id   = aa_res.getSequenceID()
                    
                    label_alt_id  = ""
                    conf = atm.getConformation()
                    if conf: label_alt_id = conf.getID()

                    self.add_atom(atm, label_comp_id, label_asym_id,
                                  label_seq_id, label_alt_id, auth_seq_id)

                label_seq_id += 1


    def molecule_non_polymer(self, mol):
        ## get a chain ID for the molecule
        label_asym_id = ""
        if len(self.available_chain_ids):
            label_asym_id = self.available_chain_ids[0]
            self.available_chain_ids = self.available_chain_ids[1:]

        for atm in mol.atomIterator():
            label_comp_id = mol.getName() or "UNK"
            label_seq_id  = 1
            auth_seq_id   = 1

            label_alt_id  = ""
            conf = atm.getConformation()
            if conf: label_alt_id = conf.getID()

            self.add_hetatm(atm, label_comp_id, label_asym_id,
                            label_seq_id, label_alt_id, auth_seq_id)


    def molecule_water(self, mol):
        for atm in mol.atomIterator():
            label_alt_id  = ""
            conf = atm.getConformation()
            if conf: label_alt_id = conf.getID()

            self.add_hetatm(atm, "HOH", self.water_chain,
                            self.water_res_seq, label_alt_id,
                            self.water_res_seq)

        self.water_res_seq += 1


    def add_hetatm(self, atm, label_comp_id, label_asym_id,
                   label_seq_id, label_alt_id, auth_seq_id):
        self.add_atom2(atm, label_comp_id, label_asym_id,
                       label_seq_id, label_alt_id, auth_seq_id, "HETATM")


    def add_atom(self, atm, label_comp_id, label_asym_id,
                 label_seq_id, label_alt_id, auth_seq_id):
        self.add_atom2(atm, label_comp_id, label_asym_id,
                       label_seq_id, label_alt_id, auth_seq_id, "ATOM")


    def add_atom2(self, atm, label_comp_id, label_asym_id,
                  label_seq_id, label_alt_id, auth_seq_id, group_PDB):

        row = mmCIF.mmCIFRow()
        self.cif_data.atom_site.addRow(row)

        ## unique serial number/id for atom
        row.id             = 0

        ## residue/entity/monomer properties
        row.group_PDB      = group_PDB
        row.label_comp_id  = label_comp_id
        row.label_asym_id  = label_asym_id
        row.label_seq_id   = label_seq_id
        row.label_alt_id   = label_alt_id
        row.auth_seq_id    = auth_seq_id

        ## atom properties
        row.type_symbol    = atm.getElement()
        row.label_atom_id  = atm.getAtomLabel()

        vec                = atm.getPosition()
        row.Cartn_x        = vec[0]
        row.Cartn_y        = vec[1]
        row.Cartn_z        = vec[2]
        
        row.occupancy      = atm.getOccupancy()
        row.B_iso_or_equiv = atm.getTemperatureFactor()
        
        if atm.getU():
            (u00, u11, u22, u01, u02, u12) = atm.getU()
            
            urow = mmCIF.mmCIFRow()
            self.cif_data.atom_site_anisotrop.addRow(urow)

            urow.id = 0
            setattr(urow, "U[1][1]", u00)
            setattr(urow, "U[2][2]", u11)
            setattr(urow, "U[3][3]", u22)
            setattr(urow, "U[1][2]", u01)
            setattr(urow, "U[1][3]", u02)
            setattr(urow, "U[2][3]", u12)



class StructureToPDBSaver:
    def save(self, fil, structure):
        self.pdb_file = PDB.PDBFile()
        self.structure = structure

        ## waters shall always be in chain S
        self.water_chain = "S"
        self.water_res_seq = 1

        ## list of chain ID's available for use
        ## "S" should not be in the list -- it's the default solvent chain
        self.available_chain_ids = ChainIDs

        ## first step removes the chainID's in use
        for chain in self.structure.chainIterator():
            if chain.getID() in self.available_chain_ids:
                i = self.available_chain_ids.index(chain.getID())
                self.available_chain_ids = self.available_chain_ids[i+1:]

        ## write out the structure
        for mol in self.structure.moleculeIterator():
            if mol.getType() == "polymer":
                self.molecule_polymer(mol)

            elif mol.getType() == "non-polymer":
                self.molecule_non_polymer(mol)

            elif mol.getType() == "water":
                self.molecule_water(mol)
                
        self.pdb_file.saveFile(fil)
    

    def molecule_polymer(self, mol):
        for chain in mol.chainIterator():
            for aa_res in chain.aminoAcidResidueIterator():
                for atm in aa_res.atomIterator():
                    resName = aa_res.getName()
                    chainID = chain.getID()
                    resSeq  = aa_res.getSequenceID()
                    iCode   = ""
                    altLoc  = ""
                    
                    conf = atm.getConformation()
                    if conf:
                        altLoc = conf.getID()

                    if resSeq and resSeq[-1] in 'ABCDEFGHIJKLMNOPQRSTIVWXYZ':
                        iCode  = resSeq[-1]
                        resSeq = resSeq[:-1]
                    
                    self.add_atom(atm, resName, chainID, resSeq, iCode, altLoc)


    def molecule_non_polymer(self, mol):
        ## get a chain ID for the molecule
        chainID = ""
        if len(self.available_chain_ids):
            chainID = self.available_chain_ids[0]
            self.available_chain_ids = self.available_chain_ids[1:]

        for atm in mol.atomIterator():
            resName = mol.getName() or "UNK"
            resSeq  = 1
            iCode   = ""
            altLoc  = ""

            conf = atm.getConformation()
            if conf:
                altLoc = conf.getID()
            
            self.add_hetatm(atm, resName, chainID, resSeq, iCode, altLoc)


    def molecule_water(self, mol):
        for atm in mol.atomIterator():
            altLoc = ""
            conf = atm.getConformation()
            if conf:
                altLoc = conf.getID()
            
            self.add_hetatm(atm, "HOH", self.water_chain,
                            self.water_res_seq, "", altLoc)
            self.water_res_seq += 1


    def add_hetatm(self, atm, resName, chainID, resSeq, iCode, altLoc):
        self.add_atom2(atm, resName, chainID, resSeq, iCode, altLoc,
                       PDB.HETATM())


    def add_atom(self, atm, resName, chainID, resSeq, iCode, altLoc):
        self.add_atom2(atm, resName, chainID, resSeq, iCode, altLoc,
                       PDB.ATOM())
        

    def add_atom2(self, atm, resName, chainID, resSeq, iCode, altLoc,
                  pdb_atom):
        
        self.pdb_file.pdbrecord_list.append(pdb_atom)

        pdb_atom.name = atm.getAtomLabel()
        pdb_atom.resName = resName
        pdb_atom.altLoc = altLoc
        pdb_atom.chainID = chainID
        pdb_atom.resSeq = resSeq
        pdb_atom.iCode = iCode
        
        vec = atm.getPosition()
        pdb_atom.x = vec[0]
        pdb_atom.y = vec[1]
        pdb_atom.z = vec[2]
        
        pdb_atom.occupancy = atm.getOccupancy()
        pdb_atom.tempFactor = atm.getTemperatureFactor()
        pdb_atom.element = atm.getElement()
        
        if atm.getU() != None:
            pdb_anisou = PDB.ANISOU()
            self.pdb_file.pdbrecord_list.append(pdb_anisou)

            pdb_anisou.name = pdb_atom.name
            pdb_anisou.resName = pdb_atom.resName
            pdb_anisou.chainID = pdb_atom.chainID
            pdb_anisou.resSeq = pdb_atom.resSeq
            pdb_anisou.iCode = pdb_atom.iCode
            pdb_anisou.element = pdb_atom.element

            (u00, u11, u22, u01, u02, u12) = atm.getU()
            u00 = int(u00 * 10000.0)
            u11 = int(u11 * 10000.0)
            u22 = int(u22 * 10000.0)
            u01 = int(u01 * 10000.0)
            u02 = int(u02 * 10000.0)
            u12 = int(u12 * 10000.0)

            setattr(pdb_anisou, "u[0][0]", u00)
            setattr(pdb_anisou, "u[1][1]", u11)
            setattr(pdb_anisou, "u[2][2]", u22)
            setattr(pdb_anisou, "u[0][1]", u01)
            setattr(pdb_anisou, "u[0][2]", u02)
            setattr(pdb_anisou, "u[1][2]", u12)



def GetAdaptor(path, mode, format):
    """Returns the correct adaptor class by inspecting the format
    string, which is the extention of the fil being saved or loaded."""
    
    ## check for compressed file
    ext = format
    
    if type(path) == types.StringType and format == "":
        (base, ext) = os.path.splitext(path)
        ext = string.lower(ext)

    elif type(path) != types.StringType and format == "":
        ## default format
        format = ".pdb"

    ## remove the compression extention if it is there
    if ext in ['.z', '.gz', '.bz2']:
        (base, ext) = os.path.splitext(base)
        ext = string.lower(ext)

    ## select adaptor based on file extention
    if ext == '.cif':
        if mode == "r": return mmCIFToStructureLoader
        if mode == "w": return StructureTommCIFSaver

    ## defualt to PDB
    else:
        if mode == "r": return PDBToStructureLoader
        if mode == "w": return StructureToPDBSaver



def LoadStructure(path_or_fil, format = ""):
    """Loads a mmCIF file(.cif) or PDB file(.pdb) into a mmPython
    Structure class and returns it.  The first argument is either a
    path string or a file object opened for reading.  The second argument
    is used to override the default file format."""
    adaptor_class = GetAdaptor(path_or_fil, "r", format)
    adaptor = adaptor_class()
    return adaptor.load(path_or_fil)



def SaveStructure(path_or_fil, structure, format = ""):
    """Saves a mmPython Structure class into a supported file type."""
    adaptor_class = GetAdaptor(path_or_fil, "w", format)
    adaptor = adaptor_class()
    adaptor.save(path_or_fil, structure)
