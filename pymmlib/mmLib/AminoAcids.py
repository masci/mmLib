## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

"""Library of amino acid residues."""

from Library import AminoAcid

AminoAcidNames = [
    "ALA", "VAL", "LEU", "ILE", "GLY", "PRO", "CYS", "MET", "HIS",
    "PHE", "TYR", "TRP", "ASN", "GLN", "SER", "THR", "LYS", "ARG",
    "ASP", "GLU", "DIS", "UNK" ]


## <ALANINE>
ALA = AminoAcid(
    name            = "ALA",
    full_name       = "Alanine",
    one_letter_name = "A",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("HB3", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "HB3"),
                       ("CA", "C"),
                       ("C", "O")])
## </ALANINE>


### <VALINE>
VAL = AminoAcid(
    name            = "VAL",
    full_name       = "Valine",
    one_letter_name = "V",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB", "H"),
                       ("CG1", "C"),
                       ("HG11", "H"),
                       ("HG12", "H"),
                       ("HG13", "H"),
                       ("CG2", "C"),
                       ("HG21", "H"),
                       ("HG22", "H"),
                       ("HG23", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB"),
                       ("CB", "CG1"),
                       ("CG1", "HG11"),
                       ("CG1", "HG12"),
                       ("CG1", "HG13"),
                       ("CB", "CG2"),
                       ("CG2", "HG21"),
                       ("CG2", "HG22"),
                       ("CG2", "HG23"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG1"))
### </VALINE>


### <LEUSINE>
LEU = AminoAcid(
    name            = "LEU",
    full_name       = "Leusine",
    one_letter_name = "L",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("HG", "H"),
                       ("CD1", "C"),
                       ("HD11", "H"),
                       ("HD12", "H"),
                       ("HD13", "H"),
                       ("CD2", "C"),
                       ("HD21", "H"),
                       ("HD22", "H"),
                       ("HD23", "H"),
                       ("C", "C"),
                       ("O", "O")],

    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "HG"),
                       ("CG", "CD1"),
                       ("CD1", "HD11"),
                       ("CD1", "HD12"),
                       ("CD1", "HD13"),
                       ("CG", "CD2"),
                       ("CD2", "HD21"),
                       ("CD2", "HD22"),
                       ("CD2", "HD23"),
                       ("CA", "C"),
                       ("C", "O")],

    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "CD1"))
### </LEUSINE>


### <ISOLEUSINE>
ILE = AminoAcid(
    name            = "ILE",
    full_name       = "Isoleusine",
    one_letter_name = "I",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB", "H"),
                       ("CG1", "C"),
                       ("HG11", "H"),
                       ("HG12", "H"),
                       ("CD1", "C"),
                       ("HD11", "H"),
                       ("HD12", "H"),
                       ("HD13", "H"),
                       ("CG2", "C"),
                       ("HG21", "H"),
                       ("HG22", "H"),
                       ("HG23", "H"),
                       ("C", "C"),
                       ("O", "O")],

    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB"),
                       ("CB", "CG1"),
                       ("CG1", "HG11"),
                       ("CG1", "HG12"),
                       ("CG1", "CD1"),
                       ("CD1", "HD11"),
                       ("CD1", "HD12"),
                       ("CD1", "HD13"),
                       ("CB", "CG2"),
                       ("CG2", "HG21"),
                       ("CG2", "HG22"),
                       ("CG2", "HG23"),
                       ("CA", "C"),
                       ("C", "O")],

    chi1_definition = ("N", "CA", "CB", "CG1"),
    chi2_definition = ("CA", "CB", "CG1", "CD1"))
### <ISOLEUSINE>


### <GLYSINE>
GLY = AminoAcid(
    name            = "GLY",
    full_name       = "Glysine",
    one_letter_name = "G",
    
    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA1", "H"),
                       ("HA2", "H"),
                       ("C", "C"),
                       ("O", "O")],

    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA1"),
                       ("CA", "HA2"),
                       ("CA", "C"),
                       ("C", "O")])
### </GLYSINE>


### <PROLINE>
PRO = AminoAcid(
    name            = "PRO",
    full_name       = "Proline",
    one_letter_name = "P",

    atom_list       = [("N", "N"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("HG1", "H"),
                       ("HG2", "H"),
                       ("CD", "C"),
                       ("HD1", "H"),
                       ("HD2", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "HG1"),
                       ("CG", "HG2"),
                       ("CG", "CD"),
                       ("CD", "HD1"),
                       ("CD", "HD2"),
                       ("CD", "N"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition   = ("N", "CA", "CB", "CG"),
    pucker_definition = ("C", "CA", "CB", "CG"))
### </PROLINE>


### <CYSTINE>
CYS = AminoAcid(
    name            = "CYS",
    full_name       = "Cystine",
    one_letter_name = "C",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("SG", "S"),
                       ("C", "C"),
                       ("O", "O")],

    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "SG"),
                       ("CA", "C"),
                       ("C", "O")],

    chi1_definition = ("N", "CA", "CB", "SG"))
### </CYSTINE>


### <METHIONINE>
MET = AminoAcid(
    name            = "MET",
    full_name       = "Methionine",
    one_letter_name = "M",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("HG1", "H"),
                       ("HG2", "H"),
                       ("SD", "S"),
                       ("CE", "C"),
                       ("HE1", "H"),
                       ("HE2", "H"),
                       ("HE3", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "HG1"),
                       ("CG", "HG2"),
                       ("CG", "SD"),
                       ("SD", "CE"),
                       ("CE", "HE1"),
                       ("CE", "HE2"),
                       ("CE", "HE3"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "SD"),
    chi3_definition = ("CB", "CG", "SD", "CE"))
### </METHIONINE>


### <HISTIDINE>
HIS = AminoAcid(
    name            = "HIS",
    full_name       = "Histidine",
    one_letter_name = "H",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("ND1", "N"),
                       ("HD1", "H"),
                       ("CE1", "C"),
                       ("HE1", "H"),
                       ("NE2", "N"),
                       ("HE2", "H"),
                       ("CD2", "C"),
                       ("HD2", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "CD2"),
                       ("CG", "ND1"),
                       ("ND1", "HD1"),
                       ("ND1", "CE1"),
                       ("CE1", "HE1"),
                       ("CE1", "NE2"),
                       ("NE2", "HE2"),
                       ("NE2", "CD2"),
                       ("CD2", "HD2"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "ND1"))
### </HISTIDINE>


### <PHENYLALALINE>
PHE = AminoAcid(
    name            = "PHE",
    full_name       = "Phenylalaline",
    one_letter_name = "F",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("CD1", "C"),
                       ("HD1", "H"),
                       ("CE1", "C"),
                       ("HE1", "H"),
                       ("CZ", "C"),
                       ("HZ", "H"),
                       ("CE2", "C"),
                       ("HE2", "H"),
                       ("CD2", "C"),
                       ("HD2", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "CD2"),
                       ("CG", "CD1"),
                       ("CD1", "HD1"),
                       ("CD1", "CE1"),
                       ("CE1", "HE1"),
                       ("CE1", "CZ"),
                       ("CZ", "HZ"),
                       ("CZ", "CE2"),
                       ("CE2", "HE2"),
                       ("CE2", "CD2"),
                       ("CD2", "HD2"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "CD1"))
### </PHENYLALALINE>


### <TYROSINE>
TYR = AminoAcid(
    name            = "TYR",
    full_name       = "Tyrosine",
    one_letter_name = "Y",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("CD1", "C"),
                       ("HD1", "H"),
                       ("CE1", "C"),
                       ("HE1", "H"),
                       ("CZ", "C"),
                       ("OH", "O"),
                       ("HH", "H"),
                       ("CE2", "C"),
                       ("HE2", "H"),
                       ("CD2", "C"),
                       ("HD2", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "CD2"),
                       ("CG", "CD1"),
                       ("CD1", "HD1"),
                       ("CD1", "CE1"),
                       ("CE1", "HE1"),
                       ("CE1", "CZ"),
                       ("CZ", "OH"),
                       ("OH", "HH"),
                       ("CZ", "CE2"),
                       ("CE2", "HE2"),
                       ("CE2", "CD2"),
                       ("CD2", "HD2"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "CD1"))
### </TYROSINE>


### <TRPTOPHAN>
TRP = AminoAcid(
    name            = "TRP",
    full_name       = "Trptophan",
    one_letter_name = "W",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("CD1", "C"),
                       ("HD1", "H"),
                       ("NE1", "N"),
                       ("HE1", "H"),
                       ("CE2", "C"),
                       ("CD2", "C"),
                       ("CE3", "C"),
                       ("HE3", "H"),
                       ("CZ3", "C"),
                       ("HZ3", "H"),
                       ("CH2", "C"),
                       ("HH2", "H"),
                       ("CZ2", "C"),
                       ("HZ2", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "CD2"),
                       ("CG", "CD1"),
                       ("CD1", "HD1"),
                       ("CD1", "NE1"),
                       ("NE1", "HE1"),
                       ("NE1", "CE2"),
                       ("CE2", "CZ2"),
                       ("CE2", "CD2"),
                       ("CD2", "CE3"),
                       ("CE3", "HE3"),
                       ("CE3", "CZ3"),
                       ("CZ3", "HZ3"),
                       ("CZ3", "CH2"),
                       ("CH2", "HH2"),
                       ("CH2", "CZ2"),
                       ("CZ2", "HZ2"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "CD1"))
### </TRPTOPHAN>


### <ASPARAGINE>
ASN = AminoAcid(
    name            = "ASN",
    full_name       = "Asparagine",
    one_letter_name = "N",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("OD1", "O"),
                       ("ND2", "N"),
                       ("HD21", "H"),
                       ("HD22", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "OD1"),
                       ("CG", "ND2"),
                       ("ND2", "HD21"),
                       ("ND2", "HD22"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "OD1"))
### </ASPARAGINE>


### <GLUTAMINE>
GLN = AminoAcid(
    name            = "GLN",
    full_name       = "Glutamine",
    one_letter_name = "Q",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("HG1", "H"),
                       ("HG2", "H"),
                       ("CD", "C"),
                       ("OE1", "O"),
                       ("NE2", "N"),
                       ("HE21", "H"),
                       ("HE22", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "HG1"),
                       ("CG", "HG2"),
                       ("CG", "CD"),
                       ("CD", "OE1"),
                       ("CD", "NE2"),
                       ("NE2", "HE21"),
                       ("NE2", "HE22"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "CD"),
    chi3_definition = ("CB", "CG", "CD", "OE1"))
### </GLUTAMINE>


### <SERINE>
SER = AminoAcid(
    name            = "SER",
    full_name       = "Serine",
    one_letter_name = "S",


    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("OG", "O"),
                       ("HG", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "OG"),
                       ("OG", "HG"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "OG"))
### <SERINE>


### <THREONINE>
THR = AminoAcid(
    name            = "THR",
    full_name       = "Threonine",
    one_letter_name = "T",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB", "H"),
                       ("OG1", "O"),
                       ("HG1", "H"),
                       ("CG2", "C"),
                       ("HG21", "H"),
                       ("HG22", "H"),
                       ("HG23", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB"),
                       ("CB", "OG1"),
                       ("OG1", "HG1"),
                       ("CB", "CG2"),
                       ("CG2", "HG21"),
                       ("CG2", "HG22"),
                       ("CG2", "HG23"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "OG1"))
### <THREONINE>


### <LYSINE>
LYS = AminoAcid(
    name            = "LYS",
    full_name       = "lysine",
    one_letter_name = "K",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("HG1", "H"),
                       ("HG2", "H"),
                       ("CD", "C"),
                       ("HD1", "H"),
                       ("HD2", "H"),
                       ("CE", "C"),
                       ("HE1", "H"),
                       ("HE2", "H"),
                       ("NZ", "N"),
                       ("HZ1", "H"),
                       ("HZ2", "H"),
                       ("HZ3", "H"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "HG1"),
                       ("CG", "HG2"),
                       ("CG", "CD"),
                       ("CD", "HD1"),
                       ("CD", "HD2"),
                       ("CD", "CE"),
                       ("CE", "HE1"),
                       ("CE", "HE2"),
                       ("CE", "NZ"),
                       ("NZ", "HZ1"),
                       ("NZ", "HZ2"),
                       ("NZ", "HZ3"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "CD"),
    chi3_definition = ("CB", "CG", "CD", "CE"),
    chi4_definition = ("CG", "CD", "CE", "NZ"))
### </LYSINE>


### <ARGININE>
ARG = AminoAcid(
    name            = "ARG",
    full_name       = "Arginine",
    one_letter_name = "R",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("HG1", "H"),
                       ("HG2", "H"),
                       ("CD", "C"),
                       ("HD1", "H"),
                       ("HD2", "H"),
                       ("NE", "N"),
                       ("HE", "H"),
                       ("CZ", "C"),
                       ("NH1", "N"),
                       ("HH11", "H"),
                       ("HH12", "H"),
                       ("NH2", "N"),
                       ("HH21", "H"),
                       ("HH22", "H"),
                       ("C", "C"),
                       ("O", "O")],

    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "HG1"),
                       ("CG", "HG2"),
                       ("CG", "CD"),
                       ("CD", "HD1"),
                       ("CD", "HD2"),
                       ("CD", "NE"),
                       ("NE", "HE"),
                       ("NE", "CZ"),
                       ("CZ", "NH1"),
                       ("NH1", "HH11"),
                       ("NH1", "HH12"),
                       ("CZ", "NH2"),
                       ("NH2", "HH21"),
                       ("NH2", "HH22"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "CD"),
    chi3_definition = ("CB", "CG", "CD", "NE"),
    chi4_definition = ("CG", "CD", "NE", "CZ"))
### </ARGININE>


### <ASPATATE>
ASP = AminoAcid(
    name            = "ASP",
    full_name       = "Aspatate",
    one_letter_name = "D",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("OD1", "O"),
                       ("OD2", "O"),
                       ("C", "C"),
                       ("O", "O")],

    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "OD1"),
                       ("CG", "OD2"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "OD1"))
### </ASPATATE>


### <GLUTAMATE>
GLU = AminoAcid(
    name            = "GLU",
    full_name       = "Glutamate",
    one_letter_name = "E",

    atom_list       = [("N", "N"),
                       ("H", "H"),
                       ("CA", "C"),
                       ("HA", "H"),
                       ("CB", "C"),
                       ("HB1", "H"),
                       ("HB2", "H"),
                       ("CG", "C"),
                       ("HG1", "H"),
                       ("HG2", "H"),
                       ("CD", "C"),
                       ("OE1", "O"),
                       ("OE2", "O"),
                       ("C", "C"),
                       ("O", "O")],
    
    bond_list       = [("N", "H"),
                       ("N", "CA"),
                       ("CA", "HA"),
                       ("CA", "CB"),
                       ("CB", "HB1"),
                       ("CB", "HB2"),
                       ("CB", "CG"),
                       ("CG", "HG1"),
                       ("CG", "HG2"),
                       ("CG", "CD"),
                       ("CD", "OE1"),
                       ("CD", "OE2"),
                       ("CA", "C"),
                       ("C", "O")],
    
    chi1_definition = ("N", "CA", "CB", "CG"),
    chi2_definition = ("CA", "CB", "CG", "CD"),
    chi3_definition = ("CA", "CB", "CG", "OE1"))
### </GLUTAMATE>


## for accessing the AminoAcid classes
AminoAcidMap = {
    "ALA" : ALA,
    "VAL" : VAL,
    "LEU" : LEU,
    "ILE" : ILE,
    "GLY" : GLY,
    "PRO" : PRO,
    "CYS" : CYS,
    "MET" : MET,
    "HIS" : HIS,
    "PHE" : PHE,
    "TYR" : TYR,
    "TRP" : TRP,
    "ASN" : ASN,
    "GLN" : GLN,
    "SER" : SER,
    "THR" : THR,
    "LYS" : LYS,
    "ARG" : ARG,
    "ASP" : ASP,
    "GLU" : GLU,
    }
