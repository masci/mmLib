## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

"""This module provides a very minimal library of the amino acids.  The
amino acid properties defined here should be properties which will never
change, and more complex properties should be loaded by other libraries."""



AminoAcidNames = [
    "ALA", "VAL", "LEU", "ILE", "GLY", "PRO", "CYS", "MET", "HIS",
    "PHE", "TYR", "TRP", "ASN", "GLN", "SER", "THR", "LYS", "ARG",
    "ASP", "GLU", "DIS", "UNK" ]


class AminoAcid:
    """Empty definition class for building a Peptide library."""
    def __init__(self):
        self.name = ""
        self.full_name = ""
        self.one_letter_name = ""
        self.chi1_definition = None
        self.chi2_definition = None
        self.chi3_definition = None
        self.chi4_definition = None
        self.penultimate_rotamer_list = []

    def __str__(self):
	return "AminoAcid::%s" % (self.name)


## Unknown
UNK = AminoAcid()
UNK.name = "UNK"
UNK.full_name = "Unknown"
UNK.one_letter_name = "U"
UNK.atom_list = [("N", "N"),
                 ("H", "H"),
                 ("CA", "C"),
                 ("HA", "H"),
                 ("C", "C"),
                 ("O", "O")]
UNK.bond_list = [("N", "H"),
                 ("N", "CA"),
                 ("CA", "HA"),
                 ("CA", "C"),
                 ("C", "O")]

## Alanine
ALA = AminoAcid()
ALA.name = "ALA"
ALA.full_name = "Alanine"
ALA.one_letter_name = "A"
ALA.atom_list = [("N", "N"),
                 ("H", "H"),
                 ("CA", "C"),
                 ("HA", "H"),
                 ("CB", "C"),
                 ("HB1", "H"),
                 ("HB2", "H"),
                 ("HB3", "H"),
                 ("C", "C"),
                 ("O", "O")]
ALA.bond_list = [("N", "H"),
                 ("N", "CA"),
                 ("CA", "HA"),
                 ("CA", "CB"),
                 ("CB", "HB1"),
                 ("CB", "HB2"),
                 ("CB", "HB3"),
                 ("CA", "C"),
                 ("C", "O")]
    
## Valine
VAL = AminoAcid()
VAL.name = "VAL"
VAL.full_name = "Valine"
VAL.one_letter_name = "V"
VAL.atom_list = [("N", "N"),
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
                 ("O", "O")]
VAL.bond_list = [("N", "H"),
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
                 ("C", "O")]
VAL.chi1_definition = ("N", "CA", "CB", "CG1") #CG2?

## Leusine
LEU = AminoAcid()
LEU.name = "LEU"
LEU.full_name = "Leusine"
LEU.one_letter_name = "L"
LEU.atom_list = [("N", "N"),
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
                 ("O", "O")]
LEU.bond_list = [("N", "H"),
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
                 ("C", "O")]
LEU.bond_definition = ( )
LEU.chi1_definition = ("N", "CA", "CB", "CG")
LEU.chi2_definition = ("CA", "CB", "CG", "CD1") #CD2?

## Isoleusine
ILE = AminoAcid()
ILE.name = "ILE"
ILE.full_name = "Isoleusine"
ILE.one_letter_name = "I"
ILE.atom_list = [("N", "N"),
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
                 ("O", "O")]
ILE.bond_list = [("N", "H"),
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
                 ("C", "O")]
ILE.bond_definition = ( )
ILE.chi1_definition = ("N", "CA", "CB", "CG1") #CG2?
ILE.chi2_definition = ("CA", "CB", "CG1", "CD1")

## Glysine
GLY = AminoAcid()
GLY.name = "GLY"
GLY.full_name = "Glysine"
GLY.one_letter_name = "G"
GLY.atom_list = [("N", "N"),
                 ("H", "H"),
                 ("CA", "C"),
                 ("HA1", "H"),
                 ("HA2", "H"),
                 ("C", "C"),
                 ("O", "O")]
GLY.bond_list = [("N", "H"),
                 ("N", "CA"),
                 ("CA", "HA1"),
                 ("CA", "HA2"),
                 ("CA", "C"),
                 ("C", "O")]

## Proline
PRO = AminoAcid()
PRO.name = "PRO"
PRO.full_name = "Proline"
PRO.one_letter_name = "P"
PRO.atom_list = [("N", "N"),
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
                 ("O", "O")]
PRO.bond_list = [("N", "CA"),
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
                 ("C", "O")]
PRO.chi1_definition = ("N", "CA", "CB", "CG")
PRO.pucker_definition = ("C", "CA", "CB", "CG")

## Cystine
CYS = AminoAcid()
CYS.name = "CYS"
CYS.full_name = "Cystine"
CYS.one_letter_name = "C"
CYS.atom_list = [("N", "N"),
                 ("H", "H"),
                 ("CA", "C"),
                 ("HA", "H"),
                 ("CB", "C"),
                 ("HB1", "H"),
                 ("HB2", "H"),
                 ("SG", "S"),
                 ("C", "C"),
                 ("O", "O")]
CYS.bond_list = [("N", "H"),
                 ("N", "CA"),
                 ("CA", "HA"),
                 ("CA", "CB"),
                 ("CB", "HB1"),
                 ("CB", "HB2"),
                 ("CB", "SG"),
                 ("CA", "C"),
                 ("C", "O")]
CYS.chi1_definition = ("N", "CA", "CB", "SG")
#CYS.chi3_definition = ("CB", "SG", "SG", "CB") disulfur bond

## Methionine
MET = AminoAcid()
MET.name = "MET"
MET.full_name = "Methionine"
MET.one_letter_name = "M"
MET.atom_list = [("N", "N"),
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
                 ("O", "O")]
MET.bond_list = [("N", "H"),
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
                 ("C", "O")]
MET.chi1_definition = ("N", "CA", "CB", "CG")
MET.chi2_definition = ("CA", "CB", "CG", "SD")
MET.chi3_definition = ("CB", "CG", "SD", "CE")
MET.chi4_definition = None

## Histidine
HIS = AminoAcid()
HIS.name = "HIS"
HIS.full_name = "Histidine"
HIS.one_letter_name = "H"
HIS.atom_list = [("N", "N"),
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
                 ("O", "O")]
HIS.bond_list = [("N", "H"),
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
                 ("C", "O")]
HIS.chi1_definition = ("N", "CA", "CB", "CG")
HIS.chi2_definition = ("CA", "CB", "CG", "ND1")

## Phenylalaline
PHE = AminoAcid()
PHE.name = "PHE"
PHE.full_name = "Phenylalaline"
PHE.one_letter_name = "F"
PHE.atom_list = [("N", "N"),
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
                 ("O", "O")]
PHE.bond_list = [("N", "H"),
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
                 ("C", "O")]
PHE.chi1_definition = ("N", "CA", "CB", "CG")
PHE.chi2_definition = ("CA", "CB", "CG", "CD1")

## Tyrosine
TYR = AminoAcid()
TYR.name = "TYR"
TYR.full_name = "Tyrosine"
TYR.one_letter_name = "Y"
TYR.atom_list = [("N", "N"),
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
                 ("O", "O")]
TYR.bond_list = [("N", "H"),
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
                 ("C", "O")]
TYR.chi1_definition = ("N", "CA", "CB", "CG")
TYR.chi2_definition = ("CA", "CB", "CG", "CD1")

## Trptophan
TRP = AminoAcid()
TRP.name = "TRP"
TRP.full_name = "Trptophan"
TRP.one_letter_name = "W"
TRP.atom_list = [("N", "N"),
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
                 ("O", "O")]
TRP.bond_list = [("N", "H"),
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
                 ("C", "O")]
TRP.chi1_definition = ("N", "CA", "CB", "CG")
TRP.chi2_definition = ("CA", "CB", "CG", "CD1")

## Asparagine
ASN = AminoAcid()
ASN.name = "ASN"
ASN.full_name = "Asparagine"
ASN.one_letter_name = "N"
ASN.atom_list = [("N", "N"),
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
                 ("O", "O")]
ASN.bond_list = [("N", "H"),
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
                 ("C", "O")]
ASN.chi1_definition = ("N", "CA", "CB", "CG")
ASN.chi2_definition = ("CA", "CB", "CG", "OD1")

## Glutamine
GLN = AminoAcid()
GLN.name = "GLN"
GLN.full_name = "Glutamine"
GLN.one_letter_name = "Q"
GLN.atom_list = [("N", "N"),
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
                 ("O", "O")]
GLN.bond_list = [("N", "H"),
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
                 ("C", "O")]
GLN.chi1_definition = ("N", "CA", "CB", "CG")
GLN.chi2_definition = ("CA", "CB", "CG", "CD")
GLN.chi3_definition = ("CB", "CG", "CD", "OE1")

## Serine
SER = AminoAcid()
SER.name = "SER"
SER.full_name = "Serine"
SER.one_letter_name = "S"
SER.atom_list = [("N", "N"),
                 ("H", "H"),
                 ("CA", "C"),
                 ("HA", "H"),
                 ("CB", "C"),
                 ("HB1", "H"),
                 ("HB2", "H"),
                 ("OG", "O"),
                 ("HG", "H"),
                 ("C", "C"),
                 ("O", "O")]
SER.bond_list = [("N", "H"),
                 ("N", "CA"),
                 ("CA", "HA"),
                 ("CA", "CB"),
                 ("CB", "HB1"),
                 ("CB", "HB2"),
                 ("CB", "OG"),
                 ("OG", "HG"),
                 ("CA", "C"),
                 ("C", "O")]
SER.chi1_definition = ("N", "CA", "CB", "OG")

## Threonine
THR = AminoAcid()
THR.name = "THR"
THR.full_name = "Threonine"
THR.one_letter_name = "T"
THR.atom_list = [("N", "N"),
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
                 ("O", "O")]
THR.bond_list = [("N", "H"),
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
                 ("C", "O")]
THR.chi1_definition = ("N", "CA", "CB", "OG1")


## Lysine
LYS = AminoAcid()
LYS.name = "LYS"
LYS.full_name = "lysine"
LYS.one_letter_name = "K"
LYS.atom_list = [("N", "N"),
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
                 ("O", "O")]
LYS.bond_list = [("N", "H"),
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
                 ("C", "O")]
LYS.chi1_definition = ("N", "CA", "CB", "CG")
LYS.chi2_definition = ("CA", "CB", "CG", "CD")
LYS.chi3_definition = ("CB", "CG", "CD", "CE")
LYS.chi4_definition = ("CG", "CD", "CE", "NZ")

## Arginine
ARG = AminoAcid()
ARG.name = "ARG"
ARG.full_name = "Arginine"
ARG.one_letter_name = "R"
ARG.atom_list = [("N", "N"),
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
                 ("O", "O")]
ARG.bond_list = [("N", "H"),
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
                 ("C", "O")]
ARG.chi1_definition = ("N", "CA", "CB", "CG")
ARG.chi2_definition = ("CA", "CB", "CG", "CD")
ARG.chi3_definition = ("CB", "CG", "CD", "NE")
ARG.chi4_definition = ("CG", "CD", "NE", "CZ")

## Aspatate
ASP = AminoAcid()
ASP.name = "ASP"
ASP.full_name = "Aspatate"
ASP.one_letter_name = "D"
ASP.atom_list = [("N", "N"),
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
                 ("O", "O")]
ASP.bond_list = [("N", "H"),
                 ("N", "CA"),
                 ("CA", "HA"),
                 ("CA", "CB"),
                 ("CB", "HB1"),
                 ("CB", "HB2"),
                 ("CB", "CG"),
                 ("CG", "OD1"),
                 ("CG", "OD2"),
                 ("CA", "C"),
                 ("C", "O")]
ASP.chi1_definition = ("N", "CA", "CB", "CG")
ASP.chi2_definition = ("CA", "CB", "CG", "OD1")

## Glutamate
GLU = AminoAcid()
GLU.name = "GLU"
GLU.full_name = "Glutamate"
GLU.one_letter_name = "E"
GLU.atom_list = [("N", "N"),
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
                 ("O", "O")]
GLU.bond_list = [("N", "H"),
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
                 ("C", "O")]
GLU.chi1_definition = ("N", "CA", "CB", "CG")
GLU.chi2_definition = ("CA", "CB", "CG", "CD")
GLU.chi3_definition = ("CA", "CB", "CG", "OE1")

## Disulfur Bond
DIS = AminoAcid()
DIS.name = "DIS"
DIS.full_name = "Disulfur Bond"
DIS.one_letter_name = None
DIS.chi1_definition = ("CA", "CB", "SG", "SG")
DIS.chi3_definition = ("CB", "SG", "SG", "CB")
DIS.chi3_definition = ("SG", "SG", "CB", "CA")


## for accessing the AminoAcid classes
AminoAcidMap = {
    "UNK" : UNK,
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
    "DIS" : DIS,
    }
