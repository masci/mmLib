## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

"""This module provides a PDB v2.2 parser.  All records in the PDB v2.2
specification have coorasponding classes defined here.  PDB files are
loaded into a list of these cassed, and also can be constrcted/modified
and written back out as PDB files."""



import string
import types
import fpformat
from   FileIO   import OpenFile


PDBError = "PDB Error"



## these are the fields in a ATOM/HETATM/ANISOU/etc... record that uniquely
## identify a record; often passed into compareRecord when searching for
## additional PDBRecord classes describing a ATOM
AtomIDFields = ["name", "altLoc", "resName", "seqRes", "iCode", "chainID"]


## these are the fields uniquely identifing a specic residue
ResidueIDFields = ["resName", "seqRes", "iCode", "chainID"]



class PDBRecord:
    """Base class for all PDB file records."""

    def __str__(self):
        return self.write()


    def compareRecord(self, rec, field_list):
        """Compares record rec to self for the fields listed in field_list."""
        for field in field_list:
            try:
                val1 = getattr(self, field)
            except AttributeError:
                val1 = ""

            try:
                val2 = getattr(rec, field)
            except AttributeError:
                val2 = ""

            if val1 != val2:
                return 0

        return 1


    def write(self):
        rec = self._name
        for (field, start, end, ftype, just) in self._field_list:
            if len(rec) > start:
                print "[ERROR] len(rec)=%d start=%d" % (i, start)

            ## add spaces to the end if necessary
            rec = rec.ljust(start - 1)
                
            ## access the namespace of this class to write the field
            try:
                s = getattr(self, field)
            except AttributeError:
                s = ""

            ## convert integer and float types
            if ftype == "string":
                pass
                
            elif ftype == "integer":
                try:
                    s = str(s)
                except ValueError:
                    raise PDBError, "field=%s %s not int" % (field, s)
                
            elif ftype[:5] == "float":
                d = int(ftype[6])
                try:
                    s = fpformat.fix(s, d)
                except ValueError:
                    raise PDBError, "field=%s %s not float" % (field, s)

            else:
                raise PDBError, "INVALID TYPE: %s" % (ftype)

            ## check for maximum length
            l = end - start + 1
            if len(s) > l:
                raise PDBError, "field=%s value=%s len=%d > max=%d" % (
                    field, s, len(s), l)

            if just == "ljust":
                s = s.ljust(l)
            else:
                s = s.rjust(l)

            rec += s.upper()

        return rec
    

    def read(self, rec):
        for (field, start, end, ftype, just) in self._field_list:
            ## adjust record reading indexes if the line doesn't contain
            ## all the fields
            if end > len(rec):
                if start > len(rec):
                    break
                end = len(rec)

            s = rec[start-1:end].strip()
            if not s:
                continue

            elif ftype == "string":
                pass

            elif ftype == "integer":
                try:
                    s = int(s)
                except ValueError:
                    print rec
                    raise PDBError, "int(\"%s\") failed" % (s)

            elif ftype[:5] == "float":
                try:
                    s = float(s)
                except ValueError:
                    print rec
                    raise PDBError, "float conversion error"

            setattr(self, field, s)


###############################################################################
## BEGIN PDB RECORD DEFINITIONS

## SECTION 2: Title Section
class HEADER(PDBRecord):
    """This section contains records used to describe the experiment and the
    biological macromolecules present in the entry: HEADER, OBSLTE, TITLE,
    CAVEAT, COMPND, SOURCE, KEYWDS, EXPDTA, AUTHOR, REVDAT, SPRSDE, JRNL,
    and REMARK records."""
    _name = "HEADER"
    _field_list = [
        ("classification", 11, 50, "string", "rjust"),
        ("depDate", 51, 59, "string", "rjust"),
        ("idCode", 63, 66, "string", "rjust")]

class OBSLTE(PDBRecord):
    """OBSLTE appears in entries which have been withdrawn from distribution.
    This record acts as a flag in an entry which has been withdrawn from the
    PDB's full release. It indicates which, if any, new entries have replaced
    the withdrawn entry.  The format allows for the case of multiple new
    entries replacing one existing entry."""
    _name = "OBSLTE"
    _field_list = [
        ("continuation", 9, 10, "string", "rjust"),
        ("repDate", 12, 20, "string", "rjust"),
        ("idCode", 22, 25, "string", "rjust"),
        ("rIdCode1", 32, 35, "string", "rjust"),
        ("rIdCode2", 37, 40, "string", "rjust"),
        ("rIdCode3", 42, 45, "string", "rjust"),
        ("rIdCode4", 47, 50, "string", "rjust"),
        ("rIdCode5", 52, 55, "string", "rjust"),
        ("rIdCode6", 57, 60, "string", "rjust"),
        ("rIdCode7", 62, 65, "string", "rjust"),
        ("rIdCode8", 67, 70, "string", "rjust")]

class TITLE(PDBRecord):
    """The TITLE record contains a title for the experiment or analysis that is
    represented in the entry. It should identify an entry in the PDB in the
    same way that a title identifies a paper."""
    _name = "TITLE "
    _field_list = [
        ("continuation", 9, 10, "string", "rjust"),
        ("title", 11, 70, "string", "ljust")]

class CAVEAT(PDBRecord):
    """CAVEAT warns of severe errors in an entry. Use caution when using an
    entry containing this record."""
    _name = "CAVEAT"
    _field_list = [
        ("continuation", 9, 10, "string", "rjust"),
        ("idCode", 12, 15, "string", "rjust"),
        ("comment", 20, 70, "string", "ljust")]

class COMPND(PDBRecord):
    """The COMPND record describes the macromolecular contents of an entry.
    Each macromolecule found in the entry is described by a set of token: value
    pairs, and is referred to as a COMPND record component. Since the concept
    of a molecule is difficult to specify exactly, PDB staff may exercise
    editorial judgment in consultation with depositors in assigning these
    names.  For each macromolecular component, the molecule name, synonyms,
    number assigned by the Enzyme Commission (EC), and other relevant details
    are specified. """ 
    _name = "COMPND"
    _field_list = [
        ("continuation", 9, 10, "string", "rjust"),
        ("compound", 11, 70, "string", "ljust")]

class SOURCE(PDBRecord):
    """The SOURCE record specifies the biological and/or chemical source of
    each biological molecule in the entry. Sources are described by both the
    common name and the scientific name, e.g., genus and species. Strain and/or
    cell-line for immortalized cells are given when they help to uniquely
    identify the biological entity studied."""
    _name = "SOURCE"
    _field_list = [
        ("continuation", 9, 10, "string", "rjust"),
        ("srcName", 11, 70, "string", "ljust")]

class KEYWDS(PDBRecord):
    """The KEYWDS record contains a set of terms relevant to the entry. Terms
    in the KEYWDS record provide a simple means of categorizing entries and may
    be used to generate index files. This record addresses some of the
    limitations found in the classification field of the HEADER record. It
    provides the opportunity to add further annotation to the entry in a concise
    and computer-searchable fashion."""
    _name = "KEYWDS"
    _field_list = [
        ("continuation", 9, 10, "string", "rjust"),
        ("keywds", 11, 70, "string", "ljust")]

class EXPDTA(PDBRecord):
    """The EXPDTA record presents information about the experiment.  The EXPDTA
    record identifies the experimental technique used. This may refer to the
    type of radiation and sample, or include the spectroscopic or modeling
    technique. Permitted values include: 
    ELECTRON DIFFRACTION
    FIBER DIFFRACTION
    FLUORESCENCE TRANSFER
    NEUTRON DIFFRACTION
    NMR
    THEORETICAL MODEL
    X-RAY DIFFRACTION""" 
    _name = "EXPDTA"
    _field_list = [
        ("continuation", 9, 10, "string", "rjust"),
        ("technique", 11, 70, "string", "ljust")]

class AUTHOR(PDBRecord):
    """The AUTHOR record contains the names of the people responsible for the
    contents of the entry."""
    _name = "AUTHOR"
    _field_list = [
        ("continuation", 9, 10, "string", "rjust"),
        ("authorList", 11, 70, "string", "ljust")]

class REVDAT(PDBRecord):
    """REVDAT records contain a history of the modifications made to an entry
    since its release."""
    _name = "REVDAT"
    _field_list = [
        ("modNum", 8, 10, "integer", "rjust"),
        ("continuation", 11, 12, "string", "rjust"),
        ("modDate", 14, 22, "string", "rjust"),
        ("modID", 24, 28, "string", "rjust"),
        ("modType", 32, 32, "integer", "rjust"),
        ("record1", 40, 45, "string", "ljust"),
        ("record2", 47, 52, "string", "ljust"),
        ("record3", 54, 59, "string", "ljust"),
        ("record4", 61, 66, "string", "ljust")]

class SPRSDE(PDBRecord):
    """The SPRSDE records contain a list of the ID codes of entries that were
    made obsolete by the given coordinate entry and withdrawn from the PDB
    release set. One entry may replace many. It is PDB policy that only the
    principal investigator of a structure has the authority to withdraw it."""
    _name = "SPRSDE"
    _field_list = [
        ("continuation", 9, 10, "string", "rjust"),
        ("sprsdeDate", 12, 20, "string", "rjust"),
        ("idCode", 22, 25, "string", "rjust"),
        ("sIdCode1", 32, 35, "string", "rjust"),
        ("sIdCode2", 37, 40, "string", "rjust"),
        ("sIdCode3", 42, 45, "string", "rjust"),
        ("sIdCode4", 47, 50, "string", "rjust"),
        ("sIdCode5", 52, 55, "string", "rjust"),
        ("sIdCode6", 57, 60, "string", "rjust"),
        ("sIdCode7", 62, 65, "string", "rjust"),
        ("sIdCode8", 67, 70, "string", "rjust")]

class JRNL(PDBRecord):
    """The JRNL record contains the primary literature citation that describes
    the experiment which resulted in the deposited coordinate set. There is at
    most one JRNL reference per entry. If there is no primary reference, then
    there is no JRNL reference. Other references are given in REMARK 1."""
    _name = "JRNL  "
    _field_list = [
        ("text", 13, 70, "string", "ljust")]

class REMARK(PDBRecord):
    """REMARK records present experimental details, annotations, comments, and
    information not included in other records. In a number of cases, REMARKs are
    used to expand the contents of other record types. A new level of structure
    is being used for some REMARK records. This is expected to facilitate
    searching and will assist in the conversion to a relational database."""
    _name = "REMARK"
    _field_list = [
        ("remarkNum", 8, 10, "integer", "rjust"),
        ("text", 12, 70, "string", "ljust")]
        
## SECTION 3: Primary Structure Section
class DBREF(PDBRecord):
    """ The DBREF record provides cross-reference links between PDB sequences
    and the corresponding database entry or entries. A cross reference to
    the sequence database is mandatory for each peptide chain with a length
    greater than ten (10) residues. For nucleic acid entries a DBREF
    record pointing to the Nucleic Acid Database (NDB) is mandatory when
    the corresponding entry exists in NDB."""
    _name = "DBREF "
    _field_list = [
        ("idCode", 8, 11, "string", "rjust"),
        ("chain_ID", 13, 13, "string", "rjust"),
        ("seqBegin", 15, 18, "integer", "rjust"),
        ("insertBegin", 19, 19, "string", "rjust"),
        ("seqEnd", 21, 24, "integer", "rjust"),
        ("insertEnd", 25, 25, "string", "rjust"),
        ("database", 27, 32, "string", "ljust"),
        ("dbAccession", 34, 41, "string", "ljust"),
        ("dbIdCode", 43, 54, "string", "ljust"),
        ("dbseqBegin", 56, 60, "integer", "rjust"),
        ("idbnsBeg", 61, 61, "string", "rjust"),
        ("dbseqEnd", 63, 67, "integer", "rjust"),
        ("dbinsEnd", 68, 68, "string", "rjust")]

class SEQADV(PDBRecord):
    """The SEQADV record identifies conflicts between sequence information
    in the ATOM records of the PDB entry and the sequence database entry
    given on DBREF. Please note that these records were designed to
    identify differences and not errors. No assumption is made as to which
    database contains the correct data. PDB may include REMARK records in
    the entry that reflect the depositor's view of which database has the
    correct sequence."""
    _name = "SEQADV"
    _field_list = [
        ("idCode", 8, 11, "string", "rjust"),
        ("resName", 13, 15, "string", "rjust"),
        ("chainID", 17, 17, "string", "rjust"),
        ("seqNum", 19, 22, "integer", "rjust"),
        ("iCode", 23, 23, "string", "rjust"),
        ("database", 25, 28, "string", "ljust"),
        ("dbIDCode", 30, 38, "string", "ljust"),
        ("dbRes", 40, 42, "string", "rjust"),
        ("dbSeq", 44, 48, "integer", "rjust"),
        ("convlict", 40, 70, "string", "ljust")]
    
class SEQRES(PDBRecord):
    """The SEQRES records contain the amino acid or nucleic acid sequence of
    residues in each chain of the macromolecule that was studied."""
    _name = "SEQRES"
    _field_list = [
        ("serNum", 9, 10, "integer", "rjust"),
        ("chainID", 12, 12, "string", "rjust"),
        ("numRes", 14, 17, "integer", "rjust"),
        ("resName1", 20, 22, "string", "rjust"),
        ("resName2", 24, 26, "string", "rjust"),
        ("resName3", 28, 30, "string", "rjust"),
        ("resName4", 32, 34, "string", "rjust"),
        ("resName5", 36, 38, "string", "rjust"),
        ("resName6", 40, 42, "string", "rjust"),
        ("resName7", 44, 46, "string", "rjust"),
        ("resName8", 48, 50, "string", "rjust"),
        ("resName9", 52, 54, "string", "rjust"),
        ("resName10", 56, 58, "string", "rjust"),
        ("resName11", 60, 62, "string", "rjust"),
        ("resName12", 64, 66, "string", "rjust"),
        ("resName13", 68, 70, "string", "rjust")]

class MODRES(PDBRecord):
    """The MODRES record provides descriptions of modifications (e.g.,
    chemical or post-translational) to protein and nucleic acid residues.
    Included are a mapping between residue names given in a PDB entry and
    standard residues."""
    _name = "MODRES"
    _field_list = [
        ("idCode", 8, 11, "string", "rjust"),
        ("resName", 13, 15, "string", "rjust"),
        ("chainID", 17, 17, "string", "rjust"),
        ("seqNum", 19, 22, "integer", "rjust"),
        ("iCode", 23, 23, "string", "rjust"),
        ("stdRes", 25, 27, "string", "rjust"),
        ("comment", 30, 70, "string", "ljust")]

## SECTION 4: Heterogen Section
class HET(PDBRecord):
    """The HET records are used to describe non-standard residues, such as
    prosthetic groups, inhibitors, solvent molecules, and ions for
    which coordinates are supplied. Groups are considered HET if they are: 
    - not one of the standard amino acids, and 
    - not one of the nucleic acids (C, G, A, T, U, and I), and 
    - not one of the modified versions of nucleic acids (+C, +G, +A,
      +T, +U, and +I), and 
    - not an unknown amino acid or nucleic acid where UNK is used to
      indicate the unknown residue name. 
    Het records also describe heterogens for which the chemical identity
    is unknown, in which case the group is assigned the hetID UNK."""
    _name = "HET   "
    _field_list = [
        ("hetID", 8, 10, "string", "rjust"),
        ("chainID", 13, 13, "string", "rjust"),
        ("seqNum", 14, 17, "integer", "rjust"),
        ("iCode", 18, 18, "string", "rjust"),
        ("numHetAtoms", 21, 25, "integer", "rjust"),
        ("text", 31, 70, "string", "ljust")]
    
class HETNAM(PDBRecord):
    """This record gives the chemical name of the compound with the
    given hetID."""
    _name = "HETNAM"
    _field_list = [
        ("continuation", 9, 10, "string", "ljust"),
        ("hetID", 12, 14, "string", "rjust"),
        ("text", 16, 70, "string", "ljust")]
    
class HETSYN(PDBRecord):
    """This record provides synonyms, if any, for the compound in the
    corresponding (i.e., same hetID) HETNAM record. This is to allow
    greater flexibility in searching for HET groups."""
    _name = "HETSYN"
    _field_list = [
        ("continuation", 9, 10, "string", "ljust"),
        ("hetID", 12, 14, "string", "rjust"),
        ("hetSynonyms", 16, 70, "string", "ljust")]

class FORMUL(PDBRecord):
    """The FORMUL record presents the chemical formula and charge of a
    non-standard group. (The formulas for the standard residues are given
    in Appendix 5.)"""
    _name = "FORMUL"
    _field_list = [
        ("compNum", 9, 10, "integer", "rjust"),
        ("hetID", 13, 15, "string", "rjust"),
        ("continuation", 17, 18, "integer", "rjust"),
        ("asterisk", 19, 19, "string", "rjust"),
        ("text", 20, 70, "string", "ljust")]

## SECTION 5: Secondary Structure Section
class HELIX(PDBRecord):
    """HELIX records are used to identify the position of helices in the
    molecule. Helices are both named and numbered. The residues where the
    helix begins and ends are noted, as well as the total length."""
    _name = "HELIX "
    _field_list = [
        ("serNum", 8, 10, "integer", "rjust"),
        ("helixID", 12, 14, "string", "rjust"),
        ("initResName", 16, 18, "string", "rjust"),
        ("initChainID", 20, 20, "string", "rjust"),
        ("initSeqNum", 22, 25, "integer", "rjust"),
        ("initICode", 26, 26, "string", "rjust"),
        ("endResName", 28, 30, "string", "rjust"),
        ("endChainID", 32, 32, "string", "rjust"),
        ("endSeqNum", 34, 37, "integer", "rjust"),
        ("endICode", 38, 38, "string", "rjust"),
        ("helixClass", 39, 40, "integer", "rjust"),
        ("comment", 41, 70, "string", "ljust"),
        ("length", 72, 76, "integer", "rjust")]

class SHEET(PDBRecord):
    """SHEET records are used to identify the position of sheets in the
    molecule. Sheets are both named and numbered. The residues where the
    sheet begins and ends are noted."""
    _name = "SHEET "
    _field_list = [
        ("strand", 8, 10, "integer", "rjust"),
        ("sheetID", 12, 14, "string", "rjust"),
        ("numStrands", 15, 16, "integer", "rjust"),
        ("initResName", 18, 20, "string", "rjust"),
        ("initChainID", 22, 22, "string", "rjust"),
        ("initSeqNum", 23, 26, "integer", "rjust"),
        ("initICode", 27, 27, "string", "rjust"),
        ("endResName", 29, 31, "string", "rjust"),
        ("endChainID", 33, 33, "string", "rjust"),
        ("endSeqNum", 34, 37, "integer", "rjust"),
        ("endICode", 38, 38, "string", "rjust"),
        ("sense", 39, 40, "integer", "rjust"),
        ("curAtom", 42, 45, "string", "rjust"),
        ("curResName", 46, 48, "string", "rjust"),
        ("curChainID", 50 ,50, "string", "rjust"),
        ("curResSeq", 51, 54, "integer", "rjust"),
        ("curICode", 55, 55, "string", "rjust"),
        ("prevAtom", 57, 60, "string", "rjust"),
        ("prevResName", 61, 63, "string", "rjust"),
        ("prevChainID", 65, 65, "string", "rjust"),
        ("prevResSeq", 66, 69, "integer", "rjust"),
        ("prevICode", 70, 70, "string", "rjust")]

class TURN(PDBRecord):
    """The TURN records identify turns and other short loop turns which
    normally connect other secondary structure segments."""
    _name = "TURN  "
    _field_list = [
        ("seq", 8, 10, "integer", "rjust"),
        ("turnID", 12, 14, "string", "rjust"),
        ("initResName", 16, 18, "string", "rjust"),
        ("initChainID", 20, 20, "string", "rjust"),
        ("initSeqNum", 21, 24, "integer", "rjust"),
        ("initICode", 25, 25, "string", "rjust"),
        ("endResName", 27, 29, "string", "rjust"),
        ("endChainID", 31, 31, "string", "rjust"),
        ("endSeqNum", 32, 35, "integer", "rjust"),
        ("endICode", 36, 36, "string", "rjust"),
        ("comment", 41, 70, "string", "ljust")]

## SECTION 6: Connectivity Annotation Section
class SSBOND(PDBRecord):
    """The SSBOND record identifies each disulfide bond in protein and
    polypeptide structures by identifying the two residues involved in the
    bond."""
    _name = "SSBOND"
    _field_list = [
        ("serNum", 8, 10, "integer", "rjust"),
        ("resName1", 12, 14, "string", "rjust"),
        ("chainID1", 16, 16, "string", "rjust"),
        ("seqNum1", 18, 21, "integer", "rjust"),
        ("iCode1", 22, 22, "string", "rjust"),
        ("resName2", 26, 28, "string", "rjust"),
        ("chainID2", 30, 30, "string", "rjust"),
        ("seqNum2", 32, 35, "integer", "rjust"),
        ("iCode2", 36, 36, "string", "rjust"),
        ("sym1", 60, 65, "string", "rjust"),
        ("sym2", 67, 72, "string", "rjust")]

class LINK(PDBRecord):
    """The LINK records specify connectivity between residues that is not
    implied by the primary structure. Connectivity is expressed in terms of
    the atom names. This record supplements information given in CONECT
    records and is provided here for convenience in searching."""
    _name = "LINK  "
    _field_list = [
        ("name1", 13, 16, "string", "rjust"),
        ("altLoc1", 17, 17, "string", "rjust"),
        ("resName1", 18, 20, "string", "rjust"),
        ("chainID1", 22, 22, "string", "rjust"),
        ("resSeq1", 23, 26, "integer", "rjust"),
        ("iCode1", 27, 27, "string", "rjust"),
        ("name2", 43, 46, "string", "rjust"),
        ("altLoc2", 47, 47, "string", "rjust"),
        ("resName2", 48, 50, "string", "rjust"),
        ("chainID2", 52, 52, "string", "rjust"),
        ("resSeq2", 53, 56, "integer", "rjust"),
        ("iCode2", 57, 57, "string", "rjust"),
        ("sym1", 60, 65, "string", "rjust"),
        ("sym2", 67, 72, "string", "rjust")]

class HYDBND(PDBRecord):
    """The HYDBND records specify hydrogen bonds in the entry."""
    _name = "HYDBND"
    _field_list = [
        ("name1", 13, 16, "string", "rjust"),
        ("altLoc1", 17, 17, "string", "rjust"),
        ("resName1", 18, 20, "string", "rjust"),
        ("chainID1", 22, 22, "string", "rjust"),
        ("resSeq1", 23, 27, "integer", "rjust"),
        ("iCode1", 28, 28, "string", "rjust"),
        ("nameH", 30, 33, "string", "rjust"),
        ("altLocH", 34, 34, "string", "rjust"),
        ("chainH", 36, 36, "string", "rjust"),
        ("resSeqH", 37, 41, "integer", "rjust"),
        ("iCodeH", 42, 42, "string", "rjust"),
        ("name2", 44, 47, "string", "rjust"),
        ("altLoc2", 48, 48, "string", "rjust"),
        ("resName2", 49, 51, "string", "rjust"),
        ("chainID2", 53, 53, "string", "rjust"),
        ("resSeq2", 54, 58, "integer", "rjust"),
        ("iCode2", 59, 59, "string", "rjust"),
        ("sym1", 60, 65, "string", "rjust"),
        ("sym2", 67, 72, "string", "rjust")]

class SLTBRG(PDBRecord):
    """The SLTBRG records specify salt bridges in the entry."""
    _name = "SLTBRG"
    _field_list = [
        ("name1", 13, 16, "string", "rjust"),
        ("altLoc1", 17, 17, "string", "rjust"),
        ("resName1", 18, 20, "string", "rjust"),
        ("chainID1", 22, 22, "string", "rjust"),
        ("resSeq1", 23, 26, "integer", "rjust"),
        ("iCode1", 27, 27, "string", "rjust"),
        ("name2", 43, 46, "string", "rjust"),
        ("altLoc2", 47, 47, "string", "rjust"),
        ("resName2", 48, 50, "string", "rjust"),
        ("chainID2", 52, 52, "string", "rjust"),
        ("resSeq2", 53, 56, "integer", "rjust"),
        ("iCode2", 57, 57, "string", "rjust"),
        ("sym1", 60, 65, "string", "rjust"),
        ("sym2", 67, 72, "string", "rjust")]

class CISPEP(PDBRecord):
    """CISPEP records specify the prolines and other peptides found to be
    in the cis conformation. This record replaces the use of footnote records
    to list cis peptides."""
    _name = "CISPEP"
    _field_list = [
        ("serial", 8, 10, "integer", "rjust"),
        ("resName1", 12, 14, "string", "rjust"),
        ("chainID1", 16, 16, "string", "rjust"),
        ("seqNum1", 18, 21, "integer", "rjust"),
        ("iCode1", 22, 22, "string", "rjust"),
        ("resName2", 26, 28, "string", "rjust"),
        ("chainID2", 30, 30, "string", "rjust"),
        ("seqNum2", 32, 35, "integer", "rjust"),
        ("iCode2", 36, 36, "string", "rjust"),
        ("modNum", 44, 46, "integer", "rjust"),
        ("measure", 54, 59, "float.2", "rjust")]
    
## SECTION 7: Miscellanious Features Section
class SITE(PDBRecord):
    """The SITE records supply the identification of groups comprising
    important sites in the macromolecule."""
    _name = "SITE  "
    _field_list = [
        ("seqNum", 8, 10, "integer", "rjust"),
        ("siteID", 12, 14, "string", "rjust"),
        ("numRes", 16, 17, "integer", "rjust"),
        ("resName1", 19, 21, "string", "rjust"),
        ("chainID1", 23, 23, "string", "rjust"),
        ("seq1", 24, 27, "integer", "rjust"),
        ("iCode1", 28, 28, "string", "rjust"),
        ("resName2", 30, 32, "string", "rjust"),
        ("chainID2", 34, 34, "string", "rjust"),
        ("seq2", 35, 38, "integer", "rjust"),
        ("iCode2", 39, 39, "string", "rjust"),
        ("resName3", 41, 43, "string", "rjust"),
        ("chainID3", 45, 45, "string", "rjust"),
        ("seq3", 46, 49, "integer", "rjust"),
        ("iCode3", 50, 50, "string", "rjust"),
        ("resName4", 52, 54, "string", "rjust"),
        ("chainID4", 56, 56, "string", "rjust"),
        ("seq4", 57, 60, "integer", "rjust"),
        ("iCode4", 61, 61, "string", "rjust")]

## SECTION 8: Crystallographic and Coordinate Transformation Section
class CRYSTn(PDBRecord):
    """The CRYSTn (n=1,2,3) record presents the unit cell parameters, space
    group, and Z value. If the structure was not determined by crystallographic
    means, CRYSTn simply defines a unit cube."""
    _field_list = [
        ("a", 7, 15, "float.3", "rjust"),
        ("b", 16, 24, "float.3", "rjust"),
        ("c", 25, 33, "float.3", "rjust"),
        ("alpha", 34, 40, "float.3", "rjust"),
        ("beta", 41, 47, "float.3", "rjust"),
        ("gamma", 48, 54, "float.3", "rjust"),
        ("sgroup", 56, 66, "string", "rjust"),
        ("z", 67, 70, "integer", "rjust")]

class CRYST1(CRYSTn):
    _name = "CRYST1"
    
class CRYST2(CRYSTn):
    _name = "CRYST2"

class CRYST3(CRYSTn):
    _name = "CRYST3"

class ORIGXn(PDBRecord):
    """The ORIGXn (n = 1, 2, or 3) records present the transformation from
    the orthogonal coordinates contained in the entry to the submitted
    coordinates."""
    _field_list = [
        ("o[n][1]", 11, 20, "float.6", "rjust"),
        ("o[n][2]", 21, 30, "float.6", "rjust"),
        ("o[n][3]", 31, 40, "float.6", "rjust"),
        ("t[n]", 46, 55, "float.5", "rjust")]

class ORIGX1(ORIGXn):
    _name = "ORIGX1"

class ORIGX2(ORIGXn):
    _name = "ORIGX2"

class ORIGX3(ORIGXn):
    _name = "ORIGX3"
    
class SCALEn(PDBRecord):
    """The SCALEn (n = 1, 2, or 3) records present the transformation from
    the orthogonal coordinates as contained in the entry to fractional
    crystallographic coordinates. Non-standard coordinate systems should
    be explained in the remarks."""
    _field_list = [
        ("s[n][1]", 11, 20, "float.6", "rjust"),
        ("s[n][2]", 21, 30, "float.6", "rjust"),
        ("s[n][3]", 31, 40, "float.6", "rjust"),
        ("u[n]", 46, 55, "float.5", "rjust")]
    
class SCALE1(SCALEn):
    _name = "SCALE1"
        
class SCALE2(SCALEn):
    _name = "SCALE2"

class SCALE3(SCALEn):
    _name = "SCALE3"

class MTRIXn(PDBRecord):
    """The MTRIXn (n = 1, 2, or 3) records present transformations expressing
    non-crystallographic symmetry."""
    _field_list = [
        ("serial", 8, 10, "integer", "rjust"),
        ("s[n][1]", 11, 20, "float.6", "rjust"),
        ("s[n][2]", 21, 30, "float.6", "rjust"),
        ("s[n][3]", 31, 40, "float.6", "rjust"),
        ("v[n]", 46, 55, "float.5", "rjust"),
        ("iGiven", 60, 60, "integer", "rjust")]
    
class MTRIX1(MTRIXn):
    _name = "MTRIX1"

class MTRIX2(MTRIXn):
    _name = "MTRIX2"

class MTRIX3(MTRIXn):
    _name = "MTRIX3"

class TVECT(PDBRecord):
    """The TVECT records present the translation vector for infinite
    covalently connected structures."""
    _name = "TVECT "
    _field_list = [
        ("serial", 8, 10, "integer", "rjust"),
        ("t[1]", 11, 20, "float.5", "rjust"),
        ("t[2]", 21, 30, "float.5", "rjust"),
        ("t[3]", 31, 40, "float.5", "rjust"),
        ("text", 41, 70, "string", "rjust")]

## SECTION 9: Coordinate Selection
class MODEL(PDBRecord):
    """The MODEL record specifies the model serial number when multiple
    structures are presented in a single coordinate entry, as is often
    the case with structures determined by NMR."""
    _name = "MODEL "
    _field_list = [
        ("serial", 11, 14, "integer", "rjust")]

class ATOM(PDBRecord):
    """The ATOM records present the atomic coordinates for standard residues.
    They also present the occupancy and temperature factor for each atom.
    Heterogen coordinates use the HETATM record type. The element symbol
    is always present on each ATOM record; segment identifier and charge
    are optional."""
    _name = "ATOM  "
    _field_list = [
        ("serial", 7, 11, "integer", "rjust"),
        ("name", 13, 16, "string", "ljust"),
        ("altLoc", 17, 17, "string", "rjust"),
        ("resName", 18, 20, "string","rjust"),
        ("chainID", 22, 22, "string", "rjust"),
        ("resSeq", 23, 26, "integer", "rjust"),
        ("iCode", 27, 27, "string", "rjust"),
        ("x", 31, 38, "float.3", "rjust"),
        ("y", 39, 46, "float.3", "rjust"),
        ("z", 47, 54, "float.3", "rjust"),
        ("occupancy", 55, 60, "float.2", "rjust"),
        ("tempFactor", 61, 66, "float.2", "rjust"),
        ("segID", 73, 76, "string", "rjust"),
        ("element", 77, 78, "string", "rjust"),
        ("charge", 79, 80, "string", "rjust")]
    
class ANISOU(PDBRecord):
    """The ANISOU records present the anisotropic temperature factors.
    Columns 7 - 27 and 73 - 80 are identical to the corresponding
    ATOM/HETATM record."""
    _name = "ANISOU"
    _field_list = [
        ("serial", 7, 11, "integer", "rjust"),
        ("name", 13, 16, "string", "ljust"),
        ("altLoc", 17, 17, "string", "rjust"),
        ("resName", 18, 20, "string","rjust"),
        ("chainID", 22, 22, "string", "rjust"),
        ("resSeq", 23, 26, "integer", "rjust"),
        ("iCode", 27, 27, "string", "rjust"),
        ("u[0][0]", 29, 35, "integer", "rjust"),
        ("u[1][1]", 36, 42, "integer", "rjust"),
        ("u[2][2]", 43, 49, "integer", "rjust"),
        ("u[0][1]", 50, 56, "integer", "rjust"),
        ("u[0][2]", 57, 63, "integer", "rjust"),
        ("u[1][2]", 64, 70, "integer", "rjust"),
        ("segID", 73, 76, "string", "rjust"),
        ("element", 77, 78, "string", "rjust"),
        ("charge", 79, 80, "string", "rjust")]

class HETATM(ATOM):
    """The HETATM records present the atomic coordinate records for atoms
    within "non-standard" groups. These records are used for water
    molecules and atoms presented in HET groups."""
    _name = "HETATM"

class SIGATM(PDBRecord):
    """The SIGATM records present the standard deviation
    of atomic parameters as they appear in ATOM and HETATM records.
    Columns 7 - 27 and 73 - 80 are identical to the corresponding
    ATOM/HETATM record."""
    _name = "SIGATM"
    _field_list = [
        ("serial", 7, 11, "integer", "rjust"),
        ("name", 13, 16, "string", "ljust"),
        ("altLoc", 17, 17, "string", "rjust"),
        ("resName", 18, 20, "string","rjust"),
        ("chainID", 22, 22, "string", "rjust"),
        ("resSeq", 23, 26, "integer", "rjust"),
        ("iCode", 27, 27, "string", "rjust"),
        ("sigX", 31, 38, "float.3", "rjust"),
        ("sigY", 39, 46, "float.3", "rjust"),
        ("sigZ", 47, 54, "float.3", "rjust"),
        ("sigOccupancy", 55, 60, "float.2", "rjust"),
        ("sigTempFactor", 61, 66, "float.2", "rjust"),
        ("segID", 73, 76, "string", "rjust"),
        ("element", 77, 78, "string", "rjust"),
        ("charge", 79, 80, "string", "rjust")]

class SIGUIJ(PDBRecord):
    """The SIGUIJ records present the standard deviations of anisotropic
    temperature factors scaled by a factor of 10**4 (Angstroms**2). 
    Columns 7 - 27 and 73 - 80 are identical to the corresponding
    ATOM/HETATM record."""
    _name = "SIGUIJ"
    _field_list = [
        ("serial", 7, 11, "integer", "rjust"),
        ("name", 13, 16, "string", "ljust"),
        ("altLoc", 17, 17, "string", "rjust"),
        ("resName", 18, 20, "string","rjust"),
        ("chainID", 22, 22, "string", "rjust"),
        ("resSeq", 23, 26, "integer", "rjust"),
        ("iCode", 27, 27, "string", "rjust"),
        ("sig[1][1]", 29, 35, "integer", "rjust"),
        ("sig[2][2]", 36, 42, "integer", "rjust"),
        ("sig[3][3]", 43, 49, "integer", "rjust"),
        ("sig[1][2]", 50, 56, "integer", "rjust"),
        ("sig[1][3]", 57, 63, "integer", "rjust"),
        ("sig[2][3]", 64, 70, "integer", "rjust"),
        ("segID", 73, 76, "string", "rjust"),
        ("element", 77, 78, "string", "rjust"),
        ("charge", 79, 80, "string", "rjust")]

class TER(PDBRecord):
    """The TER record indicates the end of a list of ATOM/HETATM records
    for a chain."""
    _name = "TER   "
    _field_list = [
        ("serial", 7, 11, "integer", "rjust"),
        ("resName", 18, 20, "string", "rjust"),
        ("chainID", 22, 22, "string", "rjust"),
        ("resSeq", 23, 26, "integer", "rjust"),
        ("iCode", 27, 27, "string", "rjust")]
    
class ENDMDL(PDBRecord):
    """The ENDMDL records are paired with MODEL records to group individual
    structures found in a coordinate entry."""
    _name = "ENDMDL"
    _field_list = []
    
## SECTION 10: Connectivity Section
class CONECT(PDBRecord):
    """The CONECT records specify connectivity between atoms for which
    coordinates are supplied. The connectivity is described using the
    atom serial number as found in the entry. CONECT records are
    mandatory for HET groups (excluding water) and for other bonds not
    specified in the standard residue connectivity table which involve
    atoms in standard residues (see Appendix 4 for the list of standard
    residues). These records are generated by the PDB."""
    _name = "CONECT"
    _field_list = [
        ("serial", 7, 11, "integer", "rjust"),
        ("serialBond1", 12, 16, "integer", "rjust"),
        ("serialBond2", 17, 21, "integer", "rjust"),
        ("serialBond3", 22, 26, "integer", "rjust"),
        ("serialBond4", 27, 31, "integer", "rjust"),
        ("serialHydBond1", 32, 36, "integer", "rjust"),
        ("serialHydBond2", 37, 41, "integer", "rjust"),
        ("serialSaltBond1", 42, 46, "integer", "rjust"),
        ("serialHydBond3", 47, 51, "integer", "rjust"),
        ("serialHydBond4", 52, 56, "integer", "rjust"),
        ("serialSaltBond2", 57, 61, "integer", "rjust")]

## SECTION 11: Bookkeeping Section
class MASTER(PDBRecord):
    """The MASTER record is a control record for bookkeeping. It lists the
    number of lines in the coordinate entry or file for selected record
    types."""
    _name = "MASTER"
    _field_list = [
        ("numRemark", 11, 15, "integer", "rjust"),
        ("O", 16, 20, "integer", "rjust"),
        ("numHet", 21, 25, "integer", "rjust"),
        ("numHelix", 26, 30, "integer", "rjust"),
        ("numSheet", 31, 35, "integer", "rjust"),
        ("numTurn", 36, 40, "integer", "rjust"),
        ("numSite", 41, 45, "integer", "rjust"),
        ("numXForm", 46, 50, "integer", "rjust"),
        ("numCoord", 51, 55, "integer", "rjust"),
        ("numTer", 56, 60, "integer", "rjust"),
        ("numConect", 61, 65, "integer", "rjust"),
        ("numSeq", 66, 70, "integer", "rjust")]

class END(PDBRecord):
    """The END record marks the end of the PDB file."""
    _name = "END   "
    _field_list = []


## PDB Record Name -> Record Class Map
PDBRecordMap = {
    HEADER._name : HEADER,
    OBSLTE._name : OBSLTE,
    TITLE._name  : TITLE,
    CAVEAT._name : CAVEAT,
    COMPND._name : COMPND,
    SOURCE._name : SOURCE,
    KEYWDS._name : KEYWDS,
    EXPDTA._name : EXPDTA,
    AUTHOR._name : AUTHOR,
    REVDAT._name : REVDAT,
    SPRSDE._name : SPRSDE,
    JRNL._name   : JRNL,
    REMARK._name : REMARK,
    DBREF._name  : DBREF,
    SEQADV._name : SEQADV,
    SEQRES._name : SEQRES,
    MODRES._name : MODRES,
    HET._name    : HET,
    HETNAM._name : HETNAM,
    HETSYN._name : HETSYN,
    FORMUL._name : FORMUL,
    HELIX._name  : HELIX,
    SHEET._name  : SHEET,
    TURN._name   : TURN,
    SSBOND._name : SSBOND,
    LINK._name   : LINK,
    HYDBND._name : HYDBND,
    SLTBRG._name : SLTBRG,
    CISPEP._name : CISPEP,
    SITE._name   : SITE,
    CRYST1._name : CRYST1,
    CRYST2._name : CRYST2,
    CRYST3._name : CRYST3,
    ORIGX1._name : ORIGX1,
    ORIGX2._name : ORIGX2,
    ORIGX3._name : ORIGX3,
    SCALE1._name : SCALE1,
    SCALE2._name : SCALE2,
    SCALE3._name : SCALE3,
    MTRIX1._name : MTRIX1,
    MTRIX2._name : MTRIX2,
    MTRIX3._name : MTRIX3,
    MODEL._name  : MODEL,
    ATOM._name   : ATOM,
    ANISOU._name : ANISOU,
    HETATM._name : HETATM,
    SIGATM._name : SIGATM,
    SIGUIJ._name : SIGUIJ,
    TER._name    : TER,
    ENDMDL._name : ENDMDL,
    CONECT._name : CONECT,
    MASTER._name : MASTER,
    END._name    : END }

## this list defines the order the records have to appear in the PDB
## file; there is also a indicator if the record is optional or mandatory
PDBRecordOrder = [
    (HEADER._name, HEADER, "mandatory"),
    (OBSLTE._name, OBSLTE, "optional"),
    (TITLE._name, TITLE, "mandatory"),
    (CAVEAT._name, CAVEAT, "optional"),
    (COMPND._name, COMPND, "mandatory"),
    (SOURCE._name, SOURCE, "mandatory"),
    (KEYWDS._name, KEYWDS, "mandatory"),
    (EXPDTA._name, EXPDTA, "mandatory"),
    (AUTHOR._name, AUTHOR, "mandatory"),
    (REVDAT._name, REVDAT, "mandatory"),
    (SPRSDE._name, SPRSDE, "optional"),
    (JRNL._name, JRNL, "optional"),
    (REMARK._name, REMARK, "optional"),
    (DBREF._name, DBREF, "optional"),
    (SEQADV._name, SEQADV, "optional"),
    (SEQRES._name, SEQRES, "optional"),
    (MODRES._name, MODRES, "optional"),
    (HET._name, HET, "optional"),
    (HETNAM._name, HETNAM, "optional"),
    (HETSYN._name, HETSYN, "optional"),
    (FORMUL._name, FORMUL, "optional"),
    (HELIX._name, HELIX, "optional"),
    (SHEET._name, SHEET, "optional"),
    (TURN._name, TURN, "optional"),
    (SSBOND._name, SSBOND, "optional"),
    (LINK._name, LINK, "optional"),
    (HYDBND._name, HYDBND, "optional"),
    (SLTBRG._name, SLTBRG, "optional"),
    (CISPEP._name, CISPEP, "optional"),
    (SITE._name, SITE, "optional"),
    (CRYST1._name, CRYST1, "mandatory"),
    (ORIGX1._name, ORIGX1, "mandatory"),
    (ORIGX2._name, ORIGX2, "mandatory"),
    (ORIGX3._name, ORIGX3, "mandatory"),
    (SCALE1._name, SCALE1, "mandatory"),
    (SCALE2._name, SCALE2, "mandatory"),
    (SCALE3._name, SCALE3, "mandatory"),
    (MTRIX1._name, MTRIX1, "optional"),
    (MTRIX2._name, MTRIX2, "optional"),
    (MTRIX3._name, MTRIX3, "optional"),
    (TVECT._name, TVECT, "optional"),
    (MODEL._name, MODEL, "optional"),
    (ATOM._name, ATOM, "optional"),
    (SIGATM._name, SIGATM, "optional"),
    (ANISOU._name, ANISOU, "optional"),
    (SIGUIJ._name, SIGUIJ, "optional"),
    (TER._name, TER, "optional"),
    (HETATM._name, HETATM, "optional"),
    (ENDMDL._name, ENDMDL, "optional"),
    (CONECT._name, CONECT, "optional"),
    (MASTER._name, MASTER, "mandatory"),
    (END._name, END, "mandatory")
    ]

## END PDB RECORD DEFINITIONS
###############################################################################


class PDBFile:
    """Class for managing a PDB file.  Load, save, edit, and create PDB
    files with this class."""

    def __init__(self):
        self.pdbrecord_list = []


    def loadFile(self, fil):
        fil = OpenFile(fil, "r")
        for ln in fil.readlines():
            ## prep line by removing newline and upper-caseing line
            ln = string.rstrip(ln)
            ln = string.upper(ln)
            ## find the record data element for the given line
            try:
                rname = ln[:6]
            except IndexError:
                rname = ln[:6].ljust(6)

            try:
                pdb_record_class = PDBRecordMap[rname]
            except KeyError:
                print "Unknown record type=%s" % (rname)
                continue

            ## create/add/parse the record
            pdb_record = pdb_record_class()
            self.pdbrecord_list.append(pdb_record)
            pdb_record.read(ln)


    def saveFile(self, fil):
        fil = OpenFile(fil, "w")

        ## this re-orders the PDB records with some basic grouping
        ## rules in the PDB v2.2 specification, it doesn't
        ## get it completely correct yet...
        self.doctorUp()

        for pdb_record in self.pdbrecord_list:
            fil.write(str(pdb_record) + "\n")
        fil.close()


    def doctorUp(self):
        ## first sort the records
        record_map = {}
        for rec in self.pdbrecord_list:
            try:
                record_map[rec._name].append(rec)
            except KeyError:
                record_map[rec._name] = []
                record_map[rec._name].append(rec)

        ## create new list in correct order, adding mandatory records
        serial_num = 1
        pdbrecord_list = []
        for (name, rec_class, requirement) in PDBRecordOrder:

            ## skip some record types that require special sort/order
            ## handling
            if name in [SIGATM._name, ANISOU._name, SIGUIJ._name, TER._name]:
                continue

            ## ATOM and HETATM records break the normal ordering rules,
            ## they have associated records SIGATM, ANISOU, and SIGUIJ
            ## which should appear in that order after them
            if name in [ATOM._name, HETATM._name]:

                ## if there are no records, continue
                if not record_map.has_key(name):
                    continue

                ## this function searches for other records describing
                ## the atom given in atom_rec from the list of records
                ## passed in rec_list
                def get_rec(atom_rec, rec_list):
                    for rec in rec_list:
                        if not atom_rec.compareRecord(rec, AtomIDFields):
                            continue
                        rec_list.remove(rec)
                        return rec
                    return None

                for atom_rec in record_map[name]:
                    atom_rec.serial = serial_num
                    
                    pdbrecord_list.append(atom_rec)

                    if record_map.has_key(SIGATM._name):
                        rec = get_rec(atom_rec, record_map[SIGATM._name])
                        if rec:
                            rec.serial = serial_num
                            pdbrecord_list.append(rec)

                    if record_map.has_key(ANISOU._name):
                        rec = get_rec(atom_rec, record_map[ANISOU._name])
                        if rec:
                            rec.serial = serial_num
                            pdbrecord_list.append(rec)

                    if record_map.has_key(SIGUIJ._name):
                        rec = get_rec(atom_rec, record_map[SIGUIJ._name])
                        if rec:
                            rec.serial = serial_num
                            pdbrecord_list.append(rec)

                    serial_num += 1

            elif requirement == "mandatory" and not record_map.has_key(name):
                pdbrecord_list.append(rec_class())
            elif record_map.has_key(name):
                pdbrecord_list += record_map[name]

        self.pdbrecord_list = pdbrecord_list


    def getRecordList(self):
        return self.pdbrecord_list


    def selectRecordList(self, *nvlist):
        """Preforms a SQL-like 'AND' select aginst all the records in the
        PDB file.  The arguments are a variable list of tuples of the form:
          (<column-name>, <column-value>)
        For example:
          selectRecordList(('_name','ATOM  '),('resName', 'LYS'))
        returns a list of ATOM records which are part of a LYS residue."""
        ## clever optimization trickies
        (attr, val) = nvlist[0]

        rec_list = []
        for rec in self.pdbrecord_list:
            try:
                if getattr(rec, attr) != val: continue
            except AttributeError:
                continue

            add = 1
            for (attr, val) in nvlist:
                try:
                    if getattr(rec, attr) != val:
                        add = 0
                        break
                except AttributeError:
                    add = 0
                    break

            if add:
                rec_list.append(rec)

        return rec_list


    def getChainList(self):
        """Returns a list of all the chain ids in the PDB file."""
        chain_list = []
        for pdb_record in self.pdbrecord_list:
            try:
                chain_id = getattr(pdb_record, "chainID")
            except AttributeError:
                continue

            if chain_id not in chain_list:
                chain_list.append(chain_id)
        return chain_list



##
## <testing>
##
    
if __name__ == "__main__":
    import sys

    try:
        path = sys.argv[1]
    except IndexError:
        print "usage: PDB.py <PDB file path>"
        sys.exit(1)

    pdbfil = PDBFile()
    pdbfil.loadFile(path)
    pdbfil.saveFile(sys.stdout)

##
## </testing>
##
