## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Brookhaven PDB v2.2 file parser.  All records in the PDB v2.2
specification have coorasponding classes defined here.  PDB files are
loaded into a list of these cassed, and also can be constrcted/modified
and written back out as PDB files.
"""
from __future__ import generators
import string
import types
import fpformat
from Structure import *


class PDBError(Exception):
    """
    """
    pass


class PDBValueError(PDBError):
    """
    """
    def __init__(self, text):
        self.text = text

    def __str__(self):
        return self.text


class PDBRecord(dict):
    """Base class for all PDB file records.
    """
    def __str__(self):
        return self.write()

    def write(self):
        """Return a properly formed PDB record string from the instance
        dictionary values.
        """
        ln = self._name
        
        for (field, start, end, ftype, just, get_func) in self._field_list:

            ## add spaces to the end if necessary
            assert len(ln) <= (start - 1)
            ln = ln.ljust(start - 1)

            ## used later
            field_char_len = end - start + 1
            

            ## access the namespace of this class to write the field
            ## if a class has a special function defined for retrieveing
            ## this record, it should use it
            if get_func:
                ln += get_func(self)
                continue

            ## get the data
            s = self.get(field, "")

            ## if the data is blank, then just add the spaces and continue
            if s == None or s == "":
                ln += " " * field_char_len
                continue

            ## convert integer and float types
            if ftype.startswith("string"):
                pass

            elif ftype.startswith("integer"):
                s = str(s)

            elif ftype.startswith("float"):
                try:
                    s = fpformat.fix(s, int(ftype[6]))
                except ValueError:
                    raise PDBValueError("field=%s %s not float" % (field, s))

            ## assert type
            try:
                assert type(s) == StringType
            except AssertionError:
                print "### s",str(type(s)), str(s), ftype, field
                print ln
                raise

            ## check for maximum length
            if len(s) > field_char_len:
                ln += s[:field_char_len]
            else:
                if just.startswith("ljust"):
                    ln += s.ljust(field_char_len)
                else:
                    ln += s.rjust(field_char_len)

        return ln

    def read(self, line):
        """Read the PDB record line and convert the fields to the appropriate
        dictionary values for this class.
        """
        for (field, start, end, ftype, just, get_func) in self._field_list:
            s = line[start-1:end]

            ## ignore blank fields
            if s.isspace():
                continue

            elif ftype.startswith("string"):
                
                if just.endswith("lstrip"):
                    s = s.lstrip()
                elif just.endswith("rstrip"):
                    s = s.rstrip()
                else:
                    s = s.strip()

            elif ftype.startswith("integer"):
                try:
                    s = int(s)
                except ValueError:
                    if s.strip() != "":
                        warning("PDB parser: int(%s) failed on record" % (s))
                        warning(str(line))
                    continue

            elif ftype.startswith("float"):
                try:
                    s = float(s)
                except ValueError:
                    if s.strip() != "":
                        warning("PDB parser: float(%s) failed on record" % (s))
                        warning(str(line))
                    continue

            self[field] = s

    def reccat(self, rec_list, field):
        """Return the concatenation of field in all the records in rec_list.
        """
        if type(rec_list) != ListType:
            rec_list = [rec_list]

        retval = ""
        for rec in rec_list:
            x = rec.get(field)
            if x != None:
                retval += x
        return retval

    def reccat_list(self, rec_list, field, sep):
        """Call reccat, then split the result by the seperator.
        """
        listx = self.reccat(rec_list, field).split(sep)
        listx = [x.strip() for x in listx]
        return listx

    def reccat_tuplelist(self, rec_list, field, sep1, sep2):
        """Call reccat_list with sep1 as the list seperator, then split
        the items into tuples by sep2.
        """
        listx = []
        for x in self.reccat_list(rec_list, field, sep1):
            i = x.find(sep2)
            if i == -1:
                continue
            key = x[:i].strip()
            val = x[i+1:].strip()
            listx.append((key, val))
        return listx  

    def reccat_dictlist(self, rec_list, field, master_key):
        listx = []
        dictx = {}
        for (key, val) in self.reccat_tuplelist(rec_list, field, ";", ":"):
            if key == master_key:
                if dictx:
                    listx.append(dictx)
                    dictx = {}
            dictx[key] = val
        if dictx:
            listx.append(dictx)
        return listx

    def reccat_multi(self, rec_list, primary_key, translations):
        """Create a list of dictionaries from a list of records.  This
        method has complex behavior to support translations of several
        PDB records into a Python format.  The primary key is used to
        seperate the dictionaries within the list, and the translation
        argument is a list of strings or 2-tuples,  If the translation is a
        string, the value from the PDB record field is copied to the return
        dictionary.  If the field is a 2-tuple==t, then t[0] is the return
        dictionary key whose value is a list formed from the list of
        PDB fields in t[1].
        """
        if type(rec_list) != ListType:
            rec_list = [rec_list]
        
        listx = []
        for rec in rec_list:

            ## XXX: add primary key generation for bad records
            try:
                pkey = rec[primary_key]
            except KeyError:
                ## if the record has no primary key, retrieve it from the
                ## last dictionary which is in the same order as the
                ## record list
                try:
                    pkey = listx[-1][primary_key]
                except KeyError:
                    continue
                except IndexError:
                    continue

            ## search for a dictionary in listx with the same primary key
            dictx = None
            for dx in listx:
                if dx[primary_key] == pkey:
                    dictx = dx
                    break

            ## new dictx if not found
            if dictx == None:
                dictx = {primary_key: pkey}
                listx.append(dictx)

            ## translate the PDB record into dictx
            for trans in translations:

                ## source is a list of fields which should be
                ## added to a list under the dest key in dictx
                if type(trans) == TupleType:
                    (dest, srcs) = trans
                    
                    for sx in srcs:
                        if dictx.has_key(dest):
                            try:
                                dictx[dest].append(rec[sx])
                            except KeyError:
                                pass
                        else:
                            try:
                                dictx[dest] = [rec[sx]]
                            except KeyError:
                                pass

                ## source is a single record field which should be
                ## added to dictx under the dest key
                else:
                    try:
                        dictx[trans] = rec[trans]
                    except KeyError:
                        pass

        return listx    

        

    
###############################################################################
## BEGIN PDB RECORD DEFINITIONS

## SECTION 2: Title Section
class HEADER(PDBRecord):
    """This section contains records used to describe the experiment and the
    biological macromolecules present in the entry: HEADER, OBSLTE, TITLE,
    CAVEAT, COMPND, SOURCE, KEYWDS, EXPDTA, AUTHOR, REVDAT, SPRSDE, JRNL,
    and REMARK records.
    """
    __slots__ = []
    
    _name = "HEADER"
    _field_list = [
        ("classification", 11, 50, "string", "rjust", None),
        ("depDate", 51, 59, "string", "rjust", None),
        ("idCode", 63, 66, "string", "rjust", None)]

class OBSLTE(PDBRecord):
    """OBSLTE appears in entries which have been withdrawn from distribution.
    This record acts as a flag in an entry which has been withdrawn from the
    PDB's full release. It indicates which, if any, new entries have replaced
    the withdrawn entry.  The format allows for the case of multiple new
    entries replacing one existing entry.
    """
    __slots__ = []
    
    _name = "OBSLTE"
    _field_list = [
        ("continuation", 9, 10, "string", "rjust", None),
        ("repDate", 12, 20, "string", "rjust", None),
        ("idCode", 22, 25, "string", "rjust", None),
        ("rIdCode1", 32, 35, "string", "rjust", None),
        ("rIdCode2", 37, 40, "string", "rjust", None),
        ("rIdCode3", 42, 45, "string", "rjust", None),
        ("rIdCode4", 47, 50, "string", "rjust", None),
        ("rIdCode5", 52, 55, "string", "rjust", None),
        ("rIdCode6", 57, 60, "string", "rjust", None),
        ("rIdCode7", 62, 65, "string", "rjust", None),
        ("rIdCode8", 67, 70, "string", "rjust", None)]

    def process(self, recs):
        """Processes continued record list to a list of dictionary objects.
        Each dictionary contains the data from one OBSLTE idCode.
        """
        return self.reccat_multi(
            recs, "idCode",
            ["repDate",
             ("rIdCodes", ["rIdCode1", "rIdCode2", "rIdCode3", "rIdCode4",
                           "rIdCode5", "rIdCode6", "rIdCode7", "rIdCode8"])])
    
class TITLE(PDBRecord):
    """The TITLE record contains a title for the experiment or analysis that is
    represented in the entry. It should identify an entry in the PDB in the
    same way that a title identifies a paper.
    """
    __slots__ = []
    
    _name = "TITLE "
    _field_list = [
        ("continuation", 9, 10, "string", "rjust", None),
        ("title", 11, 70, "string", "ljust", None)]

    def process(self, recs):
        return self.reccat(recs, "title")

class CAVEAT(PDBRecord):
    """CAVEAT warns of severe errors in an entry. Use caution when using an
    entry containing this record.
    """
    __slots__ = []
    
    _name = "CAVEAT"
    _field_list = [
        ("continuation", 9, 10, "string", "rjust", None),
        ("idCode", 12, 15, "string", "rjust", None),
        ("comment", 20, 70, "string", "ljust", None)]

    def process(self, recs):
        """Returns a list of dictionaries with keys idCode and comment.
        """
        cavet_list = []
        for rec in recs:
            idCode = rec.get("idCode")
            if idCode == None:
                continue

            ## search for cavet entry with same idCode
            cav = None
            for cavx in cavet_list:
                if cavx.get("idCode") == idCode:
                    cav = cavx
                    break

            ## create new cavet dict if necessary
            if cav == None:
                cav = {"idCode" : idCode}
                cavet_list.append(cav)

            ## add comment
            comment = rec.get("comment")
            if comment != None:
                if cav.has_key("comment"):
                    cav["comment"] += comment
                else:
                    cav["comment"] = comment
        return cavet_list

class COMPND(PDBRecord):
    """The COMPND record describes the macromolecular contents of an entry.
    Each macromolecule found in the entry is described by a set of token: value
    pairs, and is referred to as a COMPND record component. Since the concept
    of a molecule is difficult to specify exactly, PDB staff may exercise
    editorial judgment in consultation with depositors in assigning these
    names.  For each macromolecular component, the molecule name, synonyms,
    number assigned by the Enzyme Commission (EC), and other relevant details
    are specified.
    """ 
    __slots__ = []
    
    _name = "COMPND"
    _field_list = [
        ("continuation", 9, 10, "string", "rjust", None),
        ("compound", 11, 70, "string", "ljust", None)]

    def process(self, recs):
        return self.reccat_dictlist(recs, "compound", "MOL_ID")

class SOURCE(PDBRecord):
    """The SOURCE record specifies the biological and/or chemical source of
    each biological molecule in the entry. Sources are described by both the
    common name and the scientific name, e.g., genus and species. Strain and/or
    cell-line for immortalized cells are given when they help to uniquely
    identify the biological entity studied.
    """
    __slots__ = []
    
    _name = "SOURCE"
    _field_list = [
        ("continuation", 9, 10, "string", "rjust", None),
        ("srcName", 11, 70, "string", "ljust", None)]

    def process(self, recs):
        return self.reccat_dictlist(recs, "srcName", "MOL_ID")

class KEYWDS(PDBRecord):
    """The KEYWDS record contains a set of terms relevant to the entry. Terms
    in the KEYWDS record provide a simple means of categorizing entries and may
    be used to generate index files. This record addresses some of the
    limitations found in the classification field of the HEADER record. It
    provides the opportunity to add further annotation to the entry in a
    concise and computer-searchable fashion.
    """
    __slots__ = []
    
    _name = "KEYWDS"
    _field_list = [
        ("continuation", 9, 10, "string", "rjust", None),
        ("keywds", 11, 70, "string", "ljust", None)]

    def process(self, recs):
        return self.reccat_list(recs, "keywds", ",")

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
    X-RAY DIFFRACTION
    """ 
    __slots__ = []
    
    _name = "EXPDTA"
    _field_list = [
        ("continuation", 9, 10, "string", "rjust", None),
        ("technique", 11, 70, "string", "ljust", None)]
    _technique_list = [
        "ELECTRON DIFFRACTION",
        "FIBER DIFFRACTION",
        "FLUORESCENCE TRANSFER",
        "NEUTRON DIFFRACTION",
        "NMR",
        "THEORETICAL MODEL",
        "X-RAY DIFFRACTION"]

    def process(self, recs):
        """Returns a list of 2-tuples: (technique, comment) where technique
        is one of the accepted tequniques.
        """
        expdta_list = []

        for item in self.reccat_list(recs, "technique", ";"):
            tech = None
            cmnt = None

            for techx in self._technique_list:
                if item.startswith(techx):
                    tech = techx
                    cmnt = item[len(techx):].strip() or None
                    break

            if tech != None:
                expdta_list.append((tech, cmnt))

        return expdta_list

class AUTHOR(PDBRecord):
    """The AUTHOR record contains the names of the people responsible for the
    contents of the entry.
    """
    __slots__ = []
    
    _name = "AUTHOR"
    _field_list = [
        ("continuation", 9, 10, "string", "rjust", None),
        ("authorList", 11, 70, "string", "ljust", None)]

    def process(self, recs):
        return self.reccat_list(recs, "authorList", ",")

class REVDAT(PDBRecord):
    """REVDAT records contain a history of the modifications made to an entry
    since its release.
    """
    __slots__ = []
    
    _name = "REVDAT"
    _field_list = [
        ("modNum", 8, 10, "integer", "rjust", None),
        ("continuation", 11, 12, "string", "rjust", None),
        ("modDate", 14, 22, "string", "rjust", None),
        ("modID", 24, 28, "string", "rjust", None),
        ("modType", 32, 32, "integer", "rjust", None),
        ("record1", 40, 45, "string", "ljust", None),
        ("record2", 47, 52, "string", "ljust", None),
        ("record3", 54, 59, "string", "ljust", None),
        ("record4", 61, 66, "string", "ljust", None)]

    def process(self, recs):
        return self.reccat_multi(
            recs, "modNum",
            ["modDate",
             "modID",
             "modType",
             ("records", ["record1", "record2", "record3", "record4"])])

class SPRSDE(PDBRecord):
    """The SPRSDE records contain a list of the ID codes of entries that were
    made obsolete by the given coordinate entry and withdrawn from the PDB
    release set. One entry may replace many. It is PDB policy that only the
    principal investigator of a structure has the authority to withdraw it.
    """
    __slots__ = []
    
    _name = "SPRSDE"
    _field_list = [
        ("continuation", 9, 10, "string", "rjust", None),
        ("sprsdeDate", 12, 20, "string", "rjust", None),
        ("idCode", 22, 25, "string", "rjust", None),
        ("sIdCode1", 32, 35, "string", "rjust", None),
        ("sIdCode2", 37, 40, "string", "rjust", None),
        ("sIdCode3", 42, 45, "string", "rjust", None),
        ("sIdCode4", 47, 50, "string", "rjust", None),
        ("sIdCode5", 52, 55, "string", "rjust", None),
        ("sIdCode6", 57, 60, "string", "rjust", None),
        ("sIdCode7", 62, 65, "string", "rjust", None),
        ("sIdCode8", 67, 70, "string", "rjust", None)]

    def process(self, recs):
        return self.reccat_multi(
            recs, "idCode",
            ["sprsdeDate",
             ("sIdCodes", ["sIdCode1", "sIdCode2", "sIdCode3", "sIdCode4",
                           "sIdCode5", "sIdCode6", "sIdCode7", "sIdCode8"])])

class JRNL(PDBRecord):
    """The JRNL record contains the primary literature citation that describes
    the experiment which resulted in the deposited coordinate set. There is at
    most one JRNL reference per entry. If there is no primary reference, then
    there is no JRNL reference. Other references are given in REMARK 1.
    """
    __slots__ = []
    
    _name = "JRNL  "
    _field_list = [
        ("text", 13, 70, "string", "ljust", None)]

class REMARK(PDBRecord):
    """REMARK records present experimental details, annotations, comments, and
    information not included in other records. In a number of cases, REMARKs
    are used to expand the contents of other record types. A new level of
    structure is being used for some REMARK records. This is expected to
    facilitate searching and will assist in the conversion to a relational
    database.
    """
    __slots__ = []
    
    _name = "REMARK"
    _field_list = [
        ("remarkNum", 8, 10, "integer", "rjust", None),
        ("text", 12, 70, "string", "ljust", None)]
        
## SECTION 3: Primary Structure Section
class DBREF(PDBRecord):
    """ The DBREF record provides cross-reference links between PDB sequences
    and the corresponding database entry or entries. A cross reference to
    the sequence database is mandatory for each peptide chain with a length
    greater than ten (10) residues. For nucleic acid entries a DBREF
    record pointing to the Nucleic Acid Database (NDB) is mandatory when
    the corresponding entry exists in NDB.
    """
    __slots__ = []
    
    _name = "DBREF "
    _field_list = [
        ("idCode", 8, 11, "string", "rjust", None),
        ("chain_ID", 13, 13, "string", "rjust", None),
        ("seqBegin", 15, 18, "integer", "rjust", None),
        ("insertBegin", 19, 19, "string", "rjust", None),
        ("seqEnd", 21, 24, "integer", "rjust", None),
        ("insertEnd", 25, 25, "string", "rjust", None),
        ("database", 27, 32, "string", "ljust", None),
        ("dbAccession", 34, 41, "string", "ljust", None),
        ("dbIdCode", 43, 54, "string", "ljust", None),
        ("dbseqBegin", 56, 60, "integer", "rjust", None),
        ("idbnsBeg", 61, 61, "string", "rjust", None),
        ("dbseqEnd", 63, 67, "integer", "rjust", None),
        ("dbinsEnd", 68, 68, "string", "rjust", None)]

class SEQADV(PDBRecord):
    """The SEQADV record identifies conflicts between sequence information
    in the ATOM records of the PDB entry and the sequence database entry
    given on DBREF. Please note that these records were designed to
    identify differences and not errors. No assumption is made as to which
    database contains the correct data. PDB may include REMARK records in
    the entry that reflect the depositor's view of which database has the
    correct sequence.
    """
    __slots__ = []
    
    _name = "SEQADV"
    _field_list = [
        ("idCode", 8, 11, "string", "rjust", None),
        ("resName", 13, 15, "string", "rjust", None),
        ("chainID", 17, 17, "string", "rjust", None),
        ("seqNum", 19, 22, "integer", "rjust", None),
        ("iCode", 23, 23, "string", "rjust", None),
        ("database", 25, 28, "string", "ljust", None),
        ("dbIDCode", 30, 38, "string", "ljust", None),
        ("dbRes", 40, 42, "string", "rjust", None),
        ("dbSeq", 44, 48, "integer", "rjust", None),
        ("convlict", 40, 70, "string", "ljust", None)]
    
class SEQRES(PDBRecord):
    """The SEQRES records contain the amino acid or nucleic acid sequence of
    residues in each chain of the macromolecule that was studied.
    """
    __slots__ = []
    
    _name = "SEQRES"
    _field_list = [
        ("serNum", 9, 10, "integer", "rjust", None),
        ("chainID", 12, 12, "string", "rjust", None),
        ("numRes", 14, 17, "integer", "rjust", None),
        ("resName1", 20, 22, "string", "rjust", None),
        ("resName2", 24, 26, "string", "rjust", None),
        ("resName3", 28, 30, "string", "rjust", None),
        ("resName4", 32, 34, "string", "rjust", None),
        ("resName5", 36, 38, "string", "rjust", None),
        ("resName6", 40, 42, "string", "rjust", None),
        ("resName7", 44, 46, "string", "rjust", None),
        ("resName8", 48, 50, "string", "rjust", None),
        ("resName9", 52, 54, "string", "rjust", None),
        ("resName10", 56, 58, "string", "rjust", None),
        ("resName11", 60, 62, "string", "rjust", None),
        ("resName12", 64, 66, "string", "rjust", None),
        ("resName13", 68, 70, "string", "rjust", None)]

class MODRES(PDBRecord):
    """The MODRES record provides descriptions of modifications (e.g.,
    chemical or post-translational) to protein and nucleic acid residues.
    Included are a mapping between residue names given in a PDB entry and
    standard residues.
    """
    __slots__ = []
    
    _name = "MODRES"
    _field_list = [
        ("idCode", 8, 11, "string", "rjust", None),
        ("resName", 13, 15, "string", "rjust", None),
        ("chainID", 17, 17, "string", "rjust", None),
        ("seqNum", 19, 22, "integer", "rjust", None),
        ("iCode", 23, 23, "string", "rjust", None),
        ("stdRes", 25, 27, "string", "rjust", None),
        ("comment", 30, 70, "string", "ljust", None)]

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
    is unknown, in which case the group is assigned the hetID UNK.
    """
    __slots__ = []
    
    _name = "HET   "
    _field_list = [
        ("hetID", 8, 10, "string", "rjust", None),
        ("chainID", 13, 13, "string", "rjust", None),
        ("seqNum", 14, 17, "integer", "rjust", None),
        ("iCode", 18, 18, "string", "rjust", None),
        ("numHetAtoms", 21, 25, "integer", "rjust", None),
        ("text", 31, 70, "string", "ljust", None)]
    
class HETNAM(PDBRecord):
    """This record gives the chemical name of the compound with the
    given hetID.
    """
    __slots__ = []
    
    _name = "HETNAM"
    _field_list = [
        ("continuation", 9, 10, "string", "ljust", None),
        ("hetID", 12, 14, "string", "rjust", None),
        ("text", 16, 70, "string", "ljust", None)]
    
class HETSYN(PDBRecord):
    """This record provides synonyms, if any, for the compound in the
    corresponding (i.e., same hetID) HETNAM record. This is to allow
    greater flexibility in searching for HET groups.
    """
    __slots__ = []
    
    _name = "HETSYN"
    _field_list = [
        ("continuation", 9, 10, "string", "ljust", None),
        ("hetID", 12, 14, "string", "rjust", None),
        ("hetSynonyms", 16, 70, "string", "ljust", None)]

class FORMUL(PDBRecord):
    """The FORMUL record presents the chemical formula and charge of a
    non-standard group. (The formulas for the standard residues are given
    in Appendix 5.)
    """
    __slots__ = []
    
    _name = "FORMUL"
    _field_list = [
        ("compNum", 9, 10, "integer", "rjust", None),
        ("hetID", 13, 15, "string", "rjust", None),
        ("continuation", 17, 18, "integer", "rjust", None),
        ("asterisk", 19, 19, "string", "rjust", None),
        ("text", 20, 70, "string", "ljust", None)]

## SECTION 5: Secondary Structure Section
class HELIX(PDBRecord):
    """HELIX records are used to identify the position of helices in the
    molecule. Helices are both named and numbered. The residues where the
    helix begins and ends are noted, as well as the total length.
    """
    __slots__ = []
    
    _name = "HELIX "
    _field_list = [
        ("serNum", 8, 10, "integer", "rjust", None),
        ("helixID", 12, 14, "string", "rjust", None),
        ("initResName", 16, 18, "string", "rjust", None),
        ("initChainID", 20, 20, "string", "rjust", None),
        ("initSeqNum", 22, 25, "integer", "rjust", None),
        ("initICode", 26, 26, "string", "rjust", None),
        ("endResName", 28, 30, "string", "rjust", None),
        ("endChainID", 32, 32, "string", "rjust", None),
        ("endSeqNum", 34, 37, "integer", "rjust", None),
        ("endICode", 38, 38, "string", "rjust", None),
        ("helixClass", 39, 40, "integer", "rjust", None),
        ("comment", 41, 70, "string", "ljust", None),
        ("length", 72, 76, "integer", "rjust", None)]

class SHEET(PDBRecord):
    """SHEET records are used to identify the position of sheets in the
    molecule. Sheets are both named and numbered. The residues where the
    sheet begins and ends are noted.
    """
    __slots__ = []
    
    _name = "SHEET "
    _field_list = [
        ("strand", 8, 10, "integer", "rjust", None),
        ("sheetID", 12, 14, "string", "rjust", None),
        ("numStrands", 15, 16, "integer", "rjust", None),
        ("initResName", 18, 20, "string", "rjust", None),
        ("initChainID", 22, 22, "string", "rjust", None),
        ("initSeqNum", 23, 26, "integer", "rjust", None),
        ("initICode", 27, 27, "string", "rjust", None),
        ("endResName", 29, 31, "string", "rjust", None),
        ("endChainID", 33, 33, "string", "rjust", None),
        ("endSeqNum", 34, 37, "integer", "rjust", None),
        ("endICode", 38, 38, "string", "rjust", None),
        ("sense", 39, 40, "integer", "rjust", None),
        ("curAtom", 42, 45, "string", "rjust", None),
        ("curResName", 46, 48, "string", "rjust", None),
        ("curChainID", 50 ,50, "string", "rjust", None),
        ("curResSeq", 51, 54, "integer", "rjust", None),
        ("curICode", 55, 55, "string", "rjust", None),
        ("prevAtom", 57, 60, "string", "rjust", None),
        ("prevResName", 61, 63, "string", "rjust", None),
        ("prevChainID", 65, 65, "string", "rjust", None),
        ("prevResSeq", 66, 69, "integer", "rjust", None),
        ("prevICode", 70, 70, "string", "rjust", None)]

class TURN(PDBRecord):
    """The TURN records identify turns and other short loop turns which
    normally connect other secondary structure segments.
    """
    __slots__ = []
    
    _name = "TURN  "
    _field_list = [
        ("seq", 8, 10, "integer", "rjust", None),
        ("turnID", 12, 14, "string", "rjust", None),
        ("initResName", 16, 18, "string", "rjust", None),
        ("initChainID", 20, 20, "string", "rjust", None),
        ("initSeqNum", 21, 24, "integer", "rjust", None),
        ("initICode", 25, 25, "string", "rjust", None),
        ("endResName", 27, 29, "string", "rjust", None),
        ("endChainID", 31, 31, "string", "rjust", None),
        ("endSeqNum", 32, 35, "integer", "rjust", None),
        ("endICode", 36, 36, "string", "rjust", None),
        ("comment", 41, 70, "string", "ljust", None)]

## SECTION 6: Connectivity Annotation Section
class SSBOND(PDBRecord):
    """The SSBOND record identifies each disulfide bond in protein and
    polypeptide structures by identifying the two residues involved in the
    bond.
    """
    __slots__ = []
    
    _name = "SSBOND"
    _field_list = [
        ("serNum", 8, 10, "integer", "rjust", None),
        ("resName1", 12, 14, "string", "rjust", None),
        ("chainID1", 16, 16, "string", "rjust", None),
        ("seqNum1", 18, 21, "integer", "rjust", None),
        ("iCode1", 22, 22, "string", "rjust", None),
        ("resName2", 26, 28, "string", "rjust", None),
        ("chainID2", 30, 30, "string", "rjust", None),
        ("seqNum2", 32, 35, "integer", "rjust", None),
        ("iCode2", 36, 36, "string", "rjust", None),
        ("sym1", 60, 65, "string", "rjust", None),
        ("sym2", 67, 72, "string", "rjust", None)]

class LINK(PDBRecord):
    """The LINK records specify connectivity between residues that is not
    implied by the primary structure. Connectivity is expressed in terms of
    the atom names. This record supplements information given in CONECT
    records and is provided here for convenience in searching.
    """
    __slots__ = []
    
    _name = "LINK  "
    _field_list = [
        ("name1", 13, 16, "string", "rjust", None),
        ("altLoc1", 17, 17, "string", "rjust", None),
        ("resName1", 18, 20, "string", "rjust", None),
        ("chainID1", 22, 22, "string", "rjust", None),
        ("resSeq1", 23, 26, "integer", "rjust", None),
        ("iCode1", 27, 27, "string", "rjust", None),
        ("name2", 43, 46, "string", "rjust", None),
        ("altLoc2", 47, 47, "string", "rjust", None),
        ("resName2", 48, 50, "string", "rjust", None),
        ("chainID2", 52, 52, "string", "rjust", None),
        ("resSeq2", 53, 56, "integer", "rjust", None),
        ("iCode2", 57, 57, "string", "rjust", None),
        ("sym1", 60, 65, "string", "rjust", None),
        ("sym2", 67, 72, "string", "rjust", None)]

class HYDBND(PDBRecord):
    """The HYDBND records specify hydrogen bonds in the entry.
    """
    __slots__ = []
    
    _name = "HYDBND"
    _field_list = [
        ("name1", 13, 16, "string", "rjust", None),
        ("altLoc1", 17, 17, "string", "rjust", None),
        ("resName1", 18, 20, "string", "rjust", None),
        ("chainID1", 22, 22, "string", "rjust", None),
        ("resSeq1", 23, 27, "integer", "rjust", None),
        ("iCode1", 28, 28, "string", "rjust", None),
        ("nameH", 30, 33, "string", "rjust", None),
        ("altLocH", 34, 34, "string", "rjust", None),
        ("chainH", 36, 36, "string", "rjust", None),
        ("resSeqH", 37, 41, "integer", "rjust", None),
        ("iCodeH", 42, 42, "string", "rjust", None),
        ("name2", 44, 47, "string", "rjust", None),
        ("altLoc2", 48, 48, "string", "rjust", None),
        ("resName2", 49, 51, "string", "rjust", None),
        ("chainID2", 53, 53, "string", "rjust", None),
        ("resSeq2", 54, 58, "integer", "rjust", None),
        ("iCode2", 59, 59, "string", "rjust", None),
        ("sym1", 60, 65, "string", "rjust", None),
        ("sym2", 67, 72, "string", "rjust", None)]

class SLTBRG(PDBRecord):
    """The SLTBRG records specify salt bridges in the entry.
    """
    __slots__ = []
    
    _name = "SLTBRG"
    _field_list = [
        ("name1", 13, 16, "string", "rjust", None),
        ("altLoc1", 17, 17, "string", "rjust", None),
        ("resName1", 18, 20, "string", "rjust", None),
        ("chainID1", 22, 22, "string", "rjust", None),
        ("resSeq1", 23, 26, "integer", "rjust", None),
        ("iCode1", 27, 27, "string", "rjust", None),
        ("name2", 43, 46, "string", "rjust", None),
        ("altLoc2", 47, 47, "string", "rjust", None),
        ("resName2", 48, 50, "string", "rjust", None),
        ("chainID2", 52, 52, "string", "rjust", None),
        ("resSeq2", 53, 56, "integer", "rjust", None),
        ("iCode2", 57, 57, "string", "rjust", None),
        ("sym1", 60, 65, "string", "rjust", None),
        ("sym2", 67, 72, "string", "rjust", None)]

class CISPEP(PDBRecord):
    """CISPEP records specify the prolines and other peptides found to be
    in the cis conformation. This record replaces the use of footnote records
    to list cis peptides.
    """
    __slots__ = []
    
    _name = "CISPEP"
    _field_list = [
        ("serial", 8, 10, "integer", "rjust", None),
        ("resName1", 12, 14, "string", "rjust", None),
        ("chainID1", 16, 16, "string", "rjust", None),
        ("seqNum1", 18, 21, "integer", "rjust", None),
        ("iCode1", 22, 22, "string", "rjust", None),
        ("resName2", 26, 28, "string", "rjust", None),
        ("chainID2", 30, 30, "string", "rjust", None),
        ("seqNum2", 32, 35, "integer", "rjust", None),
        ("iCode2", 36, 36, "string", "rjust", None),
        ("modNum", 44, 46, "integer", "rjust", None),
        ("measure", 54, 59, "float.2", "rjust", None)]
    
## SECTION 7: Miscellanious Features Section
class SITE(PDBRecord):
    """The SITE records supply the identification of groups comprising
    important sites in the macromolecule.
    """
    __slots__ = []
    
    _name = "SITE  "
    _field_list = [
        ("seqNum", 8, 10, "integer", "rjust", None),
        ("siteID", 12, 14, "string", "rjust", None),
        ("numRes", 16, 17, "integer", "rjust", None),
        ("resName1", 19, 21, "string", "rjust", None),
        ("chainID1", 23, 23, "string", "rjust", None),
        ("seq1", 24, 27, "integer", "rjust", None),
        ("iCode1", 28, 28, "string", "rjust", None),
        ("resName2", 30, 32, "string", "rjust", None),
        ("chainID2", 34, 34, "string", "rjust", None),
        ("seq2", 35, 38, "integer", "rjust", None),
        ("iCode2", 39, 39, "string", "rjust", None),
        ("resName3", 41, 43, "string", "rjust", None),
        ("chainID3", 45, 45, "string", "rjust", None),
        ("seq3", 46, 49, "integer", "rjust", None),
        ("iCode3", 50, 50, "string", "rjust", None),
        ("resName4", 52, 54, "string", "rjust", None),
        ("chainID4", 56, 56, "string", "rjust", None),
        ("seq4", 57, 60, "integer", "rjust", None),
        ("iCode4", 61, 61, "string", "rjust", None)]

## SECTION 8: Crystallographic and Coordinate Transformation Section
class CRYSTn(PDBRecord):
    """The CRYSTn (n=1,2,3) record presents the unit cell parameters, space
    group, and Z value. If the structure was not determined by crystallographic
    means, CRYSTn simply defines a unit cube.
    """
    __slots__ = []
    
    _field_list = [
        ("a", 7, 15, "float.3", "rjust", None),
        ("b", 16, 24, "float.3", "rjust", None),
        ("c", 25, 33, "float.3", "rjust", None),
        ("alpha", 34, 40, "float.3", "rjust", None),
        ("beta", 41, 47, "float.3", "rjust", None),
        ("gamma", 48, 54, "float.3", "rjust", None),
        ("sgroup", 56, 66, "string", "rjust", None),
        ("z", 67, 70, "integer", "rjust", None)]

class CRYST1(CRYSTn):
    __slots__ = []
    
    _name = "CRYST1"
    
class CRYST2(CRYSTn):
    __slots__ = []
    
    _name = "CRYST2"

class CRYST3(CRYSTn):
    __slots__ = []
    
    _name = "CRYST3"

class ORIGXn(PDBRecord):
    """The ORIGXn (n = 1, 2, or 3) records present the transformation from
    the orthogonal coordinates contained in the entry to the submitted
    coordinates.
    """
    __slots__ = []
    
    _field_list = [
        ("o[n][1]", 11, 20, "float.6", "rjust", None),
        ("o[n][2]", 21, 30, "float.6", "rjust", None),
        ("o[n][3]", 31, 40, "float.6", "rjust", None),
        ("t[n]", 46, 55, "float.5", "rjust", None)]

class ORIGX1(ORIGXn):
    __slots__ = []
    
    _name = "ORIGX1"

class ORIGX2(ORIGXn):
    __slots__ = []
    
    _name = "ORIGX2"

class ORIGX3(ORIGXn):
    __slots__ = []
    
    _name = "ORIGX3"
    
class SCALEn(PDBRecord):
    """The SCALEn (n = 1, 2, or 3) records present the transformation from
    the orthogonal coordinates as contained in the entry to fractional
    crystallographic coordinates. Non-standard coordinate systems should
    be explained in the remarks.
    """
    __slots__ = []
    
    _field_list = [
        ("s[n][1]", 11, 20, "float.6", "rjust", None),
        ("s[n][2]", 21, 30, "float.6", "rjust", None),
        ("s[n][3]", 31, 40, "float.6", "rjust", None),
        ("u[n]", 46, 55, "float.5", "rjust", None)]
    
class SCALE1(SCALEn):
    __slots__ = []
    
    _name = "SCALE1"
        
class SCALE2(SCALEn):
    __slots__ = []
    
    _name = "SCALE2"

class SCALE3(SCALEn):
    __slots__ = []
    
    _name = "SCALE3"

class MTRIXn(PDBRecord):
    """The MTRIXn (n = 1, 2, or 3) records present transformations expressing
    non-crystallographic symmetry.
    """
    __slots__ = []
    
    _field_list = [
        ("serial", 8, 10, "integer", "rjust", None),
        ("s[n][1]", 11, 20, "float.6", "rjust", None),
        ("s[n][2]", 21, 30, "float.6", "rjust", None),
        ("s[n][3]", 31, 40, "float.6", "rjust", None),
        ("v[n]", 46, 55, "float.5", "rjust", None),
        ("iGiven", 60, 60, "integer", "rjust", None)]
    
class MTRIX1(MTRIXn):
    __slots__ = []
    
    _name = "MTRIX1"

class MTRIX2(MTRIXn):
    __slots__ = []
    
    _name = "MTRIX2"

class MTRIX3(MTRIXn):
    __slots__ = []
    
    _name = "MTRIX3"

class TVECT(PDBRecord):
    """The TVECT records present the translation vector for infinite
    covalently connected structures.
    """
    __slots__ = []
    
    _name = "TVECT "
    _field_list = [
        ("serial", 8, 10, "integer", "rjust", None),
        ("t[1]", 11, 20, "float.5", "rjust", None),
        ("t[2]", 21, 30, "float.5", "rjust", None),
        ("t[3]", 31, 40, "float.5", "rjust", None),
        ("text", 41, 70, "string", "rjust", None)]

## SECTION 9: Coordinate Selection
def ATOM_get_name(rec):
    """This should help older applications which do not use
    the element field of the ATOM record, these applications
    used column alignment to distinguish calcium (CA) from, say,
    a alpha-carbon (CA)
    """
    name    = rec.get("name")    or ""
    element = rec.get("element") or ""

    if len(element) == 2:
        name = name.ljust(4)[:4]
    else:
        l = len(name)
        if    l == 0:            name = "".ljust(4)
        elif  name[0].isdigit(): name = name.ljust(4)[:4]
        elif  l <  4:            name = " " + name.ljust(3)[:3]
        
    return name

class MODEL(PDBRecord):
    """The MODEL record specifies the model serial number when multiple
    structures are presented in a single coordinate entry, as is often
    the case with structures determined by NMR.
    """
    __slots__ = []
    
    _name = "MODEL "
    _field_list = [
        ("serial", 11, 14, "integer", "rjust", None)]

class ATOM(PDBRecord):
    """The ATOM records present the atomic coordinates for standard residues.
    They also present the occupancy and temperature factor for each atom.
    Heterogen coordinates use the HETATM record type. The element symbol
    is always present on each ATOM record; segment identifier and charge
    are optional.
    """
    __slots__ = []

    _name = "ATOM  "
    _field_list = [
        ("serial", 7, 11, "integer", "rjust", None),
        ("name", 13, 16, "string", "ljust.rstrip", ATOM_get_name),
        ("altLoc", 17, 17, "string", "rjust", None),
        ("resName", 18, 20, "string", "rjust", None),
        ("chainID", 22, 22, "string", "rjust", None),
        ("resSeq", 23, 26, "integer", "rjust", None),
        ("iCode", 27, 27, "string", "rjust", None),
        ("x", 31, 38, "float.3", "rjust", None),
        ("y", 39, 46, "float.3", "rjust", None),
        ("z", 47, 54, "float.3", "rjust", None),
        ("occupancy", 55, 60, "float.2", "rjust", None),
        ("tempFactor", 61, 66, "float.2", "rjust", None),
        ("segID", 73, 76, "string", "rjust", None),
        ("element", 77, 78, "string", "rjust", None),
        ("charge", 79, 80, "string", "rjust", None)]

class ANISOU(PDBRecord):
    """The ANISOU records present the anisotropic temperature factors.
    Columns 7 - 27 and 73 - 80 are identical to the corresponding
    ATOM/HETATM record.
    """
    __slots__ = []
    
    _name = "ANISOU"
    _field_list = [
        ("serial", 7, 11, "integer", "rjust", None),
        ("name", 13, 16, "string", "ljust", ATOM_get_name),
        ("altLoc", 17, 17, "string", "rjust", None),
        ("resName", 18, 20, "string", "rjust", None),
        ("chainID", 22, 22, "string", "rjust", None),
        ("resSeq", 23, 26, "integer", "rjust", None),
        ("iCode", 27, 27, "string", "rjust", None),
        ("u[0][0]", 29, 35, "integer", "rjust", None),
        ("u[1][1]", 36, 42, "integer", "rjust", None),
        ("u[2][2]", 43, 49, "integer", "rjust", None),
        ("u[0][1]", 50, 56, "integer", "rjust", None),
        ("u[0][2]", 57, 63, "integer", "rjust", None),
        ("u[1][2]", 64, 70, "integer", "rjust", None),
        ("segID", 73, 76, "string", "rjust", None),
        ("element", 77, 78, "string", "rjust", None),
        ("charge", 79, 80, "string", "rjust", None)]

class HETATM(ATOM):
    """The HETATM records present the atomic coordinate records for atoms
    within "non-standard" groups. These records are used for water
    molecules and atoms presented in HET groups.
    """
    __slots__ = []
    
    _name = "HETATM"

class SIGATM(PDBRecord):
    """The SIGATM records present the standard deviation
    of atomic parameters as they appear in ATOM and HETATM records.
    Columns 7 - 27 and 73 - 80 are identical to the corresponding
    ATOM/HETATM record.
    """
    _name = "SIGATM"
    _field_list = [
        ("serial", 7, 11, "integer", "rjust", None),
        ("name", 13, 16, "string", "ljust", ATOM_get_name),
        ("altLoc", 17, 17, "string", "rjust", None),
        ("resName", 18, 20, "string", "rjust", None),
        ("chainID", 22, 22, "string", "rjust", None),
        ("resSeq", 23, 26, "integer", "rjust", None),
        ("iCode", 27, 27, "string", "rjust", None),
        ("sigX", 31, 38, "float.3", "rjust", None),
        ("sigY", 39, 46, "float.3", "rjust", None),
        ("sigZ", 47, 54, "float.3", "rjust", None),
        ("sigOccupancy", 55, 60, "float.2", "rjust", None),
        ("sigTempFactor", 61, 66, "float.2", "rjust", None),
        ("segID", 73, 76, "string", "rjust", None),
        ("element", 77, 78, "string", "rjust", None),
        ("charge", 79, 80, "string", "rjust", None)]

class SIGUIJ(PDBRecord):
    """The SIGUIJ records present the standard deviations of anisotropic
    temperature factors scaled by a factor of 10**4 (Angstroms**2). 
    Columns 7 - 27 and 73 - 80 are identical to the corresponding
    ATOM/HETATM record.
    """
    __slots__ = []
    
    _name = "SIGUIJ"
    _field_list = [
        ("serial", 7, 11, "integer", "rjust", None),
        ("name", 13, 16, "string", "ljust", ATOM_get_name),
        ("altLoc", 17, 17, "string", "rjust", None),
        ("resName", 18, 20, "string","rjust", None),
        ("chainID", 22, 22, "string", "rjust", None),
        ("resSeq", 23, 26, "integer", "rjust", None),
        ("iCode", 27, 27, "string", "rjust", None),
        ("sig[1][1]", 29, 35, "integer", "rjust", None),
        ("sig[2][2]", 36, 42, "integer", "rjust", None),
        ("sig[3][3]", 43, 49, "integer", "rjust", None),
        ("sig[1][2]", 50, 56, "integer", "rjust", None),
        ("sig[1][3]", 57, 63, "integer", "rjust", None),
        ("sig[2][3]", 64, 70, "integer", "rjust", None),
        ("segID", 73, 76, "string", "rjust", None),
        ("element", 77, 78, "string", "rjust", None),
        ("charge", 79, 80, "string", "rjust", None)]

class TER(PDBRecord):
    """The TER record indicates the end of a list of ATOM/HETATM records
    for a chain.
    """
    __slots__ = []
    
    _name = "TER   "
    _field_list = [
        ("serial", 7, 11, "integer", "rjust", None),
        ("resName", 18, 20, "string", "rjust", None),
        ("chainID", 22, 22, "string", "rjust", None),
        ("resSeq", 23, 26, "integer", "rjust", None),
        ("iCode", 27, 27, "string", "rjust", None)]
    
class ENDMDL(PDBRecord):
    """The ENDMDL records are paired with MODEL records to group individual
    structures found in a coordinate entry.
    """
    __slots__ = []
    
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
    residues). These records are generated by the PDB.
    """
    __slots__ = []
    
    _name = "CONECT"
    _field_list = [
        ("serial", 7, 11, "integer", "rjust", None),
        ("serialBond1", 12, 16, "integer", "rjust", None),
        ("serialBond2", 17, 21, "integer", "rjust", None),
        ("serialBond3", 22, 26, "integer", "rjust", None),
        ("serialBond4", 27, 31, "integer", "rjust", None),
        ("serialHydBond1", 32, 36, "integer", "rjust", None),
        ("serialHydBond2", 37, 41, "integer", "rjust", None),
        ("serialSaltBond1", 42, 46, "integer", "rjust", None),
        ("serialHydBond3", 47, 51, "integer", "rjust", None),
        ("serialHydBond4", 52, 56, "integer", "rjust", None),
        ("serialSaltBond2", 57, 61, "integer", "rjust", None)]

## SECTION 11: Bookkeeping Section
class MASTER(PDBRecord):
    """The MASTER record is a control record for bookkeeping. It lists the
    number of lines in the coordinate entry or file for selected record
    types.
    """
    __slots__ = []
    
    _name = "MASTER"
    _field_list = [
        ("numRemark", 11, 15, "integer", "rjust", None),
        ("O", 16, 20, "integer", "rjust", None),
        ("numHet", 21, 25, "integer", "rjust", None),
        ("numHelix", 26, 30, "integer", "rjust", None),
        ("numSheet", 31, 35, "integer", "rjust", None),
        ("numTurn", 36, 40, "integer", "rjust", None),
        ("numSite", 41, 45, "integer", "rjust", None),
        ("numXForm", 46, 50, "integer", "rjust", None),
        ("numCoord", 51, 55, "integer", "rjust", None),
        ("numTer", 56, 60, "integer", "rjust", None),
        ("numConect", 61, 65, "integer", "rjust", None),
        ("numSeq", 66, 70, "integer", "rjust", None)]

class END(PDBRecord):
    """The END record marks the end of the PDB file.
    """
    __slots__ = []
    
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
    (TITLE._name,  TITLE,  "mandatory"),
    (CAVEAT._name, CAVEAT, "optional"),
    (COMPND._name, COMPND, "mandatory"),
    (SOURCE._name, SOURCE, "mandatory"),
    (KEYWDS._name, KEYWDS, "mandatory"),
    (EXPDTA._name, EXPDTA, "mandatory"),
    (AUTHOR._name, AUTHOR, "mandatory"),
    (REVDAT._name, REVDAT, "mandatory"),
    (SPRSDE._name, SPRSDE, "optional"),
    (JRNL._name,   JRNL,   "optional"),
    (REMARK._name, REMARK, "optional"),
    (DBREF._name,  DBREF,  "optional"),
    (SEQADV._name, SEQADV, "optional"),
    (SEQRES._name, SEQRES, "optional"),
    (MODRES._name, MODRES, "optional"),
    (HET._name,    HET,    "optional"),
    (HETNAM._name, HETNAM, "optional"),
    (HETSYN._name, HETSYN, "optional"),
    (FORMUL._name, FORMUL, "optional"),
    (HELIX._name,  HELIX,  "optional"),
    (SHEET._name,  SHEET,  "optional"),
    (TURN._name,   TURN,   "optional"),
    (SSBOND._name, SSBOND, "optional"),
    (LINK._name,   LINK,   "optional"),
    (HYDBND._name, HYDBND, "optional"),
    (SLTBRG._name, SLTBRG, "optional"),
    (CISPEP._name, CISPEP, "optional"),
    (SITE._name,   SITE,   "optional"),
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
    (TVECT._name,  TVECT,  "optional"),
    (MODEL._name,  MODEL,  "optional"),
    (ATOM._name,   ATOM,   "optional"),
    (SIGATM._name, SIGATM, "optional"),
    (ANISOU._name, ANISOU, "optional"),
    (SIGUIJ._name, SIGUIJ, "optional"),
    (TER._name,    TER,    "optional"),
    (HETATM._name, HETATM, "optional"),
    (ENDMDL._name, ENDMDL, "optional"),
    (CONECT._name, CONECT, "optional"),
    (MASTER._name, MASTER, "mandatory"),
    (END._name, END,       "mandatory")
    ]


## END PDB RECORD DEFINITIONS
###############################################################################


class PDBFile(list):
    """Class for managing a PDB file.  This class inherits from a Python
    list object, and contains a list of PDBRecord objects.
    Load, save, edit, and create PDB files with this class.
    """
    def __setattr__(self, i, rec):
        assert isinstance(rec, PDBRecord)
        list.__setattr__(self, i, rec)

    def append(self, rec):
        assert isinstance(rec, PDBRecord)
        list.append(self, rec)

    def insert(self, i, rec):
        assert isinstance(rec, PDBRecord)
        list.insert(self, i, rec)

    def load_file(self, fil, update_cb = None):
        """Loads a PDB file from File object fil.
        """
        fil = OpenFile(fil, "r")
    
        ## get file size for update callbacks
        line_number    = 0
        percent_done   = 0
        fil_read_bytes = 0

        ## some file objects do not support seek/tell
        if hasattr(fil, "seek") and hasattr(fil, "tell"):
            try:
                fil.seek(0, 2)
                fil_size_bytes = fil.tell()
                fil.seek(0, 0)
            except:
                # this is a adverage file size ;)
                fil_size_bytes = 1304189
        else:
            fil_size_bytes = 1304189

        for ln in fil.readlines():
            line_number    += 1
            fil_read_bytes += len(ln)

            ## call update callback
            if update_cb != None:
                pdone = (fil_read_bytes * 100)/fil_size_bytes
                if pdone != percent_done and pdone <= 100:
                    percent_done = pdone
                    update_cb(percent_done)

            ## find the record data element for the given line
            ln = ln.rstrip()
            rname = ln[:6].ljust(6)
            
            try:
                pdb_record_class = PDBRecordMap[rname]
            except KeyError:
                debug("PDB parser: unknown record type: %s"%(rname))
                continue

            ## create/add/parse the record
            pdb_record = pdb_record_class()
            pdb_record.read(ln)
            self.append(pdb_record)

    def save_file(self, fil):
        """Saves the PDBFile object in PDB file format to File object fil.
        """
        fil = OpenFile(fil, "w")
        for pdb_record in self:
            fil.write(str(pdb_record) + "\n")
        fil.flush()

    def record_processor(self, processor, filter_func = None):
        """Iterates the PDB records in self, and searches for handling
        methods in the processor object for reading the objects.  There
        are several choices for methods names for the processor objects.
        """
        def call_processor(recs):
            if type(recs) == ListType:
                rec = recs[0]
            else:
                rec = recs

            ## check filter function, if given, to determine if the record
            ## should be processed
            if filter_func and not filter_func(rec):
                return

            ## form method names to search for
            name = rec._name.strip()
            process = "process_%s" % (name)
            preprocess = "preprocess_%s" % (name)

            ## call process handler for records
            if hasattr(processor, process):
                getattr(processor, process)(recs)
            elif hasattr(processor, "process_default"):
                getattr(processor, "process_default")(recs)

            ## call preprocessor and processor for records
            if hasattr(rec, "process"):
                pfunc = getattr(rec, "process")
                if hasattr(processor, preprocess):
                    pfunc_cb = getattr(processor, preprocess)
                    pfunc_cb(pfunc(recs))
                elif hasattr(processor, "preprocess_default"):
                    pfunc_cb = getattr(processor, "preprocess_default")
                    pfunc_cb(name, pfunc(recs))

        cont_list = []

        for i in range(len(self)):

            ## current record
            rec = self[i]

            ## next record
            try:
                nrec = self[i+1]
            except IndexError:
                nrec = None

            ## check next record for continuation
            if nrec and rec._name == nrec._name:
                if nrec.has_key("continuation"):
                    cont_list.append(rec)
                else:
                    if cont_list:
                        call_processor(cont_list)
                        cont_list = []
                    call_processor(rec)
            else:
                if cont_list:
                    if rec._name == cont_list[0]._name:
                        cont_list.append(rec)
                        call_processor(cont_list)
                    else:
                        call_processor(cont_list)
                        call_processor(rec)
                    cont_list = []
                else:
                    call_processor(rec)
                
        if cont_list:
            call_processor(cont_list)


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
            anisou_rec["u[0][0]"] = int(atm.U[0,0] * 10000.0)
            anisou_rec["u[1][1]"] = int(atm.U[1,1] * 10000.0)
            anisou_rec["u[2][2]"] = int(atm.U[2,2] * 10000.0)
            anisou_rec["u[0][1]"] = int(atm.U[0,1] * 10000.0)
            anisou_rec["u[0][2]"] = int(atm.U[0,2] * 10000.0)
            anisou_rec["u[1][2]"] = int(atm.U[1,2] * 10000.0)

        if atm.sig_U != None:
            siguij_rec = SIGUIJ()
            self.pdb_file.append(siguij_rec)
            atom_common(atom_rec, siguij_rec)
            siguij_rec["u[0][0]"] = int(atm.U[0,0] * 10000.0)
            siguij_rec["u[1][1]"] = int(atm.U[1,1] * 10000.0)
            siguij_rec["u[2][2]"] = int(atm.U[2,2] * 10000.0)
            siguij_rec["u[0][1]"] = int(atm.U[0,1] * 10000.0)
            siguij_rec["u[0][2]"] = int(atm.U[0,2] * 10000.0)
            siguij_rec["u[1][2]"] = int(atm.U[1,2] * 10000.0)


### <testing>
def mkc():

    def prnt(field, i1, i2, ftype, just):
        print '  {'                       + \
              field.ljust(15)             + \
              ('%d, '%(i1)).rjust(4)       + \
              ('%d, '%(i2)).rjust(4)       + \
              ('%s,'%(pftype)).ljust(24)  + \
              pjust                       + \
              '},'

    for name in PDBRecordMap.keys():
        print '  {' + \
              ('"%s",'%(name)).ljust(10) + \
              "%s_fields" % (name.strip().lower()) + \
              '},'
    
    for name, rec in PDBRecordMap.items():
        print "static struct PDBFieldDef %s_fields[] = {" % (
            name.strip().lower())

        for (field, i1, i2, ftype, just, junk) in rec._field_list:

            if ftype   == "string":
                pftype = "PDB_FIELD_TYPE_STRING"
            elif ftype == "integer":
                pftype = "PDB_FIELD_TYPE_INTEGER"
            elif ftype == "float.2":
                pftype = "PDB_FIELD_TYPE_FLOAT_2"
            elif ftype == "float.3":
                pftype = "PDB_FIELD_TYPE_FLOAT_3"
            elif ftype == "float.4":
                pftype = "PDB_FIELD_TYPE_FLOAT_4"
            elif ftype == "float.5":
                pftype = "PDB_FIELD_TYPE_FLOAT_5"
            elif ftype == "float.6":
                pftype = "PDB_FIELD_TYPE_FLOAT_6"

            if just.count("ljust") > 0:
                pjust = "PDB_FIELD_RJUST"
            else:
                pjust = "PDB_FIELD_LJUST"

            prnt('"%s",'%(field), i1, i2, pftype, pjust)

        prnt("NULL,", 0, 0, "PDB_FIELD_TYPE_STRING", "PDB_FIELD_RJUST")
        print '};'
        print

if __name__ == "__main__":
    mkc()
    
##     import sys
##     try:
##         path = sys.argv[1]
##     except IndexError:
##         print "usage: PDB.py <PDB file path>"
##         sys.exit(1)
##     pdbfil = PDBFile()
##     pdbfil.load_file(path)
##     pdbfil.save_file(sys.stdout)

### </testing>

