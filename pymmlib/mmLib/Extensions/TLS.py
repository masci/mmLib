## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Utility classes for loading, manipulating, and analyzing TLS parameters.
"""
from __future__ import generators
import re
import sys
import string
import fpformat

from mmLib.mmTypes   import *
from mmLib.AtomMath  import *
from mmLib.Structure import *
from mmLib.PDB       import *
from mmLib.Viewer    import *


class TLSError(Exception):
    """Base exception class for TLS module exceptions.
    """
    pass


class TLSFileFormatError(TLSError):
    """Raised when a file format error is encountered while loading a
    TLS group description file.
    """
    pass


class TLSGroupDesc(object):
    """Description of one TLS Group.
    """
    def __init__(self):
        self.name       = ""     # text name
        self.origin     = None   # array(3)
        self.range_list = []     # (chain1, res1, chain2, res2, selection)
        self.T          = None   # array(3,3)
        self.L          = None   # array(3,3)
        self.S          = None   # array(3,3)

    def set_name(self, name):
        """Sets the TLS group name.
        """
        self.name = name

    def set_origin(self, x, y, z):
        """Sets the TLS group origin of calculations.
        """
        self.origin = array((x, y, z), Float)

    def add_range(self, chain_id1, frag_id1, chain_id2, frag_id2, selection):
        """Adds a segment of residues to the TLS group.  Not too sure how to
        handle segments which span chains, so assert on that condition.
        """
        assert chain_id1==chain_id2
        self.range_list.append(
            (chain_id1, frag_id1, chain_id2, frag_id2, selection))

    def set_T(self, t11, t22, t33, t12, t13, t23):
        """Sets the T tensor from the component arguments.  Units are in
        square Angstroms.
        """
        self.T = array( [[t11, t12, t13],
                         [t12, t22, t23],
                         [t13, t23, t33]], Float)

    def set_L(self,  l11, l22, l33, l12, l13, l23):
        """Sets the L tensor from the component arguments.  Units are in
        square Radians.
        """
        self.L = array( [[l11, l12, l13],
                         [l12, l22, l23],
                         [l13, l23, l33]], Float)

    def set_L_deg2(self,  l11, l22, l33, l12, l13, l23):
        """Sets the L tensor from the component arguments.  Units are in
        square Degrees.
        """
        self.set_L(l11*DEG2RAD2, l22*DEG2RAD2, l33*DEG2RAD2,
                   l12*DEG2RAD2, l13*DEG2RAD2, l23*DEG2RAD2)

    def set_S(self, s2211, s1133, s12, s13, s23, s21, s31, s32):
        """Sets the S tensor from the component arguments.  Units are in
        Radians*Angstroms.  The trace of S is set to 0.0.
        """
        s22 = 2.0*(s2211)/3.0 + s1133/3.0
        s11 = s22 - s2211
        s33 = s11 - s1133

        self.S = array([[s11, s12, s13],
                        [s21, s22, s23],
                        [s31, s32, s33]])
        
    def set_S_deg(self, s2211, s1133, s12, s13, s23, s21, s31, s32):
        """Sets the S tensor from the component arguments.  Units are in
        Degrees*Angstroms.  The trace of S is set to 0.0.
        """
        self.set_S(s2211*DEG2RAD, s1133*DEG2RAD, s12*DEG2RAD,
                   s13*DEG2RAD,   s23*DEG2RAD,   s21*DEG2RAD,
                   s31*DEG2RAD,   s32*DEG2RAD)

    def set_tls_group(self, tls_group):
        """Sets the TLSGroupDesc tensor values from the TLSGroup instance.
        """
        self.set_origin(
            tls_group.origin[0],
            tls_group.origin[1],
            tls_group.origin[2])

        self.set_T(
            tls_group.T[0,0],
            tls_group.T[1,1],
            tls_group.T[2,2],
            tls_group.T[0,1],
            tls_group.T[0,2],
            tls_group.T[1,2])
            
        self.set_L(
            tls_group.L[0,0],
            tls_group.L[1,1],
            tls_group.L[2,2],
            tls_group.L[0,1],
            tls_group.L[0,2],
            tls_group.L[1,2])

        self.set_S(
            tls_group.S[1,1] - tls_group.S[0,0],
            tls_group.S[0,0] - tls_group.S[2,2],
            tls_group.S[0,1],
            tls_group.S[0,2],
            tls_group.S[1,2],
            tls_group.S[1,0],
            tls_group.S[2,0],
            tls_group.S[2,1])

    def is_null(self):
        """Returns True if the T,L,S tensors are not set, or are set
        with values of zero.
        """
        if self.T==None or self.L==None or self.S==None:
            return True
        return False

    def calc_tls_name(self):
        """Creates a name for the TLS group using the selected
        residue ranges.
        """
        listx = []
        for (chain_id1, frag_id1, chain_id2, frag_id2, sel) in self.range_list:
            listx.append("%s%s-%s%s %s" % (
                chain_id1, frag_id1, chain_id2, frag_id2, sel))
        return string.join(listx,';')

    def generate_tls_group(self, struct):
        """Returns a TLSGroup containing the appropriate atoms from the
        argument Structure object, and with the origin, T, L, S tensors
        set.
        """
        tls_group = TLSGroup()

        if self.name=="":
            tls_group.name = self.calc_tls_name()
        else:
            tls_group.name = self.name

        if self.origin!=None:
            tls_group.set_origin(self.origin[0],
                                 self.origin[1],
                                 self.origin[2])

        if self.T!=None:
            tls_group.set_T(self.T[0,0], self.T[1,1], self.T[2,2],
                            self.T[0,1], self.T[0,2], self.T[1,2])

        if self.L!=None:
            tls_group.set_L(self.L[0,0], self.L[1,1], self.L[2,2],
                            self.L[0,1], self.L[0,2], self.L[1,2])

        if self.S!=None:
            ## s2211, s1133, s12, s13, s23, s21, s31, s32)
            tls_group.set_S(
                self.S[1,1] - self.S[0,0],
                self.S[0,0] - self.S[2,2],
                self.S[0,1],
                self.S[0,2],
                self.S[1,2],
                self.S[1,0],
                self.S[2,0],
                self.S[2,1])
        
        for (chain_id1, frag_id1, chain_id2, frag_id2, sel) in self.range_list:

            chain1 = struct.get_chain(chain_id1)
            if chain1 == None:
                print "[ERROR] no chain id: %s" % (chain_id1)
                return None

            try:
                seg = chain1[frag_id1:frag_id2]
            except KeyError:
                print "[ERROR] unable to find segment: %s-%s" % (
                    frag_id1, frag_id2)
                return None

            for atm in seg.iter_atoms():
                tls_group.append(atm)
            
        return tls_group


class TLSFile(object):
    """Read/Write a TLS files containing one or more TLSGroupDesc
    objects.
    """
    def __init__(self):
        self.path          = None
        self.tls_desc_list = []
        self.file_format   = None

    def set_file_format(self, file_format):
        """Set a instance of a TLSFileFormat subclass to use for reading and
        writing TLS description files.
        """
        self.file_format = file_format
    
    def load(self, fil, path=None):
        """Load TLS description file using the current TLSFileFormat instance.
        """
        assert self.file_format!=None
        self.path = path
        self.tls_desc_list = self.file_format.load(fil)
        
    def save(self, fil, path=None):
        """Save TLS description file using the curent TLSFileFormat instance.
        """
        assert self.file_format!=None
        self.path = path
        self.file_format.save(fil, self.tls_desc_list)

    def generate_tls_group_list(self, struct):
        """Returns a list of TLSGroup instances constructed from the Structure
        instance argument and the TLSGroupDesc instances contained in this
        class.
        """
        tls_group_list = []

        for tls_desc in self.tls_desc_list:
            tls_group = tls_desc.generate_tls_group(struct)
            if tls_group!=None:
                tls_group_list.append(tls_group)

        return tls_group_list

    
class TLSFileFormat(object):
    """Base class for TLS file types.
    """
    def load_supported(self):
        """Returns True if file loading is supported, otherwise returns False.
        """
        return False

    def save_supported(self):
        """Return True if file saving is supported, otherwise returns False.
        """
        return False
    
    def load(self, fil):
        """Returns a list containing one TLSGroupDesc object for each TLS group
        description found in the file.
        """
        pass

    def save(self, fil, tls_desc_list):
        """Writes the list of TLSGroupDesc object to the given file.
        """
        pass


class TLSFileFormatPDB(TLSFileFormat):
    """Reads TLS descriptions from the REMARK records in PDB files.  These
    records are only written by REFMAC5.
    """
    pdb_regex_dict = {
        "group": re.compile(
        "\s*TLS GROUP :\s+(\d+)\s*$"),

        "range": re.compile(
        "\s*RESIDUE RANGE :\s+(\w)\s+(\w+)\s+(\w)\s+(\w+)\s*$"),

        "origin": re.compile(
        "\s*ORIGIN\s+FOR\s+THE\s+GROUP\s+[(]A[)]:\s+(\S+)"\
        "\s+(\S+)\s+(\S+)\s*$"),

        "t11_t22": re.compile(
        "\s*T11:\s*(\S+)\s+T22:\s*(\S+)\s*$"),

        "t33_t12": re.compile(
        "\s*T33:\s*(\S+)\s+T12:\s*(\S+)\s*$"),

        "t13_t23": re.compile(
        "\s*T13:\s*(\S+)\s+T23:\s*(\S+)\s*$"),

        "l11_l22": re.compile(
        "\s*L11:\s*(\S+)\s+L22:\s*(\S+)\s*$"),

        "l33_l12": re.compile(
        "\s*L33:\s*(\S+)\s+L12:\s*(\S+)\s*$"),

        "l13_l23": re.compile(
        "\s*L13:\s*(\S+)\s+L23:\s*(\S+)\s*$"),

        "s11_s12_s13": re.compile(
        "\s*S11:\s*(\S+)\s+S12:\s*(\S+)\s+S13:\s*(\S+)\s*$"),

        "s21_s22_s23": re.compile(
        "\s*S21:\s*(\S+)\s+S22:\s*(\S+)\s+S23:\s*(\S+)\s*$"),

        "s31_s32_s33": re.compile(
        "\s*S31:\s*(\S+)\s+S32:\s*(\S+)\s+S33:\s*(\S+)\s*$")
        }

    def load_supported(self):
        return True

    def load(self, fil):
        self.tls_desc      = None
        self.tls_scrap     = {}
        self.tls_desc_list = []
       
        pdb_file = PDBFile()
        pdb_file.load_file(fil)
        pdb_file.record_processor(self)

        return self.tls_desc_list

    def complete_T(self):
        for key in ("t11","t22","t33","t12","t13","t23"):
            if not self.tls_scrap.has_key(key):
                return
        try:
            self.tls_desc.set_T(
                self.tls_scrap["t11"],
                self.tls_scrap["t22"],
                self.tls_scrap["t33"],
                self.tls_scrap["t12"],
                self.tls_scrap["t13"],
                self.tls_scrap["t23"])
        except AttributeError:
            raise TLSFileFormatError()

    def complete_L(self):
        for key in ("l11","l22","l33","l12","l13","l23"):
            if not self.tls_scrap.has_key(key):
                return
        try:
            self.tls_desc.set_L_deg2(
                self.tls_scrap["l11"],
                self.tls_scrap["l22"],
                self.tls_scrap["l33"],
                self.tls_scrap["l12"],
                self.tls_scrap["l13"],
                self.tls_scrap["l23"])
        except AttributeError:
            raise TLSFileFormatError()
        
    def complete_S(self):
        for key in ("s11","s22","s33","s12","s13","s23","s21","s31","s32"):
            if not self.tls_scrap.has_key(key):
                return
            
        ## s2211, s1133, s12, s13, s23, s21, s31, s32
        try:
            self.tls_desc.set_S_deg(
                self.tls_scrap["s22"] - self.tls_scrap["s11"],
                self.tls_scrap["s11"] - self.tls_scrap["s33"],
                self.tls_scrap["s12"],
                self.tls_scrap["s13"],
                self.tls_scrap["s23"],
                self.tls_scrap["s21"],
                self.tls_scrap["s31"],
                self.tls_scrap["s32"])
        except AttributeError:
            raise TLSFileFormatError()
            
    def process_REMARK(self, rec):
        """Callback for the PDBFile parser for REMARK records.  If the
        REMARK records contain TLS group information, then it is
        extracted and added to the TLSGroups list.
        """
        ## TLS REMARKs use remarkNum 3
        if rec.get("remarkNum", 0)!=3:
            return

        ## no text == no tls info
        if not rec.has_key("text"):
            return

        ## attempt to extract information from the text
        text = rec["text"]
        for (re_key, re_tls) in self.pdb_regex_dict.items():
            mx = re_tls.match(text)
            if mx != None:
                break

        ## no match
        if mx == None:
            return

        if re_key == "group":
            self.tls_desc = TLSGroupDesc()
            self.tls_desc_list.append(self.tls_desc)
            self.tls_scrap = {}
        
        elif re_key == "origin":
            (x, y, z) = mx.groups()
            
            try:
                self.tls_desc.set_origin(float(x), float(y), float(z))
            except AttributeError:
                raise TLSFileFormatError()
            except ValueError:
                raise TLSFileFormatError()

        elif re_key == "range":
            (chain_id1, frag_id1, chain_id2, frag_id2) = mx.groups()
            try:
                self.tls_desc.add_range(
                    chain_id1, frag_id1, chain_id2, frag_id2, "")
            except AttributeError:
                raise TLSFileFormatError()

        elif re_key == "t11_t22":
            (t11, t22) = mx.groups()
            try:
                self.tls_scrap["t11"] = float(t11)
                self.tls_scrap["t22"] = float(t22)
            except ValueError:
                raise TLSFileFormatError()
            self.complete_T()
            
        elif re_key == "t33_t12":
            (t33, t12) = mx.groups()
            try:
                self.tls_scrap["t33"] = float(t33)
                self.tls_scrap["t12"] = float(t12)
            except ValueError:
                raise TLSFileFormatError()
            self.complete_T()
            
        elif re_key == "t13_t23":
            (t13, t23) = mx.groups()
            try:
                self.tls_scrap["t13"] = float(t13)
                self.tls_scrap["t23"] = float(t23)
            except ValueError:
                raise TLSFileFormatError()
            self.complete_T()
            
        elif re_key == "l11_l22":
            (l11, l22) = mx.groups()
            try:
                self.tls_scrap["l11"] = float(l11)
                self.tls_scrap["l22"] = float(l22)
            except ValueError:
                raise TLSFileFormatError()
            self.complete_L()

        elif re_key == "l33_l12":
            (l33, l12) = mx.groups()
            try:
                self.tls_scrap["l33"] = float(l33)
                self.tls_scrap["l12"] = float(l12)
            except ValueError:
                raise TLSFileFormatError()
            self.complete_L()
            
        elif re_key == "l13_l23":
            (l13, l23) = mx.groups()
            try:
                self.tls_scrap["l13"] = float(l13)
                self.tls_scrap["l23"] = float(l23)
            except ValueError:
                raise TLSFileFormatError()
            self.complete_L()

        elif re_key == "s11_s12_s13":
            (s11, s12, s13) = mx.groups()
            try:
                self.tls_scrap["s11"] = float(s11)
                self.tls_scrap["s12"] = float(s12)
                self.tls_scrap["s13"] = float(s13)
            except ValueError:
                raise TLSFileFormatError()
            self.complete_S()

        elif re_key == "s21_s22_s23":
            (s21, s22, s23) = mx.groups()
            try:
                self.tls_scrap["s21"] = float(s21)
                self.tls_scrap["s22"] = float(s22)
                self.tls_scrap["s23"] = float(s23)
            except ValueError:
                raise TLSFileFormatError()
            self.complete_S()

        elif re_key == "s31_s32_s33":
            (s31, s32, s33) = mx.groups()
            try:
                self.tls_scrap["s31"] = float(s31)
                self.tls_scrap["s32"] = float(s32)
                self.tls_scrap["s33"] = float(s33)
            except ValueError:
                raise TLSFileFormatError()
            self.complete_S()


class TLSFileFormatTLSOUT(TLSFileFormat):
    """Read/Write REFMAC5 TLSIN/TLSOUT files.
    """
    tlsout_regex_dict = {
        "group": re.compile(
        "(?:^TLS\s*$)|(?:^TLS\s+(.*)$)"),
        
        "range": re.compile(
        "^RANGE\s+[']([A-Z])\s*([0-9A-Z.]+)\s*[']\s+"\
        "[']([A-Z])\s*([0-9A-Z.]+)\s*[']\s*(\w*).*$"),
        
        "origin": re.compile(
        "^ORIGIN\s+(\S+)\s+(\S+)\s+(\S+).*$"),
        
        "T": re.compile(
        "^T\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+).*$"),
        
        "L": re.compile(
        "^L\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+).*$"),

        "S": re.compile(
        "^S\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)"\
        "\s+(\S+)\s+(\S+).*$")
        }

    def convert_frag_id_load(self, frag_id):
        """Converts the residue sequence code to a mmLib fragment_id.
        """
        if len(frag_id)==0:
            return ""
        if frag_id[-1] == ".":
            frag_id = frag_id[:-1]
        return frag_id

    def convert_frag_id_save(self, frag_id):
        """Converts a mmLib fragment_id to a TLSOUT fragment id.
        """
        if len(frag_id)==0:
            return "."
        if frag_id[-1] in string.digits:
            frag_id += "."
        return frag_id
        
    def load_supported(self):
        return True
    def save_supported(self):
        return True

    def load(self, fil):
        tls_desc_list = []
        tls_desc      = None
        
        for ln in fil.readlines():
            ln = ln.rstrip()

            ## search all regular expressions for a match
            for (re_key, re_tls) in self.tlsout_regex_dict.items():
                mx = re_tls.match(ln)
                if mx != None:
                    break

            ## no match was found for the line
            if mx==None:
                continue
            
            ## do not allow a match if tls_info == None, because then
            ## somehow we've missed the TLS group begin line
            if tls_desc==None and re_key!="group":
                raise TLSFileFormatError()

            if re_key == "group":
                tls_desc = TLSGroupDesc()
                tls_desc_list.append(tls_desc)

                if mx.group(1)!=None:
                    tls_desc.set_name(mx.group(1))

            elif re_key == "origin":
                (x, y, z) = mx.groups()
                try:
                    tls_desc.set_origin(float(x), float(y), float(z))
                except ValueError:
                    raise TLSFileFormatError()

            elif re_key == "range":
                (chain_id1, frag_id1, chain_id2, frag_id2, sel) = mx.groups()

                frag_id1 = self.convert_frag_id_load(frag_id1)
                frag_id2 = self.convert_frag_id_load(frag_id2)

                tls_desc.add_range(
                    chain_id1, frag_id1, chain_id2, frag_id2, sel)

            elif re_key == "T":
                ## REFMAC ORDER: t11 t22 t33 t12 t13 t23
                (t11, t22, t33, t12, t13, t23) = mx.groups()

                try:
                    tls_desc.set_T(float(t11), float(t22), float(t33),
                                   float(t12), float(t13), float(t23))
                except ValueError:
                    raise TLSFileFormatError()

            elif re_key == "L":
                ## REFMAC ORDER: l11 l22 l33 l12 l13 l23
                (l11, l22, l33, l12, l13, l23) = mx.groups()

                try:
                    tls_desc.set_L_deg2(float(l11), float(l22), float(l33),
                                        float(l12), float(l13), float(l23))
                except ValueError:
                    raise TLSFileFormatError()
                    
            elif re_key == "S":
                ## REFMAC ORDER:
                ## <S22 - S11> <S11 - S33> <S12> <S13> <S23> <S21> <S31> <S32>
                (s2211, s1133, s12, s13, s23, s21, s31, s32) = mx.groups()

                try:
                    tls_desc.set_S_deg(float(s2211), float(s1133), float(s12),
                                       float(s13),   float(s23),   float(s21),
                                       float(s31),   float(s32))
                except ValueError:
                    raise TLSFileFormatError()

        return tls_desc_list

    def ft8(self, *fx):
        """Converts the fx float tuple to a TLSOUT tensor output string.
        """
        strx = ""
        for x in fx:
            strx += fpformat.fix(x, 4).rjust(8)
        return strx

    def tlsout_tls_desc(self, tls_desc):
        """Converts TLSGroupDesc instance to a multi-line string format
        ready to write to a TLSOUT file.
        """        
        listx = []

        if tls_desc.name!="":
            listx.append("TLS %s" % (tls_desc.name))
        else:
            listx.append("TLS")

        for (chain_id1, frag_id1,
             chain_id2, frag_id2, sel) in tls_desc.range_list:

            frag_id1 = self.convert_frag_id_save(frag_id1)
            frag_id2 = self.convert_frag_id_save(frag_id2)

            listx.append("RANGE  '%s%s' '%s%s' %s" % (
                chain_id1, frag_id1.rjust(5),
                chain_id2, frag_id2.rjust(5), sel))

        if tls_desc.origin!=None:
            listx.append("ORIGIN %8.4f %8.4f %8.4f" % (
                tls_desc.origin[0], tls_desc.origin[1], tls_desc.origin[2]))

        if tls_desc.T!=None:
            ## REFMAC ORDER: t11 t22 t33 t12 t13 t23
            x = self.ft8(tls_desc.T[0,0], tls_desc.T[1,1], tls_desc.T[2,2],
                         tls_desc.T[0,1], tls_desc.T[0,2], tls_desc.T[1,2])
            listx.append("T   %s" % (x))

        if tls_desc.L!=None:
            ## REFMAC ORDER: l11 l22 l33 l12 l13 l23
            x = self.ft8(
                tls_desc.L[0,0] * RAD2DEG2,
                tls_desc.L[1,1] * RAD2DEG2,
                tls_desc.L[2,2] * RAD2DEG2,
                tls_desc.L[0,1] * RAD2DEG2,
                tls_desc.L[0,2] * RAD2DEG2,
                tls_desc.L[1,2] * RAD2DEG2)

            listx.append("L   %s" % (x))

        if tls_desc.S!=None:
            ## REFMAC ORDER:
            ## <S22 - S11> <S11 - S33> <S12> <S13> <S23> <S21> <S31> <S32>
            x = self.ft8(
                (tls_desc.S[1,1] - tls_desc.S[0,0]) * RAD2DEG,
                (tls_desc.S[0,0] - tls_desc.S[2,2]) * RAD2DEG,
                tls_desc.S[0,1] * RAD2DEG,
                tls_desc.S[0,2] * RAD2DEG,
                tls_desc.S[1,2] * RAD2DEG,
                tls_desc.S[1,0] * RAD2DEG,
                tls_desc.S[2,0] * RAD2DEG,
                tls_desc.S[2,1] * RAD2DEG)
            
            listx.append("S   %s" % (x))

        return string.join(listx, "\n")

    def save(self, fil, tls_desc_list):
        ## with this line, the tensor components will be read by some
        ## programs in a different order, mangling the tensor values
        fil.write("REFMAC\n\n")
        
        for tls_desc in tls_desc_list:
            tls_desc_str = self.tlsout_tls_desc(tls_desc)
            fil.write(tls_desc_str + "\n")
        

class TLSGroup(AtomList):
    """A subclass of AtomList implementing methods for performing TLS
    calculations on the contained Atom instances.
    """
    def __init__(self, *args):
        AtomList.__init__(self, *args)

        self.name           = "" 
        self.origin         = zeros(3, Float)
        self.T              = array([[0.0, 0.0, 0.0],
                                     [0.0, 0.0, 0.0],
                                     [0.0, 0.0, 0.0]], Float)
        self.L              = array([[0.0, 0.0, 0.0],
                                     [0.0, 0.0, 0.0],
                                     [0.0, 0.0, 0.0]], Float)
        self.S              = array([[0.0, 0.0, 0.0],
                                     [0.0, 0.0, 0.0],
                                     [0.0, 0.0, 0.0]], Float)

    def str_old(self):
        tstr  = "TLS %s\n" % (self.name)

        tstr += "ORIGIN  %f %f %f\n" % (
            self.origin[0], self.origin[1], self.origin[2])

        tstr += "T %f %f %f %f %f %f\n" % (
            self.T[0,0], self.T[1,1], self.T[2,2], self.T[0,1], self.T[0,2],
            self.T[1,2])

        tstr += "L %f %f %f %f %f %f\n" % (
            self.L[0,0], self.L[1,1], self.L[2,2], self.L[0,1], self.L[0,2],
            self.L[1,2])

        tstr += "S %f %f %f %f %f %f %f %f\n" % (
            self.S[1,1]-self.S[0,0], self.S[0,0]-self.S[2,2], self.S[0,1],
            self.S[0,2], self.S[1,2], self.S[1,0], self.S[2,0], self.S[2,1])

        return tstr
    
    def set_origin(self, x, y, z):
        """Sets the x, y, z components of the TLS origin vector.
        """
        self.origin = array([x, y, z])

    def set_T(self, t11, t22, t33, t12, t13, t23):
        """Sets the components of the symmetric T tensor.  Units in square
        Angstroms.
        """
        self.T = array([[t11, t12, t13],
                        [t12, t22, t23],
                        [t13, t23, t33]])

    def set_L(self, l11, l22, l33, l12, l13, l23):
        """Sets the components of the symmetric L tensor from arguments.
        Units should be in square radians.
        """
        self.L = array([[l11, l12, l13],
                        [l12, l22, l23],
                        [l13, l23, l33]])

    def set_S(self, s2211, s1133, s12, s13, s23, s21, s31, s32):
        """Sets the componets of the asymmetric S tenssor.  The trace
        of the S tensor is set with the standard convention of
        the Trace(S) = 0.  Units in Radians*Angstroms.
        """
        s22 = 2.0*(s2211)/3.0 + s1133/3.0
        s11 = s22 - s2211
        s33 = s11 - s1133

        self.S = array([[s11, s12, s13],
                        [s21, s22, s23],
                        [s31, s32, s33]])

    def is_null(self):
        """Returns True if the T,L,S tensors are not set, or are set
        with values of zero.
        """
        if allclose(trace(self.T), 0.0) or allclose(trace(self.L), 0.0):
            return True
        return False

    def check_valid_model(self):
        """Returns True if the TLS model is mathmatically valid.  This
        requires the eigenvalues of the L and T tensor to be positive,
        and the Utls calculated anisotropic ADPs must must also have
        positive eigenvalues.
        """
        small = 1e-10

        if min(eigenvalues(self.L))<=small:
            return False

        if min(eigenvalues(self.T))<=small:
            return False

        T_red = self.calc_COR()["rT'"]
        if min(eigenvalues(T_red))<=small:
            return False

        for atm, Utls in self.iter_atm_Utls():
            if min(eigenvalues(Utls))<=small:
                print "Utls EEK"
                return False

        return True

    def calc_TLS_least_squares_fit(self):
        """Perform a least-squares fit of the atoms contained in self
        to the three TLS tensors: self.T, self.L, and self.S using the
        origin given by self.origin.
        """

        ## use label indexing to avoid confusion!
        T11, T22, T33, T12, T13, T23, L11, L22, L33, L12, L13, L23, \
        S11, S22, S33, S12, S13, S23, S21, S31, S32 = (
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)

        A = zeros((len(self)*6, 21), Float)
        b = zeros(len(self) * 6,  Float)

        i = -1
        for atm in self:
            i += 1

            ## set x, y, z as the vector components from the TLS origin
            x, y, z = atm.position - self.origin

            ## indecies of the components of U
            u11 =   i * 6
            u22 = u11 + 1
            u33 = u11 + 2
            u12 = u11 + 3
            u13 = u11 + 4
            u23 = u11 + 5

            ## set the B vector
            Uatm   = atm.get_U()            
            b[u11] = Uatm[0,0]
            b[u22] = Uatm[1,1]
            b[u33] = Uatm[2,2]
            b[u12] = Uatm[0,1]
            b[u13] = Uatm[0,2]
            b[u23] = Uatm[1,2]

            ## C Matrix
            xx = x*x
            yy = y*y
            zz = z*z

            xy = x*y
            xz = x*z
            yz = y*z

            A[u11, T11] = 1.0
            A[u11, L22] = zz
            A[u11, L33] = yy
            A[u11, L23] = -2.0*yz
            A[u11, S31] = -2.0*y
            A[u11, S21] = 2.0*z

            A[u22, T22] = 1.0
            A[u22, L11] = zz
            A[u22, L33] = xx
            A[u22, L13] = -2.0*xz
            A[u22, S12] = -2.0*z
            A[u22, S32] = 2.0*x

            A[u33, T33] = 1.0
            A[u33, L11] = yy
            A[u33, L22] = xx
            A[u33, L12] = -2.0*xy
            A[u33, S23] = -2.0*x
            A[u33, S13] = 2.0*y

            A[u12, T12] = 1.0
            A[u12, L33] = -xy
            A[u12, L23] = xz
            A[u12, L13] = yz
            A[u12, L12] = -zz
            A[u12, S11] = -z
            A[u12, S22] = z
            A[u12, S31] = x
            A[u12, S32] = -y

            A[u13, T13] = 1.0
            A[u13, L22] = -xz
            A[u13, L23] = xy
            A[u13, L13] = -yy
            A[u13, L12] = yz
            A[u13, S11] = y
            A[u13, S33] = -y
            A[u13, S23] = z
            A[u13, S21] = -x

            A[u23, T23] = 1.0
            A[u23, L11] = -yz
            A[u23, L23] = -xx
            A[u23, L13] = xy
            A[u23, L12] = xz
            A[u23, S22] = -x
            A[u23, S33] = x
            A[u23, S12] = y
            A[u23, S13] = -z

        (C, resids, rank, s) = linear_least_squares(A, b)

        ## verify the trase of S is zero
        assert allclose(C[S11]+C[S22]+C[S33], 0.0)

        self.T = array([ [ C[T11], C[T12], C[T13] ],
                         [ C[T12], C[T22], C[T23] ],
                         [ C[T13], C[T23], C[T33] ] ], Float)

        self.L = array([ [ C[L11], C[L12], C[L13] ],
                         [ C[L12], C[L22], C[L23] ],
                         [ C[L13], C[L23], C[L33] ] ], Float)

        self.S = array([ [ C[S11], C[S12], C[S13] ],
                         [ C[S21], C[S22], C[S23] ],
                         [ C[S31], C[S32], C[S33] ] ], Float)

    def calc_Utls(self, T, L, S, vec):
        """Returns the calculated value for the anisotropic U tensor for
        the given position.
        """
        x, y, z = vec

        xx = x*x
        yy = y*y
        zz = z*z

        xy = x*y
        yz = y*z
        xz = x*z

        u11 = T[0,0] \
              + L[1,1]*zz + L[2,2]*yy - 2.0*L[1,2]*yz \
              + 2.0*S[1,0]*z - 2.0*S[2,0]*y

        u22 = T[1,1] \
              + L[0,0]*zz + L[2,2]*xx - 2.0*L[2,0]*xz \
              - 2.0*S[0,1]*z + 2.0*S[2,1]*x

        u33 = T[2,2] \
              + L[0,0]*yy + L[1,1]*xx - 2.0*L[0,1]*xy \
              - 2.0*S[1,2]*x + 2.0*S[0,2]*y

        u12 = T[0,1] \
              - L[2,2]*xy + L[1,2]*xz + L[2,0]*yz - L[0,1]*zz \
              - S[0,0]*z + S[1,1]*z + S[2,0]*x - S[2,1]*y

        u13 = T[0,2] \
              - L[1,1]*xz + L[1,2]*xy - L[2,0]*yy + L[0,1]*yz \
              + S[0,0]*y - S[2,2]*y + S[1,2]*z - S[1,0]*x

        u23 = T[1,2] \
              - L[0,0]*yz - L[1,2]*xx + L[2,0]*xy + L[0,1]*xz \
              - S[1,1]*x + S[2,2]*x + S[0,1]*y - S[0,2]*z

        return array([[u11, u12, u13],
                      [u12, u22, u23],
                      [u13, u23, u33]])

    def iter_atm_Utls(self, T=None, L=None, S=None, o=None):
        """Iterates all the atoms in the TLS object, returning the 2-tuple
        (atm, U) where U is the calcuated U value from the current values
        of the TLS object's T,L,S, tensors and origin.
        """
        T      = T or self.T
        L      = L or self.L
        S      = S or self.S
        origin = o or self.origin
        
        for atm in self:
            Utls = self.calc_Utls(T, L, S, atm.position - origin)
            yield atm, Utls

    def shift_COR(self):
        """Shift the TLS group to the center of reaction.
        """
        calcs       = self.calc_COR()
        
        self.origin = calcs["COR"].copy()
        self.T      = calcs["T'"].copy()
        self.L      = calcs["L'"].copy()
        self.S      = calcs["S'"].copy()

        return calcs
        
    def calc_COR(self):
        """Calculate new tensors based on the center for reaction.
        This method returns a dictionary of the calculations:

        T^: T tensor in the coordinate system of L
        L^: L tensor in the coordinate system of L
        S^: S tensor in the coordinate system of L

        COR: Center of Reaction

        T',S',L': T,L,S tensors in origonal coordinate system
                  with the origin shifted to the center of reaction.
        """
        ## small_L is the smallest magnitude of L before it is considered 0.0
        small_L = 0.002 * DEG2RAD2

        tls_info = {}

        ## set the L tensor eigenvalues and eigenvectors
        (eval_L, evec_L) = eigenvectors(self.L)

        tls_info["L1_eigen_val"] = eval_L[0]
        tls_info["L2_eigen_val"] = eval_L[1]
        tls_info["L3_eigen_val"] = eval_L[2]
        
        tls_info["L1_eigen_vec"] = evec_L[0]
        tls_info["L2_eigen_vec"] = evec_L[1]
        tls_info["L3_eigen_vec"] = evec_L[2]

        ## transpose the original the evec_L so it can be used
        ## to rotate the other tensors
        evec_L = transpose(evec_L)
        
        ## carrot-L tensor (tensor WRT principal axes of L)
        cL      = zeros([3,3], Float)
        cL[0,0] = eval_L[0]
        cL[1,1] = eval_L[1]
        cL[2,2] = eval_L[2]

        tls_info["L^"] = cL

        ## carrot-T tensor (T tensor WRT principal axes of L)
        cT = matrixmultiply(
            matrixmultiply(transpose(evec_L), self.T), evec_L)

        tls_info["T^"] = cT

        ## carrot-S tensor (S tensor WRT principal axes of L)
        cS = matrixmultiply(
            matrixmultiply(transpose(evec_L), self.S), evec_L)

        ## correct for left-handed libration eigenvectors
        det = determinant(evec_L)
        if int(det) != 1:
            cS = -cS

        tls_info["S^"] = cS
        
        ## ^rho: the origin-shift vector in the coordinate system of L
        cL1122 = cL[1,1] + cL[2,2]
        cL2200 = cL[2,2] + cL[0,0]
        cL0011 = cL[0,0] + cL[1,1]

        if cL1122>small_L:
            crho0 = (cS[1,2]-cS[2,1]) / cL1122
        else:
            crho0 = 0.0

        if cL2200>small_L:
            crho1 = (cS[2,0]-cS[0,2]) / cL2200
        else:
            crho1 = 0.0

        if cL0011>small_L:
            crho2 = (cS[0,1]-cS[1,0]) / cL0011
        else:
            crho2 = 0.0

        crho = array([crho0, crho1, crho2], Float)

        tls_info["RHO^"] = crho
        
        ## rho: the origin-shift vector in orthogonal coordinates
        rho = matrixmultiply(evec_L, crho)
        
        tls_info["RHO"] = rho
        tls_info["COR"] = array(self.origin, Float) + rho

        ## set up the origin shift matrix PRHO WRT orthogonal axes
        PRHO = array([ [    0.0,  rho[2], -rho[1]],
                       [-rho[2],     0.0,  rho[0]],
                       [ rho[1], -rho[0],     0.0] ], Float)

        ## set up the origin shift matrix cPRHO WRT libration axes
        cPRHO = array([ [    0.0,  crho[2], -crho[1]],
                        [-crho[2],     0.0,  crho[0]],
                        [ crho[1], -crho[0],     0.0] ], Float)

        ## calculate tranpose of cPRHO, ans cS
        cSt = transpose(cS)
        cPRHOt = transpose(cPRHO)

        ## calculate S'^ = S^ + L^*pRHOt
        cSp = cS + matrixmultiply(cL, cPRHOt)
        tls_info["S'^"] = cSp

        ## L'^ = L^ = cL
        tls_info["L'^"] = cL

        ## calculate T'^ = cT + cPRHO*S^ + cSt*cPRHOt + cPRHO*cL*cPRHOt *
        cTp = cT + \
              matrixmultiply(cPRHO, cS) + \
              matrixmultiply(cSt, cPRHOt) + \
              matrixmultiply(matrixmultiply(cPRHO, cL), cPRHOt)
        tls_info["T'^"] = cTp

        ## transpose of PRHO and S
        PRHOt = transpose(PRHO)
        St = transpose(self.S)

        ## calculate S' = S + L*PRHOt
        Sp = self.S + matrixmultiply(self.L, PRHOt)
        tls_info["S'"] = Sp

        ## calculate T' = T + PRHO*S + St*PRHOT + PRHO*L*PRHOt
        Tp = self.T + \
             matrixmultiply(PRHO, self.S) + \
             matrixmultiply(St, PRHOt) + \
             matrixmultiply(matrixmultiply(PRHO, self.L), PRHOt)
        tls_info["T'"] = Tp

        ## L' is just L
        tls_info["L'"] = self.L.copy()

        ## now calculate the TLS motion description using 3 non
        ## intersecting screw axes, with one

        ## libration axis 1 shift in the L coordinate system        
        if cL[0,0]>small_L:
            cL1rho = array([0.0, -cSp[0,2]/cL[0,0], cSp[0,1]/cL[0,0]], Float)
        else:
            cL1rho = zeros(3, Float)

        ## libration axis 2 shift in the L coordinate system
        if cL[1,1]>small_L:
            cL2rho = array([cSp[1,2]/cL[1,1], 0.0, -cSp[1,0]/cL[1,1]], Float)
        else:
            cL2rho = zeros(3, Float)

        ## libration axis 2 shift in the L coordinate system
        if cL[2,2]>small_L:
            cL3rho = array([-cSp[2,1]/cL[2,2], cSp[2,0]/cL[2,2], 0.0], Float)
        else:
            cL3rho = zeros(3, Float)

        ## libration axes shifts in the origional orthogonal
        ## coordinate system
        tls_info["L1_rho"] = matrixmultiply(evec_L, cL1rho)
        tls_info["L2_rho"] = matrixmultiply(evec_L, cL2rho)
        tls_info["L3_rho"] = matrixmultiply(evec_L, cL3rho)

        ## calculate screw pitches (A*R / R*R) = (A/R)
        if cL[0,0]>small_L:
            tls_info["L1_pitch"] = cS[0,0]/cL[0,0]
        else:
            tls_info["L1_pitch"] = 0.0
            
        if cL[1,1]>small_L:
            tls_info["L2_pitch"] = cS[1,1]/cL[1,1]
        else:
            tls_info["L2_pitch"] = 0.0

        if cL[2,2]>small_L:
            tls_info["L3_pitch"] = cS[2,2]/cL[2,2]
        else:
            tls_info["L3_pitch"] = 0.0

        ## now calculate the reduction in T for the screw rotation axes
        cTred = cT.copy()

        for i in (0, 1, 2):
            for k in (0, 1, 2):
                if i==k:
                    continue
                if cL[k,k]>small_L:
                    cTred[i,i] -= (cS[k,i]**2) / cL[k,k]

        for i in (0, 1, 2):
            for j in (0, 1, 2):
                for k in (0, 1, 2):
                    if j==i:
                        continue
                    if cL[k,k]>small_L:
                        cTred[i,j] -= (cS[k,i]*cS[k,j]) / cL[k,k]

        ## rotate the newly calculated reduced-T tensor from the carrot
        ## coordinate system (coordinate system of L) back to the structure
        ## coordinate system
        tls_info["rT'"] = matrixmultiply(
            transpose(evec_L), matrixmultiply(cTred, evec_L))

        return tls_info

    def verify_equivalent_models(self):
        """XXX: Fixme
        """
        ### Verify that the the math we've just done is correct by
        ### comparing the original TLS calculated U tensors with the
        ### U tensors calculated from the COR
        T_cor = tls_info["T'"].copy()
        L_cor = tls_info["L'"].copy()
        S_cor = tls_info["S'"].copy()

        for atm in self:
            x  = atm.position - self.origin
            xp = atm.position - tls_info["COR"]

            Utls     = self.calc_Utls(self.T, self.L, self.S, x)
            Utls_cor = self.calc_Utls(T_cor, L_cor, S_cor, xp)

            assert allclose(Utls, Utls_cor)

    def calc_tls_info(self):
        """Calculates a number of statistics about the TLS group tensors,
        goodness of fit, various parameter averages, center of reaction
        tensors, etc...  If tls_info["valid_model"] is not None, then the TLS
        group parameters are invalid.
        """
        tls_info = self.calc_COR()

        ## EXPERIMENTAL DATA

        ## number of atoms
        tls_info["num_atoms"] = len(self)

        ## mean temp_factor/anisotropy from experimental data (PDB file)
        tls_info["exp_mean_temp_factor"] = self.calc_adv_temp_factor()
        tls_info["exp_mean_anisotropy"]  = self.calc_adv_anisotropy()

        ## TLS FIT
        tls_info["valid_model"] = self.check_valid_model()
        if not tls_info["valid_model"]:
            return tls_info

        ## goodness of fit 
        tls_info["R"] = self.calc_R()
        tls_info["mean_S"], tls_info["mean_S_sigma"] = self.calc_mean_S()
        tls_info["mean_dp2"], tls_info["mean_dp2_sigma"] = self.calc_mean_DP2()
        tls_info["sum_dp2"] = self.calc_sum_DP2()
        
        ## model temp factors
        n = 0
        mean_max_tf = 0.0
        mean_tf     = 0.0
        mean_aniso  = 0.0

        for atm, Utls in self.iter_atm_Utls():
            n += 1

            evals  = eigenvalues(Utls)
            max_ev = max(evals)
            min_ev = min(evals)

            mean_max_tf += u2b * max_ev
            mean_tf     += u2b * trace(Utls) / 3.0
            mean_aniso  += calc_anisotropy(Utls)

        tls_info["tls_mean_max_tf"] = mean_max_tf / float(n)
        tls_info["tls_mean_tf"]     = mean_tf     / float(n)
        tls_info["tls_mean_aniso"]  = mean_aniso  / float(n)        

        return tls_info

    def calc_R(self):
        """Calculate the R factor of U vs. Utls.
        """
        Rn = 0.0
        Rd = 0.0

        for atm, Utls in self.iter_atm_Utls():
            U = atm.get_U()

            Rn += abs(U[0,0] - Utls[0,0])
            Rn += abs(U[1,1] - Utls[1,1])
            Rn += abs(U[2,2] - Utls[2,2])
            Rn += abs(U[0,1] - Utls[0,1])
            Rn += abs(U[0,2] - Utls[0,2])
            Rn += abs(U[1,2] - Utls[1,2])

            Rd += abs(U[0,0])
            Rd += abs(U[1,1])
            Rd += abs(U[2,2])
            Rd += abs(U[0,1])
            Rd += abs(U[0,2])
            Rd += abs(U[1,2])
            
        return Rn / Rd

    def calc_mean_S(self):
        """Calculates the mean Suij of U vs. Utls.
        """
        num      = 0
        mean_S   = 0.0
        MSD      = 0.0
        
        S_dict = {}
        for atm, Utls in self.iter_atm_Utls():
            Uatm = atm.get_U()
            atm_S = calc_Suij(Uatm, Utls)
            S_dict[atm] = atm_S
            mean_S += atm_S
            num += 1

        mean_S = mean_S / float(num)

        for S in S_dict.values():
            MSD += (mean_S - S)**2

        MSD = MSD / float(num)

        return mean_S, math.sqrt(MSD)

    def calc_sum_DP2(self):
        """Calculates the sum of DP2(U,Utls) for all Atoms in the
        TLSGroup.
        """
        sum_DP2 = 0.0
        
        for atm, Utls in self.iter_atm_Utls():
            Uatm = atm.get_U()
            sum_DP2 += calc_DP2uij(Uatm, Utls)

        return sum_DP2
                
    def calc_mean_DP2(self):
        """Calculates the mean dP2 and standard deviation of U vs. Utls
        for every atom in the TLS group.
        Returns the 2-tuple. (mean DP2, rmsd)
        """
        num      = 0
        mean_DP2 = 0.0
        MSD      = 0.0
        
        DP2_dict = {}
        for atm, Utls in self.iter_atm_Utls():
            Uatm = atm.get_U()
            atm_DP2 = calc_DP2uij(Uatm, Utls)
            DP2_dict[atm] = atm_DP2 
            mean_DP2 += atm_DP2
            num += 1

        mean_DP2 = mean_DP2 / float(num)

        for DP2 in DP2_dict.values():
            MSD += (mean_DP2 - DP2)**2
        
        MSD = MSD / float(num)

        return mean_DP2, math.sqrt(MSD)


class TLSStructureAnalysis(object):
    """Algorithm object for rigid body searches on Structure objects.
    """
    def __init__(self, struct):
        self.struct = struct

    def iter_segments(self, chain, seg_len):
        """This iteratar yields a series of Segment objects of width
        seg_len.  The start at the beginning Fragment of the Chain,
        and the start point walks the chain one Fragment at a time
        until there are not enough Fragments left to cut Segments of
        seg_width.
        """
        chain_len = len(chain)

        for i in range(chain_len):
            start = i
            end   = i + seg_len

            if end>chain_len:
                break

            yield chain[start:end]

    def atom_filter(self, atm, **args):
        use_side_chains         = args.get("use_side_chains", True)
        include_hydrogens       = args.get("include_hydrogens", False)
        include_frac_occupancy  = args.get("include_frac_occupancy", False)
        include_single_bond     = args.get("include_single_bond", True)
        
        ## omit atoms which are not at full occupancy
        if not include_frac_occupancy and atm.occupancy<1.0:
            return False

        ## omit hydrogens
        if not include_hydrogens and atm.element=="H":
            return False

        ## omit side chain atoms
        if not use_side_chains and atm.name not in ("C", "N", "CA", "O"):
            return False

        ## omit atoms with a single bond 
        if not include_single_bond and len(atm.bond_list)<=1:
            return False
                        
        return True
        
    def iter_fit_TLS_segments(self, **args):
        """Run the algorithm to fit TLS parameters to segments of the
        structure.  This method has many options, which are outlined in
        the source code for the method.  This returns a list of dictionaries
        containing statistics on each of the fit TLS groups, the residues
        involved, and the TLS object itself.
        """
        
        ## arguments
        chain_ids               = args.get("chain_ids", None)
        origin                  = args.get("origin_of_calc")
        residue_width           = args.get("residue_width", 6)
        use_side_chains         = args.get("use_side_chains", True)
        include_hydrogens       = args.get("include_hydrogens", False)
        include_frac_occupancy  = args.get("include_frac_occupancy", False)
        include_single_bond     = args.get("include_single_bond", True)
        
        for chain in self.struct.iter_chains():

            ## skip some chains
            if chain_ids!=None and chain.chain_id not in chain_ids:
                continue

            ## don't bother with non-biopolymers and small chains
            if chain.count_amino_acids()<residue_width:
                continue
            
            for segment in self.iter_segments(chain, residue_width):
                frag_id1 = segment[0].fragment_id
                frag_id2 = segment[-1].fragment_id
                name     = "%s-%s" % (frag_id1, frag_id2)

                ## create the TLSGroup
                tls_group = TLSGroup()

                ## add atoms into the TLSGroup
                ## filter the atoms going into the TLS group                
                for atm in segment.iter_atoms():
                    if self.atom_filter(atm, **args):
                        tls_group.append(atm)

                ## check for enough atoms(parameters) after atom filtering
                if len(tls_group)<20:
                    continue
                
                ## calculate tensors and print
                tls_group.calc_TLS_least_squares_fit()
                tls_group.shift_COR()
                tls_info = tls_group.calc_tls_info()

                ## check if the TLS model is valid
                if not tls_info["valid_model"]:
                    continue

                ## add additional information
                tls_info["name"]      = name
                tls_info["chain_id"]  = chain.chain_id
                tls_info["frag_id1"]  = frag_id1
                tls_info["frag_id2"]  = frag_id2
                tls_info["tls_group"] = tls_group
                tls_info["residues"]  = segment
                tls_info["segment"]   = segment
                    
                ## this TLS group passes all our tests -- yield it
                yield tls_info


    def fit_TLS_segments(self, **args):
        """Returns the list iterated by iter_fit_TLS_segments
        """
        tls_info_list = []
        for tls_info in self.iter_fit_TLS_segments(**args):
            tls_info_list.append(tls_info)
        return tls_info_list


###############################################################################
### GLViewer Rendering components for TLS Groups
###


def goodness_color(x):
    """x in range 0.0->1.0
    """
    if x<=0.0:
        return (0.0, 0.0, 0.0)

    r = math.sqrt(x)
    g = max(0.0, x**3)
    b = max(0.0, math.sin(2.0 * math.pi * x))

    return (r, g, b)


class GLTLSAtomList(GLAtomList):
    """OpenGL visualizations of TLS group atoms.
    """
    def __init__(self, **args):
        self.tls_group = args["tls_group"]        
        GLAtomList.__init__(self, **args)
        self.glo_init_properties(**args)
    
    def glo_install_properties(self):
        GLAtomList.glo_install_properties(self)

        ## Show/Hide
        self.glo_add_property(
            { "name":        "fan_visible",
              "desc":        "Show COR-Backbone Fan",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile" })

        self.glo_add_property(
            { "name":        "L1_animation_visible",
              "desc":        "Show L<sub>1</sub> Screw Animation",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     True,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L2_animation_visible",
              "desc":        "Show L<sub>2</sub> Screw Animation",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     True,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L3_animation_visible",
              "desc":        "Show L<sub>3</sub> Screw Animation",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     True,
              "action":      "redraw" })

        ## TLS
        self.glo_add_property(
            { "name":        "tls_color",
              "desc":        "TLS Group Color",
              "catagory":    "TLS",
              "type":        "enum_string",
              "default":     "Green",
              "enum_list":   self.gldl_color_list,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "fan_opacity",
              "desc":        "COR-Backbone Fan Opacity",
              "catagory":    "TLS",
              "type":        "float",
              "range":       PROP_OPACITY_RANGE,
              "default":     1.0,
              "action":      "recompile_fan" })
        self.glo_add_property(
            { "name":        "L1_scale",
              "desc":        "Scale L<sub>1</sub> Rotation", 
              "catagory":    "TLS",
              "type":        "float",
              "default":     1.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L2_scale",
              "desc":        "Scale L<sub>2</sub> Rotation", 
              "catagory":    "TLS",
              "type":        "float",
              "default":     1.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L3_scale",
              "desc":        "Scale L<sub>3</sub> Rotation", 
              "catagory":    "TLS",
              "type":        "float",
              "default":     1.0,
              "action":      "redraw" })

        ## TLS Analysis
        self.glo_add_property(
            { "name":        "COR",
              "desc":        "TLS Center of Reaction", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "T",
              "desc":        "T<sup>COR</sup> Tensor (A<sup>2</sup>)",
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3,3)",
              "default":     zeros((3,3), Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "rT",
              "desc":        "T<sup>r</sup> Tensor (A<sup>2</sup>)",
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3,3)",
              "default":     zeros((3,3), Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L",
              "desc":        "L<sup>COR</sup> Tensor (DEG<sup>2</sup>)",
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3,3)",
              "default":     zeros((3,3), Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "S",
              "desc":        "S<sup>COR</sup> Tensor (A*DEG)",
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3,3)",
              "default":     zeros((3,3), Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_eigen_vec",
              "desc":        "L<sub>1</sub> Eigen Vector", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_eigen_vec",
              "desc":        "L<sub>2</sub> Eigen Vector", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_eigen_vec",
              "desc":        "L<sub>3</sub> Eigen Vector", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_eigen_val",
              "desc":        "L<sub>1</sub> Eigen Value", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_eigen_val",
              "desc":        "L<sub>2</sub> Eigen Value", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_eigen_val",
              "desc":        "L<sub>3</sub> Eigen Value", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_rho",
              "desc":        "L<sub>1</sub> Position from COR", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_rho",
              "desc":        "L<sub>2</sub> Position from COR", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_rho",
              "desc":        "L<sub>3</sub> Position from COR", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_pitch",
              "desc":        "L<sub>1</sub> Screw Pitch (A/DEG)", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_pitch",
              "desc":        "L<sub>2</sub> Screw Pitch (A/DEG)", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_pitch",
              "desc":        "L<sub>3</sub> Screw Pitch (A/DEG)", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })

        ## Simulation State
        self.glo_add_property(
            { "name":        "L1_rot",
              "desc":        "L<sub>1</sub> Rotation", 
              "catagory":    "TLS",
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L2_rot",
              "desc":        "L<sub>2</sub> Rotation", 
              "catagory":    "TLS",
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L3_rot",
              "desc":        "L<sub>3</sub> Rotation", 
              "catagory":    "TLS",
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })

    def gldl_install_draw_methods(self):
        GLAtomList.gldl_install_draw_methods(self)
        
        self.gldl_draw_method_install(
            { "name":                "fan",
              "func":                self.draw_fan,
              "visible_property":    "fan_visible",
              "opacity_property":    "fan_opacity",
              "recompile_action":    "recompile_fan" })

    def gldl_iter_multidraw_self(self):
        for draw_flag in GLAtomList.gldl_iter_multidraw_self(self):
            for draw_flag2 in self.gldl_iter_multidraw_animate():
                yield True
            
    def gldl_iter_multidraw_animate(self):
        """
        """
        ## optimization: if a rotation of 0.0 degrees was already
        ## drawn, then there is no need to draw it again
        zero_rot = False
        
        for Lx_axis, Lx_rho, Lx_pitch, Lx_rot, Lx_scale in (
            ("L1_eigen_vec", "L1_rho", "L1_pitch", "L1_rot", "L1_scale"),
            ("L2_eigen_vec", "L2_rho", "L2_pitch", "L2_rot", "L2_scale"),
            ("L3_eigen_vec", "L3_rho", "L3_pitch", "L3_rot", "L3_scale") ):

            if Lx_axis=="L1_eigen_vec" and \
               self.properties["L1_animation_visible"]==False:
                continue
            if Lx_axis=="L2_eigen_vec" and \
               self.properties["L2_animation_visible"]==False:
                continue
            if Lx_axis=="L3_eigen_vec" and \
               self.properties["L3_animation_visible"]==False:
                continue

            axis  = self.properties[Lx_axis]
            rho   = self.properties[Lx_rho]
            pitch = self.properties[Lx_pitch]
            rot   = self.properties[Lx_rot] * self.properties[Lx_scale]
            screw = axis * (rot * pitch)

            if allclose(rot, 0.0):
                if zero_rot:
                    continue
                zero_rot = True

            self.driver.glr_push_matrix()

            self.driver.glr_translate(-rho)
            self.driver.glr_rotate_axis(rot, axis)
            self.driver.glr_translate(rho + screw)
            yield True
            
            self.driver.glr_pop_matrix()

    def glal_iter_atoms(self):
        """
        """
        for atm in self.tls_group:
            yield atm

    def glal_calc_color(self, atom):
        """Overrides the GLAtomList coloring behavior and just
        colors using the tls_color.
        """
        return self.gldl_property_color_rgbf("tls_color")

    def glal_calc_color_U(self, atom):
        r, g, b = self.glal_calc_color(atom)
        dim     = 0.8
        return (r*dim, g*dim, b*dim)
    def glal_calc_color_Uellipse(self, atom):
        return self.glal_calc_color_U(atom)
    def glal_calc_color_Urms(self, atom):
        return self.glal_calc_color_U(atom)
    
    def glal_calc_color_trace(self):
        return self.gldl_property_color_rgbf("tls_color") 

    def glal_calc_U(self, atom):
        """Always return the reduced T tensor.
        """
        return self.properties["rT"]

    def draw_fan(self):
        """Draws a fan from the TLS group center of reaction to the
        TLS group backbone atoms.
        """
        COR     = self.properties["COR"]
        r, g, b = self.gldl_property_color_rgbf("tls_color")
        a       = self.properties["fan_opacity"]

        self.driver.glr_set_material_rgba(r, g, b, a)
        
        self.driver.glr_normalize_enable()
        self.driver.glr_light_two_sides_enable()
        self.driver.glr_begin_triangle_fan()

        ## driver optimization
        driver = self.driver
        ##

        v1 = None
        v2 = None

        for atm in self.tls_group:
            if atm.name not in ("N", "CA", "C"):
                continue

            if v1==None:
                v1 = atm.position - COR
                continue
            elif v2==None:
                v2 = atm.position - COR
                driver.glr_normal(cross(v1, v2))
                driver.glr_vertex3(0.0, 0.0, 0.0)
            else:
                v1 = v2
                v2 = atm.position - COR

            driver.glr_normal(cross(v1, v2))
            driver.glr_vertex(v1)
            driver.glr_vertex(v2)

        self.driver.glr_end()
        self.driver.glr_light_two_sides_disable()
        self.driver.glr_normalize_disable()
        

class GLTLSGroup(GLDrawList):
    """Top level visualization object for a TLS group.
    """
    def __init__(self, **args):
        self.tls_group = args["tls_group"]
        self.tls_info  = args["tls_info"]
        self.tls_name  = args["tls_name"]

        GLDrawList.__init__(self)
        self.glo_set_name(self.tls_name)

        ## add a child GLTLSAtomList for the animated atoms
        self.gl_atom_list = GLTLSAtomList(
            tls_group        = self.tls_group,
            trace            = True,
            lines            = False,
            fan_visible      = False)

        self.gl_atom_list.glo_set_name("TLS Atom Animation")
        self.gl_atom_list.glo_set_properties_id("gl_atom_list")
        self.glo_add_child(self.gl_atom_list)

        self.glo_link_child_property(
            "symmetry", "gl_atom_list", "symmetry")

        self.glo_link_child_property(
            "main_chain_visible", "gl_atom_list", "main_chain_visible")
        self.glo_link_child_property(
            "oatm_visible", "gl_atom_list", "oatm_visible")
        self.glo_link_child_property(
            "side_chain_visible", "gl_atom_list", "side_chain_visible") 
        self.glo_link_child_property(
            "hetatm_visible", "gl_atom_list", "hetatm_visible") 
        self.glo_link_child_property(
            "water_visible", "gl_atom_list", "water_visible")        
        self.glo_link_child_property(
            "hydrogen_visible", "gl_atom_list", "hydrogen_visible") 

        self.glo_link_child_property(
            "tls_color", "gl_atom_list", "tls_color")        

        self.glo_link_child_property(
            "fan_visible", "gl_atom_list", "fan_visible")
        self.glo_link_child_property(
            "fan_opacity", "gl_atom_list", "fan_opacity")
        self.glo_link_child_property(
            "axes_rT", "gl_atom_list", "U")
        self.glo_link_child_property(
            "ellipse_rT", "gl_atom_list", "ellipse")
        self.glo_link_child_property(
            "rms_rT", "gl_atom_list", "rms")

        self.glo_link_child_property(
            "adp_prob", "gl_atom_list", "adp_prob")

        self.glo_link_child_property(
            "COR", "gl_atom_list", "origin")
        self.glo_link_child_property(
            "COR", "gl_atom_list", "atom_origin")
        self.glo_link_child_property(
            "COR", "gl_atom_list", "COR")
        self.glo_link_child_property(
            "T", "gl_atom_list", "T")
        self.glo_link_child_property(
            "rT", "gl_atom_list", "rT")
        self.glo_link_child_property(
            "L", "gl_atom_list", "L")
        self.glo_link_child_property(
            "S", "gl_atom_list", "S")

        self.glo_link_child_property(
            "L1_eigen_vec", "gl_atom_list", "L1_eigen_vec")
        self.glo_link_child_property(
            "L2_eigen_vec", "gl_atom_list", "L2_eigen_vec")
        self.glo_link_child_property(
            "L3_eigen_vec", "gl_atom_list", "L3_eigen_vec")
        
        self.glo_link_child_property(
            "L1_eigen_val", "gl_atom_list", "L1_eigen_val")
        self.glo_link_child_property(
            "L2_eigen_val", "gl_atom_list", "L2_eigen_val")
        self.glo_link_child_property(
            "L3_eigen_val", "gl_atom_list", "L3_eigen_val")

        self.glo_link_child_property(
            "L1_rho", "gl_atom_list", "L1_rho")
        self.glo_link_child_property(
            "L2_rho", "gl_atom_list", "L2_rho")
        self.glo_link_child_property(
            "L3_rho", "gl_atom_list", "L3_rho")

        self.glo_link_child_property(
            "L1_pitch", "gl_atom_list", "L1_pitch")
        self.glo_link_child_property(
            "L2_pitch", "gl_atom_list", "L2_pitch")
        self.glo_link_child_property(
            "L3_pitch", "gl_atom_list", "L3_pitch")

        self.glo_link_child_property(
            "L1_rot", "gl_atom_list", "L1_rot")
        self.glo_link_child_property(
            "L2_rot", "gl_atom_list", "L2_rot")
        self.glo_link_child_property(
            "L3_rot", "gl_atom_list", "L3_rot")
 
        ## initalize properties
        self.glo_add_update_callback(self.tls_update_cb)

        if not self.tls_group.is_null():
            self.glo_init_properties(
                COR          = self.tls_info["COR"],

                T            = self.tls_info["T'"],
                rT           = self.tls_info["rT'"],
                L            = self.tls_info["L'"] * RAD2DEG2,
                S            = self.tls_info["S'"] * RAD2DEG,

                L1_eigen_vec = self.tls_info["L1_eigen_vec"],
                L2_eigen_vec = self.tls_info["L2_eigen_vec"],
                L3_eigen_vec = self.tls_info["L3_eigen_vec"],

                L1_eigen_val = self.tls_info["L1_eigen_val"] * RAD2DEG2,
                L2_eigen_val = self.tls_info["L2_eigen_val"] * RAD2DEG2,
                L3_eigen_val = self.tls_info["L3_eigen_val"] * RAD2DEG2,

                L1_rho       = self.tls_info["L1_rho"],
                L2_rho       = self.tls_info["L2_rho"],
                L3_rho       = self.tls_info["L3_rho"],

                L1_pitch     = self.tls_info["L1_pitch"] * (1.0/RAD2DEG),
                L2_pitch     = self.tls_info["L2_pitch"] * (1.0/RAD2DEG),
                L3_pitch     = self.tls_info["L3_pitch"] * (1.0/RAD2DEG),
                **args)
        else:
            self.glo_init_properties(**args)

    def set_tls_groupXXX(self, tls_group):
        """Set a new TLSGroup.
        """
        self.tls_group = tls_group

        if not self.tls_group.is_null():
            self.tls_info = self.tls_group.calc_tls_info()

            self.properties.update(
                COR          = self.tls_info["COR"],
                T            = self.tls_info["T'"],
                Tr           = self.tls_info["rT'"],
                L            = self.tls_info["L'"] * RAD2DEG2,
                S            = self.tls_info["S'"] * RAD2DEG,

                L1_eigen_vec = self.tls_info["L1_eigen_vec"],
                L2_eigen_vec = self.tls_info["L2_eigen_vec"],
                L3_eigen_vec = self.tls_info["L3_eigen_vec"],

                L1_eigen_val = self.tls_info["L1_eigen_val"] * RAD2DEG2,
                L2_eigen_val = self.tls_info["L2_eigen_val"] * RAD2DEG2,
                L3_eigen_val = self.tls_info["L3_eigen_val"] * RAD2DEG2,

                L1_rho       = self.tls_info["L1_rho"],
                L2_rho       = self.tls_info["L2_rho"],
                L3_rho       = self.tls_info["L3_rho"],

                L1_pitch     = self.tls_info["L1_pitch"] * (1.0/RAD2DEG),
                L2_pitch     = self.tls_info["L2_pitch"] * (1.0/RAD2DEG),
                L3_pitch     = self.tls_info["L3_pitch"] * (1.0/RAD2DEG) )

        else:
            self.tls_info = None

            self.properties.update(
                COR          = GLObject.PropertyDefault,
                T            = GLObject.PropertyDefault,
                Tr           = GLObject.PropertyDefault,
                L            = GLObject.PropertyDefault,
                S            = GLObject.PropertyDefault,

                L1_eigen_vec = GLObject.PropertyDefault,
                L2_eigen_vec = GLObject.PropertyDefault,
                L3_eigen_vec = GLObject.PropertyDefault,

                L1_eigen_val = GLObject.PropertyDefault,
                L2_eigen_val = GLObject.PropertyDefault,
                L3_eigen_val = GLObject.PropertyDefault,

                L1_rho       = GLObject.PropertyDefault,
                L2_rho       = GLObject.PropertyDefault,
                L3_rho       = GLObject.PropertyDefault,

                L1_pitch     = GLObject.PropertyDefault,
                L2_pitch     = GLObject.PropertyDefault,
                L3_pitch     = GLObject.PropertyDefault )

    def glo_install_properties(self):
        GLDrawList.glo_install_properties(self)

        ## TLS Analysis
        self.glo_add_property(
            { "name":        "COR",
              "desc":        "TLS Center of Reaction", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "T",
              "desc":        "T<sup>COR</sup> Tensor (A<sup>2</sup>)",
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3,3)",
              "default":     zeros((3,3), Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "rT",
              "desc":        "T<sup>r</sup> Tensor (A<sup>2</sup>)",
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3,3)",
              "default":     zeros((3,3), Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L",
              "desc":        "L<sup>COR</sup> Tensor (DEG<sup>2</sup>)",
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3,3)",
              "default":     zeros((3,3), Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "S",
              "desc":        "S<sup>COR</sup> Tensor (A*DEG)",
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3,3)",
              "default":     zeros((3,3), Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_eigen_vec",
              "desc":        "L<sub>1</sub> Eigen Vector", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_eigen_vec",
              "desc":        "L<sub>2</sub> Eigen Vector", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_eigen_vec",
              "desc":        "L<sub>3</sub> Eigen Vector", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_eigen_val",
              "desc":        "L<sub>1</sub> Eigen Value", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_eigen_val",
              "desc":        "L<sub>2</sub> Eigen Value", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_eigen_val",
              "desc":        "L<sub>3</sub> Eigen Value", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_rho",
              "desc":        "L<sub>1</sub> Position from COR", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_rho",
              "desc":        "L<sub>2</sub> Position from COR", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_rho",
              "desc":        "L<sub>3</sub> Position from COR", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_pitch",
              "desc":        "L<sub>1</sub> Screw Pitch (A/DEG)", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_pitch",
              "desc":        "L<sub>2</sub> Screw Pitch (A/DEG)", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_pitch",
              "desc":        "L<sub>3</sub> Screw Pitch (A/DEG)", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })

        ## Show/Hide
        self.glo_add_property(
            { "name":        "symmetry",
              "desc":        "Show Symmetry Equivelant",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":      "main_chain_visible",
              "desc":      "Show Main Chain Atoms",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   True,
              "action":    ["recompile", "recalc_positions"] })
        self.glo_add_property(
            { "name":      "oatm_visible",
              "desc":      "Show Main Chain Carbonyl Atoms",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   True,
              "action":    ["recompile", "recalc_positions"] })
        self.glo_add_property(
            { "name":      "side_chain_visible",
              "desc":      "Show Side Chain Atoms",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   True,
              "action":    ["recompile", "recalc_positions"] })
        self.glo_add_property(
            { "name":      "hetatm_visible",
              "desc":      "Show Hetrogen Atoms",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   True,
              "action":    ["recompile", "recalc_positions"] })
        self.glo_add_property(
            { "name":      "water_visible",
              "desc":      "Show Waters",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   False,
              "action":    ["recompile", "recalc_positions"] })
        self.glo_add_property(
            { "name":      "hydrogen_visible",
              "desc":      "Show Hydrogens",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   False,
              "action":    ["recompile", "recalc_positions"] })
        self.glo_add_property(
            { "name":        "fan_visible",
              "desc":        "Show COR-Backbone Fan",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "TLS_visible",
              "desc":        "Show TLS T<sup>r</sup> Ellipsoid/Screw Axes",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     True,
              "action":      "recompile_tensors" })
        self.glo_add_property(
            { "name":        "U",
              "desc":        "Show U<sup>TLS</sup> Thermal Axes",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile_Utls_axes" })
        self.glo_add_property(
            { "name":        "ellipse",
              "desc":        "Show U<sup>TLS</sup> Thermal Ellipsoids",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile_Utls_ellipse" })
        self.glo_add_property(
            { "name":        "rms",
              "desc":        "Show U<sup>TLS</sup> Thermal Peanuts",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile_Utls_rms" })

        self.glo_add_property(
            { "name":        "axes_rT",
              "desc":        "Show T<sup>r</sup> Thermal Axes", 
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "ellipse_rT",
              "desc":        "Show  T<sup>r</sup> Thermal Ellipsoids", 
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "rms_rT",
              "desc":        "Show  T<sup>r</sup> Thermal Peanuts",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile_Utls_rms" })
                
        self.glo_add_property(
            { "name":        "L1_visible",
              "desc":        "Show Screw L1 Displacement Surface", 
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile_surface" })
        self.glo_add_property(
            { "name":        "L2_visible",
              "desc":        "Show Screw L2 Displacement Surface", 
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile_surface" })
        self.glo_add_property(
            { "name":        "L3_visible",
              "desc":        "Show Screw L3 Displacement Surface",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile_surface" })
        
        ## TLS
        self.glo_add_property(
            { "name":        "add_biso",
              "desc":        "Add Atom B<sup>ISO</sup> to U<sup>TLS</sup>",
              "catagory":    "TLS",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "tls_color",
              "desc":        "TLS Group Visualization Color",
              "catagory":    "TLS",
              "type":        "enum_string",
              "default":     "Green",
              "enum_list":   self.gldl_color_list,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":       "adp_prob",
              "desc":       "Isoprobability Magnitude",
              "catagory":   "TLS",
              "type":       "integer",
              "range":      PROP_PROBABILTY_RANGE,
              "default":    50,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "L_axis_scale",
              "desc":       "Scale Screw Axis Length",
              "catagory":   "TLS",
              "type":       "float",
              "default":    5.00,
              "action":     "recompile_tensors" })        
        self.glo_add_property(
            { "name":       "L_axis_radius",
              "desc":       "Screw Axes Radius",
              "catagory":   "TLS",
              "type":       "float",
              "default":    0.4,
              "action":     "recompile_tensors" })
        self.glo_add_property(
            { "name":        "ellipse_opacity",
              "desc":        "U<sup>TLS</sup> Thermal Ellipsoid Opacity",
              "catagory":    "TLS",
              "type":        "float",
              "range":       PROP_OPACITY_RANGE,
              "default":     1.0,
              "action":      "recompile_Utls_ellipse" })
        self.glo_add_property(
            { "name":        "rms_opacity",
              "desc":        "U<sup>TLS</sup> Thermal Peanut Opacity",
              "catagory":    "TLS",
              "type":        "float",
              "range":       PROP_OPACITY_RANGE,
              "default":     1.0,
              "action":      "recompile_Utls_rms" })
        self.glo_add_property(
            { "name":        "surface_opacity",
              "desc":        "Screw Surface Opacity",
              "catagory":    "TLS",
              "type":        "float",
              "range":       PROP_OPACITY_RANGE,
              "default":     1.0,
              "action":      "recompile_surface" })
        self.glo_add_property(
            { "name":        "fan_opacity",
              "desc":        "COR-Backbone Fan Opacity",
              "catagory":    "TLS",
              "type":        "float",
              "range":       PROP_OPACITY_RANGE,
              "default":     1.0,
              "action":      "recompile_fan" })
        self.glo_add_property(
            { "name":        "time",
              "desc":        "Simulation Time",
              "catagory":    "TLS",
              "type":        "float",
              "range":       PROP_OPACITY_RANGE,
              "default":     0.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L1_rot",
              "desc":        "L<sub>1</sub> Simulated Rotation (DEG)",
              "catagory":    "TLS",
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L2_rot",
              "desc":        "L<sub>2</sub> Simulated Rotation (DEG)", 
              "catagory":    "TLS",
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L3_rot",
              "desc":        "L<sub>3</sub> Simulated Rotation (DEG)",
              "catagory":    "TLS",
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })

    def gldl_install_draw_methods(self):
        self.gldl_draw_method_install(
            { "name":                "tls_tensors",
              "func":                self.draw_tensors,
              "transparent":         False,
              "visible_property":    "TLS_visible",
              "recompile_action":    "recompile_tensors" })
        self.gldl_draw_method_install(
            { "name":                "Utls_axes",
              "func":                self.draw_Utls_axes,
              "transparent":         False,
              "visible_property":    "U",
              "recompile_action":    "recompile_Utls_axes" })
        self.gldl_draw_method_install(
            { "name":                "Utls_ellipse",
              "func":                self.draw_Utls_ellipse,
              "visible_property":    "ellipse",
              "opacity_property":    "ellipse_opacity",
              "recompile_action":    "recompile_Utls_ellipse" })
        self.gldl_draw_method_install(
            { "name":                "Utls_rms",
              "func":                self.draw_Utls_rms,
              "visible_property":    "rms",
              "opacity_property":    "rms_opacity",
              "recompile_action":    "recompile_Utls_rms" })
        self.gldl_draw_method_install(
            { "name":                "L1_surface",
              "func":                self.draw_L1_surface,
              "visible_property":    "L1_visible",
              "opacity_property":    "surface_opacity",
              "recompile_action":    "recompile_surface" })
        self.gldl_draw_method_install(
            { "name":                "L2_surface",
              "func":                self.draw_L2_surface,
              "visible_property":    "L2_visible",
              "opacity_property":    "surface_opacity",
              "recompile_action":    "recompile_surface" })
        self.gldl_draw_method_install(
            { "name":                "L3_surface",
              "func":                self.draw_L3_surface,
              "visible_property":    "L3_visible",
              "opacity_property":    "surface_opacity",
              "recompile_action":    "recompile_surface" })
         
    def tls_update_cb(self, updates, actions):
        if "time" in updates:
            self.update_time()

    def update_time(self):
        """Changes the time of the TLS group simulating harmonic motion.
        """
        if self.tls_group.is_null():
            return

        ## time should be in the range 0.0-0.1.0
        sin_tm = math.sin(2.0 * math.pi * self.properties["time"])

        ## calculate L eignvalue displacements at the given
        ## probability levels
        C = GAUSS3C[self.properties["adp_prob"]]

        try:
            L1_c = C * math.sqrt(self.properties["L1_eigen_val"])
        except ValueError:
            L1_c = 0.0
            L1_peak = 0.0

        try:
            L2_c = C * math.sqrt(self.properties["L2_eigen_val"])
        except ValueError:
            L2_c = 0.0
            L2_peak = 0.0

        try:
            L3_c = C * math.sqrt(self.properties["L3_eigen_val"])
        except ValueError:
            L3_c = 0.0
            L3_peak = 0.0

        L1_rot  = L1_c * sin_tm
        L2_rot  = L2_c * sin_tm
        L3_rot  = L3_c * sin_tm

        self.glo_update_properties(L1_rot=L1_rot, L2_rot=L2_rot, L3_rot=L3_rot)

    def gldl_iter_multidraw_self(self):
        """Specialized draw list invokation to recycle the draw list for
        symmetry related copies.  Cartesian versions of the symmetry rotation
        and translation operators are generated by GLStructure/UnitCell
        classes.
        """
        if self.properties["symmetry"]==False:
            yield True
            
        else:

            gl_struct = self.glo_get_glstructure()
            if gl_struct==None:
                yield True

            else:
                for symop in gl_struct.iter_orth_symops():
                    self.driver.glr_push_matrix()

                    glMultMatrixf(
                        (symop.R[0,0], symop.R[1,0], symop.R[2,0], 0.0,
                         symop.R[0,1], symop.R[1,1], symop.R[2,1], 0.0,
                         symop.R[0,2], symop.R[1,2], symop.R[2,2], 0.0,
                         symop.t[0],   symop.t[1],   symop.t[2],   1.0) )

                    yield True
                    self.driver.glr_pop_matrix()

    def gltls_iter_atoms(self):
        """Special atom iterator for the TLS drawing functions yields:
        atm, Utls
        """
        T = self.tls_group.T
        L = self.tls_group.L
        S = self.tls_group.S
        o = self.tls_group.origin
        
        for atm, visible in self.gl_atom_list.glal_iter_atoms_filtered():
            if not visible:
                continue

            Utls = self.tls_group.calc_Utls(T, L, S, atm.position - o)

            if self.properties["add_biso"]==True:
                Utls = Utls + (B2U * atm.temp_factor * identity(3, Float))
            
            yield atm, Utls
    
    def draw_tensors(self):
        """Draw tensor axis.
        """
        if self.tls_group.is_null():
            return

        self.driver.glr_push_matrix()

        ## get the TLS color
        r, g, b = self.gldl_property_color_rgbf("tls_color")
        self.driver.glr_set_material_rgb(r, g, b)
        self.driver.glr_translate(self.properties["COR"])
        ## T: units (A^2)
        self.driver.glr_Uellipse(
            (0.0,0.0,0.0), self.properties["rT"], self.properties["adp_prob"])

        ## L: units (DEG^2)
        L_scale = self.properties["L_axis_scale"]
        
        for Lx_eigen_val, Lx_eigen_vec, Lx_rho, Lx_pitch in [
            ("L1_eigen_val", "L1_eigen_vec", "L1_rho", "L1_pitch"),
            ("L2_eigen_val", "L2_eigen_vec", "L2_rho", "L2_pitch"),
            ("L3_eigen_val", "L3_eigen_vec", "L3_rho", "L3_pitch")]:

            L_eigen_vec = self.properties[Lx_eigen_vec]
            L_eigen_val = self.properties[Lx_eigen_val]
            L_rho       = self.properties[Lx_rho]
            L_pitch     = self.properties[Lx_pitch]
            
            if L_eigen_val<=0.0:
                continue

            C = GAUSS3C[self.properties["adp_prob"]]
            
            L_rot = C * (L_scale * math.sqrt(L_eigen_val))
            L_v   = L_eigen_vec * L_rot

            ## line from COR to center of screw/rotation axis
            ## draw lines from COR to the axis
            self.driver.glr_lighting_disable()
            self.driver.glr_line((0.0, 0.0, 0.0), L_rho)

            ## draw axis
            self.driver.glr_axis(
                L_rho - (0.5*L_v), L_v, self.properties["L_axis_radius"])

            ## draw disks with translational displacement
            L_screw_dis = L_eigen_vec * L_rot * L_pitch

            self.driver.glr_axis(
                L_rho - (0.5 * L_screw_dis),
                L_screw_dis,
                1.5 * self.properties["L_axis_radius"])

        self.driver.glr_pop_matrix()

    def draw_Utls_axes(self):
        """Render the anisotropic thremal axes calculated from the TLS
        model.
        """
        if self.tls_group.is_null():
            return
        
        prob = self.properties["adp_prob"]
        rgbf = self.gldl_property_color_rgbf("tls_color")

        glr_Uaxes = self.driver.glr_Uaxes

        for atm, Utls in self.gltls_iter_atoms():
            glr_Uaxes(atm.position, Utls, prob, rgbf, 1.0)

    def draw_Utls_ellipse(self):
        """Render the anisotropic thremal ellipsoids at the given probability
        contour calculated from the TLS model.
        """
        if self.tls_group.is_null():
            return
        
        prob    = self.properties["adp_prob"]
        r, g, b = self.gldl_property_color_rgbf("tls_color")
        a       = self.properties["ellipse_opacity"]
        self.driver.glr_set_material_rgba(r, g, b, a)

        glr_Uellipse = self.driver.glr_Uellipse

        for atm, Utls in self.gltls_iter_atoms():
            glr_Uellipse(atm.position, Utls, prob)

    def draw_Utls_rms(self):
        """Render the anisotropic thremal peanuts calculated from the TLS
        model.
        """
        if self.tls_group.is_null():
            return
        
        r, g, b = self.gldl_property_color_rgbf("tls_color")
        a       = self.properties["rms_opacity"]
        self.glr_set_material_rgb(r, g, b, a)

        for atm, Utls in self.gltls_iter_atoms():
            self.glr_Urms(atm.position, Utls)

    def draw_L1_surface(self):
        if self.tls_group.is_null():
            return
        self.draw_tls_surface(
            self.properties["L1_eigen_vec"],
            self.properties["L1_eigen_val"],
            self.properties["L1_rho"],
            self.properties["L1_pitch"])

    def draw_L2_surface(self):
        if self.tls_group.is_null():
            return
        self.draw_tls_surface(
            self.properties["L2_eigen_vec"],
            self.properties["L2_eigen_val"],
            self.properties["L2_rho"],
            self.properties["L2_pitch"])

    def draw_L3_surface(self):
        if self.tls_group.is_null():
            return
        self.draw_tls_surface(
            self.properties["L3_eigen_vec"],
            self.properties["L3_eigen_val"],
            self.properties["L3_rho"],
            self.properties["L3_pitch"])

    def draw_tls_surface(self, Lx_eigen_vec, Lx_eigen_val, Lx_rho, Lx_pitch):
        """Draws the TLS probability surface for a single non-intersecting
        screw axis.  Lx_eigen_val is the vaiance (mean square deviation MSD)
        of the rotation about the Lx_eigen_vec axis.
        """
        ## create a unique list of bonds which will be used to
        ## render the TLS surface; this list may be passed in a argument
        ## to avoid multiple calculations for each screw-rotation axis
        bond_list = []
        in_dict   = {}

        for atm, Utls in self.gltls_iter_atoms():
            in_dict[atm] = True

        for atm, Utls in self.gltls_iter_atoms():
            for bond in atm.iter_bonds():
                if in_dict.has_key(bond.get_partner(atm)):
                    bond_list.append(bond)
        
        ## this just won't work...
        if Lx_eigen_val==0.0:
            return

        C = GAUSS3C[self.properties["adp_prob"]]
        try:
            Lx_s = C * math.sqrt(Lx_eigen_val * DEG2RAD2)
        except ValueError:
            return

        Lx_pitch      = Lx_pitch * (1.0/DEG2RAD)
        COR           = self.properties["COR"]
        Lx_origin     = COR + Lx_rho
        steps         = 3
        rot_step      = Lx_s / float(steps)

        self.driver.glr_light_two_sides_enable()
        self.driver.glr_lighting_enable()
        self.driver.glr_normalize_enable()

        r, g, b = self.gldl_property_color_rgbf("tls_color")
        a       = self.properties["surface_opacity"]
        self.driver.glr_set_material_rgba(r, g, b, a)

        self.driver.glr_begin_quads()

        ## driver optimization
        glr_normal = self.driver.glr_normal
        glr_vertex = self.driver.glr_vertex
        ##

        for step in range(steps):
            step1 = rot_step * float(step)
            step2 = step1 + rot_step

            for sign in (-1.0, 1.0):
                rot1   = step1 * sign
                rot2   = step2 * sign

                Rstep1 = rmatrixu(Lx_eigen_vec, rot1)
                Rstep2 = rmatrixu(Lx_eigen_vec, rot2)

                screw1 = Lx_eigen_vec * (rot1 * Lx_pitch)
                screw2 = Lx_eigen_vec * (rot2 * Lx_pitch)

                for bond in bond_list:

                    pos1 = bond.atom1.position - Lx_origin
                    pos2 = bond.atom2.position - Lx_origin

                    v1 = matrixmultiply(Rstep1, pos1) + screw1
                    v2 = matrixmultiply(Rstep2, pos1) + screw2
                    v3 = matrixmultiply(Rstep2, pos2) + screw2
                    v4 = matrixmultiply(Rstep1, pos2) + screw1
                    
                    ## one normal perpendicular to the quad
                    glr_normal(cross(v2-v1, v4-v1))

                    glr_vertex(v1 + Lx_origin)
                    glr_vertex(v2 + Lx_origin)
                    glr_vertex(v3 + Lx_origin)
                    glr_vertex(v4 + Lx_origin)

        self.driver.glr_end()
        self.driver.glr_light_two_sides_disable()
        self.driver.glr_normalize_disable()


class GLTLSChain(GLDrawList):
    """Collects a list of GLTLSGroup instances which are all in the
    same chain.
    """
    def __init__(self, **args):
        GLDrawList.__init__(self)
        self.glo_set_name("TLS Chain %s" % (args["chain_id"]))
        self.glo_init_properties(**args)

    def glo_install_properties(self):
        GLDrawList.glo_install_properties(self)

        ## show/hide
        self.glo_add_property(
            { "name":        "symmetry",
              "desc":        "Show Symmetry Equivelant",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "" })
        self.glo_add_property(
            { "name":      "main_chain_visible",
              "desc":      "Show Main Chain Atoms",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   True,
              "action":    "" })
        self.glo_add_property(
            { "name":      "oatm_visible",
              "desc":      "Show Main Chain Carbonyl Atoms",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   True,
              "action":    ["recompile", "recalc_positions"] })
        self.glo_add_property(
            { "name":      "side_chain_visible",
              "desc":      "Show Side Chain Atoms",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   True,
              "action":    "" })
        self.glo_add_property(
            { "name":      "hetatm_visible",
              "desc":      "Show Hetrogen Atoms",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   True,
              "action":    "" })
        self.glo_add_property(
            { "name":      "water_visible",
              "desc":      "Show Waters",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   False,
              "action":    "" })
        self.glo_add_property(
            { "name":      "hydrogen_visible",
              "desc":      "Show Hydrogens",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   False,
              "action":    "" })
        self.glo_add_property(
            { "name":        "fan_visible",
              "desc":        "Show COR-Backbone Fan",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "TLS_visible",
              "desc":        "Show TLS T<sup>r</sup> Ellipsoid/Screw Axes",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     True,
              "action":      "" })
        self.glo_add_property(
            { "name":        "U",
              "desc":        "Show U<sup>TLS</sup> Thermal Axes",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "" })
        self.glo_add_property(
            { "name":        "ellipse",
              "desc":        "Show U<sup>TLS</sup> Thermal Ellipsoids",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "" })
        self.glo_add_property(
            { "name":        "rms",
              "desc":        "Show U<sup>TLS</sup> Thermal Peanuts",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "" })
        self.glo_add_property(
            { "name":        "axes_rT",
              "desc":        "Show T<sup>r</sup> Thermal Axes", 
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "" })
        self.glo_add_property(
            { "name":        "ellipse_rT",
              "desc":        "Show T<sup>r</sup> Thermal Ellipsoids", 
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "" })
        self.glo_add_property(
            { "name":        "rms_rT",
              "desc":        "Show T<sup>r</sup> Thermal Peanuts",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "" })
        self.glo_add_property(
            { "name":        "L1_visible",
              "desc":        "Show L<sub>1</sub> Screw Displacement Surface", 
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "" })
        self.glo_add_property(
            { "name":        "L2_visible",
              "desc":        "Show L<sub>2</sub> Screw Displacement Surface", 
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "" })
        self.glo_add_property(
            { "name":        "L3_visible",
              "desc":        "Show L<sub>3</sub> Screw Displacement Surface",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "" })
        
        ## TLS
        self.glo_add_property(
            { "name":        "add_biso",
              "desc":        "Add Atom B<sup>ISO</sup> to U<sup>TLS</sup>",
              "catagory":    "TLS",
              "type":        "boolean",
              "default":     False,
              "action":      "" })
        self.glo_add_property(
            { "name":       "adp_prob",
              "desc":       "Isoprobability Magnitude",
              "catagory":   "TLS",
              "type":       "integer",
              "range":      PROP_PROBABILTY_RANGE,
              "default":    50,
              "action":     "" })

        self.glo_add_property(
            { "name":       "L_axis_scale",
              "desc":       "Scale Screw Axis Length",
              "catagory":   "TLS",
              "type":       "float",
              "default":    5.00,
              "action":     "recompile_tensors" })
        self.glo_add_property(
            { "name":       "L_axis_radius",
              "desc":       "Screw Axes Radius",
              "catagory":   "TLS",
              "type":       "float",
              "default":    0.4,
              "action":     "" })
        self.glo_add_property(
            { "name":        "ellipse_opacity",
              "desc":        "U<sup>TLS</sup> Thermal Ellipsoid Opacity",
              "catagory":    "TLS",
              "type":        "float",
              "range":       PROP_OPACITY_RANGE,
              "default":     1.0,
              "action":      "" })
        self.glo_add_property(
            { "name":        "rms_opacity",
              "desc":        "U<sup>TLS</sup> Thermal Peanut Opacity",
              "catagory":    "TLS",
              "type":        "float",
              "range":       PROP_OPACITY_RANGE,
              "default":     1.0,
              "action":      "" })
        self.glo_add_property(
            { "name":        "surface_opacity",
              "desc":        "Screw Displacement Surface Opacity",
              "catagory":    "TLS",
              "type":        "float",
              "range":      PROP_OPACITY_RANGE,
              "default":     1.0,
              "action":      "" })
        self.glo_add_property(
            { "name":        "fan_opacity",
              "desc":        "COR-Backbone Fan Opacity",
              "catagory":    "TLS",
              "type":        "float",
              "range":       PROP_OPACITY_RANGE,
              "default":     1.0,
              "action":      "recompile_fan" })
##         ## color methods
##         self.glo_add_property(
##             { "name":        "color_method",
##               "desc":        "TLS Group Coloring Scheme",
##               "catagory":    "Color Methods",
##               "type":        "enum_string",
##               "default":     "Color by Group",
##               "enum_list":   ["Color By Group", "Color By Goodness of Fit"],
##               "action":      "" })


##     def color_by_group(self):
##         """Color TLS Groups by 
##         """
##         for gl_tls_group in self.glo_iter_children():
##             pass

    def add_gl_tls_group(self, gl_tls_group):
        self.glo_add_child(gl_tls_group)
        
        child_id = gl_tls_group.glo_get_properties_id()

        self.glo_link_child_property(
            "symmetry", child_id, "symmetry")        

        self.glo_link_child_property(
            "main_chain_visible", child_id, "main_chain_visible")
        self.glo_link_child_property(
            "oatm_visible", child_id, "oatm_visible")
        self.glo_link_child_property(
            "side_chain_visible", child_id, "side_chain_visible") 
        self.glo_link_child_property(
            "hetatm_visible", child_id, "hetatm_visible") 
        self.glo_link_child_property(
            "water_visible", child_id, "water_visible")        
        self.glo_link_child_property(
            "hydrogen_visible", child_id, "hydrogen_visible") 

        self.glo_link_child_property(
            "fan_visible", child_id, "fan_visible")
        self.glo_link_child_property(
            "fan_opacity", child_id, "fan_opacity")
        self.glo_link_child_property(
            "TLS_visible", child_id, "TLS_visible")

        self.glo_link_child_property(
            "U", child_id, "U")
        self.glo_link_child_property(
            "ellipse", child_id, "ellipse")
        self.glo_link_child_property(
            "rms", child_id, "rms")

        self.glo_link_child_property(
            "axes_rT", child_id, "axes_rT")
        self.glo_link_child_property(
            "ellipse_rT", child_id, "ellipse_rT")
        self.glo_link_child_property(
            "rms_rT", child_id, "rms_rT")

        self.glo_link_child_property(
            "L1_visible", child_id, "L1_visible")
        self.glo_link_child_property(
            "L2_visible", child_id, "L2_visible")
        self.glo_link_child_property(
            "L3_visible", child_id, "L3_visible")
        self.glo_link_child_property(
            "add_biso", child_id, "add_biso")
        self.glo_link_child_property(
            "adp_prob", child_id, "adp_prob")
        self.glo_link_child_property(
            "L_axis_scale", child_id, "L_axis_scale")
        self.glo_link_child_property(
            "L_axis_radius", child_id, "L_axis_radius")
        self.glo_link_child_property(
            "ellipse_opacity", child_id, "ellipse_opacity")
        self.glo_link_child_property(
            "rms_opacity", child_id, "rms_opacity")
        self.glo_link_child_property(
            "surface_opacity", child_id, "surface_opacity")



## <testing>
if __name__ == "__main__":
    print "==============================================="
    print "TEST CASE 1: TLS Class"
    print

    tls = TLS()
    tls.name = "All protein"
    tls.set_origin(18.885, 49.302, 13.315)
    tls.set_T(0.0263, 0.0561, 0.0048, -0.0128, 0.0065, -0.0157)
    tls.set_L(0.9730, 5.1496, 0.8488,  0.2151,-0.1296,  0.0815)
    tls.set_S(0.0007, 0.0281, 0.0336, -0.0446,-0.2288, -0.0551,
              0.0487, 0.0163)

    print tls

    print "eigenvalues(T)"
    print eigenvalues(tls.T)
    print "eigenvalues(L)"
    print eigenvalues(tls.L)

    print "==============================================="

    print 

    print "==============================================="
    print "TEST CASE 2: TLS/REFMAC file loading"

    import sys

    fil = open(sys.argv[1], "r")
    tls_list = LoadTLSGroups(fil)
    
    for tls in tls_list:
        print "----- TLS Group #%d ---" % (tls_list.index(tls))
        print tls

        print "eigenvalues(T)"
        print eigenvalues(tls.T)
        print "eigenvalues(L)"
        print eigenvalues(tls.L)

        print "-----------------------"

    print "==============================================="
## </testing>
