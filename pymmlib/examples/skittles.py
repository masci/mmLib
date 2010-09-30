#!/bin/env python
## Copyright 2002-2010 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
##
## skittles - a TLS validation tool for PDBs
##
## DESCRIPTION: This script is provided as an example on how to validate
##   structures using pymmlib.
## AUTHORS:
##   - P. Christoph Champ <champc _at_ u.washington.edu>
##   - Frank Zucker <frankz _at_ u.washington.edu>
##   - Ethan A. Merritt <merritt _at_ u.washington.edu>
## UPDATES:
##   - 2009-12-06: More documentation
##   - 2009-12-05: Ready for distribution
##   - 2009-10-30: Started script

from __future__ import generators
import sys
import os
import time
import re
import math

## Pymmlib
import pdbmodule
import tlsvld
#from mmLib import FileIO, TLS, Constants, AtomMath, Library
#from mmLib import AtomMath
from mmLib import Library, Constants

## Internal constants
AMINO_BACKBONE   = ["N", "CA", "C"]
NUCLEIC_BACKBONE = ["P", "O5'", "C5'", "C4'", "C3'", "O3'",
                    "P", "O5*", "C5*", "C4*", "C3*", "O3*"]
BACKBONE_ATOMS = AMINO_BACKBONE + NUCLEIC_BACKBONE
B2UE4 = 10000.0 / (8.0 * math.pi * math.pi)

tls_regex_dict = {
    "group":       re.compile("\s*TLS GROUP :\s+(\d+)\s*$"),
    "range":       re.compile("\s*RESIDUE RANGE :\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$"),
    "phenix_range":re.compile("\s*SELECTION:\s+CHAIN\s+(\w+)\s+AND\s+RESID\s+(-?\w+):(-?\w+)\s*$"),
    "origin":      re.compile("\s*ORIGIN\s+FOR\s+THE\s+GROUP\s+[(]A[)]:([\s\-\.0-9]+)$"),
    "t11_t22":     re.compile("\s*T11:\s*(\S+)\s+T22:\s*(\S+)\s*$"),
    "t33_t12":     re.compile("\s*T33:\s*(\S+)\s+T12:\s*(\S+)\s*$"),
    "t13_t23":     re.compile("\s*T13:\s*(\S+)\s+T23:\s*(\S+)\s*$"),
    "l11_l22":     re.compile("\s*L11:\s*(\S+)\s+L22:\s*(\S+)\s*$"),
    "l33_l12":     re.compile("\s*L33:\s*(\S+)\s+L12:\s*(\S+)\s*$"),
    "l13_l23":     re.compile("\s*L13:\s*(\S+)\s+L23:\s*(\S+)\s*$"),
    "s11_s12_s13": re.compile("\s*S11:\s*(\S+)\s+S12:\s*(\S+)\s+S13:\s*(\S+)\s*$"),
    "s21_s22_s23": re.compile("\s*S21:\s*(\S+)\s+S22:\s*(\S+)\s+S23:\s*(\S+)\s*$"),
    "s31_s32_s33": re.compile("\s*S31:\s*(\S+)\s+S32:\s*(\S+)\s+S33:\s*(\S+)\s*$")
    }
tlsout_regex_dict = {
    "group":  re.compile("(?:^\s*TLS\s*$)|(?:^\s*TLS\s+(.*)$)"),
    "range":  re.compile("^\s*RANGE\s+[']([A-Z])\s*([-0-9A-Z.]+)\s*[']\s+[']([A-Z])\s*([-0-9A-Z.]+)\s*[']\s*(\w*).*$"),
    "origin": re.compile("^\s*ORIGIN\s+(\S+)\s+(\S+)\s+(\S+).*$"),
    "T":      re.compile("^\s*T\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+).*$"),
    "L":      re.compile("^\s*L\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+).*$"),
    "S":      re.compile("^\s*S\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+).*$")
    }
#REMARK   3   OVERALL ANISOTROPIC B VALUE.
#REMARK   3    B11 (A**2) : 0.22000
#REMARK   3    B22 (A**2) : 0.03000
#REMARK   3    B33 (A**2) : -0.25000
#REMARK   3    B12 (A**2) : 0.00000
#REMARK   3    B13 (A**2) : 0.00000
#REMARK   3    B23 (A**2) : 0.00000

## SEE: http://www.wwpdb.org/documentation/format32/remark3.html
re_tls_range = re.compile("\s*RESIDUE RANGE :\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$")
re_tls_origin = re.compile("\s*ORIGIN\s+FOR\s+THE\s+GROUP\s+[(]A[)]:([\s\-\.0-9]+)$")

re_overall_b = re.compile("\s*OVERALL ANISOTROPIC B VALUE.*$")
re_b11 = re.compile("\s*B11 [(]A\*\*2[)] :\s+(\S+)*$")
re_b22 = re.compile("\s*B22 [(]A\*\*2[)] :\s+(\S+)*$")
re_b33 = re.compile("\s*B33 [(]A\*\*2[)] :\s+(\S+)*$")
re_b12 = re.compile("\s*B12 [(]A\*\*2[)] :\s+(\S+)*$")
re_b13 = re.compile("\s*B13 [(]A\*\*2[)] :\s+(\S+)*$")
re_b23 = re.compile("\s*B23 [(]A\*\*2[)] :\s+(\S+)*$")

def my_walk(path):
    """A better path walker
    """
    if os.path.isfile(path):
        yield path
    elif os.path.isdir(path):
        for x in os.listdir(path):
            for y in my_walk(os.path.join(path, x)):
                yield y

def convert_frag_id_save(frag_id):
    """Converts a mmLib fragment_id to a TLSOUT fragment id.
    """
    if len(frag_id) == 0:
        return "."
    if frag_id[-1].isdigit():
        frag_id += "."
    return frag_id

def set_origin(x, y, z):
    """Sets the TLS group origin of calculations.
    """
    self.origin = numpy.array((x, y, z), float)

def read_pdb(path):
    sec     = time.time()
    records = pdbmodule.read(path)
    sec     = time.time() - sec

    Ueq = []
    atom_obj = []
    atom_dict = {}
    listx = []
    state = 0
    overall = {}
    n = 1
    has_anisou = 0

    for rec in records:
        ## Start with empty dict
        atom_dict = {'u[1][1]': 0, 'chainID': '', 'name': '',
                    'resSeq': 0, 'u[0][1]': 0, 'u[0][0]': 0,
                    'u[0][2]': 0, 'u[1][2]': 0, 'element': '',
                    'RECORD': '', 'u[2][2]': 0, 'resName': '',
                    'serial': 0, 'has_anisou': 0,
                    'tempFactor': 0.0, 'x': 0.0, 'y': 0.0, 'z': 0.0,
                    'occupancy': 0.0}

        rec_type = rec["RECORD"]

        if rec_type == "REMARK":
            try:
                text = rec["text"]
            except KeyError:
                continue

            ## Capture the "OVERALL ANISOTROPIC B VALUE" fields
            if state == 0:
                m = re_overall_b.match(text)
                if m != None:
                    state = 1
            elif state == 1:
                m = re_b11.match(text)
                if m != None:
                    b11 = m.groups(1)
                    if b11[0] != "NULL":
                        overall[0] = float(b11[0])
                m = re_b22.match(text)
                if m != None:
                    b22 = m.groups(1)
                    if b22[0] != "NULL":
                        overall[1] = float(b22[0])
                m = re_b33.match(text)
                if m != None:
                    b33 = m.groups(1)
                    if b22[0] != "NULL":
                        overall[2] = float(b33[0])
                m = re_b12.match(text)
                if m != None:
                    b12 = m.groups(1)
                    if b22[0] != "NULL":
                        overall[3] = float(b12[0])
                m = re_b13.match(text)
                if m != None:
                    b13 = m.groups(1)
                    if b22[0] != "NULL":
                        overall[4] = float(b13[0])
                m = re_b23.match(text)
                if m != None:
                    b23 = m.groups(1)
                    if b22[0] != "NULL":
                        overall[5] = float(b23[0])
                    state = 2

            #if tls_desc.name != "":
            #    listx.append("TLS %s" % (tls_desc.name))
            #else:
            #    listx.append("TLS")

            ## Capture the TLS "RANGE"s fields
            m = re_tls_range.match(text)
            if m != None:
                (chain_id1, frag_id1, chain_id2, frag_id2) = m.groups()
                sel = "ALL"
                frag_id1 = convert_frag_id_save(frag_id1)
                frag_id2 = convert_frag_id_save(frag_id2)

                listx.append("#RANGE  '%s%s' '%s%s' %s" % (
                    chain_id1, frag_id1.rjust(5),
                    chain_id2, frag_id2.rjust(5), sel))

            ## Capture the TLS "ORIGIN"s fields
            m = re_tls_origin.match(text)
            if m != None:
                strx = m.group(1)
                ## this is nasty -- I wish I could trust the numbers
                ## to stay in fixed columns, but I can't
                ox = [0.0, 0.0, 0.0]
                for i in (0,1,2):
                    j = strx.find(".")
                    if j==-1:
                        break
                    x = strx[ max(0, j-4) : j+5]
                    strx = strx[j+5:]
                    ox[i] = float(x)

                try:
                    listx.append("#ORIGIN   %8.4f %8.4f %8.4f" % (
                        ox[0], ox[1], ox[2]))
                except:
                    print "ERROR!"
                    pass

            ## Capture the actual TLS values/fields.
            ## NOTE: These can be provided via an external TLSIN/TLSOUT file.
            #m = re_tls_T.match(text)
            #if m != None:
            #if tls_desc.T is not None:
            #    ## REFMAC ORDER: t11 t22 t33 t12 t13 t23
            #    listx.append("T   %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f" % (
            #        tls_desc.T[0,0], tls_desc.T[1,1], tls_desc.T[2,2],
            #        tls_desc.T[0,1], tls_desc.T[0,2], tls_desc.T[1,2]))
            #
            #m = re_tls_L.match(text)
            #if m != None:
            #if tls_desc.L is not None:
            #     ## REFMAC ORDER: l11 l22 l33 l12 l13 l23
            #     listx.append("L   %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f" % (
            #        tls_desc.L[0,0] * Constants.RAD2DEG2,
            #        tls_desc.L[1,1] * Constants.RAD2DEG2,
            #        tls_desc.L[2,2] * Constants.RAD2DEG2,
            #        tls_desc.L[0,1] * Constants.RAD2DEG2,
            #        tls_desc.L[0,2] * Constants.RAD2DEG2,
            #        tls_desc.L[1,2] * Constants.RAD2DEG2))
            #
            #m = re_tls_S.match(text)
            #if m != None:
            #if tls_desc.S is not None:
            #    ## REFMAC ORDER:
            #    ## <S22 - S11> <S11 - S33> <S12> <S13> <S23> <S21> <S31> <S32>
            #    listx.append(
            #        "S   %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f" % (
            #        (tls_desc.S[1,1] - tls_desc.S[0,0]) * Constants.RAD2DEG,
            #        (tls_desc.S[0,0] - tls_desc.S[2,2]) * Constants.RAD2DEG,
            #        tls_desc.S[0,1] * Constants.RAD2DEG,
            #        tls_desc.S[0,2] * Constants.RAD2DEG,
            #        tls_desc.S[1,2] * Constants.RAD2DEG,
            #        tls_desc.S[1,0] * Constants.RAD2DEG,
            #        tls_desc.S[2,0] * Constants.RAD2DEG,
            #        tls_desc.S[2,1] * Constants.RAD2DEG))

        ## Now, capture the ATOM and/or HETATM fields
        elif rec_type == "ATOM  " or rec_type == "HETATM":
            if Library.library_is_standard_residue(rec["resName"]) and \
               rec["name"].strip() in Constants.BACKBONE_ATOMS:
                atom_dict["serial"] = rec["serial"]
                atom_dict["name"] = rec["name"]
                atom_dict["resName"] = rec["resName"]
                atom_dict["chainID"] = rec["chainID"]
                atom_dict["resSeq"] = rec["resSeq"]
                atom_dict["x"] = rec["x"]
                atom_dict["y"] = rec["y"]
                atom_dict["z"] = rec["z"]
                atom_dict["occupancy"] = rec["occupancy"]
                atom_dict["tempFactor"] = rec["tempFactor"]
                #atom_dict["u[1][1]"] = 0.0
                # Initialize ANISOU values to isotropic = I3x3 * B/8*pi^2 
		#     in case this atom has no ANISOU record, e.g. 
		#     if it is next to but not in a TLS group
		Uiso = rec["tempFactor"] * B2UE4
                atom_dict["u[0][0]"] = Uiso
                atom_dict["u[1][1]"] = Uiso
                atom_dict["u[2][2]"] = Uiso
                atom_dict["u[0][1]"] = 0.0
                atom_dict["u[0][2]"] = 0.0
                atom_dict["u[1][2]"] = 0.0
                atom_dict["has_anisou"] = 0
                atom_obj.append(atom_dict)

        ## Capture the ANISOU fields, if they are present
        elif rec_type == "ANISOU":
            #print rec["serial"], rec["resName"]
            #print rec["u[0][0]"], rec["u[1][1]"], rec["u[2][2]"]

            try:
                Usum = rec["u[0][0]"] + rec["u[1][1]"] + rec["u[2][2]"]
            except KeyError:
                print "# STRANGE Uij VALUE(S) FOR ATOM=%s; RESIDUE=%s:%s:%s" % (
                    rec["serial"], rec["chainID"], rec["resName"], rec["resSeq"])
                continue

            try:
                #text = rec["text"]
                foo = 1
            except KeyError:
                continue

            if Library.library_is_standard_residue(rec["resName"]) and \
               rec["name"].strip() in Constants.BACKBONE_ATOMS:
                ## We only keep the backbone atom's Uij values.
                U11 = rec["u[0][0]"]
                U22 = rec["u[1][1]"]
                U33 = rec["u[2][2]"]
                U12 = U21 = rec["u[0][1]"]
                U13 = U31 = rec["u[0][2]"]
                U23 = U32 = rec["u[1][2]"]

                ADD_B_OVERALL = False # (default: False)
                if ADD_B_OVERALL:
                    atom_obj[-1]["u[0][0]"] = rec["u[0][0]"] + overall[0]
                    atom_obj[-1]["u[1][1]"] = rec["u[1][1]"] + overall[1]
                    atom_obj[-1]["u[2][2]"] = rec["u[2][2]"] + overall[2]
                    atom_obj[-1]["u[0][1]"] = rec["u[0][1]"] + overall[3]
                    atom_obj[-1]["u[0][2]"] = rec["u[0][2]"] + overall[4]
                    atom_obj[-1]["u[1][2]"] = rec["u[1][2]"] + overall[5]
                else:
                    atom_obj[-1]["u[0][0]"] = rec["u[0][0]"]
                    atom_obj[-1]["u[1][1]"] = rec["u[1][1]"]
                    atom_obj[-1]["u[2][2]"] = rec["u[2][2]"]
                    atom_obj[-1]["u[0][1]"] = rec["u[0][1]"]
                    atom_obj[-1]["u[0][2]"] = rec["u[0][2]"]
                    atom_obj[-1]["u[1][2]"] = rec["u[1][2]"]

                    atom_obj[-1]["has_anisou"] = 1
                has_anisou = 1

    return "\n".join(listx), overall, atom_obj

def get_chain_type(res_name):
    """Assign an integer value for residue type, identified by the actual name
    of the residue. Assigns "0" to amino acids, "1" to nucleic acids, and "9"
    to everything else.
    """
    if Library.library_is_amino_acid(res_name):
        return 0
    elif Library.library_is_nucleic_acid(res_name):
        return 1
    else:
        ## Unknown residue type
        print "# UNKOWN RESIDUE!"
        return 9

def get_res_frac(chain_type, atm_name, res_name):
    """Assigns a fractional value to the atom type within a given residue. The
    mod of an index of atom types should return fractional values of "0.1",
    "0.3", and "0.7" for most amino acid residues.
    """
    chain_type = get_chain_type(res_name)

    if chain_type == 0:
        frac = "%d" % (10*(int(Constants.AMINO_BACKBONE.index(atm_name)) % 3)/3.0)
    elif chain_type == 1:
        frac = "%d" % (10*(int(Constants.NUCLEIC_BACKBONE.index(atm_name)) % 6)/6.0)
    else:
        frac = 99

    return frac

if __name__ == "__main__":
    try:
        path = sys.argv[1]
    except IndexError:
        sys.exit(1)

    i = 0
    for pathx in my_walk(path):
        i += 1
        listx, overall, atom_obj = read_pdb(pathx)
        print "# %s" % listx

        chain_type = get_chain_type(atom_obj[0]["resName"])

        ## This is the header of the output
        print "#bond_num, bond, tls_group, cc_UV, Suij, Rosenfeld, diff_trace_UV, sum_square_diff"

        for u in range(1, len(atom_obj)):
            ## Example record/dictionary
            ## rec = {'u[1][1]': 5363, 'chainID': 'A', 'name': ' N',
            ##        'resSeq': 23, 'u[0][1]': -1790, 'u[0][0]': 5637,
            ##        'u[0][2]': -220, 'u[1][2]': -817, 'element': ' N',
            ##        'RECORD': 'ANISOU', 'u[2][2]': 5399, 'resName': 'GLU',
            ##        'serial': 1}
            cid1      = atom_obj[u-1]["chainID"]
            cid2      = atom_obj[u]["chainID"]
            res_name1 = atom_obj[u-1]["resName"].strip()
            res_name2 = atom_obj[u]["resName"].strip()
            res_num1  = atom_obj[u-1]["resSeq"]
            res_num2  = atom_obj[u]["resSeq"]
            atm_name1 = atom_obj[u-1]["name"].strip()
            atm_name2 = atom_obj[u]["name"].strip()
            x1        = atom_obj[u-1]["x"]
            y1        = atom_obj[u-1]["y"]
            z1        = atom_obj[u-1]["z"]
            x2        = atom_obj[u]["x"]
            y2        = atom_obj[u]["y"]
            z2        = atom_obj[u]["z"]
            Upos      = tlsvld.array([x1,y1,z1])
            Vpos      = tlsvld.array([x2,y2,z2])
            tls_group = 0

            ## Calculate the distance between two atoms
            d = math.sqrt( (x1-x2)**2 + (y1-y2)**2 + (x1-x2)**2 )
            #d = AtomMath.calc_distance(atom_obj[u-1], atom_obj[u])
            if ( (d > 2.0) and (cid1 == cid2) ):
                print "# BOND BREAK!"
                continue

            ## Add newlines for new chain_id
            if cid1 != cid2:
                chain_type = get_chain_type(res_name2)
                print "\n"
                continue

            ## Assign a "bond-number" to the two given atoms
            bond = "%s-%s" % (atm_name1, atm_name2)
            try:
                bond_pos = "%s.%s" % (
                    res_num1, get_res_frac(chain_type, atm_name1, res_name2))
            except:
                print "# POSSIBLE LIGAND FOUND!"
                continue

            ## Ready the data for Fortran77-style arrays
            U11 = atom_obj[u-1]["u[0][0]"]
            U22 = atom_obj[u-1]["u[1][1]"]
            U33 = atom_obj[u-1]["u[2][2]"]
            U12 = U21 = atom_obj[u-1]["u[0][1]"]
            U13 = U31 = atom_obj[u-1]["u[0][2]"]
            U23 = U32 = atom_obj[u-1]["u[1][2]"]
            U = tlsvld.array([[U11,U12,U13],[U21,U22,U23],[U31,U32,U33]])
            V11 = atom_obj[u]["u[0][0]"]
            V22 = atom_obj[u]["u[1][1]"]
            V33 = atom_obj[u]["u[2][2]"]
            V12 = V21 = atom_obj[u]["u[0][1]"]
            V13 = V31 = atom_obj[u]["u[0][2]"]
            V23 = V32 = atom_obj[u]["u[1][2]"]
            V = tlsvld.array([[V11,V12,V13],[V21,V22,V23],[V31,V32,V33]])

            ## Calculate all of the residuals
            Suij = rosenfeld = diff_trace_UV = sum_square_diff = 0.0
            cc_UV = tlsvld.calc_ccuij(U, V)
            #cc_UV = AtomMath.calc_CCuij([[U11,U12,U13],
            #                             [U21,U22,U23],
            #                             [U31,U32,U33]],
            #                            [[V11,V12,V13],
            #                             [V21,V22,V23],
            #                             [V31,V32,V33]])
            rosenfeld = tlsvld.rosenfeld(Upos, Vpos, U, V)
            sum_square_diff = tlsvld.sum_square_diff(U, V)
            diff_trace_UV = tlsvld.diff_trace_uv(U, V)

            ## The actual output of this entire script.
            print "%s\t%s:%s-%s:%s\t%s\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f" % (
                bond_pos, cid1, res_name1, res_num1, bond, tls_group,
                cc_UV, Suij, rosenfeld/10000,
                diff_trace_UV/10000, sum_square_diff/10000**2)
