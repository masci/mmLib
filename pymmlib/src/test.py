#!/usr/bin/env python
## Copyright 2002-2010 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Speed test for the C PDB reading module -- it appears to be about 450%
faster than the Python implementation.
"""

from __future__ import generators

import sys
import os
import time
import types
import re
import shelve
import pdbmodule
from Numeric import *
from LinearAlgebra import *

aniso_bin_width = 0.025
aniso_bins      = 40

res_dict = {
    "GLY": True,
    "ALA": True,
    "VAL": True,
    "LEU": True,
    "ILE": True,
    "PRO": True,
    "PHE": True,
    "TYR": True,
    "TRP": True,
    "MET": True,
    "CYS": True,
    "SER": True,
    "THR": True,
    "ASP": True,
    "GLU": True,
    "HIS": True,
    "LYS": True,
    "ARG": True,
    "ASN": True,
    "GLN": True
    }

re_program = re.compile("\s*PROGRAM[^:]*:\s*(.*)$")
re_res     = re.compile("\s*RESOLUTION\s+RANGE\s+HIGH[^:]+:\s+(\S+)\s*$")

def my_walk(path):
    """All Python's path walk functions suck.
    """
    if os.path.isfile(path):
        yield path
    elif os.path.isdir(path):
        for x in os.listdir(path):
            for y in my_walk(os.path.join(path, x)):
                yield y

def mk_aniso_dict():
    aniso_dict = {}
    for bin in range(1,aniso_bins+1):
        aniso_dict[bin] = 0
    return aniso_dict

def inc_aniso_dict(aniso_dict, aniso):
    assert aniso > 0.0

    x = aniso / aniso_bin_width
    r = x % 1.0
    if r > 0.0:
        bin = int(x) + 1
    else:
        bin = int(x)
    aniso_dict[bin] += 1

def prnt_aniso_dict(aniso_dict):
    for bin in range(1, aniso_bins+1):
        print "%d: %d" % (bin, aniso_dict[bin])

def set_aniso_dict(rec, aniso_dict):
    assert rec["RECORD"] == "ANISOU"

    U = array(
        [ [rec["u[0][0]"], rec["u[0][1]"], rec["u[0][2]"]],
          [rec["u[0][1]"], rec["u[1][1]"], rec["u[1][2]"]],
          [rec["u[0][2]"], rec["u[1][2]"], rec["u[2][2]"]] ],
        Float) / 10000.0

    evals = eigenvalues(U)

    try:
        aniso = min(evals) / max(evals)
    except ZeroDivisionError:
        pass
    else:
        if aniso > 0.0:
            inc_aniso_dict(aniso_dict, aniso)
            return aniso

    return None

def read_pdb(path):
    sec     = time.time()
    records = pdbmodule.read(path)
    sec     = time.time() - sec

    ## collect some statistics from the PDB records
    stats = {}
    stats["path"] = path

    ## adverages
    Bnum_prot = 0
    Badv_prot = 0.0
    Unum_prot = 0
    Uadv_prot = 0.0
    Anum_prot = 0
    Aadv_prot = 0.0

    for rec in records:
        rec_type = rec["RECORD"]

        if rec_type == "ATOM  " or rec_type == "HETATM":

            if res_dict.has_key(rec.get("resName", "")):
                try:
                    stats["atoms"] += 1
                except KeyError:
                    stats["atoms"] = 1

                try:
                    Badv_prot += rec["tempFactor"]
                except KeyError:
                    pass
                else:
                    Bnum_prot += 1

                cid = rec.get("chainID")
                if cid:
                    try:
                        pchains = stats["prot_chains"]
                    except KeyError:
                        stats["prot_chains"] = [cid]
                    else:
                        if cid not in pchains:
                            pchains.append(cid)

            elif rec_type == "HETATM":
                try:
                    stats["hetatoms"] += 1
                except KeyError:
                    stats["hetatoms"] = 1            

        elif rec_type == "ANISOU":
            try:
                stats["anisou"] += 1
            except KeyError:
                stats["anisou"] = 1

            if res_dict.has_key(rec.get("resName", "")):
                Usum = rec["u[0][0]"]+rec["u[1][1]"]+rec["u[2][2]"]
                Uadv_prot += float(Usum)/30000.0
                Unum_prot += 1

                try:
                    aniso = set_aniso_dict(rec, stats["aniso_dict"])
                except KeyError:
                    stats["aniso_dict"] = mk_aniso_dict()
                    aniso = set_aniso_dict(rec, stats["aniso_dict"])

                if aniso != None:
                    Anum_prot += 1
                    Aadv_prot += aniso

        elif rec_type == "REMARK":
            try:
                text = rec["text"]
            except KeyError:
                continue

            m = re_res.match(text)
            if m != None:
                try:
                    stats["res"] = float(m.group(1))
                except ValueError:
                    pass
                continue

            m = re_program.match(text)
            if m != None:
                stats["program"] = m.group(1)

        elif rec_type == "REVDAT":
            if rec.get("modType") == 0:
                stats["date"] = rec.get("modDate", "dd-mmm-yy")

        elif rec_type == "EXPDTA":
            stats["tech"] = rec["technique"]

        elif rec_type == "HEADER":
            stats["id"] = rec.get("idCode", "XXXX")

    ## finish calculating some adverages
    try:
        stats["Badv_prot"] = Badv_prot / Bnum_prot
    except ZeroDivisionError:
        pass

    try:
        stats["Uadv_prot"] = Uadv_prot / Unum_prot
    except ZeroDivisionError:
        pass

    try:
        stats["Aadv_prot"] = Aadv_prot / Anum_prot
    except ZeroDivisionError:
        pass

    return stats

def prnt_stats(stats):
    for (key, val) in stats.items():

        if type(val) == types.FloatType:
            val = "%.2f" % (val)
        else:
            val = str(val)

        print "%s: %s" % (key.rjust(15), val)

if __name__ == "__main__":
    try:
        path = sys.argv[1]
    except IndexError:
        sys.exit(1)

    i = 0
    for pathx in my_walk(path):
        i += 1
        stats = read_pdb(pathx)

        print
        prnt_stats(stats)

        if stats.has_key("Aadv_prot") and stats.has_key("id"):
            db = shelve.open("anisodb.dat")
            db[stats["id"]] = stats
            db.close()
