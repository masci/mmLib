from __future__ import generators
import os
import sys
import re
from mmLib.mmTypes import *

"""Misc. testing utility code code common to the test programs.
"""

def walk(path, start_path, *regex_args):
    """Iterate over all files rooted at path containing the substring
    in the filename, including extentions.
    """
    if os.path.isfile(path):
        yield path
        return

    re_list = [re.compile(x) for x in regex_args]

    for (dirpath, dirnames, filenames) in os.walk(path):
        for filename in filenames:
            do_yield = False

            for rex in re_list:
                match = rex.match(filename)
                if match:
                    do_yield = True
                    break

            yield_path = os.path.join(dirpath, filename)

            if start_path != None:
                if start_path == yield_path:
                    start_path = None
                else:
                    do_yield = False
                
            if do_yield:
                yield yield_path

def walk_pdb(path, start_path = None):
    return walk(path, start_path, "pdb\w+\.gz", "w+\.pdb", "w+\.pdb\.gz")

def walk_cif(path, start_path = None):
    return walk(path, start_path, "\w+\.cif", "\w+\.cif\.gz")

def walk_pdb_cif(path, start_path = None):
    return walk(path, start_path,
                "pdb\w+\.gz", "w+\.pdb", "w+\.pdb\.gz",
                "\w+\.cif", "\w+\.cif\.gz")

def pdb_stats(path):
    re_model = re.compile("^MODEL\s+(\d+).*")
    re_atom = re.compile("^(?:ATOM|HETATM)\s*(\d+).*")

    model = 1
    serial_map = {}
    stats = {"atoms" : 0}

    for ln in OpenFile(path, "r").readlines():

        ## change model
        m = re_model.match(ln)
        if m != None:
            model = m.group(1)
            continue

        ## count atoms
        m = re_atom.match(ln)
        if m != None:
            stats["atoms"] += 1

            ser = m.group(1)
            ser = "%s-%s" % (ser, model)

            if serial_map.has_key(ser):
                print "DUPLICATE ID"
                print "[1]",serial_map[ser]
                print "[2]",ln
                sys.exit(1)
            else:
                serial_map[ser] = ln

    return stats

def cif_stats(path):
    re_atom = re.compile("^(?:ATOM|HETATM)\s+(\d+)\s+.*")

    atom_site_ids = {}    
    start_counting_atoms = False
    stats = {"atoms" : 0}

    for ln in OpenFile(path, "r").readlines():

        if ln.startswith("_atom_site."):
            start_counting_atoms = True

        if not start_counting_atoms:
            continue

        ## count atoms
        m = re_atom.match(ln)
        if m != None:
            stats["atoms"] += 1
            aid = m.group(1)

            if atom_site_ids.has_key(aid):
                print "DUPLICATE ID"
                print "[1]",atom_site_ids[aid]
                print "[2]",ln
                sys.exit(1)
            else:
                atom_site_ids[aid] = ln
        
    return stats


if __name__ == "__main__":
    for x in walk("/mnt/cdrom", None, "pdb\w+\.gz"):
        print x
