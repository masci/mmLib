# speed test for the C PDB reading module -- it looks like it's about
# 450% faster than the Python implementation
from __future__ import generators

import sys
import os
import time
import pdbmodule

def my_walk(path):
    """All Python's path walk functions suck.
    """
    if os.path.isfile(path):
        yield path
    elif os.path.isdir(path):
        for x in os.listdir(path):
            for y in my_walk(os.path.join(path, x)):
                yield y


def read_pdb(path):
    sec     = time.time()
    records = pdbmodule.read(path)
    sec     = time.time() - sec


    ## collect some statistics from the PDB records
    stats = {}
    stats["time"] = sec

    for rec in records:
        rec_type = rec["RECORD"]
        
        if rec_type == "REMARK":
            try:
                text = rec["text"]
            except KeyError:
                continue

            if text.find("RESOLUTION RANGE HIGH") == 1:
                try:
                    stats["res"] = float(text[33:])
                except ValueError:
                    pass

        elif rec_type == "ATOM  " or rec_type == "HETATM":
            try:
                stats["atoms"] += 1
            except KeyError:
                stats["atoms"] = 1

        elif rec_type == "ANISOU":
            try:
                stats["anisou"] += 1
            except KeyError:
                stats["anisou"] = 1

    return stats


if __name__ == "__main__":
    try:
        path = sys.argv[1]
    except IndexError:
        sys.exit(1)

    i = 0
    for pathx in my_walk(path):
        i += 1
        stats = read_pdb(pathx)
        print "%d:%s:%s" % (i, pathx, stats)
