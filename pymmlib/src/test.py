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

    if records != None:
        print "%s: %d records in %.2f seconds" % (
            path, len(records), sec)
    else:
        print "%s: NO RECORDS" % (path)

    for rec in records:
        if rec["RECORD"] == "REMARK":
            try:
                text = rec["text"]
            except KeyError:
                pass
            else:
                if text.find("RESOLUTION RANGE HIGH") == 1:
                    print text


if __name__ == "__main__":

    try:
        path = sys.argv[1]
    except IndexError:
        sys.exit(1)

    i = 0

    for pathx in my_walk(path):
        i += 1
        print str(i)+": ",
        read_pdb(pathx)
