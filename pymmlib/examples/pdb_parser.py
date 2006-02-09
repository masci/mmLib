#!/usr/bin/env python
## This program exersizes the PDB parser by walking through a directory
## of PDB files and processing each one.  The PDBProcessor class is a
## very simple custom PDB file processor

import os, sys
from mmLib import PDB


class PDBProcessor(object):
    """Implement callbacks for PDB record types.  If you want the callback
    with the raw mmLib.PDB classes, prefix the method name with 'process_',
    if you want callback argument to be the result of the mmLib.PDB record
    class's 'process' method, then use the prefix 'preprocess_'.
    Implement only the callback you want to handle.
    """
    def process_HEADER(self, x):
        print "HEADER"
        print x
        print
    
    def preprocess_COMPND(self, x):
        print "COMPND"
        print x
        print

    def preprocess_OBSLTE(self, x):
        print "OBSLTE"
        print x
        print
        
    def preprocess_REVDAT(self, x):
        print "REVDAT"
        print x
        print

    def preprocess_SPRSDE(self, x):
        print "SPRSDE"
        print x
        print

    def preprocess_default(self, rec_name, x):
        """This method will be called, if it exists, for any PDB handler
        without its own method handler.
        """
        print "-----",rec_name
        print x
        print


def main(path):
    pdb_file = PDB.PDBFile()
    pdb_file.load_file(path)

    proc = PDBProcessor()
    pdb_file.record_processor(proc)


if __name__ == "__main__":
    import os

    try:
        path = sys.argv[1]
    except IndexError:
        print "usage: pdb_test.py <PDB file or directory of PDB files>"
        sys.exit(1)

    if os.path.isfile(path):
        main(path)
    elif os.path.isdir(path):
        for name in os.listdir(path):
            name = os.path.join(path, name)
            if not os.path.isfile(name):
                continue
            try:
                main(name)
            except:
                print "ERROR: ",name
                raise
