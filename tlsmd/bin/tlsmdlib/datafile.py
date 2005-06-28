## TLS Minimized Domains (TLSMD)
## Copyright 2002 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import cPickle
import bsddb
from misc import *

SAVEKEYS = [
    "error", "lsq_residual", "method",
    "cor_x", "cor_y", "cor_z",
    "t11", "t22", "t33", "t12", "t13", "t23",
    "l11", "l22", "l33", "l12", "l13", "l23",
    "s2211", "s1133", "s12", "s13", "s23", "s21", "s31", "s32" ]

class TLSMDFile(object):
    """Manages several files which save the state information of a TLSMD
    computation.  Current files handled are:

    GRH: Contains the TLS descriptions of all of the TLS fit subsegments of
         a structure.
    """
    def __init__(self, tlsdb_file):
        self.tlsdb_file = tlsdb_file
        ## graph-file uses Berkeley DB
        self.grh_db = bsddb.hashopen(self.tlsdb_file, "c")

    def grh_db_key(self, chain_id, frag_id1, frag_id2):
        db_key = "%s:%s:%s" % (chain_id, frag_id1, frag_id2)
        return db_key

    def grh_get_tls_record(self, chain_id, frag_id1, frag_id2):
        db_key = self.grh_db_key(chain_id, frag_id1, frag_id2)

        if not self.grh_db.has_key(db_key):
            return None

        rec = self.grh_db[db_key]
        tls_rec = cPickle.loads(rec)

        return tls_rec

    def grh_append_tls_record(self, tls):
        """Stores a tls (Python dictionary) record to the GRH file and to
        the in-memory cache.
        """
        chain_id = tls["chain_id"]
        frag_id1 = tls["frag_id1"]
        frag_id2 = tls["frag_id2"]

        db_key = self.grh_db_key(chain_id, frag_id1, frag_id2)
        
        ## create a uncluttered Python dictionary for pickeling/storage
        tls_rec = {}
        tls_rec["chain_id"] = chain_id
        tls_rec["frag_id1"] = frag_id1
        tls_rec["frag_id2"] = frag_id2

        for key in SAVEKEYS:
            if tls.has_key(key):
                tls_rec[key] = tls[key]
        
        self.grh_db[db_key] = cPickle.dumps(tls_rec)
        self.grh_db.sync()


### <testing>
if __name__=="__main__":
    import sys
    prefix = sys.argv[1]

    tlsmdfile = TLSMDFile(prefix)

    for tls in tlsmdfile.grh_iter_tls_records():
        print "%s %s %s %s" % (tls["method"], tls["chain_id"],
                               tls["frag_id1"], tls["frag_id2"])

### </testing>
