#!/usr/bin/env python
## TLS Minimized Domains (TLSMD)
## Copyright 200-20052 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import os
import sys
import xmlrpclib

PARTITION_SERVER_URL = "http://skuld.bmsc.washington.edu:10500"

class PartitionClient(object):
    def __init__(self, pserver_url):
        self.pserver_url = pserver_url
        self.pserver = xmlrpclib.ServerProxy(PARTITION_SERVER_URL)

    def calc_partitions(self, fileobj, npartition_range):
        pdb_bin = xmlrpclib.Binary(fileobj.read())
        return self.pserver.calc_partitions(
            pdb_bin, npartition_range[0], npartition_range[1])

def prnt_presult(presult):
    for chain_id, npart, partlist in presult:
        sys.stdout.write("<%s %d> %s\n" % (chain_id, npart, ", ".join(partlist)))
    

def usage():
    sys.stderr.write(
        "partition_client.py <pdb path> <num partitions begin> <num partitions end>\n\n")
    
def main():
    try:
        path = sys.argv[1]
        npart_begin = int(sys.argv[2])
        npart_end = int(sys.argv[3])
    except (IndexError, ValueError):
        usage()
        raise SystemExit

    if not os.path.isfile(path):
        sys.stderr.write("file not found: %s\n\n" % (path))
        raise SystemExit

    pclient = PartitionClient(PARTITION_SERVER_URL)
    presult = pclient.calc_partitions(open(path, "r"), (npart_begin, npart_end))
    if isinstance(presult, list):
        prnt_presult(presult)
    else:
        sys.stderr.write("server returned an error: %s\n\n" % (presult))
        raise SystemExit
    
if __name__ == "__main__":
    main()
