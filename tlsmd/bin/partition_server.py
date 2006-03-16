#!/usr/bin/env python
## TLS Minimized Domains (TLSMD)
## Copyright 200-20052 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import os
import sys
import tempfile
import xmlrpclib
import SocketServer
import SimpleXMLRPCServer

from tlsmdlib import tlsmd_analysis

SERVER_HOST_PORT = ("localhost", 10500)

class Forking_XMLRPCServer(
            SocketServer.ForkingMixIn,
            SimpleXMLRPCServer.SimpleXMLRPCServer):
    """Use customized XMLRPC server which forks for requests and uses the customized
    request handler.
    """
    def __init__(self, host_port):
        SimpleXMLRPCServer.SimpleXMLRPCServer.__init__(
            self,
            host_port,
            SimpleXMLRPCServer.SimpleXMLRPCRequestHandler,
            True)

class PartitionServer(object):
    def calc_partitions(self, pdb_bin, prange1, prange2):
        olddir = os.getcwd()

        ## create temp directory and save files
        tempdir = tempfile.mkdtemp()
        os.chdir(tempdir)
        fd, path = tempfile.mkstemp(".pdb", "struct", tempdir, "w")
        fileobj = os.fdopen(fd, "w")
        fileobj.write(pdb_bin.data)
        fileobj.close()
        sys.stdout.write("tempfile: %s\n" % (path))

        ## TLSMD analysis
        anal = tlsmd_analysis.TLSMDAnalysis(struct_file_path = path)
        anal.run_optimization()

        presult = []
        
        for chain in anal.iter_chains():
            partition_collection = chain.partition_collection
            
            for cpartition in partition_collection.iter_chain_partitions():
                if cpartition.ntls >= prange1 and cpartition.ntls <= prange2:
                    partlist = []
                    presult.append((chain.chain_id, cpartition.ntls, partlist))

                    prev_seg = None
                    for seg in cpartition.iter_tls_segments():
                        if prev_seg is None:
                            prev_seg = seg
                            continue
                        partlist.append("%s:%s" % (prev_seg.frag_id2, seg.frag_id1))
                        prev_seg = seg

        ## remove temp directory
        os.chdir(olddir)
        for path in os.listdir(tempdir):
            os.remove(os.path.join(tempdir, path))
        os.rmdir(tempdir)

        return presult

def main():
    pserver = PartitionServer()
    xmlrpc_server = Forking_XMLRPCServer(SERVER_HOST_PORT)
    xmlrpc_server.register_instance(pserver)
    xmlrpc_server.serve_forever()
    
if __name__ == "__main__":
    main()
