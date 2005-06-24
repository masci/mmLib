#!/home/jpaint/local/bin/python
## TLS Minimized Domains (TLSMD)
## Copyright 2002 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import sys

from tlsmdlib import misc
from tlsmdlib import tlsmd_analysis

def usage():
    print """\
tlsmdgridd.py - run tlsmd as a XMLRPC server daemon for

Command Line Usage:
    tlsmdgridd.py <port number>

Authors:
    Jay Painter <jpaint@u.washington.edu>

"""
    sys.exit(1)

if __name__ == "__main__":
    try:
        port = int(sys.argv[1])
    except (IndexError, ValueError):
        usage()

    print "TLSMDGridD Starting on Port %d" % (port)

    server = tlsmd_analysis.TLSGraphChainXMLRPCServer()

    try:
        server.run_server(("localhost", port))
    except KeyboardInterrupt:
        print "Killed"

    sys.exit(0)
