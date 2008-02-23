#!/usr/bin/python
## TLS Minimized Domains (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys

from tlsmdlib.datafile import *


def main():
    datafile = TLSMDFile(sys.argv[1])

    chain_id = sys.argv[2]
    frag_id1 = sys.argv[3]
    frag_id2 = sys.argv[4]

    data = datafile.grh_get_tls_record(chain_id, frag_id1, frag_id2)

    for key in data.keys():
        print "%10s = %s" % (key, str(data[key]))

if __name__=="__main__":
    main()
