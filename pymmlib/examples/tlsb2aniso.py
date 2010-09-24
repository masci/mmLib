#!/usr/bin/env python
## Copyright 2002-2010 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import sys
import getopt
import numpy

from mmLib import Constants, Structure, FileIO, TLS

def usage():
    print "tlsb2aniso.py <tlsin> <struct-in> <struct-out>"
    print
    print "description:"
    print "    Adds the TLS model ADPs to the struct-in ADPs and"
    print "    writes the result to struct-out."
    print

def main(tlsin, struct_in, struct_out):
    struct = FileIO.LoadStructure(file = struct_in)

    ## make the TLS groups
    fil = open(tlsin, "r")
    tls_file = TLS.TLSFile()
    tls_file.set_file_format(TLS.TLSFileFormatTLSOUT())
    tls_file.load(fil)
    
    tls_group_list = []
    for tls_desc in tls_file.tls_desc_list:
        tls = tls_desc.construct_tls_group_with_atoms(struct)
        tls.tls_desc = tls_desc
        tls_group_list.append(tls)
        print tls.origin

    for tls in tls_group_list:
        for atm, Utls in tls.iter_atm_Utls():
            atm.U = Utls + (Constants.B2U * atm.temp_factor * numpy.identity(3, float))
            atm.temp_factor = Constants.U2B * (numpy.trace(atm.U)/3.0)

    ## save the struct
    FileIO.SaveStructure(file = struct_out, struct = struct)


if __name__ == "__main__":
    try:
        tlsin, struct_in, struct_out = sys.argv[1:]
    except:
        usage()
        sys.exit(1)

    main(tlsin, struct_in, struct_out)
