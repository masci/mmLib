#!/home/tlsmd/local/bin/python

## TLS Motion Determination (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys
import getopt

from tlsmdlib import conf
from tlsmdlib import tlsmd_analysis
from tlsmdlib import html

def usage():
    print "tlsmd.py - search for TLS rigid domains in x-ray crystal"
    print "           structures refined with individual atomic"
    print "           B-factors (temperature factors)"
    print
    print "Command Line Usage:"
    print "  tlsmd.py  [-v ] verbose output "
    print "            [-x <URL of WebTLSMDD XMLRPC Server>]"
    print "            [-j <Job ID of WebTLSMDD Job>]"
    print "            [-r <html report dir>]"
    print "            [-d <TLS database file>]"
    print "            [-n] TLS database is complete"
    print "            [-m <tls model>] Models: ISOT(default)/ANISO/NLISOT/NLANISO"
    print "            [-w <Weighting Model>] Models: NONE(default)/IUISO"
    print "            [-a <Atoms>] ALL(default)/MAINCHAIN"
    print "            [-i <struct_id>] Override struct_id in PDB file"
    print "            [-s] output Gnuplot plots using SVG"
    print "            [-t struct.pdb:chain_id ]"
    print "            struct.pdb"
    print
    print "Authors:"
    print "  Jay Painter <jpaint@u.washington.edu>"
    print
    sys.exit(1)

###############################################################################
## main() functions for the various modes of execution
##

def analysis_main(struct_path, opt_dict):
    """Runs the TLS analysis on a structure.
    """
    ## get source filename split up for use in
    ## constructing output file names
    struct_path         = os.path.realpath(struct_path)
    
    directory, filename = os.path.split(struct_path)
    basename, ext       = os.path.splitext(filename)

    ## set option vars and defaults
    tlsdb_file        = opt_dict.get("-d")
    tlsdb_complete    = opt_dict.has_key("-n")
    chain_ids         = opt_dict.get("-c")

    if opt_dict.has_key("-t"):
        tpath, tchain_id = opt_dict["-t"].split(":")
        conf.globalconf.target_struct_path = tpath
        conf.globalconf.target_struct_chain_id = tchain_id
        
    if opt_dict.has_key("-v"):
        conf.globalconf.verbose = True

    if opt_dict.has_key("-s"):
        conf.globalconf.use_svg = True
        
    if opt_dict.has_key("-x"):
        conf.globalconf.webtlsmdd = opt_dict["-x"]

    if opt_dict.has_key("-j"):
        conf.globalconf.job_id = opt_dict["-j"]

    if opt_dict.has_key("-i"):
        conf.globalconf.struct_id = opt_dict["-i"]
        
    ## set the TLS model to use
    if opt_dict.has_key("-m"):
        tls_model = opt_dict.get("-m").upper()
        if tls_model in ["ANISO", "ISOT", "NLANISO", "NLISOT"]:
            conf.globalconf.tls_model = tls_model
        else:
            print "[ERROR] Invalid TLS Model: %s" % (tls_model)
            usage()

    ## weighting scheme
    if opt_dict.has_key("-w"):
        val = opt_dict["-w"].upper()
        val = val.upper()
        if val=="IUISO":
            conf.globalconf.weight_model = val
        else:
            print "[ERROR] Invalid Weight Model: %s" % (val)
            usage()

    ## atoms to include
    if opt_dict.has_key("-a"):
        val = opt_dict["-a"].upper()
        if val not in ["ALL", "MAINCHAIN"]:
            usage()
        conf.globalconf.include_atoms = val

    ## create the analysis processor and load the structure, select chains
    anal = tlsmd_analysis.TLSMDAnalysis(
        struct_path    = struct_path,
        sel_chain_ids  = chain_ids,
        tlsdb_file     = tlsdb_file,
        tlsdb_complete = tlsdb_complete)
    
    anal.run_optimization()

    ## generate HTML report if a directory is given
    if opt_dict.has_key("-r") and len(anal.chains)>0:
        report_dir = opt_dict["-r"]
        report = html.HTMLReport(anal)
        report.write(report_dir)

if __name__ == "__main__":
    try:
        (opts, args) = getopt.getopt(sys.argv[1:], "a:t:c:d:i:w:m:r:j:x:nvs")
    except getopt.GetoptError:
        usage()

    opt_dict = {}
    for (flag, data) in opts:
        opt_dict[flag] = data

    try:
        path = args[0]
    except IndexError:
        usage()

    try:
        analysis_main(path, opt_dict)
    except KeyboardInterrupt:
        print "Killed"

    sys.exit(0)
