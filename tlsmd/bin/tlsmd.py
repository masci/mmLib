#!/home/jpaint/local/bin/python
## TLS Minimized Domains (TLSMD)
## Copyright 2002 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys
import string
import getopt

from tlsmdlib import misc
from tlsmdlib import tlsmd_analysis
from tlsmdlib import html

def usage():
    print "tlsmd.py - search for TLS rigid domains in x-ray crystal"
    print "           structures refined with individual atomic"
    print "           B-factors (temperature factors)"
    print
    print "Command Line Usage:"
    print "  tlsmd.py  [-t <num threads>]"
    print "            [-f <grid server config file>]"
    print "            [-x <URL of WebTLSMDD XMLRPC Server>]"
    print "            [-j <Job ID of WebTLSMDD Job>]"
    print "            [-r <html report dir>]"
    print "            [-d <tls database file>]"
    print "            [-m <tls model> Models: HYBRID(default)/ANISO]"
    print "            [-w <Weighting Model> Models: NONE(default)/IUISO ]"
    print "            [-a <Atoms> ALL(default)/MAINCHAIN ]"
    print "            [-i <struct_id> Override struct_id in PDB file]"
    print "            struct.pdb"
    print
    print "To run a a tls search in grid computation mode requires a"
    print "configuration file with the URLs of the grid computation"
    print "servers. Here is a example:"
    print 
    print "## WebTLSMD Grid Server File"
    print "http://localhost:9500"
    print "http://localhost:9501"
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
    database_file     = opt_dict.get("-d", "%s.db" % (basename))
    chain_id          = opt_dict.get("-c")
    grid_config_file  = opt_dict.get("-f")
    num_threads       = int(opt_dict.get("-t", 1))
    tls_model         = opt_dict.get("-m")

    if opt_dict.has_key("-x"):
        misc.GLOBALS["WEBTLSMDD"] = opt_dict["-x"]

    if opt_dict.has_key("-j"):
        misc.GLOBALS["JOB_ID"] = opt_dict["-j"]

    if opt_dict.has_key("-i"):
        misc.GLOBALS["STRUCT_ID"] = opt_dict["-i"]
        
    ## set the TLS model to use
    if tls_model==None:
        tls_model = "HYBRID"
    tls_model = tls_model.upper()
    if tls_model in ["ANISO", "HYBRID", "ISO"]:
        tlsmd_analysis.TLS_MODEL = tls_model
    else:
        print "[ERROR] Invalid TLS Model: %s" % (tls_model)
        usage()

    ## weighting scheme
    for opt, val in opt_dict.items():
        if opt=="-w":
            val = val.upper()

            if val=="NONE":
                continue
            
            if val in ["IUISO"]:
                tlsmd_analysis.WEIGHTS.append(val)
            else:
                print "[ERROR] Invalid Weight: %s" % (val)
                usage()

    ## atoms to include
    if opt_dict.has_key("-a"):
        val = opt_dict["-a"].upper()
        if val not in ["ALL", "MAINCHAIN", "CA"]:
            usage()
        tlsmd_analysis.set_include_atoms(val)

    ## create the analysis processor and load the structure, select chains
    anal = tlsmd_analysis.TLSMDAnalysis(struct_path)
    anal.load_struct()
    anal.select_chains(chain_id)

    anal_chains = []
    for chain in anal.chains:
        anal_chains.append(chain.chain_id)
    chnx = string.join(anal_chains, ", ")

    print "PROGRAM SETTINGS"
    print "    STRUCTURE FILE ====================: %s"%(anal.struct_path)
    print "    STRUCTURE ID ----------------------: %s"%(anal.struct_id)
    print "    CHAIN ID SELECTED FOR ANALYSIS ====: %s"%(chnx)
    print "    DATABASE FILE PATH ----------------: %s"%(database_file)
    print "    GRID SERVER CONFIG FILE ===========: %s"%(str(grid_config_file))
    print

    ## database chains
    anal.calc_tls_segments(database_file, num_threads, grid_config_file)
    anal.calc_chain_minimization()

    ## generate HTML report if a directory is given
    if opt_dict.has_key("-r"):
        report_dir = opt_dict["-r"]
        report = html.HTMLReport(anal)
        report.write(report_dir)

if __name__ == "__main__":
    import getopt

    try:
        (opts, args) = getopt.getopt(sys.argv[1:], "a:c:d:i:w:m:r:f:t:j:x:")
    except getopt.GetoptError:
        usage()

    opt_dict = {}
    for (flag, data) in opts:
        opt_dict[flag] = data

    try:
        path = args[0]
    except IndexError:
        usage()
    else:
        analysis_main(path, opt_dict)

