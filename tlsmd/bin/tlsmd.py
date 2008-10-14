#!/usr/bin/python
## TLS Motion Determination (TLSMD)
## Copyright 2002-2008 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys
import getopt
import traceback

from tlsmdlib import conf, tlsmd_analysis, email


def usage():
    print "tlsmd.py - search for TLS rigid domains in x-ray crystal"
    print "           structures refined with individual atomic"
    print "           B-factors (temperature factors)"
    print
    print "Command Line Usage:"
    print "  tlsmd.py  [-v ] verbose output "
    print "            [-c <comma-seperated chain IDS>] process only given chain IDs"
    print "            [-x <URL of WebTLSMDD XMLRPC Server>]"
    print "            [-j <Job ID of WebTLSMDD Job>]"
    print "            [-r <html report dir>] write HTML report to directory"
    print "            [-m <tls model>] Models: ISOT(default)/ANISO/NLISOT/NLANISO"
    print "            [-w <Weighting Model>] Models: NONE(default)/IUISO"
    print "            [-a <Atoms>] ALL(default)/MAINCHAIN"
    print "            [-i <struct_id>] Override struct_id in PDB file"
    print "            [-s] output Gnuplot plots using SVG (default=False)"
    print "            [-k] skip generating JMOL files (default=False)"
    print "            [-h] skip generating HTML files (default=False)"
    print "            [-t struct.pdb:chain_id ] compare TLS displacments with another structure"
    print "            [-e] recombine linear TLSMD segments to find disjoint TLS groups"
    print "            [-o <num_adjecent_residues>] ADP smoothing using the given number of residues (default=0)"
    print "            [-u <min_subsegment_size>] (default=4)"
    print "            [-n <num_segments>] (default=20)"
    print "            [-b] email traceback to this address if a Python exception occurs"
    print "            struct.pdb"
    print
    print "Authors:"
    print "  Jay Painter <jpaint@u.washington.edu>"
    print
    raise SystemExit


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
    chain_ids         = opt_dict.get("-c")

    conf.globalconf.recombination = opt_dict.has_key("-e")

    if opt_dict.has_key("-o"):
        try:
            nsmooth = int(opt_dict["-o"])
        except ValueError:
            print "[ERROR] -o argument must be a integer"
            usage()
        conf.globalconf.adp_smoothing = nsmooth

    ## Added. Christoph Champ, 2008-05-13
    if opt_dict.has_key("-u"):
        print opt_dict["-u"]
        try:
            min_subseg = int(opt_dict["-u"])
        except ValueError:
            print "[ERROR] -u argument must be an integer"
            usage()
        conf.globalconf.min_subsegment_size = min_subseg

    ## Added. Christoph Champ, 2008-06-09
    if opt_dict.has_key("-n"):
        print opt_dict["-n"]
        try:
            num_segs = int(opt_dict["-n"])
        except ValueError:
            print "[ERROR] -n argument must be an integer"
            usage()
        conf.globalconf.nparts = num_segs

    if opt_dict.has_key("-t"):
        tpath, tchain_id = opt_dict["-t"].split(":")
        conf.globalconf.target_struct_path = tpath
        conf.globalconf.target_struct_chain_id = tchain_id
        
    if opt_dict.has_key("-v"):
        conf.globalconf.verbose = True

    if opt_dict.has_key("-s"):
        conf.globalconf.use_svg = True

    ## Added. Christoph Champ, 2008-05-13
    if opt_dict.has_key("-k"):
        conf.globalconf.skip_jmol = True

    ## Added. Christoph Champ, 2008-08-05
    if opt_dict.has_key("-h"):
        conf.globalconf.skip_html = True
        
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

    try:
        tlsmd_analysis.TLSMD_Main(
            struct_file_path    = struct_path,
            sel_chain_ids       = chain_ids,
            html_report_dir     = opt_dict.get("-r"))
    except:
        if opt_dict.has_key("-b"):
            email.SendTracebackEmail("tlsmd.py traceback")
        raise


if __name__ == "__main__":
    try:
        ## added option "-u" = min_subsegment_size. Christoph Champ, 2008-05-13
        ## added option "-k" = skip_jmol. Christoph Champ, 2008-05-13
        ## added option "-n" = nparts. Christoph Champ, 2008-06-09
        ## used letters: abcdeijmnorstuvwx
        ## available   : fghlpqyz
        (opts, args) = getopt.getopt(sys.argv[1:], "n:u:a:t:c:d:i:w:m:r:j:x:khvseo:b")
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
	print "completed"
    except KeyboardInterrupt:
        print "Killed"
    except:
	print "Died"

    sys.exit(0)
