#!/usr/bin/env python

## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

from __future__ import generators
import sys
import math

from mmLib.Structure      import *
from mmLib.FileLoader     import *
from mmLib.Extensions.TLS import *


class Tabulator(object):
    def __init__(self):
        self.column_dicts = []
        self.rows = []

    def add_column(self, name, desc, width, format):
        column_dict = {
            "name":   name,
            "desc":   desc,
            "width":  width,
            "format": format}
        self.column_dicts.append(column_dict)

    def prnt_cols(self):
        print "# ",

        last = self.column_dicts[-1]
        
        for cdict in self.column_dicts:
            desc  = cdict["desc"]
            width = cdict["width"]

            print desc.ljust(width)[:width],

            if cdict!=last:
                print " ",

        print

    def add_row(self, row_dict):
        self.rows.append(row_dict)

    def prnt_row(self, row_dict):
        print "  ",

        last = self.column_dicts[-1]
        
        for cdict in self.column_dicts:
            name  = cdict["name"]
            width = cdict["width"]
            
            try:
                x = cdict["format"] % (row_dict[name])
            except KeyError:
                x = ""
                
            x = x.ljust(width)[:width]
            print x,

            if cdict!=last:
                print " ",

        print


def prnt_header(args):
    print "## PATH: %s" % (args["path"])
    print "## SEGMENT LENGTH: %d" % (args["seg_len"])
    print "## MAINCHAIN ONLY: %s" % (str(args["mainchain_only"]))
    print "## OMIT SINGLE BONDED ATOMS: %s" % (str(args["omit_single_bonded"]))


def prnt_stats(stats):
        ## print out results
        print str(stats["name"]).ljust(8),

        print str(stats["segment_num"]).ljust(5),

        x = "%d" % (stats["num_atoms"])
        print x.ljust(8),

        x = "%.3f" % (stats["mean_B"])
        print x.ljust(7),

        x = "%4.2f" % (stats["mean_A"])
        print x.ljust(5),

        x = "%.3f" % (stats["R"])
        print x.ljust(6),

        x = "%.4f" % (stats["mean_dp2"])
        print x.ljust(7),

        x = "%.4f" % (stats["mean_dp2_sigma"])
        print x.ljust(7),

        x = "%5.3f" % (stats["mean_S"])
        print x.ljust(6),

        x = "%5.3f" % (stats["mean_S_sigma"])
        print x.ljust(6),

        x = "%6.4f" % (trace(stats["rT'"]))
        print x.ljust(7),

        x = "%6.4f" % (trace(stats["L'"])*RAD2DEG2)
        print x.ljust(10),

        print

def mk_tab2():
    tab2 = Tabulator()
    tab2.add_column("name",         "ResRng" , 8, "%s")
    tab2.add_column("segment_num",  "SegNo" ,  5, "%d")
    tab2.add_column("num_atoms",    "Atoms",   6, "%d")
    tab2.add_column("mean_B",       "<B>",     6, "%6.2f")
    tab2.add_column("mean_A",       "<A>",     4, "%4.2f")
    tab2.add_column("R",            "R",       6, "%.3f")
    tab2.add_column("Tr",           "t(Tr)",   6, "%6.4f")
    tab2.add_column("L",            "tr(L)",   6, "%5.2f")
    tab2.add_column("lsqr",         "LSQR",    6, "%6.4f")
    tab2.add_column("pv_lsqr",      "LSQRpv",  6, "%6.4f")
    tab2.add_column("lsq_ratio",    "Ratio",   6, "%4.2f")
    tab2.add_column("pv_Tr",        "pv_t(Tr)",   8, "%6.4f")
    tab2.add_column("pv_L",         "pv_tr(L)",   8, "%5.2f")
    tab2.add_column("pv_num_atoms", "pv_Atoms",   8, "%d")
    
    return tab2

            
def prnt_stats2(tab, stats):

    lsqr = stats["lsq_residual"]
    pv_lsqr = stats["ca_pivot"]["lsq_residual"]
    
    pv_T = stats["ca_pivot"]["T"]
    pv_L = stats["ca_pivot"]["L"]

    row = {
        "name":           stats["name"],
        "segment_num":    stats["segment_num"],
        "num_atoms":      stats["num_atoms"],
        "mean_B":         stats["mean_B"],
        "mean_A":         stats["mean_A"],
        "R":              stats["R"],
        "Tr":             trace(stats["rT'"]),
        "L":              trace(stats["L'"])*RAD2DEG2,
        "lsqr":           lsqr,
        "pv_lsqr":        pv_lsqr,
        "lsq_ratio":      pv_lsqr / lsqr,
        "pv_Tr":          trace(pv_T),
        "pv_L":           trace(pv_L)*RAD2DEG2,
        "pv_num_atoms":   stats["ca_pivot"]["num_atoms"],
        }            

    tab.prnt_row(row)


def main(**args):
    tab = mk_tab2()

    prnt_header(args)
    tab.prnt_cols()

    struct = LoadStructure(fil=args["path"])
    tls_analysis = TLSStructureAnalysis(struct)

    tls_num = 0
    for tls_info in tls_analysis.iter_fit_TLS_segments(
        residue_width       = args["seg_len"],
        use_side_chains     = not args["mainchain_only"],
        include_single_bond = not args["omit_single_bonded"]):

        tls_num += 1
        tls_info["segment_num"] = tls_num
        
        ## calculate adverage temp factor and anisotropy
        tls_group = tls_info["tls_group"]

        Umean = 0.0
        Amean = 0.0
        for atm in tls_group:
            Umean += trace(atm.get_U())/3.0
            Amean += atm.calc_anisotropy()

        Umean = Umean / len(tls_group)
        Amean = Amean / len(tls_group)

        tls_info["mean_U"] = Umean 
        tls_info["mean_B"] = Umean * U2B
        tls_info["mean_A"] = Amean

        prnt_stats2(tab, tls_info)


def usage():
    print "search_tls.py - A utility to fit TLS groups to anisotropically"
    print "                or isotropically refined protein structures for"
    print "                motion analysis."
    print
    print "DESCRIPTION:"
    print "    Compute a range of TLS tensors by walking the amino"
    print "    acid backbone one residue at a time, spanning a continous"
    print "    sequence segment of a given length.  Each TLS calculation"
    print "    produces one line of output with some interesting statistics."
    print "    Please do not assume this is a scientifically useful"
    print "    thing to do!"
    print
    print "OPTIONS:"
    print "    -l <length>"
    print "        Set the length, in sequential amino acids, of the"
    print "        TLS groups which will be fit to the protein. The"
    print "        default is 3."
    print "    -m  Search mainchain atoms only"
    print "    -s  Omit atoms with only one bond"
    print "    -n  Omit TLS groups with negitive L/T Eigenvalues"
    print "    -c <chain_id>"
    print "        Only search the given chain."
    print


if __name__ == "__main__":
    import getopt

    try:
        (opts, args) = getopt.getopt(sys.argv[1:], "l:msnc:dx")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    ## program defaults
    seg_len            = 3
    mainchain_only     = False
    omit_single_bonded = False
    omit_neg_eigen     = False
    chain              = None
    use_dP2_min        = False
    cluster            = False

    ## read program options
    for (opt, item) in opts:
        if opt=="-l":
            try:
                seg_len = int(item)
            except ValueError:
                usage()
                sys.exit(1)

        if opt=="-m":
            mainchain_only = True

        if opt=="-s":
            omit_single_bonded = True

        if opt=="-n":
            omit_neg_eigen = True

        if opt=="-c":
            chain = item

        if opt=="-d":
            use_dP2_min = True

        if opt=="-x":
            cluster = True

    ## make sure a file name was entered
    if len(args)!=1:
        usage()
        sys.exit(1)

    main(path               = args[0],
         seg_len            = seg_len,
         mainchain_only     = mainchain_only,
         omit_single_bonded = omit_single_bonded,
         omit_neg_eigen     = omit_neg_eigen,
         chain              = chain,
         use_dP2_min        = use_dP2_min,
         cluster            = cluster)
