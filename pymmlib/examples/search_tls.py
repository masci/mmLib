#!/usr/bin/env python

## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

from __future__ import generators
import sys
import math

from mmLib.Structure import *
from mmLib.FileLoader import *
from mmLib.Extensions.TLS import *


def main(**args):
    print "## PATH: %s" % (args["path"])
    print "## SEGMENT LENGTH: %d" % (args["seg_len"])
    print "## MAINCHAIN ONLY: %s" % (str(args["mainchain_only"]))
    print "## OMIT SINGLE BONDED ATOMS: %s" % (str(args["omit_single_bonded"]))
    print "## RES   NUM   Atoms    <B>     <A>   R      <DP2>   s<DP2>  <DP2N>  "\
          "s<DP2N> <S>    s<S>   t(T)    t(L)"

    struct = LoadStructure(fil = args["path"])

    tls_analysis = TLSStructureAnalysis(struct)

    stats_list = tls_analysis.fit_TLS_segments(
        residue_width       = args["seg_len"],
        use_side_chains     = not args["mainchain_only"],
        include_single_bond = not args["omit_single_bonded"])

    
    for stats in stats_list:

        tls = stats["tls"]

        ## if a chain is specified, then skip all other chains
##         if args.get("chain")!=None and args.get("chain")!=chain.chain_id:
##             continue

        ## calculate adverage temp factor and anisotropy
        Umean = 0.0
        Amean = 0.0
        for atm in stats["tls"]:
            Umean += trace(atm.get_U())/3.0
            Amean += atm.calc_anisotropy()

        Umean = Umean / float(len(tls))
        Bmean = Umean * 8.0 * math.pi**2
        Amean = Amean / float(len(tls))

        ## print out results
        print str(stats["name"]).ljust(8),

        print str(stats_list.index(stats)).ljust(5),

        x = "%d" % (len(tls))
        print x.ljust(8),

        x = "%.3f" % (Bmean)
        print x.ljust(7),

        x = "%4.2f" % (Amean)
        print x.ljust(5),

        x = "%.3f" % (stats["R"])
        print x.ljust(6),

        x = "%.4f" % (stats["mean_DP2"])
        print x.ljust(7),

        x = "%.4f" % (stats["sigma_DP2"])
        print x.ljust(7),

        x = "%.4f" % (stats["mean_DP2N"])
        print x.ljust(7),

        x = "%.4f" % (stats["sigma_DP2N"])
        print x.ljust(7),
        
        x = "%5.3f" % (stats["mean_S"])
        print x.ljust(6),

        x = "%5.3f" % (stats["sigma_S"])
        print x.ljust(6),

        x = "%6.4f" % (trace(tls.T))
        print x.ljust(7),

        x = "%6.4f" % (trace(tls.L)*rad2deg2)
        print x.ljust(10),

        print
            

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
        (opts, args) = getopt.getopt(sys.argv[1:], "l:msnc:d")
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
         use_dP2_min        = use_dP2_min)
