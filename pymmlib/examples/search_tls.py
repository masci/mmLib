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


def iter_segments(chain, seg_len):
    """
    """
    segment = []
    for res in chain.iter_amino_acids():
        segment.append(res)

        if len(segment)<seg_len:
            continue
        
        if len(segment)>seg_len:
            segment = segment[1:]

        atom_list = AtomList()
        for rx in segment:
            for atm in rx.iter_atoms():
                atom_list.append(atm)

        yield atom_list


def main(**args):
    print "## PATH: %s" % (args["path"])
    print "## SEGMENT LENGTH: %d" % (args["seg_len"])
    print "## MAINCHAIN ONLY: %s" % (str(args["mainchain_only"]))
    print "## OMIT SINGLE BONDED ATOMS: %s" % (str(args["omit_single_bonded"]))
    print "## res range::group num::num atoms::Badv::Aadv::R::dP2::Suij"

    struct = LoadStructure(fil = args["path"])

    ## list of all TLS groups
    tls_list = []

    for chain in struct.iter_chains():

        ## if a chain is specified, then skip all other chains
        if args.get("chain")!=None and args.get("chain")!=chain.chain_id:
            continue

        for seg_atom_list in iter_segments(chain, args["seg_len"]):

            atm0 = seg_atom_list[0]
            atmX = seg_atom_list[-1]
            name = "%s-%s" % (atm0.fragment_id, atmX.fragment_id)

            ## new tls group for segment
            tls = TLSGroup()

            ## filter atoms being added to the group
            ## add all atoms to the TLS group which are at full occupancy
            for atm in seg_atom_list:
                if atm.element=="H":
                    continue

                if atm.occupancy<1.0:
                    continue

                if args["mainchain_only"]==True and\
                   atm.name not in ("N", "CA", "C", "O"):
                    continue

                if args["omit_single_bonded"]==True and\
                   len(atm.bond_list)<=1:
                    continue
                
                tls.append(atm)
            ##

            if len(tls)==0:
                print "## empty group"
                continue

            tls_list.append(tls)
            tls.origin = tls.calc_centroid()

            ## calculate tensors and print
            if args.get("use_dP2_fit")==True:
                tls.calc_TLS_least_squares_fit()
            else:
                tls.calc_TLS_dP2_fit()

            if args["omit_neg_eigen"]==True:
                if min(eigenvalues(tls.T))<=0.0:
                    print "## Negitive T Eigenvalue"
                    continue
                if min(eigenvalues(tls.L))<=0.0:
                    print "## Negitive L Eigenvalue"
                    continue
            
            calcs = tls.shift_COR()
            Rfact = tls.calc_R()
            dP2   = tls.calc_adv_DP2uij()
            Suij  = tls.calc_adv_Suij()

            ## calculate adverage temp factor and anisotropy
            Uadv = 0.0
            Aadv = 0.0
            for atm in tls:
                Uadv += trace(atm.get_U())/3.0
                Aadv += atm.calc_anisotropy()

            Uadv = Uadv / float(len(tls))
            Badv = Uadv * 8.0 * math.pi**2
            Aadv = Aadv / float(len(tls))

            ## print out results
            print str(name).ljust(8),

            print str(tls_list.index(tls)).ljust(5),

            x = "%d" % (len(tls))
            print x.ljust(8),

            x = "%.3f" % (Badv)
            print x.ljust(10),

            x = "%4.2f" % (Aadv)
            print x.ljust(10),
            
            x = "%.3f" % (Rfact)
            print x.ljust(8),

            x = "%.4f" % (dP2)
            print x.ljust(10),

            x = "%5.3f" % (Suij)
            print x.ljust(10),
            
            x = "%6.4f" % (trace(tls.T))
            print x.ljust(10),

            x = "%6.4f" % (trace(tls.L)*rad2deg2)
            print x.ljust(10),


            eval = eigenvalues(tls.L)

            x = "%6.4f" % ((math.sqrt(eval[0]) +
                            math.sqrt(eval[1]) +
                            math.sqrt(eval[2]))*rad2deg)

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
