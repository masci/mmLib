#!/usr/bin/env python

## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import sys
import math

from mmLib.Structure import *
from mmLib.FileLoader import *
from mmLib.Extensions.TLS import *

def usage():
    print "search_tls.py <file path>"
    print
    print "description:"
    print "    Compute a range of TLS tensors by walking the amino"
    print "    acid backbone one atom at a time, spanning a continous"
    print "    segment of 9 atoms (3 residues).  Each TLS calculation"
    print "    produces one line of output with some useful statistics."
    print "    Please do not assume this is a scientifically useful"
    print "    thing to do!"
    print


def iter_mainchain(chain, atom_len):
    """Yields a list of continous mainchain atoms length of atom_len, walking
    the chain one atom at a time.
    """
    extra_dict = {}
    atom_list  = []    

    for res in chain.iter_amino_acids():

        for name in ["N", "CA", "C"]:
            try:
                atm = res[name]
            except KeyError:
                raise StopIteration

            atom_list.append(atm)

            ## add side chain atoms
            if atm.name == "CA":
                extra_dict[atm] = []
                for atmx in res:
                    if atmx.name not in ["N", "CA", "C", "O"]:
                        extra_dict[atm].append(atmx)

            ## add C-O oxygen atoms
            if atm.name == "C":
                try:
                    atm0 = res["O"]
                except KeyError:
                    pass
                else:
                    extra_dict[atm] = [atm0]
                    
            if len(atom_list) >= atom_len:
                if  len(atom_list) > atom_len:
                    atm0 = atom_list[0]
                    atom_list = atom_list[1:]
                    try:
                        del extra_dict[atm0]
                    except KeyError:
                        pass

                ## add extra atoms to the yield list
                yield_list = []
                for atmx in atom_list:
                    yield_list.append(atmx)
                    try:
                        yield_list += extra_dict[atmx]
                    except KeyError:
                        pass

                yield yield_list


def iter_mainchain2(chain, atom_len):
    """Yields a list of continous mainchain atoms length of atom_len, walking
    the chain one atom at a time.
    """
    extra_dict = {}
    atom_list  = []    

    for res in chain.iter_amino_acids():

        for name in ["N", "CA", "C"]:
            try:
                atm = res[name]
            except KeyError:
                raise StopIteration

            atom_list.append(atm)

            ## add C-O oxygen atoms
            if atm.name == "C":
                try:
                    atm0 = res["O"]
                except KeyError:
                    pass
                else:
                    extra_dict[atm] = [atm0]
                    
            if len(atom_list) >= atom_len:
                if  len(atom_list) > atom_len:
                    atm0 = atom_list[0]
                    atom_list = atom_list[1:]
                    try:
                        del extra_dict[atm0]
                    except KeyError:
                        pass

                ## add extra atoms to the yield list
                yield_list = []
                for atmx in atom_list:
                    yield_list.append(atmx)
                    try:
                        yield_list += extra_dict[atmx]
                    except KeyError:
                        pass

                yield yield_list

def segment_list(chain, seg_len):
    """Simple amino acid segment generator.  Given a chain and segment length,
    this returns a list of all continuos segments of amino acids of that length
    in the chain, starting at the beginning of the chain.
    """
    segment = []
    segment_list = []

    for res in chain.iter_amino_acids():
        segment.append(res)

        if len(segment) < seg_len:
            continue

        segment_list.append(segment)
        segment = segment[1:]

    return segment_list


def main(path):
    print """
    Calculating TLS parameters for a single rigid body group composed of
    all the amino acids
    """
    print "# <group num> <num atoms> <Badv> <R> <R/Uadv> <trT> <trL> <trS>"

    struct = LoadStructure(fil = path)

    ## list of all TLS groups
    tls_list = []

    for chain in struct.iter_chains():
        for seg in iter_mainchain2(chain, 3*6):

            atm0 = seg[0]
            atmX = seg[-1]

            name = str(atm0.fragment) + atm0.name + \
                   str(atmX.fragment) + atmX.name 

            ## new tls group for segment
            tls = TLSGroup()
            tls_list.append(tls)

            tls.seg = seg            

            ## add segment atoms
            for atm in seg:
                tls.append(atm)

            tls.origin = tls.calc_centroid()

            ## calculate tensors and print
            tls.calc_tls_tensors()

            Rfact = tls.calc_R()
            calcs = tls.calc_COR()

            trT   = trace(calcs["T'"])/3.0
            trL   = trace(calcs["L'"])/3.0

            trS   = ( abs(calcs["S'"][0,0]) + abs(calcs["S'"][1,1]) +\
                      abs(calcs["S'"][2,2]) ) / 3.0

            Uadv = 0.0
            for atm in tls:
                if atm.U != None:
                    Uadv += trace(atm.U)/3.0

            Uadv = Uadv / float(len(tls))

            ## print out results
            #print str(name).ljust(35),

            print str(tls_list.index(tls)).ljust(5),

            x = "%d" % (len(tls))
            print x.ljust(8),

            x = "%.3f" % (Uadv * 8 * math.pi * math.pi)
            print x.ljust(10),

            x = "%.3f" % (Rfact)
            print x.ljust(8),

            x = "%.3f" % (Rfact/Uadv)
            print x.ljust(10),
            
            x = "%.3f" % (trT)
            print x.ljust(12),
            x = "%.3f" % (trL * rad2deg2)
            print x.ljust(12),
            x = "%.3f" % (trS * rad2deg)
            print x.ljust(12),

            print
            

if __name__ == "__main__":
    try:
        path = sys.argv[1]
    except IndexError:
        usage()
        sys.exit(1)

    main(path)
