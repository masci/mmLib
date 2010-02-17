## TLS Motion Determination (TLSMD)
## Copyright 2002-2010 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

## Python modules
import sys
import getopt
import math
import random
import numpy

## pymmlib
from mmLib import FileIO, Structure, TLS, AtomMath, Constants

def usage():
    print "%s [-t <tlsin>] [-p <pdbout>] <structure file>"%(sys.argv[0])
    print
    print "description:"
    print "    Performs a least squares fit of TLS tensors to the"
    print "    tempature factors of the given structure file.  If"
    print "    no TLS groups are defined by the TLSIN file, then"
    print "    one group is created per chain."
    print


T1 = numpy.array([[10.0,  0.0,  0.0],
                  [ 0.0, 25.0,  0.0],
                  [ 0.0,  0.0, 15.0]], float) * Constants.B2U

L1 = numpy.array([[ 5.0,  0.0,  0.0],
                  [ 0.0,  6.0,  0.0],
                  [ 0.0,  0.0,  3.0]], float) * Constants.DEG2RAD2

S1 = numpy.array([[ 1.0,  0.0, -0.2],
                  [ 0.0, -2.0,  0.0],
                  [ 0.2,  0.0,  1.0]], float) * Constants.DEG2RAD

T2 = numpy.array([[10.0,  0.0,  0.0],
                  [ 0.0, 25.0,  0.0],
                  [ 0.0,  0.0, 15.0]], float) * Constants.B2U

L2 = numpy.array([[ 5.0,  0.0,  0.0],
                  [ 0.0,  8.0,  0.0],
                  [ 0.0,  0.0,  3.0]], float) * Constants.DEG2RAD2


T3 = numpy.array([[10.0,  0.0,  0.0],
                  [ 0.0, 25.0,  0.0],
                  [ 0.0,  0.0, 15.0]], float) * Constants.B2U

L3 = numpy.array([[ 5.0,  0.0,  0.0],
                  [ 0.0, 20.0,  0.0],
                  [ 0.0,  0.0,  3.0]], float) * Constants.DEG2RAD2



def random_vec():
    rcoord = lambda: random.random() - 0.5
    return AtomMath.normalize(numpy.array([rcoord(), rcoord(), rcoord()], float))


def rt_random(T):
    theta = (random.random() - 0.5) * math.pi
    R = AtomMath.rmatrixu(random_vec(), theta)
    return numpy.dot(numpy.dot(R, T), numpy.transpose(R))
    
def rt_random2(T1, T2):
    theta = (random.random() - 0.5) * math.pi
    R = AtomMath.rmatrixu(random_vec(), theta)
    T1R = numpy.dot(numpy.dot(R, T1), numpy.transpose(R))
    T2R = numpy.dot(numpy.dot(R, T2), numpy.transpose(R))
    return T1, T2
    
def main(path, opt_dict):
    struct = FileIO.LoadStructure(file = path)

    tls_group_list = []

    ## make the TLS groups
    if opt_dict.has_key("-t"):

        try:
            fil = open(opt_dict["-t"], "r")
        except IOError:
            print "ERROR: TLSIN File not found %s" % (opt_dict["-t"])
            sys.exit(-1)
        
        tls_file = TLS.TLSFile()
        tls_file.set_file_format(TLS.TLSFileFormatTLSOUT())
        tls_file.load(fil)
        
        for tls_desc in tls_file.tls_desc_list:
            #tls = tls_desc.generate_tls_group(struct) ## old def
            tls = tls_desc.construct_tls_group_with_atoms(struct)
            tls.tls_desc = tls_desc
            tls_group_list.append(tls)
            print tls.name

    else:
        tls_desc = TLS.TLSGroupDesc()
        tls_desc.add_range("A", "1", "A", "25", "ALL")
        tls_group = tls_desc.construct_tls_group_with_atoms(struct)
        tls_group_list.append(tls_group)
        tls_group.tls_desc = tls_desc

        tls_group.origin = tls_group.calc_centroid()
        tls_group.T = rt_random(T1)
        tls_group.L, tls_group.S = rt_random2(L1, S1)

        tls_desc = TLS.TLSGroupDesc()
        tls_desc.add_range("A", "26", "A", "35", "ALL")
        tls_group = tls_desc.construct_tls_group_with_atoms(struct)
        tls_group_list.append(tls_group)
        tls_group.tls_desc = tls_desc

        tls_group.origin = tls_group.calc_centroid()
        tls_group.T = rt_random(T2)
        tls_group.L, tls_group.S = rt_random2(L2, S1)
        
        tls_desc = TLS.TLSGroupDesc()
        tls_desc.add_range("A", "36", "A", "52", "ALL")
        tls_group = tls_desc.construct_tls_group_with_atoms(struct)
        tls_group_list.append(tls_group)
        tls_group.tls_desc = tls_desc

        tls_group.origin = tls_group.calc_centroid()
        tls_group.T = rt_random(T3)
        tls_group.L, tls_group.S = rt_random2(L3, S1)

    for atm in struct.iter_all_atoms():
        atm.U = None

    ## write out ideal TLS ANISOU records
    if opt_dict.has_key("-p"):
        for tls_group in tls_group_list:
            for atm, Utls in tls_group.iter_atm_Utls():
                atm.temp_factor = Constants.U2B * (numpy.trace(Utls) / 3.0)
                atm.U = Utls
                
        ## save the struct
        FileIO.SaveStructure(file = opt_dict["-p"], struct = struct)



if __name__ == "__main__":
    try:
        (opts, args) = getopt.getopt(sys.argv[1:], "t:p:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    opt_dict = {}
    for (flag, data) in opts:
        opt_dict[flag] = data

    try:
        path = args[0]
    except IndexError:
        usage()
        sys.exit(1)

    main(path, opt_dict)
