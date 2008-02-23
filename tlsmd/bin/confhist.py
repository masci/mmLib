#!/usr/bin/python
## TLS Motion Determination (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys
import copy

from mmLib.Structure   import *
from mmLib.FileLoader  import *
from tlsmdlib.datafile import *

def prnt_sep():
    print "======================================================================================="

class SubSegConfs(object):
    def __init__(self, n, m, p):
        """
        number of residues: n
        smallest segment length: m
        number of partitions: p
        """
        self.n = n
        self.m = m
        self.p = p

        ## keep track of the total number of configurations
        self.total = 0

        ## initalize partition indexing
        self.I = []
        for i in xrange(p - 1):
            self.I.append(m * (i + 1))

    def calc_configurations(self):
        """
        calculate the number of configurations there should be
        XXX: this is wrong!
        """
        n = self.n
        m = self.m
        p = self.p
        return ((n - p*m + 1) * (1 + ( (p-2)*(n - p*m) + 1))) / 2

    def check_i(self):
        """
        make sure there are no partitions smaller than m
        """
        for i in xrange(1, len(self.I)):
            assert self.I[i]-self.I[i-1] >= self.m

    def go(self):
        self.xxx2(0)

    def xxx2(self, i):
        """
        recurse over all possible configurations
        """
        nvertex = self.p - 1
        istart  = self.I[i]
        iend    = self.n - (self.m * (self.p - 1 - i))

        ## recurse
        if i < (self.p-2):
            icur = istart

            while icur <= iend:
                self.total += 1
                self.I[i] = icur
                self.process_configuration()

                self.xxx2(i+1)

                ## incriment position of current vertex and all
                ## vetexes which come after it
                icur += 1
                
                segm = 1
                j = i + 1
                while j < nvertex:
                    self.I[j] = icur + segm*self.m 
                    segm += 1
                    j += 1

        ## last moving vertex
        else:
            icur = istart
            while icur <= iend:
                self.total += 1
                self.I[i] = icur
                self.process_configuration()
                icur += 1
                
    def process_configuration(self):
        """
        override for more interesting behavior
        """
        print "I ",self.I

class ChainPartitonPointHistogram(object):
    def __init__(self, chain):
        self.chain = chain
        self.num_partitions = 0
        self.pp_names = ["" for i in xrange(len(chain)+1)]
        self.pp_freqs = [0  for i in xrange(len(chain)+1)]

        for i in xrange(1, len(chain)):
            self.pp_names[i] = "%s:%s" % (chain[i-1], chain[i])

    def count(self, i):
        self.pp_freqs[i] += 1

    def prnt(self):
        fil = open("cphistogram.txt","w")

        print "ChainPartitonPointHistogram"
        for i in xrange(len(self.chain)+1):
            print "[%s]: %10d" % (self.pp_names[i], self.pp_freqs[i])
            fil.write("%d %d %f\n" % (i, int(self.chain[i-1].fragment_id), self.pp_freqs[i]))

        fil.close()

class Histogram(object):
    def __init__(self, nbins, hmax, hmin = 0.0):
        self.nbins = nbins
        self.hmax = hmax
        self.hmin = hmin
        self.bsize = (self.hmax - self.hmin) / self.nbins
        self.bins = [0 for x in xrange(self.nbins)]
        
    def count(self, val):
        bin = int((val - self.hmin) / self.bsize)
        self.bins[bin] += 1

    def prnt(self):
        fil = open("histogram.txt","w")

        print "Histogram"
        for i in xrange(len(self.bins)):
            bin_min = self.hmin + (i * self.bsize)
            bin_max = self.hmin + ((i+1) * self.bsize)
            print "[%f-%f]: %d" % (bin_min, bin_max, self.bins[i])
            fil.write("%d %f\n" % (i, self.bins[i]))

        fil.close()

class TopList(list):
    def __init__(self, ntop):
        list.__init__(self)
        self.ntop = ntop
        self.tmax = None

    def add(self, val):
        if val < self.tmax:
            list.append(self, val)
            list.sort(self)
            self.tmax = list.__getitem__(self, -1)

class ConfResidHistorgram(SubSegConfs):
    def __init__(self, dbfile, chain, m, p):
        self.datafile = TLSMDFile(dbfile)
        self.chain = chain
        SubSegConfs.__init__(self, len(self.chain), m, p)

        self.clen = len(self.chain)
        self.ntop = 25
        self.top  = []

        self.rmat, self.rmatflag = self.calc_rmat()

        self.histogram    = Histogram(100, self.rmat[0,self.clen-1])
        self.cp_histogram = ChainPartitonPointHistogram(self.chain)

    def calc_rmat(self):
        """
        create a residual matrix
        """
        prnt_sep()
        print "Generating a residual matrix..."

        rmat = zeros((self.clen, self.clen), Float)
        rmatflag = zeros((self.clen, self.clen), Int)

        for i in xrange(self.clen):
            for j in xrange(i + self.m - 1, self.clen):

                frag_id1 = self.chain.fragment_list[i].fragment_id
                frag_id2 = self.chain.fragment_list[j].fragment_id

                data = self.datafile.grh_get_tls_record(self.chain.chain_id, frag_id1, frag_id2)
                if data == None or data.has_key("lsq_residual") == False:
                    print "No Database Record: %s-%s" % (frag_id1, frag_id2)
                else:
                    rmat[i,j] = data["lsq_residual"]
                    rmatflag[i,j] = 1

        print "Done."

        return rmat, rmatflag

    def process_configuration(self):
        segments = []

        i = 0
        j = 0        
        for v in self.I:
            i = j
            j = v
            segments.append((i,j-1))

        i = j
        j = self.n
        segments.append((i,j-1))

        if self.total%10000==0:
            print "[%d] %s" % (self.total, str(segments))

        self.calc_configuration_residual(segments)

    def calc_configuration_residual(self, segments):
        lsqr = 0.0

        for i, j in segments:
            if self.rmatflag[i,j]==0:
                return
            lsqr += self.rmat[i,j]

        self.histogram.count(lsqr)

        if len(self.top) < self.ntop:
            self.top.append((lsqr, copy.copy(self.I)))
        elif lsqr < self.top[-1][0]:
            del self.top[-1]
            self.top.append((lsqr, copy.copy(self.I)))
            self.top.sort()

    def prnt(self):
        prnt_sep()
        self.histogram.prnt()
        prnt_sep()

        for lsqr, I in self.top:
            for i in I:
                self.cp_histogram.count(i)

        self.cp_histogram.prnt()
        
    def prnt_top(self):
        for lsqr, desc in self.top:
            print desc

def protein_segment(chain):
    i = None
    j = None

    seg = Segment(chain_id = chain.chain_id)

    for frag in chain.iter_amino_acids():
        seg.add_fragment(frag)

    return seg

def main():
    dbfile = sys.argv[1]

    prnt_sep()
    print "Loading Structure..."
    struct = LoadStructure(fil=sys.argv[2])


    chain_id = sys.argv[3]
    seg = protein_segment(struct.get_chain(chain_id))

    crh = ConfResidHistorgram(dbfile, seg, 15, 5)

    prnt_sep()
    print "Analyizing all configurations..."
    crh.go()

    crh.prnt()

if __name__=="__main__":
    main()

