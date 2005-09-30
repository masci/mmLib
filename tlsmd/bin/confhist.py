#!/usr/bin/env python
## TLS Minimized Domains (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys

from mmLib.Structure   import *
from mmLib.FileLoader  import *
from tlsmdlib.datafile import *

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
        for i in range(p - 1):
            self.I.append(m * (i + 1))

    def calc_configurations(self):
        """
        calculate the number of configurations there should be
        """
        n = self.n
        m = self.m
        p = self.p
        return ((n - p*m + 1) * (1 + ( (p-2)*(n - p*m) + 1))) / 2

    def check_i(self):
        """
        make sure there are no partitions smaller than m
        """
        for i in range(1, len(self.I)):
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
        

class ConfResidHistorgram(SubSegConfs):
    def __init__(self, dbfile, chain, m, p):
        self.datafile = TLSMDFile(dbfile)
        self.chain = chain
        SubSegConfs.__init__(self, len(self.chain), m, p)

        self.clen = len(self.chain)

        self.ntop = 40
        self.top  = []

        self.rmat = self.calc_rmat()
        self.init_histogram()

    def calc_rmat(self):
        """
        create a residual matrix
        """
        print "Generating a residual matrix..."

        rmat = zeros((self.clen, self.clen), Float)

        for i in range(self.clen):
            for j in range(i + self.m - 1, self.clen):

                frag_id1 = self.chain.fragment_list[i].fragment_id
                frag_id2 = self.chain.fragment_list[j].fragment_id

                data = self.datafile.grh_get_tls_record(self.chain.chain_id, frag_id1, frag_id2)
                if data == None:
                    print "EEK! %s-%s" % (frag_id1, frag_id2)

                try:
                    rmat[i,j] = data["lsq_residual"]
                except KeyError:
                    pass

        print "Done."

        return rmat

    def init_histogram(self):
        """
        initalize histogram
        """
        self.nbins = 30
        self.hrange_max = self.rmat[0,self.clen-1]
        self.hrange_min = 0.7
        self.bsize = (self.hrange_max - self.hrange_min) / self.nbins
        self.bins = [0 for x in range(self.nbins)]

    def bin_value(self, val):
        bin = int((val - self.hrange_min) / self.bsize)
        self.bins[bin] += 1

    def prnt_histogram(self):
        print "Histogram"
        for i in range(len(self.bins)):
            bin_min = self.hrange_min + (i * self.bsize)
            bin_max = self.hrange_min + ((i+1) * self.bsize)
            print "[%f-%f]: %d" % (bin_min, bin_max, self.bins[i])

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

        self.calc_configuration_residual(segments)

    def calc_configuration_residual(self, segments):
        lsqr = 0.0

        for i, j in segments:
            lsqr += self.rmat[i,j]

        self.bin_value(lsqr)

        if len(self.top) < self.ntop:
            desc = "%s: %f" % (str(segments), lsqr)
            self.top.append((lsqr, desc))
        elif lsqr < self.top[-1][0]:
            del self.top[-1]
            desc = "%s: %f" % (str(segments), lsqr)
            self.top.append((lsqr, desc))
            self.top.sort()

    def prnt_top(self):
        for lsqr, desc in self.top:
            print desc

def protein_segment(chain):
    i = None
    j = None

    for frag in chain.iter_amino_acids():
        if i==None: i = chain.index(frag)
        j = chain.index(frag)

    return chain[i:j]

def main():
    dbfile = sys.argv[1]
    struct = LoadStructure(fil=sys.argv[2])
    chain_id = sys.argv[3]

    seg = protein_segment(struct.get_chain(chain_id))

    crh = ConfResidHistorgram(dbfile, seg, 5, 5)
    crh.go()
    crh.prnt_top()
    crh.prnt_histogram()


if __name__=="__main__":
    main()

