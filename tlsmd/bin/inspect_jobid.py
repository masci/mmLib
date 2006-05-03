#!/usr/bin/env python

import os
import sys
import time

import xmlrpclib
from tlsmdlib import const, conf

webtlsmdd = xmlrpclib.ServerProxy(conf.WEBTLSMDD)

def main():
    if len(sys.argv) != 2:
        print "Usage: ./inspect_jobid.py <jobid>"
        sys.exit()

    job_list = webtlsmdd.job_list()
    for jdict in job_list:
        if sys.argv[1] in jdict.values():
            for k,v in jdict.iteritems():
                print k, ":", v

if __name__ == "__main__":
    main()


