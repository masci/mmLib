#!/usr/bin/python

import os
import sys
import time

import xmlrpclib
from tlsmdlib import const, conf

webtlsmdd = xmlrpclib.ServerProxy(conf.WEBTLSMDD)
DUMP_ALL = True

def main():
    if len(sys.argv) != 2:
        print "Usage: ./inspect_jobid.py <jobid>"
        sys.exit()

    job_list = webtlsmdd.job_list()
    if DUMP_ALL:
        for jdict in job_list:
            for k,v in jdict.iteritems():
                print k, ":", v
    elif DUMP_ALL == False:
        for jdict in job_list:
            if sys.argv[1].strip("/") in jdict.values():
                for k,v in jdict.iteritems():
                    print k, ":", v

if __name__ == "__main__":
    main()


