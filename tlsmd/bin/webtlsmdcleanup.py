#!/usr/bin/python
## TLS Motion Determination (TLSMD)
## Copyright 2002-2009 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

## Pyton modules
import os
import sys
import time
import socket
import xmlrpclib

## TLSMD
from tlsmdlib import const, conf, misc, mysql_support

SECS_IN_DAY = float(60*60*24)
DELETE_DAYS = 14

## GLOBALS
webtlsmdd = xmlrpclib.ServerProxy(conf.WEBTLSMDD)
mysql = mysql_support.MySQLConnect()

def check_remove(jdict):
    state = jdict.get("state")
    if state == None:
        return True

    submit_time = jdict.get("submit_time")
    if submit_time == None:
        return False

    email = jdict.get("email")
    if email: pass

    days = round((time.time() - float(submit_time)) / SECS_IN_DAY)
    jdict["days"] = days

    if days > DELETE_DAYS:
        return True

    if state == "submit1" and days > 1:
        return True

    if state == "lost_directory" and days > 1:
        return True

    return False

def main():
    t = misc.timestamp()
    print "[%s] WebTLSMD Job Cleanup: Checking database for old jobs to remove." % t

    job_list = mysql.job_list()
    jdict_remove_list = []

    for jdict in job_list:
        if check_remove(jdict):
            jdict_remove_list.append(jdict)

    for jdict in jdict_remove_list:
        try:
            print "[%s] %10s  %40s  %d Days Old" % (
                t, jdict["job_id"], jdict.get("email", "No Email"), 
                jdict["days"])
        except:
            print "[%s] %10s  Bad submission status" % (t, jdict["job_id"])

        mysql.archive_old_jobs(jdict["job_id"])
        webtlsmdd.remove_job(jdict["job_id"])


if __name__=="__main__":
    main()
