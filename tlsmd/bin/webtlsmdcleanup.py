#!/usr/bin/python
## TLS Motion Determination (TLSMD)
## Copyright 2002-2010 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

## Pyton modules
import os
import sys
from time import time
from shutil import rmtree

## TLSMD
from tlsmdlib import conf, misc, mysql_support

SECS_IN_DAY = float(60*60*24)
DELETE_DAYS = 120

## GLOBALS
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

    days = round((time() - float(submit_time)) / SECS_IN_DAY)
    jdict["days"] = days

    if days > DELETE_DAYS:
        return True

    ## This is a temporary setting to clear out jobs submitted during
    ## development. DEBUG is normally set to False
    DEBUG = False
    if DEBUG:
        if state in ["died", "errors", "submit1"]:
            return True

    if state in ["submit1", "lost_directory", "syserror"] and days > 1:
        return True

    return False

def remove_job(job_id):
    """Removes the job from both the database and working directory.
    If job is still running when this function is called, it will first call
    KillJob(), then remove the associated data and files.
    """
    if not mysql.job_exists(job_id):
        return False

    try:
        job_dir = os.path.join(conf.TLSMD_WORK_DIR, job_id)
        rmtree(job_dir)
    except:
        return False

    mysql.delete_jdict(job_id)
    return True

def main():
    t = misc.timestamp()
    msg  = "[%s] WebTLSMD Job Cleanup: " % t
    msg += "Checking database for old jobs to remove."
    print msg

    job_list = mysql.job_list()
    jdict_remove_list = []

    for jdict in job_list:
        if jdict["via_pdb"] and jdict["state"] != ("running" or "queued"):
            print "[%s] saving to database PDB: %s" % (
                t, jdict["structure_id"])
            mysql.archive_pdb_jobs(jdict["job_id"])

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

        try:
            remove_job(jdict["job_id"])
        except:
            msg  = "[%s] ERROR: Could not remove files/directories " % t
            msg += "associated with job_id: %s" % jdict["job_id"]
            print msg


if __name__=="__main__":
    main()

