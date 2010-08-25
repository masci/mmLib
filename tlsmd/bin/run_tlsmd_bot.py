#!/bin/env python
# -*- coding: utf-8 -*-
## TLS Minimized Domains (TLSMD)
## Copyright 2002-2010 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

## Python modules
import os
import sys
import time
import socket
import string
import random
import math
import numpy
import re
import xmlrpclib
import subprocess

## Pymmlib
from mmLib import Library ## checks if is_{amino,nucleic}_acid()

## TLSMD
from tlsmdlib import conf, console, const, misc, mysql_support

## GLOBALS
webtlsmdd = xmlrpclib.ServerProxy(conf.WEBTLSMDD)
mysql = mysql_support.MySQLConnect()

## JOB SELECT (example):
# mysql -B -N -e 'select id from pdb.remarks where tlsmd IS NULL order by rand() limit 1;'
# for i in `mysql -B -N -e 'select pdb_id from tlsmddb.via_pdb;'`; do mysql -e "UPDATE pdb.remarks SET tlsmd='1' WHERE id='$i';"; done

def timestring(secs):
    tm_struct = time.localtime(secs)
    return time.strftime("%Y-%m-%d %H:%M %Z", tm_struct)

def secdiffstring(secs):
    secs = int(secs)

    hours = secs / 3600
    secs = secs - (hours * 3600)

    min = secs / 60
    secs = secs - (min * 60)

    x = "%1d:%2d.%2d" % (hours, min, secs)
    return x.replace(" ", "0")

def timediffstring(begin, end):
    secs = int(end - begin)
    return secdiffstring(secs)

def left_justify_string(keyword, value):
    """Returns a string with dotted separation.
    """
    return '%s' % keyword .ljust(40, ".") + ": " + '%s\n' % value

def check_job_id(form):
    """Retrieves and confirms the job_id from a incomming form. Returns
    None on error, or the job_id on success.
    """
    if form.has_key("job_id"):
        job_id = form["job_id"].value
        if len(job_id) < conf.MAX_JOB_ID_LEN:
            if job_id.startswith("TLSMD"):
                if mysql.job_exists(job_id):
                    return job_id
    return None

def vet_struct_id(data, max_len):
    if isinstance(data, unicode):
        return False
    if len(data) > max_len:
        return False
    if not data.isalnum():
        return False
    return True

def start_job(pdbid):
    pdbid = pdbid.upper()

    if mysql.pdb_exists(pdbid) != None:
        return "PDB: %s was already run" % pdbid

    pdbfile_bin = webtlsmdd.fetch_pdb(pdbid)
    pdbfile = pdbfile_bin.data
    if len(pdbfile) == 0:
        return "FAILED: Could not download PDB %s from RCSB." % pdbid

    job_id = prepare_submission(pdbfile)

    try:
        mysql.set_pdb_db(pdbid)
    except:
        return "ERROR: Could not write to internal PDB DB"

    mysql.job_set_via_pdb(job_id, "1")
    mysql.job_set_jmol_view(job_id, "0")
    mysql.job_set_jmol_animate(job_id, "0")
    mysql.job_set_histogram(job_id, "0")
    mysql.job_set_private_job(job_id, "0")

    ip_addr = os.environ.get("REMOTE_ADDR", "Unknown")
    mysql.job_set_remote_addr(job_id, ip_addr)

    mysql.job_set_state(job_id, "queued")

    return "NOTE: Starting PDB %s with job_id %s" % (pdbid, job_id)

def prepare_submission(pdbfile):
    """class SubmitPDBPage
    """
    job_id = mysql.job_new()

    ## basic sanity checks
    ## If check_upload returns anything but a empty string, the server will
    ## inform the user of the problem and not proceed any further.
    ln = pdbfile.split("\n")
    r = check_upload(job_id, ln)
    if r != '':
        console.stdoutln("WARNING: %s" % str(r))
        sys.exit(0)

    result = webtlsmdd.set_structure_file(job_id, xmlrpclib.Binary(pdbfile))
    if result != "":
        console.stdoutln("ERROR: Failed to submit structure. %s. Please try again." % result)
        sys.exit(0)

    return job_id

def redirect_page(self, pdbid):
    return ""

def running_stddev(atomnum, restype, resnum, chain, tfactor):
    """Calculates a running standard deviation for the average B-factors
    of a given set of residues (controlled by the 'window' variable).
    """
    tmpfile = misc.generate_security_code()
    n = atm = res_tfac = 0
    avg_tfac = []
    res_id = []
    prevrestype = restype[0]
    prevresnum = resnum[0]
    prevchain = chain[0]
    ## Save B_{mean} per residue for each chain
    while n < len(tfactor):
        if( (prevresnum == resnum[n]) and (prevrestype == restype[n]) ):
            res_tfac = res_tfac + tfactor[n]
            atm = atm + 1
        else:
            avg_tfac.append(res_tfac/atm) # store previous guy
            res_id.append(resnum[n-1])    # store previous guy
            res_tfac = tfactor[n]
            atm = 1
            prevrestype = restype[n]
            prevresnum = resnum[n]
            if(prevchain != chain[n]):
                prevchain = chain[n]
        n = n + 1
    avg_tfac.append(res_tfac/atm)        # store last guy
    res_id.append(resnum[n-1])           # store last guy

    ## Save RMSD(B) +/-5 residues
    ## FIXME EAM
    ## Not correct, because it crosses chain boundaries and because the wrong 
    ## value is calculated (std of mean, rather than the std of the atoms)
    nbad = 0
    for s in range(5, len(avg_tfac)-5):
        stddev11 = numpy.std(avg_tfac[s-5:s+5])
        if stddev11 < conf.MIN_STDDEV_BFACT or stddev11 > conf.MAX_STDDEV_BFACT:
            nbad = nbad + 1

    return nbad, tmpfile

def check_upload(job_id, file):
    """Runs sanity checks on uploaded file
    """
    ## Checks if PDB contains valids aa/na residues
    ## PDB must have at least 30 ATOMs
    ## PDB can not have lowercase alt. res. numbers
    ## Check Standard deviation of temp. factors
    ## Check that not all occupancies are 0.00
    atom_num = []
    res_type = []
    res_num = []
    chain = []
    temp_factors = []
    bad_std = -1
    num_total = 0
    num_good = 0
    occupancy = 0.0
    ignore = 0
    line_num = 0
    for line in file:
        line_num += 1
        if line.startswith('HEADER'):
            header_id = re.sub(r"^HEADER.{56}(....)", '\\1', line).strip()

        elif line.startswith('EXPDTA    NMR'):
            return "NMR structure! Skipping: %s [%s]" % (job_id, header_id)

        elif re.match(r'^REMARK   2 RESOLUTION\. ([0-9\.]{1,}) ANGSTROMS.*', line):
            resolution = re.sub(r'^REMARK   2 RESOLUTION\. ([0-9\.]{1,}) ANGSTROMS.*', '\\1', line).strip()

        elif re.match('^ATOM.....................[0-9][a-z]', line):
            ## E.g., Don't allow "100b". Force it to be "100B"
            return "Lowercase alternate residue names: %s [%s]" % (job_id, header_id)

        elif line.startswith('ATOM') and (
            Library.library_is_standard_residue(line[17:20].strip())):
            num_total += 1
            if float(line[56:60].strip()) < 1.00:
                ## ignore occupancies < 1.00
                ignore += 1
                continue
            else:
                num_good += 1
                atom_num.append(int(line[7:11].strip()))
                res_type.append(line[17:20].strip())
                res_num.append(int(line[23:26].strip()))
                chain.append(line[21:22])
                occupancy += float(line[56:60].strip())
                temp_factors.append(float(line[60:65].strip()))
        else:
            continue

    if(len(atom_num) < 30):
        return "Not a PDB structure or has unrecognized residue names: %s [%s]" % (
            job_id, header_id)

    if(occupancy / num_good == 0.0):
        return "All occupancies are 0.0. TLSMD won't run on this structure: %s [%s]" % (
            job_id, header_id)

    bad_std, tmpfile = running_stddev(atom_num, res_type, res_num, chain, temp_factors)
    if bad_std > 0:
        ## If there are a string of "bad" B-factors, return a plot showing the
        ## "bad" regions and do not proceed any further in the analysis.
        return_string = "STDDEV %s > Bfact < %s for job_id: %s [%s]" % (
            conf.MAX_STDDEV_BFACT, conf.MIN_STDDEV_BFACT, job_id, header_id)
        return return_string

    return ''


def main():
    try:
        pdbid = sys.argv[1]
    except IndexError:
        sys.exit(1)

    r = start_job(pdbid.upper())
    console.stdoutln("%s" % r)

if __name__=="__main__":
    main()
    sys.exit(0)
