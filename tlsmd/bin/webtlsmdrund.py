#!/home/jpaint/local/bin/python
## TLS Minimized Domains (TLSMD)
## Copyright 2002 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys
import time
import string
import re
import fcntl
import popen2
import xmlrpclib


## CONFIGURATION
VERSION   = "0.0.1"
WEBTLSMDD = "http://localhost:10100"

TLSMD_CMD = [
    "%s/bin/tlsmd.py" % (os.environ["TLSMD_ROOT"]),
    "-rANALYSIS" ]

if os.environ.has_key("TLSMD_GRID_FILE"):
    TLSMD_CMD.append("-f%s" % (os.environ["TLSMD_GRID_FILE"]))
    

## regular expression for parsing the output of TLSMD while
## it is running -- this is used for updating the calculation
## progress information
RE_PROCESS_CHAIN = re.compile(
    "process_chain\(chain_id=(\w+),\s*frag_id=\{(\w+)\.\.(\w+)\},.*$")

def log_write(x):
    sys.stdout.write(x+"\n")
    sys.stdout.flush()

def log_job_start(jdict):
    tlsmd = jdict["tlsmd"]
    
    ln  = ""
    ln += "[%s]: " % (time.asctime(time.localtime(time.time())))
    ln += string.join(tlsmd, " ")
    log_write(ln)
    
def log_job_end(jdict):
    ln  = ""
    ln += "[%s]: " % (time.asctime(time.localtime(time.time())))
    ln += "Finished Job %s" % (jdict["job_id"])
    log_write(ln)

def log_error(jdict, err):
    ln  = ""
    ln += "[%s]: " % (time.asctime(time.localtime(time.time())))
    ln += "ERROR: %s" % (err)
    log_write(ln)

def get_cdict(chains, chain_id):
    for cdict in chains:
        if chains["chain_id"]==chain_id:
            return cdict
    return None

def run_tlsmd(webtlsmdd, jdict):
    job_id = jdict["job_id"]
    tlsmd  = jdict["tlsmd"]

    stdout, stdin = popen2.popen4(tlsmd)
    stdin.close()

    logfil = open("log.txt", "w")

    chain_time_dict = {}
    time_chain_id   = None
    time_begin      = None
    
    segments = 0

    while True:
        ln = stdout.readline()
	if len(ln)==0:
            break

        ## for recording processing time...
        if ln.startswith("BEGIN TIMEING CHAIN ID "):
            if time_chain_id==None:
                time_chain_id = ln[23]
                time_begin = time.time()

        if ln.startswith("END TIMEING CHAIN ID "):
            chain_id = ln[21]
            if time_chain_id==chain_id:
                time_taken = time.time() - time_begin

                if chain_time_dict.has_key(chain_id):
                    chain_time_dict[chain_id] += time_taken
                else:
                    chain_time_dict[chain_id] = time_taken

            time_chain_id = None
            time_begin    = None

        m = RE_PROCESS_CHAIN.match(ln)
        if m!=None:
            segments += 1
            chain_id, frag_id1, frag_id2 = m.groups()

            ## only update the run-time database periodically
            if segments%100==0:
                chain_id, frag_id1, frag_id2 = m.groups()
                webtlsmdd.job_data_set(job_id, "run_chain_id", chain_id)
                webtlsmdd.job_data_set(job_id, "run_frag_id1", frag_id1)
                webtlsmdd.job_data_set(job_id, "run_frag_id2", frag_id2)
                
        logfil.write(ln)
	logfil.flush()

    stdout.close()
    logfil.close()

    ## record the time used for each chain
    chains = webtlsmdd.job_data_get(job_id, "chains")
    for cdict in chains:
        chain_id = cdict["chain_id"]
        if chain_time_dict.has_key(chain_id):
            cdict["processing_time"] = chain_time_dict[chain_id]
        webtlsmdd.job_data_set(job_id, "chains", chains)

def run_job(webtlsmdd, jdict):
    job_id = jdict["job_id"]
    log_job_start(jdict)

    webtlsmdd.job_data_set(job_id, "run_start_time", time.time())

    old_dir = os.getcwd()

    ## change to the job directory, and run TLSMD
    try:
        os.chdir(jdict["job_dir"])
    except os.error, err:
        log_error(jdict, str(err))
    else:
        run_tlsmd(webtlsmdd, jdict)

    os.chdir(old_dir)
    log_job_end(jdict)

    job_id = jdict["job_id"]
    webtlsmdd.job_data_set(job_id, "state", "completed")
    webtlsmdd.job_data_set(job_id, "run_end_time", time.time())

def get_job(webtlsmdd):
    """Remove the top job from the queue file and return it.
    """
    ## retrieve the top queued job
    jdict = False
    i = 0
    while True:
        jdict = webtlsmdd.job_get_dict_index(i)
        if jdict==False:
            break
        if jdict.get("state")=="queued":
            break
        i += 1
    if jdict==False:
        return None

    job_id = jdict["job_id"]

    ## change state of the job and re-load to catch
    ## any updates which may have happened
    if webtlsmdd.job_data_set(job_id, "state", "running")==False:
        return None
    jdict = webtlsmdd.job_get_dict(job_id)

    tlsmd = TLSMD_CMD[:]

    ## Job ID and webtlsmdd URL
    tlsmd.append("-j%s" % (job_id))
    tlsmd.append("-x%s" % (WEBTLSMDD))

    ## select TLS model
    tls_model = jdict["tls_model"]
    if tls_model=="ISOT":
        tlsmd.append("-mHYBRID")
    else:
        tlsmd.append("-mANISO")

    ## select LSQ weighting
    tlsmd.append("-w%s" % (jdict["weight"]))

    ## select chain IDs to analyize
    cids = []
    for cdict in jdict["chains"]:
        if cdict["selected"]==True:
            cids.append(cdict["chain_id"])
    tlsmd.append("-c%s" % (string.join(cids,",")))

    ## included atoms
    include_atoms = jdict["include_atoms"]
    tlsmd.append("-a%s" % (include_atoms))

    ## input PDB file
    tlsmd.append(jdict["pdb_filename"])

    jdict["tlsmd"] = tlsmd
    return jdict
    
def main():
    log_write("Starting WebTLSMDRunD v%s" % (VERSION))
    log_write("using xmlrpc server webtlsmdd.py at %s" % (WEBTLSMDD))

    webtlsmdd = xmlrpclib.ServerProxy(WEBTLSMDD, allow_none=1)

    while True:
        jdict = get_job(webtlsmdd)
        if jdict==None:
            time.sleep(1.0)
            continue

        run_job(webtlsmdd, jdict)

if __name__=="__main__":
    main()

