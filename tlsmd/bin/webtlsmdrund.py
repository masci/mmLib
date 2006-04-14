#!/home/tlsmd/local/bin/python
## TLS Motion Determination (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys
import time
import re
import fcntl
import subprocess
import socket
import xmlrpclib


from tlsmdlib import const, conf, email


def log_write(x):
    sys.stdout.write(x + "\n")
    sys.stdout.flush()

def log_job_start(jdict):
    tlsmd = jdict["tlsmd"]
    
    ln  = ""
    ln += "[%s]: " % (time.asctime(time.localtime(time.time())))
    ln += " ".join(tlsmd)
    log_write(ln)

def chain_size_string(jdict):
    if jdict.has_key("chains") == False:
        return "---"
    listx = []
    for cdict in jdict["chains"]:
        listx.append("%s:%d" % (cdict["chain_id"], cdict["length"]))
    return ";".join(listx)
    
def log_job_end(jdict):
    ln  = ""
    ln += "[%s]: " % (time.asctime(time.localtime(time.time())))
    ln += "Finished Job %s" % (jdict["job_id"])
    log_write(ln)
 
    ## write to a special log file
    if jdict.get("private_job", True):
        private_text = "private"
    else:
        private_text = "public"
    
    l = [time.asctime(time.localtime(time.time())),
         jdict.get("ip_addr", "000.000.000.000"),
	 jdict.get("email", "nobody@nowhere.com"),
	 private_text,
	 jdict.get("job_id", "EEK!!"),
	 jdict.get("structure_id", "----"),
         chain_size_string(jdict) ]

    try:
        open(conf.LOG_PATH, "a").write(" ".join(l) + "\n")
    except IOError:
        log_write("ERROR: cannot open logfile %s" % (conf.LOG_PATH))

def log_error(jdict, err):
    ln  = ""
    ln += "[%s]: " % (time.asctime(time.localtime(time.time())))
    ln += "ERROR: %s" % (err)
    log_write(ln)


def get_cdict(chains, chain_id):
    for cdict in chains:
        if chains["chain_id"] == chain_id:
            return cdict
    return None


def run_tlsmd(webtlsmdd, jdict):
    job_id = jdict["job_id"]
    tlsmd  = jdict["tlsmd"]

    ## write the tlsmd execution command out to a file
    open("tlsmdcmd.txt", "w").write(" ".join(tlsmd) + '\n')

    pobj = subprocess.Popen(tlsmd,
                            stdin = subprocess.PIPE,
                            stdout = subprocess.PIPE,
                            stderr = subprocess.STDOUT,
                            close_fds = True,
                            bufsize = 0)

    logfil = open("log.txt", "w")

    chain_time_dict = {}
    time_chain_id   = None
    time_begin      = None
    
    while True:
        ln = pobj.stdout.readline()
	if len(ln) == 0:
            break

        logfil.write(ln)
	logfil.flush()

    pobj.wait()
    logfil.close()


def run_job(webtlsmdd, jdict):
    job_id = jdict["job_id"]
    log_job_start(jdict)

    webtlsmdd.job_set_run_time_begin(job_id, time.time())

    old_dir = os.getcwd()

    ## change to the job directory, and run TLSMD
    job_dir = webtlsmdd.job_get_job_dir(job_id)
    try:
        os.chdir(job_dir)
    except os.error, err:
        log_error(jdict, str(err))
    else:
        run_tlsmd(webtlsmdd, jdict)

    os.chdir(old_dir)
    log_job_end(jdict)

    webtlsmdd.job_set_state(job_id, "completed")
    webtlsmdd.job_set_run_time_end(job_id, time.time())

    ## send email now that the job is complete
    send_mail(job_id)


def get_job(webtlsmdd):
    """Remove the top job from the queue file and return it.
    """
    try:
        job_id = webtlsmdd.get_next_queued_job_id()
    except socket.error:
        log_write("[ERROR] unable to connect to webtlsmdd.py")
        raise SystemExit

    if job_id == "":
        return None

    ## change state of the job and re-load to catch
    ## any updates which may have happened
    if webtlsmdd.job_set_state(job_id, "running") == False:
        return None

    jdict = webtlsmdd.job_get_dict(job_id)
    if jdict == None:
        return None

    tlsmd = [conf.TLSMD_PROGRAM_PATH, "-b", "-rANALYSIS" ]

    ## Job ID and webtlsmdd URL
    tlsmd.append("-j%s" % (job_id))
    tlsmd.append("-x%s" % (conf.WEBTLSMDD))

    ## override PDB ID
    tlsmd.append("-i%s" % (jdict["structure_id"]))

    ## plot style
    if jdict.get("plot_format") == "SVG":
        tlsmd.append("-s")

    ## select TLS model
    tls_model = jdict["tls_model"]
    if tls_model == "ISOT":
        tlsmd.append("-mISOT")
    else:
        tlsmd.append("-mANISO")

    ## select LSQ weighting
    #tlsmd.append("-w%s" % (jdict["weight"]))

    ## select chain IDs to analyze
    cids = []
    for cdict in jdict["chains"]:
        if cdict["selected"] == True:
            cids.append(cdict["chain_id"])
    tlsmd.append("-c%s" % (",".join(cids)))

    ## included atoms
    include_atoms = jdict["include_atoms"]
    tlsmd.append("-a%s" % (include_atoms))

    ## input PDB file
    tlsmd.append(jdict["pdb_filename"])

    jdict["tlsmd"] = tlsmd
    return jdict


MAIL_MESSAGE = """\
This is a automated message sent to you by the TLS Motion 
Determination (TLSMD) Server to inform you the analysis of the structure
you submitted is complete.  The link below will take you directly
to the completed analysis:

http://skuld.bmsc.washington.edu<ANALYSIS_URL>

If you experience any problems, please send us some email.

Regards,
Jay Painter <jpaint@u.washington.edu>
Ethan Merritt <merritt@u.washington.edu>

"""

def send_mail(job_id):
    if not os.path.isfile(conf.MSMTP):
        log_write("mail client not found: %s" % (conf.MSMTP))
        return
    
    webtlsmdd = xmlrpclib.ServerProxy(conf.WEBTLSMDD)

    jdict = webtlsmdd.job_get_dict(job_id)
    if jdict == False:
        log_write("job_id not found: %s" % (job_id))
        return

    address = jdict.get("email", "")
    if len(address) == 0:
        log_write("no email address")
        return

    analysis_url = jdict.get("analysis_url", "")
    if len(analysis_url)==0:
        log_write("no analysis_url")
        return

    ## send mail using msmtp
    email.SendEmail(
        address, 
        "Your TLSMD Job %s is Complete" % (job_id),
        MAIL_MESSAGE.replace("<ANALYSIS_URL>", analysis_url))
    
    log_write("sent mail to: %s" % (address))

def fetch_and_run_jobs_forever():
    log_write("starting webtlsmdrund.py  version %s" % (const.VERSION))
    log_write("using xmlrpc server webtlsmdd.py at URL.........: %s" % (conf.WEBTLSMDD))

    webtlsmdd = xmlrpclib.ServerProxy(conf.WEBTLSMDD)

    while True:
        jdict = get_job(webtlsmdd)
        if jdict == None:
            time.sleep(5.0)
            continue

        run_job(webtlsmdd, jdict)

def main():
    try:
        fetch_and_run_jobs_forever()
    except:
        email.SendTracebackEmail("webtlsmdrund.py exception")
        raise

if __name__=="__main__":
    if len(sys.argv) == 2:
        send_mail(sys.argv[1])
    else:
        main()

