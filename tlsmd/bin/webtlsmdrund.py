#!/usr/bin/python
# coding=UTF-8
## TLS Motion Determination (TLSMD)
## Copyright 2002-2009 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

## Python modules
import os
import sys
import time
import signal
import re
import fcntl
from subprocess import Popen, call, PIPE
import subprocess
import socket
import xmlrpclib
import shutil
import tarfile

## TLSMD
from tlsmdlib import const, conf, email, misc, mysql_support

## XXX: This is never used, 2009-05-29
## if there are no lines of output from the tlsmd process
## in TIMEOUT_SECS, then kill the process
TIMEOUT_SECS = 1 * (60 * 60)

def log_write(x):
    sys.stdout.write("[%s] " % misc.timestamp() + x + "\n")
    sys.stdout.flush()

def log_job_start(tlsmd):
    ln  = ""
    ln += " ".join(tlsmd)
    log_write(ln)

def log_job_died(job_id):
    ln  = ""
    ln += "Killed Job %s from CLI" % (job_id)
    log_write(ln)

def chain_size_string(chain_sizes):
    listx = []
    for c in chain_sizes.split(';'):
        chid, length, selected, type = misc.parse_chains(c)
        if selected == "1":
            listx.append("%s:%s" % (chid, length))
    return ";".join(listx)

def chain_type_string(chain_sizes):
    listx = []
    for c in chain_sizes.split(';'):
        chid, length, selected, type = misc.parse_chains(c)
        if selected == "1":
            listx.append("%s:%s" % (chid, type))
    return ";".join(listx)

def log_job_end(jdict):
    ln  = ""
    ln += "Finished Job %s" % (jdict["job_id"])
    log_write(ln)

    ## write to a special log file
    if jdict.get("private_job", True):
        private_text = "private"
    else:
        private_text = "public"

    submit_time = jdict.get('submit_time', 0.0)
    run_time_begin = jdict.get('run_time_begin', 0.0)
    run_time_end = jdict.get('run_time_end', 0.0)
    processing_time = timediff(run_time_begin, run_time_end)

    l = ["[Submit time: %s]"  % (timestring(submit_time)), ## 2009-04-21 17:31 PDT
         "[Start time: %s] " % (timestring(run_time_begin)),
         "[End time: %s] " % (timestring(run_time_end)),
         "[Processing time: %s] " % (processing_time),
         "[IP : %s] " % (jdict.get("ip_address", "000.000.000.000")),
         "[Email: %s] " % (jdict.get("email", "nobody@nowhere.com")),
         "[Privacy: %s] " % (private_text),
         "[Job ID: %s] " % (jdict.get("job_id", "EEK!!")),
         "[Structure ID: %s] " % (jdict.get("structure_id", "----")),
         "[Header ID: %s] " % (jdict.get("header_id", "----")),
         "[Chain sizes: %s] " % (chain_size_string(jdict["chain_sizes"])),
         "[Chain types: %s] " % (chain_type_string(jdict["chain_sizes"])),
         "[TLS Model: %s] " % (jdict.get('tls_model', 'None')),
         "[Weight: %s] " % (jdict.get('weight', 'None')),
         "[Atoms: %s] " % (jdict.get('include_atoms', 'None')),
         "[Nparts: %s] " % (jdict.get('nparts', 'None')),
         "[Max segs: %s] " % (jdict.get('chain_max_segs', 'None')),
         "[Resolution: %s] " % (jdict.get('resolution', 'None')),
         "[Initial residuals: %s] " % (jdict.get('initial_residuals', 'None')),
         "[Final residuals: %s] " % (jdict.get('final_residuals', 'None')),
         "[STDDEV Bfact: %s] " % (jdict.get('stddev_bfact', 'None')),
         "[Jmol animate: %s] " % (jdict.get('generate_jmol_animate', 'None')),
         "[Jmol viewer: %s] " % (jdict.get('generate_jmol_view', 'None')),
         "[Histogram: %s] " % (jdict.get('generate_histogram', 'None')),
         "[State: %s] " % (jdict.get('state', 'None'))]

    try:
        open(conf.LOG_FILE, "a").write(" ".join(l) + "\n")
    except IOError:
        log_write("ERROR: cannot open logfile %s" % (conf.LOG_FILE))

def timestring(raw_time):
    t = time.asctime(time.localtime(raw_time))
    return t

def itimestring(secs):
    tm_struct = time.localtime(secs)
    return time.strftime("%Y-%m-%d %H:%M:%S", tm_struct)

def timediff(begin, end):
    secs = int(end - begin)
    hours = secs / 3600
    secs = secs - (hours * 3600)
    min = secs / 60
    secs = secs - (min * 60)
    x = "%1d:%2d.%2d" % (hours, min, secs)
    return x.replace(" ", "0")

def log_error(err):
    ln  = ""
    ln += "ERROR: %s" % (err)
    log_write(ln)

def get_cdict(chains, chain_id):
    for cdict in chains:
        if chains["chain_id"] == chain_id:
            return cdict
    return None

def run_tlsmd(mysql, jdict):
    """main tlsmd fork/exec routine"""
    tlsmd = jdict["tlsmd"]

    ## write the tlsmd execution command out to a file
    open("tlsmdcmd.txt", "w").write(" ".join(tlsmd) + '\n')

    ## FORK SECTION: fork/execvp; 2008-01-22
    ## e.g., args="python tlsmd.py -b -rANALYSIS -jTLSMD9370_iNMLMncN \
    ##               -i1FIN -mISOT -cA -aALL struct.pdb"
    args = ["python"] + tlsmd
    pid = os.fork()
    if not pid:
        ## We are child

        redirect = os.open("log.txt", os.O_WRONLY | os.O_CREAT)
        os.dup2(redirect, 1) ## standard output (1)
        os.dup2(redirect, 2) ## standard error  (2)
        os.close(redirect)

        save_pid = os.getpid()  ## capture child pid
        mysql.job_set_pid(jdict["job_id"], save_pid)

        ## Set up the environment
        pwd = {}
        pwd["PWD"] = os.path.join(conf.TLSMD_WORK_DIR, jdict["job_id"])
        os.environ.update(pwd)

        os.execvp("python", args)
        ## NOTE: We never come back from 'execvp()'

    time.sleep(5.0)
    return

def run_job(mysql, jdict):
    job_id = jdict["job_id"]
    log_job_start(jdict["tlsmd"])

    mysql.job_set_run_time_begin(job_id, time.time())

    old_dir = os.getcwd()

    ## change to the job directory, and run TLSMD
    job_dir = os.path.join(conf.TLSMD_WORK_DIR, jdict["job_id"])
    try:
        os.chdir(job_dir)
    except os.error, err:
        log_error(str(err))
    else:
        run_tlsmd(mysql, jdict)

    return

def check_logfile_for_errors(file):
    """Searches through the log.txt file for warnings and errors
    """
    warnings = False
    errors = False
    completed = False
    infil = open(file,'r').readlines()
    for line in infil:
        if re.match(r'^\s*Warning:', line):
            warnings = True
        elif re.match(r'^.*[Ee][Rr][Rr][Oo][Rr]', line):
            errors = True
        elif line.startswith('completed'):
            completed = True
        else:
            continue

    if not completed:
        return "died"
    if errors:
        return "errors"
    if warnings:
        return "warnings"
    return "success"

def cleanup_job(mysql, jdict):
    """Cleanup job directory upon completion and email user
    """
    old_dir = os.getcwd()
    job_id = jdict["job_id"]
    job_dir = os.path.join(conf.TLSMD_WORK_DIR, job_id)
    pdb_dir = os.path.join(conf.WEBTLSMDD_PDB_DIR, jdict["structure_id"])
    try:
        os.chdir(job_dir)
    except os.error, err:
        log_error(str(err))
        mysql.job_set_state(job_id, "lost_directory")
        return

    ## Check if user submitted job via PDB code (i.e., from pdb.org)
    if int(jdict["via_pdb"]) == 1:
        if os.path.exists(pdb_dir):
            shutil.rmtree(pdb_dir)
        try:
            shutil.copytree(job_dir, pdb_dir)
        except OSError:
            raise

    ## create tarball
    tar = tarfile.open("%s.tar.gz" % job_id, "w:gz")
    tar.add("ANALYSIS")
    tar.close()

    os.chdir(old_dir)
    ## check 'log.txt' for warnings
    if check_logfile_for_errors(job_dir + "/log.txt") == "warnings":
        mysql.job_set_state(job_id, "warnings")

    mysql.job_set_run_time_end(job_id, time.time())
    log_job_end(mysql.job_get_dict(job_id))

    ## send email now that the job is complete
    send_mail(jdict)

def get_job(mysql):
    """Remove the top job from the queue file and return it.
    """
    ## NOTE: This is the very first step in the entire webtlsmdrund.py process
    try:
        job_id = mysql.get_next_queued_job_id()
    except socket.error:
        log_write("[ERROR] unable to connect to MySQL")
        raise SystemExit

    if job_id == "":
        return None

    ## change state of the job and re-load to catch any updates which may
    ## have occurred
    if mysql.job_set_state(job_id, "running") == False:
        return None
    mysql.job_set_state(job_id, "running")

    jdict = mysql.job_get_dict(job_id)
    if jdict == None:
        return None

    ## Construct the run command for this job
    tlsmd = [conf.TLSMD_PROGRAM_PATH, "-b", "-rANALYSIS" ]

    ## Job ID
    tlsmd.append("-j%s" % (job_id))

    ## override PDB ID
    tlsmd.append("-i%s" % (jdict["structure_id"]))

    ## Jmol/Histogram features
    if jdict["generate_jmol_view"] == True:
        tlsmd.append("--generate-jmol-viewer")
    if jdict["generate_jmol_animate"] == True:
        tlsmd.append("--generate-jmol-animate")
    if jdict["generate_histogram"] == True:
        tlsmd.append("--generate-histogram")

    ## plot style
    if jdict["plot_format"] == "SVG":
        tlsmd.append("-s")

    ## select TLS model
    tls_model = jdict["tls_model"]
    if tls_model == "ISOT":
        tlsmd.append("-mISOT")
    else:
        tlsmd.append("-mANISO")

    ## set maximum number of segments
    tlsmd.append("-n%s" % (jdict["nparts"]))

    ## select chain IDs to analyze
    cids = []
    chains = jdict["chain_sizes"].rstrip(";")
    for c in chains.split(';'):
        chid, length, selected, type = misc.parse_chains(c)
        if selected == "1":
            cids.append(chid)
    tlsmd.append("-c%s" % (",".join(cids)))

    ## included atoms
    tlsmd.append("-a%s" % (jdict["include_atoms"]))

    ## input PDB file
    tlsmd.append(conf.PDB_FILENAME)

    jdict["tlsmd"] = tlsmd
    return jdict


MAIL_MESSAGE = """\
This is an automated message sent to you by the TLS Motion 
Determination (TLSMD) Server to inform you the analysis of the structure
you submitted is complete.  The link below will take you directly
to the completed analysis:

<BASE_URL><ANALYSIS_URL>

having the following user comments:
"<USER_COMMENT>"

If you experience any problems, please send us some email.

Regards,
Christoph Champ <champc@u.washington.edu>
Ethan Merritt <merritt@u.washington.edu>

"""

def send_mail(jdict):
    if not os.path.isfile(conf.MAIL):
        log_write("ERROR: mail client not found: %s" % (conf.MAIL))
        return

    job_id = jdict["job_id"]
    if jdict == False:
        log_write("WARNING: Trying to email but job_id not found: %s" % (job_id))
        return

    address = jdict.get("email", "")
    if len(address) == 0:
        return

    job_url = "%s/%s" % (conf.TLSMD_WORK_URL, job_id)
    analysis_url = "%s/ANALYSIS/index.html" % (job_url)
    if len(analysis_url) == 0:
        log_write("WARNING: no analysis_url: %s" % job_id)
        return

    user_comment = jdict.get("user_comment", "")
    if len(user_comment) == 0:
        user_comment = "no comment"

    ## send mail using msmtp
    mail_message = MAIL_MESSAGE
    mail_message = mail_message.replace("<BASE_URL>", conf.BASE_PUBLIC_URL)
    mail_message = mail_message.replace("<ANALYSIS_URL>", analysis_url)
    mail_message = mail_message.replace("<USER_COMMENT>", user_comment)
    email.SendEmail(
        address,
        "Your TLSMD Job %s is Complete" % (job_id),
        mail_message)

def is_pid_running(full_cmd):
    """Checks if a given PID is still running. Returns "True" if job is still
    running.
    """
    try:
        p = Popen(full_cmd, shell=True, stdout=PIPE)
        output = p.communicate()[0]
        if re.match('.*defunct.*', output):
            return False ## Sometimes python goes <defunct> on a PID. Job no longer running.
        elif output:
            return True ## If anything besides "defunct" is returned from ps, job is still running
    except Exception, e:
        print >>sys.stderr, "Execution failed:", e
        return False

    return False

def check_for_pid(pid):
    """Function for checking if process is running. Returns "True" if job is
    still running.
    """
    try:
        cmd = "ps -p %s --no-heading" % pid
        res = is_pid_running(cmd)

        if res:
           return True ## PID still running
        else:
           return False ## PID _not_ running

    except Exception, e:
        ## Something wrong happened with the ps command
        print "ERROR: Can't run is_pid_running(): %s" % e

    return False

def job_completed(mysql, jdict):
    """Checks if "running" job is finished and sets state according to log.txt
    messages. If there is a match, return False (job not completed).
    Otherwise return True (job completed)
    """
    job_id = jdict["job_id"]
    job_dir = os.path.join(conf.TLSMD_WORK_DIR, job_id)
    try:
        os.chdir(job_dir)
    except:
        mysql.job_set_state(job_id, "lost_directory")
        return True

    pid = mysql.job_get_pid(job_id)
    if pid == None:
        ## job PID wasn't stored in database; something must have gone wrong
        mysql.job_set_state(job_id, "syserror")
        return True
    else:
        pid = int(pid)

    if check_for_pid(pid):
        ## Job is still running
        return False

    logfile = os.path.join(job_dir, "log.txt")
    mysql.job_set_state(job_id, check_logfile_for_errors(logfile))
    return True

def fetch_and_run_jobs_forever():
    log_write("starting webtlsmdrund.py  version %s" % (const.VERSION))
    log_write("using xmlrpc server webtlsmdd.py at URL.........: %s" % (conf.WEBTLSMDD))

    mysql = mysql_support.MySQLConnect()

    running_list = [] ## New array to hold list of jobs currently running
    while True:
        ## First clear out any jobs in run queue that have finished
        ## Get information on jobs in queue
        job_list = mysql.job_list()
        for myjdict in job_list:
            if myjdict.get("state") == None:
                ## This is a "Partially Submitted" job, so go to next in list
                continue
            ## Comment out the following three lines to refresh the states in
            ## the database
            if myjdict.get("state") != "running":
                ## Check for any other states besides None and "running"
                continue
            if job_completed(mysql, myjdict):
                for n in range(len(running_list)):
                    if myjdict["job_id"] == running_list[n]:
                        if mysql.job_get_state(myjdict["job_id"]) != "died":
                            cleanup_job(mysql, myjdict)

                        ## remove completed job_id from array
                        running_list.pop(n)
                        break

        ## Check whether there is a slot free to start a new run
        if len(running_list) >= conf.MAX_PARALLEL_JOBS:
            time.sleep(5.0)
            continue

        ## If there is a run slot free, get next job in queue
        jdict = get_job(mysql)
        if jdict == None:
            time.sleep(5.0)
            continue

        ## Start the new job
        running_list.append(jdict["job_id"])
        run_job(mysql, jdict)

def main():
    try:
        fetch_and_run_jobs_forever()
    except KeyboardInterrupt:
        sys.exit(1)
    except:
        email.SendTracebackEmail("webtlsmdrund.py exception")
        raise

if __name__=="__main__":
    if len(sys.argv) == 2:
        send_mail(sys.argv[1])
    else:
        main()

