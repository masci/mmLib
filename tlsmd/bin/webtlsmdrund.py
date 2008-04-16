#!/usr/bin/python
# coding=UTF-8
## TLS Motion Determination (TLSMD)
## Copyright 2002-2008 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys
import time
import signal
import re
import fcntl
from subprocess import Popen, call, PIPE  # Christoph Champ, 2008-03-17
import subprocess
import socket
import xmlrpclib
import shutil
import tarfile	# Christoph Champ, 2007-12-03
import datetime # Christoph Champ, 2008-01-29

from tlsmdlib import const, conf, email

## if there are no lines of output from the tlsmd process
## in TIMEOUT_SECS, then kill the process
TIMEOUT_SECS = 1 * (60 * 60)

def log_write(x):
    sys.stdout.write(x + "\n")
    sys.stdout.flush()

def log_job_start(jdict):
    tlsmd = jdict["tlsmd"]
    
    ln  = ""
    ## Changed date format to international: YYYY-MM-DD HH:MM:SS; Christoph Champ, 2008-01-29
    ln += "[%s]: " % (datetime.datetime.fromtimestamp(time.time()).isoformat(' ')[:-7])
    ln += " ".join(tlsmd)
    log_write(ln)

def chain_size_string(jdict):
    if jdict.has_key("chains") == False:
        return "---"
    listx = []
    for cdict in jdict["chains"]:
        listx.append("%s:%d" % (cdict["chain_id"], cdict["length"]))
    return ";".join(listx)
    
def log_job_died(job_id):
    ln  = ""
    ln += "[%s]: " % (datetime.datetime.fromtimestamp(time.time()).isoformat(' ')[:-7])
    ln += "Killed Job %s from CLI" % (job_id)
    log_write(ln)

def log_job_end(jdict):
    ln  = ""
    ln += "[%s]: " % (datetime.datetime.fromtimestamp(time.time()).isoformat(' ')[:-7])
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
    l = ["[Submit time: %s]"  % (timestring(submit_time)),
         "[Start time: %s] " % (timestring(run_time_begin)),
         "[End time: %s] " % (timestring(run_time_end)),
         "[Processing time: %s] " % (processing_time),
         "[IP : %s] " % (jdict.get("ip_addr", "000.000.000.000")),
         "[Email: %s] " % (jdict.get("email", "nobody@nowhere.com")),
         "[Privacy: %s] " % (private_text),
         "[Job ID: %s] " % (jdict.get("job_id", "EEK!!")),
         "[Structure ID: %s] " % (jdict.get("structure_id", "----")),
         "[Chain sizes: %s] " % (chain_size_string(jdict)),
         "[TLS Model: %s] " % (jdict.get('tls_model', 'None')),
         "[Weight: %s] " % (jdict.get('weight', 'None')),
         "[Atoms: %s] " % (jdict.get('include_atoms', 'None')),
         "[State: %s] " % (jdict.get('state', 'None'))]

    try:
        open(conf.LOG_PATH, "a").write(" ".join(l) + "\n")
    except IOError:
        log_write("ERROR: cannot open logfile %s" % (conf.LOG_PATH))

def timestring(raw_time):
    t = time.asctime(time.localtime(raw_time))
    return t

def timediff(begin, end):
    secs = int(end - begin)
    hours = secs / 3600
    secs = secs - (hours * 3600)
    min = secs / 60
    secs = secs - (min * 60)
    x = "%1d:%2d.%2d" % (hours, min, secs)
    return x.replace(" ", "0")
    
    
def log_error(jdict, err):
    ln  = ""
    ln += "[%s]: " % (datetime.datetime.fromtimestamp(time.time()).isoformat(' ')[:-7])
    ln += "ERROR: %s" % (err)
    log_write(ln)


def get_cdict(chains, chain_id):
    for cdict in chains:
        if chains["chain_id"] == chain_id:
            return cdict
    return None

def run_tlsmd(webtlsmdd, jdict):
    """main tlsmd fork/exec routine"""
    tlsmd  = jdict["tlsmd"]

    ## write the tlsmd execution command out to a file
    open("tlsmdcmd.txt", "w").write(" ".join(tlsmd) + '\n')

    ### FORK SECTION: fork/execvp; 2008-01-22
    ## e.g., args="python tlsmd.py -b -rANALYSIS -jTLSMD9370_iNMLMncN -xhttp://localhost:10100 -i1FIN -mISOT -cA -aALL struct.pdb"
    args=["python"] + tlsmd
    pid=os.fork()
    if not pid:
	## We are child

	redirect = os.open("log.txt", os.O_WRONLY | os.O_CREAT)
	os.dup2(redirect, 1)                        # standard output (1)
	os.dup2(redirect, 2)                        # standard error (2)
	os.close(redirect)

	## Capture child pid. Christoph Champ, 2008-02-03
        save_pid = os.getpid()
	## Switched to using database field instead of file. Christoph Champ, 2008-03-14
	webtlsmdd.job_set_pid(jdict["job_id"],save_pid)

	## Set up the environment
	## TODO Find out which keys to clear. Christoph Champ, 2008-02-12
	## os.environ.clear()  ## probably too drastic, but some cleaning might be good
	pwd = {}
	pwd["PWD"] = jdict["job_dir"]
	os.environ.update(pwd)
	#print os.environ.values
	
	os.execvp("python", args)
	## WE NEVER COME BACK FROM execvp

    time.sleep(5.0)
    return

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
 
    return

def check_logfile_for_errors(file):
    """Searches through the log.txt file for warnings and errors"""
    ## Started by Christoph Champ, 2008-03-02
    
    warnings = False
    errors = False
    completed = False
    infil=open(file,'r').readlines()
    for line in infil:
	## Switched to re.match() for regex capability. Christoph Champ, 2008-03-11
	if re.match(r'^\s*Warning:',line):
	    warnings = True
	#elif re.match(r'^\s*[Ee][Rr][Rr][Oo][Rr]',line):
	elif re.match(r'^.*[Ee][Rr][Rr][Oo][Rr]',line):
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

def cleanup_job(webtlsmdd, jdict):
    """Cleanup job directory upon completion and email user
    """

    ## Force it to use correct path. Christoph Champ, 2008-02-01
    old_dir=os.getcwd()
    job_dir=conf.TLSMD_WORK_DIR+"/"+jdict["job_id"]
    try:
        os.chdir(job_dir)
    except os.error, err:
        log_error(jdict, str(err))
	webtlsmdd.job_set_state(jdict["job_id"], "lost_directory")
	return

    ### Check if user submitted job via PDB code (i.e., from pdb.org). Christoph Champ, 2008-02-12
    if jdict.get("via_pdb", False) and "pdb_dir" in jdict and len(jdict["pdb_dir"]) != 0:
        #print "%s: Archiving job..." % (time.ctime())
        pdb_dir = jdict['pdb_dir']
    
        if os.path.exists(pdb_dir):
            shutil.rmtree(pdb_dir)
        try:
            shutil.copytree(job_dir, pdb_dir)
        except OSError:
            raise

    ## create tarball. Christoph Champ, 2007-12-03
    tar=tarfile.open("%s.tar.gz"%jdict["job_id"], "w:gz")
    tar.add("ANALYSIS")
    tar.close()

    os.chdir(old_dir)
    ## check 'log.txt' for warnings. Christoph Champ, 2008-03-02
    if check_logfile_for_errors(job_dir+"/log.txt")=="warnings":
       webtlsmdd.job_set_state(jdict["job_id"], "warnings")
    webtlsmdd.job_set_run_time_end(jdict["job_id"], time.time())
    log_job_end(webtlsmdd.job_get_dict(jdict["job_id"]))

    ## send email now that the job is complete
    send_mail(jdict["job_id"])

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

    ## Construct the run command for this job
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
This is an automated message sent to you by the TLS Motion 
Determination (TLSMD) Server to inform you the analysis of the structure
you submitted is complete.  The link below will take you directly
to the completed analysis:

http://verdandi.bmsc.washington.edu<ANALYSIS_URL>

having the following user comments:
<USER_COMMENT>

If you experience any problems, please send us some email.

Regards,
Christoph Champ <champc@u.washington.edu>
Ethan Merritt <merritt@u.washington.edu>

"""

def send_mail(job_id):
    if not os.path.isfile(conf.MAIL):
        log_write("ERROR: mail client not found: %s" % (conf.MAIL))
        return

    webtlsmdd = xmlrpclib.ServerProxy(conf.WEBTLSMDD)

    jdict = webtlsmdd.job_get_dict(job_id)
    if jdict == False:
        log_write("WARNING: job_id not found: %s" % (job_id))
        return

    address = jdict.get("email", "")
    if len(address) == 0:
        log_write("NOTE: no email address")
        return

    analysis_url = jdict.get("analysis_url", "")
    if len(analysis_url)==0:
        log_write("WARNING: no analysis_url")
        return

    ## user_comment added. Christoph Champ, 2007-12-18
    user_comment = jdict.get("user_comment", "")
    if len(user_comment)==0:
        log_write("NOTE: no user_comment")
        #return # this was stopping "via_pdb" from emailing the user. Christoph Champ, 2008-02-12

    ## send mail using msmtp
    ## user_comment added. Christoph Champ, 2007-12-18
    mail_message = MAIL_MESSAGE
    mail_message = mail_message.replace("<ANALYSIS_URL>", analysis_url)
    mail_message = mail_message.replace("<USER_COMMENT>", user_comment)
    email.SendEmail(
        address,
        "Your TLSMD Job %s is Complete" % (job_id),
        mail_message)

    log_write("NOTE: sent mail to: %s" % (address))

def is_pid_running(full_cmd):
    ## Check if a given PID is still running. Christoph Champ, 2008-03-17
    try:
        p = Popen(full_cmd, shell=True, stdout=PIPE)
        output = p.communicate()[0]
        if re.match('.*defunct.*',output):
            return False ## Sometimes python goes <defunct> on a PID. Job no longer running.
        elif output:
            return True ## If anything besides "defunct" is returned from ps, job is still running
    except Exception, e:
        print >>sys.stderr, "Execution failed:", e
        return False

    return False

def check_for_pid(webtlsmdd,jdict):
    """Function for checking if process is running.
    """
    try:
        tmp_pid=webtlsmdd.job_get_pid(jdict["job_id"])
        pid=int(tmp_pid)
        cmd = "ps -p %s --no-heading" % pid
        res = is_pid_running(cmd)

        if res:
           return True ## PID still running
        else:
           return False ## PID _not_ running

    except Exception, e:
        # Something wrong happened with the ps command
        print "ERROR: Can't run is_pid_running(): %s"%e

    return False

def job_completed(webtlsmdd,jdict):
    ## Checks if "running" job is finished and sets state according to log.txt messages. Christoph Champ, 2008-03-21
    ## If there is a match, return False (job not completed). Otherwise return True (job completed); 2008-02-01
    try:
	os.chdir(conf.TLSMD_WORK_DIR+"/"+jdict["job_id"])
    except:
	webtlsmdd.job_set_state(jdict["job_id"], "lost_directory") ## There must be zero spaces in "state". Christoph Champ, 2008-02-07
	return True

    if check_for_pid(webtlsmdd,jdict):
	return False
    
    webtlsmdd.job_set_state(jdict["job_id"], check_logfile_for_errors(jdict["job_dir"]+"/log.txt"))
    return True

def fetch_and_run_jobs_forever():
    log_write("starting webtlsmdrund.py  version %s" % (const.VERSION))
    log_write("using xmlrpc server webtlsmdd.py at URL.........: %s" % (conf.WEBTLSMDD))

    webtlsmdd = xmlrpclib.ServerProxy(conf.WEBTLSMDD)

    running_list=[] # New array to hold list of jobs currently running
    while True:
        ### First clear out any jobs in run queue that have finished
        ### Get information on jobs in queue; 2008-01-22
        job_list = webtlsmdd.job_list()
        for myjdict in job_list:
	    if myjdict.get("state") == None:
		## This is a "Partially Submitted" job, so go to next in list
		continue
	    ## Comment out the following three lines to refresh the states in the database
	    if (myjdict.get("state") != "running"):
		## Check for any other states besides None and "running"
		continue
	    if job_completed(webtlsmdd,myjdict):
		for n in range(len(running_list)):
		    if myjdict["job_id"]==running_list[n]:
			if webtlsmdd.job_get_state(myjdict["job_id"]) != "died":
			   cleanup_job(webtlsmdd,myjdict) # Note: passing full jdict (not just job_id). Christoph Champ, 2008-02-12
			running_list.pop(n)  # remove completed job_id from array
			break

	### EAM 7-Feb-2008 revamp the logic limiting us to 4 run slots

	### Check whether there is a slot free to start a new run
	if len(running_list) >= 4:
	    time.sleep(5.0)
	    continue

	### If there is a run slot free, get next job in queue
	jdict = get_job(webtlsmdd)
	if jdict == None:
	    time.sleep(5.0)
	    continue

	### Start the new job
	running_list.append(jdict["job_id"])
	run_job(webtlsmdd,jdict)

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

