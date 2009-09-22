#!/usr/bin/python
# coding=UTF-8 
## TLS Motion Determination (TLSMD)
## Copyright 2002-2009 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

## Python
import os
import sys
import time
import string
import random
import traceback
## NOTE: The order of these signals is important!
from signal import SIG_IGN ## Needed for daemon_main()
from signal import SIGUSR1 ## Needed for SignalJob()
from signal import SIGHUP  ## Needed for KillJob()
import signal
import cPickle
import bsddb
import socket
import xmlrpclib
import SocketServer
import SimpleXMLRPCServer
import urllib
import gzip       ## for fetching PDBs from pdb.org
import StringIO
import subprocess ## for render of 'struct.png'
import re         ## for "raw grey"-backbone rendering. 2008-12-03

## pymmlib
from mmLib import FileIO

## TLSMD
from tlsmdlib import conf, const, tls_calcs, email, misc, mysql_support

mysql = mysql_support.MySQLConnect()

def fatal(text):
    sys.stderr.write("[FATAL ERROR] %s\n" % (text))
    raise SystemExit


def parse_molauto(infile, outfile):
    """Parses the molauto output to force each chain to have its own unique
       colour. 2008-12-03
    """
    file = open(outfile, "w")
    for line in open(infile).readlines():
        file.write("%s" % line)
        if(re.match(r'^  set segments', line)):
            file.write("  set segments 10;\n")
            file.write("  set planecolour hsb 0.6667 1 1;")
        elif(re.match(r'^  set planecolour', line)):
            colour = line
        elif(re.match(r'^  .* from ', line)):
            chain1 = chain2 = line
            chain1 = re.sub(r'^  .* from ([A-Z])[0-9]{1,} to ([A-Z])[0-9].*$', '\\1', line).strip()
            chain2 = re.sub(r'^  .* from ([A-Z])[0-9]{1,} to ([A-Z])[0-9].*$', '\\2', line).strip()
            file.write("%s" % line)
            if(chain1 != chain2):
                file.write("%s" % colour)
        else:
            file.write("%s" % line)
    file.close()
    return

def render_struct(job_dir):
    """Generate struct.png via molauto/parse_molauto/molscript/render
    """
    ## cmd: molauto smallAB.pdb|parse_molauto.pl|molscript -r |render -bg white -size 200x200 -png mymol.png
    cmdlist = ["%s %s/struct.pdb | %s | %s -r | %s -bg white -size %s -png %s/struct.png 1>&2" % (
              conf.MOLAUTO, job_dir, conf.PARSE_MOLAUTO_SCRIPT,
              conf.MOLSCRIPT, conf.RENDER,
              conf.RENDER_SIZE, job_dir)]
    proc = subprocess.Popen(cmdlist,
                            shell = True,
                            stdin = subprocess.PIPE,
                            stdout = subprocess.PIPE,
                            stderr = subprocess.PIPE,
                            close_fds = True,
                            bufsize = 32768)

def extract_raw_backbone(infile, outfile):
    """Extracts all backbone atoms (amino acid or nucleic acid) 
       for use in tlsanim2r3d
    """
    file = open(outfile, "w")
    for line in open(infile).readlines():
        if(re.match(r'^ATOM........ CA .*', line) or \
           re.match(r'^ATOM........ P  .*', line) or \
           re.match(r'^ATOM........ C[543]\'.*', line) or \
           re.match(r'^ATOM........ O[53]\'.*', line)):
            chain_id = line[21:22]
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())

            file.write("1 0 %s 0 0 %.3f %.3f %.3f\n" % (chain_id, x, y, z))
    file.close()
    return

def generate_raw_grey_struct(job_dir):
    """Generate 'raw' input for tlsanim2r3d, but only for the non-animated
       sections, for _all_ chains
    """
    ## cmd: ./extract_raw_chains.pl <smallAB.pdb |./tlsanim2r3d - >ANALYSIS/struct.r3d
    ## cmd: ./tlsanim2r3d < struct.raw >ANALYSIS/struct.r3d

    extract_raw_backbone("%s/struct.pdb" % job_dir, "%s/struct.raw" % job_dir)

    cmdlist = ["%s < %s/struct.raw > %s/struct.r3d 2>/dev/null" % (
              conf.TLSANIM2R3D, job_dir, job_dir)]
    proc = subprocess.Popen(cmdlist,
                            shell = True,
                            stdin = subprocess.PIPE,
                            stdout = subprocess.PIPE,
                            stderr = subprocess.PIPE,
                            close_fds = True,
                            bufsize = 32768)

def generate_bases_r3d(job_dir, chain_id):
    """Generate 'raw' input for tlsanim2r3d, but only for the non-animated
       sections, for _all_ chains
    """
    ## cmd: grep '^ATOM.................B.*' | rings3d -bases >>bases.r3d
    cmdlist = ["grep '^ATOM.................%s.*' %s/struct.pdb | %s -bases >>%s/bases.r3d 2>/dev/null" % (
              chain_id, job_dir, conf.RINGS3D, job_dir)]
    proc = subprocess.Popen(cmdlist,
                            shell = True,
                            stdin = subprocess.PIPE,
                            stdout = subprocess.PIPE,
                            stderr = subprocess.PIPE,
                            close_fds = True,
                            bufsize = 32768)

def generate_sugars_r3d(job_dir, chain_id):
    """Generate 'raw' input for tlsanim2r3d, but only for the non-animated
       sections, for _all_ chains
    """
    ## cmd: grep '^ATOM.................B.*' | rings3d -bases >>bases.r3d
    cmdlist = ["grep '^ATOM.................%s.*' %s/struct.pdb | %s -ribose >>%s/sugars.r3d 2>/dev/null" % (
              chain_id, job_dir, conf.RINGS3D, job_dir)]
    proc = subprocess.Popen(cmdlist,
                            shell = True,
                            stdin = subprocess.PIPE,
                            stdout = subprocess.PIPE,
                            stderr = subprocess.PIPE,
                            close_fds = True,
                            bufsize = 32768)

def SetStructureFile(webtlsmdd, job_id, struct_bin):
    """Creates job directory, saves structure file to the job directory,
    and sets all jdict defaults.
    """
    if not mysql.job_exists(job_id):
        return False

    try:
        os.chdir(conf.TLSMD_WORK_DIR)
    except OSError:
        return "Unable to change to conf.TLSMD_WORK_DIR = '%s'" % (
            conf.TLSMD_WORK_DIR)

    try:
        os.mkdir(job_id)
    except OSError:
        return "Unable to make job directory %s" % (job_id)

    job_dir = os.path.join(conf.TLSMD_WORK_DIR, job_id)
    os.chdir(job_dir)

    ## save PDB file
    pdb_filename = conf.PDB_FILENAME
    filobj = open(pdb_filename, "w")
    filobj.write(struct_bin.data)
    filobj.close()

    ## Generate summary/thumb 'struct.png' image
    if conf.THUMBNAIL:
        render_struct(job_dir)

    ## Generate 'struct.r3d' for Raster3D
    if conf.GEN_RAW_GREY:
        generate_raw_grey_struct(job_dir)

    ## set basic properties of the job
    job_url = "%s/%s" % (conf.TLSMD_WORK_URL, job_id)

    log_url = "%s/log.txt" % (job_url)
    log_file = "%s/log.txt" % (job_dir)
    if not os.path.exists(log_file):
        open(log_file, 'w').close() ## touch log.txt

    tarball_url       = "%s/%s.tar.gz" % (job_url, job_id)
    analysis_dir      = "%s/ANALYSIS" % (job_dir)
    analysis_base_url = "%s/ANALYSIS" % (job_url)
    analysis_url      = "%s/ANALYSIS/index.html" % (job_url)

    ## submission time and initial state
    submit_time = time.time()
    mysql.job_set_state(job_id, "submit1")
    mysql.job_set_submit_time(job_id, submit_time)

    ## This is for internal use only
    tm_struct = time.localtime(submit_time)
    submit_date = time.strftime("%Y-%m-%d %H:%M:%S", tm_struct)
    mysql.job_set_submit_date(job_id, submit_date)

    ## now load the structure and build the submission form
    try:
        struct = FileIO.LoadStructure(fil = pdb_filename)
    except:
        return "The Python Macromolecular Library was unable to load your structure file."

    if not struct.structure_id:
        struct.structure_id = "XXXX"
    mysql.job_set_structure_id(job_id, struct.structure_id)

    ## Select Chains for Analysis
    num_atoms = 0
    num_aniso_atoms = 0
    largest_chain_seen = 0

    mycb_desc = ""
    chains = []
    for chain in struct.iter_chains():
        naa = chain.count_amino_acids()
        nna = chain.count_nucleic_acids()
        num_frags = 0
        if naa > 0:
            num_frags = naa
        elif nna > 0:
            num_frags = nna

            ## this chain has nucleic acids in it, so generate r3d file for
            ## just the sugars
            generate_bases_r3d(job_dir, chain.chain_id)
            generate_sugars_r3d(job_dir, chain.chain_id)

        ## minimum number of residues (amino/nucleic) per chain
        if (naa > 0 and naa < conf.MIN_AMINO_PER_CHAIN) or\
           (nna > 0 and nna < conf.MIN_NUCLEIC_PER_CHAIN):
            continue

        largest_chain_seen = max(num_frags, largest_chain_seen)

        ## form name
        cb_name = 'CHAIN%s' % (chain.chain_id)
        mycb_desc = mycb_desc + chain.chain_id + ":"

        ## chains = "A:10:0:aa;B:20:1:na;C:30:0:na;"
        ## create chain description label cb_desc
        if naa > 0:
            mycb_desc = mycb_desc + str(num_frags) + ":1:aa;"
        elif nna > 0:
            mycb_desc = mycb_desc + str(num_frags) + ":1:na;"
        else:
            continue

        for atm in chain.iter_all_atoms():
            num_atoms += 1
            if atm.U is not None:
                num_aniso_atoms += 1

    if num_atoms < 1:
        webtlsmdd.remove_job(job_id)
        return 'Your submitted structure contained no atoms'

    if largest_chain_seen > conf.LARGEST_CHAIN_ALLOWED:
        webtlsmdd.remove_job(job_id)
        return 'Your submitted structure contained a chain exceeding the 1700 residue limit'

    mysql.job_set_chain_sizes(job_id, mycb_desc)

    ## set defaults
    mysql.job_set_user_name(job_id, "")
    mysql.job_set_email(job_id, "")
    mysql.job_set_user_comment(job_id, "")
    mysql.job_set_plot_format(job_id, "PNG")
    mysql.job_set_nparts(job_id, conf.globalconf.nparts)
    mysql.job_set_via_pdb(job_id, "0")

    mysql.job_set_private_job(job_id, "0")
    mysql.job_set_jmol_view(job_id, "0")
    mysql.job_set_jmol_animate(job_id, "0")
    mysql.job_set_histogram(job_id, "0")
    if conf.PRIVATE_JOBS:
        mysql.job_set_private_job(job_id, "1")
    if conf.globalconf.generate_jmol_view:
        mysql.job_set_jmol_view(job_id, "1")
    if conf.globalconf.generate_jmol_animate:
        mysql.job_set_jmol_animate(job_id, "1")
    if conf.globalconf.generate_histogram:
        mysql.job_set_histogram(job_id, "1")

    try:
        aniso_ratio = float(num_aniso_atoms) / float(num_atoms)
    except ZeroDivisionError:
        return 'Your submitted structure contained no atoms'

    if aniso_ratio > 0.90:
        mysql.job_set_tls_model(job_id, "ANISO")
    else:
        mysql.job_set_tls_model(job_id, "ISOT")

    mysql.job_set_weight_model(job_id, "NONE")
    mysql.job_set_include_atoms(job_id, "ALL")

    return ""

def RequeueJob(webtlsmdd, job_id):
    """Pushes job to the end of the list
    """
    if (webtlsmdd.jobdb.job_data_get(job_id,'state') == 'running'):
        return False
    else:
        gdict = webtlsmdd.jobdb.retrieve_globals()
        job_num = gdict['next_job_num'] 
        gdict['next_job_num'] = job_num + 1
        webtlsmdd.jobdb.store_globals(gdict)
        webtlsmdd.jobdb.job_data_set(job_id, 'job_num', job_num)
        return True

def RemoveJob(webtlsmdd, job_id):
    """Removes the job from both the database and working directory.
    If job is still running when this function is called, it will first call
    KillJob(), then remove the associated data and files.
    """
    if not mysql.job_exists(job_id):
        return False

    job_dir = os.path.join(conf.TLSMD_WORK_DIR, job_id)
    if job_dir and os.path.isdir(job_dir):
        for root, dirs, files in os.walk(job_dir, topdown = False):
            for name in files:
                os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))
        os.rmdir(job_dir)

    mysql.delete_jdict(job_id)
    return True

def SignalJob(webtlsmdd, job_id):
    """Causes a job stuck on a certain task to skip that step and move on to
       the next step. It will eventually have a state "warnings"
    """
    ## FIXME: Doesn't seem to work, 2009-06-12
    if not mysql.job_exists(job_id):
        return False

    job_dir = os.path.join(conf.TLSMD_WORK_DIR, job_id)
    if job_dir and os.path.isdir(job_dir):
        try:
            pid = int(mysql.job_get_pid(job_id))
        except:
            return False
        try:
            ## Send signal SIGUSR1 and try to continue to job process.
            os.kill(pid, SIGUSR1)
        except:
            return False

    return True

def KillJob(webtlsmdd, job_id):
    """Kills jobs in state "running" by pid and moves them to the 
       "Completed Jobs" section as "killed" state
    """
    if not mysql.job_exists(job_id):
        return False

    job_dir = os.path.join(conf.TLSMD_WORK_DIR, job_id)
    if job_dir and os.path.isdir(job_dir):
        try:
            if mysql.job_get_pid(job_id) == None:
                return False
            else:
                pid = int(mysql.job_get_pid(job_id))
        except:
            return False
        try:
            os.kill(pid, SIGHUP)
        except:
            return False

    return True

def Refmac5RefinementPrep(job_id, chain_ntls):
    """Called with a list of tuples (chain_id, ntls).
    Generates PDB and TLSIN files for refinement with REFMAC5.
    Returns a single string if there is an error, otherwise a
    dictionary of results is returned.
    """
    struct_id = mysql.job_get_structure_id(job_id)
    job_dir = os.path.join(conf.TLSMD_WORK_DIR, job_id)
    analysis_dir = os.path.join(job_dir, "ANALYSIS")
    job_url = os.path.join(conf.TLSMD_PUBLIC_URL, "jobs", job_id)
    analysis_base_url = "%s/ANALYSIS" % (job_url)

    if not os.path.isdir(analysis_dir):
        return "Job analysis directory does not exist"

    old_dir = os.getcwd()
    os.chdir(analysis_dir)

    ## input structure
    pdbin = "%s.pdb" % (struct_id)
    if not os.path.isfile(pdbin):
        pdbin = None
        for pdbx in glob.glob("*.pdb"):
            if len(pdbx) == 8:
                struct_id = pdbx[:4]
                pdbin = pdbx
                break
        if pdbin is None:
            os.chdir(old_dir)
            return "Input PDB File %s Not Found" % (pdbin)

    ## the per-chain TLSOUT files from TLSMD must be merged
    tlsins = []
    for chain_id, ntls in chain_ntls:
        tlsin = "%s_CHAIN%s_NTLS%d.tlsout" % (struct_id, chain_id, ntls)
        if not os.path.isfile(tlsin):
            os.chdir(old_dir)
            return "Input TLSIN File %s Not Found" % (tlsin)
        tlsins.append(tlsin)

    ## form unique pdbout/tlsout filenames
    outbase = job_id
    pdbout = "%s.pdb" % (outbase)

    ## the tlsout from this program is going to be the tlsin
    ## for refinement, so it's important for the filename to have
    ## the tlsin extension so the user is not confused
    tlsout    = "%s.tlsin"  % (outbase)
    phenixout = "%s.phenix" % (outbase)

    ## make urls for linking
    pdbout_url    = "%s/%s" % (analysis_base_url, pdbout)
    tlsout_url    = "%s/%s" % (analysis_base_url, tlsout)
    phenixout_url = "%s/%s" % (analysis_base_url, phenixout)

    ## create the REFMAC/PHENIX files
    tls_calcs.refmac5_prep(pdbin, tlsins, pdbout, tlsout)
    tls_calcs.phenix_prep(pdbin, tlsins, phenixout)

    os.chdir(old_dir)
    return dict(pdbout = pdbout,
                pdbout_url = pdbout_url,
                tlsout = tlsout,
                tlsout_url = tlsout_url,
                phenixout = phenixout,
                phenixout_url = phenixout_url)


class WebTLSMDDaemon(object):
    def __init__(self, db_file):
        self.db_file = db_file
        self.jobdb = None

    def job_list(self):
        """Returns a ordered list of all jdicts in the database
        """
        return self.jobdb.jdict_list()

    def job_pdb_list(self):
        """Returns a ordered list of all jdicts in the database
           containing 'via_pdb'
        """
        return self.jobdb.jdict_pdb_list()

    def job_new(self):
        return self.jobdb.job_new()

    def job_exists(self, job_id):
        return self.jobdb.job_exists(job_id)

    def job_get_dict(self, job_id):
        jdict = self.jobdb.retrieve_jdict(job_id)
        if jdict is None:
            return False
        return jdict

    def get_next_queued_job_id(self):
        job_list = self.job_list()
        for jdict in job_list:
            if jdict.get("state") == "queued":
                return jdict["job_id"]
        return ""

    def set_structure_file(self, job_id, struct_bin):
        """Creates job directory, saves structure file to the job directory,
        and sets all jdict defaults.
        """
        return SetStructureFile(self, job_id, struct_bin)

    def remove_job(self, job_id):
        """Removes the job from both the database and working directory.
        If job is still running when this function is called, it will first call
        KillJob(), then remove the associated data and files.
        """
        try:
            KillJob(self, job_id)
        except:
            pass
        return RemoveJob(self, job_id)

    def signal_job(self, job_id):
        """Signals a job stuck on a certain task to skip that step and move on
           to the next step. It will eventually have a state "warnings".
        """
        return SignalJob(self, job_id)

    def kill_job(self, job_id):
        """Kills jobs in state "running" by pid and moves them to the
           "Completed Jobs" section as "killed" state.
        """
        return KillJob(self, job_id)

    def job_set_remote_addr(self, job_id, remote_addr):
        self.jobdb.job_data_set(job_id, "ip_addr", remote_addr)
        return remote_addr
    def job_get_remote_addr(self, job_id):
        return self.jobdb.job_data_get(job_id, "ip_addr")

    def job_set_state(self, job_id, state):
        self.jobdb.job_data_set(job_id, "state", state)
        return state
    def job_get_state(self, job_id):
        return self.jobdb.job_data_get(job_id, "state")

    def job_get_analysis_dir(self, job_id):
        return self.jobdb.job_data_get(job_id, "analysis_dir")

    def job_get_analysis_url(self, job_id):
        return self.jobdb.job_data_get(job_id, "analysis_url")

    def job_get_analysis_base_url(self, job_id):
        return self.jobdb.job_data_get(job_id, "analysis_base_url")

    def job_get_job_url(self, job_id):
        return self.jobdb.job_data_get(job_id, "job_url")

    def job_get_job_dir(self, job_id):
        return self.jobdb.job_data_get(job_id, "job_dir")

    def job_get_pdb_dir(self, job_id):
        return self.jobdb.job_data_get(job_id, "pdb_dir")

    def job_set_pdb_dir(self, job_id, pdb_id):
        directory = os.path.join(conf.WEBTLSMDD_PDB_DIR, pdb_id)
        self.jobdb.job_data_set(job_id, "pdb_dir", directory)
        return directory

    def job_set_pid(self, job_id, os_pid):
        self.jobdb.job_data_set(job_id, "pid", os_pid)
        return os_pid
    def job_get_pid(self, job_id):
        return self.jobdb.job_data_get(job_id, "pid")

    def job_get_log_url(self, job_id):
        return self.jobdb.job_data_get(job_id, "log_url")

    def job_get_tarball_url(self, job_id):
        return self.jobdb.job_data_get(job_id, "tarball_url")

    def job_set_chains(self, job_id, chains):
        self.jobdb.job_data_set(job_id, "chains", chains)
        return chains
    def job_get_chains(self, job_id):
        return self.jobdb.job_data_get(job_id, "chains")

    def job_set_user(self, job_id, user):
        self.jobdb.job_data_set(job_id, "user", user)
        return user
    def job_get_user(self, job_id):
        return self.jobdb.job_data_get(job_id, "user")

    def job_set_user_name(self, job_id, user_name):
        self.jobdb.job_data_set(job_id, "user_name", user_name)
        return user_name
    def job_get_user_name(self, job_id):
        return self.jobdb.job_data_get(job_id, "user_name")

    def job_set_user_comment(self, job_id, user_comment):
        self.jobdb.job_data_set(job_id, "user_comment", user_comment)
        return user_comment
    def job_get_user_comment(self, job_id):
        return self.jobdb.job_data_get(job_id, "user_comment")

    def job_set_email(self, job_id, email_address):
        self.jobdb.job_data_set(job_id, "email", email_address)
        return email_address
    def job_get_email(self, job_id):
        return self.jobdb.job_data_get(job_id, "email")

    def job_set_private_job(self, job_id, private_job):
        self.jobdb.job_data_set(job_id, "private_job", private_job)
        return private_job
    def job_get_private_job(self, job_id):
        return self.jobdb.job_data_get(job_id, "private_job")

    def job_set_structure_id(self, job_id, structure_id):
        self.jobdb.job_data_set(job_id, "structure_id", structure_id)
        return structure_id
    def job_get_structure_id(self, job_id):
        return self.jobdb.job_data_get(job_id, "structure_id")

    ## This data comes from the HEADER line
    ## It will be a four character string and will be "xxxx" for RefMac.
    def job_set_header_id(self, job_id, header_id):
        self.jobdb.job_data_set(job_id, "header_id", header_id)
        return header_id
    def job_get_header_id(self, job_id):
        return self.jobdb.job_data_get(job_id, "header_id")

    ## This data comes from the REMARK line
    ## It will be a string something like "1.80"
    def job_set_resolution(self, job_id, resolution):
        self.jobdb.job_data_set(job_id, "resolution", resolution)
        return resolution
    def job_get_resolution(self, job_id):
        return self.jobdb.job_data_get(job_id, "resolution")

    def job_set_initial_residuals(self, job_id, initial_residuals):
        self.jobdb.job_data_set(job_id, "initial_residuals", initial_residuals)
        return initial_residuals
    def job_get_initial_residuals(self, job_id):
        return self.jobdb.job_data_get(job_id, "initial_residuals")

    def job_set_final_residuals(self, job_id, final_residuals):
        self.jobdb.job_data_set(job_id, "final_residuals", final_residuals)
        return final_residuals
    def job_get_final_residuals(self, job_id):
        return self.jobdb.job_data_get(job_id, "final_residuals")

    def job_set_stddev_bfact(self, job_id, stddev_bfact):
        self.jobdb.job_data_set(job_id, "stddev_bfact", stddev_bfact)
        return stddev_bfact
    def job_get_stddev_bfact(self, job_id):
        return self.jobdb.job_data_get(job_id, "stddev_bfact")

    def job_set_chain_max_segs(self, job_id, chain_max_segs):
        self.jobdb.job_data_set(job_id, "chain_max_segs", chain_max_segs)
        return chain_max_segs
    def job_get_chain_max_segs(self, job_id):
        return self.jobdb.job_data_get(job_id, "chain_max_segs")

    ## FIXME: Why are there two of these? #1
    def job_set_tls_model(self, job_id, tls_model):
        self.jobdb.job_data_set(job_id, "tls_model", tls_model)
        return tls_model
    def job_get_tls_model(self, job_id):
        return self.jobdb.job_data_get(job_id, "tls_model")

    ## FIXME: Why are there two of these? #2
    def job_set_tls_model(self, job_id, tls_model):
        self.jobdb.job_data_set(job_id, "tls_model", tls_model)
        return tls_model
    def job_get_tls_model(self, job_id):
        return self.jobdb.job_data_get(job_id, "tls_model")

    def job_set_weight_model(self, job_id, weight_model):
        self.jobdb.job_data_set(job_id, "weight", weight_model)
        return weight_model
    def job_get_weight_model(self, job_id):
        return self.jobdb.job_data_get(job_id, "weight")

    def job_set_include_atoms(self, job_id, include_atoms):
        self.jobdb.job_data_set(job_id, "include_atoms", include_atoms)
        return include_atoms
    def job_get_include_atoms(self, job_id):
        return self.jobdb.job_data_get(job_id, "include_atoms")

    def job_set_plot_format(self, job_id, plot_format):
        self.jobdb.job_data_set(job_id, "plot_format", plot_format)
        return plot_format
    def job_get_plot_format(self, job_id):
        return self.jobdb.job_data_get(job_id, "plot_format")

    ## Generate JMol view feature (default=True)
    def job_set_jmol_view(self, job_id, generate_jmol_view):
        self.jobdb.job_data_set(job_id, "generate_jmol_view", generate_jmol_view)
        return generate_jmol_view
    def job_get_jmol_view(self, job_id):
        return self.jobdb.job_data_get(job_id, "generate_jmol_view")

    ## Generate JMol animate feature (default=True)
    def job_set_jmol_animate(self, job_id, generate_jmol_animate):
        self.jobdb.job_data_set(job_id, "generate_jmol_animate", generate_jmol_animate)
        return generate_jmol_animate
    def job_get_jmol_animate(self, job_id):
        return self.jobdb.job_data_get(job_id, "generate_jmol_animate")

    ## Generate Histogram plot (default=False)
    def job_set_histogram(self, job_id, generate_histogram):
        self.jobdb.job_data_set(job_id, "generate_histogram", generate_histogram)
        return generate_histogram
    def job_get_histogram(self, job_id):
        return self.jobdb.job_data_get(job_id, "generate_histogram")

    ## Set/get number of partitions/chain
    def job_set_nparts(self, job_id, nparts):
        self.jobdb.job_data_set(job_id, "nparts", nparts)
        return nparts
    def job_get_nparts(self, job_id):
        return self.jobdb.job_data_get(job_id, "nparts")

    def job_set_run_time_begin(self, job_id, run_time_begin):
        self.jobdb.job_data_set(job_id, "run_time_begin", run_time_begin)
        return run_time_begin
    def job_get_run_time_begin(self, job_id):
        return self.jobdb.job_data_get(job_id, "run_time_begin")

    def job_set_run_time_end(self, job_id, run_time_end):
        self.jobdb.job_data_set(job_id, "run_time_end", run_time_end)
        return run_time_end
    def job_get_run_time_end(self, job_id):
        return self.jobdb.job_data_get(job_id, "run_time_end")

    def job_set_via_pdb(self, job_id, bool):
        self.jobdb.job_data_set(job_id, "via_pdb", bool)
        return bool

    def requeue_job(self, job_id):
        """Pushes the job to the back of the queue"""
        return RequeueJob(self, job_id)

    def refmac5_refinement_prep(self, job_id, chain_ntls):
        """Called with a list of tuples (chain_id, ntls).
        Generates PDB and TLSIN files for refinement with REFMAC5.
        Returns a single string if there is an error, otherwise a
        dictionary of results is returned.
        """
        return Refmac5RefinementPrep(job_id, chain_ntls)

    def pdb_exists(self, pdbid):
        try:
            f = open(conf.WEBTLSMDD_PDBID_FILE, 'r')
        except IOError:
            ## if it doesn't exist create it
            f = open(conf.WEBTLSMDD_PDBID_FILE, 'w+')

        for line in f:
            if line == pdbid + '\n':
                f.close()
                return True

        f.close()

        return False

    def fetch_pdb(self, pdbid):
        """Retrives the PDB file from RCSB"""
        try:
            cdata = urllib.urlopen("%s/%s.pdb.gz" % (conf.GET_PDB_URL,pdbid)).read()
            sys.stdout.write("FOUND PDB: %s" % pdbid)
            data = gzip.GzipFile(fileobj = StringIO.StringIO(cdata)).read()
        except IOError:
            return xmlrpclib.Binary("")
        return xmlrpclib.Binary(data)

    def set_pdb_db(self, pdbid):
        try:
            f = open(conf.WEBTLSMDD_PDBID_FILE, 'a')
        except IOError:
            return False
        f.write(pdbid + '\n')
        f.close()
        return True

class WebTLSMD_XMLRPCRequestHandler(SimpleXMLRPCServer.SimpleXMLRPCRequestHandler):
    """Override the standard XMLRPC request handler to open the database before
    calling the method.
    """
    ## TODO: Can this be removed? 2009-06-01
    def handle(self):
        #self.server.webtlsmdd.jobdb = JobDatabase(self.server.webtlsmdd.db_file)
        return SimpleXMLRPCServer.SimpleXMLRPCRequestHandler.handle(self)


class WebTLSMD_XMLRPCServer(
            SocketServer.ForkingMixIn,
            SimpleXMLRPCServer.SimpleXMLRPCServer):
    """Use customized XMLRPC server which forks for requests and uses the
       customized request handler.
    """
    def __init__(self, host_port):
        SimpleXMLRPCServer.SimpleXMLRPCServer.__init__(
            self,
            host_port,
            WebTLSMD_XMLRPCRequestHandler,
            False)


def daemon_main():
    rtype, baseurl, port = conf.WEBTLSMDD.split(":")
    host_port = ("localhost", int(port))

    sys.stdout.write("STARTING webtlsmdd.py DAEMON..................: %s\n" % misc.timestamp())
    sys.stdout.write("webtlsmdd.py xmlrpc server version %s\n" % (const.VERSION))
    sys.stdout.write("using database file...........................: %s\n" % (conf.WEBTLSMDD_DATABASE))
    sys.stdout.write("listening for incoming connections at URL.....: %s\n" % (conf.WEBTLSMDD))
    sys.stdout.write("job (working) directory.......................: %s\n" % (conf.TLSMD_WORK_DIR))

    os.chdir(conf.TLSMD_WORK_DIR)

    ## Switched from handle_SIGCHLD to SIG_IGN. Christoph Champ, 2008-03-10
    signal.signal(signal.SIGCHLD, SIG_IGN)

    webtlsmdd = WebTLSMDDaemon(conf.WEBTLSMDD_DATABASE)    

    try:
        xmlrpc_server = WebTLSMD_XMLRPCServer(host_port)
    except socket.error:
        sys.stderr.write("[ERROR] unable to bind to host,port: %s\n" % (str(host_port)))
        raise SystemExit

    xmlrpc_server.webtlsmdd = webtlsmdd
    xmlrpc_server.register_instance(webtlsmdd)
    xmlrpc_server.serve_forever()

def main():
    try:
        daemon_main()
    except:
        email.SendTracebackEmail("webtlsmdd.py exception")
        raise

def inspect():
    database_path = sys.argv[2]

    webtlsmdd = WebTLSMDDaemon(database_path)

    if sys.argv[1] == "list":
        for dbkey in webtlsmdd.db.keys():
            jdict = webtlsmdd.job_get_dict(dbkey)
            print dbkey, jdict.get("email")

    if sys.argv[1] == "remove":
        dbkey = sys.argv[2]
        del webtlsmdd.db[dbkey]
        webtlsmdd.db.sync()

def usage():
    print "webtlsmdd.py [list | remove] args..."

if __name__=="__main__":
    if len(sys.argv) == 1:
        try:
            main()
        except KeyboardInterrupt:
            raise SystemExit
    else:
        inspect()
    sys.exit(0)
