#!/usr/bin/python
# coding=UTF-8 
## TLS Motion Determination (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys
import time
import string
import random
import traceback
from signal import SIG_IGN # Needed for daemon_main(). Christoph Champ, 2008-03-10
from signal import SIGHUP # Needed for KillJob(). Christoph Champ, 2008-03-18
import signal
import cPickle
import bsddb
import socket
import xmlrpclib
import SocketServer
import SimpleXMLRPCServer
import urllib
import gzip
import StringIO

from mmLib import FileIO
from tlsmdlib import conf, const, tls_calcs, email


def fatal(text):
    sys.stderr.write("[FATAL ERROR] %s\n" % (text))
    raise SystemExit


def generate_security_code(code_length = 8):
    """Generates a random 8
    """
    random.seed()
    codelist = list(5 * string.ascii_letters) 
    random.shuffle(codelist)
    code = "".join(random.sample(codelist, code_length)) 
    return code


class JobDatabase(object):
    def __init__(self, db_file):
        self.db_file = db_file
        self.db = bsddb.hashopen(self.db_file, "c")

        if self.retrieve_globals() is None:
            self.init_globals()

    def __store_dict(self, dbkey, dictx):
        self.db[dbkey] = cPickle.dumps(dictx)
        self.db.sync()
    
    def __retrieve_dict(self, dbkey):
        try:
            data = self.db[dbkey]
        except KeyError:
            return None
        rdict = cPickle.loads(data)
        for key, value in rdict.iteritems():
            if isinstance(rdict, unicode):
                rdict[key] = ""
        return rdict

    def init_globals(self):
        gdict = dict(next_job_num = 1)
        self.store_globals(gdict)

    def retrieve_globals(self):
        """Retrieves the globals dictionary from the database or
        creates one and returns it if one does not exist.
        """
        return self.__retrieve_dict("GLOBALS")

    def store_globals(self, gdict):
        self.__store_dict("GLOBALS", gdict)

    def store_jdict(self, jdict):
        job_id = jdict["job_id"]
        assert job_id.startswith("TLSMD")
        self.__store_dict(job_id, jdict)
    
    def retrieve_jdict(self, job_id):
        assert job_id.startswith("TLSMD")
        return self.__retrieve_dict(job_id)

    def delete_jdict(self, job_id):
        assert job_id.startswith("TLSMD")
        if not self.db.has_key(job_id):
            return False
        del self.db[job_id]
        self.db.sync()
        return True

    def job_exists(self, job_id):
        return self.db.has_key(job_id)

    def jdict_list(self):
        """Returns a sorted list of all jdicts in the database
        """
        ## retrieve all jdicts from database and 
        listx = []
        for dbkey in self.db.keys():
            if dbkey.startswith("TLSMD"):
                jdict = self.retrieve_jdict(dbkey)

                job_id = jdict["job_id"]

                if jdict.has_key("job_num"):
                    job_num = jdict["job_num"]
                else:
                    job_nums = job_id[5:]
                    j = job_nums.find("_")
                    job_nums = job_nums[:j]
                    job_num = int(job_nums)
		    
                listx.append((job_num, jdict))

        listx.sort()

        job_list = []
        for job_num, jdict in listx:
            job_list.append(jdict)

        return job_list

    def job_new(self):
        gdict = self.retrieve_globals()
        job_num = gdict["next_job_num"]
        gdict["next_job_num"] =  job_num + 1
        self.store_globals(gdict)

        ## assign job_id
        security_code = generate_security_code()
        job_id = "TLSMD%d_%s" % (job_num, security_code)
	    
        ## create job dictionary
        jdict = {}
        jdict["job_id"] = job_id
        jdict["job_num"] = job_num
        self.store_jdict(jdict)
        return job_id

    def job_data_set(self, job_id, key, value):
        jdict = self.retrieve_jdict(job_id)
        if jdict is None:
            return False
        jdict[key] = value
        self.store_jdict(jdict)
        return True

    def job_data_get(self, job_id, key):
        jdict = self.retrieve_jdict(job_id)
        if jdict is None:
            return ""
        value = jdict.get(key)
	if isinstance(value, unicode):
            return ""
        if value is None:
            return ""
        return value


def SetStructureFile(webtlsmdd, job_id, struct_bin):
    """Creates job directory, saves structure file to the job directory,
    and sets all jdict defaults.
    """
    if not webtlsmdd.job_exists(job_id):
        return False

    try:
        os.chdir(conf.TLSMD_WORK_DIR)
    except OSError:
        return "Unable to change to conf.TLSMD_WORK_DIR = '%s'" % (conf.TLSMD_WORK_DIR)

    try:
        os.mkdir(job_id)
    except OSError:
        return "Unable to make job directory %s" % (job_id)

    job_dir = os.path.join(conf.TLSMD_WORK_DIR, job_id)
    os.chdir(job_dir)
    webtlsmdd.jobdb.job_data_set(job_id, "job_dir", job_dir)

    ## save PDB file
    pdb_filename = "struct.pdb"
    webtlsmdd.jobdb.job_data_set(job_id, "pdb_filename", pdb_filename)
    filobj = open(pdb_filename, "w")
    filobj.write(struct_bin.data)
    filobj.close()

    ## set basic properties of the job
    job_url = "%s/%s" % (conf.TLSMD_WORK_URL, job_id)
    webtlsmdd.jobdb.job_data_set(job_id, "job_url", job_url)

    log_url = "%s/log.txt" % (job_url)
    webtlsmdd.jobdb.job_data_set(job_id, "log_url", log_url)

    ## create tarball path and url. Christoph Champ, 2007-12-03
    tarball_url = "%s/%s.tar.gz" % (job_url,job_id)
    webtlsmdd.jobdb.job_data_set(job_id, "tarball_url", tarball_url)

    analysis_dir = "%s/ANALYSIS" % (job_dir)
    webtlsmdd.jobdb.job_data_set(job_id, "analysis_dir", analysis_dir)

    analysis_base_url = "%s/ANALYSIS" % (job_url)
    webtlsmdd.jobdb.job_data_set(job_id, "analysis_base_url", analysis_base_url)

    analysis_url = "%s/ANALYSIS/index.html" % (job_url)
    webtlsmdd.jobdb.job_data_set(job_id, "analysis_url", analysis_url)

    webtlsmdd.jobdb.job_data_set(job_id, "version", const.VERSION)

    ## submission time and initial state
    webtlsmdd.jobdb.job_data_set(job_id, "state", "submit1")
    webtlsmdd.jobdb.job_data_set(job_id, "submit_time", time.time())

    ## now load the structure and build the submission form
    try:
        struct = FileIO.LoadStructure(fil = pdb_filename)
    except:
        return "The Python Macromolecular Library was unable to load your structure file."

    if not struct.structure_id:
        struct.structure_id = "XXXX"
    webtlsmdd.jobdb.job_data_set(job_id, "structure_id", struct.structure_id)

    ## Select Chains for Analysis
    num_atoms = 0
    num_aniso_atoms = 0
    largest_chain_seen = 0

    chains = []
    for chain in struct.iter_chains():
        naa = chain.count_amino_acids()
        nna = chain.count_nucleic_acids()
        num_frags = 0
        if naa > 0:
            num_frags = naa
        elif nna > 0:
            num_frags = nna

        if num_frags < 10:
            continue

        largest_chain_seen = max(num_frags, largest_chain_seen)

        ## form name
        cb_name = 'CHAIN%s' % (chain.chain_id)

        ## create chain description label cb_desc
        if naa > 0:
            cb_desc = 'Chain %s (%d Amino Acid Residues)' % (chain.chain_id, num_frags)
        elif nna > 0:
            cb_desc = 'Chain %s (%d Nucleic Acid Residues)' % (chain.chain_id, num_frags)
        else:
            continue
            
        for atm in chain.iter_all_atoms():
            num_atoms += 1
            if atm.U is not None:
                num_aniso_atoms += 1

        listx = []
        i = 0
        for frag in chain.iter_fragments():
            i += 1
            if i > 5:
                break
            listx.append(frag.res_name)
        cb_preview = string.join(listx, " ")

        cdict = {}
        chains.append(cdict)
        cdict["chain_id"] = chain.chain_id
        cdict["length"] = num_frags
        cdict["name"] = cb_name
        cdict["desc"] = cb_desc
        cdict["preview"] = cb_preview
        cdict["selected"] = True

    if num_atoms < 1:
        webtlsmdd.remove_job(job_id)
        return 'Your submitted structure contained no atoms'

    if largest_chain_seen > 1700:
        webtlsmdd.remove_job(job_id)
        return 'Your submitted structure contained a chain exceeding the 1700 residue limit'

    webtlsmdd.jobdb.job_data_set(job_id, "chains", chains)

    ## defaults
    webtlsmdd.jobdb.job_data_set(job_id, "user", "")
    webtlsmdd.jobdb.job_data_set(job_id, "passwd", "")
    webtlsmdd.jobdb.job_data_set(job_id, "email", "")
    webtlsmdd.jobdb.job_data_set(job_id, "comment", "")
    webtlsmdd.jobdb.job_data_set(job_id, "private_job", False) ## This is overwritten in html.py. Christoph Champ, 2008-02-09
    webtlsmdd.jobdb.job_data_set(job_id, "plot_format", "PNG")

    try:
        aniso_ratio = float(num_aniso_atoms) / float(num_atoms)
    except ZeroDivisionError:
        return 'Your submitted structure contained no atoms'

    if aniso_ratio > 0.90:
        webtlsmdd.jobdb.job_data_set(job_id, "tls_model", "ANISO")
    else:
        webtlsmdd.jobdb.job_data_set(job_id, "tls_model", "ISOT")

    webtlsmdd.jobdb.job_data_set(job_id, "weight", "")
    webtlsmdd.jobdb.job_data_set(job_id, "include_atoms", "ALL")

    return ""

def RequeueJob(webtlsmdd, job_id):
    """Pushes job to the end of the list"""

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
    """
    if not webtlsmdd.jobdb.job_exists(job_id):
        return False

    job_dir = webtlsmdd.job_get_job_dir(job_id)
    if job_dir and job_dir.startswith(conf.TLSMD_WORK_DIR) and os.path.isdir(job_dir):
        for root, dirs, files in os.walk(job_dir, topdown = False):
            for name in files:
                os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))
        os.rmdir(job_dir)

    webtlsmdd.jobdb.delete_jdict(job_id)
    return True

def KillJob(webtlsmdd, job_id):
    """Kills jobs in state "running" by pid and moves them to the "Completed Jobs" section
       as "killed" state
    """
    ## Christoph Champ, 2008-03-13
    ## DEBUG
    #debug=open("/tmp/debug.log","a+")
    #debug.write("KillJob: DEBUG1\n")
    #debug.close()

    if not webtlsmdd.jobdb.job_exists(job_id):
        return False

    job_dir = webtlsmdd.job_get_job_dir(job_id)
    if job_dir and job_dir.startswith(conf.TLSMD_WORK_DIR) and os.path.isdir(job_dir):
        try:
	    ## Switched to storing pid in database. Christoph Champ, 2008-03-14
	    tmp_pid=webtlsmdd.job_get_pid(job_id)
	    pid=int(tmp_pid)
        except:
	    return False
        try:
	    ## Switched to SIGHUP because SIGUSR1 was not working. Christoph Champ, 2008-03-18
	    os.kill(pid,SIGHUP)
        except:
	    return False

    ## We want to keep the job_id around in order to inform the user that their job has been "killed". Christoph Champ, 2008-03-13
    webtlsmdd.job_set_state(job_id,"killed")
    return True

def Refmac5RefinementPrep(webtlsmdd, job_id, chain_ntls):
    """Called with a list of tuples (chain_id, ntls).
    Generates PDB and TLSIN files for refinement with REFMAC5.
    Returns a single string if there is an error, otherwise a
    dictionary of results is returned.
    """
    struct_id = webtlsmdd.job_get_structure_id(job_id)
    analysis_dir = webtlsmdd.job_get_analysis_dir(job_id)
    analysis_base_url = webtlsmdd.job_get_analysis_base_url(job_id)

    if not os.path.isdir(analysis_dir):
        return "Job analysis directory does not exist"

    old_dir = os.getcwd()
    os.chdir(analysis_dir)

    ## input structure
    pdbin  = "%s.pdb" % (struct_id)
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
    listx = [struct_id]
    for chain_id, ntls in chain_ntls:
        listx.append("CHAIN%s" % (chain_id))
        listx.append("NTLS%d" % (ntls))
    outbase ="_".join(listx)
    pdbout = "%s.pdb" % (outbase)

    ## the tlsout from this program is going to be the tlsin
    ## for refinement, so it's important for the filename to have
    ## the tlsin extension so the user is not confused
    tlsout = "%s.tlsin" % (outbase)
    phenixout = "%s.phenix" % (outbase) ## PHENIX, Christoph Champ, 2007-11-06

    ## make urls for linking
    pdbout_url = "%s/%s" % (analysis_base_url, pdbout)
    tlsout_url = "%s/%s" % (analysis_base_url, tlsout)
    phenixout_url = "%s/%s" % (analysis_base_url, phenixout) ## PHENIX, Christoph Champ, 2007-11-06

    ## create the files
    tls_calcs.refmac5_prep(pdbin, tlsins, pdbout, tlsout)
    tls_calcs.phenix_prep(pdbin, tlsins, phenixout) ## PHENIX, Christoph Champ, 2007-11-06

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
	    ## Changed. Christoph Champ, 2008-01-29
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
        """
        return RemoveJob(self, job_id)

    def kill_job(self, job_id):
	"""Kills jobs in state "running" by pid and moves them to the "Completed Jobs" section
	   as "killed" state
	"""
	## Christoph Champ, 2008-03-07
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

    def job_get_job_dir(self, job_id):
        return self.jobdb.job_data_get(job_id, "job_dir")
    
    def job_get_pdb_dir(self, job_id):
        return self.jobdb.job_data_get(job_id, "pdb_dir")
    
    def job_set_pdb_dir(self, job_id, pdb_id):
        directory = os.path.join(conf.WEBTLSMDD_PDB_DIR, pdb_id)
        self.jobdb.job_data_set(job_id, "pdb_dir", directory)
        return directory

    ## New pid field added. Christoph Champ, 2008-03-14
    def job_set_pid(self, job_id, os_pid):
        self.jobdb.job_data_set(job_id, "pid", os_pid)
        return os_pid
    def job_get_pid(self, job_id):
        return self.jobdb.job_data_get(job_id, "pid")

    def job_get_log_url(self, job_id):
        return self.jobdb.job_data_get(job_id, "log_url")

    ## tarball url. Christoph Champ, 2007-12-03
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

    ## New user_comment field added. Christoph Champ, 2007-12-18
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

    def job_set_tls_model(self, job_id, tls_model):
        self.jobdb.job_data_set(job_id, "tls_model", tls_model)
        return tls_model
    def job_get_tls_model(self, job_id):
        return self.jobdb.job_data_get(job_id, "tls_model")

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
        return Refmac5RefinementPrep(self, job_id, chain_ntls)

    def pdb_exists(self, pdbid):
        try:
            f = open(conf.WEBTLSMDD_PDBID_FILE, 'r')
        except IOError:
            # if it doesn't exist create it
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
	    ## Changed to global variable. Christoph Champ, 2008-03-10
            #cdata = urllib.urlopen("http://www.rcsb.org/pdb/files/%s.pdb.gz" % (pdbid)).read()
            cdata = urllib.urlopen("%s/%s.pdb.gz" % (conf.GET_PDB_URL,pdbid)).read()
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
    def handle(self):
        self.server.webtlsmdd.jobdb = JobDatabase(self.server.webtlsmdd.db_file)
        return SimpleXMLRPCServer.SimpleXMLRPCRequestHandler.handle(self)


class WebTLSMD_XMLRPCServer(
            SocketServer.ForkingMixIn,
            SimpleXMLRPCServer.SimpleXMLRPCServer):
    """Use customized XMLRPC server which forks for requests and uses the customized
    request handler.
    """
    def __init__(self, host_port):
        SimpleXMLRPCServer.SimpleXMLRPCServer.__init__(
            self,
            host_port,
            WebTLSMD_XMLRPCRequestHandler,
            False)
                
## Making sure this is never used/needed. Christoph Champ, 2008-03-17
#def handle_SIGCHLD(signum, frame):
#    try:
#        os.waitpid(-1, os.WNOHANG)
#    except:
#	pass

def daemon_main():
    rtype, baseurl, port = conf.WEBTLSMDD.split(":")
    host_port = ("localhost", int(port))

    sys.stdout.write("webtlsmdd.py xmlrpc server version %s\n" % (const.VERSION))
    sys.stdout.write("using database file...........................: %s\n" % (conf.WEBTLSMDD_DATABASE))
    sys.stdout.write("listening for incoming connections at URL.....: %s\n" % (conf.WEBTLSMDD))
    sys.stdout.write("job (working) directory.......................: %s\n" % (conf.TLSMD_WORK_DIR))

    os.chdir(conf.TLSMD_WORK_DIR)

    ## Switching to a different signal handler. Christoph Champ, 2008-03-10
    #signal.signal(signal.SIGCHLD, handle_SIGCHLD)
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
