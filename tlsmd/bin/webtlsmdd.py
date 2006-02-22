#!/usr/bin/env python

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

import cPickle
import bsddb

import xmlrpclib
import SimpleXMLRPCServer

from mmLib import FileIO

from tlsmdlib import conf, const


def debug(text):
    return
    open("debug.txt", "a").write(text)


def generate_security_code(code_length = 8):
    """Generates a random 8
    """
    codelist = list(5 * string.ascii_letters) 
    random.shuffle(codelist)
    code = "".join(random.sample(codelist, code_length)) 
    return code


class WebTLSMD_XMLRPCServer(SimpleXMLRPCServer.SimpleXMLRPCServer):
    def handle_error(self, request, client_address):
        print "error!"
        TCPServer.handle_error(self, request, client_address)


class WebTLSMDDaemon2(object):
    def __init__(self, db_file):
        self.db_file = db_file
        self.db = bsddb.hashopen(self.db_file, "c")
        self.retrieve_globals()

    def store_dict(self, dbkey, dictx):
        self.db[dbkey] = cPickle.dumps(dictx)
        self.db.sync()
    
    def retrieve_dict(self, dbkey):
        try:
            data = self.db[dbkey]
        except KeyError:
            return None
        return cPickle.loads(data)

    def retrieve_globals(self):
        """Retrieves the globals dictionary from the database or
        creates one and returns it if one does not exist.
        """
        gdict = self.retrieve_dict("GLOBALS")

        ## no global dictionary; create one which is consistent with
        ## the TLSMD jobs in the database
        if gdict == None:
            gdict = {}

            job_list = self.job_list()
            if len(job_list) == 0:
                gdict["next_job_num"] = 1
            else:
                jdict = job_list[-1]
                job_id = jdict["job_id"]
                job_num = int(job_id[5:])
                gdict["next_job_num"] = job_num + 1

            self.store_globals(gdict)

        debug("DATABASE GLOBALS\n")
	for key, value in gdict.items():
	    debug("%20s : %s\n" % (key, value))

        return gdict

    def store_globals(self, gdict):
        self.store_dict("GLOBALS", gdict)

    def store_jdict(self, jdict):
        dbkey = jdict["job_id"]
        self.store_dict(jdict["job_id"], jdict)
    
    def retrieve_jdict(self, job_id):
        jdict = self.retrieve_dict(job_id)
        for key, val in jdict.items():
            if val==None:
                debug("ERROR: JDICT for %s has key/value %s:None" % (jdict["job_id"], str(key)))
        return jdict

    def job_list(self):
        """Returns a ordered list of all jdicts in the database
        """
	debug("calling job_list\n")

        ## retrieve all jdicts from database and 
        listx = []
        for dbkey in self.db.keys():
            if dbkey.startswith("TLSMD"):
                debug("retrieving %s\n" % (dbkey))

                jdict = self.retrieve_jdict(dbkey)
                debug(str(jdict))

                job_id = jdict["job_id"]
                job_nums = job_id[5:]
		j = job_nums.find("_")
		if j>0:
                    job_nums = job_nums[:j]
	
                debug("job num = %s\n" % (job_nums))
	
		try:
                    job_num = int(job_nums)
                except ValueError:
                    debug("ERROR: unable to determine job number for JOBID %s\n" % (dbkey))
		    continue
		    
                listx.append((job_num, jdict))

        listx.sort()

        job_list = []
        for job_num, jdict in listx:
            job_list.append(jdict)

	debug(str(job_list))
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
	self.store_jdict(jdict)
        return job_id

    def job_exists(self, job_id):
        return self.db.has_key(job_id)

    def job_get_dict(self, job_id):
        jdict = self.retrieve_jdict(job_id)
        if jdict==None:
            return False
        return jdict
    
    def job_get_dict_index(self, i):
        job_list = self.job_list()
        try:
            jdict = job_list[i]
        except IndexError:
            return False
        return jdict

    def job_delete(self, job_id):
        if not self.db.has_key(job_id):
            return False

        del self.db[job_id]
        self.db.sync()
        return True

    def job_data_set(self, job_id, key, value):
        jdict = self.retrieve_jdict(job_id)
        if jdict == None:
            return False
        jdict[key] = value
        self.store_jdict(jdict)
        return True

    def job_data_get(self, job_id, key):
        jdict = self.retrieve_jdict(job_id)
        if jdict == None:
            return False
        return jdict.get(key, False)

    def get_next_queued_job_id(self):
        job_list = self.job_list()
	for jdict in job_list:
            if jdict.get("state") == "running" or jdict.get("state") == "queued":
                return jdict["job_id"]
	return False

    def set_structure_file(self, job_id, struct_bin):
        """Creates job directory, saves structure file to the job directory,
        and sets all jdict defaults.
        """
        if not self.job_exists(job_id):
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
        self.job_data_set(job_id, "job_dir", job_dir)

        ## save PDB file
        pdb_filename = "struct.pdb"
        self.job_data_set(job_id, "pdb_filename", pdb_filename)
        filobj = open(pdb_filename, "w")
        filobj.write(struct_bin.data)
        filobj.close()
        
        job_url = "%s/%s" % (conf.TLSMD_WORK_URL, job_id)
        self.job_data_set(job_id, "job_url", job_url)

        log_url = "%s/log.txt" % (job_url)
        self.job_data_set(job_id, "log_url", log_url)

        analysis_dir = "%s/ANALYSIS" % (job_dir)
        self.job_data_set(job_id, "analysis_dir", analysis_dir)

        analysis_base_url = "%s/ANALYSIS" % (job_url)
        self.job_data_set(job_id, "analysis_base_url", analysis_base_url)

        analysis_url = "%s/ANALYSIS/index.html" % (job_url)
        self.job_data_set(job_id, "analysis_url", analysis_url)

        ## submission time and initial state
        self.job_data_set(job_id, "state", "submit1")
        self.job_data_set(job_id, "submit_time", time.time())

        ## now load the structure and build the submission form
        try:
            struct = FileIO.LoadStructure(fil = pdb_filename)
        except:
            return "The Python Macromolecular Library was unable to load your structure file."
            
	if not struct.structure_id:
	    struct.structure_id = "XXXX"
        self.job_data_set(job_id, "structure_id", struct.structure_id)

        ## Select Chains for Analysis
        num_atoms          = 0
        num_aniso_atoms    = 0
        largest_chain_seen = 0

        chains = []
        for chain in struct.iter_chains():
            naa = chain.count_amino_acids()
            if naa < 10:
                continue

            largest_chain_seen = max(naa, largest_chain_seen)

            for atm in chain.iter_all_atoms():
                num_atoms += 1
                if atm.U != None:
                    num_aniso_atoms += 1

            ## form name
            cb_name = 'CHAIN%s' % (chain.chain_id)
            
            ## create chain description label cb_desc
            cb_desc = 'Chain %s (%d Amino Acid Residues)' % (chain.chain_id, chain.count_amino_acids())

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
            cdict["length"] = naa
            cdict["name"] = cb_name
            cdict["desc"] = cb_desc
            cdict["preview"] = cb_preview
            cdict["selected"] = True

        if largest_chain_seen > 1700:
            self.remove_job(job_id)
	    return 'Your submitted structure contained a chain exceeding the 1700 residue limit'

        self.job_data_set(job_id, "chains", chains)

        ## defaults
        self.job_data_set(job_id, "user", "")
        self.job_data_set(job_id, "passwd", "")
        self.job_data_set(job_id, "email", "")
        self.job_data_set(job_id, "comment", "")
        self.job_data_set(job_id, "private_job", False)
        self.job_data_set(job_id, "plot_format", "PNG")

        aniso_ratio = float(num_aniso_atoms) / float(num_atoms)
        if aniso_ratio > 0.90:
            self.job_data_set(job_id, "tls_model", "ANISO")
        else:
            self.job_data_set(job_id, "tls_model", "ISOT")
            
        self.job_data_set(job_id, "weight", "")
        self.job_data_set(job_id, "include_atoms", "ALL")

	return ""

    def remove_job(self, job_id):
        """Removes the job from both the database and working directory.
        """
        if not self.job_exists(job_id):
            return False
        
        job_dir = self.job_data_get(job_id, "job_dir")

        if job_dir and job_dir.startswith(conf.TLSMD_WORK_DIR) and os.path.isdir(job_dir):

            for root, dirs, files in os.walk(job_dir, topdown = False):
                for name in files:
                    os.remove(os.path.join(root, name))
                for name in dirs:
                    os.rmdir(os.path.join(root, name))

            os.rmdir(job_dir)

        self.job_delete(job_id)
        return True
        
    def run_server(self, host, port):
        xmlrpc_server = WebTLSMD_XMLRPCServer(
            (host, port),
            SimpleXMLRPCServer.SimpleXMLRPCRequestHandler,
            True)

        xmlrpc_server.register_function(self.job_list,               "job_list")
        xmlrpc_server.register_function(self.job_new,                "job_new")
        xmlrpc_server.register_function(self.job_exists,             "job_exists")
        xmlrpc_server.register_function(self.job_get_dict,           "job_get_dict")
        xmlrpc_server.register_function(self.job_get_dict_index,     "job_get_dict_index")
        xmlrpc_server.register_function(self.job_delete,             "job_delete")
        xmlrpc_server.register_function(self.job_data_set,           "job_data_set")
        xmlrpc_server.register_function(self.job_data_get,           "job_data_get")
        xmlrpc_server.register_function(self.get_next_queued_job_id, "get_next_queued_job_id")
        xmlrpc_server.register_function(self.set_structure_file,     "set_structure_file")
        xmlrpc_server.register_function(self.remove_job,             "remove_job")
        
        xmlrpc_server.serve_forever()


def main():
    rtype, baseurl, port = conf.WEBTLSMDD.split(":")

    print "webtlsmdd.py xmlrpc server version %s" % (const.VERSION)
    print "using database file %s" % (conf.WEBTLSMDD_DATABASE)
    print "listening for incoming connections at %s" % (conf.WEBTLSMDD)
    
    webtlsmdd = WebTLSMDDaemon2(conf.WEBTLSMDD_DATABASE)    
    webtlsmdd.run_server("localhost", int(port))


def inspect():
    database_path = sys.argv[2]
    
    webtlsmdd = WebTLSMDDaemon2(database_path)

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
    if len(sys.argv)==1:
        try:
            main()
        except KeyboardInterrupt:
            pass
    else:
        inspect()
    sys.exit(0)
