#!/home/jpaint/local/bin/python
## TLS Minimized Domains (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys
import time

import cPickle
import bsddb

import xmlrpclib
import SimpleXMLRPCServer

###############################################################################
## CONFIGURATION

HOST = "localhost"
PORT = 10100

###############################################################################

class WebTLSMDDaemon2(object):
    def __init__(self, db_file):
        self.db_file = db_file
        self.db = bsddb.hashopen(self.db_file, "c")

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

        ## no global dictionary; create one which is consistant with
        ## the TLSMD jobs in the database
        if gdict==None:
            gdict = {}

            job_list = self.job_list()
            if len(job_list)==0:
                gdict["next_job_num"] = 1
            else:
                jdict = job_list[-1]
                job_id = jdict["job_id"]
                job_num = int(job_id[5:])
                gdict["next_job_num"] = job_num + 1

            self.store_globals(gdict)

        return gdict

    def store_globals(self, gdict):
        self.store_dict("GLOBALS", gdict)

    def store_jdict(self, jdict):
        dbkey = jdict["job_id"]
        self.store_dict(jdict["job_id"], jdict)
    
    def retrieve_jdict(self, job_id):
        return self.retrieve_dict(job_id)

    def job_list(self):
        """Returns a ordered list of all jdicts in the database
        """
        ## retrieve all jdicts from database and 
        listx = []
        for dbkey in self.db.keys():
            if dbkey.startswith("TLSMD"):
                jdict = self.retrieve_jdict(dbkey)
                job_num = int(dbkey[5:])
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
        job_id = "TLSMD%d" % (job_num)

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
        if jdict==None:
            return False
        jdict[key] = value
        self.store_jdict(jdict)
        return True

    def job_data_get(self, job_id, key):
        jdict = self.retrieve_jdict(job_id)
        if jdict==None:
            return False
        return jdict.get(key, False)

    def get_next_queued_job_id(self):
        job_list = self.job_list()
	for jdict in job_list:
            if jdict.get("state")=="running" or jdict.get("state")=="queued":
                return jdict["job_id"]
	return False

    def run_server(self, host, port):
        xmlrpc_server = SimpleXMLRPCServer.SimpleXMLRPCServer(
            (host, port), SimpleXMLRPCServer.SimpleXMLRPCRequestHandler, True)

        xmlrpc_server.register_function(self.job_list,     "job_list")
        xmlrpc_server.register_function(self.job_new,      "job_new")
        xmlrpc_server.register_function(self.job_exists,   "job_exists")
        xmlrpc_server.register_function(self.job_get_dict, "job_get_dict")
        xmlrpc_server.register_function(
            self.job_get_dict_index, "job_get_dict_index")
        xmlrpc_server.register_function(self.job_delete,   "job_delete")
        xmlrpc_server.register_function(self.job_data_set, "job_data_set")
        xmlrpc_server.register_function(self.job_data_get, "job_data_get")
        xmlrpc_server.register_function(self.get_next_queued_job_id, "get_next_queued_job_id")

        xmlrpc_server.serve_forever()

def main():
    database_file = os.environ["TLSMD_DATABASE"]
    webtlsmdd = WebTLSMDDaemon2(database_file)
    webtlsmdd.run_server(HOST, PORT)

def inspect():
    database_file = os.environ["TLSMD_DATABASE"]
    webtlsmdd = WebTLSMDDaemon2(database_file)

    if sys.argv[1]=="list":
        for dbkey in webtlsmdd.db.keys():
            print dbkey

    if sys.argv[1]=="remove":
        dbkey = sys.argv[2]
	del webtlsmdd.db[dbkey]
        webtlsmdd.db.sync()

if __name__=="__main__":

    if len(sys.argv)==1:
        try:
            main()
        except KeyboardInterrupt:
            pass
    
    else:
        inspect()

    sys.exit(0)
