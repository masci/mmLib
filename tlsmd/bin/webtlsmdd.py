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
PORT = 20100
FILE = os.path.join(os.environ["TLSMD_ROOT"], "data", "webtlsmdd.db")

###############################################################################

class WebTLSMDData(object):
    def __init__(self):
        self.job_num  = 1
        self.job_list = []

    def get_jdict(self, job_id):
        for jdict in self.job_list:
            if job_id==jdict["job_id"]:
                return jdict
        return None


class WebTLSMDDaemon(object):
    def save(self, data):
        cPickle.dump(data, open(FILE, "wb"))

    def load(self):
        if os.path.exists(FILE):
            return cPickle.load(open(FILE, "rb"))
        else:
            return WebTLSMDData()

    def job_list(self):
        data = self.load()
        return data.job_list

    def job_new(self):
        data = self.load()
        job_id = "TLSMD%d" % (data.job_num)
        data.job_num += 1
        
        jdict = {}
        jdict["job_id"] = job_id
        data.job_list.append(jdict)
        self.save(data)

        return job_id

    def job_exists(self, job_id):
        data = self.load()
        jdict = data.get_jdict(job_id)
        if jdict==None:
            return False
        return True

    def job_get_dict(self, job_id):
        data = self.load()
        jdict = data.get_jdict(job_id)
        if jdict==None:
            return False
        return jdict
    
    def job_get_dict_index(self, i):
        data = self.load()
        try:
            jdict = data.job_list[i]
        except IndexError:
            return False
        return jdict

    def job_delete(self, job_id):
        data = self.load()
        jdict = data.get_jdict(job_id)
        if jdict==None:
            return False
        data.job_list.remove(jdict)
        self.save(data)
        return True

    def job_data_set(self, job_id, key, value):
        data = self.load()
        jdict = data.get_jdict(job_id)
        if jdict==None:
            return False
        jdict[key] = value
        self.save(data)
        return True

    def job_data_get(self, job_id, key):
        data = self.load()
        jdict = data.get_jdict(job_id)
        if jdict==None:
            return False
        return jdict.get(key, False)

    def run_server(self, host, port):
        xmlrpc_server = SimpleXMLRPCServer.SimpleXMLRPCServer(
            (host, port), SimpleXMLRPCServer.SimpleXMLRPCRequestHandler, False)

        xmlrpc_server.register_function(self.job_list,     "job_list")
        xmlrpc_server.register_function(self.job_new,      "job_new")
        xmlrpc_server.register_function(self.job_exists,   "job_exists")
        xmlrpc_server.register_function(self.job_get_dict, "job_get_dict")
        xmlrpc_server.register_function(
            self.job_get_dict_index, "job_get_dict_index")
        xmlrpc_server.register_function(self.job_delete,   "job_delete")
        xmlrpc_server.register_function(self.job_data_set, "job_data_set")
        xmlrpc_server.register_function(self.job_data_get, "job_data_get")

        xmlrpc_server.serve_forever()


class WebTLSMDDaemon2(object):
    def __init__(self, db_file):
        self.db_file = db_file
        self.db = bsddb.hashopen(self.db_file, "c")

    def store_dict(self, dbkey, dictx):
        self.db[dbkey] = cPickle.dumps(dictx)
        self.grh_db.sync()
    
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
        gdict = retrieve_dict("GLOBALS")

        ## no global dictionary; create one which is consistant with
        ## the TLSMD jobs in the database
        if gdict==None:
            gdict = {}

            job_list = self.job_list()
            if len(job_list)==0:
                gdict["job_num"] = 1
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
                jdict = self.retrieve_jdict(db_key)
                job_num = int(dbkey[5:])
                listx.append((job_num, jdict))

        listx.sort()

        job_list = []
        for job_num, jdict in listx:
            job_list.append(jdict)

        return job_list

    def job_new(self):
        gdict = retrieve_globals()

        ## assign job_id
        job_id = "TLSMD%d" % (gdict["next_job_num"])
        gdict["next_job_num"] += 1

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
        self.grh_db.sync()
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

    def run_server(self, host, port):
        xmlrpc_server = SimpleXMLRPCServer.SimpleXMLRPCServer(
            (host, port), SimpleXMLRPCServer.SimpleXMLRPCRequestHandler, False)

        xmlrpc_server.register_function(self.job_list,     "job_list")
        xmlrpc_server.register_function(self.job_new,      "job_new")
        xmlrpc_server.register_function(self.job_exists,   "job_exists")
        xmlrpc_server.register_function(self.job_get_dict, "job_get_dict")
        xmlrpc_server.register_function(
            self.job_get_dict_index, "job_get_dict_index")
        xmlrpc_server.register_function(self.job_delete,   "job_delete")
        xmlrpc_server.register_function(self.job_data_set, "job_data_set")
        xmlrpc_server.register_function(self.job_data_get, "job_data_get")

        xmlrpc_server.serve_forever()



def main():
    webtlsmdd = WebTLSMDDaemon()
    webtlsmdd.run_server(HOST, PORT)

def convert():
    wd1 = WebTLSMDDaemon()
    wd2 = WebTLSMDDaemon("convert.db")

    job_list = wd1.job_list()

    for jdict in job_list:
        wd2.store_jdict(jdict)


if __name__=="__main__":
    try:
        main()
    except KeyboardInterrupt:
        pass

    sys.exit(0)
