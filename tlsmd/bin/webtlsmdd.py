#!/home/tlsmd/local/bin/python
## TLS Minimized Domains (TLSMD)
## Copyright 2002 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys
import time

import cPickle
import xmlrpclib
import SimpleXMLRPCServer

###############################################################################
## CONFIGURATION

HOST = "localhost"
PORT = 10100
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

def main():
    webtlsmdd = WebTLSMDDaemon()
    webtlsmdd.run_server(HOST, PORT)

if __name__=="__main__":
    try:
        main()
    except KeyboardInterrupt:
        pass

    sys.exit(0)
