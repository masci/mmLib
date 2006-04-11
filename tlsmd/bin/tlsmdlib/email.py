## TLS Motion Determination (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys
import subprocess
import traceback

import conf

def SendEmail(address, subject, body):
    if not os.path.isfile(conf.MSMTP):
        sys.stderr.write("mail client not found: %s" % (conf.MSMTP))
        return
    
    mlist = ["To: %s" % (address),
             "Subject: %s" % (subject),
             "",
             body]

    ## send mail using msmtp
    try:
        pobj = subprocess.Popen([conf.MSMTP, address],
                                stdin = subprocess.PIPE,
                                stdout = subprocess.PIPE,
                                stderr = subprocess.STDOUT,
                                close_fds = True,
                                bufsize = 8192)
    except OSError:
        sys.stderr.write("[ERROR] mail client failed to execute: %s" % (conf.MSMTP))
        return
        
    pobj.stdin.write("\n".join(mlist))
    pobj.stdin.close()
    pobj.wait()
    

def SendTracebackEmail(context):
    SendEmail(conf.TRACEBACK_EMAIL, context, traceback.format_exc())
