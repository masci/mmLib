## TLS Motion Determination (TLSMD)
## Copyright 2002-2008 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys
import subprocess
import traceback

import conf

def SendEmail(address, subject, body):
    if not os.path.isfile(conf.MAIL):
        sys.stderr.write("mail client not found: %s" % (conf.MAIL))
        return
    
    ## send mail using /usr/bin/mail
    try:
        pobj = subprocess.Popen([conf.MAIL, "-s", subject, address],
                                stdin = subprocess.PIPE,
                                stdout = subprocess.PIPE,
                                stderr = subprocess.STDOUT,
                                close_fds = True,
                                bufsize = 8192)
    except OSError:
        sys.stderr.write("[ERROR] mail client failed to execute: %s" % (conf.MAIL))
        return
        
    pobj.stdin.write(body)
    pobj.stdin.close()
    pobj.wait()
    

def SendTracebackEmail(context):
    SendEmail(conf.TRACEBACK_EMAIL, context, traceback.format_exc())
