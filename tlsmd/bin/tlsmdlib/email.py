## TLS Motion Determination (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys
import popen2
import traceback

import conf

def SendEmail(address, subject, body):
    if not os.path.isfile(conf.MSMTP):
        sys.stderr.write("Mail Client %s Not Found" % (conf.MSMTP))
        return
    
    mlist = ["To: %s" % (address),
             "Subject: %s" % (subject),
             "",
             body]

    ## send mail using msmtp
    pobj = popen2.Popen4([conf.MSMTP, address])
    pobj.tochild.write("\n".join(mlist))
    pobj.wait()
    

def SendTracebackEmail(context):
    SendEmail(conf.TRACEBACK_EMAIL, context, traceback.format_exc())
