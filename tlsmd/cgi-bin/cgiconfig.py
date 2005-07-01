## TLS Minimized Domains (TLSMD)
## Copyright 2002 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import re
import string

WEBTLSMDD              = "http://localhost:10100"
VERSION                = "0.0.1"
LAST_MODIFIED_BY       = "Jay Painter"
LAST_MODIFIED_BY_EMAIL = "jpaint@u.washington.edu"
LAST_MODIFIED_DATE     = "Jume 11, 2005"

TLSMD_WORK_DIR         = "/home/jpaint/public_html/webtlsmd/run"
TLSMD_WORK_URL         = "/~jpaint/webtlsmd/run"

LINK_SPACE             = '&nbsp;&nbsp;&nbsp;&nbsp;'

DOCUMENTATION_URL      = "/~jpaint/cgi-bin/documentation.html"
DOCUMENTATION_PATH     = "/home/jpaint/tlsmd/doc/documentation.html"


def get_documentation_block(block_name):
    reblock = re.compile("\s*<!--\s*BLOCK:\s*%s\s*-->.*" % (block_name))

    fil = open(DOCUMENTATION_PATH, "r")
    inblock = False
    listx = []

    for ln in fil.readlines():
        if inblock==False:
            if reblock.match(ln):
                inblock = True
            continue
        if ln.count("<!-- ENDBLOCK -->")>0:
            break
        listx.append(ln)

    return string.join(listx, "")
            
        
