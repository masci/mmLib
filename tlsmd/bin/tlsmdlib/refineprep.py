## TLS Minimized Domains (TLSMD)
## Copyright 2002-2008 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys
import time
import socket
import xmlrpclib
import cgitb; cgitb.enable()
import cgi

import const
import conf

## GLOBALS
webtlsmdd = xmlrpclib.ServerProxy(conf.WEBTLSMDD)

CAPTION = """\
<b>Refmac5:</b> Download both the modified PDBIN file for your structure and the corresponding
TLSIN file. Feed these to REFMAC5 as a starting point for multi-TLS group refinement.
See the TLSMD documentation for detailed instructions.
<p>
<b>PHENIX:</b> The PHENIX file contains a description of the TLS groups you selected.
This file is intended to be read by the PHENIX.refine input scripts. 
"""

class Page(object):
    def __init__(self, form):
        self.form = form

    def html_title(self, title):
        return '<center><h1>%s</h1></center>' % (title)

    def html_head_nocgi(self, title):
        x  = ''
        x += '<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" '
        x += '"http://www.w3.org/TR/html4/loose.dtd">\n\n'
        x += '<html>'
        x += '<head>'
        x += '  <title>%s</title>' % (title)
        x += '  <link rel="stylesheet" href="../tlsmd.css" type="text/css" media="screen">'
        x += '</head>'
        x += '<body>'
        x += '<div id="page">'
        return x

    def html_head(self, title):
        x = ''
        x += 'Content-Type: text/html\n\n'
        x += self.html_head_nocgi(title)
        return x

    def html_foot(self):
        x = ''
        x += '<center>'
        x += '<p><small><b>Version %s</b> Released %s' % (const.VERSION, const.RELEASE_DATE)
        x += ' by %s ' % (const.AUTHOR)
        x += '<i>%s</i></small></p>' % (const.EMAIL)
        x += '</center>'
        x += '</div>'
        x += '</body></html>'
        return x


class ErrorPage(Page):
    def __init__(self, form, text):
        Page.__init__(self, form)
        self.text = text
    
    def html_page(self):
        title = 'TLSMD: Error'
        
        x  = ''
        x += self.html_head(title)
        x += self.html_title(title)
        x += '<br>'
        x += '<center><h3>An Error Occured</h3></center>'

        if self.text!=None:
            x += self.text
            
        x += self.html_foot()
        return x


class RefinePrepError(Exception):
    def __init__(self, text):
        Exception.__init__(self)
        self.text = text


class RefinePrepPage(Page):
    def html_page(self):
        
        job_id = check_job_id(self.form)
    
        title = 'Input Files for TLS Refinement'
        
        x  = ''
        x += self.html_head(title)
        x += self.html_title(title)

        x += '<center>'
        x += '<h3>'
        x += 'Step 2: Download the generated XYZIN(PDBIN), TLSIN, and PHENIX files below'
        x += '</h3>'
        x += '</center>'
        
        ## extract ntls selections from CGI form
        chain_ntls = []
        for key in self.form.keys():
            if key.startswith("NTLS_CHAIN"):
                chain_id = key[-1]
                try:
                    ntls = int(self.form[key].value)
                except ValueError:
                    continue
                chain_ntls.append((chain_id, ntls))

        chain_ntls.sort()

        ## make sure there were selections
        if len(chain_ntls) == 0:
            raise RefinePrepError("Form Processing Error: No Chains Selected")

        # EAM DEBUG 1
        #fault_html = "DEBUG in RefinePrepPage at #1: job_id = %s chain_ntls = %s" % (job_id, chain_ntls)
        #dpage = ErrorPage(self.form, fault_html)
        #print dpage.html_page()
        # EAM DEBUG

        ## call webtlsmdd to generate files
        result = webtlsmdd.refmac5_refinement_prep(job_id, chain_ntls)
        if isinstance(result, str):
            raise RefinePrepError(result)

        x += '<p>%s</p>' % (CAPTION)

        # EAM DEBUG 2
        #fault_html = "DEBUG in RefinePrepPage at #2: job_id = %s" % (job_id)
        #dpage = ErrorPage(self.form, fault_html)
        #print dpage.html_page()
        # EAM DEBUG

        ## success -- make download links
        x += '<center>'
        x += '<table border="0" style="background-color:#eeeeee">'
        x += '<tr>'
        x += '<td align="right"><b>PDBIN File</b></td>'
        x += '<td><a href="%s" type="text/plain">%s</a></td>' % (result["pdbout_url"], result["pdbout"])
        x += '</tr><tr>'
        x += '<td align="right"><b>TLSIN File</b></td>'
        x += '<td><a href="%s" type="text/plain">%s</a></td>' % (result["tlsout_url"], result["tlsout"])
        x += '</tr><tr>'
        x += '<td align="right"><b>PHENIX File</b></td>'
        x += '<td><a href="%s" type="text/plain">%s</a></td>' % (result["phenixout_url"], result["phenixout"])
        x += '</table>'

        x += '<br>'

        x += '<center>'
        x += '<h3>'
        x += 'Step 3: Read this '
        x += '<a href="/~tlsmd/documentation.html#refmac5">How-To</a>'
        x += '</h3>'
        x += '</center>'

        x += self.html_foot()
        return x


def check_job_id(form):
    """Retrieves and confirms the job_id from a incomming form.  Returns
    None on error, or the job_id on success.
    """
    if form.has_key("job_id"):
        job_id = form["job_id"].value
        if len(job_id) < 20:
            if job_id.startswith("TLSMD"):
                if webtlsmdd.job_exists(job_id):
                    return job_id
    return None


def main():
    form = cgi.FieldStorage()

    page = None
    job_id = check_job_id(form)
    if job_id == None:
        page = ErrorPage(form, "<center>The Job ID seems to be expired or invalid.</center>")
    else:
        page = RefinePrepPage(form)

    try:
        print page.html_page()

    except RefinePrepError, err:
        text = '<center><p>%s</p></center>' % (err.text)
        page = ErrorPage(form, text)
        print page.html_page()

    except xmlrpclib.Fault, fault:
        fault_html = "xmlrpclib.Fault from refineprep.py:<br>fault code: %s<br>fault string: %s" % (fault.faultCode, fault.faultString)
        page = ErrorPage(form, fault_html)
        print page.html_page()

    except socket.error, err:
        page = ErrorPage(form, "socket.error: " + str(err))
        print page.html_page()


if __name__=="__main__":
    main()
    sys.exit(0)
