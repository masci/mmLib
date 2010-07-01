## TLS Minimized Domains (TLSMD)
## Copyright 2002-2010 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

## Python modules
import os
import sys
import time
import socket
import xmlrpclib
import cgitb; cgitb.enable()
import cgi

## TLSMD
import captions
import conf
import const
import misc
import mysql_support

## GLOBALS
webtlsmdd = xmlrpclib.ServerProxy(conf.WEBTLSMDD)
mysql = mysql_support.MySQLConnect()


class Page(object):
    """Formats HTML in the refinement preparation page.
    """

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
        x += '  <link rel="stylesheet" href="%s/%s" type="text/css" media="screen">' % (
                 conf.TLSMD_PUBLIC_URL, "tlsmd.css")
        x += '</head>\n'
        x += '<body>'
        x += '<div id="page">'
        return x

    def html_head(self, title):
        x  = ''
        x += 'Content-Type: text/html\n\n'
        x += self.html_head_nocgi(title)
        return x

    def html_foot(self):
        x  = ''
        x += '<center>'
        x += '<p><small><b>Version %s</b> Released %s' % (
             const.VERSION, const.RELEASE_DATE)
        x += ' by %s ' % (const.AUTHOR)
        x += '<i>%s</i></small></p>' % (const.EMAIL)
        x += '</center>'
        x += '</div>'
        x += '</body></html>'
        return x


class ErrorPage(Page):
    """Returns a description of any errors generated by user inputs.
    """
    def __init__(self, form, text):
        Page.__init__(self, form)
        self.text = text

    def html_page(self):
        title = 'TLSMD: Error'

        x  = ''
        x += self.html_head(title)
        x += self.html_title(title)
        x += '<br/>'
        x += '<center><h3>An Error Occured</h3></center>'

        if self.text != None:
            x += '<center>'
            x += self.text
            x += '</center>'

        x += self.html_foot()
        return x


class RefinePrepError(Exception):
    """Returns a description of any refinement preparation errors.
    """
    def __init__(self, text):
        Exception.__init__(self)
        self.text = text


class RefinePrepPage(Page):
    """Main class to generate the necessary files needed for a refinement of a
    given structure and it TLS partitions.
    """

    def __init__(self, form):
        self.form = form

    def html_page(self):
        """Creates the input files for needed for TLS Refinement.
        """

        ## initialize values (some are just dummy values)
        job_id = "TLSMD0000_xxxxxxxx"
        struct_id = "XXXX"
        chain_ntls = []
        wilson = float(conf.MAX_WILSON_B)

        ## now fill in the initialized values with those in the web form.
        for key in self.form.keys():
            if key.startswith("job_id"):
                job_id = self.form[key].value

            elif key.startswith("struct_id"):
                struct_id = self.form[key].value

            ## extract ntls selections from CGI form
            elif key.startswith("NTLS_CHAIN"):
                chain_id = key[-1]
                try:
                    ntls = int(self.form[key].value)
                except ValueError:
                    continue
                chain_ntls.append((chain_id, ntls))

            ## Extract the "Constant B for pure TLS model (e.g., Wilson B)"
            ## value from the refinement prep form
            elif key.startswith("wilson"):
                try:
                    if misc.is_float(self.form[key].value):
                        wilson = float(self.form[key].value)
                        if wilson > conf.MAX_WILSON_B:
                            wilson = conf.MAX_WILSON_B
                    else:
                        msg  = "Value for constant B '%s' is not valid. " % (
                            self.form[key].value)
                        msg += "Please enter a floating point number."
                        raise RefinePrepError(msg)
                except ValueError:
                    continue

        chain_ntls.sort()

        ## make sure there were selections
        if len(chain_ntls) == 0:
            raise RefinePrepError("Form Processing Error: No Chains Selected")

        ## call webtlsmdd to generate files (Refmac5 + PHENIX)
        result = webtlsmdd.refmac5_refinement_prep(job_id, struct_id, 
                                                   chain_ntls, wilson)
        if isinstance(result, str):
            raise RefinePrepError(result)

        ## Success! Now provide a description and make download links for the
        ## input files needed for refinement.
        title = 'Input Files for TLS Refinement'

        x  = ''
        x += self.html_head(title)
        x += self.html_title(title)

        x += '<center>\n'
        x += '<h3>'
        x += 'Step 2: Download the generated XYZIN(PDBIN), TLSIN, and PHENIX files below'
        x += '</h3>'
        x += '</center>'

        x += '<p>%s</p>' % (captions.REFINEMENT_FILES_DOWNLOAD_INFO)

        x += '<center>\n'
        x += '<table class="submit_table">'

        ## "REFMAC (TLS + Biso) files"
        x += '<th colspan="2">REFMAC (TLS + Biso) files</th>'
        x += '<tr>'
        x += '<td align="right"><b>PDBIN File: </b></td>'
        x += '<td><a href="%s" type="text/plain">%s</a></td>' % (
            result["pdbout_url1"], result["pdbout1"].split("/")[1])
        x += '</tr><tr>'
        x += '<td align="right"><b>TLSIN File: </b></td>'
        x += '<td><a href="%s" type="text/plain">%s</a></td>' % (
            result["tlsout_url1"], result["tlsout1"].split("/")[1])
        x += '</tr>'

        ## "REFMAC (Pure TLS) files"
        x += '<tr><td colspan="2"><hr></td></tr>'
        x += '<th colspan="2">REFMAC (Pure TLS) files</th>'
        x += '<tr>'
        x += '<td align="right"><b>PDBIN File: </b></td>'
        x += '<td><a href="%s" type="text/plain">%s</a></td>' % (
            result["pdbout_url2"], result["pdbout2"].split("/")[1])
        x += '</tr><tr>'
        x += '<td align="right"><b>TLSIN File: </b></td>'
        x += '<td><a href="%s" type="text/plain">%s</a></td>' % (
            result["tlsout_url2"], result["tlsout2"].split("/")[1])
        x += '</tr>'

        ## "PHENIX files"
        x += '<tr><td colspan="2"><hr></td></tr>'
        x += '<th colspan="2">PHENIX files</th>'
        x += '<tr>'
        x += '<td align="right"><b>PHENIX File: </b></td>'
        x += '<td><a href="%s" type="text/plain">%s</a></td>' % (
            result["phenix_url"], result["phenix"].split("/")[1])
        x += '</tr>'
        x += '</table>'

        x += '<br/>'

        ## Documentation link
        x += '<center>\n'
        x += '<h3>'
        x += 'Step 3: Read this '
        x += '<a href="/~tlsmd/documentation.html#refmac5">How-To</a>'
        x += '</h3>'
        x += '</center>'

        x += self.html_foot()

        return x


def check_job_id(form):
    """Retrieves and confirms the job_id from a incoming form.
    Returns None on error, or the job_id on success.
    """
    if form.has_key("job_id"):
        job_id = form["job_id"].value
        if len(job_id) < 20:
            if job_id.startswith("TLSMD"):
                if mysql.job_exists(job_id):
                    return job_id
            elif mysql.pdb_exists(job_id):
                return job_id
    return None


def main():
    form = cgi.FieldStorage()

    page = RefinePrepPage(form)
    try:
        print page.html_page()

    except RefinePrepError, err:
        text = '<p>%s</p>' % (err.text)
        page = ErrorPage(form, text)
        print page.html_page()

    ## XXX: Might be used for debugging.
    #except xmlrpclib.Fault, fault:
    #    fault_html  = "xmlrpclib.Fault from refineprep.py:<br/>"
    #    fault_html += "fault code: %s<br/>fault string: '%s'" % (
    #        fault.faultCode, fault.faultString)
    #    page = ErrorPage(form, fault_html)
    #    print page.html_page()

    except socket.error, err:
        page = ErrorPage(form, "socket.error: " + str(err))
        print page.html_page()


if __name__=="__main__":
    main()
    sys.exit(0)
