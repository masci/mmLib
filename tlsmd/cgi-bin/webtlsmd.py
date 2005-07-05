## TLS Minimized Domains (TLSMD)
## Copyright 2002 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys
import time
import socket
import string

import xmlrpclib
import cgitb; cgitb.enable()
import cgi

from mmLib.Structure  import *
from mmLib.FileLoader import *


## GLOBALS
from cgiconfig import *
webtlsmdd = xmlrpclib.ServerProxy(WEBTLSMDD, allow_none=True)

def timestring(secs):
    tm_struct = time.localtime(secs)
    return time.strftime("%m-%d-%y %H:%M %Z" ,tm_struct)


def timediffstring(begin, end):
    secs = end - begin
    hours = secs / float(3600)
    return "%.2f" % (hours)


def html_title(title):
    """Title
    """
    x  = ''
    x += '<center><h1>%s</h1></center>' % (title)
    return x


def html_nav_bar(page_name=None):
    """Site navigation bar.
    """
    x  = ''
    x += '<center>'

    if page_name=="index":
        x += 'Home'
    else:
        x += '<a href="webtlsmd.cgi">Home</a>'

    x += LINK_SPACE

    if page_name=="queue":
        x += 'Job Queue'
    else:
        x += '<a href="webtlsmd.cgi?page=queue">Job Queue</a>'

    x += LINK_SPACE

    if page_name=="edit":
        x += 'Edit Job'
    else:
        x += '<a href="webtlsmd.cgi?page=edit">Edit Job</a>'

    x += LINK_SPACE

    if page_name=="completed":
        x += 'Completed Jobs'
    else:
        x += '<a href="webtlsmd.cgi?page=completed">Completed Jobs</a>'

    x += LINK_SPACE

    x += '<a href="../webtlsmd/examples">Examples</a>'
    x += LINK_SPACE
    x += '<a href="../webtlsmd/doc/documentation.html">Documentation</a>'
    x += '</center>'
    x += '<br>'
    return x


def html_job_nav_bar(webtlsmdd, job_id):
    """Navigation bar to the TLSMD output files.
    """
    analysis_dir = webtlsmdd.job_data_get(job_id, "analysis_dir")
    analysis_index = os.path.join(analysis_dir, "index.html")
    analysis_url = webtlsmdd.job_data_get(job_id, "analysis_url")

    job_dir = webtlsmdd.job_data_get(job_id, "job_dir")
    logfile = os.path.joing(job_dir, "log.txt")
    log_url = webtlsmdd.job_data_get(job_id, "log_url")

    x  = ''
    x += '<center>'
    x += '<h3>'

    if os.path.isfile(analysis_index):
        x += '<a href="%s">Completed Analysis</a>' % (analysis_url)

    x += LINK_SPACE

    if os.path.isfile(logfile):
        x += '<a href="%s">Logfile</a>' % (log_url)

    x += '</h3>'
    x += '</center>'
    x += '<br>'
    return x


def html_job_edit_form(fdict):
    x  = ''
    x += '<center>'

    x += '<form '\
         'enctype="multipart/form-data" '\
         'action="webtlsmd.cgi" '\
         'method="get">'

    x += '<input type="hidden" name="page" value="%s">' % (
        fdict.get("page", "index"))
    x += '<input type="hidden" name="edit_form" value="TRUE">'
    x += '<input type="hidden" name="job_id" value="%s">' % (
        fdict["job_id"])

    x += '<table border="1" width="100%">'

    ## user/email/passcode/structure name
    x += '<tr>'
    x += '<th colspan="2">User Information</th>'
    x += '<th>Session Information</th>'
    x += '</tr>'

    x += '<tr><td colspan="2">'
    x += '<table>'

    ## user
    x += '<tr>'
    x += '<td align="right"><label>User name:</td><td>'

    user = fdict.get("user", "")
    if len(user)>0:
        x += '<input type="hidden" name="user" value="%s">' % (user)
        x += '<b>%s</b>' % (user)
    else:
        x += '<input '\
             'type="text" '\
             'name="user" '\
             'value="%s" '\
             'size="10" '\
             'maxlength="10">' % (user)
    x += '</label></td>'
    x += '</tr>'

    ## password
    x += '<tr>'
    x += '<td align="right"><label>Password:</td><td>'
    passwd = fdict.get("passwd", "")
    if len(passwd)>0:
        x += '<input type="hidden" name="passwd" value="%s">' % (passwd)
        x += '<b>%s</b>' % (passwd)
    else:
        x += '<input '\
             'type="text" '\
             'name="passwd" '\
             'value="%s" '\
             'size="10" '\
             'maxlength="10">' % (passwd)
    x += '</label></td>'
    x += '</tr>'

    ## email address
    x += '<tr>'
    x += '<td align="right"><label>EMail Address:</td><td>'
    x += '<input '\
         'type="text" '\
         'name="email" '\
         'value="%s" '\
         'size="25" '\
         'maxlength="40">' % (fdict.get("email", ""))
    x += '</label></td>'
    x += '</tr>'

    ## structure code
    x += '<tr>'
    x += '<td align="right"><label>Structure Code:</td><td>'
    x += '<input '\
         'type="text" '\
         'name="structure_id" '\
         'value="%s" '\
         'size="4" '\
         'maxlength="4">' % (fdict.get("structure_id", ""))
    x += '</label></td>'
    x += '</tr>'

    ## comment
    x += '<tr>'
    x += '<td align="right"><label>Comment:</td><td>'
    x += '<input '\
         'type="text" '\
         'name="comment" '\
         'value="%s" '\
         'size="25" '\
         'maxlength="40">' % (fdict.get("comment", ""))
    x += '</label></td>'
    x += '</tr>'

    x += '</table>'
    x += '</td>'

    ## session info
    x += '<td valign="top"><table>'

    x += '<tr><td align="right">TLSMD Job ID:</td>'
    x += '<td><b>%s</b></td></tr>' % (fdict["job_id"])

    x += '<tr><td align="right">Job State:</td>'
    x += '<td><b>%s</b></td></tr>' % (fdict["state"])
    
    x += '<tr><td align="right">Submission IP Address: </td>'
    x += '<td><b>%s</b></td></tr>' % (fdict.get("ip_addr", ""))

    x += '<tr><td align="right">Submission Date: </td>'

    if fdict.has_key("submit_time"):
        date = timestring(fdict["submit_time"])
    else:
        date = "No Time"
    x += '<td><b>%s</b></td></tr>' % (date)

    x += '</table></td>'

    x += '</tr>'

    ## Select Chains for Analysis
    x += '<tr><th colspan="3">Select Chains for Analysis</th></tr>'

    x += '<tr><td colspan="3">'
    x += '<table>'
    for cdict in fdict["chains"]:
        x += '<tr><td>'
        x += '<label>'
        if cdict["selected"]==True:
            x += '<input type="checkbox" name="%s" value="TRUE" checked>' % (
                cdict["name"])
        else:
            x += '<input type="checkbox" name="%s" value="TRUE">' % (
                cdict["name"])
        x += '%s: <small>%s...</small>' % (cdict["desc"], cdict["preview"])
        x += '</label>'

        x += '</td></tr>'

    x += '</table></td></tr>'

    x += '<tr>'
    x += '<th>TLS Model</th>'\
         '<th>Least Squares Weighting</th>'\
         '<th>Include Atoms</th>'
    x += '</tr>'

    x += '<tr>'

    ## TLS Model
    x += '<td><table>'

    x += '<tr><td><label>'
    if fdict.get("tls_model")==None or fdict.get("tls_model")=="ISOT":
        x += '<input type="radio" name="tls_model" value="ISOT" checked>'
    else:
        x += '<input type="radio" name="tls_model" value="ISOT">'

    x += 'Isotropic'
    x += '</label></td></tr>'

    x += '<tr><td><label>'
    if fdict.get("tls_model")=="ANISO":
        x += '<input type="radio" name="tls_model" value="ANISO" checked>'
    else:
        x += '<input type="radio" name="tls_model" value="ANISO">'
    x += 'Anisotropic'
    x += '</label></td></tr>'

    x += '</table></td>'

    ## Least Squares Weighting
    x += '<td><table>'

    x += '<tr><td><label>'
    if fdict.get("weight")==None or fdict.get("weight")=="IUISO":
        x += '<input type="radio" name="weight" value="IUISO" checked>'
    else:
        x += '<input type="radio" name="weight" value="IUISO">'
    x += 'Inverse Atomic B<sub>iso</sub>'
    x += '</label></td></tr>'

    x += '<tr><td><label>'
    if fdict.get("weight")=="NONE":
        x += '<input type="radio" name="weight" value="NONE" checked>'
    else:
        x += '<input type="radio" name="weight" value="NONE">'
    x += 'No Weighting'
    x += '</label></td></tr>'

    x += '</table></td>'

    ## Include Atoms
    x += '<td><table>'

    x += '<tr><td><label>'
    if fdict.get("include_atoms")==None or fdict.get("include_atoms")=="ALL":
        x += '<input type="radio" name="include_atoms" value="ALL" checked>'
    else:
        x += '<input type="radio" name="include_atoms" value="ALL">'
    x += 'Include All Atoms'
    x += '</label></td></tr>'

    x += '<tr><td><label>'
    if fdict.get("include_atoms")=="MAINCHAIN":
        x += '<input type="radio" name="include_atoms" '\
             'value="MAINCHAIN" checked>'
    else:
        x += '<input type="radio" name="include_atoms" '\
             'value="MAINCHAIN">'
    x += 'Main Chain Atoms'
    x += '</label></td></tr>'

    x += '<tr><td><label>'
    if fdict.get("include_atoms")=="CA":
        x += '<input type="radio" name="include_atoms" '\
             'value="CA" checked>'
    else:
        x += '<input type="radio" name="include_atoms" '\
             'value="CA">'
    x += 'C-Alpha Atoms'
    x += '</label></td></tr>'

    x += '</table></td>'

    x += '</tr>'

    ## end form

    x += '<tr><td colspan="3">'

    x += '<table width="100%">'
    x += '<tr>'
    
    x += '<td align="left">'
    if fdict.has_key("removebutton"):
        x += '<input type="submit" name="submit" value="Remove Job">'
    x += '</td>'
    
    x += '<td align="right">'
    x += '<input type="submit" name="submit" value="OK">'
    x += LINK_SPACE
    x += '<input type="submit" name="submit" value="Cancel">'
    x += '</tr>'
    x += '</table>'

    x += '</td></tr>'
    x += '</table>'
    x += '</form>'
    return x


def html_job_info_table(fdict):
    x  = ''
    x += '<center>'

    x += '<table border="1" width="100%">'

    ## user/email/passcode/structure name
    x += '<tr>'
    x += '<th colspan="2">User Information</th>'
    x += '<th>Session Information</th>'
    x += '</tr>'

    x += '<tr><td colspan="2">'
    x += '<table>'

    ## user
    x += '<tr>'
    x += '<td align="right"><label>User name:</td><td>'
    x += '<b>%s</b>' % (fdict.get("user",""))
    x += '</label></td>'
    x += '</tr>'

    ## password
    x += '<tr>'
    x += '<td align="right"><label>Password:</td><td>'
    x += '<b>%s</b>' % (fdict.get("passwd", ""))
    x += '</label></td>'
    x += '</tr>'

    ## email address
    x += '<tr>'
    x += '<td align="right"><label>EMail Address:</td><td>'
    x += '<b>%s</b>' % (fdict.get("email", ""))
    x += '</label></td>'
    x += '</tr>'

    ## structure code
    x += '<tr>'
    x += '<td align="right"><label>Structure Code:</td><td>'
    x += '<b>%s</b>' % (fdict.get("structure_id", ""))
    x += '</label></td>'
    x += '</tr>'

    ## comment
    x += '<tr>'
    x += '<td align="right"><label>Comment:</td><td>'
    x += '<b>%s</b>' % (fdict.get("comment", ""))
    x += '</label></td>'
    x += '</tr>'

    x += '</table>'
    x += '</td>'

    ## session info
    x += '<td valign="top"><table>'

    x += '<tr><td align="right">TLSMD Job ID:</td>'
    x += '<td><b>%s</b></td></tr>' % (fdict["job_id"])

    x += '<tr><td align="right">Job State:</td>'
    x += '<td><b>%s</b></td></tr>' % (fdict["state"])
    
    x += '<tr><td align="right">Submission IP Address: </td>'
    x += '<td><b>%s</b></td></tr>' % (fdict.get("ip_addr", ""))

    x += '<tr><td align="right">Submission Date: </td>'

    if fdict.has_key("submit_time"):
        date = timestring(fdict["submit_time"])
    else:
        date = "No Time"
    x += '<td><b>%s</b></td></tr>' % (date)

    x += '</table></td>'

    x += '</tr>'

    ## Select Chains for Analysis
    x += '<tr><th colspan="3">Select Chains for Analysis</th></tr>'

    x += '<tr><td colspan="3">'
    x += '<table>'
    for cdict in fdict["chains"]:
        x += '<tr><td>'
        if cdict["selected"]==True:
            x += '<tr><td>'
            x += '%s: <small>%s...</small>' % (cdict["desc"], cdict["preview"])
            x += '</td></tr>'

    x += '</table></td></tr>'

    ## column titles
    x += '<tr>'
    x += '<th>TLS Model</th>'\
         '<th>Least Squares Weighting</th>'\
         '<th>Include Atoms</th>'
    x += '</tr>'

    x += '<tr>'

    ## TLS Model
    x += '<td>'
    if fdict.get("tls_model")==None or fdict.get("tls_model")=="ISOT":
        x += 'Isotropic'
    elif fdict.get("tls_model")=="ANISO":
        x += 'Anisotropic'
    x += '</td>'

    ## Least Squares Weighting
    x += '<td>'
    if fdict.get("weight")==None or fdict.get("weight")=="IUISO":
        x += 'Inverse Atomic B<sub>iso</sub>'
    elif fdict.get("weight")=="NONE":
        x += 'No Weighting'
    x += '</td>'

    ## Include Atoms
    x += '<td>'
    if fdict.get("include_atoms")==None or fdict.get("include_atoms")=="ALL":
        x += 'Include All Atoms'
    elif fdict.get("include_atoms")=="MAINCHAIN":
        x += 'Main Chain Atoms'
    elif fdict.get("include_atoms")=="CA":
        x += 'C-Alpha Atoms'
        
    x += '</td>'

    x += '</tr>'

    ## end form
    if fdict.has_key("removebutton"):
        x += '<form '\
             'enctype="multipart/form-data" '\
             'action="webtlsmd.cgi" '\
             'method="get">'

        ## Job ID, user, passwd
        x += '<input type="hidden" name="page" value="%s">' % (
            fdict.get("page", "index"))
        x += '<input type="hidden" name="edit_form" value="TRUE">'
        x += '<input type="hidden" name="job_id" value="%s">' % (
            fdict["job_id"])
        x += '<input type="hidden" name="user" value="%s">' % (
            fdict["user"])
        x += '<input type="hidden" name="passwd" value="%s">' % (
            fdict["passwd"])

        x += '<tr>'
        x += '<td colspan="3" align="left">'
        x += '<input type="submit" name="submit" value="Remove Job">'
        x += '</td>'
        x += '</form>'

    x += '</tr>'
    x += '</table>'
    return x


def check_job_id(form, webtlsmdd):
    """Retrieves and confirms the job_id from a incomming form.  Returns
    None on error, or the job_id on success.
    """
    if form.has_key("job_id"):
        job_id = form["job_id"].value
        if len(job_id)<20:
            if job_id.startswith("TLSMD"):
                if webtlsmdd.job_exists(job_id):
                    return job_id
    return None


def vet_data(data, max_len):
    if len(data)>max_len:
        return False
    if not data.isalnum():
        return False
    return True

def vet_email(email):
    if len(email)>45:
        return False
    return True

def vet_comment(comment):
    if len(comment)>50:
        return False
    return True


def extract_job_edit_form(form, webtlsmdd):
    """Extract the input from the Job Edit Form and update the webtlsmdd
    database with the information.
    """
    if not form.has_key("edit_form"):
        return False

    job_id = check_job_id(form, webtlsmdd)
    if job_id==None:
        return False

    if form.has_key("user"):
        user = form["user"].value.strip()
        if vet_data(user, 10):
            webtlsmdd.job_data_set(job_id, "user", user)

    if form.has_key("passwd"):
        passwd = form["passwd"].value
        if len(passwd)<10:
            webtlsmdd.job_data_set(job_id, "passwd", passwd)
        
    if form.has_key("email"):
        email = form["email"].value.strip()
        if vet_email(email):
            webtlsmdd.job_data_set(job_id, "email", email)

    if form.has_key("structure_id"):
        structure_id = form["structure_id"].value.strip()
        if vet_data(structure_id, 4):
            webtlsmdd.job_data_set(job_id, "structure_id", structure_id)

    if form.has_key("comment"):
        comment = form["comment"].value.strip()
        if vet_comment(comment):
            webtlsmdd.job_data_set(job_id, "comment", comment)

    chains = webtlsmdd.job_data_get(job_id, "chains")
    for cdict in chains:
        if form.has_key(cdict["name"]):
            cdict["selected"] = True
        else:
            cdict["selected"] = False
    webtlsmdd.job_data_set(job_id, "chains", chains)

    if form.has_key("tls_model"):
        tls_model = form["tls_model"].value.strip()
        if tls_model in ["ISOT", "ANISO"]:
            webtlsmdd.job_data_set(job_id, "tls_model", tls_model)

    if form.has_key("weight"):
        weight = form["weight"].value.strip()
        if weight in ["NONE", "IUISO"]:
            webtlsmdd.job_data_set(job_id, "weight", weight)

    if form.has_key("include_atoms"):
        include_atoms = form["include_atoms"].value.strip()
        if include_atoms in ["ALL", "MAINCHAIN", "CA"]:
            webtlsmdd.job_data_set(job_id, "include_atoms", include_atoms)

    return True

def remove_job(webtlsmdd, job_id):
    """Removes job from database and deletes working directory and
    contents.
    """
    job_dir = webtlsmdd.job_data_get(job_id, "job_dir")

    if job_dir and \
       job_dir.startswith(TLSMD_WORK_DIR) and \
       os.path.isdir(job_dir):

        for root, dirs, files in os.walk(job_dir, topdown=False):
            for name in files:
                os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))

        os.rmdir(job_dir)

    webtlsmdd.job_delete(job_id)
         

class Page(object):
    def __init__(self, form):
        self.form = form

    def html_head_nocgi(self, title):
        x  = ''
        x += '<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" '
        x += '"http://www.w3.org/TR/html4/loose.dtd">\n\n'
        x += '<html>'
        x += '<head>'
        x += '  <title>%s</title>' % (title)
        x += '  <style type="text/css" media=screen>'
        x += '  <!-- '
        x += '  BODY { background-color: white;'
        x += '         margin-left: 5%; margin-right: 5%;'
        x += '         border-left: 5%; border-right: 5%;'
        x += '         margin-top: 2%; border-top: 2%;}'
        x += '  -->'
        x += '  </style>'
        x += '</head>'
        x += '<body>'
        return x

    def html_head(self, title):
        x = ''
        x += 'Content-Type: text/html\n\n'
        x += self.html_head_nocgi(title)
        return x

    def html_foot(self):
        x = ''
        x += '<center>'
        x += '<p><small><b>Version %s</b> Last Modified %s' % (
            VERSION, LAST_MODIFIED_DATE)
        x += ' by %s ' % (LAST_MODIFIED_BY)
        x += '<i>%s</i></small></p>' % (LAST_MODIFIED_BY_EMAIL)
        x += '</center>'
        x += '</body></html>'
        return x


class ErrorPage(Page):
    def __init__(self, form, text=None):
        Page.__init__(self, form)
        self.text = text
    
    def html_page(self):
        title = 'TLSMD: Error'
        
        x  = ''
        x += self.html_head(title)
        x += html_title(title)
        x += html_nav_bar()
        x += '<br>'
        x += '<center><h3>A Error Occured</h3></center>'
        if self.text!=None:
            x += '<pre>%s</pre>' % (self.text)
        x += self.html_foot()
        return x


class IndexPage(Page):
    def html_page(self):
        title = 'TLSMD: TLS Minimized Domains'

        x  = ''
        x += self.html_head(title)
        x += html_title(title)
        x += html_nav_bar("index")

        x += '<form enctype="multipart/form-data" '\
             'action="webtlsmd.cgi" '\
             'method="post">'
        
        x += '<input type="hidden" name="page" value="submission_form">'
        x += '<center>'
        x += '<table>'
        x += '<tr>'
        x += '<td align="left">Upload PDB File:</td>'
        x += '<td><input name="pdbfile" size="50" type="file"></td>'
        x += '</tr>'
        x += '<tr>'
        x += '<td colspan="2" align="center">'
        x += '<input value="Submit" type="submit">'
        x += '</td>'
        x += '</tr>'
        x += '</table>'
        x += '</center>'
        x += '</form>'

        x += get_documentation_block("INDEX")

        x += self.html_foot()
        return x


class QueuePage(Page):
    def html_head_nocgi(self, title):
        x  = ''
        x += '<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" '
        x += '"http://www.w3.org/TR/html4/loose.dtd">\n\n'
        x += '<html>'
        x += '<head>'
        x += '  <meta http-equiv="refresh" content="60">'
        x += '  <title>%s</title>' % (title)
        x += '  <style type="text/css" media=screen>'
        x += '  <!-- '
        x += '  BODY { background-color: white;'
        x += '         margin-left: 5%; margin-right: 5%;'
        x += '         border-left: 5%; border-right: 5%;'
        x += '         margin-top: 2%; border-top: 2%;}'
        x += '  -->'
        x += '  </style>'
        x += '</head>'
        x += '<body>'
        return x
    
    def html_foot(self):
        x = ''
        x += '<center>'
        x += '<p><small><b>Version %s</b> Last Updated %s</p>' % (
            VERSION,
            timestring(time.time()))
        x += '</center>'
        x += '</body></html>'
        return x
    
    def get_job_list(self):
        """Get a list of all the jobs in the job queue file.
        """
        jl = webtlsmdd.job_list()

        job_list = []
        for jdict in jl:
            if jdict["state"]!="completed":
                job_list.append(jdict)

        return job_list

    def html_running_job_table(self, job_list):
        ## get the job dictionary of the running job
        run_jdict = None
        for jdict in job_list:
            if jdict.get("state")=="running":
                run_jdict = jdict
                break

        ## if there is no running job
        if run_jdict==None:
            x  = '<center>'
            x += '<h3>No Running Job</h3>'
            x += '</center>'
            return x

        jdict = run_jdict

        x  = '<center>'
	x += '<h3>Running Jobs</h3>'
        x += '<table border="1" width="100%">'
        x += '<tr>'
        x += '<th><font size="-5">Job ID</font></th>'
        x += '<th><font size="-5">Struct ID</font></th>'
        x += '<th><font size="-5">User</font></th>'
        x += '<th><font size="-5">Submission Date</font></th>'
        x += '<th><font size="-5">Currently Processing</font></th>'
        x += '<th><font size="-5">Processing Time<br>Used (Hours)</font></th>'
        x += '</tr>'

        x += '<tr>'
        x += '<td><a href="webtlsmd.cgi?page=edit&job_id=%s">%s</a>' % (
            jdict["job_id"] ,jdict["job_id"])
        x += '<td>%s</td>' % (jdict["structure_id"])
        x += '<td>%s</td>' % (jdict["user"])
        x += '<td>%s</td>' % (timestring(jdict["submit_time"]))

        tls_seg = 'Chain <b>%s</b> Residues <b>%s-%s</b>' % (
            jdict.get("run_chain_id", ""),
            jdict.get("run_frag_id1", ""),
            jdict.get("run_frag_id2", ""))
        x += '<td>%s</td>' % (tls_seg)

        begin = jdict.get("run_start_time")
        if begin==None:
            begin = jdict["submit_time"]
        hours = timediffstring(begin, time.time())
        x += '<td>%s</td>' % (hours)

        x += '</tr>'
        x += '</table>'
        x += '</center>'
        return x

    def html_queued_job_table(self, job_list):
        queued_list = []
        for jdict in job_list:
            if jdict.get("state")=="queued":
                queued_list.append(jdict)

        if len(queued_list)==0:
            x  = '<center>'
            x += '<h3>No Queued Jobs</h3>'
            x += '</center>'
            return x

        x  = ''
        x += '<center>'
	x += '<h3>Queued Jobs</h3>'
        x += '<table border="1" width="100%">'
        x += '<tr>'
        x += '<th><font size="-5">Job ID</font></th>'
        x += '<th><font size="-5">Struct ID</font></th>'
        x += '<th><font size="-5">User</font></th>'
        x += '<th><font size="-5">Submission Date</font></th>'
        x += '</tr>'

        for jdict in queued_list:
            x += '<tr>'

            x += '<td><a href="webtlsmd.cgi?page=edit&job_id=%s">%s</a>' % (
                jdict["job_id"] ,jdict["job_id"])
            x += '<td>%s</td>' % (jdict["structure_id"])
            x += '<td>%s</td>' % (jdict["user"])
            x += '<td>%s</td>' % (timestring(jdict["submit_time"]))
                                  
            x += '</tr>'

        x += '</table>'
        x += '</center>'
        return x

    def html_limbo_job_table(self, job_list):
        limbo_list = []
        for jdict in job_list:
            if jdict.get("state") not in ["queued", "running", "completed"]:
                limbo_list.append(jdict)

        if len(limbo_list)==0:
            return None

        x  = ''
        x += '<center>'
	x += '<h3>Partially Completed Jobs</h3>'
        x += '<table border="1" width="100%">'
        x += '<tr>'
        x += '<th><font size="-5">Job ID</font></th>'
        x += '<th><font size="-5">Struct ID</font></th>'
        x += '<th><font size="-5">State</font></th>'
        x += '<th><font size="-5">User</font></th>'
        x += '<th><font size="-5">Submission Date</font></th>'
        x += '</tr>'

        for jdict in limbo_list:
            x += '<tr>'

            x += '<td><a href="webtlsmd.cgi?page=edit&job_id=%s">%s</a>' % (
                jdict["job_id"] ,jdict["job_id"])
            x += '<td>%s</td>' % (jdict["structure_id"])
            x += '<td>%s</td>' % (jdict["state"])
            x += '<td>%s</td>' % (jdict["user"])
            x += '<td>%s</td>' % (timestring(jdict["submit_time"]))
                                  
            x += '</tr>'

        x += '</table>'
        x += '</center>'
        return x

    def html_page(self):
        title = 'TLSMD: Job Queue'
        
        x  = ''
        x += self.html_head(title)
        x += html_title(title)
        x += html_nav_bar("queue")

        job_list = self.get_job_list()
        x += self.html_running_job_table(job_list)
	x += '<br>'	
        x += self.html_queued_job_table(job_list)

        limbo = self.html_limbo_job_table(job_list)
        if limbo!=None:
            x += '<br>'
            x += limbo

        x += self.html_foot()
        return x


class CompletedPage(Page):
    def html_foot(self):
        x = ''
        x += '<center>'
        x += '<p><small><b>Version %s</b> Last Updated %s</p>' % (
            VERSION,
            timestring(time.time()))
        x += '</center>'
        x += '</body></html>'
        return x
    
    def get_job_list(self):
        """Get a list of all the jobs in the job queue file.
        """
        jl = webtlsmdd.job_list()

        job_list = []
        for jdict in jl:
            if jdict["state"] in ["completed", "defunct"]:
                job_list.append(jdict)
        job_list.reverse()

        return job_list

    def html_job_table(self, job_list):
        x  = ''
        x += '<center>'
        x += '<table border="1" width="100%">'
        x += '<tr>'
        x += '<th><font size="-5">Job ID</font></th>'
        x += '<th><font size="-5">Struct ID</font></th>'
        x += '<th><font size="-5">User</font></th>'
        x += '<th><font size="-5">Status</font></th>'
        x += '<th><font size="-5">Submission Date</font></th>'
	x += '<th><font size="-5">Processing Time<br> Used (Hours)</font></th>'
        x += '</tr>'

        for jdict in job_list:
            x += '<tr>'
            x += '<td><a href="webtlsmd.cgi?page=edit&job_id=%s">%s</a>' % (
                jdict["job_id"] ,jdict["job_id"])
            x += '<td>%s</td>' % (jdict["structure_id"])
            x += '<td>%s</td>' % (jdict["user"])
            x += '<td>%s</td>' % (jdict["state"])
            x += '<td>%s</td>' % (timestring(jdict["submit_time"]))

            begin = jdict.get("run_start_time")
            if begin==None:
	        begin = jdict["submit_time"]
            end = jdict.get("run_end_time")
	    if end==None:
                hours = "Unknown"
	    else:
	        hours = timediffstring(begin, time.time())
            x += '<td>%s</td>' % (hours)

            x += '</tr>'

        x += '</table>'
        x += '</center>'
        return x

    def html_page(self):
        title = 'TLSMD: Completed Jobs'
        
        x  = ''
        x += self.html_head(title)
        x += html_title(title)
        x += html_nav_bar("completed")

        job_list = self.get_job_list()
        if len(job_list)>0:
            x += self.html_job_table(job_list)
        else:
            x += '<center><h3>No Completed Jobs</h3></center>'

        x += self.html_foot()
        return x


class EditPage(Page):
    def html_page(self):
        title = 'TLSMD: Edit Job'
        
        x  = ''
        x += self.html_head(title)
        x += html_title(title)

        if self.check_auth():
            x += html_nav_bar()

            if self.form.has_key("submit") and \
               self.form["submit"].value=="Remove Job":
                x += self.remove()
            else:
                x += self.edit()

        else:
            x += html_nav_bar("edit")
            x += self.auth()
        
        x += self.html_foot()
        return x

    def check_auth(self):
        job_id = check_job_id(self.form, webtlsmdd)
        if job_id==None:
            return False
        
        if self.form.has_key("submit"):
            if self.form["submit"].value=="Cancel":
                return False
    
        user = ""
        if self.form.has_key("user"):
            user = self.form["user"].value.strip()

        passwd = ""
        if self.form.has_key("passwd"):
            passwd = self.form["passwd"].value

        user_cmp = webtlsmdd.job_data_get(job_id, "user")
        if len(user_cmp)>0 and user_cmp!=user:
            return False

        passwd_cmp = webtlsmdd.job_data_get(job_id, "passwd")
        if len(passwd_cmp)>0 and passwd_cmp!=passwd:
            return False
        
        return True

    def auth(self):
        ## try to retrieve the job_id for authorization
        job_id = check_job_id(self.form, webtlsmdd)

        ## html generation
        x  = ''
        x += '<center>'

        x += '<form enctype="multipart/form-data" '\
             'action="webtlsmd.cgi" '\
             'method="post">'
        x += '<input type="hidden" name="page" value="edit">'

        x += '<table width="40%" border="1">'

        x += '<tr>'
        x += '<td align="center">'
        x += '<table>'

        if job_id==None:
            x += '<tr>'
            x += '<td align="right">Job ID:</td>'
            x += '<td>'

            x += '<select name="job_id">'

            for jdict in webtlsmdd.job_list():
                x += '<option value="%s">%s</option>' % (
                    jdict["job_id"], jdict["job_id"])
            
            x += '</select>'

            x += '</td>'
            x += '</tr>'
        else:
            x += '<tr>'
            x += '<td align="right">Job ID:</td>'
            x += '<td><b>%s</b></td>' % (job_id)
            x += '<input type="hidden" name="job_id" value="%s">' % (job_id)
            x += '</tr>'

        x += '<tr>'
        x += '<td align="right">Username:</td>'
        x += '<td><input type="text" name="user" '\
             'size="10" maxlength="10"></td>'
        x += '</tr>'

        x += '<tr>'
        x += '<td align="right">Password:</td>'
        x += '<td><input type="text" name="passwd" '\
             'size="10" maxlength="10"></td>'
        x += '</tr>'
        
        x += '</table>'

        x += '<tr>'
        x += '<td align="right">'
        x += '<input type="submit" name="submit" value="OK">'
        x += LINK_SPACE
        x += '<input type="submit" name="submit" value="Cancel">'
        x += '</td>'
        x += '</tr>'
        
        x += '</table>'
        x += '</center>'
        x += '</form>'
        return x

    def edit(self):
        x = ''
        job_id = check_job_id(self.form, webtlsmdd)
        if job_id==None:
            x += '<center>'
            x += '<h3><b>ERROR:</b> no job_id in form</h3>'
            x += '</center>'
            return x

        ## if the job is not in the "queued" state, then it is not
        ## safe to edit
        state = webtlsmdd.job_data_get(job_id, "state")
        if state=="queued":
            extract_job_edit_form(self.form, webtlsmdd)

        ## get the state dictionary for the entire job
        fdict = webtlsmdd.job_get_dict(job_id)
        fdict["page"] = "edit"
        fdict["removebutton"] = True
            
        if state=="running" or state=="completed":
            x += html_job_nav_bar(webtlsmdd, job_id)
            x += html_job_info_table(fdict)
        else:
            x += html_job_edit_form(fdict)
        
        return x

    def remove(self):
        job_id = check_job_id(self.form, webtlsmdd)
        if job_id==None:
            return '<center><h3>'\
                   '<b>ERROR:</b> no job_id in form</h3></center>'
        
        remove_job(webtlsmdd, job_id)

        x  = ''
        x += '<center>'
        x += '<h3>Job %s has been removed.</h3>' % (job_id)
        x += '</center>'
        return x


class SubmissionException(Exception):
    def __init__(self, html):
        Exception.__init__(self)
        self.html = html


class SubmissionFormPage(Page):
    def html_page(self):        
        title = 'TLSMD: Fill Out Job Submission Form'
        
        x  = ''
        x += self.html_head(title)
        x += html_title(title)

        try:
            job_id = self.prepare_submission()
        except SubmissionException, err:
             x += html_nav_bar()
             x += '<center><h3>ERROR:%s</h3></center>' % (err.html)
        else:
            x += '<center>'
            x += '<h3>You must read and fill out this form to complete'
            x += ' your submission!</h3>'
            x += '</center>'
            
            x += get_documentation_block("SUBMIT1")
            x += self.job_edit_form(job_id)

        x += self.html_foot()
        return x

    def job_edit_form(self, job_id, show_warnings=False):
        fdict = webtlsmdd.job_get_dict(job_id)
        fdict["page"] = "submission"
        x = html_job_edit_form(fdict)
        return x

    def prepare_submission(self):
        if self.form.has_key("pdbfile")==False or \
           self.form["pdbfile"].file==None:
            raise SubmissionException("No PDB file uploaded")
	
        ## make working directory
        try:
            os.chdir(TLSMD_WORK_DIR)
        except os.error, err:
            raise SubmissionException(
                '<p>Cannot change to working directory: %s</p>' % (str(err)))

        job_id = webtlsmdd.job_new()
        os.umask(022)
        try:
            os.mkdir(job_id)
        except os.error, err:
            webtlsmdd.job_delete(job_id)
            raise SubmissionException(
                '<p>Cannot make directory: %s</p>' % (str(err)))

        job_dir = os.path.join(TLSMD_WORK_DIR, job_id)
        os.chdir(job_dir)
        webtlsmdd.job_data_set(job_id, "job_dir", job_dir)

        ## save PDB file
        pdb_filename = "struct.pdb"

        num_lines = 0
        infil = self.form["pdbfile"].file
        outfil = open(pdb_filename,"w")
        while True:
            ln = infil.readline()
            if not ln:
                break
            outfil.write(ln)
            num_lines += 1
        outfil.close()    

        ## error out if there weren't many lines
        if num_lines<10:
            webtlsmdd.job_delete(job_id)
            raise SubmissionException(
                '<p>Only Recieved %d lines</p>' % (num_lines))

        webtlsmdd.job_data_set(job_id, "pdb_filename", pdb_filename)

        job_url = "%s/%s" % (TLSMD_WORK_URL, job_id)
        webtlsmdd.job_data_set(job_id, "job_url", job_url)

        log_url = "%s/log.txt" % (job_url)
        webtlsmdd.job_data_set(job_id, "log_url", log_url)

        analysis_dir      = "%s/ANALYSIS" % (job_dir)
        analysis_base_url = "%s/ANALYSIS" % (job_url)
        analysis_url      = "%s/ANALYSIS/index.html" % (job_url)
        webtlsmdd.job_data_set(job_id, "analysis_dir", analysis_dir)
        webtlsmdd.job_data_set(
            job_id, "analysis_base_url", analysis_base_url)
        webtlsmdd.job_data_set(job_id, "analysis_url", analysis_url) 

        ip_addr = os.environ.get("REMOTE_ADDR", "Unknown")
        webtlsmdd.job_data_set(job_id, "ip_addr", ip_addr)

        ## submission time and initial state
        webtlsmdd.job_data_set(job_id, "state", "submit1")
        webtlsmdd.job_data_set(job_id, "submit_time", time.time())

        ## now load the structure and build the submission form
        struct = LoadStructure(fil=pdb_filename)

        webtlsmdd.job_data_set(
            job_id, "structure_id", struct.structure_id)

        ## Select Chains for Analysis
        chains = []
        for chain in struct.iter_chains():
            if chain.count_amino_acids()<10:
                continue

            ## form name
            cb_name = 'CHAIN%s' % (chain.chain_id)
            
            ## create chain description label cb_desc
            cb_desc = 'Chain %s, %d Residues' % (
                chain.chain_id,
                chain.count_amino_acids())

            listx = []
            i = 0
            for frag in chain.iter_fragments():
                i += 1
                if i>5: break
                listx.append(frag.res_name)
            cb_preview = string.join(listx, " ")
            
            cdict = {}
            chains.append(cdict)
            cdict["chain_id"] = chain.chain_id
            cdict["name"]     = cb_name
            cdict["desc"]     = cb_desc
            cdict["preview"]  = cb_preview
            cdict["selected"] = True

        webtlsmdd.job_data_set(job_id, "chains", chains)

        ## defaults
        webtlsmdd.job_data_set(job_id, "user", "")
        webtlsmdd.job_data_set(job_id, "passwd", "")
        webtlsmdd.job_data_set(job_id, "email", "")
        webtlsmdd.job_data_set(job_id, "comment", "")

        webtlsmdd.job_data_set(job_id, "tls_model", "ISOT")
        webtlsmdd.job_data_set(job_id, "weight", "IUISO")
        webtlsmdd.job_data_set(job_id, "include_atoms", "ALL")

	return job_id


class SubmissionPage(Page):
    def html_page(self):
        try:
            job_id = self.complete_submission()
	except SubmissionException, err:
	    title = 'TLSMD: Job Submission Failed'
            html  = '<center><h3>ERROR: %s</h3></center>' % (err.html)
	else:
            title = 'TLSMD: Job Submission Succeeded'
            html  = '<center>'
	    html += '<h3>SUCCESS! Your job ID is %s</h3>' % (job_id)
	    html += '</center>'
	    
        x  = self.html_head(title)
        x += html_title(title)
        x += html_nav_bar()
        x += html
        x += self.html_foot()
        return x

    def complete_submission(self):
        ## check for submission key
        if not self.form.has_key("submit"):
            raise SubmissionException('Submission Error')

        ## get job_id; verify job exists
        job_id = check_job_id(self.form, webtlsmdd)
        if job_id==None:
            raise SubmissionException('Submission Error')

        ## make sure the job is in the right state to be submitted
	state = webtlsmdd.job_data_get(job_id, "state")
	if state=="queued":
	    raise SubmissionException("Your job is already queued")
    	elif state=="running":
	    raise SubmissionException("Your job is already running")

        ## verify the submission IP address
        ip_addr = os.environ.get("REMOTE_ADDR", "Unknown")
        ip_addr_verify = webtlsmdd.job_data_get(job_id, "ip_addr")
        if ip_addr!=ip_addr_verify:
            raise SubmissionException('Submission IP Address Mismatch')

        ## completely remove the job
        if self.form["submit"].value=="Cancel":
            remove_job(webtlsmdd, job_id)
            raise SubmissionException('You cancelled the job')

        extract_job_edit_form(self.form, webtlsmdd)

        ## if everything with the form is okay, then change
        ## the job state to queued
        webtlsmdd.job_data_set(job_id, "state", "queued")

        return job_id


def main():
    page = None
    form = cgi.FieldStorage()

    ## mode of app
    if not form.has_key("page"):
        page = IndexPage(form)

    elif form["page"].value=="edit":
        page = EditPage(form)

    elif form["page"].value=="queue":
        page = QueuePage(form)

    elif form["page"].value=="submission_form":
        page = SubmissionFormPage(form)

    elif form["page"].value=="submission":
        page = SubmissionPage(form)

    elif form["page"].value=="completed":
        page = CompletedPage(form)

    if page==None:
        page = IndexPage(form)

    try:
        print page.html_page()

    except xmlrpclib.Fault, err:
        page = ErrorPage(form, "xmlrpclib.Fault: " +str(err))
        print page.html_page()

    except socket.error, err:
        page = ErrorPage(form, "socket.error: " + str(err))
        print page.html_page()


if __name__=="__main__":
    main()
    sys.exit(0)
