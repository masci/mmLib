## TLS Minimized Domains (TLSMD)
## Copyright 200-20052 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys
import time
import socket
import string
import random

import xmlrpclib
import cgitb; cgitb.enable()
import cgi

from mmLib import FileLoader

import const
import conf

## GLOBALS
webtlsmdd = xmlrpclib.ServerProxy(conf.WEBTLSMDD)

STYLE_SHEET = """
BODY { background-color:white;margin-left:5%;margin-right:5%;border-left:5%;border-right:5%;margin-top:2%;border-top:2%; }
.step_title { line-height:2.0em;font-size:large }
.submit_table { padding:10px;border-width:thin;border-color:blue;border-style:solid;background-color:#eeeeee; }
.inner_table { font-size:small;width:90%;padding:10px;border-width:thin;border-color:black;border-style:solid;background-color:#dddddd; }
.inner_title { background-color:#ccccff; line-height:2.0em; }
.ninner_table { font-size:small;width:95%;padding:2%;border-width:thin;
                border-color:#aaaaaa;border-style:solid;background-color:#eeeeee }
"""


def timestring(secs):
    tm_struct = time.localtime(secs)
    return time.strftime("%m-%d-%y %H:%M %Z" ,tm_struct)

def secdiffstring(secs):
    secs = int(secs)
    
    hours = secs / 3600
    secs = secs - (hours * 3600)
 
    min = secs / 60
    secs = secs - (min * 60)

    x = "%1d:%2d.%2d" % (hours, min, secs)
    return x.replace(" ", "0")

def timediffstring(begin, end):
    secs  = int(end - begin)
    return secdiffstring(secs)

def html_title(title):
    """Title
    """
    return '<center><h1>%s</h1></center>' % (title)

def html_nav_bar(page_name=None):
    """Site navigation bar.
    """
    l = ['<center>']

    if page_name=="home":
        l.append('Home')
    else:
        l.append('<a href="%s/index.html">Home</a>' % (conf.TLSMD_BASE_URL))

    l.append(LINK_SPACE)

    if page_name=="":
        l.append('Start a New Job')
    else:
        l.append('<a href="webtlsmd.cgi?page=submit1">Start a New Job</a>')

    l.append(LINK_SPACE) 

    if page_name=="queue":
        l.append('Job Status')
    else:
        l.append('<a href="webtlsmd.cgi">Job Status</a>')

    l.append(LINK_SPACE)

    l.append('<a href="%s/examples/index.html">Examples</a>' % (conf.TLSMD_BASE_URL))

    l.append(LINK_SPACE)
    
    l.append('<a href="%s/documentation.html">Documentation</a>' % (conf.TLSMD_BASE_URL))
    l.append('</center>')
    l.append('<br>')

    return "".join(l)

def html_job_nav_bar(webtlsmdd, job_id):
    """Navigation bar to the TLSMD output files.
    """
    analysis_dir = webtlsmdd.job_data_get(job_id, "analysis_dir")
    analysis_index = os.path.join(analysis_dir, "index.html")
    analysis_url = webtlsmdd.job_data_get(job_id, "analysis_url")

    job_dir = webtlsmdd.job_data_get(job_id, "job_dir")
    logfile = os.path.join(job_dir, "log.txt")
    log_url = webtlsmdd.job_data_get(job_id, "log_url")

    if not os.path.isfile(analysis_index) and not os.path.isfile(logfile):
        return ''

    x  = ''
    x += '<center>'
    x += '<h3>'

    if os.path.isfile(analysis_index):
	x += '<a href="%s">Click Here: View Completed TLSMD Analysis</a>' % (analysis_url)

    x += LINK_SPACE

    if os.path.isfile(logfile):
        x += '<a href="%s">View TLSMD Logfile</a>' % (log_url)

    x += '</h3>'
    x += '</center>'
    x += '<br>'
    return x


def html_job_edit_form(fdict):
    x  = ''
    x += '<center>'

    x += '<form enctype="multipart/form-data" action="webtlsmd.cgi" method="post">'
    x += '<input type="hidden" name="page" value="%s">' % (fdict.get("page", "index"))
    x += '<input type="hidden" name="edit_form" value="TRUE">'
    x += '<input type="hidden" name="job_id" value="%s">' % (fdict["job_id"])

    x += '<table border="1" width="100%">'

    ## user/email/passcode/structure name
    x += '<tr>'
    x += '<th colspan="2">User Information</th>'
    x += '<th>Session Information</th>'
    x += '</tr>'

    x += '<tr><td colspan="2">'
    x += '<table>'

    ## keep job private
    x += '<tr><td></td>'
    x += '<td>'
    x += '<label>'
    x += '<input type="checkbox" name="private_job" value="TRUE">'
    x += 'Keep Job Private'
    x += '</label>'
    x += '</td>'
    x += '</tr>'

    ## email address
    x += '<tr>'
    x += '<td align="right"><label>EMail Address:</td><td>'
    x += '<input type="text" name="email" value="%s" size="25" maxlength="40">' % (fdict.get("email", ""))
    x += '</label></td>'
    x += '</tr>'

    ## structure code
    x += '<tr>'
    x += '<td align="right"><label>Structure Code:</td><td>'
    x += '<input type="text" name="structure_id" value="%s" size="4" maxlength="4">' % (fdict.get("structure_id", ""))
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
    for cdict in fdict.get("chains", []):
        x += '<tr><td>'
        x += '<label>'
        if cdict["selected"]==True:
            x += '<input type="checkbox" name="%s" value="TRUE" checked>' % (cdict["name"])
        else:
            x += '<input type="checkbox" name="%s" value="TRUE">' % (cdict["name"])
        x += '%s' % (cdict["desc"])
        x += '</label>'

        x += '</td></tr>'

    x += '</table></td></tr>'

    ## end form

    x += '<tr><td colspan="3">'

    x += '<table width="100%">'
    x += '<tr>'
    
    x += '<td align="left">'
    if fdict.has_key("removebutton"):
        x += '<input type="submit" name="submit" value="Remove Job">'
    x += '</td>'
    
    x += '<td align="right">'
    x += '<input type="submit" name="submit" value="Next">'
    x += '</tr>'
    x += '</table>'

    x += '</td></tr>'
    x += '</table>'
    x += '</form>'
    return x

def html_session_info_table(fdict):
    if fdict.has_key("submit_time"):
        date = timestring(fdict["submit_time"])
    else:
        date = ""

    l = ['<table class="inner_table">',

         '<tr class="inner_title"><th>',
         '<a id="cid2" href="javascript:',
         "ToggleDivVisibility('cid2','id2','Show Session Information','Hide Session Information')",
         '">Show Session Information</a>',
         '</th></tr>',

         '<tr><td align="center">',

         '<div id="id2" style="display:none"><table class="ninner_table">',
         '<tr><td align="right">TLSMD Job ID:</td>',
         '<td><b>%s</b></td></tr>' % (fdict["job_id"]),
         
         '<tr><td align="right">Job State:</td>',
         '<td><b>%s</b></td></tr>' % (fdict["state"]),
         
         '<tr><td align="right">Submission IP Address: </td>',
         '<td><b>%s</b></td></tr>' % (fdict.get("ip_addr", "")),
         
         '<tr><td align="right">Submission Date: </td>',
         '<td><b>%s</b></td></tr>' % (date),
         '</table></div>',
         
         '</table>']

    return "".join(l)

def html_user_info_table(fdict):
    l = ['<table class="inner_table">',

         '<tr class="inner_title"><th colspan="2">User Information</th></tr>',

         '<tr><td align="center">',
         '<table class="ninner_table">',

         '<tr>',
         '<td align="right"><label for="user_name">Your Name</label></td>',
         '<td><input type="text" id="user_name" name="user_name" value="%s" size="25" maxlength="40"></td>' % (fdict.get("user_name","")),
         '</tr>',

         '<tr>',
         '<td align="right"><label for="email">EMail Address</label></td>',
         '<td><input type="text" id="email" name="email" value="%s" size="25" maxlength="40"></td>' % (fdict.get("email", "")),
         '</tr>',

         '</table>',
         '</td></tr></table>']

    return "".join(l)

def html_program_settings_table(fdict):

    opt_plot_svg = ""
    opt_plot_png = ""

    opt_atoms_all = ""
    opt_atoms_mnchn = ""

    l = ['<table class="inner_table">',
         '<tr class="inner_title"><th>TLSMD Program Options</th></tr>',

         '<tr><td align="center">',
         '<table width="100%">',
         '<tr><td align="center" valign="top">',

         ## left table
         '<table class="ninner_table">',

         '<tr><td>',
         '<label><input type="checkbox" id="private_job" name="private_job" value="TRUE">Keep Job Private</label>',
         '</td></tr>',

         '<tr><td>',
         '<label for="structure_id">4-Letter Structure ID </label>',
         '<input type="text" id="structure_id" name="structure_id" value="%s" size="4" maxlength="4">' % (fdict.get("structure_id", "")),
         '</td></tr>',

         '</table>',

         '</td><td align="center" valign="top">',

         ## right table
         '<table class="ninner_table">',
         '<tr style="line-height:2em"><th>Select Chains for Analysis</th></tr>']
         
    for cdict in fdict.get("chains", []):
        if cdict["selected"]==True:
            x = '<label><input type="checkbox" id="%s" name="%s" value="TRUE" checked>' % (cdict["name"], cdict["name"])
        else:
            x = '<label><input type="checkbox" id="%s" name="%s" value="TRUE">' % (cdict["name"], cdict["name"])
            
        l +=['<tr><td>', x, cdict["desc"], '</label></td></tr>' ]

    l +=['</table>',
         
         '</td></tr>',
         '</table>',

         ## advanced options
         '<tr class="inner_title"><th>',
         '<a id="cid1" href="javascript:',
         "ToggleDivVisibility('cid1','id1','Show Advanced Program Options','Hide Advanced Program Options')",
         '">Show Advanced Program Options</a>',
         '</th></tr>',
         
         '<tr><td align="center">',
         '<div id="id1" style="display:none">',
         '<table class="ninner_table">',
         '<tr>',

         '<td valign="top">',
         '<fieldset><legend>Plot Output Format</legend>',
         '<div style="font-size:xx-small">Select the output format for plots.<br>SVG works with the Adobe plugin and Firefox 1.5.</div><br>',
         '<label><input name="plot_format" type="radio" value="PNG" tabindex="35" checked>PNG Images</label><br>',
         '<label><input name="plot_format" type="radio" value="SVG" tabindex="35">SVG</label>',
         '</fieldset>',
         '</td>',

         '<td valign="top">',
         '<fieldset><legend>Atom Class Selection</legend>',
         '<div style="font-size:xx-small">Analyze all protein atoms, or just the main chain atoms.</div><br>',
         '<label><input name="include_atoms" type="radio" value="ALL" tabindex="35" checked>All Atoms</label><br>',
         '<label><input name="include_atoms" type="radio" value="MAINCHAIN" tabindex="35">Mainchain Atoms (N,CA,C,O,CB)</label>',
         '</fieldset>',
         '</td>',

         '</tr>',
         '</table>',

         '</div>',
         '</td></tr>',
         '</table>']
         
    return "".join(l)


def html_job_edit_form2(fdict, title=""):
    if fdict.has_key("removebutton"):
        remove_button = '<input type="submit" name="submit" value="Remove Job">'
    else:
        remove_button = ''
    
    l = ['<script language=javascript type="text/javascript">',
         'function ToggleDivVisibility(control_id, target_id, show_val, hide_val) {',
         '  var ctrl_element = document.getElementById(control_id);',
         '  var target_element = document.getElementById(target_id);',
         
         '  if (target_element.style.display != "none") {',
         '    target_element.style.display = "none";',
         '    ctrl_element.firstChild.nodeValue = show_val;',
         '  } else {',
         '    target_element.style.display = "inline";',
         '    ctrl_element.firstChild.nodeValue = hide_val;',
         '  }',
         '}',
         '</script>',

         '<center>',

         '<form enctype="multipart/form-data" action="webtlsmd.cgi" method="post">',
         
         '<input type="hidden" name="page" value="%s">' % (fdict.get("page", "index")),
         '<input type="hidden" name="edit_form" value="TRUE">',
         '<input type="hidden" name="job_id" value="%s">' % (fdict["job_id"]),
         
         '<table width="100%" class="submit_table">',
         '<tr><th class="step_title">%s</th></tr>' % (title),

         '<tr><td align="center">', html_user_info_table(fdict), '</td></tr>',
         '<tr><td align="center">', html_program_settings_table(fdict), '</td></tr>',
         '<tr><td align="center">', html_session_info_table(fdict), '</td></tr>',

         '<tr><td align="center"><input type="submit" name="submit" value="Submit Job"></td></tr>',

         '</table>',
         '</form>',
         '</center>']

    return "".join(l)

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
        date = "---"
    x += '<td><b>%s</b></td></tr>' % (date)

    x += '<tr><td align="right">Processing Start Date: </td>'
    if fdict.has_key("run_start_time"):
        date = timestring(fdict["run_start_time"])
    else:
        date = "---"
    x += '<td><b>%s</b></td></tr>' % (date)

    x += '<tr><td align="right">Processing End Date: </td>'
    if fdict.has_key("run_end_time"):
        date = timestring(fdict["run_end_time"])
    else:
        date = "---"
    x += '<td><b>%s</b></td></tr>' % (date)

    x += '<tr><td align="right">Processing Time(HH:MM): </td>'
    if fdict.has_key("run_end_time") and fdict.has_key("run_start_time"):
        hours = timediffstring(fdict["run_start_time"], fdict["run_end_time"])
    else:
        hours = "---"
    x += '<td><b>%s</b></td></tr>' % (hours)

    x += '</table></td>'

    x += '</tr>'

    ## Select Chains for Analysis
    x += '<tr><th colspan="3">Selected Chains</th></tr>'

    x += '<tr><td colspan="3">'
    x += '<table cellpadding="5">'
    x += '<tr><th><font size="-5">Chain</font></th><th><font size="-5">Processing Time (HH:MM.SS)</font></th></tr>'
    for cdict in fdict.get("chains", []):
        x += '<tr><td>'
        if cdict["selected"]==True:
            x += '<tr>'
	    
            x += '<td>%s</td>' % (cdict["desc"])

            if cdict.has_key("processing_time"):
                hours = secdiffstring(cdict["processing_time"])
	    else:
		hours = "---"
            x += '<td>%s</td>' % (hours)

	    x += '</tr>' 

    x += '</table></td></tr>'

    ## column titles
    x += '<tr>'
    x += '<th>TLS Model</th><th>Least Squares Weighting</th><th>Include Atoms</th>'
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
        x += '<form enctype="multipart/form-data" action="webtlsmd.cgi" method="post">'

        ## Job ID, user, passwd
        x += '<input type="hidden" name="page" value="%s">' % (fdict.get("page", "index"))
        x += '<input type="hidden" name="edit_form" value="TRUE">'
        x += '<input type="hidden" name="job_id" value="%s">' % (fdict["job_id"])
        x += '<input type="hidden" name="user" value="%s">' % (fdict["user"])
        x += '<input type="hidden" name="passwd" value="%s">' % (fdict["passwd"])

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

    if form.has_key("private_job"):
        webtlsmdd.job_data_set(job_id, "private_job", True)

    if form.has_key("user_name"):
        user_name = form["user_name"].value.strip()
        user_name = user_name[:100]
        webtlsmdd.job_data_set(job_id, "user_name", user_name)
            
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
        if include_atoms in ["ALL", "MAINCHAIN"]:
            webtlsmdd.job_data_set(job_id, "include_atoms", include_atoms)

    if form.has_key("plot_format"):
        plot_format = form["plot_format"].value.strip()
        if plot_format in ["PNG", "SVG"]:
            webtlsmdd.job_data_set(job_id, "plot_format", plot_format)


    return True


def remove_job(webtlsmdd, job_id):
    """Removes job from database and deletes working directory and
    contents.
    """
    job_dir = webtlsmdd.job_data_get(job_id, "job_dir")

    if job_dir and job_dir.startswith(conf.TLSMD_WORK_DIR) and os.path.isdir(job_dir):

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
        x += '<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">'
        x += '<html>'
        x += '<head>'
        x += '  <title>%s</title>' % (title)
        x += '  <style type="text/css" media=screen>'
        x += '  <!-- '
        x += STYLE_SHEET
        x += '  -->'
        x += '  </style>'
        x += '</head>'
        x += '<body>'
        return x

    def html_head(self, title):
        return 'Content-Type: text/html\n\n' + self.html_head_nocgi(title)

    def html_foot(self):
        l = ['<center>',
             '<p><small><b>Version %s</b> Last Modified %s' % (const.VERSION, const.RELEASE_DATE),
             '</small></p>',
             '</center>',
             '</body></html>']
        
        return "".join(l)


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


class QueuePage(Page):
    def __init__(self, form):
        Page.__init__(self, form)
        if self.form.has_key("admin"):
            self.admin = self.verify_admin(self.form["admin"].value)
        else:
            self.admin = False

    def verify_admin(self, passcode):
        try:
            code = open(conf.ADMIN_PASSWORD_FILE, "r").read().strip()
	except IOError:
            return False
        return code==passcode

    def rcsb_href(self, struct_id):
        if struct_id.lower()=="xxxx":
            return struct_id
        return '<a href="http://pdbbeta.rcsb.org/pdb/explore.do?structureId=%s">%s</a>' % (struct_id, struct_id)

    def html_head_nocgi(self, title):
        l = ['<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">',
             '<html>',
             '<head>',
             '  <meta http-equiv="refresh" content="120">',
             '  <title>%s</title>' % (title),
             '  <style type="text/css" media=screen>',
             '  <!-- ',
             STYLE_SHEET,
             '  -->',
             '  </style>',
             '</head>',
             '<body>']
        
        return "".join(l)
    
    def html_foot(self):
        l = ['<center>',
             '<p><small><b>Version %s</b> Last Updated %s</p>' % (conf.VERSION, timestring(time.time())),
             '</center>',
             '</body></html>']
        
        return "".join(l)

    def html_page(self):
        title = 'TLSMD: Job Status'
        job_list = self.get_job_list()
	
        l = [self.html_head(title),
             html_title(title),
             html_nav_bar("queue"),
             self.html_private_form(),
             '<center><b>',
             'Or click on the Job ID you wish to view',
             '</b></center>',
             '<br>',
             self.html_running_job_table(job_list),
             '<br>',
             self.html_queued_job_table(job_list),
             '<br>',
             self.html_completed_job_table(job_list)]

        limbo = self.html_limbo_job_table(job_list)
        if limbo!=None:
            l.append('<br>')
            l.append(limbo)

        l.append(self.html_foot())
        
        return "".join(l)

    def html_private_form(self):
        l = ['<form action="webtlsmd.cgi" method="post">',
             '<input type="hidden" name="page" value="explore">',
             
             '<center>',
             '<b>To access a private job, enter its Job ID below</b>',
             '</center>',
             
             '<center>',
             '<input type="text" name="job_id" size="50">',
             '</center>',
             
             '</form>']

        return "".join(l)

    def explore_href(self, jdict):
        if self.admin==True:
            page = "admin"
        else:
            page = "explore"
 
        if self.admin==True:
            l = ['<a href="webtlsmd.cgi?page=%s&amp;job_id=%s">%s</a>' % (page, jdict["job_id"] ,jdict["job_id"])]

            if jdict.get("user_name", "")!="":
                l.append('<br>%s' % (jdict["user_name"]))

            if jdict.get("email", "")!="":
                l.append('<br>%s' % (jdict["email"]))

            return "".join(l)

        if jdict.get("private_job", False)==True:
            return 'private'
    
        return '<a href="webtlsmd.cgi?page=%s&amp;job_id=%s">%s</a>' % (page, jdict["job_id"] ,jdict["job_id"])
    
    def chain_size_string(self, jdict):
        if jdict.has_key("chains")==False:
            return "---"

        listx = []
        for cdict in jdict["chains"]:
            listx.append("%s:%d" % (cdict["chain_id"], cdict["length"]))

	strx = ''
	while len(listx)>0:
            l3 = listx[:5]
	    listx = listx[5:]

	    strx += " ".join(l3)
	    if len(listx)>0:
                strx += '<br>'
	
	return '%s' % (strx)
	
    def get_job_list(self):
        """Get a list of all the jobs in the job queue file.
        """
        return webtlsmdd.job_list()

    def html_running_job_table(self, job_list):
        ## get the job dictionary of the running job
        run_jdict = None
        for jdict in job_list:
            if jdict.get("state")=="running":
                run_jdict = jdict
                break

        ## if there is no running job
        jdict = run_jdict

        x  = '<center>'
	x += '<b>Running Jobs</b>'
        x += '<table border="0" cellpadding="3" width="100%" style="background-color:#eeeeee; font-size:x-small">'
        x += '<tr style="background-color:#bbbbbb">'
        x += '<th>Job ID</th>'
        x += '<th>Struct ID</th>'
	x += '<th>Chain:Num Res</th>'
        x += '<th>Submission Date</th>'
        x += '<th>Currently Processing</th>'
        x += '<th>Processing Time<br>Used (HH:MM.SS)</th>'
        x += '</tr>'

        if jdict!=None:
            x += '<tr>'

            x += '<td>%s</td>' % (self.explore_href(jdict))
            x += '<td>%s</td>' % (self.rcsb_href(jdict.get("structure_id", "XXXX")))
	    x += '<td>%s</td>' % (self.chain_size_string(jdict))
            x += '<td>%s</td>' % (timestring(jdict["submit_time"]))

            tls_seg = 'Chain <b>%s</b> %s%% Complete' % (jdict.get("run_chain_id", ""), jdict.get("pcomplete", ""))
            x += '<td>%s</td>' % (tls_seg)

            if jdict.has_key("run_start_time"):
                 hours = timediffstring(jdict["run_start_time"], time.time())
            else:
                 hours = "---"
            x += '<td align="right">%s</td>' % (hours)

            x += '</tr>'
        else:
	    x += '<tr>'
	    x += '<td colspan="6" align="center">'
	    x += 'No Jobs Running'
	    x += '</td>'
	    x += '</tr>'

        x += '</table>'
        x += '</center>'
        return x

    def html_queued_job_table(self, job_list):
        queued_list = []
        for jdict in job_list:
            if jdict.get("state")=="queued":
                queued_list.append(jdict)

        x  = ''
        x += '<center>'
	x += '<b>%d Queued Jobs</b>' % (len(queued_list))
        x += '<table border="0" cellpadding="3" width="100%" style="background-color:#eeeeee; font-size:x-small">'
        x += '<tr style="background-color:#bbbbbb">'
        x += '<th>Job ID</th>'
        x += '<th>Struct ID</th>'
	x += '<th>Chain:Num Res</th>'
        x += '<th>Submission Date</th>'
        x += '</tr>'

        for jdict in queued_list:
            x += '<tr>'
            
            x += '<td>%s</td>' % (self.explore_href(jdict))
            x += '<td>%s</td>' % (self.rcsb_href(jdict.get("structure_id", "XXXX")))
            x += '<td>%s</td>' % (self.chain_size_string(jdict))	
            x += '<td>%s</td>' % (timestring(jdict["submit_time"]))
                                  
            x += '</tr>'

        if len(queued_list)==0:
	    x += '<tr>'
	    x += '<td colspan="4" align="center">'
	    x += 'No Jobs Queued'
	    x += '</td>'
	    x += '</tr>'

        x += '</table>'
        x += '</center>'
        return x

    def html_completed_job_table(self, job_list):
        completed_list = []
        for jdict in job_list:
            if jdict["state"] in ["completed", "defunct"]:
                completed_list.append(jdict)

        completed_list.reverse()

        l = ['<center><b>%d Completed Jobs</b></center>' % (len(completed_list)),
             '<center>',
             '<table border="0" cellpadding="3" width="100%" style="background-color:#eeeeee; font-size:x-small">',
             '<tr style="background-color:#bbbbbb">',
             '<th>Job ID</th>',
             '<th>Struct ID</th>',
             '<th>Status</th>',
             '<th>Submission Date</th>',
             '<th>Processing Time<br> Used (HH:MM.SS)</th>',
             '</tr>']

        alt_color = True

        for jdict in completed_list:
            if alt_color:
                l.append('<tr style="background-color:#ddddee">')
            else:
                l.append('<tr>')
            alt_color = not alt_color
                                
            l.append('<td>%s</td>' % (self.explore_href(jdict)))
            l.append('<td>%s</td>' % (self.rcsb_href(jdict.get("structure_id", "XXXX"))))
            l.append('<td>%s</td>' % (jdict["state"]))
            l.append('<td>%s</td>' % (timestring(jdict["submit_time"])))

            if jdict.has_key("run_end_time") and jdict.has_key("run_start_time"):
                hours = timediffstring(jdict["run_start_time"], jdict["run_end_time"])
	    else:
		hours = "---"
            l.append('<td align="right">%s</td>' % (hours))

            l.append('</tr>')

        l.append('</table>')
        l.append('</center>')
        return "".join(l)
    
    def html_limbo_job_table(self, job_list):
        limbo_list = []
        for jdict in job_list:
            if jdict.get("state") not in ["queued", "running", "completed"]:
                limbo_list.append(jdict)

        if len(limbo_list)==0:
            return None

        x  = ''
        x += '<center>'
	x += '<b>Partially Submitted Jobs</b>'
        x += '<table border="0" width="100%" style="background-color:#eeeeee; font-size:x-small">'
        x += '<tr style="background-color:#bbbbbb">'
        x += '<th>Job ID</th>'
        x += '<th>Struct ID</th>'
        x += '<th>State</th>'
        x += '<th>Submission Date</th>'
        x += '</tr>'

        for jdict in limbo_list:
            x += '<tr>'

            x += '<td>%s</td>' % (self.explore_href(jdict))
            x += '<td>%s</td>' % (self.rcsb_href(jdict.get("structure_id", "XXXX")))
            x += '<td>%s</td>' % (jdict["state"])
            x += '<td>%s</td>' % (timestring(jdict["submit_time"]))
                                  
            x += '</tr>'

        x += '</table>'
        x += '</center>'
        return x
    

class ExploreJobPage(Page):
    def html_page(self):
        job_id = check_job_id(self.form, webtlsmdd)
	if job_id==None:
	    title = 'TLSMD: Explore Job'
	    x  = self.html_head(title)
	    x += html_title(title)
	    x += '<center><h3>ERROR: Invalid Job ID</h3></center>'
	    x += self.html_foot()
	    return x
	    
	title = 'TLSMD: Explore Job ID %s' % (job_id)
        x  = ''
	x += self.html_head(title)
	x += html_title(title)
	x += html_nav_bar()
	x += html_job_nav_bar(webtlsmdd, job_id)
        jdict = webtlsmdd.job_get_dict(job_id)
	x += html_job_info_table(jdict)
        x += self.html_foot()
	return x

       
class AdminJobPage(Page):
    def html_page(self):
        job_id = check_job_id(self.form, webtlsmdd)
	if job_id==None:
	    title = 'TLSMD: View Job'
	    x  = self.html_head(title)
	    x += html_title(title)
	    x += '<center><h3>ERROR: Invalid Job ID</h3></center>'
	    x += self.html_foot()
	    return x
        
        title = 'TLSMD: Administrate Job %s' % (job_id)

        x  = ''
        x += self.html_head(title)
        x += html_title(title)

        x += html_nav_bar()

        if self.form.has_key("submit") and self.form["submit"].value=="Remove Job":
            x += self.remove(job_id)
        else:
            x += self.edit(job_id)
        
        x += self.html_foot()
        return x

    def edit(self, job_id):
        x = ''

        ## if the job is not in the "queued" state, then it is not
        ## safe to edit
        state = webtlsmdd.job_data_get(job_id, "state")
        if state=="queued":
            extract_job_edit_form(self.form, webtlsmdd)

        ## get the state dictionary for the entire job
        fdict = webtlsmdd.job_get_dict(job_id)
        fdict["page"] = "admin"
        fdict["removebutton"] = True
            
        if state=="running" or state=="completed":
            x += html_job_nav_bar(webtlsmdd, job_id)
            x += html_job_info_table(fdict)
        else:
            x += html_job_edit_form(fdict)
        
        return x

    def remove(self, job_id):
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


SUBMIT1_NOTE = """\
Analysis of large structures is
computationally expensive, so you may have to wait hours to days for
the server to generate a complete analysis depending on how
heavily it is loaded.<br><br>
"""

class Submit1Page(Page):
    def html_page(self):
        title = 'TLSMD: Start a New Job'

        l = [self.html_head(title),
             html_title(title),
             '<center>',

             '<form enctype="multipart/form-data" action="webtlsmd.cgi" method="post">',
             '<input type="hidden" name="page" value="submit2">',

             '<table class="submit_table">',
             '<tr><th colspan="2" class="step_title">Step 1: Select your PDB file to upload</th></tr>',

             '<tr>',
             '<td align="left">Upload PDB File:</td>',
             '<td><input name="pdbfile" size="50" type="file"></td>',
             '</tr>',

             '<tr><td colspan="2" align="center">',
             '<input value="Upload File and Proceed to Step 2" type="submit">',
             '</td></tr>',

             '</table>',
             '</form>',
             
             '</center>',

             self.html_foot()]

        return "".join(l)


class Submit2Page(Page):
    def html_page(self):        
        title = 'TLSMD: Start a New Job'
        
        x  = ''
        x += self.html_head(title)
        x += html_title(title)

        try:
            job_id = self.prepare_submission()
        except SubmissionException, err:
             x += html_nav_bar()
             x += '<center><h3>ERROR:%s</h3></center>' % (err.html)
        else:            
            x += self.job_edit_form(job_id)

        x += self.html_foot()
        return x

    def job_edit_form(self, job_id, show_warnings=False):
        fdict = webtlsmdd.job_get_dict(job_id)
        fdict["page"] = "submit3"
        x = html_job_edit_form2(fdict, "Step 2: Fill out Submission Form, then Submit Job")
        return x

    def generate_security_code(self):
        l = list(5 * string.ascii_letters) 
        random.shuffle(l)
	code = "".join(random.sample(l,8)) 
	return code

    def prepare_submission(self):
        if self.form.has_key("pdbfile")==False or self.form["pdbfile"].file==None:
            raise SubmissionException("No PDB file uploaded")
	
        ## make working directory
        try:
            os.chdir(conf.TLSMD_WORK_DIR)
        except os.error, err:
            raise SubmissionException('<p>Cannot change to working directory: %s</p>' % (str(err)))

        job_id = webtlsmdd.job_new(self.generate_security_code())
        os.umask(022)
        try:
            os.mkdir(job_id)
        except os.error, err:
            webtlsmdd.job_delete(job_id)
            raise SubmissionException('<p>Cannot make directory: %s</p>' % (str(err)))

        job_dir = os.path.join(conf.TLSMD_WORK_DIR, job_id)
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
            raise SubmissionException('<p>Only Recieved %d lines</p>' % (num_lines))

        webtlsmdd.job_data_set(job_id, "pdb_filename", pdb_filename)

        job_url = "%s/%s" % (conf.TLSMD_WORK_URL, job_id)
        webtlsmdd.job_data_set(job_id, "job_url", job_url)

        log_url = "%s/log.txt" % (job_url)
        webtlsmdd.job_data_set(job_id, "log_url", log_url)

        analysis_dir      = "%s/ANALYSIS" % (job_dir)
        analysis_base_url = "%s/ANALYSIS" % (job_url)
        analysis_url      = "%s/ANALYSIS/index.html" % (job_url)
        webtlsmdd.job_data_set(job_id, "analysis_dir", analysis_dir)
        webtlsmdd.job_data_set(job_id, "analysis_base_url", analysis_base_url)
        webtlsmdd.job_data_set(job_id, "analysis_url", analysis_url) 

        ip_addr = os.environ.get("REMOTE_ADDR", "Unknown")
        webtlsmdd.job_data_set(job_id, "ip_addr", ip_addr)

        ## submission time and initial state
        webtlsmdd.job_data_set(job_id, "state", "submit1")
        webtlsmdd.job_data_set(job_id, "submit_time", time.time())

        ## now load the structure and build the submission form
        struct = FileLoader.LoadStructure(fil=pdb_filename)
	if not struct.structure_id:
	    struct.structure_id = "XXXX"

        webtlsmdd.job_data_set(job_id, "structure_id", struct.structure_id)

        ## Select Chains for Analysis
        num_atoms          = 0
        num_aniso_atoms    = 0
        largest_chain_seen = 0

        chains = []

        for chain in struct.iter_chains():
            naa = chain.count_amino_acids()
            if naa<10:
                continue

            largest_chain_seen = max(naa, largest_chain_seen)

            for atm in chain.iter_all_atoms():
                num_atoms += 1
                if atm.U!=None:
                    num_aniso_atoms += 1

            ## form name
            cb_name = 'CHAIN%s' % (chain.chain_id)
            
            ## create chain description label cb_desc
            cb_desc = 'Chain %s (%d Amino Acid Residues)' % (chain.chain_id, chain.count_amino_acids())

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
            cdict["length"]   = naa
            cdict["name"]     = cb_name
            cdict["desc"]     = cb_desc
            cdict["preview"]  = cb_preview
            cdict["selected"] = True

        if largest_chain_seen>1700:
            webtlsmdd.job_delete(job_id)
	    strx = '<p>Your submitted structure contained a chain exceeding the 1700 residue limit</p>'
            raise SubmissionException(strx)

        webtlsmdd.job_data_set(job_id, "chains", chains)

        ## defaults
        webtlsmdd.job_data_set(job_id, "user", "")
        webtlsmdd.job_data_set(job_id, "passwd", "")
        webtlsmdd.job_data_set(job_id, "email", "")
        webtlsmdd.job_data_set(job_id, "comment", "")
        webtlsmdd.job_data_set(job_id, "private_job", False)
        webtlsmdd.job_data_set(job_id, "plot_format", "PNG")

        aniso_ratio = float(num_aniso_atoms)/float(num_atoms)
        if aniso_ratio>0.90:
            webtlsmdd.job_data_set(job_id, "tls_model", "ANISO")
        else:
            webtlsmdd.job_data_set(job_id, "tls_model", "ISOT")
            
        webtlsmdd.job_data_set(job_id, "weight", "")
        webtlsmdd.job_data_set(job_id, "include_atoms", "ALL")

	return job_id

SUBMIT3_CAP1 = """\
You may monitor the progress of your TLSMD submission by its Job ID
on the Job Status page, available by clicking the link on the top
of this page.  All queued, running and completed jobs are listed on
the Job Status page.  Through this page you may explore the output
of your job, and lookup your job by its Job ID if you have chosen to
keep your job private.
"""

class Submit3Page(Page):
    def html_page(self):
        try:
            job_id = self.complete_submission()
	except SubmissionException, err:
	    title = 'TLSMD: Job Submission Failed'
            html  = '<center><h3>ERROR: %s</h3></center>' % (err.html)
	else:
            title = 'TLSMD: Job Submission Succeeded'

            l = ['<center>',
                 
                 '<table class="submit_table">',
                 '<tr><th class="step_title">Step 3: Finished!  Job successfully submitted.</th></tr>',

                 '<tr><td align="center">Your job ID is %s</td></tr>' % (job_id),

                 '<tr><td>',
                 '<p>Visit and bookmark your ',
                 '<a href="webtlsmd.cgi?page=explore&amp;job_id=%s">Explore Job %s</a> ' % (job_id, job_id),
                 'page, this page is the status page of your job, and it is ',
                 'updated as your job progresses through the queue.  Once your ',
                 'job is complete, a link to the completed TLSMD analysis will appear ',
                 'on it.',
                 '</p>',
                 
                 '<p>%s</p>' % (SUBMIT3_CAP1),

                 '</td></tr>',
                 '</table>',
                 '</center>']
            
            html = "".join(l)
	    
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
        if self.form["submit"].value=="Cancel Job Submission":
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

    if form.has_key("page"):

        if form["page"].value=="explore":
            page = ExploreJobPage(form)

        elif form["page"].value=="admin":
            page = AdminJobPage(form)

        elif form["page"].value=="submit1":
            page = Submit1Page(form)

        elif form["page"].value=="submit2":
            page = Submit2Page(form)

        elif form["page"].value=="submit3":
            page = Submit3Page(form)

    if page==None:
        page = QueuePage(form)

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
