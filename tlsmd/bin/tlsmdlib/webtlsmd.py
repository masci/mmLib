# -*- coding: utf-8 -*-
## TLS Minimized Domains (TLSMD)
## Copyright 2002-2008 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
##
## NOTE: This file contains changes and additions by Christoph Champ, 2007-10-24

import os
import sys
import time
import socket
import string
import random
import math

import xmlrpclib
import cgitb; cgitb.enable()
import cgi

import const
import conf
import numpy		# Added by Christoph Champ, 2007-10-23
import subprocess	# Added by Christoph Champ, 2007-11-20
import datetime		# Added by Christoph Champ, 2008-02-01
import re		# Added by Christoph Champ, 2008-02-07

## GLOBALS
webtlsmdd = xmlrpclib.ServerProxy(conf.WEBTLSMDD)

def timestring(secs):
    tm_struct = time.localtime(secs)
    ## Switched to international time format. Christoph Champ, 2008-02-07
    #return time.strftime("%m-%d-%y %H:%M %Z" ,tm_struct)
    return time.strftime("%Y-%m-%d %H:%M %Z" ,tm_struct)

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
    l = ['<div id="navcontainer">',
         '  <ul>',
	 '    <li><a href="%s/index.html">Home</a></li>' % (conf.TLSMD_BASE_URL),
         '    <li><a href="webtlsmd.cgi?page=submit1">Start a New Job</a></li>', 
	 '    <li><a href="webtlsmd.cgi">Job Status</a></li>',
	 '    <li><a href="%s/examples/index.html">Examples</a></li>' % (conf.TLSMD_BASE_URL),
	 '    <li><a href="%s/documentation.html">Documentation</a></li>' % (conf.TLSMD_BASE_URL),
	 '  </ul>',
	 '</div>'
         ]
    return "\n".join(l)

def html_job_nav_bar(webtlsmdd, job_id):
    """Navigation bar to the TLSMD output files.
    """
    analysis_dir = webtlsmdd.job_get_analysis_dir(job_id)
    analysis_index = os.path.join(analysis_dir, "index.html")
    analysis_url = webtlsmdd.job_get_analysis_url(job_id)

    job_dir = webtlsmdd.job_get_job_dir(job_id)
    logfile = os.path.join(job_dir, "log.txt")
    log_url = webtlsmdd.job_get_log_url(job_id)

    ## tarball. Christoph Champ, 2007-12-03
    #tarball_url = os.path.join(job_dir,"%s.tar.bz2"%job_id)
    tarball_url = webtlsmdd.job_get_tarball_url(job_id)

    if not os.path.isfile(analysis_index) and not os.path.isfile(logfile):
        return ''

    x  = ''
    x += '<center>'

    if os.path.isfile(analysis_index):
        x += '<h3>View <a href="%s">Completed TLSMD Analysis</a></h3>' % (analysis_url)

    if os.path.isfile(logfile):
        x += '<h3>View <a href="%s">TLSMD Logfile</a></h3>' % (log_url)

    ## tarball link. Christoph Champ, 2007-12-03
    if os.path.isfile(analysis_index):
        x += '<h3>Download <a href="%s">Local Copy of TLSMD Analysis output (tarball)</a></h3>' % (tarball_url)

    x += '</center>'
    x += '<br>'
    return x


def html_job_edit_form(fdict, pdb=False):
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
    if not pdb:
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
    if not pdb:
        x += '<tr>'
        x += '<td align="right"><label>Structure Code:</td><td>'
        x += '<input disabled type="text" name="structure_id" value="%s" size="4" maxlength="4">' % (fdict.get("structure_id", ""))
        x += '</label></td>'
        x += '</tr>'

        x += '</td>'
    
    x += '</table>'

    ## session info
    x += '<td valign="top"><table>'

    x += '<tr><td align="right">TLSMD Job ID:</td>'
    x += '<td><b>%s</b></td></tr>' % (fdict["job_id"])

    x += '<tr><td align="right">Job State:</td>'
    try:
	x += '<td><b>%s</b></td></tr>' % (fdict["state"])
    except:
	x += '<td><b>None</b></td></tr>'
    
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
    if not pdb:
        x += '<tr><th colspan="3">Select Chains for Analysis</th></tr>'

        x += '<tr><td colspan="3">'
        x += '<table>'
        for cdict in fdict.get("chains", []):
            x += '<tr><td>'
            x += '<label>'
            if cdict["selected"]:
                x += '<input type="checkbox" name="%s" value="TRUE" checked>' % (cdict["name"])
            else:
                x += '<input type="checkbox" name="%s" value="TRUE">' % (cdict["name"])
            x += '%s' % (cdict["desc"])
            x += '</label>'

            x += '</td></tr>'

        x += '</table></td></tr>'
    else:
        # select all the chains by default
        for cdict in fdict.get("chains", []):
            x += '<input type="hidden" name="%s" value="TRUE">' % (cdict['name'])

    x += '</table>'
    ## end form

    x += '<tr><td colspan="3">'

    x += '<table width="100%">'
    x += '<tr>'
    
    x += '<td align="left">'
    if fdict.has_key("removebutton"):
        x += '<input type="submit" name="submit" value="Remove Job">'
    if fdict.has_key("requeuebutton"):
        x += '<input type="submit" name="submit" value="Requeue Job">'
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

	 ## New user_comment added. Christoph Champ, 2007-12-18
         '<tr>',
         '<td align="right"><label for="user_comment">Associated Notes</label></td>',
         '<td><input type="text" id="user_comment" name="user_comment" value="%s" size="40" maxlength="128"></td>' % (fdict.get("user_comment","")),
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

         ## Changed default to 'private'. Christoph Champ, 2007-12-18
         '<tr><td>',
         '<label><input type="checkbox" id="private_job" name="private_job" value="TRUE" checked>Keep Job Private</label>',
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
        if cdict["selected"]:
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

    ## user comments. Christoph Champ, 2007-12-18
    x += '<tr>'
    x += '<td align="right"><label>Associated Notes:</td><td>'
    x += '<b>%s</b>' % (fdict.get("user_comment", ""))
    x += '</label></td>'
    x += '</tr>'

    x += '</table>'
    x += '</td>'

    ## session info
    x += '<td valign="top"><table>'

    x += '<tr><td align="right">TLSMD Job ID:</td>'
    x += '<td><b>%s</b></td></tr>' % (fdict["job_id"])

    x += '<tr><td align="right">Job State:</td>'
    if fdict.has_key("state"):
        jobstate = (fdict["state"])
    else:
        jobstate = "unknown"
    x += '<td><b>%s</b></td></tr>' % (jobstate)
    
    x += '<tr><td align="right">Submission IP Address: </td>'
    x += '<td><b>%s</b></td></tr>' % (fdict.get("ip_addr", ""))

    x += '<tr><td align="right">Submission Date: </td>'
    if fdict.has_key("submit_time"):
        date = timestring(fdict["submit_time"])
    else:
        date = "---"
    x += '<td><b>%s</b></td></tr>' % (date)

    x += '<tr><td align="right">Processing Start Date: </td>'
    if fdict.has_key("run_time_begin"):
        date = timestring(fdict["run_time_begin"])
    else:
        date = "---"
    x += '<td><b>%s</b></td></tr>' % (date)

    x += '<tr><td align="right">Processing End Date: </td>'
    if fdict.has_key("run_time_end"):
        date = timestring(fdict["run_time_end"])
    else:
        date = "---"
    x += '<td><b>%s</b></td></tr>' % (date)

    x += '<tr><td align="right">Processing Time(HH:MM): </td>'
    if fdict.has_key("run_time_end") and fdict.has_key("run_time_begin"):
        hours = timediffstring(fdict["run_time_begin"], fdict["run_time_end"])
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
        if cdict["selected"]:
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
    if fdict.get("tls_model") is None or fdict.get("tls_model")=="ISOT":
        x += 'Isotropic'
    elif fdict.get("tls_model")=="ANISO":
        x += 'Anisotropic'
    x += '</td>'

    ## Least Squares Weighting
    x += '<td>'
    if fdict.get("weight") is None or fdict.get("weight")=="IUISO":
        x += 'Inverse Atomic B<sub>iso</sub>'
    elif fdict.get("weight")=="NONE":
        x += 'No Weighting'
    x += '</td>'

    ## Include Atoms
    x += '<td>'
    if fdict.get("include_atoms") is None or fdict.get("include_atoms")=="ALL":
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
    if isinstance(data, unicode):
        return False
    if len(data) > max_len:
        return False
    if not data.isalnum():
        return False
    return True

def vet_email(email_address):
    if len(email_address) > 45:
        return False
    return True


def extract_job_edit_form(form, webtlsmdd):
    """Extract the input from the Job Edit Form and update the webtlsmdd
    database with the information.
    """
    if not form.has_key("edit_form"):
        return False

    job_id = check_job_id(form, webtlsmdd)
    if job_id is None:
        return False

    if form.has_key("user"):
        user = form["user"].value.strip()
        if vet_data(user, 10):
            webtlsmdd.job_set_user(job_id, user)

    if form.has_key("private_job"):
        webtlsmdd.job_set_private_job(job_id, True)

    if form.has_key("user_name"):
        user_name = form["user_name"].value.strip()
        user_name = user_name[:100]
        webtlsmdd.job_set_user_name(job_id, user_name)
            
    if form.has_key("email"):
        email_address = form["email"].value.strip()
        if vet_email(email_address):
            webtlsmdd.job_set_email(job_id, email_address)

    if form.has_key("structure_id"):
        structure_id = form["structure_id"].value.strip()
        if vet_data(structure_id, 4):
            webtlsmdd.job_set_structure_id(job_id, structure_id)

    ## New user_comment field added. Christoph Champ, 2007-12-18
    if form.has_key("user_comment"):
        user_comment = form["user_comment"].value.strip()
        user_comment = user_comment[:128]
        webtlsmdd.job_set_user_comment(job_id, user_comment)

    chains = webtlsmdd.job_get_chains(job_id)
    for cdict in chains:
        if form.has_key(cdict["name"]):
            cdict["selected"] = True
        else:
            cdict["selected"] = False
    webtlsmdd.job_set_chains(job_id, chains)

    if form.has_key("tls_model"):
        tls_model = form["tls_model"].value.strip()
        if tls_model in ["ISOT", "ANISO"]:
            webtlsmdd.job_set_tls_model(job_id, tls_model)

    if form.has_key("weight"):
        weight = form["weight"].value.strip()
        if weight in ["NONE", "IUISO"]:
            webtlsmdd.job_set_weight_model(job_id, weight)

    if form.has_key("include_atoms"):
        include_atoms = form["include_atoms"].value.strip()
        if include_atoms in ["ALL", "MAINCHAIN"]:
            webtlsmdd.job_set_include_atoms(job_id, include_atoms)

    if form.has_key("plot_format"):
        plot_format = form["plot_format"].value.strip()
        if plot_format in ["PNG", "SVG"]:
            webtlsmdd.job_set_plot_format(job_id, plot_format)

    return True



class Page(object):
    def __init__(self, form):
        self.form = form

    def html_head_nocgi(self, title, redirect=None):
        x  = ''
        x += '<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">'
        x += '<html>'
        x += '<head>'
        x += '  <title>%s</title>' % (title)
        x += '  <link rel="stylesheet" href="../tlsmd.css" type="text/css" media="screen">'
        x += '  <link rel="stylesheet" href="../tlsmd_print.css" type="text/css" media="print">'
        if redirect != None:
            x += '<meta http-equiv="REFRESH" content="3; URL=%s">' % (redirect)
        x += '</head>'
        x += '<body><div id="page">'
        return x

    def html_head(self, title, redirect=None):
        if redirect == None:
            return 'Content-Type: text/html\n\n' + self.html_head_nocgi(title)
        else:
            return 'Content-Type: text/html\n\n' + self.html_head_nocgi(title, redirect)

    def html_foot(self):
        l = ['<center>',
             '<p><small><b>Version %s</b> Last Modified %s' % (const.VERSION, const.RELEASE_DATE),
             '</small></p>',
             '</center>',
             '</div></body></html>']
        
        return "".join(l)


class ErrorPage(Page):
    def __init__(self, form, text=None):
        Page.__init__(self, form)
        self.text = text
    
    def html_page(self):
        title = 'TLSMD: Error'

        l = [self.html_head(title, None),
             html_title(title),
             html_nav_bar(),
             '<br>',
             '<center><p class="perror">Error<br>' ]
        
        if self.text is not None:
            l.append(self.text)

        l.append('</p></center>')

        l.append(self.html_foot())
        return "".join(l)


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
        return code == passcode

    def rcsb_href(self, jdict):
        if jdict.get("private_job", False):
            return "----"
        struct_id = jdict.get("structure_id", "xxxx")
        if struct_id.lower() == "xxxx":
            return struct_id
        return '<a href="http://www.pdb.org/pdb/explore.do?structureId=%s">%s</a>' % (struct_id, struct_id)

    def html_head_nocgi(self, title):
        l = ['<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">',
             '<html>',
             '<head>',
             '  <title>%s</title>' % (title),
	     '  <link rel="stylesheet" href="../tlsmd.css" type="text/css" media="screen">',
	     '  <link rel="stylesheet" href="../tlsmd_print.css" type="text/css" media="print">',
             '</head>',
             '<body><div id="page">']
        
        return "".join(l)
    
    #'<p><small><b>Version %s</b> Last Updated %s</p>' % (const.VERSION, timestring(time.time())),
    def html_foot(self):
        l = ['<center>',
             '<p><small><b>Version %s</b> Last Updated %s PST</p>' % (const.VERSION, (datetime.datetime.fromtimestamp(time.time()).isoformat(' ')[:-10])),
             '</center>',
             '</div></body></html>']
        
        return "".join(l)

    def html_page(self):
        title = 'TLSMD: Job Status'
        job_list = self.get_job_list()
	
        l = [self.html_head(title, None),
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
        if self.admin:
            page = "admin"
        else:
            page = "explore"
 
        if self.admin:
            job_id = jdict.get("job_id", "")

            user_name = jdict.get("user_name", "")
            if isinstance(user_name, unicode):
                user_name = ""

            ## New user_comment added. Christoph Champ, 2007-12-18
            user_comment = jdict.get("user_comment", "")
            if isinstance(user_comment, unicode):
                user_comment = ""

	    email_address = jdict.get("email")
            if isinstance(email_address, unicode):
                email_address = ""
            
            l = ['<a href="webtlsmd.cgi?page=%s&amp;job_id=%s">%s</a>' % (page, job_id, job_id)]

            if user_name != "":
                l.append('<br>%s' % (user_name))

            ## New user_comment added. Christoph Champ, 2007-12-18
            if user_comment != "":
                l.append('<br>%s' % (user_comment))

            if email_address != "":
                l.append('<br>%s' % (email_address))

            return "".join(l)

        if jdict.get("private_job", False):
	    ## Return job number only (non-clickable)
            job_number = re.match(r'[^_]*', jdict["job_id"])
            if job_number: return job_number.group(0)
	    return 'private'
    
        return '<a href="webtlsmd.cgi?page=%s&amp;job_id=%s">%s</a>' % (page, jdict["job_id"] ,jdict["job_id"])
    
    def chain_size_string(self, jdict):
        if jdict.has_key("chains") == False:
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

#============================================================================================
    def get_progress_status(self, job_id):
	### FIXME Too basic right now. Christoph Champ, 2008-02-07 ###
        ### Open directory, loop over filenames, and check for certain files ###
	prog="0" ## Start out at 0%
	try:
	     os.chdir(conf.TLSMD_WORK_DIR + '/' + job_id + '/ANALYSIS')
	except:
	#    ## Immediately fails; ANALYSIS dir is not written right away
	#    webtlsmdd.job_set_state(job_id,"lost_directory") ## cannot have spaces in state
	     return prog
	pngfile=re.compile('\S+NTLS1_BMEAN.png').match		## match first *NTLS1_BMEAN.png file found in the job_dir
	reffile=re.compile('\S+REFINEMENT_PREP.html').match	## match first *REFINEMENT_PREP.html file found in the job_dir
	for fname in os.listdir(os.getcwd()):
            if pngfile(fname): prog="25" ## arbitrary 60% done
	    elif reffile(fname): prog="95"
	    #else: prog="0"
        return prog
#============================================================================================

    def html_running_job_table(self, job_list):

	## get an array of "running" jobs from the job dictionary. Christoph Champ, 2008-01-30
        run_jdict = []
        for jdict in job_list:
            if jdict.get("state") == "running":
		run_jdict.append(jdict)

	x  = ['<center>',
	      '<b>%d Running Jobs</b>' % (len(run_jdict)),
	      '<table border="0" cellpadding="3" width="100%" class="status_table">',
	      '<tr class="status_table_head">',
	      '<th>Job ID</th>',
	      '<th>Structure ID</th>',
	      '<th>Chain:Num Res</th>',
	      '<th>Submission Date</th>',
	      '<th colspan="2">Running Time (HH:MM.SS)</th>',
	      '</tr>']

	## creates mutiple rows, _if_ there are multiple "running" jobs. Christoph Champ, 2008-01-30
	## TODO Update the progress bar to reflect where we are in the "running" state
        row1 = True
        #progress = 0 ## Isn't this overwritten immediately?
        for jdict in run_jdict:
            if row1:
                x.append('<tr class="status_table_row1">')
            else:
                x.append('<tr class="status_table_row2">')
            row1 = not row1

            x += ['<td>%s</td>' % (self.explore_href(jdict)),
                  '<td>%s</td>' % (self.rcsb_href(jdict)),
                  '<td>%s</td>' % (self.chain_size_string(jdict)),
                  '<td>%s</td>' % (timestring(jdict["submit_time"]))]

            if jdict.has_key("run_time_begin"):
                 hours = timediffstring(jdict["run_time_begin"], time.time())
            else:
                 hours = "---"
            ## progress bar should be something like (number of files created) / (20*number of chains)
            #progress = 25
	    progress = self.get_progress_status(jdict["job_id"]) ## Returns a integer value
	    ## Testing with "%"-done indicator
	    ## TODO "% done" Doesn't look to good in style/format. Christoph Champ, 2008-02-12
            x += '<td align=left><div class="prog-border"><div class="prog-bar" style="width: %s%%;"></div></div></td>' %progress
            #x += '<td><div class="prog-border"><div class="prog-bar" style="width: %s%%;"></div></div> %s%%</td>' %(progress,progress)
	    x += '<td align="right">%s</td></tr>' % (hours)

	## for zero running jobs
        if len(run_jdict) == 0:
            x += ['<tr>',
                  '<td colspan="6" align="center">',
                  'No Jobs Running',
                  '</td>',
                  '</tr>']

	#x += '</table></center>' ## Old
	x.append('</table></center>')
	return "".join(x)

    def html_queued_job_table(self, job_list):
        queued_list = []
        for jdict in job_list:
            if jdict.get("state") == "queued":
                queued_list.append(jdict)

        l = ['<center>',
             '<b>%d Queued Jobs</b>' % (len(queued_list)),
             '<table border="0" cellpadding="3" width="100%" class="status_table">',
             '<tr class="status_table_head">',
             '<th>Job ID</th>',
             '<th>Struct ID</th>',
             '<th>Chain:Num Res</th>',
             '<th>Submission Date</th>',
             '</tr>']

        row1 = True
        for jdict in queued_list:
            if row1:
                l.append('<tr class="status_table_row1">')
            else:
                l.append('<tr class="status_table_row2">')
            row1 = not row1

            l += ['<td>%s</td>' % (self.explore_href(jdict)),
                  '<td>%s</td>' % (self.rcsb_href(jdict)),
                  '<td>%s</td>' % (self.chain_size_string(jdict)),
                  '<td>%s</td>' % (timestring(jdict["submit_time"])),
                  '</tr>' ]

        if len(queued_list) == 0:
	    l += ['<tr>',
                  '<td colspan="4" align="center">',
                  'No Jobs Queued',
                  '</td>',
                  '</tr>']

        l.append('</table></center>')
        
        return "".join(l)

    def html_completed_job_table(self, job_list):
        completed_list = []
        for jdict in job_list:
            if jdict.get("state") in ["completed", "defunct"]:
                completed_list.append(jdict)

        completed_list.reverse()

        l = ['<center><b>%d Completed Jobs</b></center>' % (len(completed_list)),
             '<center>',
             '<table border="0" cellpadding="3" width="100%" class="status_table">',
             '<tr class="status_table_head">',
             '<th>Job ID</th>',
             '<th>Struct ID</th>',
             '<th>Status</th>',
             '<th>Submission Date</th>',
             '<th>Processing Time (HH:MM.SS)</th>',
             '</tr>']

        row1 = True
        for jdict in completed_list:
            if row1:
                l.append('<tr class="status_table_row1">')
            else:
                l.append('<tr class="status_table_row2">')
            row1 = not row1
                                
            l.append('<td>%s</td>' % (self.explore_href(jdict)))
            l.append('<td>%s</td>' % (self.rcsb_href(jdict)))
            l.append('<td>%s</td>' % (jdict.get("state")))
            l.append('<td>%s</td>' % (timestring(jdict["submit_time"])))

            if jdict.has_key("run_time_begin") and jdict.has_key("run_time_end"):
                hours = timediffstring(jdict["run_time_begin"], jdict["run_time_end"])
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

        if len(limbo_list) == 0:
            return None

        x  = ''
        x += '<center>'
	x += '<b>Partially Submitted Jobs</b>'
        x += '<table border="0" width="100%" class="status_table">'
        x += '<tr class="status_table_head">'
        x += '<th>Job ID</th>'
        x += '<th>Struct ID</th>'
        x += '<th>State</th>'
        x += '<th>Submission Date</th>'
        x += '</tr>'

        for jdict in limbo_list:
            x += '<tr>'

            x += '<td>%s</td>' % (self.explore_href(jdict))
            x += '<td>%s</td>' % (self.rcsb_href(jdict))
            x += '<td>%s</td>' % (jdict.get("state"))
            x += '<td>%s</td>' % (timestring(jdict.get("submit_time")))
                                  
            x += '</tr>'

        x += '</table>'
        x += '</center>'
        return x
    

class ExploreJobPage(Page):
    def html_page(self):
        job_id = check_job_id(self.form, webtlsmdd)
	if job_id is None:
	    title = 'TLSMD: Explore Job'
	    x  = self.html_head(title, None)
	    x += html_title(title)
	    x += '<center><p class="perror">ERROR: Invalid Job ID</p></center>'
	    x += self.html_foot()
	    return x
	    
	title = 'TLSMD: Explore Job ID %s' % (job_id)
        x  = ''
	x += self.html_head(title, None)
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
        jdict = webtlsmdd.job_get_dict(job_id)
        pdb = jdict.get('via_pdb', False)

	if job_id is None:
	    title = 'TLSMD: View Job'
	    x  = self.html_head(title, None)
	    x += html_title(title)
	    x += '<center><p class="perror">ERROR: Invalid Job ID</p></center>'
	    x += self.html_foot()
	    return x
        
        title = 'TLSMD: Administrate Job %s' % (job_id)

        x  = ''
        x += self.html_head(title, None)
        x += html_title(title)

        x += html_nav_bar()

        if self.form.has_key("submit") and self.form["submit"].value == "Remove Job":
            x += self.remove(job_id)

        elif self.form.has_key("submit") and self.form["submit"].value == "Requeue Job":
            x += self.requeue(job_id)
        else:
            x += self.edit(job_id, pdb)
        
        x += self.html_foot()
        return x

    def edit(self, job_id, pdb):
        x = ''

        ## if the job is not in the "queued" state, then it is not safe to edit
        state = webtlsmdd.job_get_state(job_id)
        if state == "queued":
            extract_job_edit_form(self.form, webtlsmdd)

        ## get the state dictionary for the entire job
        fdict = webtlsmdd.job_get_dict(job_id)
        fdict["page"] = "admin"
        fdict["removebutton"] = True
        fdict["requeuebutton"] = True
            
        if state == "running" or state == "completed":
            x += html_job_nav_bar(webtlsmdd, job_id)
            x += html_job_info_table(fdict)
        else:
            x += html_job_edit_form(fdict, pdb)
        
        return x

    def remove(self, job_id):
        webtlsmdd.remove_job(job_id)
        x  = ''
        x += '<center>'
        x += '<h3>Job %s has been removed.</h3>' % (job_id)
        x += '</center>'
        return x

    def requeue(self, job_id):
        result = webtlsmdd.requeue_job(job_id)
        x  = ''
        x += '<center>'
        if result:
            x += "<h3>Job %s has been pushed to the back.</h3>" % (job_id)
        else:
            x += "<h3>Job %s could not be requeued because it is running.</h3>" % (job_id)
        x += '</center>'
        return x


class SubmissionException(Exception):
    def __init__(self, err):
        Exception.__init__(self)
        self.err = err

    def __str__(self):
        return self.err


SUBMIT1_NOTE = """\
Analysis of large structures is
computationally expensive, so you may have to wait hours to days for
the server to generate a complete analysis depending on how
heavily it is loaded.<br><br>
"""

class Submit1Page(Page):
    def html_page(self):
        title = 'TLSMD: Start a New Job'

        ## Added html_nav_bar(). Christoph Champ, 2007-12-03
        l = [self.html_head(title, None),
             html_title(title),
             html_nav_bar(),
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
             
             '<center><h4>OR</h4></center>',

             '<center>',

             '<form action="webtlsmd.cgi" method="post">',
             '<input type="hidden" name="page" value="submit_pdb">',
             '<table class="submit_table">',
             '<tr><th colspan="2" class="step_title">Enter a PDB ID:</th>',
             '<td><input name="pdbid" size="4" maxlength="4" type="text"></td>',
             '<td><input value="Submit" type="submit"</td>',
             '</tr>',
             '</center>',
             '</table>',
             '<br><i><font color=red>TLSMD requires crystallographically refined B factors.',
             '<br>Please do not submit NMR structures, theoretical models, ',
             '<br>or any PDB file with unrefined Bs',
             '</font></i>',

             self.html_foot()]

        return "".join(l)


class Submit2Page(Page):

    def html_page(self):        
        title = 'TLSMD: Start a New Job'
        
        ## Added html_nav_bar(). Christoph Champ, 2007-12-03
        l = [self.html_head(title, None),
             html_title(title),html_nav_bar() ]

        try:
            job_id = self.prepare_submission()
        except SubmissionException, err:
             l.append('<center><p class="perror">ERROR:<br>%s</p></center>' % (err))
        else:            
            l.append(self.job_edit_form(job_id))

        l.append(self.html_foot())
        return "".join(l)

    def job_edit_form(self, job_id, show_warnings = False):
        fdict = webtlsmdd.job_get_dict(job_id)
        fdict["page"] = "submit3"
        return html_job_edit_form2(fdict, "Step 2: Fill out Submission Form, then Submit Job")

    def prepare_submission(self):
        if self.form.has_key("pdbfile") == False or self.form["pdbfile"].file is None:
            raise SubmissionException("No PDB file uploaded")

        ## allocate a new JobID
        job_id = webtlsmdd.job_new()
        ip_addr = os.environ.get("REMOTE_ADDR", "Unknown")
        webtlsmdd.job_set_remote_addr(job_id, ip_addr)
        
        ## save PDB file
        infil = self.form["pdbfile"].file
        line_list = []
        while True:
            ln = infil.readline()
            if not ln:
                break
            line_list.append(ln)

        ## error out if there weren't many lines
        if len(line_list) < 10:
            webtlsmdd.remove_job(job_id)
            raise SubmissionException('Only Recieved %d lines of upload' % (len(line_list)))

        ## basic sanity checks
        r = check_upload(line_list)
        if r != '':
            raise SubmissionException(str(r))

        ## pass the PDB file to the application server
        result = webtlsmdd.set_structure_file(job_id, xmlrpclib.Binary("".join(line_list)))
        if result != "":
            raise SubmissionException(result)

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
            html  = '<center><p class="perror">ERROR:<br>%s</p></center>' % (err)
	else:
            title = 'TLSMD: Job Submission Succeeded'

            l = ['<center>',
                 
                 '<table class="submit_table">',
                 '<tr><th class="step_title">Step 3: Finished!  Job successfully submitted.</th></tr>',

                 '<tr><td align="center">Your job ID is <B>%s</B></td></tr>' % (job_id),

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
	    
        x  = self.html_head(title, None)
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
        if job_id is None:
            raise SubmissionException('Submission Error')

        ## make sure the job is in the right state to be submitted
	state = webtlsmdd.job_get_state(job_id)
	if state == "queued":
	    raise SubmissionException("Your job is already queued")
    	elif state == "running":
	    raise SubmissionException("Your job is already running")

        ## verify the submission IP address
        ip_addr = os.environ.get("REMOTE_ADDR", "Unknown")
        ip_addr_verify = webtlsmdd.job_get_remote_addr(job_id)
        if ip_addr != ip_addr_verify:
            raise SubmissionException('Submission IP Address Mismatch')

        ## completely remove the job
        if self.form["submit"].value == "Cancel Job Submission":
            webtlsmdd.remove_job(job_id)
            raise SubmissionException('You cancelled the job')

        extract_job_edit_form(self.form, webtlsmdd)

        ## if everything with the form is okay, then change
        ## the job state to queued
        webtlsmdd.job_set_state(job_id, "queued")

        return job_id

class SubmitPDBPage(Page):
    """Handles requests submitted via a PDB ID"""

    def html_page(self):
        if "pdbid" not in self.form:
            raise SubmissionException("Please enter a PDB ID")
        elif len(self.form["pdbid"].value) < 4 or not self.form["pdbid"].value.isalnum():
            raise SubmissionException("Invalid PDB ID.  Please try again.")

        pdbid = self.form["pdbid"].value.upper()
         
        if webtlsmdd.pdb_exists(pdbid):
            return self.redirect_page(pdbid)

        pdbfile_bin = webtlsmdd.fetch_pdb(pdbid)
        pdbfile = pdbfile_bin.data
        
        if len(pdbfile) == 0:
            raise SubmissionException("Could not download PDB File from RCSB.")

        ## basic sanity checks. Christoph Champ, 2007-10-24
        ln=pdbfile.split("\n")
        r = check_upload(ln)
        if r != '':
            raise SubmissionException(str(r))
 
        job_id = self.prepare_submission(pdbfile)

        if not webtlsmdd.set_pdb_db(pdbid):
            raise SubmissionException("Could not write to internal PDB DB")

        webtlsmdd.job_set_pdb_dir(job_id, pdbid)
        
        fdict = webtlsmdd.job_get_dict(job_id)
        fdict["page"] = "submit3"
        
        title = "Enter contact info:"
        l = [self.html_head(title, None), html_title(title)]
        l.append(html_job_edit_form(fdict, pdb=True))
        l.append(self.html_foot())

        return "".join(l)

    def prepare_submission(self, pdbfile):
        job_id = webtlsmdd.job_new()
        ip_addr = os.environ.get("REMOTE_ADDR", "Unknown")
        webtlsmdd.job_set_remote_addr(job_id, ip_addr)
        webtlsmdd.job_set_via_pdb(job_id, True)
        result = webtlsmdd.set_structure_file(job_id, xmlrpclib.Binary(pdbfile))
        if result != "":
            return SubmissionException("Failed to submit structure.")
        return job_id

    def redirect_page(self, pdbid):
        # check to see if this job is still running
        try:
            os.chdir(conf.WEBTLSMDD_PDB_DIR + '/' + pdbid)
        except OSError:
            title = "This structure is currently being analyzed, please check back later."
            page = [self.html_head(title),
                    html_title(title),
                    self.html_foot()]
            return "".join(page)


        title = "This protein has already been analyzed"
        analysis_url = "http://skuld.bmsc.washington.edu/~tlsmd/pdb/%s/ANALYSIS" % (pdbid)
        analysis_title = "Analysis of %s" % (pdbid)
        redirect = [self.html_head(title, redirect=analysis_url), 
                    html_title(title),
                    '<center>',
                    '<br><h2>Click below to see the results:</h2>',
                    '<h3><a href="%s">%s</a>' % (analysis_url, analysis_title),
                    '<br><br>',
                    '<font size=-2>You will be redirected automatically in 3 seconds</font>'
                    '</center>'
                    ]
        redirect.append(self.html_foot())
        return "".join(redirect)

def generate_random_filename(code_length = 8):
    """Generates a random 8 character string. Christoph Champ, 2007-12-03
    """
    random.seed()
    codelist = list(5 * string.ascii_letters)
    random.shuffle(codelist)
    code = "".join(random.sample(codelist, code_length))
    return code

def running_stddev(atomnum,restype,resnum,chain,tfactor):
    """Calculates a running standard deviation"""
    ######### EAM 3-Dec-2007 ##########
    tmpfile=generate_random_filename()
    n=nres=atm=nbad=res_tfac=0
    avg_tfac=[]
    std=[]
    res_id=[]
    prevrestype=restype[0]
    prevresnum=resnum[0]
    while n<len(tfactor):
        if( (prevresnum == resnum[n]) and (prevrestype == restype[n]) ):
	   res_tfac=res_tfac+tfactor[n]
	   atm=atm+1
        else:
	   avg_tfac.append(res_tfac/atm) # store previous guy
           res_id.append(resnum[n-1])    # store previous guy
	   res_tfac=tfactor[n]
	   atm=1
           prevrestype=restype[n]
           prevresnum=resnum[n]
        n=n+1
    avg_tfac.append(res_tfac/atm)        # store last guy
    res_id.append(resnum[n-1])           # store last guy

    ## Save B_{mean} per residue for each chain
    fdat=open('%s/%s.dat'%(conf.WEBTMP_PATH,tmpfile),'w')
    r=0
    for b in avg_tfac:
	fdat.write("%s\t%s\n"%(res_id[r],b))
	r=r+1
        if (r < len(res_id)) and (res_id[r] < res_id[r-1]):
           fdat.write("\n")
    fdat.close()

    ## Save RMSD(B) +/-5 residues
    ### FIXME EAM
    ### Not correct, because it crosses chain boundaries
    ### and because the wrong value is calculated (std of mean, 
    ### rather than the std of the atoms)
    fstd=open('%s/%s.std'%(conf.WEBTMP_PATH,tmpfile),'w')
    for s in range(5,len(avg_tfac)-5):
        stddev11 = numpy.std(avg_tfac[s-5:s+5])
	fstd.write("%s\t%s\n"%(res_id[s],stddev11))
	if stddev11 < 0.05:
	   nbad=nbad+1
        if (s < len(res_id)) and (res_id[s+1] < res_id[s]):
           fstd.write("\n")
    fstd.close()

    return nbad, tmpfile

def check_upload(file):
    """Runs sanity checks on uploaded file"""
    # Standard deviation of temperature factors check
    atom_num=[]
    res_type=[]
    res_num=[]
    chain=[]
    temp_factors=[]
    bad_std= -1
    for line in file:
        if line.startswith('EXPDTA    NMR'):
            return "NMR structure! Please do not submit NMR structures, theoretical models, or any PDB file with unrefined Bs."
        elif line.startswith('ATOM'):
	    atomnum=line[7:11]
            restype=line[17:20]
            ch_id=line[22:22]
            resnum=line[23:26]
            tfactor=line[61:66]
            atomnum.strip()
            restype.strip()
            resnum.strip()
            tfactor.strip()
            atomnum=int(atomnum)
            resnum=int(resnum)
            tfactor=float(tfactor)
            atom_num.append(atomnum)
            res_type.append(restype)
            res_num.append(resnum)
            chain.append(ch_id)
            temp_factors.append(tfactor)
        else:
            continue
    if(len(atom_num)<30):
        return "Not a PDB structure"

    bad_std,tmpfile=running_stddev(atom_num,res_type,res_num,chain,temp_factors)
    if bad_std > 0:
	f=open('%s/%s.gnu'%(conf.WEBTMP_PATH,tmpfile),'w')
	f.write("set style fill solid 0.15 noborder\n")
	f.write("set style data linespoints\n")
	f.write("set output '%s/%s.png'\n" %(conf.WEBTMP_PATH,tmpfile))
	f.write("set yrange [0:*]\n")
	## f.write("set ylabel 'Å^2' norotate tc rgb 'blue'\n")
	f.write("set ytics nomirror tc rgb 'blue'\n")
	f.write("set y2range [0:1]\n")
	## f.write("set y2label 'Å^2' norotate tc rgb 'red'\n")
	f.write("set y2tics nomirror tc rgb 'red'\n")
	f.write("set format y2 '%.1f'\n")
	f.write("set xlabel 'residue number'\n")
	f.write("set grid\n")
	f.write("set title 'Distribution of B factors in submitted structure (Å^2)'\n")
	#f.write("set label 1 'bad_std = %d' at graph 0.1, 0.9 noenh\n" % bad_std)   # EAM DEBUG
	#f.write("set term png font '%s' enhanced size 800,400\n" % gnuplot_font)
	f.write("set term png font '%s' enhanced truecolor\n" % conf.GNUPLOT_FONT)
	f.write("plot '%s/%s.std' using 1:($2<0.1 ? 999 : 0) axes x1y2 w filledcurve lt -1 notitle, \\\n"%(conf.WEBTMP_PATH,tmpfile))
	f.write("     '%s/%s.dat' using 1:2 axes x1y1 lt 3 pt 1 title 'B_{mean} per residue', \\\n"%(conf.WEBTMP_PATH,tmpfile))
	f.write("     '%s/%s.std' using 1:2 axes x1y2 lt 1 pt 1 title 'RMSD(B) +/-5 residues', \\\n"%(conf.WEBTMP_PATH,tmpfile))
	f.write("     0.05 axes x1y2 with lines lc rgb 'red' notitle\\\n")
	f.close()
	subprocess.Popen([r"%s"%conf.GNUPLOT,"%s/%s.gnu"%(conf.WEBTMP_PATH,tmpfile)]).wait()
	return "Standard deviation of temperature factors is less than 0.05 for those residues in the shaded regions below:<br/><img src='%s/%s/%s.png'/>" % (conf.BASE_PUBLIC_URL,"webtmp",tmpfile)

    return ''


def main():
    page = None
    form = cgi.FieldStorage()

    if form.has_key("page"):

        if form["page"].value == "explore":
            page = ExploreJobPage(form)

        elif form["page"].value == "admin":
            page = AdminJobPage(form)

        elif form["page"].value == "submit1":
            page = Submit1Page(form)

        elif form["page"].value == "submit2":
            page = Submit2Page(form)

        elif form["page"].value == "submit3":
            page = Submit3Page(form)

        elif form["page"].value == "submit_pdb":
            page = SubmitPDBPage(form)

    if page is None:
        page = QueuePage(form)

    try:
        print page.html_page()
    except xmlrpclib.Fault, fault:
        fault_html = "xmlrpclib.Fault:<br>fault code: %s<br>fault string: %s" % (
            fault.faultCode, fault.faultString.replace("\n","<br>"))

        page = ErrorPage(form, fault_html)
        print page.html_page()

    except socket.error, err:
        page = ErrorPage(form, "socket.error: " + str(err))
        print page.html_page()

    except SubmissionException, err:
        page = ErrorPage(form, str(err))
        print page.html_page()


if __name__=="__main__":
    main()
    sys.exit(0)
