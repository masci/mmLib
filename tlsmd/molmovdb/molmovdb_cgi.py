## TLS Motion Determination (TLSMD)
## Copyright 2002-2009 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import sys
import os
import time
import traceback
import cgitb
import cgi

import numpy
import kid

from mmLib import ConsoleOutput
from tlsmdlib import conf, tlsmd_analysis, table, console

## <CONFIGURE>
BASE_URL = "http://localhost/~tlsmd/molmovdb/"
BASE_DIR = "/home/tlsmd/public_html/molmovdb/"
STYLESHEET_URL = "molmovdb.css"
## </CONFIGURE>


NOID_HTML_TEMPLATE = """\
<HTML>
  <HEAD>
    <TITLE>ERROR: No Morph ID</TITLE>
    <LINK rel="stylesheet" href="${stylesheet_url}" type="text/css" media="screen">
  </HEAD>
  <BODY><DIV id="page">
    <H1>ERROR: No Morph ID</H1>
    <P>The TLSMD CGI script must be invoked with a valid morph ID<br>
      http://<scripturl>?ID=<idstring> Argument required
    </P>
  </DIV></BODY>
</BODY>
</HTML>
"""

INVALIDID_HTML_TEMPLATE = """\
<HTML>
  <HEAD>
    <TITLE>ERROR: Invalid Morph ID ${morph_id}</TITLE>
    <LINK rel="stylesheet" href="${stylesheet_url}" type="text/css" media="screen">
  </HEAD>
  <BODY><DIV id="page">
    <CENTER>
      <H1>ERROR: Invalid Morph ID ${morph_id}</H1>
      <P>The source structure file <i>${struct_path}</i> cannot be found.</P>
    </CENTER>
  </DIV></BODY>
</BODY>
</HTML>
"""

REDIRECT_HTML_TEMPLATE = """\
<HTML>
  <HEAD>
    <TITLE>Redirect</TITLE>
    <LINK rel="stylesheet" href="${stylesheet_url}" type="text/css" media="screen">
  </HEAD>
  <META http-equiv="refresh" content="0; Url=${tlsmd}">
  <BODY><DIV id="page">
    <CENTER>
      <P>You should have been redirected to <A HREF="${tlsmd}">${tlsmd}</A></P>
    </CENTER>
  </DIV></BODY>
</BODY>
</HTML>
"""

RELOAD_UNTIL_COMPLETE_HTML_TEMPLATE = """\
<HTML>
  <HEAD>
    <TITLE>TLSMD Processing</TITLE>
    <LINK rel="stylesheet" href="${stylesheet_url}" type="text/css" media="screen">
  </HEAD>
  <META http-equiv="refresh" content="20; Url=${tlsmd}">
  <BODY><DIV id="page">
    <CENTER>
      <P>This page will reload itself every 20 seconds until the TLSMD job is complete.<BR>
      <A HREF="${tlsmd}"><${tlsmd}></A>
      </P>
    </DIV></CENTER>
  </BODY>
</BODY>
</HTML>
"""

SUBPROCESS_TRACEBACK_HTML_TEMPLATE = """\
<HTML>
  <HEAD>
    <TITLE>TLSMD Processing</TITLE>
    <LINK rel="stylesheet" href="${stylesheet_url}" type="text/css" media="screen">
  </HEAD>
  <BODY><DIV id="page">
    <H1>TLSMD Error on Morph ID ${morph_id}</H1>
    <PRE id="traceback">${tracebk}</PRE>
  </DIV></BODY>
</BODY>
</HTML>
"""


def SimpleTemplateReplace(template, **args):
    """Simple key/value replacement of template string with
    the dictionary values of args.
    """
    for key, value in args.iteritems():
        keytag = "${%s}" % (key)
        template = template.replace(keytag, value)
    return template


def CGIResponse(html):
    """Write a CGI repsonse with the argument html
    to standard output.
    """
    sys.stdout.write("Content-type: text/html\n\n")
    sys.stdout.write(html)
    

class MorphLocations(object):
    def __init__(self, morph_id):
        self.morph_id = morph_id
        self.base_url = os.path.join(BASE_URL, morph_id)
        self.base_dir = os.path.join(BASE_DIR, morph_id)

    def directory(self):
        return self.base_dir

    def source_path(self):
        return "ff0.pdb", os.path.join(self.directory(), "ff0.pdb")

    def source_url(self):
        return os.path.join(self.base_url, "ff0.pdb")

    def chain_datafile_index_path(self):
        filename = "tlsmd_chain_index.txt"
        path = os.path.join(self.directory(), filename)
        return filename, path

    def chain_datafile_path(self, chain_id):
        filename = "tlsmd_chain_%s.txt" % (chain_id)
        path = os.path.join(self.directory(), filename)
        return filename, path

    def html_path(self):
        filename = "tlsmd.html"
        path = os.path.join(self.directory(), filename)
        return filename, path

    def html_url(self):
        return os.path.join(self.base_url, "tlsmd.html")
    

def HingePredictionTable(chain):
    """Adds a table.StringTable instance to chain.tbl contining one row per
    residue from chain. The table's first column contains the residue
    sequence number + insertion code, and the remaining columns from left
    to right contain the TLSMD inspired hinge residues from the
    2-TLS, 3-TLS, ..., N-TLS group model.
    """
    cols = 3 + conf.globalconf.nparts - 1
    tbl = table.StringTable(
        len(chain), cols, "0",
        title = "TLSMD Predicted Hinge Residues for Chain %s" % (
            chain.chain_id),
        column_titles = ["SeqNum", "ResName", "MolMovDB Hinge Prediction"] + ["%d TLS Groups" % (ntls) for ntls in xrange(2,conf.globalconf.nparts + 1)])

    chain.tbl = tbl

    ## initalize table
    first_last_frag_id = (chain[0].fragment_id, chain[-1].fragment_id)
    tbl_row_index = {}
    for i, frag in enumerate(chain.iter_fragments()):
        tbl[i, 0] = frag.fragment_id
        tbl[i, 1] = frag.res_name
        tbl[i, 2] = 1
        tbl_row_index[frag.fragment_id] = i

    partition_collection = chain.partition_collection

    ## set table values from TLSMD results
    col = 3
    for cpartition in partition_collection.iter_chain_partitions():
        if cpartition.num_tls_segments() == 1:
            continue
        for tls in cpartition.iter_tls_segments():
            for frag_id1, frag_id2 in tls.iter_segment_ranges():
                for frag_id in (frag_id1, frag_id2):
                    if frag_id not in first_last_frag_id:
                        row = tbl_row_index[frag_id]
                        tbl[row, col] = 1
                        tbl[row, 2] = 0
        col += 1


def TLSMDHingePredictor(struct_fobj):
    """Runs the TLSMD algorithm on the input structure (passed to
    this function as the file object of the PDB file).
    """
    tlsmd = tlsmd_analysis.TLSMDAnalysis(struct_file_object = struct_fobj)
    conf.globalconf.nparts = 5
    tlsmd_analysis.IndependentTLSSegmentOptimization(tlsmd)
    
    for chain in tlsmd.iter_chains():
        HingePredictionTable(chain)

    return tlsmd

def SubProcess(morph):
    """The CGI Script forks off this child function as a child process to
    perform the TLSMD calculations and write out results files.
    """
    ## close all file descriptors
    for fd in xrange(0, 25):
        try:
            os.close(fd)
        except OSError:
            pass

    ## run TLSMD
    struct_basename, struct_path = morph.source_path()
    tlsmd = TLSMDHingePredictor(open(struct_path, "r"))

    ## output tabulated data files
    idx_filename, idx_path = morph.chain_datafile_index_path()
    idx_fobj = open(idx_path, "w")
    for chain in tlsmd.iter_chains():
        filename, path = morph.chain_datafile_path(chain.chain_id)
        idx_fobj.write(filename + "\n")
        open(path, "w").write(str(chain.tbl))
    idx_fobj.close()

    ## write out HTML file
    html_filename, html_path = morph.html_path()
    template = kid.Template(file = "molmovdb.kid")
    template.tlsmd = tlsmd
    template.stylesheet_url = STYLESHEET_URL
    
    open(html_path, "w").writelines(template.generate(output = "html"))


def LaunchSubProcess(morph):
    """Catch any tracebacks in the TLSMD computational code and write
    out the traceback in HTML.
    """
    try:
        SubProcess(morph)
    except:
        html_filename, html_path = morph.html_path()
        open(html_path, "w").write(
            SimpleTemplateReplace(SUBPROCESS_TRACEBACK_HTML_TEMPLATE,
                                  tracebk = traceback.format_exc(),
                                  morph_id = morph.morph_id,
                                  stylesheet_url = STYLESHEET_URL))

    ## child processes of a fork() should always exit like this
    sys.exit(n). _exit()

def CGIProcess(morph):
    ## write temporary reload template
    html_filename, html_path = morph.html_path()
    open(html_path, "w").write(
        SimpleTemplateReplace(RELOAD_UNTIL_COMPLETE_HTML_TEMPLATE,
                              tlsmd = morph.html_url(),
                              stylesheet_url = STYLESHEET_URL))
    
    ## redirect browser
    CGIResponse(
        SimpleTemplateReplace(REDIRECT_HTML_TEMPLATE,
                              tlsmd = morph.html_url(),
                              stylesheet_url = STYLESHEET_URL))


def CGI_Main(morph_id):
    """main() when run as a CGI script.
    """
    ## shut off mmLib and tlsmd console output
    ConsoleOutput.disable()
    console.disable()
    
    ## enable CGI exception catching
    cgitb.enable()

    ## don't let an attacker be tricky
    assert morph_id.find("/") == -1 and morph_id.find("..") == -1 and len(morph_id) < 20
    morph = MorphLocations(morph_id)

    struct_filename, struct_path = morph.source_path()
    if os.path.exists(struct_path):
        pid = os.fork()
        if pid == 0:
            LaunchSubProcess(morph)
        else:
            CGIProcess(morph)
    else:
        CGIResponse(
            SimpleTemplateReplace(INVALIDID_HTML_TEMPLATE,
                                  morph_id = morph_id,
                                  struct_path = struct_path,
                                  stylesheet_url = STYLESHEET_URL))


def main():
    """Choose between command-line main and CGI main.
    """
    form = cgi.FieldStorage()
    morph_id = None

    if form.has_key("ID"):
        morph_id = form["ID"].value
    elif not os.environ.has_key("SERVER_SOFTWARE"):
        try:
            morph_id = sys.argv[1]
        except IndexError:
            sys.stderr.write("usage: tlsmd.cgi ID\n\n")
            raise SystemExit

    if morph_id is not None:
        CGI_Main(morph_id)
    else:
        CGIResponse(
            SimpleTemplateReplace(NOID_HTML_TEMPLATE,
                                  stylesheet_url = STYLESHEET_URL))
            
    sys.exit(0)


if __name__ == "__main__":
    main()
