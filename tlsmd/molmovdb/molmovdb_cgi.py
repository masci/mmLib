## TLS Motion Determination (TLSMD)
## Copyright 2002-2006 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import cgitb
import cgi

import sys
import os
import math
import urllib
import cStringIO
import numpy

import kid

from mmLib import Constants, FileIO, Structure, TLS, AtomMath, ConsoleOutput
from tlsmdlib import conf, atom_selection, tlsmdmodule, \
     tlsmd_analysis, table, tls_calcs, console, gnuplots


BASE_URL = "http://dev.molmovdb.org/uploads/"
BASE_DIR = "/usr/local/server/uploads"

NOID_HTML = """\
Content-type: text/html\n\n
<HTML>
  <HEAD><TITLE>TLSMD CGI Script Error</TITLE></HEAD>
  <BODY>
    <CENTER>
      <H1>TLSMD CGI Script Error</H1>
      <P>[ERROR] http://<scripturl>?ID=<idstring> Argument required</P>
    </CENTER>
  </BODY>
</BODY>
</HTML>
"""

REDIRECT_HTML_TEMPLATE = """\
Content-type: text/html\n\n
<HTML>
  <HEAD><TITLE>Redirect</TITLE></HEAD>
  <META http-equiv="refresh" content="0; Url=<TLSMD>">
  <BODY>
    <CENTER>
      <P>You should have been redirected to <A HREF="<TLSMD>"><TLSMD></A></P>
    </CENTER>
  </BODY>
</BODY>
</HTML>
"""


class MorphLocations(object):
    def __init__(self, morph_id):
        self.morph_id = morph_id
        self.base_url = os.path.join(BASE_URL, morph_id)
        self.base_dir = os.path.join(BASE_DIR, morph_id)

    def source_url(self):
        return os.path.join(self.base_url, "ff0.pdb")

    def directory(self):
        return self.base_dir

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
    
    def stylesheet_url(self):
        return "http://www.drizzle.com/~jpaint/molmovdb.css"


def HingePredictionTable(chain):
    """Adds a table.StringTable instance to chain.tbl contining one row per
    residue from chain.  The table's first column contains the residue
    sequence number + insertion code, and the remaining columns from left
    to right contain the TLSMD inspired hinge residues from the
    2-TLS, 3-TLS, ..., N-TLS group model.
    """
    cols = 3 + conf.globalconf.nparts - 1
    tbl = table.StringTable(
        len(chain), cols, "0",
        title = "TLSMD Predicted Hinge Residues for Chain %s" % (chain.chain_id),
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


def CGI_Main(form):
    """main() when run as a CGI script.
    """
    ## shut off mmLib and tlsmd console output
    ConsoleOutput.disable()
    console.disable()
    
    ## enable CGI exception catching
    cgitb.enable()

    ## get the morph ID
    morph_id = form["ID"].value
    ## don't let an attacker be tricky
    assert morph_id.find("/") == -1 and morph_id.find("..") == -1 and len(morph_id) < 20
    morph = MorphLocations(morph_id)

    ## get the PDB file
    pdb_buff = urllib.urlopen(morph.source_url()).read()

    ## run TLSMD
    fobj = cStringIO.StringIO(pdb_buff)
    tlsmd = TLSMDHingePredictor(fobj)

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
    template.stylesheet_url = morph.stylesheet_url()
    
    open(html_path, "w").writelines(template.generate(output = "html"))
    
    ## redirect browser
    sys.stdout.write(REDIRECT_HTML_TEMPLATE.replace("<TLSMD>", morph.html_url()))


def main():
    """Choose between command-line main and CGI main.
    """
    form = cgi.FieldStorage()
    if form.has_key("ID"):
        CGI_Main(form)
    else:
        sys.stdout.write(NOID_HTML)

    
if __name__ == "__main__":
    main()
