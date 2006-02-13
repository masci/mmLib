## TLS Motion Determination (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Installation settings.  When installed to run as a web application, the layout
of the directories from the tlsmd root directory is as follows:

$TLSMD_ROOT/cgi-bin 
$TLSMD_ROOT/jobs
$TLSMD_ROOT/jmol
$TLSMD_ROOT/examples

"""
import os
import time

## BEGIN: CONFIGURATION PATHS AND URLS
TLSMD_ROOT             = os.environ.get("TLSMD_ROOT", "/home/jpaint/tlsmd")
TLSMD_WWW_ROOT         = "/home/jpaint/public_html"
TLSMD_BASE_URL         = "/~jpaint"
WEBTLSMDD              = "http://localhost:10200"
ADMIN_PASSWORD_FILE    = "/home/tlsmd/database/cgi-admin"
## END: CONFIGURATION PATHS AND URLS



## derived paths
GNUPLOT_FONT           = os.path.join(TLSMD_ROOT, "fonts/LucidaSansOblique.ttf")
REFINEPREP_URL         = os.path.join(TLSMD_BASE_URL, "cgi-bin/refineprep.cgi")
TLSMD_WORK_DIR         = os.path.join(TLSMD_WWW_ROOT, "jobs")
TLSMD_WORK_URL         = "%s/jobs" % (TLSMD_BASE_URL)

## the isoprobability contour level for all visualizations
ADP_PROB = 85

## number of TLS partitons for each chain
NPARTS = 20

class GlobalConfiguration(object):
    def __init__(self):
        self.tls_model = "ISOT"
        self.weight_model = "UNIT"
        self.include_atoms = "ALL"
        self.min_subsegment_size = 4
        self.adp_prob = ADP_PROB
        self.nparts = NPARTS
        self.verbose = False
        self.use_svg = False
        self.webtlsmdd = None
        self.job_id = None
        self.struct_id = None
        self.start_time = time.time()
        self.target_struct_path = None
        self.target_struct_chain_id = None

    def prnt(self):
        print "TLSMD GLOBAL CONFIGURATION"
        print "    TLS PARAMETER FIT ENGINE.....: %s" % (self.tls_model)
        print "    MIN_SUBSEGMENT_SIZE..........: %d" % (self.min_subsegment_size)
        print "    ATOM B-FACTOR WEIGHT_MODEL...: %s" % (self.weight_model)
        print "    PROTEIN ATOMS CONSIDERED.....: %s" % (self.include_atoms)
        print

    def verify(self):
        assert self.tls_model     in ["ANISO", "ISOT", "NLANISO", "NLISOT"]
        assert self.weight_model  in ["UNIT", "IUISO"]
        assert self.include_atoms in ["ALL", "MAINCHAIN", "CA"]
        
globalconf = GlobalConfiguration()

