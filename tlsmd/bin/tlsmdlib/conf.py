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
import const

## BEGIN: CONFIGURATION PATHS AND URLS
TLSMD_ROOT             = os.environ.get("TLSMD_ROOT", "/home/jpaint/tlsmd")
TLSMD_WWW_ROOT         = "/home/jpaint/public_html"
TLSMD_BASE_URL         = "/~jpaint"
WEBTLSMDD              = "http://localhost:10200"
WEBTLSMDD_DATABASE     = "/home/jpaint/database/webtlsmd.db"
ADMIN_PASSWORD_FILE    = "/home/tlsmd/database/cgi-admin"
MSMTP                  = "/usr/bin/msmtp"
## END: CONFIGURATION PATHS AND URLS


## derived paths
TLSMD_PROGRAM_PATH     = os.path.join(TLSMD_ROOT, "bin", "tlsmd.py")
GNUPLOT_FONT           = os.path.join(TLSMD_ROOT, "fonts/LucidaSansOblique.ttf")
REFINEPREP_URL         = "%s/cgi-bin/refineprep.cgi" % (TLSMD_BASE_URL)
TLSMD_WORK_DIR         = os.path.join(TLSMD_WWW_ROOT, "jobs")
TLSMD_WORK_URL         = "%s/jobs" % (TLSMD_BASE_URL)
JMOL_DIR               = "../../../jmol"

## the isoprobability contour level for all visualizations
ADP_PROB = 50

## number of TLS partitons for each chain
NPARTS = 20

## the pixel width of the TLS visualization rendered ray traces
VIS_WIDTH = 800

## the JMol viewer is a square window, generated with this pixel size
JMOL_SIZE = 600

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
        print "TLS Motion Determination (TLSMD) Version %s" % (const.VERSION)
        print
        print "TLS PARAMETER FIT ENGINE...........: %s" % (self.tls_model)
        print "MIN_SUBSEGMENT_SIZE................: %d" % (self.min_subsegment_size)
        print "ATOM B-FACTOR WEIGHT_MODEL.........: %s" % (self.weight_model)
        print "PROTEIN ATOMS CONSIDERED...........: %s" % (self.include_atoms)
        print

    def verify(self):
        assert self.tls_model     in ["ANISO", "ISOT", "NLANISO", "NLISOT"]
        assert self.weight_model  in ["UNIT", "IUISO"]
        assert self.include_atoms in ["ALL", "MAINCHAIN", "CA"]
        
globalconf = GlobalConfiguration()

