## TLS Motion Determination (TLSMD)
## Copyright 2002-2008 by TLSMD Development Group (see AUTHORS file)
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
import console

## BEGIN: CONFIGURATION PATHS AND URLS
BASE_PUBLIC_URL        = "http://verdandi.bmsc.washington.edu"  ## Added by Christoph Champ, 2008-02-07
TLSMD_ROOT             = os.environ.get("TLSMD_ROOT", "/home/tlsmd/tlsmd")
TLSMD_WWW_ROOT         = "/home/tlsmd/public_html"
TLSMD_PUBLIC_URL       = "http://verdandi.bmsc.washington.edu/~tlsmd" # Added by Christoph Champ, 2007-12-13
TLSMD_BASE_URL         = "/~tlsmd"
WEBTLSMDD              = "http://localhost:10100"
WEBTLSMDD_DATABASE     = "/home/tlsmd/database/webtlsmd.db"
ADMIN_PASSWORD_FILE    = "/home/tlsmd/database/cgi-admin"
MAIL                   = "/bin/mail"
TRACEBACK_EMAIL        = "tlsmdtraceback"
LOG_PATH               = "/home/tlsmd/log/tlsmd_runlog.txt"
PDB_URL                = "http://www.pdb.org/pdb/explore/explore.do?structureId=" ## Added by Christoph Champ, 2008-02-20
GET_PDB_URL            = "http://www.rcsb.org/pdb/files" ## Added by Christoph Champ, 2008-03-10
## END: CONFIGURATION PATHS AND URLS

## override default configuration
ALTCONF = "/etc/tlsmd.conf"
if os.path.exists(ALTCONF):
    execfile(ALTCONF)

## derived paths
TLSMD_PROGRAM_PATH     = os.path.join(TLSMD_ROOT, "bin", "tlsmd.py")
WEBTMP_PATH	       = "/var/www/webtmp"		# Added by Christoph Champ, 2007-12-03
GNUPLOT                = "/usr/local/bin/gnuplot"	# Added by Christoph Champ, 2007-12-03
GNUPLOT_FONT           = os.path.join(TLSMD_ROOT, "fonts/LucidaSansOblique.ttf")
#REFINEPREP_URL         = "%s/cgi-bin/refineprep.cgi" % (TLSMD_BASE_URL)
REFINEPREP_URL         = "%s/cgi-bin/refineprep.cgi" % (TLSMD_PUBLIC_URL) # Changed. Christoph Champ, 2008-03-17
TLSMD_WORK_DIR         = os.path.join(TLSMD_WWW_ROOT, "jobs")
TLSMD_WORK_URL         = "%s/jobs" % (TLSMD_BASE_URL)
JMOL_DIR               = "../../../jmol"		# Directory path must be relative, not an absolute URL. Christoph Champ, 2008-03-17
#JMOL_DIR               = "%s/jmol" % (TLSMD_PUBLIC_URL) # Added by Christoph Champ, 2007-12-13
WEBTLSMDD_PDB_DIR      = os.path.join(TLSMD_WWW_ROOT,"pdb")
WEBTLSMDD_PDBID_FILE   = os.path.join(WEBTLSMDD_PDB_DIR,"pdbids.txt")

## the isoprobability contour level for all visualizations
ADP_PROB = 50

## number of TLS partitons for each chain
NPARTS = 20

## the pixel width of the TLS visualization rendered ray traces
VIS_WIDTH  = 640
VIS_HEIGHT = 400

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
        self.recombination = False
        self.adp_smoothing = 0

    def prnt(self):
        console.stdoutln("TLS Motion Determination (TLSMD) Version %s" % (const.VERSION))
        console.endln()
        console.kvformat("TLS PARAMETER FIT ENGINE", self.tls_model)
        console.kvformat("MIN_SUBSEGMENT_SIZE", self.min_subsegment_size)
        console.kvformat("ATOM B-FACTOR WEIGHT_MODEL", self.weight_model)
        console.kvformat("PROTEIN ATOMS CONSIDERED", self.include_atoms)
        console.endln()

    def verify(self):
        assert self.tls_model     in ["ANISO", "ISOT", "NLANISO", "NLISOT"]
        assert self.weight_model  in ["UNIT", "IUISO"]
        assert self.include_atoms in ["ALL", "MAINCHAIN", "CA"]
        
globalconf = GlobalConfiguration()

