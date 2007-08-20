## TLS Motion Determination (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

"""Installation settings.  When installed to run as a web application, the layout
of the directories from the tlsmd root directory is as follows:
"""

import os
import time
import console

from django.conf import settings

class GlobalConfiguration(object):
    def __init__(self):
        self.tls_model = "ISOT"
        self.weight_model = "UNIT"
        self.include_atoms = "ALL"
        self.min_subsegment_size = 4
        self.adp_prob = settings.ADP_PROB
        self.nparts = settings.NPARTS
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
        console.stdoutln("TLS Motion Determination (TLSMD) Version %s" % (settings.VERSION))
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

