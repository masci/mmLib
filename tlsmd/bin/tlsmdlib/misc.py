## TLS Minimized Domains (TLSMD)
## Copyright 2002 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

###############################################################################
## Imports
##

## standard Python library
import os
import sys
import math
import string
import time
import cPickle
import thread
import xmlrpclib
import SimpleXMLRPCServer

from threading          import *
from Queue              import *

## Numeric Python imports
from Numeric       import *
from LinearAlgebra import *

## mmLib imports
from mmLib.FileLoader     import *
from mmLib.Structure      import *
from mmLib.Extensions.TLS import *

###############################################################################
## Constants
##

## unreasonably small T/L eigenvalues
## B<=0.01
TSMALL = 0.1 * B2U

## L RMSD <= 0.1 DEG
LSMALL = (0.1)**2 * DEG2RAD2

## threashold for small/negitive Eigenvalues for T,L,U tensors
SMALL_EIGENVALUE = 1.0e-10

###############################################################################
## Globals
##
GLOBALS = {
    "TLS_MODEL":           "HYBRID",
    "WEIGHT_MODEL":        "UNIT",
    "INCLUDE_ATOMS":       "ALL",
    "MIN_SUBSEGMENT_SIZE": 6,

    "VERSION":        "0.0.1",
    "RELEASE_DATE":   "June 21, 2005",
    "AUTHOR":         "Jay Painter",
    "EMAIL":          "jpaint@u.washington.edu",

    "WEBTLSMDD":      None,
    "JOB_ID":         None,
    "REFINEPREP_URL": "/~jpaint/cgi-bin/refineprep.cgi"
    }

###############################################################################
## Utility Funcs
##
_STIME = 0.0

def start_timing():
    global _STIME
    _STIME = time.time()

def end_timing():
    global _STIME
    tm = time.time() - _STIME
    return "Computation Time: %5.2f sec" % (tm)

def begin_chain_timing(chain_id):
    print "BEGIN TIMING CHAIN %s %f" % (chain_id, time.time())

def end_chain_timing(chain_id):
    print "END TIMING CHAIN %s %f" % (chain_id, time.time())
