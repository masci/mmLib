## TLS Motion Determination (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Program constants
"""
from mmLib.mmTypes import *

## unreasonably small T/L eigenvalues B<=0.01
TSMALL = 0.1 * B2U

## L RMSD <= 0.1 DEG
LSMALL = (0.1)**2 * DEG2RAD2

## program info
VERSION = "0.5.7"
RELEASE_DATE = "CVS Version"
AUTHOR = "Jay Painter"
EMAIL = "jpaint@u.washington.edu"

## mainchain atom definitions
MAINCHAIN_ATOMS = ["N","CA","C","O","CB"]
