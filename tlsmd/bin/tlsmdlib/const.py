## TLS Motion Determination (TLSMD)
## Copyright 2002-2010 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""TLSMD program constants
"""
## program info
## NOTE: 3rd decimal place is for internal development (on Verdandi)
LINK_SPACE    = '&nbsp;&nbsp;&nbsp;&nbsp;'
VERSION       = "1.4.0"
RELEASE_DATE  = "11 February 2010"
AUTHOR        = "Ethan Merritt"
EMAIL         = "merritt@u.washington.edu"

## pymmlib
from mmLib import Constants

## unreasonably small T/L eigenvalues B <= 0.01
TSMALL = 0.1 * Constants.B2U

## L RMSD <= 0.1 DEG
LSMALL = (0.1)**2 * Constants.DEG2RAD2

## mainchain atom definitions
AMINO_ACID_MAINCHAIN_ATOMS   = ["N", "CA", "C", "O", "CB"]
NUCLEIC_ACID_MAINCHAIN_ATOMS = ["P", "O5*", "C5*", "C4*", "C3*", "O3*",
                                     "O5'", "C5'", "C4'", "C3'", "O3'"]
MAINCHAIN_ATOMS = AMINO_ACID_MAINCHAIN_ATOMS + NUCLEIC_ACID_MAINCHAIN_ATOMS
MAINCHAIN_ATOM_DICT = {}
for res_name in MAINCHAIN_ATOMS:
    MAINCHAIN_ATOM_DICT[res_name] = True
import re ## used when needing a regular expressions match
RE_MAINCHAIN_ATOMS = re.compile("^ATOM.*  (N|CA|C|O|CB|P|O5\*|O5'|C5\*|C5'|C4\*|C4'|C3\*|C3'|O3\*|O3') (.){20}.*$", re.I)

## used in structcmp.TLSConformationPredctionHypothosis to calc_superposition
SUPER_ATOMS = ["N","CA","C"]

## used with Skittles
JUNCTION_ATOMS = ["N", "C", "O3'", "O3*", "P"]

## used in tls_animate.py to filter selected atoms
DISPLACE_ATOM_NAME_DICT = {
    "CA": True, "P": True, 
    "O5*": True, "C5*": True, "C4*": True, "C3*": True, "O3*": True,
    "O5'": True, "C5'": True, "C4'": True, "C3'": True, "O3'": True
    }

## used for index position
ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
