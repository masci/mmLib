## Copyright 2002-2010 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import math

PI       = math.pi
PI2      = math.pi**2
PI3      = math.pi**3

RAD2DEG  = 180.0 / PI
DEG2RAD  = PI / 180.0
RAD2DEG2 = RAD2DEG**2
DEG2RAD2 = DEG2RAD**2

## converts between U (angstrom^2) temperature factor values and B temperature
## factor values.
U2B = 8.0 * PI2
B2U = 1.0 / (8.0 * PI2)
B2UE4 = B2U * 10000.0

## Path to render (Raster3D) binary
RENDER = "/usr/local/bin/render"

AMINO_BACKBONE   = ["N", "CA", "C"]
NUCLEIC_BACKBONE = ["P", "O5'", "C5'", "C4'", "C3'", "O3'",
                    "P", "O5*", "C5*", "C4*", "C3*", "O3*"]
BACKBONE_ATOMS = AMINO_BACKBONE + NUCLEIC_BACKBONE
