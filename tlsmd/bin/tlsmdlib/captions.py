## TLS Minimized Domains (TLSMD)
## Copyright 2002 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

MOTION_ANALYSIS_TEXT = """\
This analysis explicity shows how each protein chain can be split into
multiple TLS groups using one to twenty adjecent, continous groups
along the chain sequence.  It goes on to analyize the implied rigid body
translational and rotatational motion of each group, as well as its
quality of fit to the refined atomic displacement parameters (B-factors).
"""

MULTI_CHAIN_ALIGNMENT_TEXT = """\
When multiple chains are present in your input structure, a side-by-side
sequence alignment is generated to show how the TLS group selection of one
chain aligns with the selection the other chains.  This analysis is only
meaningful when there are multiple chains of the same sequence in the
asymetric unit.
"""

REFINEMENT_PREP_TEXT = """\
The macromolecular refinement program Refmac5 from CCP4 implements a
refinement mode where a TLS description can be added to individual isotropic
temperature factors.  Traditionally, one TLS group is assigned for each
chain because there has been no technique for selecting multiple TLS groups
from the crystallographic data.  However, this is the calculation TLSMD
performs.  By using this TLSMD refinement preperation, you can choose
the number of TLS groups to use per chain, and generate a input PDB
and TLSIN file for refinement using Refmac5.
"""

LSQR_CAPTION = """\
TLSMD selects TLS groups by the minimization of a residual function.
This plot shows the value of the residual as a function of the number
of TLS groups allowed to be used in in the minimization.  Using a given
number of TLS groups, the residual value in the plot above is the lowest
found when all possible choices of protein chain continous TLS groups are
considered.  The details of these TLS groups are analyzed below.
"""

SEG_ALIGN_CAPTION = """\
This plot shows the location along the protein sequence of the
optimal TLS group segments, and how those segments align with
the optimal TLS group segments as the number of TLS groups used
increases.
"""

TRANSLATION_GRAPH_CAPTION = """\
This graph shows the TLS group translational displacement magnitude
of the three principal components of the reduced T tensor at a
isoprobability magnitude of 85%.  The line colors are the same as
those used for the TLS groups in the structure visualization.
"""

LIBRATION_GRAPH_CAPTION = """\
This graph shows the displacement caused by the three TLS group screw axes
on the mainchain atoms of the protein.  The screw displacement axes are
calculated in terms of a Gaussian variance-covariance tensor, and displacment
magnituce is shown at a 85% isoprobability magnitude like the translational
displacement.  Protein segments with hinge-like flexibility show up as peaks in
this graph.
"""

FIT_GRAPH_CAPTION = """\
This graph assesses the quality of the TLS prediction for each TLS group
spanning the residue chain by graphing the difference in the refined (input)
mainchain atom B factors from the TLS model predicted B factors.  If the
TLS model was a perfect fit to the input structure data, this would be
a line at 0.0.
"""
