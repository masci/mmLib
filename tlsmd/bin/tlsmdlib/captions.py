## TLS Minimized Domains (TLSMD)
## Copyright 2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

MOTION_ANALYSIS_TEXT = """\
For each protein chain, this analysis show the optimal division into
1 TLS group, 2 TLS groups, 3 TLS groups, etc, up to 20 groups.
It goes on to analyze the implied rigid body translational and rotatational
motion of each group, as well as its quality of fit to the refined atomic
displacement parameters (B-factors). If you have multiple chains in your
structure, they are treated independently.
"""

MULTI_CHAIN_ALIGNMENT_TEXT = """\
When multiple chains are present in your input structure, a side-by-side
sequence alignment is generated to show how the TLS group selection for one
chain aligns with the selection for the other chains.  This analysis is only
meaningful when there are multiple chains of the same sequence in the
asymmetric unit.
"""

REFINEMENT_PREP_TEXT = """\
The CCP4 macromolecular refinement program Refmac5 implements a
refinement mode in which refinement of individual isotropic temperature factors
is supplemented by refinement of one or more TLS groups.
It is usually plausible to try assigning a single TLS group to each protein
molecule, but how one would subdivide the molecules into multiple groups
has been problematic.  This is one important function served by TLSMD analysis.
TLSMD will help you choose how many TLS groups to subdivide your protein chains
into, and it will generate a PDB and TLSIN files for input to further refinement
using Refmac5.
"""

REFINEMENT_PREP_INFO = """\
The optimized TLS groups calculated by TLSMD from a isotropically refined structure may be used to further refine the structure with the TLS + restrained refinement mode of Refmac5.  Given the number of TLS groups you would like to use for each chain, TLSMD will generate a special structure model file (PDBIN) and TLS tensor file(TLSIN) you can use as input files to Refmac5.  These files are generated specifically for Refmac5 refinement by splitting the temperature factor magnitude of each atom between the TLS model and individual atomic temperature factors.  Because of this, the TLS model for refinement is different than the one used for motion analysis.  Select the number of TLS groups to use for each chain and click <b>OK</b> to generate these input files.
"""

LSQR_CAPTION = """\
TLSMD selects the optimal partition of a chain into <N> TLS groups by minimizing
an overall residual function.
This plot shows the value of the residual as a function of the number
of TLS groups allowed in the partition.
Adding additional TLS groups will always make this residual lower, but there
is an issue of diminishing returns as you go to larger numbers of groups.
The details of each optimal <N> group partition are analyzed below.
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
those used for the TLS groups in the various structure visualizations.
"""

LIBRATION_GRAPH_CAPTION = """\
This graph shows the displacement of mainchain atoms implied by the three
screw axes of the TLS group to which they belong.  The screw displacement axes are
calculated in terms of a Gaussian variance-covariance tensor, and displacment
magnitude is shown at a 85% isoprobability magnitude like the translational
displacement.  Protein segments undergoing hinge-like motion show up as peaks in
this graph.
"""

FIT_GRAPH_CAPTION = """\
This graph assesses the quality of the TLS prediction for each residue
by graphing the difference between the refined (input) mainchain atom B factors
and the corresponding B factors implied by the TLS model alone.  If the
TLS model were a perfect description of the observed thermal motion
described by the input structural model, this plot would consist of a
line at 0. <b>Warning:</b> If the input structural model was itself generated
by refinement of a similar TLS model, this will introduce a strong bias towards
low residual differences.
"""
