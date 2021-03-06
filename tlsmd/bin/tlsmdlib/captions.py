## TLS Minimized Domains (TLSMD)
## Copyright 2002-2010 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

OPTIMIZATION_PARAMS_TEXT = """\
The optimization algorithm at the core of TLSMD has several options which are automatically chosen for you by analyzing your input structure during submission.  Most of these options require detailed explanation which we will present later in publication.  The <i>Included Atoms</i> option is, however, pretty straight forward.  It controls whether to include all atoms, or just the main-chain atoms in the optimization.  Since TLSMD can take quite a long time to run, structures containing chains longer than 400 residues are automatically set to use only the main-chain protein atoms <b>N</b>, <b>CA</b>, and <b>C</b>.  If no chain in the structure exceeds 400 residues in length, the all protein atoms are used in the optimization.
"""

MOTION_ANALYSIS_TEXT = """\
For each protein chain, this analysis show the optimal division into
1 TLS group, 2 TLS groups, 3 TLS groups, etc, up to 20 groups.
It goes on to analyze the implied rigid body translational and rotational
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
The CCP4 refinement program Refmac5
and the PHENIX program PHENIX.refine
both implement a refinement mode in which individual isotropic temperature factors
are supplemented by refinement of one or more TLS groups.
TLSMD will help you choose how many TLS groups to subdivide your protein chains
into, and it will generate a input files for you to use in further refinement
using Refmac5 or PHENIX.
"""

REFINEMENT_PREP_INFO = """\
You may find the graph below useful in selecting the number of groups.  Often there is a 'dogleg' or 'elbow' point in the curve to the right of which adding more TLS groups does not significantly reduce the least squares residual.  The number of TLS groups at that point is a reasonable number to select.
<p><b>Refmac5</b>:
Given the number of TLS groups you would like to use for each chain, TLSMD will generate a special coordinate file (PDBIN) and TLS tensor file (TLSIN) that you can use for input to Refmac5.  These files are generated specifically for Refmac5 refinement by splitting the temperature factor magnitude of each atom between the TLS model and individual atomic temperature factors.  Because of this, the starting B factors for further refinement are different from the ones used for motion analysis.  
<p><b>PHENIX</b>:
Given the number of TLS groups you would like to use for each chain, TLSMD will generate a single file for input to PHENIX that describes the residues making up each group.  PHENIX does not currently allow you to input a starting model for the TLS parameters themselves.
"""

REFINEMENT_FILES_DOWNLOAD_INFO = """\
<p><b>Refmac5:</b> Download both the modified PDBIN file for your structure and
the corresponding TLSIN file. Feed these to REFMAC5 as a starting point for
multi-TLS group refinement. See the TLSMD documentation for detailed instructions.</p>
<p><b>PHENIX:</b> The PHENIX file contains a description of the TLS groups you
selected. This file is intended to be read by the PHENIX.refine input scripts.</p>
"""

LSQR_CAPTION = """\
TLSMD selects the optimal partition of a chain into 1 to 20 TLS groups by minimizing
an overall residual function.
This plot shows the value of the residual as a function of the number
of TLS groups allowed in the partition.
Adding additional TLS groups will always make this residual lower, but there
is an issue of diminishing returns as you go to larger numbers of groups.
The individual TLS segments are analyzed below.
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
This graph shows the displacement of main chain atoms implied by the three
screw axes of the TLS group to which they belong.  The screw displacement axes are
calculated in terms of a Gaussian variance-covariance tensor, and displacement
magnitude is shown at a 85% isoprobability magnitude like the translational
displacement.  Protein segments undergoing hinge-like motion show up as peaks in
this graph.
"""

FIT_GRAPH_CAPTION = """\
This graph assesses the quality of the TLS prediction for each residue
by graphing the difference between the refined (input) main chain atom B factors
and the corresponding B factors implied by the TLS model alone.  If the
TLS model were a perfect description of the observed thermal motion
described by the input structural model, this plot would consist of a
line at 0.
"""

RMSD_DEVIATION_GRAPH_CAPTION = """\
This plot is similar to the plot above, except that instead of plotting
(Bpred - Bobs) for the CA atom only, the colored line plots the
RMSD of |Bpred-Bobs| for all atoms in the residue. If the TLSMD model
were a perfect description of the input B values, then the colored line
would be a horizontal line at 0.
"""

NO_VALID_CONFIGURATIONS = """\
TLSMD was unable to optimize this chain.  This is due to known limitation
in this beta version of TLSMD which will be fixed in a upcoming version.
"""

TLS_GROUP_RECOMBINATION = """\
The TLSMD optimization algorithm models TLS groups as sequential segments
of a protein or DNA/RNA chain.  This matrix shows the RMSD B values of the
individual groups on the diagonal, and the RMSD B values of combined groups as
off-diagonal elements.  This helps identify non-contiguous protein
segments which may be combined into a single TLS group.
"""

## TODO: Add captions for animated GIFs, 2009-07-08
