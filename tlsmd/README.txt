TLS Minimized Domains 0.2.0
------------------------------------------------------------------------------

TLS Minimal Domains (TLSMD) is a program which computes optimal 
TLS (rigid-body) descriptions from the isotropically or anisotropically 
refined ADPs of a protein crystal structure.   This is accomplished in a
somewhat brute-force mannor.  First, a single TLS model is fit by least
squares to every possible subsegment of a protein chain and its its 
goodness of fit assessed by the least squares fit residual.  The number of
possible subsegments of a protein depends on its numer of residues, n, and
the minimum segment width, w, is given by the equation below:

S(n,w) = (n*n + m*m - 2*n*m - n + m) / 2

Once the quality of fit to the TLS model has been assessed, a graph is
constructed by placing a vertex before the first residue in the chain,
after the last residue in the chain, and one vertex between each residue
in the protein chain.  This results in a total of n+1 vertexes.  Edges
are constructed using the TLS group subsegments spanning the residues 
the TLS groups were fit to, which a edge cost equal the the least squares
residual of the TLS model fit to the input thermal parameters (ADPs).
Constrained Bellman-Ford is then used to determine the least cost path
from the source vertex at the beginning of the chain, to the destination
vertex at the end of the protein chain.

------------------------------------------------------------------------------
Jay Painter <jpaint@u.washington.edu> <jpaint@gmail.com>
Aug 10, 2005
