## TLS Minimized Domains (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import numpy

class HCSSSP(object):
    """Hop Constrained Single Source Shortest Path graph(V,E) minimization
    based on the Bellman-Ford Algorithm but modified to work with a
    2-dimensional cost(D) matrix, path(P) matrix, and travel(T) matrix.
    """
    def HCSSSP_minimize(self, V, E, hops):
        """Hop-Constrained Single Source Shorted Path minimization,
        loosely based on the Bellman-Ford SSSP algorithm using
        Dynamic Programming.  Returns the D, P, and T matrixes.
        """
        assert len(V)>0
        assert len(E)>0

        num_vertex = len(V)

        ## initialize D/P
        infinity = 1e10

        ## a 2D cost matrix; the value at Dij describes the minimum
        ## cost to reach vertex j by traversing i edges
        D = numpy.zeros((hops+1, num_vertex), float) + infinity

        ## like BellmanFord, initialize the source vertex distance to 0.0
        for i in xrange(hops+1):
            D[i,0] = 0.0

        ## a 2D previous vertex matrix; the value at Pij is the
        ## previous vertex of the path used to achieve cost Dij,
        ## except the previous vertex it describes is not the one
        ## in row i, but the one in row i-1 (the previous row)
        P = numpy.zeros((hops+1, num_vertex), int) - 1

        ## a 2D "travel" matrix containing the edge used by the path
        ## through the previous matrix -- this is a Python 2D matrix and
        ## not a Numerical Python 2d array
        T = []
        for i in xrange(hops+1):
            T.append([None for j in xrange(num_vertex)])

        ## now run the minimization
        for h in xrange(1, hops+1):
            for edge in E:
                self.HCSSSP_minimize_relax(D, P, T, edge, h)

        ## now the matrix Dij and Pij are complete
        return D, P, T
            
    def HCSSSP_minimize_relax(self, D, P, T, edge, hop_constraint):
        """Relax vertices for the current number of hops using the cost array
        from the costs calculated using the previous number of hops.

        Current D for the given number of hops h is D[h], the D
        array for the previous number of hops is D[h-1]
        """
        vertex_i = edge[0]
        vertex_j = edge[1]
        weight   = edge[2]

        ## get the cost vector for the current hop constraint (which we are
        ## in the process of calculating), and the cost vector for
        ## the previous hop constraint (which we assume has been calculated
        ## previously)
        Dp = D[hop_constraint - 1]
        Dc = D[hop_constraint]

        ## perform relaxation for the current number of hops aginst the
        ## cost vector for the previous number of hops; this results
        ## in the current cost vector being the minimum cost using at most
        ## one more hop(edge)
        if Dc[vertex_j] > (Dp[vertex_i] + weight):
            Dc[vertex_j]               = Dp[vertex_i] + weight
            P[hop_constraint,vertex_j] = vertex_i
            T[hop_constraint][vertex_j]= edge
            
    def HCSSSP_maximize(self, V, E, hops):
        """Hop-Constrained Single Source Shorted Path minimization,
        loosely based on the Bellman-Ford SSSP algorithm using
        Dynamic Programming.  Returns the D, P, and T matrixes.
        """
        assert len(V)>0
        assert len(E)>0

        num_vertex = len(V)

        ## a 2D cost matrix; the value at Dij describes the minimum
        ## cost to reach vertex j by traversing i edges
        D = numpy.zeros((hops+1, num_vertex), float)

        ## like BellmanFord, initialize the source vertex distance to 0.0
        for i in xrange(hops+1):
            D[i,0] = 0.0

        ## a 2D previous vertex matrix; the value at Pij is the
        ## previous vertex of the path used to achieve cost Dij,
        ## except the previous vertex it describes is not the one
        ## in row i, but the one in row i-1 (the previous row)
        P = numpy.zeros((hops+1, num_vertex), int) - 1

        ## a 2D "travel" matrix containing the edge used by the path
        ## through the previous matrix -- this is a Python 2D matrix and
        ## not a Numerical Python 2d array
        T = []
        for i in xrange(hops+1):
            T.append([None for j in xrange(num_vertex)])

        ## now run the minimization
        for h in xrange(1, hops+1):
            for edge in E:
                self.HCSSSP_maximize_relax(D, P, T, edge, h)

        ## now the matrix Dij and Pij are complete
        return D, P, T
            
    def HCSSSP_maximize_relax(self, D, P, T, edge, hop_constraint):
        """Relax vertices for the current number of hops using the cost array
        from the costs calculated using the previous number of hops.

        Current D for the given number of hops h is D[h], the D
        array for the previous number of hops is D[h-1]
        """
        vertex_i = edge[0]
        vertex_j = edge[1]
        weight   = edge[2]

        ## get the cost vector for the current hop constraint (which we are
        ## in the process of calculating), and the cost vector for
        ## the previous hop constraint (which we assume has been calculated
        ## previously)
        Dp = D[hop_constraint - 1]
        Dc = D[hop_constraint]

        ## perform relaxation for the current number of hops against the
        ## cost vector for the previous number of hops; this results
        ## in the current cost vector being the minimum cost using at most
        ## one more hop(edge)
        if Dc[vertex_j] < (Dp[vertex_i] + weight):
            Dc[vertex_j]               = Dp[vertex_i] + weight
            P[hop_constraint,vertex_j] = vertex_i
            T[hop_constraint][vertex_j]= edge

    def HCSSSP_path_iter(self, V, D, P, T, hop_constraint):
        """Iterate over the path from beginning to end yielding the tuple:
        (hi, hj, edge) where hi is the row index (for D,P,T) of vertex
        i in edge, and hj is the row index (should be hi+1) of vertex j
        in edge.
        """
        edge_list  = []
        num_vertex = len(D[0])
        
        ## start at the destination vertex
        curr_v = num_vertex - 1
        h      = hop_constraint

        while curr_v>0:
            prev_vertex  = P[h,curr_v]
            edge         = T[h][curr_v]
            curr_v       = prev_vertex
            h            -= 1

            edge_list.append((h, h+1, edge))
            
        edge_list.reverse()
        for edge in edge_list:
            yield edge

