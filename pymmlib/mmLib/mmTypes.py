## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

"""Provides data types used across mmLib.  In some cases the types are
custom, and the code for those is here.  Inother cases the types are imported
from other Python packages."""



## use Vector classes from the Scientific Python package
from   Scientific.Geometry  import Vector



## use Numeric Python (Numeric, LinearAlgebra) for array(matrix)
from   Numeric              import array
from   LinearAlgebra        import eigenvalues

