## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

class XRay(object):
    """
    """
    def __init__(self):
        self.unit_cell = None
        
        self.resolution_high = None
        self.resolution_low  = None

        self.R_free = None
        self.R_work = None
        self.R_obs  = None

        self.number_reflections_obs  = None
        self.number_reflections_work = None
        self.number_reflections_free = None
