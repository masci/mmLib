## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

"""Open file class used by mmLib."""

import os
import types
import gzip


def OpenFile(path, mode):
    """Right now this only supports opening GZip'ed files, in the future
    it might be extended for URLs."""

    ## if path is not a string, assume it is a file object and
    ## return it
    if type(path) != types.StringType:
        return path

    (base, ext) = os.path.splitext(path)

    if ext == ".gz":
        return gzip.open(path, mode)

    return open(path, mode)
    

