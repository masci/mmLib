## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Specialized types used by mmLib.  In some cases the types are
custom, and the code for those is here.  Inother cases the types are
imported from other Python packages.
"""
from __future__ import generators
import os
import sys
import string
import gzip
from types import *
from Scientific.Geometry  import Vector
from Numeric import *
from LinearAlgebra import *


def OpenFile(path, mode):
    """Right now this only supports opening GZip'ed files, in the future
    it might be extended for URLs.
    """
    ## if path is not a string, assume it is a file object and
    ## return it
    if type(path) != StringType:
        return path

    (base, ext) = os.path.splitext(path)

    if ext == ".gz":
        return gzip.open(path, mode)

    return open(path, mode)


class FragmentID(object):
    """A fragment ID class acts a lot like a string, but separates the
    res_seq and icode internally.
    """
    def __init__(self, frag_id):
        self.res_seq = 1
        self.icode = ""
        try:
            self.res_seq = int(frag_id)
        except ValueError:
            try:
                self.res_seq = int(frag_id[:-1])
            except ValueError:
                pass
            else:
                self.icode = frag_id[-1]
    def __str__(self):
        return str(self.res_seq) + self.icode
    def __lt__(self, other):
        assert isinstance(other, FragmentID)
        return (self.res_seq, self.icode) < (other.res_seq, other.icode)
    def __le__(self, other):
        assert isinstance(other, FragmentID)
        return (self.res_seq, self.icode) <= (other.res_seq, other.icode)
    def __eq__(self, other):
        assert isinstance(other, FragmentID)
        return (self.res_seq, self.icode) == (other.res_seq, other.icode)
    def __ne__(self, other):
        assert isinstance(other, FragmentID)
        return (self.res_seq, self.icode) != (other.res_seq, other.icode)
    def __gt__(self, other):
        assert isinstance(other, FragmentID)
        return (self.res_seq, self.icode) > (other.res_seq, other.icode)
    def __ge__(self, other):
        assert isinstance(other, FragmentID)
        return (self.res_seq, self.icode) >= (other.res_seq, other.icode)


def warning(x):
    """Writes warnings out to the file given in the environment variable
    MMLIB_WARNING, or mmlib_warning.txt by default.
    """
    path = os.environ.get("MMLIB_WARNING", "mmlib_warning.txt")
    open(path, "a").write(x+"\n")


def debug(x):
    """Writes debugging output to the file given in the environment variable
    MMLIB_DEBUG, or mmlib_debug.txt by default.
    """
    path = os.environ.get("MMLIB_DEBUG", "mmlib_debug.txt")
    open(path, "a").write(x+"\n")


### <TESTING>
if __name__ == "__main__":
    fid = FragmentID("163A")
    print "fid",fid
    print "fid.res_seq", fid.res_seq
    print "fid.icode", fid.icode
### </TESTING>
    
