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
import math
import string
import gzip
from types import *
from Scientific.Geometry  import Vector
from Numeric import *
from LinearAlgebra import *

## useful constents
rad2deg  = 180.0 / math.pi
deg2rad  = math.pi / 180.0
rad2deg2 = rad2deg * rad2deg
deg2rad2 = deg2rad * deg2rad


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


def setmaps(smap, skey, dmap, dkey, default = None):
    """Sets the dmap/dkey with the string value from smap/skey or
    default if not smap/skey value is found.
    """
    try:
        dmap[dkey] = str(smap[skey])
    except ValueError:
        pass
    except KeyError:
        pass

    if default == None:
        return False

    dmap[dkey] = default
    return True


def setmapi(smap, skey, dmap, dkey, default = None):
    """Sets the dmap/dkey with the integer value from smap/skey or
    default if not smap/skey value is found.
    """
    try:
        dmap[dkey] = int(smap[skey])
    except ValueError:
        pass
    except KeyError:
        pass

    if default == None:
        return False

    dmap[dkey] = default
    return True


def setmapf(smap, skey, dmap, dkey, default = None):
    """Sets the dmap/dkey with the float value from smap/skey or
    default if not smap/skey value is found.
    """
    try:
        dmap[dkey] = float(smap[skey])
    except ValueError:
        pass
    except KeyError:
        pass

    if default == None:
        return False

    dmap[dkey] = default
    return True


def warning(x):
    """Writes warnings out to the file given in the environment variable
    MMLIB_WARNING.  This can be set to a file path, "stdout", "stderr",
    or a empty string for no action.  It writes to the file
    mmlib_warning.txt by default.  
    """
    x    = x + "\n"
    path = os.environ.get("MMLIB_WARNING", "mmlib_warning.txt")

    try:
        if path == "":
            return
        elif path == "stdout":
            sys.stdout.write(x)
        elif path == "stderr":
            sys.stderr.write(x)
        else:
            open(path, "a").write(x)
    except IOError:
        pass


def debug(x):
    """Writes debugging output to the file given in the environment variable
    MMLIB_DEBUG.  This can be set to a file path, "stdout", "stderr",
    or a empty string for no action.  It writes to the file
    mmlib_warning.txt by default.
    """
    x    = x + "\n"
    path = os.environ.get("MMLIB_DEBUG", "mmlib_debug.txt")

    try:
        if path == "":
            return
        elif path == "stdout":
            sys.stdout.write(x)
        elif path == "stderr":
            sys.stderr.write(x)
        else:
            open(path, "a").write(x)
    except IOError:
        pass


### <TESTING>
if __name__ == "__main__":
    fid = FragmentID("163A")
    print "fid",fid
    print "fid.res_seq", fid.res_seq
    print "fid.icode", fid.icode
### </TESTING>
    
