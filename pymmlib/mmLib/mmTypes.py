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

## turn on debugging
_DEBUG = False


## useful constents
PI       = math.pi
PI2      = math.pi**2
PI3      = math.pi**3

RAD2DEG  = 180.0 / PI
DEG2RAD  = PI / 180.0
RAD2DEG2 = RAD2DEG**2
DEG2RAD2 = DEG2RAD**2

## converting between U (angstrom^2) temp factor
## values and B temp factor values
U2B = 8.0 * PI2
B2U = 1.0 / (8.0 * PI2)


## types, functions, things which are useful and difficult to
## classify...

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


def setmap(smap, skey, dmap, dkey):
    """Sets the dmap/dkey with the value from smap/skey/
    """
    if smap.has_key(skey):
        dmap[dkey] = str(smap[skey])
        return True
    return False


def setmaps(smap, skey, dmap, dkey):
    """Sets the dmap/dkey with the string value from smap/skey/
    """
    if smap.has_key(skey):
        try:
            dmap[dkey] = str(smap[skey])
        except ValueError:
            print "setmaps(): ValueError"
            return False
        return True
    return False


def setmapi(smap, skey, dmap, dkey):
    """Sets the dmap/dkey with the integer value from smap/skey.
    """
    if smap.has_key(skey) and smap[skey]!="":
        try:
            dmap[dkey] = int(smap[skey])
        except ValueError:
            print "setmapi(): ValueError"
            return False
        return True
    return False


def setmapf(smap, skey, dmap, dkey):
    """Sets the dmap/dkey with the float value from smap/skey or
    default if not smap/skey value is found.
    """
    if smap.has_key(skey) and smap[skey]!="":
        try:
            dmap[dkey] = float(smap[skey])
        except ValueError:
            print "setmapf(): ValueError dmap[%s]=smap[%s]=%s" % (
                dkey, skey, smap[skey])
            return False
        return True
    return False


def setmapsd(smap, skey, dmap, dkey):
    """Sets the dmap/dkey with the string value from smap/skey or
    default if not smap/skey value is found.
    """
    try:
        dmap[dkey] = str(smap[skey])
    except ValueError:
        pass
    except KeyError:
        pass

    if default==None:
        return False

    dmap[dkey] = default
    return True


def setmapid(smap, skey, dmap, dkey, default=None):
    """Sets the dmap/dkey with the integer value from smap/skey or
    default if not smap/skey value is found.
    """
    try:
        dmap[dkey] = int(smap[skey])
        return
    except ValueError:
        pass
    except KeyError:
        pass

    if default==None:
        return False

    dmap[dkey] = default
    return True


def setmapfd(smap, skey, dmap, dkey, default=None):
    """Sets the dmap/dkey with the float value from smap/skey or
    default if not smap/skey value is found.
    """
    try:
        dmap[dkey] = float(smap[skey])
        return
    except ValueError:
        pass
    except KeyError:
        pass

    if default==None:
        return False

    dmap[dkey] = default
    return True


def fatal(x):
    """Fatal errors.
    """
    sys.stderr.write("[MMLIB:FATAL] %s\n" % (x))
    sys.exit(-1)
    

def warning(x):
    """Writes warnings out to the file given in the environment variable
    MMLIB_WARNING.  This can be set to a file path, "stdout", "stderr",
    or a empty string for no action.  It writes to the file
    mmlib_warning.txt by default.  
    """
    x = "[MMLIB:WARNING] %s\n" % (x)
    path = os.environ.get("MMLIB_WARNING", "stderr")

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
    if _DEBUG==False:
        return
    
    x    = "[MMLIB:DEBUG] %s\n" % (x)
    path = os.environ.get("MMLIB_DEBUG", "stderr")

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
    
