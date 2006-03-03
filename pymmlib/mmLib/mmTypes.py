## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Specialized types used by mmLib.
"""
from __future__ import generators

import os
import sys

## turn on debugging
_DEBUG = False

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
        return "%d%s" % (self.res_seq, self.icode)
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
