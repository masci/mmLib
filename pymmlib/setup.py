#!/usr/bin/env python

## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

## Distutils script to ease installation
## python setup.py doc             : runs happydoc creating "doc" dir
## python setup.py sdist           : creates a .tar.gz distribution
## python setup.py bdist_wininst   : creates a window binary installer
## python setup.py bdist_rpm       : creates a rpm distribution


import os
import sys
from distutils.core import setup



def run_setup():
    setup( name         = "PyMMLib",
           version      = "0.2",
           author       = "Jay Painter",
           author_email = "jpaint@u.washington.edu",
           url          = "http://pymmlib.sourceforge.net/",
           
           packages     = ["mmLib", "mmLib/Data", "mmLib/Extensions"] )


def make_doc():
    os.system("happydoc -d doc -t 'PyMMLib Documentation' --no-comments "\
              "--no-private-names mmLib")

    

if __name__ == "__main__":
    if sys.argv[1] == "doc":
        make_doc()
    else:
        run_setup()
