#!/usr/bin/env python

## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

## Distutils script to ease installation
##
## python setup.py sdist           : creates a .tar.gz distribution
## python setup.py bdist_wininst   : creates a window binary installer
## python setup.py bdist_rpm       : creates a rpm distribution


from distutils.core import setup

setup(   name         = "PyMMLib",
         version      = "0.1",
         author       = "Jay Painter",
         author_email = "jpaint@u.washington.edu",
         url          = "http://pymmlib.sourceforge.net/",

         packages     = ["mmLib"] )
