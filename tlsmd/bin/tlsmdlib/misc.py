## TLS Motion Determination (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import time

_STIME = 0.0

def start_timing():
    global _STIME
    _STIME = time.time()

def end_timing():
    global _STIME
    tm = time.time() - _STIME
    return "Computation Time: %5.2f sec" % (tm)

def begin_chain_timing(chain_id):
    print "BEGIN TIMING CHAIN %s %f" % (chain_id, time.time())

def end_chain_timing(chain_id):
    print "END TIMING CHAIN %s %f" % (chain_id, time.time())

def rgb_f2i(rgb):
    """Transforms the float 0.0-1.0 RGB color values to
    integer 0-255 RGB values.
    """
    r, g, b = rgb
    ri = int(255.0 * r)
    gi = int(255.0 * g)
    bi = int(255.0 * b)
    return (ri, gi, bi)

def rgb_f2s(rgbf):
    rgbs = "#%2x%2x%2x" % rgb_f2i(rgbf)
    rgbs = rgbs.replace(" ", "0")
    return rgbs
