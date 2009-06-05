## TLS Motion Determination (TLSMD)
## Copyright 2002-2009 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

## Python modules
import time
import datetime
import random
import string
import re

## TLSMD
import console
import conf

_STIME = 0.0

def start_timing():
    global _STIME
    _STIME = time.time()

def end_timing():
    global _STIME
    tm = time.time() - _STIME
    return "Computation Time: %5.2f sec" % (tm)

def timestamp():
    ## TODO: Allow for either "-7" or "-10". 2009-03-25
    ## Also, allow for "Y-M-D" format, etc.
    return datetime.datetime.fromtimestamp(time.time()).isoformat(' ')[:-7]

def start_time():
    ## returns format e.g., "26 May 2009"
    return time.strftime("%d %b %Y", time.localtime(conf.globalconf.start_time))

def begin_chain_timing(chain_id):
    console.stdoutln("BEGIN TIMING CHAIN %s %f" % (chain_id, time.time()))

def end_chain_timing(chain_id):
    console.stdoutln("END TIMING CHAIN %s %f" % (chain_id, time.time()))

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

def generate_security_code(code_length = 8):
    """Generates a random 8-character string
    """
    random.seed()
    codelist = list(5 * string.ascii_letters)
    random.shuffle(codelist)
    code = "".join(random.sample(codelist, code_length))
    return code

def parse_chains(chains):
    """Parses a given chains string returned from the MySQL DB
    """
    ## Turns "A:100:1:aa" -> "A", "100", "1", "aa"
    ## I.e.: chain_id = A, num_residues = 100,
    ##       selected = 1/True, type = aa/amino acids
    chid     = re.sub(r'([A-Za-z0-9]):([0-9]{1,}):([01]):([na]{2});?', '\\1;', chains).rstrip(";")
    length   = re.sub(r'([A-Za-z0-9]):([0-9]{1,}):([01]):([na]{2});?', '\\2;', chains).rstrip(";")
    selected = re.sub(r'([A-Za-z0-9]):([0-9]{1,}):([01]):([na]{2});?', '\\3;', chains).rstrip(";")
    type     = re.sub(r'([A-Za-z0-9]):([0-9]{1,}):([01]):([na]{2});?', '\\4;', chains).rstrip(";")
    return chid, length, selected, type
