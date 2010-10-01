## TLS Motion Determination (TLSMD)
## Copyright 2002-2010 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

## Python modules
import sys

## TLSMD
import misc

console_output_enabled = True
debug_output_enabled = False
cpu_time_output_enabled = True

import traceback
def formatExceptionInfo(maxTBlevel=7):
    """Takes the three-element tuple returned by sys.exc_info() and transforms
    each element into a more convenient form, a string. cla.__name__ gives
    the name of the exception class, while exc.__dict__["args"] gives other
    details about the exception.
    In the case of socket exceptions, these details will be in a two-element
    tuple, like:
        ("error", (32, 'Broken pipe').
    Lastly, traceback.format_tb() formats the traceback information into a
    string.
    The optional argument (maxTBlevel> in the sample code) allows users to
    control the depth of the traceback that will be formatted.
    The traceback information is not essential to identify or categorize
    exceptions, but if you want to log all the spurious unknown exceptions
    your program encounters, it is useful to write that traceback string in
    the log.
    """
    cla, exc, trbk = sys.exc_info()
    excName = cla.__name__
    try:
        excArgs = exc.__dict__["args"]
    except KeyError:
        excArgs = "<no args>"
    excTb = traceback.format_tb(trbk, maxTBlevel)
    return (excName, excArgs, excTb) ## <- Use this one for _very_ verbose!
    #return excTb ## <- Less verbose.

def stdout(text):
    """Simple STDOUT printing redirected to log.txt
    """
    if console_output_enabled:
        sys.stdout.write(text)
        sys.stdout.flush()

def stderr(text):
    """Simple STDERR printing redirected to log.txt
    """
    if console_output_enabled:
        sys.stdout.write(text)
        sys.stdout.flush()

def enable():
    global console_output_enabled
    console_output_enabled = True

def disable():
    global console_output_enabled
    console_output_enabled = False

def kvformat(key, value):
    """TOC-style of printing to stdout/log.txt file
    """
    ## E.g., "TLS PARAMETER FIT ENGINE................: ISOT"
    stdoutln(key.ljust(40, ".") + ": " + str(value))

## The following three defs are to force Unix/Linux newline usage.
def endln():
    stdout("\n")

def stdoutln(line):
    stdout("[%s] %s\n" % (misc.timestamp(), line))

def debug_stdoutln(line):
    if debug_output_enabled:
        stdout("[%s] %s\n" % (misc.timestamp(), line))

## TODO: Add warning function, 2009-01-08
def warning_stdoutln(line):
    stdout("[%s]     Warning: %s\n" % (misc.timestamp(), line))

def cpu_time_stdoutln(line):
    if cpu_time_output_enabled:
        stdout("[%s] CPU_TIME %s\n" % (misc.timestamp(), line))

def stderrln(line):
    stderr(line + "\n")
