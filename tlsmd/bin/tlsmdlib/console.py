## TLS Motion Determination (TLSMD)
## Copyright 2002-2008 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import sys

console_output_enabled = True

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
    stdout(line + "\n")

def stderrln(line):
    stderr(line + "\n")
