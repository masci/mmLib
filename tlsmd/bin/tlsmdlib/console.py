## TLS Motion Determination (TLSMD)
## Copyright 2002-2006 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import sys

def kvformat(key, value):
    stdoutln(key.ljust(40, ".") + ": " + str(value))

def endln():
    stdout("\n")

def stdoutln(line):
    stdout(line + "\n")

def stdout(text):
    sys.stdout.write(text)
    sys.stdout.flush()

def stderrln(line):
    stderr(line + "\n")

def stderr(text):
    sys.stdout.write(text)
    sys.stdout.flush()
