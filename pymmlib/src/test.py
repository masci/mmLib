# speed test for the C PDB reading module -- it looks like it's about
# 450% faster than the Python implementation
import sys
import os
import time
import pdbmodule

try:
    path = sys.argv[1]
except IndexError:
    sys.exit(1)

sec = time.time()
listx = pdbmodule.read(path)
sec = time.time() - sec

print "%s: %d records in %.2f seconds" % (path, len(listx), sec)

#for x in listx:
#    print x
