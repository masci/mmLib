## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import struct


class MTZFile:
    """Reads a CCP4 MTZ file."""
    size_of_int = struct.calcsize("i")
    size_of_float = struct.calcsize("f")
    size_of_word = struct.calcsize("f")
    def __init__(self, path):
        self.path = path

    def read_int(self, fil):
        i = fil.read(self.size_of_int)
        (i,) = struct.unpack("i", i)
        return i

    def read_float(self, fil):
        f = fil.read(self.size_of_float)
        (f,) = struct.unpack("f", f)
        return f

    def isMTZFile(self):
        """Checks header to make sure the file is a MTZ file."""
        fil = open(self.path, "rb")
        header = fil.read(4)
        if header != "MTZ ":
            return 0
        return 1

    def getHeaderOffset(self):
        fil = open(self.path, "rb")
        fil.seek(4)
        offset = fil.read(self.size_of_int)
        (offset,) = struct.unpack("i", offset)
        return (offset * self.size_of_int) - (1 * self.size_of_int)

    def readHeader(self):
        ho = self.getHeaderOffset()
        fil = open(self.path, "rb")
        fil.seek(ho)

        vers = fil.read(15)
        print 'VERS="%s"' %(vers)

        title = fil.read(70)
        print 'TITLE="%s"'%(title)




if __name__ == "__main__":
    mtz = MTZFile("/home/jpaint/shadow/data/caa_unique.mtz")
    print mtz.isMTZFile()
    mtz.readHeader()
