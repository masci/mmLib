## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

"""MTZ file parser.  It has a few flaws I'm going to fix soon.  It will not
correctly convert files saved on machines with different endian formats or
integer/float sizes.  It cannot save MTZ files.  If anyone needs these
features, email jpaint@u.washington.edu and I'll finish them. ;)"""

import types
import struct
from   FileIO  import OpenFile


## Error
MTZError = "MTZError"



## host system variable sizes 
INT_SIZE   = struct.calcsize("i")
FLOAT_SIZE = struct.calcsize("f")
WORD_SIZE  = struct.calcsize("f")



## MTZ header keywords
MTZKeywords = [
    "VERS",
    "TITLE",
    "NCOL",
    "CELL",
    "SORT",
    "SYMINF",
    "SYMM",
    "RESO",
    "VALM",
    "COL",
    "NDIF",
    "PROJECT",
    "CRYSTAL",
    "DATASET",
    "DCELL",
    "DWAVEL",
    "BATCH",
    "END"]


class MTZFile:
    """Reads a CCP4 MTZ file."""

    def __init__(self):
        self.int_size     = INT_SIZE
        self.float_size   = FLOAT_SIZE
        self.word_size    = WORD_SIZE

        self.header1_list = []
        self.header2_list = []

        self.column_list  = []
        self.data_list    = []


    def __str__(self):
        return str(self.column_list)
    

    def get_data_list(self):
        return self.data_list


    def load_file(self, fil):
        self.fil = OpenFile(fil, "rb")

        ## check MTZ file header
        self.fil.seek(0)
        if self.fil.read(4) != "MTZ ":
            raise MTZError, "File is not a MTZ file"

        ## read the offset for the file header and seek to
        ## the begging of the header and read it
        self.fil.seek(4)
        hoffset = self.get_offset(self.read_int() - 1)
        self.fil.seek(hoffset)

        header = self.fil.read()

        ## build primary/secondary header lists
        self.header1_list = []
        self.header2_list = []

        hflag = 0

        while header:
            line   = header[:80]
            header = header[80:]
            ltmp   = line.split()

            if not ltmp:
                raise MTZError, "Blank header record"

            if hflag == 0: self.header1_list.append(ltmp)
            else:          self.header2_list.append(ltmp)
            
            if ltmp[0] == "END": hflag = 1

        ## build a list of the data columns
        self.column_list = []
        
        for rlist in self.header1_list:
            rname = rlist[0]
            if rname[:3] != "COL": continue

            try:
                cname = rlist[1]
            except IndexError:
                raise MTZError, "Malformed MTZ COL header"
            
            self.column_list.append(cname)

        ## data rows start at 20*sizeof(int) into the file
        self.data_list = []

        doffset    = self.get_offset(20)
        ncol       = len(self.column_list)
        row_format = "f" * ncol ## This needs to be equiv to REAL*4
        row_bytes  = 4 * ncol
        nrow       = (hoffset - doffset) / row_bytes

        self.fil.seek(doffset)
        
        for i in range(nrow):
            row = struct.unpack(row_format, self.fil.read(row_bytes))
            self.data_list.append(row)

    
    def save_file(self, fil):
        pass


    def get_offset(self, x):
        """Seeking in a MTZ file is a painful experience.  The documentation
        for the MTZ files gives a completely incorrect description of them,
        so it's best to ignore it.  Seeks are actually in multiples of the
        size of integers on the machine the MTZ file was written in."""
        return x * self.int_size


    def read_int(self):
        i = self.fil.read(self.int_size)
        (i,) = struct.unpack("i", i)
        return i


    def write_int(self, i):
        s = struct.pack("i", i)
        self.fil.write(s)


    def read_float(self):
        f = self.fil.read(self.float_size)
        (f,) = struct.unpack("f", f)
        return f



if __name__ == "__main__":
    import sys

    try:
        path = sys.argv[1]
    except IndexError:
        print "usage: MTZ.py <MTZ path>"
        sys.exit(1)
    
    mtz = MTZFile()
    mtz.load_file(path)
    print mtz

    for row in mtz.get_data_list():
        print row
