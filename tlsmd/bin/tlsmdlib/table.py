## TLS Motion Determination (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import itertools

class StringTable(object):
    def __init__(self, m, n, default = "", spacing = 2, min_column_width = 2):
        assert m >= 0 and n >= 0
        self.__data = [[default for j in xrange(n)] for i in xrange(m)]
        self.__spacing = spacing
        self.__min_column_width = min_column_width
        self.__column_width = [min_column_width for i in xrange(n)]

    def __setitem__(self, location, value):
        i, j = location
        value = str(value)
        self.__data[i][j] = value
        self.__column_width[j] = max(self.__column_width[j], len(value))

    def __getitem__(self, location):
        i, j = location
        return self.__data[i][j]
        
    def __str__(self):
        field_width = [w + self.__spacing for w in self.__column_width]
        l = []
        for row in self.__data:
            for value, width in itertools.izip(row, field_width):
                l.append(value.rjust(width))
            l.append("\n")
        return "".join(l)

    def set_column(self, i, j, iterable):
        dataiter = iter(self.__data)
        for x in xrange(i): dataiter.next()
        width = self.__min_column_width
        for row, value in itertools.izip(self.__data, iterable):
            value = str(value)
            row[j] = value
            width = max(width, len(value))
        self.__column_width[j] = max(self.__column_width[j], width)
            
    def append_row(self, *row_tuple):
        row = []
        self.__data.append(row)
        for j, value in enumerate(row_tuple):
            value = str(value)
            row.append(value)
            vlen = len(value)
            try:
                self.__column_width[j] = max(self.__column_width[j], vlen)
            except IndexError:
                nadd = j - len(self.__column_width) + 1
                for x in xrange(nadd):
                    width = max(self.__min_column_width, vlen)
                    self.__column_width.append(width)
        
    def size(self):
        return len(self.__data), max(itertools.imap(len, self.__data))


def StringTableFromMatrix(matrix):
    """Return a StringTable instance with the dimentions and
    values of the argument matrix.
    """
    m, n = matrix.shape
    tbl = StringTable(m, n, "0.0")
    for i in range(m):
        for j in range(n):
            tbl[i,j] = matrix[i,j]
    return tbl


## testing
def testmain():
    t = StringTable(0, 0)

    for i in xrange(10):
        t.append_row(i, "stuff")

    t.set_column(0, 1, itertools.imap(lambda x: str(x)+"*", xrange(10)))

    print "size", t.size()
    print str(t)
    print "============="
    
    t = StringTable(5, 5, "*")
    t[3,4] = 5.0
    print str(t)

if __name__ == "__main__":
    testmain()

