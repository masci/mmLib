## TLS Motion Determination (TLSMD)
## Copyright 2002-2010 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import itertools

class StringTable(object):
    """Converts a given string (as an object of nrows, ncols, and actual
    string) to a tabular format. Intended for use in plotting programs (e.g.,
    gnuplot).
    """
    def __init__(self, m, n,
                 default = "",
                 spacing = 5,
                 min_column_width = 2,
                 title = None,
                 column_titles = None):

        assert m >= 0 and n >= 0
        assert column_titles is None or len(column_titles) <= n

        self.__data = [[default for j in xrange(n)] for i in xrange(m)]
        self.__spacing = spacing
        self.__min_column_width = min_column_width
        self.__title = title
        self.__column_titles = column_titles
        self.__column_width = [min_column_width for i in xrange(n)]
        if column_titles is not None:
            for i in xrange(len(column_titles)):
                title_width = len(column_titles[i]) + 1
                self.__column_width[i] = max(title_width, self.__column_width[i])

    def __setitem__(self, location, value):
        i, j = location
        value = str(value)
        self.__data[i][j] = value
        self.__column_width[j] = max(self.__column_width[j], len(value))

    def __getitem__(self, location):
        i, j = location
        return self.__data[i][j]

    def __str__(self):
        l = []
        if self.__title is not None:
            for ln in self.__title.split("\n"):
                l.append("#")
                l.append(ln)
                l.append("\n")

        field_width = [w for w in self.__column_width]
        space = " " * self.__spacing

        if self.__column_titles is not None:
            first = True
            for value, width in itertools.izip(self.__column_titles, field_width):
                if first:
                    l.append("#")
                    l.append(value.center(width-1, "."))
                    first = False
                else:
                    l.append(space)
                    l.append(value.center(width, "."))
            l.append("\n")

        field_width = [w  for w in self.__column_width]
        for row in self.__data:
            first = True
            for value, width in itertools.izip(row, field_width):
                if first:
                    first = False
                else:
                    l.append(space)
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

    def iter_rows(self):
        for row in self.__data:
            yield row

    def iter_column_titles(self):
        if self.__column_titles is not None:
            return iter(self.__column_titles)
        else:
            return iter(list())

    def title(self):
        if self.__title is not None:
            return self.__title
        else:
            return ""

    def size(self):
        return len(self.__data), max(itertools.imap(len, self.__data))


def StringTableFromMatrix(matrix):
    """Return a StringTable instance with the dimentions and values of the 
    argument matrix.
    """
    m, n = matrix.shape
    tbl = StringTable(m, n, "0.0")
    for i in range(m):
        for j in range(n):
            tbl[i,j] = matrix[i,j]
    return tbl


## testing
def testmain():
    t = StringTable(0, 2, title = "My Test Table", 
                    column_titles = ["Column 1", "Column 2"])

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
