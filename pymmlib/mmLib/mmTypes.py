## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

"""Specialized types used by mmLib.  In some cases the types are
custom, and the code for those is here.  Inother cases the types are
imported from other Python packages."""

from   __future__           import generators
import string

## we need weak references
import weakref

## get all the Python types
from   types                import *

## use Vector classes from the Scientific Python package
from   Scientific.Geometry  import Vector

## use Numeric Python (Numeric, LinearAlgebra) for array(matrix)
from   Numeric              import array, transpose, matrixmultiply
from   LinearAlgebra        import eigenvalues


class OrderedList(object):
    """Implements a ordered Python list."""
    def __init__(self):
        self.list = []
    def __len__(self):
        return len(self.list)
    def __getitem__(self, i):
        assert type(i) == IntType
        return self.list[i]
    def __delitem__(self, i):
        assert type(i) == IntType
        del self.list[i]
    def __iter__(self):
        for val in self.list: yield val
    def __contains__(self, val):
        return val in self.list
    def count(self, val):
        return self.list.count(val)
    def index(self, val):
        return self.list.index(val)
    def remove(self, val):
        self.list.remove(val)
    def add(self, val):
        self.list.append(val)

class OrderedTupleList(object):
    """Implements a ordered Python list using tuples for list elements.  The
    list is sorted according to Python tuple sorting rules, and the last
    element of the tuple is always the list data.  The add method must be
    overridden to add the proper sorting tuple."""
    def __init__(self):
        self.list = []
    def __len__(self):
        return len(self.list)
    def __getitem__(self, i):
        assert type(i) == IntType
        return self.list[i][-1]
    def __delitem__(self, i):
        assert type(i) == IntType
        del self.list[i]
    def __iter__(self):
        for item in self.list: yield item[-1]
    def __contains__(self, val):
        for item in self.list:
            if val == item[-1]: return True
        return False
    def count(self, val):
        cnt = 0
        for item in self.list:
            if val == item[-1]: cnt += 1
        return cnt
    def index(self, val):
        for item in self.list:
            if val == item[-1]: return self.list.index(item)
        raise ValueError, "list.index(x): x not in list"
    def remove(self, val):
        if item in self.list:
            if val == item[-1]: self.list.remove(item)
        raise ValueError, "list.index(x): x not in list"
    def add(self, val):
        self.list.append((val,))

class WeakrefList:
    """Implements a Python list, but it keeps weak references to the objects
    it contains.  Otherwise, it's exactly like the native list."""
    def __init__(self):
        self.__list = []
    def __len__(self):
        return len(self.__list)
    def __getitem__(self, i):
        return self.__list[i]()
    def __setitem__(self, i, val):
        self.__list[i] = weakref.ref(val, self.__del_cb)
    def __iter__(self):
        for ref in self.__list: yield ref()
    def __contains__(self, val):
        return weakref.ref(val) in self.__list
    def __del_cb(self, ref):
        self.__list.remove(ref)
    def len(self):
        return len(self.__list)
    def append(self, val):
        self.__list.append(weakref.ref(val, self.__del_cb))
    def count(self, val):
        return self.__list.count(weakref.ref(val))
    def insert(self, i, val):
        self.__list.insert(i, weakref.ref(val, self.__del_cb))
    def index(self, val):
        return self.__list.index(weakref.ref(val))
    def pop(self):
        return self.__list.pop(self)()
    def remove(self, val):
       self.__list.remove(weakref.ref(val))
    def sort(self):
        def wcmp(a, b):
            return cmp(a(), b())
        self.__list.sort(wcmp)

## Fragment IDs are tricky things
class FragmentID(object):
    """A fragment ID class acts a lot like a string, but separates the
    res_seq and icode internally."""
    def __init__(self, frag_id):
        try:
            self.res_seq  = int(frag_id)
        except ValueError:
            self.res_seq  = int(frag_id[:-1])
            self.icode    = frag_id[-1]
        else:
            self.icode    = ""
    def __str__(self):
        return str(self.res_seq) + self.icode
    def __lt__(self, other):
        assert isinstance(other, FragmentID)
        return (self.res_seq, self.icode) < (other.res_seq, other.icode)
    def __le__(self, other):
        assert isinstance(other, FragmentID)
        return (self.res_seq, self.icode) <= (other.res_seq, other.icode)
    def __eq__(self, other):
        assert isinstance(other, FragmentID)
        return (self.res_seq, self.icode) == (other.res_seq, other.icode)
    def __ne__(self, other):
        assert isinstance(other, FragmentID)
        return (self.res_seq, self.icode) != (other.res_seq, other.icode)
    def __gt__(self, other):
        assert isinstance(other, FragmentID)
        return (self.res_seq, self.icode) > (other.res_seq, other.icode)
    def __ge__(self, other):
        assert isinstance(other, FragmentID)
        return (self.res_seq, self.icode) >= (other.res_seq, other.icode)
    
### <TESTING>
def wrl_test():
    class C:
        def __init__(self, i):
            self.i = i
        def __del__(self):
            print "C:%d DELETED" % (self.i)
        def __str__(self):
            return "CLASS C:%d" % (self.i)

    l2 = []
    l = WeakrefList()

    for x in range(10):
        c = C(x)
        l.append(c)
        l2.append(c)

    print "# l.index(c)"
    print l.index(c)

    print "# l.pop()"
    print l.pop()

    print "---"
    for x in l: print x
    print "---"

    print "# del l"
    del l

if __name__ == "__main__":
    wrl_test()
    print "# exit"
### </TESTING>
    
