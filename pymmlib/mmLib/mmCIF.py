## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""mmCIF file and mmCIF dictionary parser.  Files are parsed
into a set of data structures where they can be further processed.  The
data structures can also be constructed and written back out as mmCIF.
A CIF dictionary parser is also included as a specialized version of the
mmCIF parser.
"""
from __future__ import generators
import string
import re
import types
from mmTypes import *

##
## DATA STRUCTURES FOR HOLDING CIF INFORMATION
##
## mmCIF files are parsed into:
##         mmCIFFile -> [mmCIFData] -> [mmCIFTable] -> [mmCIFRow]
##
## mmCIF dictionaries are parsed into:
##         mmCIFDictionary -> [mmCIFData] -> [mmCIFTable] -> [mmCIFRow]
##

MAX_LINE = 80


class mmCIFError(Exception):
    """Base class of errors raised by Structure objects.
    """
    pass


class mmCIFSyntaxError(Exception):
    """Base class of errors raised by Structure objects.
    """
    def __init__(self, line_num, text):
        self.line_num = line_num
        self.text = text

    def __str__(self):
        return "[line: %d] %s" % (self.line_num, self.text)


class mmCIFRow(dict):
    """Contains one row of data.  In a mmCIF file, this is one complete
    set of data found under a section.  The data can be accessed by using
    the column names as class attributes.
    """
    __slots__ = ["table"]
    
    def __eq__(self, other):
        return id(self) == id(other)

    def __deepcopy__(self, memo):
        return mmCIFRow(self)

    def __contains__(self, column):
        return dict.__contains__(self, column.lower())

    def __setitem__(self, column, value):
        dict.__setitem__(self, column.lower(), value)

    def __getitem__(self, column):
        return dict.__getitem__(self, column.lower())

    def __delitem__(self, column):
        dict.__delitem__(self, column.lower())

    def get(self, column, default = None):
        return dict.get(self, column.lower(), default)

    def has_key(self, column):
        return dict.has_key(self, column.lower())

    def mget(self, *keys):
        """Return the fist value found for the given keys in the argument
        list.
        """
        for key in keys:
            try:
                return self[key]
            except KeyError:
                continue


class mmCIFTable(list):
    """Contains columns and rows of data for a mmCIF section.  Rows of data
    are stored as mmCIFRow classes.
    """
    __slots__ = ["name", "columns", "data"]

    def __init__(self, name, columns = None):
        list.__init__(self)
        self.name = name
        self.columns = columns or []

    def __deepcopy__(self, memo):
        table = mmCIFTable(self.name, self.columns[:])
        for row in self:
            table.append(copy.deepcopy(row, memo))
        return table

    def __eq__(self, other):
        return id(self) == id(other)

    def is_single(self):
        """Return true if the table is not a _loop table with multiple
        rows of data.
        """
        return len(self) <= 1
    
    def __getitem__(self, x):
        """Retrieves mmCIFRow at index x from the table if the argument is
        a integer.  If the argument is a string, then the data from the
        first row is returned.
        """
        if type(x) == IntType:
            return list.__getitem__(self, x)

        elif type(x) == StringType:
            try:
                return self[0][x]
            except IndexError:
                raise KeyError
            except KeyError:
                raise KeyError

        raise TypeError, x
    
    def __setitem__(self, x, value):
        if type(x) == IntType and isinstance(value, mmCIFRow):
            value.table = self
            list.__setitem__(self, x, value)

        elif type(x) == StringType:
            try:
                self[0][x] = value
            except IndexError:
                row = mmCIFRow()
                row[x] = value
                self.append(row)

    def __delitem__(self, i):
        assert isinstance(row, mmCIFRow)
        self.remove(self[i])

    def append(self, row):
        assert isinstance(row, mmCIFRow)
        row.table = self
        list.append(self, row)

    def insert(self, i, row):
        assert isinstance(row, mmCIFRow)
        row.table = self
        list.insert(self, i, row)

    def remove(self, row):
        assert isinstance(row, mmCIFRow)
        del row.table
        list.remove(self, row)

    def autoset_columns(self):
        """Automaticly sets the mmCITable column names by inspecting all
        mmCIFRow objects it contains.
        """
        column_used = {}

        for cif_row in self:
            for column in cif_row.keys():
                column_used[column] = True                
                if column not in self.columns:
                    self.columns.append(column)

        for column in self.columns[:]:
            if not column_used.has_key(column):
                self.columns.remove(column)

    def get_row(self, *args):
        """Preforms a SQL-like 'AND' select aginst all the rows in the table,
        and returns the first matching row found.  The arguments are a
        variable list of tuples of the form:
          (<column-name>, <column-value>)
        For example:
          ger_row(('atom_id','CA'),('entity_id', '1'))
        returns the first matching row with atom_id==1 and entity_id==1.
        """
        for cif_row in self:
            try:
                for (column, value) in args:
                    if cif_row[column] != value:
                        continue
            except KeyError:
                continue

            return cif_row

        return None

    def iter_rows(self, *args):
        """This is the same as get_row, but it iterates over all matching
        rows in the table.
        """
        for cif_row in self:
            try:
                for (column, value) in args:
                    if cif_row[column] != value:
                        continue
            except KeyError:
                continue

            yield cif_row

    def row_index_dict(self, key):
        """Return a dictionary mapping the value of the row's value in
        column 'key' to the row itself.  If there are multiple rows with
        the same key value, they will be overwritten with the last found
        row.
        """
        dictx = {}
        for row in self:
            try:
                dictx[row[key]] = row
            except KeyError:
                pass
        return dictx

    def debug(self):
        print "mmCIFTable::%s" % (self.name)
        for row in self:
            for col in self.columns:
                print "%s=%s" % (col, row.get(col))[:80]
            print "---"


class mmCIFData(list):
    """Contains all information found under a data_ block in a mmCIF file.
    mmCIF files are represented differently here than their file format
    would suggest.  Since a mmCIF file is more-or-less a SQL database dump,
    the files are represented here with their sections as "Tables" and
    their subsections as "Columns".  The data is stored in "Rows".
    """
    __slots__ = ["name", "file"]
    
    def __init__(self, name):
        list.__init__(self)
        self.name = name

    def __deepcopy__(self, memo):
        data = mmCIFData(self.name)
        for table in self:
            data.append(copy.deepcopy(table, memo))
        return data

    def __eq__(self, other):
        return id(self) == id(other)
    
    def __getitem__(self, x):
        if type(x) == IntType:
            return list.__getitem__(self, x)

        elif type(x) == StringType:
            name = x.lower()
            for ctable in self:
                if ctable.name.lower() == name:
                    return ctable
            raise KeyError, x

        raise TypeError, x

    def __setitem__(self, x, table):
        """
        """
        assert isinstance(table, mmCIFTable)

        try:
            old_table = self[x]
        except (KeyError, IndexError):
            pass
        else:
            self.remove(old_table)

        if type(x) == IntType:
            table.data = self
            list.__setitem__(self, x, table)

        elif type(x) == StringType:
            self.append(table)

    def __delitem__(self, x):
        """Remove a mmCIFTable by index or table name.
        """
        self.remove(self[x])

    def append(self, table):
        """Append a mmCIFTable.  This will trigger the removal of any
        table with the same name.
        """
        assert isinstance(table, mmCIFTable)
        try:
            del self[table.name]
        except KeyError:
            pass
        table.data = self
        list.append(self, table)

    def insert(self, i, table):
        assert isinstance(table, mmCIFTable)
        try:
            del self[table.name]
        except KeyError:
            pass
        table.data = self
        list.insert(self, i, table)

    def remove(self, table):
        assert isinstance(table, mmCIFTable)
        del table.data
        list.remove(self, table)

    def has_key(self, x):
        try:
            self[x]
        except KeyError:
            return False
        else:
            return True

    def get(self, x, default = None):
        try:
            return self[x]
        except KeyError:
            return default

    def get_table(self, name):
        """Looks up and returns a stored mmCIFTable class by it's name.  This
        name is the section key in the mmCIF file.
        """
        try:
            return self[name]
        except KeyError:
            return None
        except IndexError:
            return None

    def split_tag(self, tag):
        cif_table_name, cif_column_name = tag[1:].split(".")
        return cif_table_name.lower(), cif_column_name.lower()

    def join_tag(self, cif_table_name, cif_column_name):
        return "_%s.%s" % (cif_table_name, cif_column_name)

    def get_tag(self, tag):
        """Get.
        """
        table_name, column = self.split_tag(tag)
        try:
            return self[table_name][column]
        except KeyError:
            return None

    def set_tag(self, tag, value):
        """Set.x
        """
        table_name, column = self.split_tag(tag)
        self[table_name][column] = value
        
    def debug(self):
        print "mmCIFData::%s" % (self._name)
        for ctable in self:
            ctable.debug()


class mmCIFSave(mmCIFData):
    """Class to store data from mmCIF dictionary save_ blocks.  I treat
    them as non-nested sections along with data_ sections.  This may
    not be correct.
    """
    pass


class mmCIFFile(list):
    """Class representing a mmCIF files.
    """
    def __deepcopy__(self, memo):
        file = mmCIFFile()
        for data in self:
            file.append(copy.deepcopy(data, memo))
        return file

    def __eq__(self, other):
        return id(self) == id(other)

    def __getitem__(self, x):
        """Retrieve a mmCIFData object by index or name.
        """
        if type(x) == IntType:
            return list.__getitem__(self, x)

        elif type(x) == StringType:
            name = x.lower()
            for cdata in self:
                if cdata.name.lower() == name:
                    return cdata
            raise KeyError, x

        raise TypeError, x
    
    def __delitem__(self, x):
        """Remove a mmCIFData by index or data name.  Raises IndexError
        or KeyError if the mmCIFData object is not found, the error raised
        depends on the argument type.
        """
        self.remove(self[x])

    def append(self, cdata):
        """Append a mmCIFData object.  This will trigger the removal of any
        mmCIFData object in the file with the same name.
        """
        assert isinstance(cdata, mmCIFData)
        try:
            del self[cdata.name]
        except KeyError:
            pass
        cdata.file = self
        list.append(self, cdata)

    def insert(self, i, cdata):
        assert isinstance(cdata, mmCIFData)
        try:
            del self[cdata.name]
        except KeyError:
            pass
        cdata.file = self
        list.insert(self, i, cdata)

    def has_key(self, x):
        for cdata in self:
            if cdata.name == x:
                return True
        return False

    def get(self, x, default = None):
        try:
            return self[x]
        except KeyError:
            return default
        
    def load_file(self, fil, update_cb = None):
        """Load and append the mmCIF data from file object fil into self.
        """
        fil = OpenFile(fil, "r")
        mmCIFFileParser().parse_file(fil, self, update_cb)

    def save_file(self, fil):
        fil = OpenFile(fil, "w")
        mmCIFFileWriter().write_file(fil, self)

    def get_data(self, name):
        try:
            return self[name]
        except KeyError:
            return None
        except IndexError:
            return None

    def debug(self):
        print "mmCIFFile"
        for cdata in self:
            cdata.debug()


class mmCIFDictionary(mmCIFFile):
    """Class representing a mmCIF dictionary.  The constructor of this class
    takes two arguments.  The first is the string path for the file, or
    alternativly a file object.
    """
    pass


##
## FILE PARSERS/WRITERS
##


class mmCIFFileParser(object):
    """Stateful parser which uses the mmCIFElementFile tokenizer to read
    a mmCIF file and convert it into the mmCIFData/mmCIFTable/mmCIFRow
    data hierarchy.
    """
    def parse_file(self, fil, cif_file, update_cb = None):
        self.update_cb = update_cb
        self.line_number = 0

        token_iter = self.gen_token_iter(fil)

        try:
            self.parse(token_iter, cif_file)
        except StopIteration:
            pass
        else:
            raise mmCIFError()

    def syntax_error(self, err):
        raise mmCIFSyntaxError(self.line_number, err)
 
    def parse(self, token_iter, cif_file):
        cif_table_cache = {}
        cif_data        = None
        cif_table       = None
        cif_row         = None
        state           = ""

        tblx,colx,strx,tokx = token_iter.next()
        
        while 1:
            if tblx != None:
                state = "RD_SINGLE"
            elif tokx != None:
                if tokx == "loop_":
                    state = "RD_LOOP"
                elif tokx.startswith("data_"):
                    state = "RD_DATA"
                elif tokx.startswith("save_"):
                    state = "RD_SAVE"
                else:
                    self.syntax_error("bad token #1: "+str(tokx))
            else:
                self.syntax_error("bad token #2")
                return            
 

            if state == "RD_SINGLE":
                try:
                    cif_table = cif_table_cache[tblx]
                except KeyError:
                    cif_table = cif_table_cache[tblx] = mmCIFTable(tblx)

                    try:
                        cif_data.append(cif_table)
                    except AttributeError:
                        self.syntax_error(
                            "section not contained in data_ block")
                        return

                    cif_row = mmCIFRow()
                    cif_table.append(cif_row)
                else:
                    try:
                        cif_row = cif_table[0] 
                    except IndexError:
                        self.syntax_error("bad token #3")
                        return

                ## check for duplicate entries
                if colx in cif_table.columns:
                    self.syntax_error("redefined subsection (column)")
                    return
                else:
                    cif_table.columns.append(colx)

                x,x,strx,tokx = token_iter.next()

                if tokx != None:
                    if tokx == ".":
                        cif_row[colx] = ""
                    else:
                        cif_row[colx] = tokx
                elif strx != None:
                    cif_row[colx] = strx
                else:
                    self.syntax_error("bad token #4")

                tblx,colx,strx,tokx = token_iter.next()
                continue

            elif state == "RD_LOOP":
                tblx,colx,strx,tokx = token_iter.next()

                if tblx == None or colx == None:
                    self.syntax_error("bad token #5")
                    return
                
                if cif_table_cache.has_key(tblx):
                    self.syntax_error("_loop section duplication")
                    return

                cif_table = mmCIFTable(tblx)

                try:
                    cif_data.append(cif_table)
                except AttributeError:
                    self.syntax_error(
                        "_loop section not contained in data_ block")
                    return

                cif_table.columns.append(colx)

                while 1:
                    tblx,colx,strx,tokx = token_iter.next()
                    if tblx == None:
                        break
                    if tblx != cif_table.name:
                        self.syntax_error("changed section names in _loop")
                        return
                    cif_table.columns.append(colx)
                    
                while 1:
                    cif_row = mmCIFRow()
                    cif_table.append(cif_row)

                    for col in cif_table.columns:
                        if tokx != None:
                            if tokx == ".":
                                cif_row[col] = ""
                            else:
                                cif_row[col] = tokx
                        elif strx != None:
                            cif_row[col] = strx

                        tblx,colx,strx,tokx = token_iter.next()

                    ## the loop ends when one of these conditions is met
                    if tblx != None:
                        break
                    if tokx == None:
                        continue
                    if tokx == "loop_":
                        break
                    if tokx.startswith("data_"):
                        break
                    if tokx.startswith("save_"):
                        break

                continue

            elif state == "RD_DATA":
                cif_data = mmCIFData(tokx[5:])
                cif_file.append(cif_data)
                cif_table_cache = {}
                cif_table = None

                tblx,colx,strx,tokx = token_iter.next()

            elif state == "RD_SAVE":
                cif_data = mmCIFSave(tokx[5:])
                cif_file.append(cif_data)
                cif_table_cache = {}
                cif_table = None

                tblx,colx,strx,tokx = token_iter.next()
                

    def gen_token_iter(self, fil):
        re_tok = re.compile(
            r"(?:"

             "(?:_(.+?)[.](\S+))"               "|"  # _section.subsection

             "(?:['\"](.*?)(?:['\"]\s|['\"]$))" "|"  # quoted strings

             "(\S+)"                                 # unquoted tokens

             ")")

        ## get file size for update callbacks
        percent_done = 0
        fil_read_bytes = 0

        ## some file objects do not support seek/tell
        if hasattr(fil, "seek") and hasattr(fil, "tell"):
            try:
                fil.seek(0, 2)
                fil_size_bytes = fil.tell()
                fil.seek(0, 0)
            except:
                # this is a adverage file size ;)
                fil_size_bytes = 1304189
        else:
            fil_size_bytes = 1304189

        ## do we have fast Python, or slow Python?
        try:
            file_iter = iter(fil)
        except TypeError:
            class FileIter(object):
                def __init__(self, fil):
                    self.fil = fil
                def next(self):
                    ln = self.fil.readline()
                    if ln=="":
                        raise StopIteration
                    else:
                        return ln

            file_iter = FileIter(fil)

        while 1:
            try:
                ln = file_iter.next()
            except StopIteration:
                break
            else:
                self.line_number += 1
                fil_read_bytes   += len(ln)

            ## call update callback
            if self.update_cb != None:
                pdone = (fil_read_bytes * 100)/fil_size_bytes
                if pdone != percent_done and pdone <= 100:
                    percent_done = pdone
                    self.update_cb(percent_done)

            ## skip comments
            if ln.startswith("#"):
                continue

            ## semi-colen multi-line strings
            if ln.startswith(";"):
                x = ln[1:]
                while 1:
                    try:
                        ln = file_iter.next()
                    except StopIteration:
                        break
                    else:
                        self.line_number += 1
                    
                    if ln.startswith(";"):
                        break
                    x += ln

                x = x.rstrip()
                yield (None, None, x, None)
                continue

            ## split line into tokens
            tok_iter = re_tok.finditer(ln)
            for tokm in tok_iter:
                yield tokm.groups()


class mmCIFFileWriter(object):
    """Writes out a mmCIF file using the data in the mmCIFData list.
    """  
    def write_file(self, fil, cif_data_list):
        self.fil = fil

        ## constant controlls the spacing between columns
        self.SPACING = 2

        ## iterate through the data sections and write them
        ## out to the file
        for cif_data in cif_data_list:
            self.cif_data = cif_data
            self.write_cif_data()

    def write(self, x):
        self.fil.write(x)

    def writeln(self, x = ""):
        self.fil.write(x+"\n")

    def write_mstring(self, mstring):
        self.write(self.form_mstring(mstring))

    def form_mstring(self, mstring):
        strx = ";"

        lw = MAX_LINE - 2

        for x in mstring.split("\n"):
            if x == "":
                strx += "\n"
                continue
            
            while len(x) > 0:
                x1 = x[:lw]
                x  = x[lw:]
                strx += x1 + "\n"

        strx += ";\n"
        return strx

    def fix(self, x):
        if type(x) != StringType:
            return str(x)
        if x == "":
            return "."
        return x

    def data_type(self, x):
        """Analyze x and return its type: token, qstring, mstring
        """
        if type(x) != StringType:
            x = str(x)
            return x, "token"

        if x == "" or x == ".":
            return ".", "token"

        if x.find("\n") != -1:
            return x, "mstring"
        
        if x.count(" ")!=0 or x.count("\t")!=0:
            if len(x) > MAX_LINE-2:
                return x, "mstring"
            if x.count("' ")!=0:
                return x, "mstring"
            return x, "qstring"

        if len(x) < MAX_LINE:
            return x, "token"
        else:
            return x, "mstring"

    def write_cif_data(self):
        if isinstance(self.cif_data, mmCIFSave):
            self.writeln("save_%s" % self.cif_data.name)
        else:
            self.writeln("data_%s" % self.cif_data.name)
        self.writeln("#")
        
        for cif_table in self.cif_data:
            ## ignore tables without data rows
            if len(cif_table) == 0:
                continue

            ## special handling for tables with one row of data
            elif len(cif_table) == 1:
                self.write_one_row_table(cif_table)

            ## _loop tables
            elif len(cif_table) > 1 and len(cif_table.columns) > 0:
                self.write_multi_row_table(cif_table)

            else:
                print "wtf?",cif_table
                sys.exit(1)

            self.writeln("#")

    def write_one_row_table(self, cif_table):
        row = cif_table[0]

        ## determine max key length for formatting output
        kmax  = 0
        table_len = len(cif_table.name) + 2
        for col in cif_table.columns:
            klen = table_len + len(col)
            assert klen < MAX_LINE
            kmax = max(kmax, klen)

        ## we need a space after the tag
        kmax += self.SPACING
        vmax  = MAX_LINE - kmax - 1

        ## write out the keys and values
        strx = ""
        
        for col in cif_table.columns:
            strx = "_%s.%s" % (cif_table.name, col)
            strx = strx.ljust(kmax)

            try:
                x0 = row[col]
            except KeyError:
                x = "?"
                dtype = "token"
            else:
                x, dtype = self.data_type(x0)

            if dtype == "token":
                if len(x) > vmax:
                    strx += "\n"
                strx += "%s\n" % (x)
                self.write(strx)

            elif dtype == "qstring":
                if len(x) > vmax:
                    strx += "\n"
                strx += "'%s'\n" % (x)
                self.write(strx)

            elif dtype == "mstring":
                strx += "\n"
                self.write(strx)
                self.write_mstring(x)

    def write_multi_row_table(self, cif_table):
        ## write the key description for the _loop
        self.writeln("loop_")
        for col in cif_table.columns:
            key = "_%s.%s" % (cif_table.name, col)
            assert len(key) < MAX_LINE
            self.writeln(key)

        col_len_map   = {}
        col_dtype_map = {}

        for row in cif_table:
            for col in cif_table.columns:
                ## get data and data type
                try:
                    x0 = row[col]
                except KeyError:
                    lenx  = 1
                    dtype = "token"
                else:
                    x, dtype = self.data_type(x0)

                    ## determine write length of data
                    if dtype == "token":
                        lenx = len(x)
                    elif dtype == "qstring":
                        lenx = len(x) + 2
                    else:
                        lenx = 0

                try:
                    col_dtype = col_dtype_map[col]
                except KeyError:
                    col_dtype_map[col] = dtype
                    col_len_map[col]   = lenx
                    continue

                ## update the column charactor width if necessary
                if col_len_map[col] < lenx:
                    col_len_map[col] = lenx

                ## modify column data type if necessary
                if col_dtype != dtype:
                    if dtype == "mstring":
                        col_dtype_map[col] = "mstring"
                    elif col_dtype == "token" and dtype == "qstring":
                        col_dtype_map[col] = "qstring"


        ## form a write list of the column names with values of None to
        ## indicate a newline
        wlist = []
        llen = 0
        for col in cif_table.columns:
            dtype = col_dtype_map[col]

            if dtype == "mstring":
                llen = 0
                wlist.append((None, None, None))
                wlist.append((col, dtype, None))
                continue

            lenx  = col_len_map[col]
            if llen == 0:
                llen = lenx
            else:
                llen += self.SPACING + lenx

            if llen > MAX_LINE-1:
                wlist.append((None, None, None))
                llen = lenx

            wlist.append((col, dtype, lenx))
            
        ## write out the data
        spacing   = " " * self.SPACING
        add_space = False
        listx     = []

        for row in cif_table:
            for (col, dtype, lenx) in wlist:

                if col == None:
                    add_space = False
                    listx.append("\n")
                    continue

                if add_space == True:
                    add_space = False
                    listx.append(spacing)

                if dtype == "token":
                    x = str(row.get(col, "."))
                    if x == "":
                        x = "."
                    x = x.ljust(lenx)
                    listx.append(x)
                    add_space = True
                    
                elif dtype == "qstring":
                    x = row.get(col, ".")
                    if x == "":
                        x = "."
                    elif x != "." and x != "?":
                        x = "'%s'" % (x)
                    x = x.ljust(lenx)
                    listx.append(x)
                    add_space = True

                elif dtype == "mstring":
                    try:
                        listx.append(self.form_mstring(row[col]))
                    except KeyError:
                        listx.append(".\n")
                    add_space = False


            add_space = False
            listx.append("\n")

            ## write out strx if it gets big to avoid using a lot of
            ## memory
            if len(listx) > 1024:
                self.write("".join(listx))
                listx = []

        ## write out the _loop section
        self.write("".join(listx))


### <testing>
if __name__ == '__main__':
    import sys

    try:
        path = sys.argv[1]
    except IndexError:
        print "usage: mmCIF.py <mmCIF file path>"
        sys.exit(1)

    cif = mmCIFDictionary()
    cif.load_file(path)
    cif.save_file(sys.stdout)
### </testing>
