## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""This module provides a primitive mmCIF file parser.  Files are parsed
into a set of data structures where they can be further processed.  The
data structures can also be constructed and written back out as mmCIF.
A CIF dictionary parser is also included as a specialized version of the
mmCIF parser.
"""
from __future__ import generators
import string
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

mmCIFError = "mmCIFError"


class mmCIFRow(dict):
    """Contains one row of data.  In a mmCIF file, this is one complete
    set of data found under a section.  The data can be accessed by using
    the column names as class attributes.
    """
    def __eq__(self, other):
        return id(self) == id(other)

    def __deepcopy__(self, memo):
        return mmCIFRow(self)

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
    
    def __setitem__(self, i, row):
        assert isinstance(row, mmCIFRow)
        row.table = self
        list.__setitem__(self, i, row)

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
        """Iterates through all rows in self, and forms a list of all
        unique column names, then sets the self.columns to that list.
        """
        self.columns = []
        for row in self:
            for col in row.keys():
                if col not in self.columns:
                    self.columns.append(col)
        self.columns.sort()

    def get_row(self, *args):
        """Returns the first row found matching the argument list.
        """
        def chk(row):
            for (key, val) in args:
                if row[key] != val:
                    return False
            return True
        
        for row in self:
            if chk(row):
                return row
        return None

    def iter_rows(self, *args):
        """Iterate over all rows matching the argument list.
        """
        def chk(row):
            for (key, val) in args:
                if row[key] != val:
                    return False
            return True

        for row in self:
            if chk(row):
                yield row

    def select_row_list(self, *nvlist):
        """Preforms a SQL-like 'AND' select aginst all the rows in the table.
        The arguments are a variable list of tuples of the form:
          (<column-name>, <column-value>)
        For example:
          select_row_list(('atom_id','CA'),('entity_id', '1'))
        returns a list of rows with atom_id==1 and entity_id==1.
        """
        ## clever optimization tricks
        (attr, val) = nvlist[0]
        
        row_list = []
        for row in self:
            if row[attr] != val: continue
            
            add = 1
            for (attr, val) in nvlist:
                if row[attr] != val:
                    add = 0
                    break

            if add: row_list.append(row)

        return row_list

    def debug(self):
        print "mmCIFTable::%s" % (self.name)
        for row in self:
            for col in self.columns:
                print "%s=%s" % (col, row[col])[:80]
            print "---"


class mmCIFData(list):
    """Contains all information found under a data_ block in a mmCIF file.
    mmCIF files are represented differently here than their file format
    would suggest.  Since a mmCIF file is more-or-less a SQL database dump,
    the files are represented here with their sections as "Tables" and
    their subsections as "Columns".  The data is stored in "Rows".
    """
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
            for ctable in self:
                if x == ctable.name:
                    return ctable
            raise KeyError, x

        raise TypeError, x

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
            for cdata in self:
                if cdata.name == x:
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
        
    def load_file(self, fil):
        """Load and append the mmCIF data from file object fil into self.
        """
        fil = OpenFile(fil, "r")
        mmCIFFileParser().parse_file(fil, self)

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
    def load_file(self, fil):
        """Load the mmCIF dictionary into self.
        """
        mmCIFDictionaryParser().parse_file(OpenFile(fil, "r"), self)


##
## FILE PARSERS/WRITERS
##

## charactors considered quotes
QUOTES = ["'", '"']

## maximum line width
MAX_LINE = 80


class mmCIFElementFile:
    """Tokenizes a mmCIF file for the state parser.
    """
    def __init__(self, file):
        self.file = file
        self.line_number = 0
        self.elements = []

    def escaped(self, line, i):
        """Given a line and a charactor index, determine if the charactor
        at index i is escaped.
        """
        j = i - 1
        while 1:
            try:
                if line[j] != "\\": break
            except IndexError:
                break
            j -= 1
        return not ((i - j) % 2)

    def split(self, line):
        tok_list = []
        j = 0
        state = "whitespace"

        for i in range(len(line)):
            if state == "whitespace":
                if line[i] not in string.whitespace:
                    j = i
                    if line[i] in QUOTES:
                        state = "quote"
                    else:
                        state = "element"
                continue

            elif state == "element":
                if line[i] in string.whitespace:
                    state = "whitespace"
                    tok = line[j:i]
                    if tok == ".":
                        tok = ""
                    tok_list.append((tok, "token"))
                continue

            elif state == "quote":
                if line[i] in QUOTES:
                    try:
                        if line[i+1] not in string.whitespace: continue
                    except IndexError:
                        pass
                    
                    state = "whitespace"
                    tok_list.append((line[j+1:i],"string")) ## strip quotes
                continue

        if  state == "element":
            tok = line[j:]
            if tok == ".":
                tok = ""
            tok_list.append((tok,"token"))

        elif state == "quote":
            tok_list.append((line[j+1:-1],"string")) ## strip out the quotes
            
        return tok_list
 
    def read_elements(self):
        state = "read data"
        temp = ""
        
        while 1:
            line = self.file.readline()
            self.line_number += 1
            if not line:
                raise EOFError, "mmCIFElementFile"

            if state == "read data":
                if line.startswith("#"):
                    continue
                if line.startswith(";"):
                    state = "read ;string"
                    temp = line[1:]
                    continue
                elements = self.split(line)
                if elements:
                    self.elements += elements
                    break
                continue

            if state == "read ;string":
                if line.startswith(";"):
                    state = "read data"
                    self.elements.append((temp.rstrip(), "string"))
                    break
                temp += line
                continue
        
    def get_next_element(self):
        if not self.elements:
            self.read_elements()
        return self.elements.pop(0)

    def peek_next_element(self):
        if not self.elements:
            self.read_elements()
        return self.elements[0]

    def replace_element(self, element):
        self.elements.insert(0, element)


class mmCIFFileParser:
    """Stateful parser which uses the mmCIFElementFile tokenizer to read
    a mmCIF file and convert it into the mmCIFData/mmCIFTable/mmCIFRow
    data hierarchy.
    """
    error = "mmCIF Syntax Error"
    def syntax_error(self, err):
        err = "[line %d] %s" % (self.cife.line_number, err)
        raise self.error, err

    def parse_file(self, fil, cif_file):
        self.done = 0
        self.cife = mmCIFElementFile(fil)

        while not self.done:
            try:
                (s,t) = self.cife.get_next_element()
            except EOFError:
                self.done = 1
                break

            if s.startswith("data_"):
                data = mmCIFData(s[5:])
                cif_file.append(data)
                self.read_data(data)
            else:
                self.syntax_error('unexpected element="%s"' % (s))

    def read_data(self, data):
        while not self.done:
            try:
                (s,t) = self.cife.get_next_element()
            except EOFError:
                self.done = 1
                break

            if s.startswith("_"):
                self.read_single(data, s)
                
            elif s.startswith("loop_"):
                self.read_loop(data, s)

            elif s.startswith("data_"):
                self.cife.replace_element((s,t))
                break

            else:
                self.syntax_error('bad element=%s' % (s))

    def read_single(self, data, s):
        (tname, cname) = s[1:].split(".")

        table = data.get_table(tname)
        if not table:
            table = mmCIFTable(tname)
            data.append(table)
            table.append(mmCIFRow())
        else:
            if cname in table.columns:
                self.syntax_error('redefined single column %s.%s' % (
                    tname, cname))
    
        table.columns.append(cname)
        row = table[0]

        try:
            (s,t) = self.cife.get_next_element()
        except EOFError:
            self.syntax_error('premature end of file')

        if s.startswith("_") and t == "token":
            self.syntax_error('expected data got section=%s' % (s))

        row[cname] = s

    def read_loop(self, data, s):
        table_name = None
        cname_list = []

        ## read in table column names
        while 1:
            try:
                (s,t) = self.cife.get_next_element()
            except EOFError:
                self.syntax_error('premature end of file')

            ## no more column names, this is data
            if t == "string" or not s.startswith("_"):
                self.cife.replace_element((s,t))
                break

            (tname, cname) = s[1:].split(".")

            if table_name != tname:
                if table_name:
                    self.syntax_error('section in loop do not match')
                table_name = tname

            cname_list.append(cname)

        table = mmCIFTable(table_name, cname_list)
        data.append(table)

        ## read in table column data
        while 1:
            row = mmCIFRow()
            table.append(row)
            
            for cname in table.columns:
                try:
                    (s,t) = self.cife.get_next_element()
                except EOFError:
                    self.syntax_error('loop values incomplete')

                ## freak out
                if s.startswith("_") and t == "token":
                    self.syntax_error('expected data got section=%s' % (s))

                row[cname] = s

            try:
                (s,t) = self.cife.peek_next_element()
            except EOFError:
                self.done = 1
                break

            ## if the next element is a mmCIF control token, then break
            if t == "token":
                if s.startswith("_") or s.startswith("data_") or \
                   s.startswith("loop_") or s.startswith("save_"):
                   break


class mmCIFDictionaryParser(mmCIFFileParser):
    """Subclassed from mmCIFFileParser and extended to support the additional
    syntax encountered in the dictionary files.  I wrote this quite a while
    ago, and now that I look at it again, I suspect it's not complete.
    """
    def parse_file(self, fil, cif_file):
        self.done = 0
        self.cife = mmCIFElementFile(fil)

        while not self.done:
            try:
                (s,t) = self.cife.get_next_element()
            except EOFError:
                self.done = 1
                break

            if s.startswith("data_"):
                data = mmCIFData(s[5:])
                data.save_list = []
            else:
                self.syntax_error('unexpected element="%s"' % (s))

            cif_file.append(data)
            self.read_data(data) 

    def read_data(self, data):
        while not self.done:
            try:
                (s,t) = self.cife.get_next_element()
            except EOFError:
                self.done = 1
                break

            if s.startswith("_"):
                self.read_single(data, s)
                
            elif s.startswith("loop_"):
                self.read_loop(data, s)

            elif s.startswith("data_"):
                self.syntax_error("two data_ subsections in dictionary")
                break

            elif s.startswith("save_"):
                save = mmCIFSave(s[5:])
                data.save_list.append(save)
                self.read_save(save)
                
            else:
                self.syntax_error('bad element=%s' % (s))

    def read_save(self, save):
        while not self.done:
            try:
                (s,t) = self.cife.get_next_element()
            except EOFError:
                self.done = 1
                break

            if s.startswith("_"):
                self.read_single(save, s)
                
            elif s.startswith("loop_"):
                self.read_loop(save, s)

            elif s.startswith("save_"):
                break

            else:
                self.syntax_error('bad element=%s' % (s))



class mmCIFFileWriter:
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

    def write(self, line):
        self.fil.write(line)

    def writeln(self, line = ""):
        self.fil.write(line + "\n")

    def fix_value(self, val):
        ## make sure the value is a string
        val = str(val)

        ## blank values are not allowd -- replace them with
        ## question marks
        if not len(val): return "?"

        ## quote strings with spaces
        if val.find(" ") != -1: return "'%s'" % (val)

        ## don't do anything to the value if we make it here
        return val

    def fix_big_string(self, val):
        ## large strings -- break them up into newlined chunks < size
        if len(val) > MAX_LINE-1:
            tmp = []
            for x in val.split("\n"):
                while len(x) > MAX_LINE-1:
                    tmp.append(x[:MAX_LINE-1])
                    x = x[MAX_LINE-1:]
                tmp.append(x)
            return string.join(tmp, "\n")

        return val

    def write_cif_data(self):
        self.writeln("data_%s" % self.cif_data.name)
        self.writeln("")
        
        for cif_table in self.cif_data:
            if   len(cif_table) == 0: continue
            elif len(cif_table) == 1: self.write_one_row_table(cif_table)
            else:                     self.write_multi_row_table(cif_table)
            self.writeln("#")

    def write_one_row_table(self, cif_table):
        row = cif_table[0]

        kmax  = 0
        for col in cif_table.columns:
            key = "_%s.%s" % (cif_table.name, col)
            kmax = max(kmax, len(key))

        ## we need a space after the tag
        kmax += self.SPACING

        for col in cif_table.columns:
            key = "_%s.%s" % (cif_table.name, col)
            self.write(key.ljust(kmax))

            try:
                val = row[col]
            except KeyError:
                val = "?"

            fval = self.fix_value(val)
            
            if   len(fval) > MAX_LINE-1:
                fval = self.fix_big_string(val)
                self.writeln()
                self.writeln(";" + fval)
                self.writeln(";")

            elif len(fval) > MAX_LINE-kmax-1: 
                self.writeln()
                self.writeln(fval)

            else:
                self.writeln(fval)

    def write_multi_row_table(self, cif_table):
        self.writeln("loop_")
        for col in cif_table.columns:
            self.writeln("_%s.%s" % (cif_table.name, col))

        ## initalize a vmax list, the list of the maximum width
        ## of that column in all the rows
        vmax = [0 for x in cif_table.columns]

        ## optimization
        col_tuple = [(cif_table.columns[i], i) \
                     for i in range(len(cif_table.columns))]

        for row in cif_table:
            for (col, i) in col_tuple:
                try:
                    val = row[col]
                except KeyError:
                    val = "?"
                val = self.fix_value(val)
                vmax[i] = max(vmax[i], len(val))
                
        ## now we know the maximum width of each field, so we can start
        ## writing out rows
        for row in cif_table:
            wlen = MAX_LINE-1
            
            for (col, i) in col_tuple:
                try:
                    val = row[col]
                except KeyError:
                    val = "?"

                ## we don't have enough space left on this line
                if vmax[i]+self.SPACING > wlen and wlen < MAX_LINE-1:
                    self.writeln()
                    wlen = MAX_LINE-1

                ## write out the value
                if vmax[i] > MAX_LINE-1:
                    fval = self.fix_big_string(val)
                    self.writeln(";" + fval)
                    self.writeln(";")

                else:
                    if wlen < MAX_LINE-1:
                        self.write("".ljust(self.SPACING))
                        wlen -= self.SPACING

                    fval = self.fix_value(val)
                    self.write(fval.ljust(vmax[i]))
                    wlen -= vmax[i]

            if wlen < MAX_LINE-1:
                self.writeln()


class mmCIFFileBuilder:
    """Builds a mmCIF file from a Structure object.
    """
    entry_columns = ["id"]

    audit_author_columns = ["name"]
    
    entity_columns = ["id", "type", "details"]

    cell_columns = [
        "entry_id", "length_a", "length_b", "length_c", "angle_alpha",
        "angle_beta", "angle_gamma", "PDB_Z"]

    symmetry_columns = [
        "entry_id", "space_group_name_H-M", "cell_setting",
        "Int_Tables_number"]
    
    atom_site_columns = [
        "group_PDB", "id", "type_symbol", "label_entity_id",
        "Cartn_x", "Cartn_y", "Cartn_z", 
        "occupancy", "B_iso_or_equiv", "Cartn_x_esd", "Cartn_y_esd",
        "Cartn_z_esd", "occupancy_esd", "B_iso_or_equiv_esd", "auth_seq_id",
        "auth_comp_id", "auth_asym_id", "auth_atom_id"]

    atom_site_anisotrop_columns = [
        "id", "type_symbol", "label_entity_id",
        "U[1][1]", "U[1][2]", "U[1][3]", "U[2][2]",
        "U[2][3]", "U[3][3]", "U[1][1]_esd", "U[1][2]_esd", "U[1][3]_esd",
        "U[2][2]_esd", "U[2][3]_esd", "U[3][3]_esd", "pdbx_auth_seq_id",
        "pdbx_auth_comp_id", "pdbx_auth_asym_id", "pdbx_auth_atom_id"]

    copy_cifdb_tables = [
        "entry","audit_author","struct","struct_keywords","struct_conf",
        "struct_siite_gen"]

    def __init__(self, struct, cif_file):
        self.struct = struct
        self.cif_data = mmCIFData("XXX")
        cif_file.append(self.cif_data)

        ## maps fragment -> entity_id
        self.entity_id_map = {}

        ## these tables can be copied directly from the structure's
        ## cif database
        for tbl in self.copy_cifdb_tables:
            try:
                cpy = copy.deepcopy(self.struct.cifdb[tbl])
            except KeyError:
                continue
            self.cif_data.append(cpy)

        ## these tables need to be formed from the atom structure
        self.add_entity()
        self.add_cell()
        self.add_symmetry()
        self.add_atom_site()

    def get_table(self, name, columns = None):
        try:
            return self.cif_data[name]
        except KeyError:
            pass
        
        table = mmCIFTable(name, columns[:])
        self.cif_data.append(table)
        return table

    def add_entity(self):
        ## maps fragment -> entity_id
        entity = self.get_table("entity", self.entity_columns)

        ## I. entity.type == "polymer"
        
        ## first detect polymer chains
        ## map of entity::mmCIFRow -> sequence list
        es_list = []
        
        for chain in self.struct.iter_chains():

            ## if the chain is a bio-polymer, it is one entity; come up
            ## with a name from its sequence and add it to the
            ## entity map
            if not chain.has_standard_residues():
                continue

            sequence = chain.sequence or chain.calc_sequence()

            ## compare aginst previously calculated sequences to
            ## determine the correct entity_id
            entity_id = None

            for (row, seq) in es_list:
                if seq == sequence:
                    entity_id = row["id"]
                    break

            if entity_id == None:
                row = mmCIFRow()
                entity.append(row)
                row["id"] = entity.index(row) + 1
                row["type"] = "polymer"
                if self.struct.library.is_amino_acid(sequence[0]):
                    row["details"] = "%d residue polypeptide"%(len(sequence))
                elif self.struct.library.is_nucleic_acid(sequence[0]):
                    row["details"] = "%d residue DNA/RNA"%(len(sequence))

                entity_id = row["id"]
                es_list.append((row, sequence))

            for res in chain.iter_standard_residues():
                self.entity_id_map[res] = entity_id


        ## II. entity.type == "non-polymer" or "water"
        er_map = {}

        for chain in self.struct.iter_chains():
            for frag in chain.iter_non_standard_residues():

                ## already assigned a entity_id for this fragment_id
                if er_map.has_key(frag.res_name):
                    self.entity_id_map[frag] = er_map[frag.res_name]

                ## we need to assign a entity_id for this fragment_id
                ## and add a row for it in the entity table
                else:
                    row = mmCIFRow()
                    entity.append(row)
                    entity_id = row["id"] = entity.index(row) + 1
                    if frag.is_water():
                        row["type"] = "water"
                        row["details"] = ""
                    else:
                        row["type"] = "non-polymer"
                        row["details"] = frag.res_name

                    er_map[frag.res_name] = entity_id
                    self.entity_id_map[frag] = entity_id

    def add_cell(self):
        if self.struct.unit_cell:
            unit_cell = self.struct.unit_cell
        else:
            return
        
        cell = self.get_table("cell", self.cell_columns)
        row = mmCIFRow()
        cell.append(row)

        row["entry_id"] = self.cif_data["entry"]["id"]
        row["length_a"] = unit_cell.a
        row["length_b"] = unit_cell.b
        row["length_c"] = unit_cell.c
        row["angle_alpha"] = unit_cell.calc_alpha_deg()
        row["angle_beta"] = unit_cell.calc_beta_deg()
        row["angle_gamma"] = unit_cell.calc_gamma_deg()

    def add_symmetry(self):
        if self.struct.unit_cell and self.struct.unit_cell.space_group:
            space_group = self.struct.unit_cell.space_group
        else:
            return

        cell = self.get_table("symmetry", self.symmetry_columns)
        row = mmCIFRow()
        cell.append(row)

        row["entry_id"] = self.cif_data["entry"]["id"]
        row["space_group_name_H-M"] = space_group.pdb_name
        row["Int_Tables_number"] = space_group.number

    def add_atom_site(self):
        atom_site = self.get_table("atom_site", self.atom_site_columns)        

        for atm1 in self.struct.iter_atoms():
            for atm2 in atm1.alt_list:
                asrow = mmCIFRow()
                atom_site.append(asrow)
                self.add_atom_site_row(asrow, atm2)

    def add_atom_site_row(self, asrow, atm):
        if atm.get_fragment().is_standard_residue():
            asrow["group_PDB"] = "ATOM"
        else:
            asrow["group_PDB"] = "HETATM"

        asrow["id"] = atom_site.index(asrow) + 1
        asrow["label_entity_id"] = self.entity_id_map[atm.get_fragment()]
        asrow["auth_atom_id"] = atm.name
        asrow["auth_comp_id"] = atm.res_name
        asrow["auth_seq_id"] = atm.fragment_id
        asrow["auth_asym_id"] = atm.chain_id
        asrow["type_symbol"] = atm.element
        asrow["Cartn_x"] = atm.position[0]
        asrow["Cartn_y"] = atm.position[1]
        asrow["Cartn_z"] = atm.position[2]
        asrow["occupancy"] = atm.occupancy
        asrow["B_iso_or_equiv"] = atm.temp_factor

        if atm.sig_position:
            asrow["Cartn_x_esd"] = atm.sig_position[0]
            asrow["Cartn_y_esd"] = atm.sig_position[1]
            asrow["Cartn_z_esd"] = atm.sig_position[2]
            asrow["occupancy_esd"] = atm.sig_occupancy
            asrow["B_iso_or_equiv_esd"] = atm.sig_temp_factor

        if atm.U:
            aniso = self.get_table("atom_site_anisotrop",
                                   self.atom_site_anisotrop_columns)

            anrow = mmCIFRow()
            aniso.append(anrow)
            anrow["id"] = asrow["id"]
            anrow["type_symbol"] = asrow["type_symbol"]
            anrow["label_entity_id"] = asrow["label_entity_id"]
            anrow["pdbx_auth_seq_id"] = asrow["auth_seq_id"]
            anrow["pdbx_auth_comp_id"] = asrow["auth_comp_id"]
            anrow["pdbx_auth_asym_id"] = asrow["auth_asym_id"]
            anrow["pdbx_auth_atom_id"] = asrow["auth_atom_id"]
            anrow["U[1][1]"] = atm.U[0,0]
            anrow["U[2][2]"] = atm.U[1,1]
            anrow["U[3][3]"] = atm.U[2,2]
            anrow["U[1][2]"] = atm.U[0,1]
            anrow["U[1][3]"] = atm.U[0,2]
            anrow["U[2][3]"] = atm.U[1,2]

            if atm.sig_U:
                anrow["U[1][1]_esd"] = atm.sig_U[0,0]
                anrow["U[2][2]_esd"] = atm.sig_U[1,1]
                anrow["U[3][3]_esd"] = atm.sig_U[2,2]
                anrow["U[1][2]_esd"] = atm.sig_U[0,1]
                anrow["U[1][3]_esd"] = atm.sig_U[0,2]
                anrow["U[2][3]_esd"] = atm.sig_U[1,2]



### <testing>
if __name__ == '__main__':
    import sys

    try:
        path = sys.argv[1]
    except IndexError:
        print "usage: mmCIF.py <mmCIF file path>"
        sys.exit(1)

    cif = mmCIFFile()
    cif.load_file(path)
    cif.save_file(sys.stdout)
### </testing>
