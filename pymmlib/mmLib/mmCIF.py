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
    pass


class mmCIFTable(list):
    """Contains columns and rows of data for a mmCIF section.  Rows of data
    are stored as mmCIFRow classes.
    """
    def __init__(self, name, columns = None):
        list.__init__(self)
        self.name = name
        self.columns = columns or []

    def __setitem__(self, i, row):
        assert isinstance(row, mmCIFRow)
        list.__setitem__(self, i, row)

    def append(self, row):
        assert isinstance(row, mmCIFRow)
        list.append(self, row)

    def insert(self, i, row):
        assert isinstance(row, mmCIFRow)
        list.insert(i, row)

    def autoset_columns(self):
        """Iterates through all rows in self, and forms a list of all
        unique column names, then sets the self.columns to that list.
        """
        self.columns = []
        for row in self:
            for col in row.keys():
                if key not in columns:
                    self.columns.append(key)
        self.columns.sort()

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

    def __getitem__(self, x):
        if type(x) == IntType:
            return list.__getitem__(self, x)

        elif type(x) == StringType:
            for ctable in self:
                if ctable.name == x:
                    return ctable
            raise KeyError, x

        raise TypeError, x

    def __delitem__(self, x):
        list.__delitem__(self, self[x])

    def __setitem__(self, i, table):
        assert isinstance(table, mmCIFTable)
        list.__setitem__(i, table)

    def append(self, table):
        assert isinstance(table, mmCIFTable)
        list.append(self, table)

    def insert(self, i, table):
        assert isinstance(table, mmCIFTable)
        list.insert(self, i, table)

    def has_key(self, x):
        for ctable in self:
            if ctable.name == x:
                return True
        return False

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


class mmCIFFile(list):
    """Class representing a mmCIF files.
    """
    def __getitem__(self, x):
        if type(x) == IntType:
            return list.__getitem__(self, x)

        elif type(x) == StringType:
            for cdata in self:
                if cdata.name == x:
                    return cdata
            raise KeyError, x

        raise TypeError, x

    def __delitem__(self, x):
        list.remove(self, self[x])

    def __setitem__(self, i, cdata):
        assert isinstance(table, mmCIFData)
        list.__setitem__(i, cdata)
        
    def append(self, cdata):
        assert isinstance(cdata, mmCIFData)
        list.append(self, cdata)

    def insert(self, i, cdata):
        assert isinstance(cdata, mmCIFData)
        list.insert(self, i, cdata)

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
    def __init__(self, path_or_fil):
        self._path = path_or_fil
        self._data_list = []

    def load(self):
        if type(self._path) == types.StringType:
            fil = OpenFile(self._path, "r")
        else:
            fil = self._path

        for cif_data in mmCIFDictionaryParser().parse_file(fil):
            self.append(cif_data)

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
        list = []
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
                    list.append((tok, "token"))
                continue

            elif state == "quote":
                if line[i] in QUOTES:
                    try:
                        if line[i+1] not in string.whitespace: continue
                    except IndexError:
                        pass
                    
                    state = "whitespace"
                    list.append((line[j+1:i],"string")) ## strip quotes
                continue

        if  state == "element":
            tok = line[j:]
            if tok == ".":
                tok = ""
            list.append((tok,"token"))

        elif state == "quote":
            list.append((line[j+1:-1],"string")) ## strip out the quotes
            
        return list
 
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
    def parse_file(self, fil):
        self.done = 0
        self.cife = mmCIFElementFile(fil)
        try:
            (s,t) = self.cife.get_next_element()
        except EOFError:
            return None

        if not s.startswith("data_"):
            self.syntax_error('unexpected element=%s' % (s))
            
        return self.read_save_list(dict)

    def read_save_list(self, dict):
        save_list = []
        
        while not self.done:
            try:
                (s,t) = self.cife.get_next_element()
            except EOFError:
                self.done = 1
                break

            if s.startswith("_"):
                self.read_single(data, s)

            elif s.startswith("save_"):
                save = self.read_save(mmCIFData(), s)
                save_list.append(save)
                
            elif s.startswith("loop_"):
                self.read_loop(s)

            elif s.startswith("data_"):
                self.syntax_error('two data_ subsections in dictionary')

            else:
                self.syntax_error('bad element=%s' % (s))

        return save_list

    def read_save(self, data, s):
        if s.startswith("save__"):
            name = s[6:]
        else:
            name = s[5:]

        if not name:
            self.syntax_error('expected opening save_ block')

        self.read_data(data)
        try:
            (s,t) = self.cife.get_next_element()
        except EOFError:
            self.done = 1
        else:
            if s != 'save_':
                self.sytax_error('expected closing save_ block')
        
        return data

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
                self.syntax_error('two data_ subsections in dictionary')
                break

            elif s.startswith("save_"):
                self.cife.replace_element((s,t))
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

    def __init__(self, struct, cif_file):
        self.struct = struct
        self.cif_data = mmCIFData("XXX")
        cif_file.append(self.cif_data)

        ## maps fragment -> entity_id
        self.entity_id_map = {}
        
        self.add_entry()
        self.add_audit_author()
        self.add_entity()
        self.add_cell()
        self.add_atom_site()

    def get_table(self, name, columns = None):
        try:
            return self.cif_data[name]
        except KeyError:
            pass
        
        table = mmCIFTable(name, columns[:])
        self.cif_data.append(table)
        return table

    def add_entry(self):
        entry = self.get_table("entry", self.entry_columns)
        row = mmCIFRow()
        entry.append(row)
        row["id"] = self.struct.exp_data["id"]

    def add_audit_author(self):
        if not self.struct.exp_data.has_key("author_list"):
            return
        
        audit_author = self.get_table("audit_author",
                                      self.audit_author_columns)

        for auth in self.struct.exp_data["author_list"]:
            aurow = mmCIFRow()
            audit_author.append(aurow)
            aurow["name"] = auth

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
        if not self.struct.unit_cell:
            return
        else:
            unit_cell = self.struct.unit_cell
        
        cell = self.get_table("cell", self.cell_columns)
        row = mmCIFRow()
        cell.append(row)

        row["entry_id"] = self.cif_data["entry"][0]["id"]
        row["length_a"] = unit_cell.a
        row["length_b"] = unit_cell.b
        row["length_c"] = unit_cell.c
        row["angle_alpha"] = unit_cell.alpha
        row["angle_beta"] = unit_cell.beta
        row["angle_gamma"] = unit_cell.gamma

    def add_atom_site(self):
        atom_site = self.get_table("atom_site", self.atom_site_columns)
        
        for atm in self.struct.iter_atoms():
            asrow = mmCIFRow()
            atom_site.append(asrow)

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
