#!/usr/bin/env python
## Copyright 2002-2010 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

from __future__ import generators
import sys
import getopt
import copy

from mmLib.mmCIF import *


class mmCIFValidator(object):
    def __init__(self):
        self.dict_list = []
        self.cif_save_cache = {}
        self.parent_tag_cache = {}

    def load_dictionary(self, path):
        """Loads a mmCIF dictionary into the manager
        """
        cif_dict = mmCIFDictionary()
        cif_dict.path = path
        cif_dict.load_file(path)
        self.dict_list.append(cif_dict)

    def iter_cif_saves(self):
        """Iterates all mmCIFSave objects in all dictionaries
        """
        for cif_dict in self.dict_list:
            for cif_save in cif_dict:
                yield cif_save

    def split_tag(self, tag):
        table_name, column = tag[1:].split(".")
        return table_name.lower(), column.lower()

    def join_tag(self, table_name, column):
        return "_%s.%s" % (table_name, column)

    def lookup_cif_save(self, tag):
        """Returns the first save block found in the list of dictionaries.
        """
        try:
            return self.cif_save_cache[tag]
        except KeyError:
            for cif_dict in self.dict_list:
                try:
                    cif_save = cif_dict[tag]
                except KeyError:
                    pass
                else:
                    self.cif_save_cache[tag] = cif_save
                    return cif_save
        return None

    def iter_column_saves(self, table_name):
        """Iterates the mmCIFSave objects of all subsection of the given
        section.
        """
        for cif_save in self.iter_cif_saves():
            try:
                tag = cif_save["item"]["name"]
            except KeyError:
                continue

            table_namex, columnx = self.split_tag(tag)
            if table_name == table_namex:
                yield cif_save

    def iter_mandatory_columns(self, table_name):
        """Iterates the mandatory subsection names of a section.
        """
        for cif_save in self.iter_column_saves(table_name):
            if cif_save["item"]["mandatory_code"] == "yes":
                tag = cif_save["item"]["name"]
                table_namex, columnx = self.split_tag(tag)
                yield columnx

    def lookup_parent(self, tag):
        """Finds the parent tag of the given tag.  Recursive method.
        """
        for cif_save in self.iter_cif_saves():
            try:
                item_linked_table = cif_save["item_linked"]
            except KeyError:
                pass
            else:
                for cif_row in item_linked_table:
                    if cif_row["child_name"] == tag:
                        return self.lookup_parent(cif_row["parent_name"])

        return tag

    def lookup_table_primary_tag(self, table_name):
        """Returns the name of the primary key column for a table, or returns
        None if there is no primary key column.
        """
        cif_save = self.lookup_cif_save(table_name)
        if cif_save == None:
            return None
        try:
            return cif_save["category_key"]["name"]
        except KeyError:
            pass
        return None

    def lookup_children(self, tag):
        """Returns a list of the tag's children.
        """
        child_tag_list = []

        cif_save = self.lookup_cif_save(tag)
        if cif_save != None:
            try:
                item_linked_table = cif_save["item_linked"]
            except KeyError:
                pass
            else:
                for cif_row in item_linked_table:
                    child_tag_list.append(cif_row["child_name"])

        return child_tag_list

    def is_root_tag(self, tag):
        """Returns True if the tag has linked children and no parent tags.
        """
        children = self.lookup_children(tag)
        if len(children) == 0:
            return False
        parent = self.lookup_parent(tag)
        if parent != tag:
            return False
        return True


class Any(object):
    """Returns True when compaired with any other object.
    """
    def __eq__(self, other):
        return True


class mmCIFMerge(object):
    def __init__(self, name, validator):
        self.cif_data = mmCIFData(name) 
        self.validator = validator

        self.merge_accel_dict = {}
        self.merge_accel_list = []

    def log(self, text):
        """Log what is happening.
        """
        sys.stderr.write("[mmCIFMerge] %s\n" % (text))

    def save_file(self, fil):
        cif_file = mmCIFFile()
        cif_file.append(copy.deepcopy(self.cif_data))
        cif_file.save_file(fil)


    ## <general methods>: Utility calls can be made with arguments with
    ##                    any CIF elements

    def get_column_values(self, cif_table, column):
        """Return a list of all the values found in the table under the
        given column.
        """
        value_list = []
        for cif_row in cif_table:
            try:
                value_list.append(cif_row[column])
            except KeyError:
                pass
        return value_list

    def sort_columns(self, cif_table):
        """Sort columns so they look nice ;)
        """
        if "id" in cif_table.columns:
            cif_table.columns.remove("id")
            cif_table.columns.insert(0, "id")

    def remove_blank_values(self, cif_data):
        """Removes any blank [.?] from the cif_data block.
        """
        for cif_table in cif_data:
            for cif_row in cif_table:
                for (column, value) in cif_row.items():
                    if value == "." or value == "?" or value == "":
                        del cif_row[column]

    def add_parent_values(self, cif_data):
        """Checks all linked values in the cif_data block for validity by 
        checking all child values aginst parent values and ensureing the 
        parent values exist. If the parent values do not exist, they are added.
        """
        for cif_table in cif_data:
            for column in cif_table.columns:
                tag        = self.validator.join_tag(cif_table.name, column)
                parent_tag = self.validator.lookup_parent(tag)

                if parent_tag == tag:
                    continue

                ## if we get here, then this table/column has a parent
                (par_tbl, par_col) = self.validator.split_tag(parent_tag)
                try:
                    parent_values = self.get_column_values(
                        cif_data[par_tbl], par_col)
                except KeyError:
                    parent_values = []

                child_values = self.get_column_values(cif_table, column)

                ## fill parent_add_list with values not found in the
                ## parent table
                parent_add_list = []
                for child_value in child_values:
                    if child_value in parent_values:
                        continue
                    if child_value in parent_add_list:
                        continue

                    parent_add_list.append(child_value)

                ## add any unaccounted for parent values
                for value in parent_add_list:
                    try:
                        parent_cif_table = cif_data[par_tbl]
                    except KeyError:
                        parent_cif_table = mmCIFTable(par_tbl)
                        cif_data.append(parent_cif_table)

                    if par_col not in parent_cif_table.columns:
                        parent_cif_table.columns.append(par_col)

                    cif_row = mmCIFRow()
                    parent_cif_table.append(cif_row)
                    cif_row[par_col] = value

                    self.log("adding parent %s=%s for %s=%s" % (
                        parent_tag, value, tag, value))

    def change_root_value(
        self, cif_data, table_name, column, old_value, new_value):
        """Changes the value of the table_name.column from old value to 
        new_value, and propagates that change throughout all linked items
        in the cif_data. If the new_value matches a row in the root table,
        then the row containing old_value is removed instead of having its
        value changed.
        """
        try:
            cif_table = cif_data[table_name]

        except KeyError:
            ## what? no table?  make one
            cif_table = mmCIFTable(table_name, [column])
            cif_data.append(cif_table)

            cif_row = mmCIFRow()
            cif_table.append(cif_row)
            cif_row[column] = new_value

        else:
            old_cif_row = None
            new_cif_row = None

            ## check for a row in the table matching the new value
            for cif_row in cif_table:
                if cif_row.get(column) == new_value:
                    new_cif_row = cif_row
                elif cif_row.get(column) == old_value:
                    old_cif_row = cif_row

            ## if there was a row in the table already matching the new
            ## value, then use it and remove the old row
            if old_cif_row == None and new_cif_row == None:
                cif_row = mmCIFRow()
                cif_table.append(cif_row)
                cif_row[column] = new_value

            elif old_cif_row != None and new_cif_row == None:
                old_cif_row[column] = new_value

            elif old_cif_row == None and new_cif_row != None:
                cif_table.remove(old_cif_row)

            elif old_cif_row != None and new_cif_row != None:
                cif_table.remove(old_cif_row)

        ## set the child values
        num_changed = 0
        root_tag    = self.validator.join_tag(table_name, column)
        child_tags  = self.validator.lookup_children(root_tag)

        for tag in child_tags:
            table_name, column = self.validator.split_tag(tag)

            try:
                cif_table = cif_data[table_name]
            except KeyError:
                continue

            for cif_row in cif_table:
                try:
                    value = cif_row[column]
                except KeyError:
                    pass
                else:
                    if value == old_value:
                        num_changed += 1
                        cif_row[column] = new_value

                        self.log("changed value %s=%s->%s" % (
                            tag, old_value, new_value))

    def change_root_value_uber_alles(
        self, cif_data, table_name, column, old_value, new_value):
        """For the given root table/column, this function keeps only one
        row (the one matching old_value), changes the value to new_value,
        and sets every possible linked table/row in cif_data with new_value.
        """
        ## STEP 1:
        ## first the root table
        try:
            cif_table = cif_data[table_name]

        except KeyError:
            cif_table = mmCIFTable(table_name, [column])
            cif_data.append(cif_table)

            cif_row = mmCIFRow()
            cif_table.append(cif_row)
            cif_row[column] = new_value

        else:
            ## find the cif_row with the old value, and the cif_row with
            ## the new value if they exist
            old_cif_row = None
            new_cif_row = None

            for cif_row in cif_table:
                if cif_row.get(column) == old_value:
                    old_cif_row = cif_row
                if cif_row.get(column) == new_value:
                    new_cif_row = cif_row

            keep_cif_row = None

            if new_cif_row != None:
                keep_cif_row = new_cif_row

            elif old_cif_row != None:
                old_cif_row[column] = new_value
                keep_cif_row = old_cif_row

            else:
                if column not in cif_table.columns:
                    cif_table.columns.append(column)

                keep_cif_row = cif_row = mmCIFRow()
                cif_table.append(cif_row)
                cif_row[column] = new_value

            ## remove all the rows except for the one matching the new_value
            remove_list = []
            for cif_row in cif_table:
                if cif_row != keep_cif_row:
                    remove_list.append(cif_row)

            for cif_row in remove_list:
                cif_table.remove(cif_row)

        ## STEP 2:
        ## Now get their little children too! If the child tables exist
        ## in the cif_data at all, then set the appropriate column of
        ## every row to new_value.
        root_tag = self.validator.join_tag(table_name, column)
        child_tags = self.validator.lookup_children(root_tag)

        for child_tag in child_tags:
            (chld_tbl, chld_col) = self.validator.split_tag(child_tag)

            try:
                cif_table = cif_data[chld_tbl]
            except KeyError:
                continue

            if chld_col not in cif_table.columns:
                cif_table.columns.append(chld_col)

            for cif_row in cif_table:
                cif_row[chld_col] = new_value

    ## </general methods>

    ## <merge utility methods>: These are part of the merge process and
    ##                          work on self.cif_data, and the parts
    ##                          of other CIF files being merged.

    def resolve_cif_data_conflicts(self, cif_data):
        """Compares primary column values of the cif_data with values
        in the merged data and changes them to new, non conflicting
        values.  This is done before the cif_data can be merged.
        """
        for cif_table in cif_data:
            tag = self.validator.lookup_table_primary_tag(cif_table.name)
            if tag == None:
                continue

            table_name, column = self.validator.split_tag(tag)

            ## we are only interested in primary tags which are also
            ## root tags
            if self.validator.is_root_tag(tag) == False:
                continue

            if column not in cif_table.columns:
                continue

            try:
                ctx = self.cif_data[table_name]
            except KeyError:
                continue

            ## dictionary of [primary_key_value] = cif_row
            primary_values = {}
            for crx in ctx:
                try:
                    value = crx[column]
                except KeyError:
                    pass
                else:
                    primary_values[value] = crx

            for cif_row in cif_table:
                try:
                    value = cif_row[column]
                    crx = primary_values[value]
                except KeyError:
                    continue

                self.log("%s: file primary value %s=%s conflicts with "\
                         "merged values" % (cif_data.file.path, tag, value))

                ## check cif_row with crx, if cif_row has no conflicting
                ## column values, then the row will cleanly merge in
                ## and nothing needs to be done
                will_merge = True
                for (col, val) in cif_row.items():
                    if crx.has_key(col) == False:
                        continue
                    if crx[col] != val:
                        will_merge = False
                        break
                if will_merge == True:
                    continue

                ## after a lot of thought, this seems to make sense
                prefix = cif_table.data.file.path.replace(".","_")
                new_value = "%s%s" % (prefix, value)

                self.change_root_value(
                    cif_data, table_name, column, value, new_value)

    def create_unused_primary_key(self, table_name, primary_column):
        """Assume primary keys are numeric, return the lowest unused
        primary key in the cif_table.
        """
        try:
            ctx = self.cif_data[table_name]
        except KeyError:
            return "1"
        else:
            primary_keys = self.get_column_values(ctx, primary_column)

            new_key = 1
            while str(new_key) in primary_keys:
                new_key += 1
            return str(new_key)        

    def make_comare_accelerator(self, cif_row, columns):
        """Returns a tuple of the values in cif_row in the order of columns.
        For rows which do not have data for a column, the special class
        Any() is created which compares True against any value.

        merge_accel: (cif_row[column1], cif_row[column2], ...)
        """
        merge_accel = []

        for column in columns:
            merge_accel.append(cif_row.get(column, Any()))

        ## last item is a special Any() with a reference to the row
        any = Any()
        any.cif_row = cif_row
        merge_accel.append(any)

        return tuple(merge_accel)

    def do_cif_row_merge(self, crx, cif_row, columns):
        """Merge the given columns of cif_row2 into cif_row1.
        """
        log_merged_columns = []

        for column in columns:
            value = crx.get(column)

            if value == None or value == "" or value == "?" or value == ".":
                try:
                    crx[column] = cif_row[column]
                except KeyError:
                    pass
                else:
                    log_merged_columns.append(column)

                    if column not in crx.table.columns:
                        crx.table.columns.append(column)

        ## log what happened, return
        if len(log_merged_columns) > 0:
            i = crx.table.index(crx)
            self.log("%s: merged columns=%s into existing row=%d" % (
                crx.table.name, string.join(log_merged_columns, ","), i))
        else:
            i1 = crx.table.index(crx)
            i2 = cif_row.table.index(cif_row)
            self.log("%s: table %s duplicate row=%d file row=%d" % (
                cif_row.table.data.file.path, crx.table.name, i1, i2))

    def merge_cif_row(self, ctx, cif_row, columns):
        """Adds the cif_row to the cif_table after checking to see if the
        row is a duplicate, or a near duplicate with extra columns. If
        the cif_row matches a cif_row already in the cif_table, only the
        new columns are merged into the existing row. Returns 1 if the
        row was appended to the table, and returns 0 if the row values were
        merged, or if the row was an exact match with an existing row.
        """
        merge_accel = self.make_comare_accelerator(cif_row, columns)

        match_found = False
        for accel in self.merge_accel_list:
            if accel == merge_accel:
                match_found = True
                break

        if match_found == True:
            ## CASE 1
            ## there is a row which matches this row already in the table
            ## merge any extra columns of data this row may have into
            ## the matching row
            crx = accel[-1].cif_row
            self.do_cif_row_merge(crx, cif_row, columns)

            self.merge_accel_list.remove(accel)
            accel = self.make_comare_accelerator(crx, columns)
            self.merge_accel_list.append(accel)

            return 0
        else:
            ## CASE 2:
            ## if we get here, then the row does not match any
            ## other row and should be appended
            crx = copy.deepcopy(cif_row)
            ctx.append(crx)

            accel = self.make_comare_accelerator(crx, columns)
            self.merge_accel_list.append(accel)

            return 1

    def merge_cif_table(self, cif_table, columns):
        """Merge the given columns of the cif_table.
        """
        ## the accelerator must be blanked for each new table
        ## to be merged
        self.merge_accel_list = []

        try:
            ctx = self.cif_data[cif_table.name]

        except KeyError:
            ## create new table if necessary
            ctx = mmCIFTable(cif_table.name, copy.deepcopy(cif_table.columns))
            self.cif_data.append(ctx)
            self.log("%s: adding new table %s" % (
                cif_table.data.file.path, cif_table.name))

        else:
            ## fill the merge accelerator list and dictionary
            ## for the row merge
            for crx in ctx:
                accel = self.make_comare_accelerator(crx, columns)
                self.merge_accel_list.append(accel)

        num = 0
        for cif_row in cif_table:
            num += self.merge_cif_row(ctx, cif_row, columns)

        self.log("%s: added %d rows into table %s" % (
            cif_table.data.file.path, num, cif_table.name))

    def merge_cif_data(self, cif_data):
        """Merge the cif_data block.
        """
        ## Before merging, remove any blank values from the cif_data
        ## block
        self.remove_blank_values(cif_data)

        ## Before merging a cif_data block, it needs have tables
        ## with linked data fully filled in to avoid accidental
        ## data overlaps when merging with other files
        self.add_parent_values(cif_data)

        ## Before merging a cif_data block, it must be scanned for
        ## tables with primary key conflicts with tables already in
        ## the merged files. Any primary keys which confict must be
        ## changed to a new value before the merge.
        self.resolve_cif_data_conflicts(cif_data)

        ## merge the tables
        for cif_table in cif_data:
            columns = copy.deepcopy(cif_table.columns)
            self.merge_cif_table(cif_table, columns)

    ## </merge utility methods>

    ## <API CALLS>

    def merge_cif_file(self, cif_file, data_list = None):
        """Merge the cif_file. If a data_list is given, only merge
        the data sections in the list, otherwise merge all data
        sections in the cif_file.
        """
        for cif_data in cif_file:
            if data_list and cif_data.name not in data_list:
                continue
            self.merge_cif_data(cif_data)

    def post_process(self):
        """Called after all cif_files have been merged, this method
        post-processes the merged file.
        """
        for cif_table in self.cif_data:
            cif_table.autoset_columns()
            self.sort_columns(cif_table)

    ## </API CALLS>


class UsageException(Exception):
    def __init__(self, text):
        self.text = text
    def __str__(self):
        return "Invalid Argument: %s" % (self.text)


def usage():
    print 'cifmerge.py - A utility for intelligent merging and manipulation of'
    print '              mmCIF files by utilizing linked data definitions'
    print '              found in  mmCIF dictionaries.'
    print
    print 'usage: cifmerge.py [-s "table.column=old:new"]'
    print '                   [-u "table.column=old:new"]'
    print '                   [-d mmCIF dictionary_filename]'
    print '                   [-f merge_filename]'
    print '                   [-n data_name]'
    print '                   OUTPUT_CIF_FILENAME'
    print
    print '    -h,-?'
    print '        Prints this help page.'
    print
    print '    -s table.column=old:new'
    print '        Changes the value of a row in table.column from old to'
    print '        new, and then changes all linked child values of that'
    print '        table.column according to the loaded mmCIF dictionaries.'
    print
    print '    -u table.column=old:new'
    print '        Uber-alles value change. This reduces table to a one row'
    print '        table with only the value new. If an old value is given,'
    print '        any other columns of data in the row will be preserved.'
    print '        All linked data items found in the mmCIF dictionaries'
    print '        in the entire file are then changed to value new.'
    print
    print '    -d dictionary_filename'
    print '        Specifies a mmCIF dictionary to use. Dictionaries are'
    print '        necessary for the correct merging of tables with linked'
    print '        values.'
    print
    print '    -f merge_filename'
    print '        A mmCIF file to merge. Use multiple times to merge'
    print '        multiple files.'
    print 
    print '    -n name'
    print '        Give the output mmCIF file data the argument name. If'
    print '        no name is specified, XXXX is used.'
    print

    raise SystemExit


def decode_change_arg(change_arg):
    """Decodes the string:table.column=old:new into:(table, column, old, new).
    """
    ieq = change_arg.find("=")
    if ieq == -1:
        return None

    table_name_column = change_arg[:ieq].strip()
    try:
        table_name, column = table_name_column.split(".")
    except ValueError:
        raise UsageException(change_arg)

    old_new = change_arg[ieq+1:].strip()

    if old_new[0] == "'" or old_new[0] == '"':
        quote = old_new[0]
    else:
        quote = None

    if quote == None:
        try:
            old, new = old_new.split(":")
        except ValueError:
            return UsageException(change_arg)
        old = old.strip()
        new = new.strip()
    else:
        old_new = old_new[1:-1]
        try:
            old, new = old_new.split("%s:%s" % (quote))
        except ValueError:
            raise UsageException(change_arg)

    return table_name, column, old, new


def main():
    ## parse options
    (opts, args) = getopt.getopt(sys.argv[1:], "h?u:s:d:f:n:")

    s_arg_list  = []
    u_arg_list  = []
    d_arg_list  = []
    f_arg_list  = []
    n_arg       = None
    output_file = None

    for (opt, item) in opts:
        if opt == "-h" or opt == "-?":
            usage()

        elif opt == "-s":
            s_arg_list.append(decode_change_arg(item))

        elif opt == "-u":
            u_arg_list.append(decode_change_arg(item))

        elif opt == "-d":
            d_arg_list.append(item)

        elif opt == "-f":
            f_arg_list.append(item)

        elif opt == "-n":
            n_arg = item

    if len(args) != 1:
        raise UsageException("One ouput file path required.")
    output_path = args[0]

    ## create validator and load dictionaries
    validator = mmCIFValidator()

    for path in d_arg_list:
        sys.stderr.write("[LOADING DICTIONARY] %s\n" % (path))
        validator.load_dictionary(path)

    ## load mmCIF files
    cif_list = []
    for path in f_arg_list:
        sys.stderr.write("[LOADING MMCIF] %s\n" % (path))

        cif_file = mmCIFFile()
        cif_file.path = path
        cif_list.append(cif_file)
        cif_file.load_file(path)

    ## create merge object and merge files
    name = n_arg or "XXXX"
    merge = mmCIFMerge(name, validator)
    for cif_file in cif_list:
        merge.merge_cif_file(cif_file)

    ## apply set/uber set
    for (table_name, column, old, new) in u_arg_list:
        merge.change_root_value_uber_alles(
            merge.cif_data, table_name, column, old, new)

    for (table_name, column, old, new) in s_arg_list:
        merge.change_root_value(
            merge.cif_data, table_name, column, old, new)

    ## merge and write
    merge.post_process()
    if output_path != None:
        sys.stderr.write("[WRITING MMCIF] %s\n" % (output_path))
        merge.save_file(output_path)


if __name__ == "__main__":

    try:
        main()
    except UsageException, uerr:
        print "[USAGE ERROR]: %s" % (uerr.text)
        print

        usage()
