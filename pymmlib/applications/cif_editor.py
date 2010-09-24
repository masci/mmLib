#!/usr/bin/env python
## Copyright 2002-2010 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import os
import sys
import math
import string
import copy
import cPickle

import pygtk
pygtk.require("2.0")
import gobject
import gtk

from mmLib import mmCIF


## constants
__version__ = "0.99.0"

CIF_EDITOR_CONF = ".cif_editor"
LARGE_DIALOG_SIZE = (400, 300)

## /constants



## large text edit window

class EditDialog(gtk.Dialog):
    def __init__(self, context, cif_row, col_name):
        self.context = context
        self.cif_row = cif_row
        self.col_name = col_name

        title = "[%s] %s.%s.%s" % (self.context.path,
                                   cif_row.table.data.name,
                                   cif_row.table.name,
                                   col_name)

        gtk.Dialog.__init__(self,
                            title,
                            self.context.mw,
                            gtk.DIALOG_DESTROY_WITH_PARENT)

        self.add_button(gtk.STOCK_OK, gtk.RESPONSE_OK)
        self.add_button(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL)
        self.set_default_size(*LARGE_DIALOG_SIZE)
        self.connect("destroy", self.destroy_cb)
        self.connect("response", self.response_cb)

        ## label
        frame = gtk.Frame()
        self.vbox.add(frame)
        frame.set_label("_data_%s _%s.%s" % (cif_row.table.data.name,
                                             cif_row.table.name,
                                             col_name))
        frame.set_border_width(5)

        ## text window
        sw = gtk.ScrolledWindow()
        frame.add(sw)
        sw.set_border_width(5)
        sw.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        self.text_view = gtk.TextView()
        sw.add(self.text_view)        

        buffer = self.text_view.get_buffer()
        iter = buffer.get_start_iter()
        buffer.insert(iter, cif_row.get(col_name, ""))

        self.show_all()

    def destroy_cb(self, dialog):
        assert dialog == self
        self.context.close_edit_dialog(self)

    def response_cb(self, dialog, code):
        assert dialog == self        
        if code == gtk.RESPONSE_OK:
            buffer = self.text_view.get_buffer()
            new_text = buffer.get_text(
                buffer.get_start_iter(), buffer.get_end_iter(), False)
            self.context.cif_row_set_value(
                self.cif_row, self.col_name, new_text)

        self.destroy()

## /large text edit window


class FileControlEditDialog(gtk.Dialog):
    """Model dialog for editing the names of cif_data, cif_table, and
    (cif_table, col_name) entries.
    """
    def __init__(self, context, cif):
        self.context = context
        self.cif = cif

        if isinstance(self.cif, mmCIF.mmCIFData):
            msg = "Enter mmCIF Data Block Name"
            entry_text = self.cif.name

        elif isinstance(self.cif, mmCIF.mmCIFTable):
            msg = "Enter mmCIF Section Name"
            entry_text = self.cif.name

        elif isinstance(self.cif, tuple):
            (cif_table, col_name) = self.cif
            msg = "Enter mmCIF Subsection Name"
            entry_text = col_name

        gtk.Dialog.__init__(self,
                            msg,
                            self.context.mw,
                            gtk.DIALOG_DESTROY_WITH_PARENT|gtk.DIALOG_MODAL)

        self.add_button(gtk.STOCK_OK, gtk.RESPONSE_OK)
        self.add_button(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL)
        self.connect("response", self.response_cb)

        self.hbox = gtk.HBox()
        self.vbox.pack_start(self.hbox, True, True, 0)
        self.hbox.set_border_width(5)

        label = gtk.Label("Name:")
        self.hbox.pack_start(label, False, False, 0)

        self.entry = gtk.Entry()
        self.hbox.pack_start(self.entry, True, True, 0)
        self.entry.set_text(entry_text)

        self.show_all()

    def response_cb(self, dialog, code):
        assert dialog == self

        if code != gtk.RESPONSE_OK:
            self.destroy()
            return

        name = self.entry.get_text()

        if isinstance(self.cif, mmCIF.mmCIFData):
            self.context.cif_data_set_name(self.cif, name)

        elif isinstance(self.cif, mmCIF.mmCIFTable):
            self.context.cif_table_set_name(self.cif, name)

        elif isinstance(self.cif, tuple):
            (cif_table, col_name) = self.cif
            i = cif_table.columns.index(col_name)
            self.context.cif_table_column_set_name(cif_table, i, name)

        self.destroy()


## <mmCIF EDITOR PANEL>


class TableTreeControl(gtk.TreeView):
    def __init__(self, context):
        self.context = context
        self.fade_list = []
        self.cif_table = None
        self.disable_edited_cb = False

        gtk.TreeView.__init__(self)
        self.set_rules_hint(True)
        self.connect("row-activated", self.row_activated_cb)
        self.connect("button-release-event", self.button_release_event_cb)

    def row_activated_cb(self, tree_view, path, column):
        """Retrieve selected node, then call the correct set method for the
        type.
        """
        assert tree_view == self
        cif_row = self.get_cif_row_from_path(path)
        col_name = column.get_data("col_name")
        self.context.open_edit_dialog(cif_row, col_name)

    def button_release_event_cb(self, tree_view, bevent):
        """
        """
        assert tree_view == self
        x = int(bevent.x)
        y = int(bevent.y)

        retval = self.get_path_at_pos(x, y)
        if retval == None:
            return False

        (path, column, x, y) = retval

        cif_row = self.get_cif_row_from_path(path)
        col_name = column.get_data("col_name")
        APP.set_help_window("_%s.%s" % (cif_row.table.name, col_name))
        return False

    def cell_edited_cb(self, cell, row, new_text):
        """Called when a cell has been edited, and a new value entered.
        """
        if self.disable_edited_cb:
            return

        col_name = cell.get_data("col_name")
        cif_row_index = int(row)
        cif_row = self.cif_table[cif_row_index]

        ## only request the update if the value changed
        if not cif_row.has_key(col_name) or cif_row[col_name] != new_text:
            self.context.cif_row_set_value(cif_row, col_name, new_text)

    def get_row(self, path):
        """Returns the cif_row at the given model path.
        """
        (i, ) = path
        return self.cif_table[i]

    def cif_row_get_path(self, cif_row):
        """Returns a model path for locating the cif_row in the model.
        """
        return (self.cif_table.index(cif_row), )

    def cif_row_get_iter(self, cif_row):
        """Returns a model iterator for changing cif_row's model values.
        """
        path = self.cif_row_get_path(cif_row)
        return self.model.get_iter(path)

    def get_cif_row_from_path(self, path):
        """Returns the cif_row object presented in the given path of the
        model.
        """
        return self.get_row(path)

    def get_selected_cif_row(self):
        """Returns the selected cif_row.
        """
        if self.cif_table == None:
            return None

        selection = self.get_selection()
        if selection == None:
            return None

        (model, iter) =  selection.get_selected()
        assert model == self.model
        if iter == None:
            return None

        path = self.model.get_path(iter)
        return self.get_row(path)

    def scroll_to_cif_row(self, cif_row):
        path = self.cif_row_get_path(cif_row)
        self.scroll_to_cell(path, None, False, 0.0, 0.0)

    def set_cif_row(self, cif_row):
        """Sets the cif_table for this control, if needed, and then scrolls to
        and selects the given cif_row.
        """
        if cif_row.table != self.cif_table:
            self.set_cif_table(cif_row.table)
        self.scroll_to_cif_row(cif_row)

    def set_cif_table(self, cif_table):
        """Resets the cif table displayed in the control.
        """
        self.cif_table = cif_table

        ## remove any existing columns
        for col in self.get_columns():
            self.remove_column(col)

        ## this sets the control blank if there is no model
        if self.cif_table == None:
            self.context.mw.table_frame.set_label("")
            self.set_model(None)
            return

        ## update the decorative frame in the window
        self.context.mw.table_frame.set_label(
            "%s.%s" % (cif_table.data.name, cif_table.name))

        ## does somthing helpful if there are no columns in the set table
        if len(self.cif_table.columns) == 0:
            cell = gtk.CellRendererText()
            col = gtk.TreeViewColumn("No data columns in table...", cell)
            self.append_column(col)
            self.model = gtk.ListStore(gobject.TYPE_STRING)
            self.set_model(self.model)
            return

        ## create the new self.model
        num_cols = len(self.cif_table.columns)
        gtk_cols = [gobject.TYPE_STRING, gobject.TYPE_BOOLEAN] * num_cols

        self.model = gtk.ListStore(*gtk_cols)
        self.set_model(self.model)

        ## add columns based on self.cif_table.columns
        for col_name in self.cif_table.columns:
            cell = gtk.CellRendererText()
            cell.set_data("col_name", col_name)
            cell.connect("edited", self.cell_edited_cb)

            model_index = self.cif_table.columns.index(col_name) * 2

            col = gtk.TreeViewColumn(col_name.replace("_","__"), cell)
            col.set_data("col_name", col_name)
            col.add_attribute(cell, "markup", model_index)
            col.add_attribute(cell, "editable", model_index+1)

            self.append_column(col)

        self.reset_rows()

    def reset_rows(self):
        """Removes all rows from the model and reloads from self.cif_table.
        """
        self.model.clear()
        for cif_row in self.cif_table:
            self.set_cif_row_values(cif_row, self.model.append())

    def fade_cb(self, data_tuple):
        """Callback for the ultra cool highlight then fade feature after
        a cif_row insert or update.
        """
        (cif_row, iter, color) = data_tuple

        call_again = False
        for i in xrange(len(color)):
            if color[i] > 0:
                color[i] = max(0, color[i] - 5)
                call_again = True

        if call_again == False:
            self.set_cif_row_values(cif_row, iter)
        else:
            self.set_cif_row_values(cif_row, iter, color)

        return call_again

    def set_cif_row_values(self, cif_row, iter = None, color = None):
        """Called to update a row of data in the model.  If no iterator is
        given, then use the position of the cif_row in the self.cif_table
        to determine the row index in the model and obtain the iterator.
        """
        if iter == None:
            iter = self.cif_row_get_iter(cif_row)

        row_list = []

        for i in xrange(len(self.cif_table.columns)):
            col_name = self.cif_table.columns[i]

            try:
                text = cif_row[col_name]
            except KeyError:
                text = "."
            else:
                if text == None:
                    text = "."

            model_index = i * 2
            editible = True

            if len(text) > 40:
                text = text[:40] + '<span foreground="red">[more]</span>'
                editible = False
            if "\n" in text:
                text = text.replace("\n", '<span foreground="blue">|</span>')
                editible = False
            if color != None:
                text = '<span foreground="#%02x%02x%02x">%s</span>' % (
                    color[0],color[1],color[2], text)

            row_list += [model_index, text, model_index + 1, editible]

        self.model.set(iter, *row_list)        

    def insert_cif_row(self, cif_row):
        """Insert a row of data into the model.  The cif_row should already
        be inserted into the cif_table.
        """
        if cif_row.table == self.cif_table:
            pos = self.cif_table.index(cif_row)
            iter = self.model.insert(pos)
            self.set_cif_row_values(cif_row, iter)
            self.scroll_to_cif_row(cif_row)
        else:
            self.set_cif_row(cif_row)

    def remove_cif_row(self, cif_row):
        """Removes one row frm the model.  The cif_row must still be in
        cif_table.
        """
        if cif_row.table == self.cif_table:
            self.scroll_to_cif_row(cif_row)
        else:
            self.set_cif_row(cif_row)

        iter = self.cif_row_get_iter(cif_row)

        self.disable_edited_cb = True
        iter = self.model.remove(iter)
        self.disable_edited_cb = False

    def update_cif_row(self, cif_row, fade = False):
        """Updates the values of the cif_row in the model.
        """
        if cif_row.table == self.cif_table:
            self.scroll_to_cif_row(cif_row)
        else:
            self.set_cif_row(cif_row)

        iter = self.cif_row_get_iter(cif_row)
        if fade:
            colors = [0x00, 0xff, 0x00]
            gtk.timeout_add(50, self.fade_cb, (cif_row, iter, colors))
        else:
            colors = None
        self.set_cif_row_values(cif_row, iter, colors)

    def do_insert_at_selected(self, ins_cif_row = None, before = False):
        cif_row = self.get_selected_cif_row()

        ## if nothing is selected, just append
        if cif_row == None:
            self.do_append(ins_cif_row)
            return

        i = cif_row.table.index(cif_row)

        if not before:
            i += 1

        if ins_cif_row == None:
            ins_cif_row = mmCIF.mmCIFRow()
        else:
            assert isinstance(ins_cif_row, mmCIF.mmCIFRow)

        self.context.cif_table_insert_row(
            cif_row.table, i, ins_cif_row)

    def do_append(self, ins_cif_row = None):
        if ins_cif_row == None:
            ins_cif_row = mmCIF.mmCIFRow()
        else:
            assert isinstance(ins_cif_row, mmCIF.mmCIFRow)

        self.context.cif_table_insert_row(
            self.cif_table, -1, ins_cif_row)

    def do_delete_selected(self):
        cif_row = self.get_selected_cif_row()
        if cif_row == None:
            return
        self.context.cif_row_remove(cif_row)



class FileTreeControl(gtk.TreeView):
    def __init__(self, context):
        self.context = context

        gtk.TreeView.__init__(self)
        self.get_selection().set_mode(gtk.SELECTION_BROWSE)

        self.connect("row-activated", self.row_activated_cb)
        self.connect("button-release-event", self.button_release_event_cb)

        self.model = gtk.TreeStore(
            gobject.TYPE_STRING,    # 0: data name
            gobject.TYPE_BOOLEAN)   # 1: editable
        self.set_model(self.model)

        cell = gtk.CellRendererText()
        column = gtk.TreeViewColumn("Data/Section/Subsection", cell)
        column.add_attribute(cell, "text", 0)
        column.add_attribute(cell, "editable", 1)
        self.append_column(column)

    def set_cif_file(self):
        """Sets the cif file from the current context.  After this is called
        this control assumes it is kept up to dat by the various update
        functions.
        """
        self.model.clear()
        for cif_data in self.context.cif_file:
            iter1 = self.model.append(None)
            self.add_model_cif_data(iter1, cif_data)

    def add_model_cif_data(self, iter1, cif_data):
        """Adds cif_data and children to the model starting at iter1.
        """
        self.model.set(iter1, 0, cif_data.name, 1, False)
        for cif_table in cif_data:
            iter2 = self.model.append(iter1)
            self.add_model_cif_table(iter2, cif_table)

    def add_model_cif_table(self, iter2, cif_table):
        """Adds cif_table and children to the model starting at iter2.
        """
        self.model.set(iter2, 0, cif_table.name, 1, False)
        for col_name in cif_table.columns:
            iter3 = self.model.append(iter2)
            self.model.set(iter3, 0, col_name, 1, False)

    def get_cif(self, path):
        """Returns the cif object selected from the given gtk.Model style path.
        """
        assert len(path) > 0
        assert len(path) <= 3

        if len(path) == 1:
            (i, ) = path
            return self.context.cif_file[i]

        elif len(path) == 2:
            (i, j) = path
            return self.context.cif_file[i][j]

        elif len(path) == 3:
            (i, j, k) = path
            cif_table = self.context.cif_file[i][j]
            return (cif_table, cif_table.columns[k])

    def get_cif_path(self, cif):
        """Returns the model path for the cif instance: cif_data, cif_table,
        or (cif_table, col_name).
        """
        if isinstance(cif, mmCIF.mmCIFData):
            i = cif.file.index(cif)
            return (i, )

        elif isinstance(cif, mmCIF.mmCIFTable):
            i = cif.data.file.index(cif.data)
            j = cif.data.index(cif)
            return (i, j)

        elif isinstance(cif, tuple):
            (cif_table, col_name) = cif
            i = cif_table.data.file.index(cif_table.data)
            j = cif_table.data.index(cif_table)
            k = cif_table.columns.index(col_name)
            return (i, j, k)

    def get_cif_iter(self, cif):
        path = self.get_cif_path(cif)
        return self.model.get_iter(path)

    def get_selected_cif(self):
        """Retrieves the selected cif_data, cif_table, or (cif_table, col_name)
        instance.  If nothing is selected, returns None.
        """
        sel = self.get_selection()
        if sel == None:
            return None

        (model, iter) = sel.get_selected()
        if iter == None:
            return None

        assert model == self.model
        path = self.model.get_path(iter)
        return self.get_cif(path)

    def row_activated_cb(self, tree_view, path, column):
        """Retrieve selected node, then call the correct set method for the
        type.
        """
        assert tree_view == self

        cif = self.get_cif(path)

        if isinstance(cif, mmCIF.mmCIFData):
            pass

        elif isinstance(cif, mmCIF.mmCIFTable):
            self.context.mw.table_ctrl.set_cif_table(cif)

        elif isinstance(cif, tuple):
            (cif_table, col_name) = cif
            self.context.mw.table_ctrl.set_cif_table(cif_table)

    def button_release_event_cb(self, tree_view, bevent):
        assert tree_view == self
        x = int(bevent.x)
        y = int(bevent.y)

        retval = self.get_path_at_pos(x, y)
        if retval == None:
            return False

        (path, column, x, y) = retval

        cif = self.get_cif(path)

        if isinstance(cif, mmCIF.mmCIFTable):
            APP.set_help_window(cif.name.upper())

        elif isinstance(cif, tuple):
            (cif_table, col_name) = cif
            APP.set_help_window("_%s.%s" % (cif_table.name, col_name))

        return False

    def select_cif(self, cif):
        """Selects one of the cif items in the tree view.  Will scroll if
        necessary.
        """
        path = self.get_cif_path(cif_data)
        self.scroll_to_cell(path, None, False, 0.0, 0.0)

    def update_table_ctrl(self, cif, col_name = None):
        """This is a little hokey and inefficent, oh well.
        """
        do_table_reset = False

        cur_table = self.context.mw.table_ctrl.cif_table
        if cur_table:
            if cif == cur_table.data or cif == cur_table:
                do_table_reset = True

        ## update table_ctrl view
        if do_table_reset:
            self.context.mw.table_ctrl.set_cif_table(cur_table)

    def cif_data_changed(self, cif_data):
        """Update because the name of cif_data changed
        """
        iter = self.get_cif_iter(cif_data)
        self.model.set(iter, 0, cif_data.name)
        self.update_table_ctrl(cif_data)

    def cif_data_insert(self, cif_data):
        """Update because cif_data was inserted.
        """
        i = cif_data.file.index(cif_data)

        if i == 0:
            iter = self.model.prepend(None)
        elif i > 0:
            path = (i-1, )
            sibling = self.model.get_iter(path)
            iter = self.model.insert_after(None, sibling)

        self.add_model_cif_data(iter, cif_data)
        self.update_table_ctrl(cif_data)

    def cif_data_remove(self, cif_data):
        """Update because cif_data is about to be removed.
        """
        iter = self.get_cif_iter(cif_data)
        self.model.remove(iter)

    def cif_table_changed(self, cif_table):
        """Update because the name of cif_table changed.
        """
        iter = self.get_cif_iter(cif_table)
        self.model.set(iter, 0, cif_table.name)
        self.update_table_ctrl(cif_table)

    def cif_table_insert(self, cif_table):
        """Update because cif_table was inserted.
        """
        path = self.get_cif_path(cif_table)
        (i, j) = path

        if j == 0:
            parent = self.model.get_iter((i, ))
            iter = self.model.prepend(parent)
        elif j > 0:
            path = (i, j-1)
            sibling = self.model.get_iter(path)
            iter = self.model.insert_after(None, sibling)

        self.add_model_cif_table(iter, cif_table)
        self.update_table_ctrl(cif_table)

    def cif_table_remove(self, cif_table):
        """Update because cif_table is about to be removed.
        """
        iter = self.get_cif_iter(cif_table)
        self.model.remove(iter)

    def cif_table_column_changed(self, cif_table, col_name):
        """Update because the name of a column in cif_table changed from
        old_col_name to new_col_name.
        """
        iter = self.get_cif_iter((cif_table, col_name))
        self.model.set(iter, 0, col_name)
        self.update_table_ctrl(cif_table, col_name)

    def cif_table_column_insert(self, cif_table, col_name):
        """Update because column col_name was inserted in cif_table.
        """
        path = self.get_cif_path((cif_table, col_name))
        (i, j, k) = path

        if k == 0:
            parent = self.model.get_iter((i, j))
            iter = self.model.prepend(parent)
        elif k > 0:
            sibling = self.model.get_iter((i, j, k-1))
            iter = self.model.insert_after(None, sibling)

        self.model.set(iter, 0, col_name)
        self.update_table_ctrl(cif_table, col_name)

    def cif_table_column_remove(self, cif_table, col_name):
        """Update because column col_name is about to be removed from
        cif_table.
        """
        iter = self.get_cif_iter((cif_table, col_name))
        self.model.remove(iter)

    def new_cif_data(self):
        """Create a new cif_data->cif_table->column->cif_row
        """
        cif_data = mmCIF.mmCIFData("New_Data")
        cif_data.append(self.new_cif_table())
        return cif_data

    def new_cif_table(self):
        """Create a new cif_table->column->cif_row
        """
        cif_table = mmCIF.mmCIFTable("New_Section")
        cif_table.columns.append("New_Subsection")
        cif_table.append(mmCIF.mmCIFRow())
        return cif_table

    def do_insert_sibling_at_selected(
        self, ins_cif = None, before = False, sel_cif = None,
        col_val_list = None):
        """Insert a sibling item before or after the selected item in the
        control.
        """
        if sel_cif == None:
            cif = self.get_selected_cif()
        else:
            cif = sel_cif

        ## if nothing is selected, perform check if the cif_file is empty
        ## if it is, then at least add a new mmCIFData
        if cif == None:
            if len(self.context.cif_file) == 0:
                if ins_cif and isinstance(ins_cif, mmCIF.mmCIFData):
                    self.context.cif_file_insert_data(
                        self.context.cif_file, -1, ins_cif)
                elif ins_cif == None:
                    self.context.cif_file_insert_data(
                        self.context.cif_file, -1, self.new_cif_data()) 

        elif isinstance(cif, mmCIF.mmCIFData):
            i = cif.file.index(cif)
            if ins_cif == None:
                ins_cif = self.new_cif_data()
            else:
                assert isinstance(ins_cif, mmCIF.mmCIFData)
            if not before:
                i += 1
            self.context.cif_file_insert_data(cif.file, i, ins_cif)

        elif isinstance(cif, mmCIF.mmCIFTable):
            i = cif.data.index(cif)
            if ins_cif == None:
                ins_cif = self.new_cif_table()
            else:
                assert isinstance(ins_cif, mmCIF.mmCIFTable)
            if not before:
                i += 1
            self.context.cif_data_insert_table(cif.data, i, ins_cif)

        elif isinstance(cif, tuple):
            (cif_table, col_name) = cif
            i = cif_table.columns.index(col_name)            
            if ins_cif == None:
                ins_cif = "New_Subsection"
            if not before:
                i += 1
            self.context.cif_table_column_insert(
                cif_table, i, ins_cif, col_val_list)

    def do_append_child_at_selected(
        self, ins_cif = None, sel_cif = None, col_val_list = None):
        """Insert a child item before or after the selected item in the
        control.
        """
        if sel_cif == None:
            cif = self.get_selected_cif()
        else:
            cif = sel_cif

        if cif == None:
            if len(self.context.cif_file) == 0:
                if ins_cif and isinstance(ins_cif, mmCIF.mmCIFData):
                    self.context.cif_file_insert_data(
                        self.context.cif_file, -1, ins_cif)
                elif ins_cif == None:
                    self.context.cif_file_insert_data(
                        self.context.cif_file, -1, self.new_cif_data()) 

        elif isinstance(cif, mmCIF.mmCIFData):
            if ins_cif == None:
                ins_cif = self.new_cif_table()
            else:
                assert isinstance(ins_cif, mmCIF.mmCIFTable)
            self.context.cif_data_insert_table(cif, -1, ins_cif)

        elif isinstance(cif, mmCIF.mmCIFTable):
            if ins_cif == None:
                self.context.cif_table_column_insert(cif, -1, "New_Subsection")
            else:
                (col_name, col_val_list) = ins_cif
                self.context.cif_table_column_insert(
                    cif, -1, col_name, col_val_list)

    def do_delete_selected(self, sel_cif = None):
        """Delete the currently selected item in the control.
        """
        if sel_cif == None:
            cif = self.get_selected_cif()
            if cif == None:
                return
        else:
            cif = sel_cif

        if isinstance(cif, mmmCIF.mCIFData):
            self.context.cif_data_remove(cif)

        elif isinstance(cif, mmCIF.mmCIFTable):
            self.context.cif_table_remove(cif)

        elif isinstance(cif, tuple):
            (cif_table, col_name) = cif
            self.context.cif_table_column_remove(cif_table, col_name)

    def do_copy_selected(self, cut = False):
        """Copy the currently selected item in the copy buffer,
        and optionally delete.
        """
        cif = self.get_selected_cif()
        if cif == None:
            return

        cif = copy.deepcopy(cif)

        if isinstance(cif, mmCIF.mmCIFData):
            APP.buffer = cif
            if cut:
                self.context.cif_data_remove(cif)

        elif isinstance(cif, mmCIF.mmCIFTable):
            APP.buffer = cif
            if cut:
                self.context.cif_table_remove(cif)

        elif isinstance(cif, tuple):
            (cif_table, col_name) = cif
            col_val_list = []
            for cif_row in cif_table:
                col_val_list.append(cif_row[col_name])
            APP.buffer = (col_name, col_val_list)
            if cut:
                self.context.cif_table_column_remove(cif_table, col_name)

    def do_paste_selected(self, before = False):
        """Insert the current application paste buffer.
        """
        ins_cif = copy.deepcopy(APP.buffer)
        if ins_cif == None:
            return

        sel_cif = self.get_selected_cif()

        if sel_cif == None:
            if isinstance(ins_cif, mmCIF.mmCIFData):
                self.do_append_child_at_selected(ins_cif)

        elif isinstance(sel_cif, mmCIF.mmCIFData):
            if isinstance(ins_cif, mmCIF.mmCIFData):
                self.do_insert_sibling_at_selected(ins_cif, before = before)
            elif isinstance(ins_cif, mmCIF.mmCIFTable):
                self.do_append_child_at_selected(ins_cif)

        elif isinstance(sel_cif, mmCIF.mmCIFTable):
            if isinstance(ins_cif, mmCIF.mmCIFData):
                self.do_insert_sibling_at_selected(
                    ins_cif, before = before, sel_cif = sel_cif.data)
            elif isinstance(ins_cif, mmCIF.mmCIFTable):
                self.do_insert_sibling_at_selected(ins_cif, before = before)
            elif isinstance(ins_cif, tuple):
                self.do_append_child_at_selected(ins_cif)

        elif isinstance(sel_cif, tuple):
            (cif_table, col_name) = sel_cif

            if isinstance(ins_cif, mmCIF.mmCIFData):
                 self.do_insert_sibling_at_selected(
                    ins_cif, before = before, sel_cif = cif_table.data)
            elif isinstance(ins_cif, mmCIF.mmCIFTable):
                self.do_insert_sibling_at_selected(
                    ins_cif, before = before, sel_cif = cif_table)
            elif isinstance(ins_cif, tuple):
                (col_name, col_val_list) = ins_cif
                self.do_insert_sibling_at_selected(
                    ins_cif = col_name,
                    col_val_list = col_val_list,
                    before = before)

    def do_change_name_selected(self):
        """Change the name of the currently selected item in the control.
        """
        cif = self.get_selected_cif()
        if cif == None:
            return
        FileControlEditDialog(self.context, cif)


class mmCIFEditorMainWindow(gtk.Window):
    """Main window class used as the center of action for a
    mmCIFEditorWindowContext.
    """
    def __init__(self, context):
        self.context = context
        gtk.Window.__init__(self)

        ## Create the toplevel window
        self.set_default_size(600, 400)
        self.connect("delete-event", self.quit_cb, self)

        table = gtk.Table(1, 3, False)
        self.add(table)

        ## Create the menubar
        MENU_ITEMS = (
            ('/_File', None, None, 0, '<Branch>' ),
            ('/File/_New', '', self.new_cb, 0,
             '<StockItem>', gtk.STOCK_NEW),
            ('/File/_Open', '', self.open_cb, 0,
             '<StockItem>', gtk.STOCK_OPEN),
            ('/File/_Save', '<control>S', self.save_cb, 0,
             '<StockItem>', gtk.STOCK_SAVE),
            ('/File/Save _As...', None, self.save_as_cb, 0,
             '<StockItem>', gtk.STOCK_SAVE),
            ('/File/sep1', None, None, 0, '<Separator>'),
            ('/File/CIF Dictionaries...', None,
             APP.open_dictionary_manager_window, 0, ''),
            ('/File/sep1', None, None, 0, '<Separator>'),
            ('/File/_Quit', '<control>Q', self.quit_cb, 0,
             '<StockItem>', gtk.STOCK_QUIT),

            ('/Edit/Undo', '<control>Z', self.undo_cb, 0,
             '<StockItem>', gtk.STOCK_UNDO),
            ('/Edit/sep1', None, None, 0, '<Separator>'),
            ('/Edit/Cut', '<control>X', self.cut_cb, 0,
             '<StockItem>', gtk.STOCK_CUT),
            ('/Edit/Copy', '<control>C', self.copy_cb, 0,
             '<StockItem>', gtk.STOCK_COPY),
            ('/Edit/Paste', '<control>V', self.paste_cb, 0,
             '<StockItem>', gtk.STOCK_PASTE),

            ('/Structure/Insert Before', '<control>u',
             self.insert_cif_before_cb, 0, '<StockItem>', gtk.STOCK_ADD),
            ('/Structure/Insert After', '<control>i',
             self.insert_cif_after_cb, 0, '<StockItem>', gtk.STOCK_ADD),
            ('/Structure/Append Child', '<control>y',
             self.append_cif_child_cb, 0, '<StockItem>', gtk.STOCK_ADD),
            ('/Structure/Delete', '<control>p',
             self.delete_cif_cb, 0, '<StockItem>', gtk.STOCK_REMOVE),
            ('/Structure/Change Item Name', '<control>o',
             self.change_cif_cb, 0, '<StockItem>', gtk.STOCK_REFRESH),

            ('/Data/Insert Row Before', '<control>h',
             self.insert_before_cif_row_cb, 0, '<StockItem>', gtk.STOCK_ADD),
            ('/Data/Insert Row After', '<control>j',
             self.insert_after_cif_row_cb, 0, '<StockItem>', gtk.STOCK_ADD),
            ('/Data/Delete', '<control>l',
             self.delete_cif_row_cb, 0, '<StockItem>', gtk.STOCK_REMOVE),

            ('/_Help', None, None, 0, '<Branch>'),
            ('/Help/_About', None, APP.open_about_window, 0, ''),
            ('/Help/Dictionary Help Browser', None,
             APP.open_help_window, 0, ''))

        self.accel_group = gtk.AccelGroup()
        self.add_accel_group(self.accel_group)

        self.item_factory = gtk.ItemFactory(
            gtk.MenuBar, '<main>', self.accel_group)
        self.item_factory.create_items(MENU_ITEMS, self)

        table.attach(self.item_factory.get_widget('<main>'),
                     # X direction              Y direction
                     0, 1,                      0, 1,
                     gtk.EXPAND | gtk.FILL,     0,
                     0,                         0)

        ## CIF display pandel widget
        self.cif_panel_hpaned = gtk.HPaned()
        table.attach(self.cif_panel_hpaned,
                     # X direction           Y direction
                     0, 1,                   1, 2,
                     gtk.EXPAND | gtk.FILL,  gtk.EXPAND | gtk.FILL,
                     0,                      0)

        ## Create statusbar 
        self.hbox = gtk.HBox()
        table.attach(self.hbox,
                     # X direction           Y direction
                     0, 1,                   2, 3,
                     gtk.EXPAND | gtk.FILL,  0,
                     0,                      0)

        self.statusbar = gtk.Statusbar()
        self.hbox.pack_start(self.statusbar, True, True, 0)
        self.statusbar.set_has_resize_grip(False)

        self.pg_frame = gtk.Frame()
        self.hbox.pack_start(self.pg_frame, False, False, 0)
        self.pg_frame.set_border_width(3)

        self.progress = gtk.ProgressBar()
        self.pg_frame.add(self.progress)
        self.progress.set_size_request(100, 0)

        ## the two editor controls: FileTreeControl, TableTreeControl
        self.file_ctrl = FileTreeControl(self.context)
        self.table_ctrl = TableTreeControl(self.context)

        self.sw1 = gtk.ScrolledWindow()
        self.cif_panel_hpaned.add1(self.sw1)
        self.sw1.set_border_width(3)
        self.sw1.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        self.sw1.add(self.file_ctrl)

        self.table_frame = gtk.Frame()
        self.cif_panel_hpaned.add2(self.table_frame)
        self.table_frame.set_border_width(3)

        self.sw2 = gtk.ScrolledWindow()
        self.table_frame.add(self.sw2)
        self.sw2.set_border_width(3)
        self.sw2.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        self.sw2.add(self.table_ctrl)

        self.show_all()

    def new_cb(self, *args):
        """/File/New callback. Opens a new edit window.
        """
        self.context.open_file(None)

    def open_cb(self, *args):
        """/File/Open callback. Open a FileSeletor window for
        opening new files.
        """
        if not hasattr(self, "open_file_selector"):
            self.open_file_selector = gtk.FileSelection(
                "Select mmCIF file (.cif)");
            self.open_file_selector.connect(
                "delete-event", self.open_destroy_cb)
            self.open_file_selector.ok_button.connect(
                "clicked", self.open_ok_cb)
            self.open_file_selector.cancel_button.connect(
                "clicked", self.open_cancel_cb)
            self.open_file_selector.show()
        else:
            self.open_file_selector.present()

    def open_destroy_cb(self, *args):
        self.open_file_selector.hide()
        return True

    def open_ok_cb(self, ok_button):
        """Called by the [OK] button on the FileSelector.
        """
        path = self.open_file_selector.get_filename()
        self.open_file_selector.hide()
        self.context.open_file(path)

    def open_cancel_cb(self, cancel_button):
        """Called by the [Cancel] button on the FileSelector.
        """
        self.open_file_selector.hide()

    def save_cb(self, *args):
        """/File/Save callback.
        """
        self.context.save_file()

    def save_as_cb(self, *args):
        """/File/Save As callback. Open a FileSeletor window for Save As...
        """
        if not hasattr(self, "save_as_file_selector"):
            self.save_as_file_selector = gtk.FileSelection(
                "Save mmCIF file as...");
            self.save_as_file_selector.connect(
                "delete-event", self.save_as_delete_cb)
            self.save_as_file_selector.ok_button.connect(
                "clicked", self.save_as_ok_cb)
            self.save_as_file_selector.cancel_button.connect(
                "clicked", self.save_as_cancel_cb)
            self.save_as_file_selector.show()
        else:
            self.save_as_file_selector.present()

    def save_as_delete_cb(self, *args):
        self.save_as_file_selector.hide()
        return True

    def save_as_ok_cb(self, ok_button):
        """Called by the [OK] button on the FileSelector.
        """
        path = self.save_as_file_selector.get_filename()
        self.save_as_file_selector.hide()
        self.context.save_file(path)

    def save_as_cancel_cb(self, cancel_button):
        """Called by the [Cancel] button on the FileSelector.
        """
        self.save_as_file_selector.hide()

    def quit_cb(self, *args):
        """/File/Quit or the window's delete-event callback.
        """
        self.destroy()
        self.context.quit()
        return True

    def undo_cb(self, *args):
        """/Edit/Undo callback.
        """
        self.context.undo()

    def cut_cb(self, *args):
        """/Edit/Cut callback.
        """
        self.file_ctrl.do_copy_selected(cut = True)

    def copy_cb(self, *args):
        """/Edit/Copy callback.
        """
        self.file_ctrl.do_copy_selected()

    def paste_cb(self, *args):
        """/Edit/Paste callback.
        """
        self.file_ctrl.do_paste_selected()

    def insert_cif_before_cb(self, *args):
        self.file_ctrl.do_insert_sibling_at_selected(before = True)

    def insert_cif_after_cb(self, *args):
        self.file_ctrl.do_insert_sibling_at_selected(before = False)

    def append_cif_child_cb(self, *args):
        self.file_ctrl.do_append_child_at_selected()

    def delete_cif_cb(self, *args):
        self.file_ctrl.do_delete_selected()

    def change_cif_cb(self, *args):
        self.file_ctrl.do_change_name_selected()

    def insert_before_cif_row_cb(self, *args):
        self.table_ctrl.do_insert_at_selected(before = True)

    def insert_after_cif_row_cb(self, *args):
        self.table_ctrl.do_insert_at_selected(before = False)

    def delete_cif_row_cb(self, *args):
        self.table_ctrl.do_delete_selected()

    def update_menuitem_enable_status(self):
        """Sets enabled status of all menuitems.
        """
        mi_save = self.item_factory.get_item('/File/Save')
        mi_save_as = self.item_factory.get_item('/File/Save As...')
        if self.context.saved:
            mi_save.set_sensitive(False)
            mi_save_as.set_sensitive(False)
        else:
            mi_save_as.set_sensitive(True)
            if self.context.path != None:
                mi_save.set_sensitive(True)

        ## turn on undo if there is an undo stack
        mi = self.item_factory.get_item('/Edit/Undo')
        if self.context.undo_list:
            mi.set_sensitive(True)
        else:
            mi.set_sensitive(False)

    def set_status_bar(self, text):
        """Sets the text on the windows statusbar.
        """
        self.statusbar.pop(0)
        self.statusbar.push(0, text)
        while gtk.events_pending():
            gtk.main_iteration(True)

    def update_cb(self, percent):
        """Callback for file loading code to inform the GUI of how
        of the file has been read
        """
        self.progress.set_fraction(percent/100.0)
        while gtk.events_pending():
            gtk.main_iteration(True)


class mmCIFEditor:
    """Implements a mmCIF file editor using logged transactions and
    notifications.  Transactions are listed below.

    Operations on mmCIFRow:
      cif_row_set_value(cif_row, col_name, value)
      cif_row_remove(cif_row)

    Operations on mmCIFTable:
      cif_table_insert_row(cif_table, i, cif_row)
      cif_table_set_name(cif_table, name)
      cif_table_remove(cif_table)

      cif_table_column_set_name(cif_table, old_col_name, new_col_name)
      cif_table_column_insert(cif_table, i, col_name)
      cif_table_column_remove(cif_table, col_name)

    Operations on mmCIFData:
      cif_data_insert_table(cif_data, i, cif_table)
      cif_data_set_name(cif_data, name)
      cif_data_remove(cif_data)

    Operations on mmCIFFile:
      cif_file_insert_data(cif_data, i, cif_data)
    """
    def __init__(self):
        self.undo_list = []

    def clear_undo_stack(self):
        """Clears out the undo stack.
        """
        self.undo_list = []

    def undo(self):
        """Pop one undo operation from the undo stack and apply.
        """
        if self.undo_list:
            undo_op = self.undo_list.pop()
            func = undo_op[0]
            args = undo_op[1:]
            func(*args)

    def cif_row_set_value(self, cif_row, col_name, value, save_undo = True):
        """Sets the value of the cif_row[col_name] to the given value.
        """
        assert isinstance(cif_row, mmCIF.mmCIFRow)

        if save_undo:
            undo = (self.cif_row_set_value_undo,
                    cif_row, col_name, cif_row.get(col_name))
            self.undo_list.append(undo)

        cif_row[col_name] = value
        self.cif_row_set_value_notify(cif_row, col_name)

    def cif_row_set_value_undo(self, cif_row, col_name, value):
        assert isinstance(cif_row, mmCIF.mmCIFRow)
        self.cif_row_set_value(cif_row, col_name, value, save_undo = False)

    def cif_row_set_value_notify(self, cif_row, col_name):
        pass

    def cif_row_remove(self, cif_row, save_undo = True):
        """Removes the cif_row from self.cif_fil.
        """
        assert isinstance(cif_row, mmCIF.mmCIFRow)

        if save_undo:
            undo = (self.cif_row_remove_undo,
                    cif_row.table.index(cif_row),
                    cif_row)
            self.undo_list.append(undo)

        self.cif_row_remove_notify(cif_row)
        cif_row.table.remove(cif_row)
        self.cif_row_remove_complete_notify(cif_row)

    def cif_row_remove_notify(self, cif_row):
        pass
    def cif_row_remove_complete_notify(self, cif_row):
        pass

    def cif_row_remove_undo(self, i, cif_row):
        assert isinstance(cif_row, mmCIF.mmCIFRow)
        self.cif_table_insert_row(cif_row.table, i, cif_row, save_undo = False)

    def cif_table_insert_row(self, cif_table, i, cif_row, save_undo = True):
        """Inserts a new cif_row into cif_table at position i.
        """
        assert isinstance(cif_table, mmCIF.mmCIFTable)
        assert isinstance(cif_row, mmCIF.mmCIFRow)

        if save_undo:
            undo = (self.cif_table_insert_row_undo, cif_row)
            self.undo_list.append(undo)

        if i == -1:
            cif_table.append(cif_row)
        else:
            cif_table.insert(i, cif_row)
        self.cif_table_insert_row_notify(cif_row)

    def cif_table_insert_row_undo(self, cif_row):
        assert isinstance(cif_row, mmCIF.mmCIFRow)
        self.cif_row_remove(cif_row, save_undo = False)

    def cif_table_insert_row_notify(self, cif_row):
        pass

    def cif_data_insert_table(self, cif_data, i, cif_table, save_undo = True):
        """Inserts a new cif_table into cif_data at position i.
        """
        assert isinstance(cif_data, mmCIF.mmCIFData)
        assert isinstance(cif_table, mmCIF.mmCIFTable)

        if cif_data.has_key(cif_table.name):
            self.error("Insert Error: Table names must be unique.")
            return

        if save_undo:
            undo = (self.cif_data_insert_table_undo, cif_table)
            self.undo_list.append(undo)

        if i == -1:
            cif_data.append(cif_table)
        else:
            cif_data.insert(i, cif_table)
        self.cif_data_insert_table_notify(cif_table)

    def cif_data_insert_table_undo(self, cif_table):
        assert isinstance(cif_table, mmCIF.mmCIFTable)
        self.cif_table_remove(cif_table, save_undo = False)

    def cif_data_insert_table_notify(self, cif_table):
        pass

    def cif_table_set_name(self, cif_table, name, save_undo = True):
        """Sets the name of a cif_table.
        """
        assert isinstance(cif_table, mmCIF.mmCIFTable)

        if cif_table.data.has_key(name):
            self.error(
                "Name Change Error: Table names must be unique.")
            return

        if save_undo:
            undo = (self.cif_table_set_name_undo, cif_table, cif_table.name)
            self.undo_list.append(undo)

        cif_table.name = name
        self.cif_table_set_name_notify(cif_table)

    def cif_table_remove(self, cif_table, save_undo = True):
        """Removes the cif_table from the cif_file.
        """
        assert isinstance(cif_table, mmCIF.mmCIFTable)

        print "\n\n## cif_table: ",dir(cif_table),"\n\n"

        if save_undo:
            undo = (self.cif_table_remove_undo,
                    cif_table.data.index(cif_table),
                    cif_table)
            self.undo_list.append(undo)

        self.cif_table_remove_notify(cif_table)
        cif_table.data.remove(cif_table)
        self.cif_table_remove_complete_notify(cif_table)

    def cif_table_remove_notify(self, cif_table):
        pass
    def cif_table_remove_complete_notify(self, cif_table):
        pass

    def cif_table_remove_undo(self, i, cif_table):
        assert isinstance(cif_table, mmCIF.mmCIFTable)

        self.cif_data_insert_table(
            cif_table.data, i, cif_table, save_undo = False)

    def cif_table_set_name_undo(self, cif_table, old_name):
        assert isinstance(cif_table, mmCIF.mmCIFTable)
        self.cif_table_set_name(cif_table, old_name, save_undo = False)

    def cif_table_set_name_notify(self, cif_table):
        pass

    def cif_table_column_set_name(
        self, cif_table, i, col_name, save_undo = True):
        """Sets the cif_table.columns[i] = col_name
        """
        assert isinstance(cif_table, mmCIF.mmCIFTable)

        ## do not allow duplicate column names
        for j in xrange(len(cif_table.columns)):
            if i != j and cif_table.columns[j] == col_name:
                self.error(
                    "Name Change Error: Column names must be unique.")
                return

        if save_undo:
            undo = (self.cif_table_column_set_name_undo,
                    cif_table, i, cif_table.columns[i])
            self.undo_list.append(undo)

        orig_col_name = cif_table.columns[i]
        cif_table.columns[i] = col_name
        for cif_row in cif_table:
            try:
                cif_row[col_name] = cif_row[orig_col_name]
            except KeyError:
                pass
            else:
                del cif_row[orig_col_name]

        self.cif_table_column_set_name_notify(cif_table, col_name)

    def cif_table_column_set_name_undo(self, cif_table, i, old_name):
        assert isinstance(cif_table, mmCIF.mmCIFTable)
        self.cif_table_column_set_name(
            cif_table, i, old_name, save_undo = False)

    def cif_table_column_set_name_notify(self, cif_table, col_name):
        pass

    def cif_table_column_insert(
        self, cif_table, i, col_name, col_val_list = None, save_undo = True):
        """Calls cif_table.insert(i, col_name)
        """
        assert isinstance(cif_table, mmCIF.mmCIFTable)

        ## do not allow duplicate column names
        if col_name in cif_table.columns:
            self.error("Insert Error: Column names must be unique.")
            return

        if save_undo:
            undo = (self.cif_table_column_remove,
                    cif_table, col_name)
            self.undo_list.append(undo)

        if i == -1:
            cif_table.columns.append(col_name)
        else:
            cif_table.columns.insert(i, col_name)

        if col_val_list:
            for i in xrange(len(cif_table)):
                cif_table[i][col_name] = col_val_list[i]
        self.cif_table_column_insert_notify(cif_table, col_name)

    def cif_table_column_insert_undo(self, cif_table, col_name):
        assert isinstance(cif_table, mmCIF.mmCIFTable)
        self.cif_table_column_remove(cif_table, col_name, save_undo = False)

    def cif_table_column_insert_notify(self, cif_table, col_name):
        pass

    def cif_table_column_remove(self, cif_table, col_name, save_undo = True):
        """Removes a column from the cif_table.
        """
        assert isinstance(cif_table, mmCIF.mmCIFTable)

        if save_undo:
            rebuild_list = []
            for cif_row in cif_table:
                rebuild_list.append(cif_row.get(col_name))
            undo = (self.cif_table_column_remove_undo,
                    cif_table,
                    cif_table.columns.index(col_name),
                    col_name,
                    rebuild_list)
            self.undo_list.append(undo)

        self.cif_table_column_remove_notify(cif_table, col_name)
        cif_table.columns.remove(col_name)
        for cif_row in cif_table:
            try:
                del cif_row[col_name]
            except KeyError:
                pass
        self.cif_table_column_remove_complete_notify(cif_table, col_name)

    def cif_table_column_remove_notify(self, cif_table, col_name):
        pass
    def cif_table_column_remove_complete_notify(self, cif_table, col_name):
        pass

    def cif_table_column_remove_undo(
        self, cif_table, i, col_name, rebuild_list):
        assert isinstance(cif_table, mmCIF.mmCIFTable)
        self.cif_table_column_insert(cif_table, i, col_name, save_undo = False)
        for i in xrange(len(cif_table)):
            cif_table[i][col_name] = rebuild_list[i]

    def cif_file_insert_data(self, cif_file, i, cif_data, save_undo = True):
        """Inserts a new cif_data into cif_file at position i.
        """
        assert isinstance(cif_file, mmCIF.mmCIFFile)
        assert isinstance(cif_data, mmCIF.mmCIFData)

        if cif_file.has_key(cif_data.name):
            self.error(
                "Insert Error: Data block names must be unique.")
            return

        if save_undo:
            undo = (self.cif_file_insert_data_undo, cif_data)
            self.undo_list.append(undo)

        if i == -1:
            cif_file.append(cif_data)
        else:
            cif_file.insert(i, cif_data)
        self.cif_file_insert_data_notify(cif_data)

    def cif_file_insert_data_undo(self, cif_data):
        assert isinstance(cif_data, mmCIF.mmCIFData)
        self.cif_data_remove(cif_data, save_undo = False)

    def cif_file_insert_data_notify(self, cif_data):
        pass

    def cif_data_remove(self, cif_data, save_undo = True):
        """Removes the cif_row from self.cif_fil.
        """
        assert isinstance(cif_data, mmCIF.mmCIFData)
        if save_undo:
            undo = (self.cif_data_remove_undo,
                    cif_data.file.index(cif_data),
                    cif_data)
            self.undo_list.append(undo)

        self.cif_data_remove_notify(cif_data)
        cif_data.file.remove(cif_data)
        self.cif_data_remove_complete_notify(cif_data)

    def cif_data_remove_notify(self, cif_data):
        pass
    def cif_data_remove_complete_notify(self, cif_data):
        pass

    def cif_data_remove_undo(self, i, cif_data):
        assert isinstance(cif_data, mmCIF.mmCIFData)
        self.cif_file_insert_data(
            cif_data.file, i, cif_data, save_undo = False)

    def cif_data_set_name(self, cif_data, name, save_undo = True):
        """Sets the name of a cif_data
        """
        assert isinstance(cif_data, mmCIF.mmCIFData)
        if cif_data.file.has_key(name):
            self.error(
                "Name Change Error: Data block names must be unique.")
            return

        if save_undo:
            undo = (self.cif_data_set_name_undo, cif_data, cif_data.name)
            self.undo_list.append(undo)

        cif_data.name = name
        self.cif_data_set_name_notify(cif_data)

    def cif_data_set_name_undo(self, cif_data, old_name):
        assert isinstance(cif_data, mmCIF.mmCIFData)
        self.cif_data_set_name(cif_data, old_name, save_undo = False)

    def cif_data_set_name_notify(self, cif_data):
        pass



class mmCIFEditorWindowContext(mmCIFEditor):
    """Manages a single mmCIF editor window.
    """
    def __init__(self, path = None):
        mmCIFEditor.__init__(self)
        self.mw = mmCIFEditorMainWindow(self)

        ## the path and current cif_file insetance
        ## the open_into_context flag is set when a new window is opened,
        ## so that Open... will open into the current context if
        ## no editing has been done
        self.open_into_context = True
        self.path = path
        self.cif_file = mmCIF.mmCIFFile()
        ## set to False when the cif_file has been edited and not saved
        self.saved = True
        ## a list of edit dialogs open
        self.edit_dialog_list = []

        ## active initialization
        if self.path != None:
            self.open_file(path)

        self.update_enabled_menuitems()

    def open_edit_dialog(self, cif_row, col_name):
        """Opens a dialog with a text widget for editing larger entries.
        These dialogs are note model, so only allow one edit window for a
        row/col_name.
        """
        for edit_dialog in self.edit_dialog_list:
            if edit_dialog.cif_row == cif_row and \
               edit_dialog.col_name == col_name:
                edit_dialog.present()
                return

        edit_dialog = EditDialog(self, cif_row, col_name)
        self.edit_dialog_list.append(edit_dialog)

    def close_edit_dialog(self, edit_dialog):
        """Called to finalize the close of an edit dialog.
        """
        self.edit_dialog_list.remove(edit_dialog)

    def open_file(self, path):
        """Opens a CIF file for editing.
        """
        if self.open_into_context:

            if path == None or path == "":
                return
            try:
                fil = open(path, "r")
            except IOError:
                self.error("Unable to open file:\n%s" % (path))
                return

            try:
                self.cif_file.load_file(fil)
            except mmCIF.mmCIFError, e:
                self.error("mmCIF parse error:\n%s" % (e))
                return

            self.open_into_context = False
            self.path = path
            self.clear_undo_stack()

            self.mw.set_title("mmCIF Editor: %s" % (self.path))
            self.mw.set_status_bar("Loading File: %s" % (self.path))
            self.mw.file_ctrl.set_cif_file()
            self.mw.progress.set_fraction(0.0)
            self.mw.set_status_bar("")

            self.update_enabled_menuitems()
        else:
            APP.new_editor_context(path)

    def save_file(self, path = None):
        """Save the mmCIF file, with path implements Save As...
        """
        assert path != ""
        ## case 1: Save As... sets self.path
        if self.path == None and path != None:
            save_path = path
        ## case 2: Save
        elif self.path != None and path == None:
            save_path = self.path
        ## case 3: Save As...
        elif self.path != None and path != None:
            save_path = path

        self.mw.set_status_bar("Saving File...please wait.")
        try:
            self.cif_file.save_file(save_path)
        except IOError, err:
            self.error("ERROR: Save file %s failed!" % (save_path))
        else:
            ## back to case 1: Save As... sets filename
            if self.path == None:
                self.path = save_path
                self.mw.set_title("mmCIF Editor: %s" % (self.path))
            if save_path == self.path:
                self.saved = True
            self.mw.set_status_bar("")
            self.update_enabled_menuitems()

    def quit(self):
        """Called when the editors main window is destroyed.
        """
        APP.editor_context_closed(self)

    def error(self, text):
        """Display a modeal dialog containing an error message.
        """

        dialog = gtk.MessageDialog(
            self.mw,
            gtk.DIALOG_DESTROY_WITH_PARENT,
            gtk.MESSAGE_ERROR,
            gtk.BUTTONS_CLOSE,
            text)

        dialog.run()
        dialog.destroy()

    def changed(self):
        self.open_into_context = False
        self.saved = False
        self.update_enabled_menuitems()

    def update_enabled_menuitems(self):
        self.mw.update_menuitem_enable_status()


    def cif_table_insert_row_notify(self, cif_row):
        self.mw.table_ctrl.insert_cif_row(cif_row)
        self.changed()
    def cif_row_set_value_notify(self, cif_row, col_name):
        self.mw.table_ctrl.update_cif_row(cif_row)
        self.changed()
    def cif_row_remove_notify(self, cif_row):
        self.mw.table_ctrl.remove_cif_row(cif_row)
        self.changed()


    def cif_data_insert_table_notify(self, cif_table):
        self.mw.file_ctrl.cif_table_insert(cif_table)
        self.changed()
    def cif_table_set_name_notify(self, cif_table):
        self.mw.file_ctrl.cif_table_changed(cif_table)
        self.changed()
    def cif_table_remove_notify(self, cif_table):
        self.mw.file_ctrl.cif_table_remove(cif_table)
        self.changed()
    def cif_table_remove_complete_notify(self, cif_table):
        if self.mw.table_ctrl.cif_table == cif_table:
            self.mw.table_ctrl.set_cif_table(None)

    def cif_table_column_set_name_notify(self, cif_table, col_name):
        self.mw.file_ctrl.cif_table_column_changed(cif_table, col_name)
        self.changed()
    def cif_table_column_insert_notify(self, cif_table, col_name):
        self.mw.file_ctrl.cif_table_column_insert(cif_table, col_name)
        self.changed()
    def cif_table_column_remove_notify(self, cif_table, col_name):
        self.mw.file_ctrl.cif_table_column_remove(cif_table, col_name)
        self.changed()
    def cif_table_column_remove_complete_notify(self, cif_table, col_name):
        self.mw.file_ctrl.update_table_ctrl(cif_table, col_name)


    def cif_file_insert_data_notify(self, cif_data):
        self.mw.file_ctrl.cif_data_insert(cif_data)
        self.changed()
    def cif_data_set_name_notify(self, cif_data):
        self.mw.file_ctrl.cif_data_changed(cif_data)
        self.changed()
    def cif_data_remove_notify(self, cif_data):
        self.mw.file_ctrl.cif_data_remove(cif_data)
        self.changed()
    def cif_data_remove_complete_notify(self, cif_data):
        if self.mw.table_ctrl.cif_table in cif_data:
            self.mw.table_ctrl.set_cif_table(None)

    def error(self, text):
        dialog = gtk.MessageDialog(
            self.mw,
            gtk.DIALOG_DESTROY_WITH_PARENT,
            gtk.MESSAGE_ERROR,
            gtk.BUTTONS_CLOSE,
            text)

        dialog.run()
        dialog.destroy()

## </mmCIF EDITOR PANEL>


## <APPLICATION GLOABAL WINDOWS>

class AboutWindow(gtk.Dialog):
    def __init__(self):
        gtk.Dialog.__init__(self, "About mmCIF Editor", None, 0)
        self.add_button(gtk.STOCK_CLOSE, gtk.RESPONSE_CLOSE)
        self.connect("delete-event", self.delete_event_cb)
        self.connect("response", self.delete_event_cb)

        self.label = gtk.Label()
        self.vbox.add(self.label)

        self.label.set_markup(
            '<span size="large">mmCIF Editor v%s</span>\n' % (__version__) +
            "(c)2003 PyMMLib Development Group\n"
            "http://pymmlib.sourceforge.net/\n"
            "\n"
            "University of Washington\n"
            "Jay Painter\n"
            "Ethan Merritt")

    def delete_event_cb(self, *args):
        self.hide()
        return True

    def present(self):
        self.show_all()
        gtk.Dialog.present(self)


class HelpWindow(gtk.Window):
    def __init__(self):
        gtk.Window.__init__(self)
        self.connect("delete_event", self.delete_notify_cb)
        self.set_title("mmCIF Dictionary Help Browser")
        self.set_default_size(400, 300)

        self.sw = gtk.ScrolledWindow()
        self.add(self.sw)
        self.sw.set_border_width(5)
        self.sw.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)

        self.text_view = gtk.TextView()
        self.sw.add(self.text_view)        
        self.text_view.set_editable(False)
        self.text_view.set_justification(gtk.JUSTIFY_LEFT)

        self.set_text_message(
            "Click on any mmCIF file item\n"
            "for help on that item.")

        self.current_cif_save = None

    def delete_notify_cb(self, *args):
        self.hide()
        return True

    def present(self):
        self.show_all()
        gtk.Window.present(self)

    def set_text_message(self, text):
        """Sets a simple text message in the window.
        """
        buffer = gtk.TextBuffer()
        tag = buffer.create_tag("title")
        tag.set_property("size-points", 14.0)
        iter = buffer.get_start_iter()
        buffer.insert_with_tags_by_name(iter, text, "title")
        self.text_view.set_buffer(buffer)

    def set_cif_save(self, cif_save):
        """Sets the Help text from the cif_save.
        """
        if self.current_cif_save == cif_save:
            return

        title_text = "No help available."
        text = None
        default_value = None
        range = "RANGE (max, min):     "
        data_type = None
        type_conditions = None
        data_units = None
        mandatory = None
        dependent = ""
        detail = "EXAMPLE(S):\n"

        if cif_save.has_key("category"):
            try:
                title_text = cif_save["category"][0]["id"]
            except KeyError, IndexError:
                pass
            try:
                text = cif_save["category"][0]["description"]
            except KeyError, IndexError:
                pass

            for row in cif_save.get("category_examples", []):
                try:
                    detail += "\n" + row["detail"] + "\n"
                except KeyError:
                    pass
                try:
                    detail += row["case"]
                except KeyError:
                    pass

        elif cif_save.has_key("item"):
            ## get title text
            try:
                title_text = cif_save["item"][0]["name"]
            except KeyError, IndexError:
                pass
            ## get description
            try:
                text = cif_save["item_description"][0]["description"]
            except KeyError, IndexError:
                pass
            ## get example details and cases
            for row in cif_save.get("item_examples", []):
                try:
                    detail += row["detail"] + "\n"
                except KeyError:
                    pass
                try:
                    detail += row["case"] + "\n"
                except KeyError:
                    pass
            ## get default value
            try:
                default_value = cif_save["item_default"][0]["value"]
            except KeyError, IndexError:
                pass
            ## get data type
            try:
               data_type = cif_save["item_type"][0]["code"]
            except KeyError, IndexError:
                pass
            ## get data type conditions
            try:
                type_conditions = cif_save["item_type_conditions"][0]["code"]
            except KeyError, IndexError:
                pass
            ## get data units
            try:
                data_units = cif_save["item_units"][0]["code"]
            except KeyError, IndexError:
                pass
            ## get mandatory item code
            try:
                mandatory = cif_save["item"][0]["mandatory_code"]
            except KeyError, IndexError:
                pass
            ## get dependents
            for row in cif_save.get("item_dependent", []):
                 try:
                     dependent += "\n          " + row["dependent_name"]
                 except KeyError:
                     pass
            ## get range values
            for row in cif_save.get("item_range", []):
                maximum = None
                minimum = None
                try:
                    maximum = row["maximum"]
                    minimum = row["minimum"]
                except KeyError:
                    pass
                if((maximum != None) & (minimum != None)):
                    if(maximum == ""):
                        maximum = " . "
                    range += "(" + maximum + ", " + minimum + ")  "

        ## take out alignment spaces
        ltmp = text.split("\n")
        ndiscard = 80

        for ln in ltmp:
            if len(ln) == 0:
                continue
            n = len(ln) - len(ln.lstrip())
            ndiscard = min(ndiscard, n)
        ltmp = [x[ndiscard:] for x in ltmp]
        text = string.join(ltmp, "\n")

        ltmp1 = detail.split("\n")
        ndiscard = 80

        for ln in ltmp1:
            if len(ln) == 0:
                continue
            n = len(ln) - len(ln.lstrip())
            ndiscard = min(ndiscard, n)
        ltmp1 = [x[ndiscard:] for x in ltmp1]
        detail = string.join(ltmp1, "\n")

        buffer = gtk.TextBuffer()

        fontsize = buffer.create_tag("title")
        fontsize.set_property("size-points", 14.0)

        wordwrap = buffer.create_tag("wrap")
        wordwrap.set_property("wrap-mode", gtk.WRAP_WORD)

        iter = buffer.get_start_iter()

        buffer.insert_with_tags_by_name(
            iter, title_text + "\n", "title")

        if((mandatory != None) & (mandatory != "no")):
            buffer.insert_with_tags_by_name(
                iter, "***Item is mandatory***\n", "wrap")
        buffer.insert_with_tags_by_name(
            iter, "\n\n" + text + "\n\n\n", "wrap")        

        if(data_type != None):
            buffer.insert_with_tags_by_name(
                iter, "DATA TYPE:     " + data_type + "\n", "wrap")

        if(type_conditions != None):
            buffer.insert_with_tags_by_name(
                iter, "DATA TYPE CONDINTIONS:     " + \
                type_conditions + "\n", "wrap")

        if(data_units != None):
            buffer.insert_with_tags_by_name(
                iter, "DATA UNITS:     " + data_units + "\n", "wrap")

        if(default_value != None):            
            buffer.insert_with_tags_by_name(
                iter, "DEFAULT VALUE:     " + default_value + "\n", "wrap")

        if(range != "RANGE (max, min):     "):
            buffer.insert_with_tags_by_name(
                iter, range + "\n", "wrap")

        if(dependent != ""):
            buffer.insert_with_tags_by_name(
                iter, "DEPENDENT(S):     " + dependent + "\n", "wrap")

        if(detail != "EXAMPLE(S):\n"):
             buffer.insert_with_tags_by_name(iter, detail, "wrap")

        self.text_view.set_buffer(buffer)


class mmCIFDictionaryManager(list):
    def __init__(self):
        list.__init__(self)

        if APP.pref.has_key("dictionary path list"):
            for path in APP.pref["dictionary path list"]:
                self.load_dictionary(path)
        else:
            APP.pref["dictionary path list"] = []

    def load_dictionary(self, path):
        """Loads a mmCIF dictionary into the manager
        """
        try:
            fil = open(path, "r")
        except IOError:
            return None

        cif_dict = mmCIF.mmCIFDictionary()
        cif_dict.path = path
        cif_dict.load_file(path)

        self.append(cif_dict)
        return cif_dict

    def add_dictionary(self, path):
        """Add a dictionary to the manager, updates preferences.
        """
        cif_dict = self.load_dictionary(path)
        if cif_dict:
            APP.pref["dictionary path list"].append(cif_dict.path)
            APP.save_preferences()

    def remove_dictionary(self, path):
        """Remove a dictionary from the manager by its path, updates
        preferences.
        """
        for cif_dict in self:
            if path == cif_dict.path:
                self.remove(cif_dict)
                APP.pref["dictionary path list"].remove(cif_dict.path)
                APP.save_preferences()
                break

    def lookup_cif_save(self, tag):
        """Returns the first save block found in the list of dictionaries
        """
        for cif_save in self:
            try:
                return cif_save[tag]
            except KeyError:
                pass
        return None


class DictionaryManagerWindow(gtk.Window):
    def __init__(self, dict_manager):
        self.dict_manager = dict_manager

        gtk.Window.__init__(self)
        self.set_default_size(*LARGE_DIALOG_SIZE)
        self.connect("delete_event", self.delete_notify_cb)
        self.set_title("mmCIF Dictionary Manager")

        self.vbox = gtk.VBox()
        self.add(self.vbox)
        self.vbox.set_border_width(5)

        self.toolbar = gtk.Toolbar()
        self.vbox.pack_start(self.toolbar, False, False, 0)
        self.toolbar.insert_stock(
            gtk.STOCK_ADD,
            "Add mmCIF Dictionary",
            None,
            self.add_dict_cb,
            None,
            -1)
        self.toolbar.append_space()
        self.toolbar.insert_stock(
            gtk.STOCK_REMOVE,
            "Remove mmCIF Dictionary",
            None,
            self.remove_dict_cb,
            None,
            -1)

        sw = gtk.ScrolledWindow()
        self.vbox.pack_start(sw, True, True, 0)
        sw.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)

        self.tree_view = gtk.TreeView()
        sw.add(self.tree_view)
        self.tree_view.get_selection().set_mode(gtk.SELECTION_SINGLE)

        cell = gtk.CellRendererText()
        self.column = gtk.TreeViewColumn("CIF Dictionaries", cell, markup=0)
        self.tree_view.append_column(self.column)

        self.model = gtk.ListStore(gobject.TYPE_STRING)
        self.tree_view.set_model(self.model)

        self.update_tree_view()

    def delete_notify_cb(self, *args):
        self.hide()
        return True

    def get_selected_dictionary_index(self):
        sel = self.tree_view.get_selection()
        try:
            (model, iter) = sel.get_selected()
            (i, ) = self.model.get_path(iter)
        except ValueError:
            return None
        return i

    def add_dict_cb(self, *args):
        if not hasattr(self, "file_selector"):
            self.file_selector = gtk.FileSelection(
                "Select mmCIF Dictionary file (.dic)");
            self.file_selector.connect(
                "delete-event", self.fsel_destroy_cb)
            self.file_selector.ok_button.connect(
                "clicked", self.fsel_ok_cb)
            self.file_selector.cancel_button.connect(
                "clicked", self.fsel_cancel_cb)
            self.file_selector.show()
        else:
            self.file_selector.present()

    def fsel_destroy_cb(self, *args):
        self.file_selector.hide()
        return True

    def fsel_ok_cb(self, *args):
        """Called by the [OK] button on the FileSelector.
        """
        path = self.file_selector.get_filename()
        self.file_selector.hide()

        ### XXX: add try/except
        self.dict_manager.add_dictionary(path)
        self.update_tree_view()

    def fsel_cancel_cb(self, *args):
        """Called by the [Cancel] button on the FileSelector.
        """
        self.file_selector.hide()

    def remove_dict_cb(self, *args):
        i = self.get_selected_dictionary_index()
        if i != None:
            del self.dict_manager[i]
            self.update_tree_view()

    def present(self):
        self.show_all()
        gtk.Window.present(self)

    def update_tree_view(self):
        self.model.clear()
        for cif_dict in self.dict_manager:
            text = ""

            for cif_data in cif_dict:
                try:
                    text +=\
                    '<span underline="single" size="large">%s</span>\n' % (
                        cif_data["dictionary"][0]["title"])
                except KeyError, IndexError:
                    pass
                try:
                    text +=\
           'Version: <span style="italic" weight="bold">%s</span>\n' % (
                        cif_data["dictionary"][0]["version"])
                except KeyError, IndexError:
                    pass
                try:
                    text +=\
           'Datablock ID: <span style="italic" weight="bold">%s</span>\n' % (
                        cif_data["dictionary"][0]["datablock_id"])
                except KeyError, IndexError:
                    pass

            path = getattr(cif_dict, "path", "")
            text += 'Path: <span style="italic" weight="bold">%s</span>' % (
                path)

            self.model.set(self.model.append(), 0, text)


class CIFEditorApplication(list):
    """Once instance of this class manages the application process.
    """
    def initalize(self):
        ## load preferences
        self.load_preferences()

        ## CIF dictionary manager
        self.dict_manager = mmCIFDictionaryManager()

        ## windows which are global to the application session
        self.help_window = HelpWindow()
        self.about_window = AboutWindow()
        self.dict_manager_window = DictionaryManagerWindow(self.dict_manager)

    def user_conf_path(self):
        """Generate the configuration file path.
        """
        return os.path.expanduser(os.path.join("~", CIF_EDITOR_CONF))

    def load_preferences(self):
        """Loads the preferences file.
        """
        try:
            self.pref = cPickle.load(open(self.user_conf_path(), "r"))
        except:
            self.pref = {}

    def save_preferences(self):
        """Save all the preferences.
        """
        try:
            cPickle.dump(self.pref, open(self.user_conf_path(), "w"))
        except IOError:
            ##XXX: model dialog warning
            pass

    def new_editor_context(self, path = None):
        """Opens a new CIF editor context/window.
        """
        context = mmCIFEditorWindowContext(path)
        self.append(context)
        return context

    def editor_context_closed(self, context):
        """Callback whever a CIF editor context is closed.
        """
        self.remove(context)
        if not len(self):
            self.save_preferences()
            gtk.main_quit()

    def open_dictionary_manager_window(self, *args):
        self.dict_manager_window.present()

    def close_dictionary_manager_window(self, *args):
        self.dict_manager_window.hide()

    def open_about_window(self, *args):
        self.about_window.present()

    def close_about_window(self, *args):
        self.about_window.hide()

    def open_help_window(self, *args):
        self.help_window.present()

    def close_help_window(self, *args):
        self.help_window.hide()

    def set_help_window(self, tag):
        cif_save = self.dict_manager.lookup_cif_save(tag)
        if cif_save:
            self.help_window.set_cif_save(cif_save)
        else:
            self.help_window.set_text_message(
                "Sorry, the descrption for the item\n"
                "could not be found in any of your listed\n"
                "mmCIF dictionary files.\n"
                "\n"
                "Item: " + tag)



## <MAIN>
if __name__=="__main__":
    global APP
    APP = CIFEditorApplication()
    APP.initalize()

    ## open windows for each path in the path_list
    if len(sys.argv[1:]):
        for path in sys.argv[1:]:
            APP.new_editor_context(path)
    else:
        APP.new_editor_context()

    gtk.main()
    sys.exit(0)
## </MAIN>

