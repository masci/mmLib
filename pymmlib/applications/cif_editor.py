#!/usr/bin/env python
## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
import sys
import math
import pygtk
pygtk.require("2.0")
import gobject
import gtk

from mmLib.mmCIF import *




## <mmCIF EDITOR PANEL>

class FileTreeModel(gtk.GenericTreeModel):
    """This shows the mmCIFData/mmCIFTable tree.

    Selections:
      -- Selection signals from this widget select the current viewd
         mmCIFTable.

      -- Drag/Drops of tables copy sections of tables from one
         CIF file to another.
    """
    def __init__(self, context, tree_view):
        self.context = context
        self.tree_view = tree_view
        gtk.GenericTreeModel.__init__(self)

    def get_iter_root(self):
        return self.context.cif_file

    def on_get_flags(self):
	"""returns the GtkTreeModelFlags for this particular type of model"""
	return 0

    def on_get_n_columns(self):
	"""returns the number of columns in the model"""
	return 1
    
    def on_get_column_type(self, index):
	"""returns the type of a column in the model"""
	return gobject.TYPE_STRING

    def on_get_path(self, node):
	"""returns the tree path (a tuple of indices at the various
	levels) for a particular node."""
        if isinstance(node, mmCIFData):
            return (node.file.index(node), )
        elif isinstance(node, mmCIFTable):
            return (node.data.file.index(node.data), node.data.index(node))

    def on_get_iter(self, path):
        """returns the node corresponding to the given path."""
        if len(path) == 1:
            return self.context.cif_file[path[0]]

        elif len(path) == 2:
            return self.context.cif_file[path[0]][path[1]]
        
    def on_get_value(self, node, column):
	"""returns the value stored in a particular column for the node"""
	return node.name

    def on_iter_next(self, node):
	"""returns the next node at this level of the tree"""
        if isinstance(node, mmCIFFile):
            return None
        elif isinstance(node, mmCIFData):
            try:
                return node.file[node.file.index(node)+1]
            except IndexError:
                return None
        elif isinstance(node, mmCIFTable):
            try:
                return node.data[node.data.index(node)+1]
            except IndexError:
                return None

    def on_iter_children(self, node):
	"""returns the first child of this node"""
        if isinstance(node, mmCIFFile) or isinstance(node, mmCIFData):
            try:
                return node[0]
            except IndexError:
                return None
        else:
            return None

    def on_iter_has_child(self, node):
	"""returns true if this node has children"""
        if isinstance(node, mmCIFFile) or isinstance(node, mmCIFData):
            rval = len(node) > 0
        else:
            rval = False
        return rval

    def on_iter_n_children(self, node):
	"""returns the number of children of this node"""
        if isinstance(node, mmCIFFile) or isinstance(node, mmCIFData):
            rval = len(node)
        else:
            rval = 0
        return rval

    def on_iter_nth_child(self, node, n):
	"""returns the nth child of this node"""
        if isinstance(node, mmCIFFile) or isinstance(node, mmCIFData):
            return node[n]
        else:
            return None

    def on_iter_parent(self, node):
	"""returns the parent of this node"""
        if isinstance(node, mmCIFFile):
            return None
        elif isinstance(node, mmCIFData):
            return node.file
        elif isinstance(node, mmCIFTable):
            return node.data


class TableListModel(gtk.ListStore):
    """Editing spreadsheet for a mmCIFTable.  For columns in the table which
    have a small number of charactors, edits are done within the spreadsheet,
    for columns with large amount of text, a model dialog is invoked for
    editing.
    """
    def __init__(self, context, cif_table):
        self.context = context
        self.cif_table = cif_table
        self.fade_list = []
        ## figure out how many columns will be needed by the ListStore,
        ## one for each column of data, pluse two for a gtk.TRUE and
        ## gtk.FALSE which control editable columns
        num_cols = len(self.cif_table.columns)
        gtk_cols = [gobject.TYPE_STRING, gobject.TYPE_BOOLEAN] * num_cols
        gtk.ListStore.__init__(self, *gtk_cols)
        self.reset_rows()

    def reset_tree_view_columns(self, tree_view):
        """Prepares the tree_view widget with the appropriate columns for
        displaying this model.
        """
        ## remove any existing columns
        for c in tree_view.get_columns():
            tree_view.remove_column(c)

        ## add columns based on self.cif_table.columns
        for col_name in self.cif_table.columns:
            cell = gtk.CellRendererText()
            cell.set_data("col_name", col_name)
            cell.connect("edited", self.cell_edited_cb)

            model_index = self.cif_table.columns.index(col_name) * 2

            column = gtk.TreeViewColumn(col_name.replace("_","__"), cell)
            column.add_attribute(cell, "markup", model_index)
            column.add_attribute(cell, "editable", model_index+1)

            tree_view.append_column(column)

    def reset_rows(self):
        """Removes all rows from the model and reloads from self.cif_table.
        """
        self.clear()
        for cif_row in self.cif_table:
            self.set_cif_row_values(cif_row, self.append())

    def on_get_path(self, cif_row):
        """Returns a model path for locating the cif_row in the model.
        """
        return (self.cif_table.index(cif_row), )

    def get_cif_row_from_path(self, path):
        """Returns the cif_row object presented in the given path of the
        model.
        """
        return self.cif_table[path[0]]

    def cell_edited_cb(self, cell, row, new_text):
        """Called when a cell has been edited, and a new value entered.
        """
        col_name = cell.get_data("col_name")
        cif_row_index = int(row)
        cif_row = self.cif_table[cif_row_index]

        ## only request the update if the value changed
        if not cif_row.has_key(col_name) or cif_row[col_name] != new_text:
            self.context.cif_row_set_value(cif_row, col_name, new_text)

    def fade_cb(self, data_tuple):
        """Callback for the ultra cool highlight then fade feature after
        a cif_row insert or update.
        """
        (cif_row, iter, color) = data_tuple

        call_again = gtk.FALSE
        for i in range(len(color)):
            if color[i] > 0:
                color[i] = max(0, color[i] - 5)
                call_again = gtk.TRUE

        if call_again == gtk.FALSE:
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
            iter = self.get_iter(self.on_get_path(cif_row))
        
        row_list = []

        for i in range(len(self.cif_table.columns)):
            col_name = self.cif_table.columns[i]
            text = cif_row.get(col_name, ".")
            model_index = i * 2
            editible = gtk.TRUE

            if len(text) > 40:
                text = text[:40] + '<span foreground="red">[more]</span>'
                editible = gtk.FALSE
            if "\n" in text:
                text = text.replace("\n", '<span foreground="blue">|</span>')
                editible = gtk.FALSE
            if color != None:
                text = '<span foreground="#%02x%02x%02x">%s</span>' % (
                    color[0],color[1],color[2], text)

            row_list += [model_index, text, model_index + 1, editible]

        self.set(iter, *row_list)        

    def insert_cif_row(self, cif_row):
        """Insert a row of data into the model.  The cif_row should already
        be inserted into the cif_table.
        """
        assert cif_row.table == self.cif_table
        pos = self.cif_table.index(cif_row)
        iter = self.insert(pos)
        self.set_cif_row_values(cif_row, iter)

    def remove_cif_row(self, cif_row):
        """Removes one row frm the model.  The cif_row must still be in
        cif_table.
        """
        assert cif_row.table == self.cif_table
        pos = self.cif_table.index(cif_row)
        iter = self.get_iter(pos)
        self.remove(iter)

    def update_cif_row(self, cif_row):
        """Updates the values of the cif_row in the model.
        """
        assert cif_row.table == self.cif_table
        pos = self.cif_table.index(cif_row)
        iter = self.get_iter(pos)
        colors = [0x00, 0xff, 0x00]
        self.set_cif_row_values(cif_row, iter, colors)
        gtk.timeout_add(50, self.fade_cb, (cif_row, iter, colors))


class mmCIFPanel(gtk.HPaned):
    """Paned widget to group together the FileView and TableView of the
    CIF editor.

    Behavior:

    This panel does not act on the cif_file, it only requests to the
    context.  The context then issues updates to the GUI.  The requests
    are made by the mmCIFEditor's interface.

    What I need to be able to do:
    -- last selected item: mmCIFData, mmCIFTable, or mmCIFRow (used for help)

    -- last selected item + focus used for delete

    -- checks
      -- current mmCIFData
      -- current mmCIFTable
      -- current mmCIFRow

    -- selections: upon selection, hierarchy should scroll to correct
                   position and become expanded/selected/visible
      -- mmCIFData
      -- mmCIFTable
      -- mmCIFRow

    -- updates: look at all cases in mmCIFEditor and handle by
                re-generating the viewer information
      -- mmCIFData:
        --rename: update file view, select renamed item
        --delete: update file view, deselect selected children
        --insert: update file view, select inserted item
      -- mmCIFTable:
        --rename: update file view, update table view label
        --delete: update file view, dselect selected children
        --insert: update file view, select inserted item
        --remove column: update table view
      -- mmCIFRow:
        --delete: update table view, dselect
        --insert: update table view, select
        --update(value changed): update table view
    """

    def get_last_selected(self):
        """Return the last mmCIF object selected.
        """
        pass
    
    def get_selections(self):
        """Returns the selected items?
        """
        pass

    def __init__(self, context):
        self.context = context

        gtk.HPaned.__init__(self)
        self.set_border_width(2)

        ## models used by the two TreeView widgets
        ## file model is the right hand tree widget for displaying the
        ## cif_data and cif_table sections of the cif_file, and
        ## table_model is used for displaying the cif_rows for a cif_table
        self.file_model = None
        self.table_model = None

        ## LEFT HALF
        self.sw1 = gtk.ScrolledWindow()
        self.sw1.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        self.add1(self.sw1)

        self.tv1 = gtk.TreeView()
        self.sw1.add(self.tv1)
        self.tv1.set_model(self.file_model)
        self.tv1.connect("row_activated", self.file_view_row_activated)
        self.tv1.connect("button-release-event", self.file_view_button_release)

        cell = gtk.CellRendererText()
        column = gtk.TreeViewColumn("data", cell, text=0)
        self.tv1.append_column(column)
        
        ## RIGHT HALF
        self.table_frame = gtk.Frame()
        self.add2(self.table_frame)

        self.sw2 = gtk.ScrolledWindow()
        self.sw2.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        self.table_frame.add(self.sw2)
        
        self.tv2 = gtk.TreeView()
        self.sw2.add(self.tv2)
        self.tv2.set_rules_hint(gtk.TRUE)
        self.tv2.connect("row_activated", self.table_view_row_activated)
        self.tv2.connect("key_press_event", self.table_view_key_press_event)

    def file_view_button_release(self, tree_view, bevent):
        (tpath, column, x, y) = tree_view.get_path_at_pos(bevent.x, bevent.y)
        return gtk.FALSE

    def file_view_row_activated(self, tree_view, path, column):
        """Retrieve selected node, then call the correct set method for the
        type.
        """
        model = tree_view.get_model()
        node = model.on_get_iter(path)
        if isinstance(node, mmCIFData):
            self.select_cif_data(node)
        elif isinstance(node, mmCIFTable):
            self.select_cif_table(node)

    def update_file_model(self):
        self.file_model = FileTreeModel(self.context, self.tv1)
        self.table_model = None
        self.tv1.set_model(self.file_model)
        self.tv2.set_model(self.table_model)

    def insert_cif_row(self, cif_row):
        """Notification of the insert of cif_row.
        """
        if cif_row.table == self.table_model.cif_table:
            self.table_model.insert_cif_row(cif_row)
        self.select_cif_row(cif_row)

    def remove_cif_row(self, cif_row):
        self.select_cif_row(cif_row)
        self.tabel_model.remove_cif_row(cif_row)

    def update_cif_row(self, cif_row):
        self.select_cif_row(cif_row)
        self.table_model.update_cif_row(cif_row)

    def select_cif_data(self, cif_data):
        """Select a cif_data instance.  If it is not visible, then make it
        visible.
        """
        file_path = self.file_model.on_get_path(cif_data)
        self.tv1.scroll_to_cell(file_path, None, gtk.FALSE, 0.0, 0.0)

    def select_cif_table(self, cif_table):
        """Select a cif_table instance.  If not visible, then make it
        visible.
        """
        if self.table_model and self.table_model.cif_table == cif_table:
            return
        
        file_path = self.file_model.on_get_path(cif_table)
        self.tv1.scroll_to_cell(file_path, None, gtk.FALSE, 0.0, 0.0)

        self.table_frame.set_label(
            "%s.%s" % (cif_table.data.name,cif_table.name))

        self.table_model = TableListModel(self.context, cif_table)
        self.table_model.reset_tree_view_columns(self.tv2)
        self.tv2.set_model(self.table_model)

    def select_cif_row(self, cif_row):
        """Select a cif_row instance.  If not visible, then make it visible.
        """
        self.select_cif_table(cif_row.table)
        row_path = self.table_model.on_get_path(cif_row)
        self.tv2.scroll_to_cell(row_path, None, gtk.FALSE, 0.0, 0.0)

    def table_view_row_activated(self, tree_view, path, column):
        assert tree_view == self.tv2

    def table_view_key_press_event(self, tree_view, kevent):
        assert tree_view == self.tv2

        selection = self.tv2.get_selection()
        (model, iter) =  selection.get_selected()
        path = model.get_path(iter)
        cif_row = model.get_cif_row_from_path(path)
        i = model.cif_table.index(cif_row) + 1
        self.context.cif_table_insert_row(model.cif_table, i, mmCIFRow())
        
        print "key", kevent.keyval
        return gtk.FALSE


class mmCIFEditorMainWindow(gtk.Window):
    """Main window class used as the center of action for a
    mmCIFEditorWindowContext.
    """
    def __init__(self, context):
        self.context = context
        gtk.Window.__init__(self)

        ## Create the toplevel window
        self.set_default_size(500, 400)
        self.connect("delete-event", self.quit_cb, self)

        table = gtk.Table(1, 3, gtk.FALSE)
        self.add(table)

        ## Create the menubar
        MENU_ITEMS = (
            ('/_File', None, None, 0, '<Branch>' ),
            ('/File/_New', '<control>N', self.new_cb, 0,
             '<StockItem>', gtk.STOCK_NEW),
            ('/File/_Open', '<control>O', self.open_cb, 0,
             '<StockItem>', gtk.STOCK_OPEN),
            ('/File/_Save', '<control>S', self.save_cb, 0,
             '<StockItem>', gtk.STOCK_SAVE),
            ('/File/Save _As...', None, self.save_as_cb, 0,
             '<StockItem>', gtk.STOCK_SAVE),
            ('/File/sep1', None, None, 0, '<Separator>'),
            ('/File/_Quit', '<control>Q', self.quit_cb, 0,
             '<StockItem>', gtk.STOCK_QUIT),

            ('/Edit/_Undo', '<control>X', self.undo_cb, 0,
             '<StockItem>', gtk.STOCK_UNDO),
            
            ('/_Help', None, None, 0, '<Branch>'),
            ('/Help/_About', None, APP.open_about_window, 0, ''),
            ('/Help/Help Browser', None, APP.open_help_window, 0, ''))

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
        self.cif_panel = mmCIFPanel(self.context)
        table.attach(self.cif_panel,
                     # X direction           Y direction
                     0, 1,                   1, 2,
                     gtk.EXPAND | gtk.FILL,  gtk.EXPAND | gtk.FILL,
                     0,                      0)

        ## Create statusbar 
        self.statusbar = gtk.Statusbar();
        table.attach(self.statusbar,
                     # X direction           Y direction
                     0, 1,                   2, 3,
                     gtk.EXPAND | gtk.FILL,  0,
                     0,                      0)

        self.show_all()

    def new_cb(self, *args):
        """/File/New callback.  Opens a new edit window.
        """
        self.context.open_file(None)

    def open_cb(self, *args):
        """/File/Open callback.  Open a FileSeletor window for
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
        return gtk.TRUE

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
        """/File/Save As callback.  Open a FileSeletor window for Save As...
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
        return gtk.TRUE

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
        return gtk.TRUE

    def undo_cb(self, *args):
        """/Edit/Undo callback.
        """
        self.context.undo()

    def enable_save(self, enable):
        """Enables all save menuitems.
        """
        menuitem = self.item_factory.get_item('/File/Save')
        if enable:
            menuitem.set_sensitive(gtk.TRUE)
        else:
            menuitem.set_sensitive(gtk.FALSE)

    def enable_undo(self, enable):
        """Enables undo menuitem.
        """
        menuitem = self.item_factory.get_item('/Edit/Undo')
        if enable:
            menuitem.set_sensitive(gtk.TRUE)
        else:
            menuitem.set_sensitive(gtk.FALSE)

    def set_status_bar(self, text):
        """Sets the text on the windows statusbar.
        """
        self.statusbar.pop(0)
        self.statusbar.push(0, text)
        while gtk.events_pending():
            gtk.main_iteration(gtk.TRUE)


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
        if save_undo:
            undo = (self.cif_row_set_value_undo,
                    cif_row, col_name, cif_row.get(col_name))
            self.undo_list.append(undo)

        cif_row[col_name] = value
        self.cif_row_set_value_notify(cif_row, col_name)

    def cif_row_set_value_undo(self, cif_row, col_name, value):
        self.cif_row_set_value(cif_row, col_name, value, False)

    def cif_row_set_value_notify(self, cif_row, col_name):
        pass
        
    def cif_row_remove(self, cif_row, save_undo = True):
        """Removes the cif_row from self.cif_fil.
        """
        if save_undo:
            undo = (self.cif_row_remove_undo,
                    cif_row.table.index(cif_row),
                    cif_row)
            self.undo_list.append(undo)

        self.cif_row_remove_notify(cif_row)
        cif_row.table.remove(cif_row)

    def cif_row_remove_undo(self, i, cif_row):
        self.insert_row(i, cif_row, False)
        
    def cif_row_remove_notify(self, cif_row):
        pass

    def cif_table_insert_row(self, cif_table, i, cif_row, save_undo = False):
        """Inserts a new cif_row into cif_table at position i.
        """
        if save_undo:
            undo = (self.cif_table_insert_row_undo, cif_row)
            self.undo_list.append(undo)

        cif_table.insert(i, cif_row)
        self.cif_table_insert_row_notify(cif_row)

    def cif_table_insert_row_undo(self, cif_row):
        self.cif_row_remove(cif_row, False)

    def cif_table_insert_row_notify(self, cif_row):
        pass

    def cif_table_set_name(self, cif_table, name, save_undo = True):
        """Sets the name of a cif_table.
        """
        if save_undo:
            undo = (self.cif_table_set_name, cif_table, cif_table.name)
            self.undo_list.append(undo)

        cif_table.name = name
        self.cif_table_set_name_notify(cif_table)

    def cif_table_set_name_notify(self, cif_table):
        pass


class mmCIFEditorWindowContext(mmCIFEditor):
    """Manages a single mmCIF editor window.
    """
    def __init__(self, path = None):
        mmCIFEditor.__init__(self)
        self.mw = mmCIFEditorMainWindow(self)

        ## set to False when the cif_file has been edited and not saved
        self.saved = True
        ## the path of the cif_file
        self.path = None

        ## active initialization
        if path != None:
            self.open_file(path)

        self.update_enabled_menuitems()
            
    def open_file(self, path):
        """Opens a CIF file for editing.
        """
        if self.path == None:
            self.path = path

            self.clear_undo_stack()

            self.mw.set_title("mmCIF Editor: %s" % (self.path))
            self.mw.set_status_bar("Loading File...please wait.")

            self.cif_file = mmCIFFile()
            self.cif_file.load_file(OpenFile(self.path, "r"))
            
            self.mw.cif_panel.update_file_model()
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
            self.mw.set_status_bar("ERROR: Save file %s failed!" % (save_path))
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

    def update_enabled_menuitems(self):
        """Called whenever a edit is made to the CIF file.  This function
        updates the enabled/disabled state of menuitems in the GUI.
        """
        self.mw.enable_save(not self.saved)
        self.mw.enable_undo(len(self.undo_list) > 0)

    def cif_row_set_value_notify(self, cif_row, col_name):
        self.mw.cif_panel.update_cif_row(cif_row)
        self.saved = False
        self.update_enabled_menuitems()

    def cif_row_remove_notify(self, cif_row):
        self.mw.cif_panel.remove_cif_row(cif_row)

    def cif_table_insert_row_notify(self, cif_row):
        self.mw.cif_panel.insert_cif_row(cif_row)


## </mmCIF EDITOR PANEL>




## <APPLICATION GLOABAL WINDOWS>

class AboutWindow(gtk.Window):
    def __init__(self, quit_notify):
        gtk.Window.__init__(self)
        self.connect("destroy", quit_notify)
        self.set_title("About mmCIF Editor")
        self.set_border_width(5)

        self.label = gtk.Label()
        self.add(self.label)

        self.label.set_text(
            """\
            mmCIF Editor v1.0
            (c)2003 PyMMLib Development Group
            http://pymmlib.sourceforge.net/
            
            University of Washington
            Jay Painter
            Ethan Merrit""")

        self.show_all()


class HelpWindow(gtk.Window):
    def __init__(self, close_notify):
        gtk.Window.__init__(self)
        self.connect("destroy", close_notify)
        self.set_title("mmCIF Editor Help Browser")
        self.set_default_size(250, 400)
        self.show_all()


class mmCIFDicionaryManager(list):
    def find_save(self):
        pass


class mmCIFDictionaryManagerWindow(gtk.Window):
    def __init__(self, close_notify):
        gtk.Window.__init__(self)
        self.connect("destroy", close_notiry)
        self.set_title("mmCIF Dictionary Manager")
        
        self.dict_manager = mmCIFDictionaryManager()

        self.show_all()
    

class CIFEditorApplication(list):
    """Once instance of this class manages the application process.
    """
    def __init__(self):
        ## windows which are global to the application session
        self.help_window = None
        self.about_window = None
        self.dict_manager_window = None

        ## random application properties, these should be saved and
        ## loaded
        self.prop = {}

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
            gtk.main_quit()

    def open_about_window(self, widget, junk1 = None, junk2 = None):
        if self.about_window == None:
            self.about_window = AboutWindow(self.close_about_window)
        self.about_window.present()

    def close_about_window(self, window):
        self.about_window.hide()
        
    def open_help_window(self, widget, junk1 = None, junk2 = None):
        if self.help_window == None:
            self.help_window = HelpWindow(self.close_help_window)
        self.help_window.present()

    def close_help_window(self, window):
        self.help_window.hide()

    def open_dict_manager_window(self, widget, junk1 = None, junk2 = None):
        if self.dict_manager_window == None:
            self.dict_manager_window = mmCIFDictionaryManagerWindow(
                self.close_dict_manager_window)
        self.dict_manager_window.present()

    def close_dict_manager_window(self, window):
        self.dict_manager_window.hide()

## <MAIN>

global APP
APP = CIFEditorApplication()

## open windows for each path in the path_list
if len(sys.argv[1:]):
    for path in sys.argv[1:]:
        APP.new_editor_context(path)
else:
    APP.new_editor_context()

gtk.main()
sys.exit(0)

## </MAIN>

