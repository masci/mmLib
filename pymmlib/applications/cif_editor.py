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


class TableTreeModel(gtk.ListStore):
    def __init__(self, ctable):
       gtk.ListStore.__init__(self,*[gobject.TYPE_STRING]*len(ctable.columns))
       self.ctable = ctable

       rl = []
       for i in range(len(self.ctable.columns)):
           rl += [i, self.ctable.columns[i]]

       for row in self.ctable:
           srl = rl[:]

           for i in range(len(srl)/2):
               srl[(i*2)+1] = row[srl[(i*2)+1]]

           iter = self.append()
           self.set(iter, *srl)


class FileTreeModel(gtk.GenericTreeModel):
    def __init__(self, cif_file):
        gtk.GenericTreeModel.__init__(self)
        self.cif_file = cif_file

    def get_iter_root(self):
        return self.cif_file

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
            return (self.cif_file.index(node), )

        elif isinstance(node, mmCIFTable):
            for cif_data in self.cif_file:
                if node == cif_data.get(node.name):
                    di = self.cif_file.index(cif_data)
                    ti = cif_data.index(node)
                    return (di, ti)

    def on_get_iter(self, path):
        """returns the node corresponding to the given path."""
        if len(path) == 1:
            return self.cif_file[path[0]]

        elif len(path) == 2:
            return self.cif_file[path[0]][path[1]]
        
    def on_get_value(self, node, column):
	"""returns the value stored in a particular column for the node"""
	return node.name

    def on_iter_next(self, node):
	"""returns the next node at this level of the tree"""
        if isinstance(node, mmCIFFile):
            return None

        elif isinstance(node, mmCIFData):
            i = self.cif_file.index(node)
            try:
                return self.cif_file[i+1]
            except IndexError:
                return None

        elif isinstance(node, mmCIFTable):
            (di, ti) = self.on_get_path(node)
            try:
                return self.cif_file[di][ti+1]
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
            return self.cif_file

        elif isinstance(node, mmCIFTable):
            (di, ti) = self.on_get_path(node)
            return self.cif_file[di]


class CIFPanel(gtk.HPaned):
    def __init__(self):
        gtk.HPaned.__init__(self)
        self.set_border_width(2)

        ## LEFT HALF
        self.sw1 = gtk.ScrolledWindow()
        self.sw1.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        self.add1(self.sw1)

        self.tv1 = gtk.TreeView()
        self.sw1.add(self.tv1)
        self.tv1.connect("row_activated", self.tview_row_activated)
        self.tv1.set_model(gtk.GenericTreeModel())

        cell = gtk.CellRendererText()
        column = gtk.TreeViewColumn("data", cell, text=0)
        self.tv1.append_column(column)
        
        ## RIGHT HALF
        self.sw2 = gtk.ScrolledWindow()
        self.sw2.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        self.add2(self.sw2)
        
        self.tv2 = gtk.TreeView()
        self.tv2.set_rules_hint(gtk.TRUE)
        self.sw2.add(self.tv2)
        self.tv2.connect("row_activated", self.rview_row_activated)

    def set_cif_file(self, cif_file):
        self.cif_file = cif_file
        model = FileTreeModel(self.cif_file)
        self.tv1.set_model(model)

    def tview_row_activated(self, tree_view, path, column):
        ## find selected node
        model = tree_view.get_model()
        node = model.on_get_iter(path)
        
        if isinstance(node, mmCIFTable):
            model = TableTreeModel(node)
            self.tv2.set_model(model)

            for c in self.tv2.get_columns():
                self.tv2.remove_column(c)

            for cn in node.columns:
                cell = gtk.CellRendererText()
                column = gtk.TreeViewColumn(
                    cn.replace("_","__"), cell, text=node.columns.index(cn))
                self.tv2.append_column(column)


    def rview_row_activated(self, tree_view, path, column):
        pass
    

class MainWindow:
    def __init__(self, quit_notify_cb):
        self.quit_notify_cb = quit_notify_cb

        ## Create the toplevel window
        self.window = gtk.Window()
        self.set_title("")
        self.window.set_default_size(500, 400)
        self.window.connect('destroy', self.quit_notify_cb, self)

        table = gtk.Table(1, 3, gtk.FALSE)
        self.window.add(table)

        ## Create the menubar
        menubar = gtk.MenuBar()

        file = gtk.MenuItem("File")
        menu = gtk.Menu()

        mi   = gtk.MenuItem("Open")
        mi.connect("activate", self.open_cb)
        menu.add(mi)

        mi   = gtk.MenuItem("Quit")
        mi.connect("activate", self.quit_notify_cb, self)
        menu.add(mi)

        file.set_submenu(menu)
        menubar.add(file)

        table.attach(menubar,
                     # X direction              Y direction
                     0, 1,                      0, 1,
                     gtk.EXPAND | gtk.FILL,     0,
                     0,                         0)


        ## CIF display pandel widget
        self.cif_panel = CIFPanel()
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

        self.window.show_all()

    def open_cb(self, widget):
        file_selector = gtk.FileSelection("Select mmCIF file (.cif)");

        ok_button = file_selector.ok_button
        ok_button.connect("clicked", self.open_ok_cb, file_selector)

        cancel_button = file_selector.cancel_button
        cancel_button.connect("clicked", self.open_cancel_cb, file_selector)
        
        file_selector.show()

    def open_ok_cb(self, ok_button, file_selector):
        path = file_selector.get_filename()
        file_selector.destroy()
        APP.new_window(path)
        
    def open_cancel_cb(self, cancel_button, file_selector):
        file_selector.destroy()

    def set_title(self, title):
        title = "mmCIF Editor: " + title
        title = title[:50]
        self.window.set_title(title)

    def set_status_bar(self, text):
        self.statusbar.pop(0)
        self.statusbar.push(0, text)

    def load_file(self, path):
        self.set_title(path)
        self.set_status_bar("Loading %s, please wait..." % (path))
        while gtk.events_pending():
            gtk.main_iteration(gtk.TRUE)

        self.cif_file = mmCIFFile()
        self.cif_file.load_file(open(path, "r"))

        self.cif_panel.set_cif_file(self.cif_file)

        self.set_status_bar("")
        return gtk.FALSE


class CIFEditorApplication(list):
    """Once instance of this class manages the application process.
    """
    def __init__(self, path_list):
        if path_list:
            for path in path_list:
                self.new_window(path)
        else:
            self.new_window()

    def new_window(self, path = None):
        mw = MainWindow(self.quit_notify_cb)
        self.append(mw)
        if path:
            mw.load_file(path)
        return mw
    
    def quit_notify_cb(self, window, mw):
        """Callback whever a MainWindow is closed.
        """
        self.remove(mw)
        if not len(self):
            gtk.main_quit()


if __name__ == "__main__":
    global APP
    APP = CIFEditorApplication(sys.argv[1:])
    gtk.main()
    sys.exit(0)
                               
