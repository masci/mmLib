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


class TableTreeModel(gtk.GenericTreeModel):
    def __init__(self, cif_file):
        self.cif_file = cif_file
	gtk.GenericTreeModel.__init__(self)

    def get_iter_root(self):
        return self.cif_file

    def on_get_flags(self):
	'''returns the GtkTreeModelFlags for this particular type of model
        '''
	return 0

    def on_get_n_columns(self):
	'''returns the number of columns in the model
        '''
	return 1

    def on_get_column_type(self, index):
	'''returns the type of a column in the model
        '''
	return gobject.TYPE_STRING

    def on_get_path(self, node):
	'''returns the tree path (a tuple of indices at the various
	levels) for a particular node.
        '''
        if isinstance(node, mmCIFData):
            return (self.cif_file.index(node), )

        elif isinstance(node, mmCIFTable):
            for cif_data in self.cif_file:
                if cif_data.get(node.name) == node:
                    return (self.cif_file.index(cif_data),
                            cif_data.index(node))

    def on_get_iter(self, path):
        '''returns the node corresponding to the given path.
        '''
        if len(path) == 1:
            return self.cif_file[path[0]]

        elif len(path) == 2:
            return self.cif_file[path[0]][path[1]]
        
    def on_get_value(self, node, column):
	'''returns the value stored in a particular column for the node
        '''
	return node.name

    def on_iter_next(self, node):
	'''returns the next node at this level of the tree
        '''
        if isinstance(node, mmCIFData):
            try:
                return self.cif_file[self.cif_file.index(node)+1]
            except IndexError:
                return None

        elif isinstance(node, mmCIFTable):
            (di, ti) = self.on_get_path(node)
            try:
                self.cif_file[di][ti+1]
            except IndexError:
                return None

    def on_iter_children(self, node):
	'''returns the first child of this node
        '''
        if node == self.cif_file:
            try:
                return self.cif_file[0]
            except IndexError:
                return None

        elif isinstance(node, mmCIFData):
            try:
                return node[0]
            except IndexError:
                return None

        else:
            return None

    def on_iter_has_child(self, node):
	'''returns true if this node has children
        '''
        return len(node) != 0

    def on_iter_n_children(self, node):
	'''returns the number of children of this node
        '''
        return len(node)

    def on_iter_nth_child(self, node, n):
	'''returns the nth child of this node
        '''
        return node[n]

    def on_iter_parent(self, node):
	'''returns the parent of this node
        '''
        if node == self.cif_file:
            return None

        elif isinstance(node, mmCIFData):
            return self.cif_file

        elif isinstance(node, mmCIFTable):
            (di, ti) = self.on_get_path(node)
            return


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
        self.tv1.connect("row_activated", self.row_activated)
        cell = gtk.CellRendererText()
        column = gtk.TreeViewColumn("Tables", cell, text=0)
        self.tv1.append_column(column)

        ## RIGHT HALF
        self.tv2 = None

    def set_cif_file(self, cif_file):
        self.cif_file = cif_file
        
        model = TableTreeModel(self.cif_file)
        self.tv1.set_model(model)
        self.struct_panel.set_structure(self.struct)

    def row_activated(self, tree_view, path, column):
        ## find selected node
        model = tree_view.get_model()
        node = model.on_get_iter(path)
        self.tv2.set_structure(self.struct, node)


class MainWindow:
    def __init__(self, quit_notify_cb):
        self.quit_notify_cb = quit_notify_cb

        ## Create the toplevel window
        self.window = gtk.Window()
        self.setTitle("")
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


        ##
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
        file_selector = gtk.FileSelection("Select file to view");

        ok_button = file_selector.ok_button
        ok_button.connect("clicked", self.open_ok_cb, file_selector)

        cancel_button = file_selector.cancel_button
        cancel_button.connect("clicked", self.open_cancel_cb, file_selector)
        
        file_selector.show()


    def open_ok_cb(self, ok_button, file_selector):
        path = file_selector.get_filename()
        file_selector.destroy()
        self.load_file(path)
        

    def open_cancel_cb(self, cancel_button, file_selector):
        file_selector.destroy()


    def setTitle(self, title):
        title = "mmCIF Editor: " + title
        title = title[:50]
        self.window.set_title(title)


    def setStatusBar(self, text):
        self.statusbar.pop(0)
        self.statusbar.push(0, text)


    def load_file(self, path):
        self.setTitle(path)
        self.setStatusBar("Loading %s, please wait..." % (path))
        while gtk.events_pending():
            gtk.main_iteration(gtk.TRUE)

        ## LOAD FILE HERE!

        self.setStatusBar("")
        return gtk.FALSE


class CIFEditorApplication(list):
    """Once instance of this class manages the application process.
    """
    def __init__(self, path_list):
        if path_list:
            for path in path_list:
                mw = MainWindow(self.quit_notify_cb)
                self.append(mw)
                mw.load_file(path)
        else:
            mw = MainWindow(self.quit_notify_cb)
            self.append(mw)
            
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
                               
