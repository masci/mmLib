#!/usr/bin/env python

## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

## GtkCIFView.py
## a utility for viewing mmCIF files

## use Python-GTK bindings for GTK+-1.2
try:
    import pygtk
    pygtk.require('1.2')
except:
    pass

from gtk import *

import os
import sys
import string

from mmLib import mmCIF



class MainWindow:
    def __init__(self):
        self.window = GtkWindow(WINDOW_TOPLEVEL)
        self.window.connect("delete_event", self.window_delete_event)
        self.window.set_usize(200, 200)
	self.window.set_uposition(20, 20)

        vbox = GtkVBox()
        self.window.add(vbox)
        vbox.show()

        self.menubar = GtkMenuBar()
        vbox.pack_start(self.menubar, expand = FALSE)
        self.menubar.show()

        ## [File] Menu
        file_mi = GtkMenuItem("File")
        self.menubar.append(file_mi)
        file_mi.show()

        menu = GtkMenu()
        file_mi.set_submenu(menu)
        menu.show()
        
        mi = GtkMenuItem("Open")
        mi.connect("activate", self.file_open)
        menu.append(mi)
        mi.show()
        
        mi = GtkMenuItem("Quit")
        mi.connect("activate", self.file_quit)
        menu.append(mi)
        mi.show()


        ## [View] Menu
        file_mi = GtkMenuItem("View")
        self.menubar.append(file_mi)
        file_mi.show()

        menu = GtkMenu()
        file_mi.set_submenu(menu)
        menu.show()
        
        mi = GtkMenuItem("Structure Window")
        mi.connect("activate", self.view_structure)
        menu.append(mi)
        mi.show()
        

        ## table list area
        vpaned = GtkVPaned()
        vbox.add(vpaned, expand = TRUE)
        vpaned.set_border_width(5)
        vpaned.show()
        
        self.hpaned = GtkHPaned()
        vpaned.add1(self.hpaned)
        self.hpaned.set_border_width(5)
        self.hpaned.show()

        swin = GtkScrolledWindow()
        swin.set_policy(POLICY_AUTOMATIC, POLICY_AUTOMATIC)
        self.hpaned.add1(swin)
        swin.show()

        self.table_clist = GtkCList(1, ["mmCIF Section"])
        swin.add(self.table_clist)
        self.table_clist.column_titles_passive()
        self.table_clist.set_column_resizeable(1, FALSE)
        self.table_clist.set_column_width(0, 100)
        self.table_clist.set_selection_mode(SELECTION_SINGLE)
        self.table_clist.connect("select_row", self.table_clist_select_row)
        self.table_clist.show()

        ## row list area
        self.rl_sw = GtkScrolledWindow()
        self.rl_sw.set_policy(POLICY_AUTOMATIC, POLICY_AUTOMATIC)
        self.hpaned.add2(self.rl_sw)
        self.rl_sw.show()

        ## help text window
        self.help_text = GtkText()
        self.help_text.set_editable(FALSE)
        vpaned.add2(self.help_text)
        self.help_text.show()

        self.window.show()

    def window_delete_event(self, win, event = None):
        mainquit()

    def file_open(self, mi):
        pass

    def file_quit(self, mi):
        mainquit()

    def view_structure(self, mi):
        pass

    def table_clist_select_row(self, clist, row, col, event):
        data = self.cif.getData(0)
        self.selected_table_name = self.table_clist.get_text(row, col)
        table = data.getTable(self.selected_table_name)

        if hasattr(self, "row_clist"):
            self.row_clist.destroy()
            self.row_clist = None

        clist = table.columnList()

        self.row_clist = GtkCList(len(table.columnList()),table.columnList())
        self.row_clist.set_column_resizeable(1, TRUE)
        self.rl_sw.add(self.row_clist)
        self.row_clist.set_selection_mode(SELECTION_BROWSE)
        self.row_clist.connect("select_row", self.row_clist_select_row)
        self.row_clist.show()

        self.row_clist.freeze()
        for row in table.rowList():
            l = [getattr(row, x)[:50] for x in table.columnList()]
            self.row_clist.append(l)
            
        self.row_clist.thaw()
        self.row_clist.show()

    def row_clist_select_row(self, clist, row, col, event):
        data = self.cif.getData(0)
        self.selected_column_name = self.row_clist.get_column_title(col)

        s = ""
        dname = "%s.%s" % (self.selected_table_name, self.selected_column_name)

        ddata = mmCIF.CIF_Dictionnary.getData(dname)
        if ddata:
            dtable = ddata.getTable("item_description")
            if dtable:
                robj = dtable.rowList()[0]
                s = robj.description

        if not s:
            s = "No help in dictionary for %s" % (dname)
        else:
            s = dname + "\n----------------\n" + string.strip(s)

        self.help_text.freeze()
        self.help_text.set_point(0)
        self.help_text.forward_delete(self.help_text.get_length())
        self.help_text.insert_defaults(s)
        self.help_text.thaw()
        

    def setCIF(self, cif):
        self.cif = cif
        self.window.set_title('GtkCIFView: %s' % (self.cif.getPath()))
        data = self.cif.getData(0)

        list = [table.name() for table in data.tableList()]
        list.sort()

        self.table_clist.freeze()
        for s in list:
            self.table_clist.append([s])
        self.table_clist.thaw()
            

def main():
    try:
        path = sys.argv[1]
    except IndexError:
        print "usage: %s <mmCIF file path>" % (sys.argv[0])
        sys.exit(1)

    mmCIF.LoadCIFDictionary("data/cif_mm.dic")
    cif = mmCIF.LoadCIFFile(path)
#    cif.debug()

    mw = MainWindow()
    mw.setCIF(cif)
    mainloop()

if __name__ == '__main__':
    main()
