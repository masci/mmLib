#!/usr/bin/env python
## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import math

import pygtk
pygtk.require("2.0")

import gobject
import gtk
import gtk.gtkgl

from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

from mmLib.Structure     import *
from mmLib.FileLoader    import LoadStructure, SaveStructure
from mmLib.GLViewer      import *
from mmLib.Extensions.PenultimateRotamers import FindBestRotamer

try:
    # try double-buffered
    glconfig = gtk.gdkgl.Config(mode = gtk.gdkgl.MODE_RGB |
                                gtk.gdkgl.MODE_DOUBLE |
                                gtk.gdkgl.MODE_DEPTH)
except gtk.gdkgl.NoMatches:
    # try single-buffered
    glconfig = gtk.gdkgl.Config(mode = gtk.gdkgl.MODE_RGB |
                                gtk.gdkgl.MODE_DEPTH)

class GtkGLViewer(gtk.gtkgl.DrawingArea):
    def __init__(self):
        gtk.gtkgl.DrawingArea.__init__(self)
        gtk.gtkgl.widget_set_gl_capability(self, glconfig)
        
        self.set_size_request(300, 300)

        self.connect('button_press_event',  self.button_press_event)
        self.connect('motion_notify_event', self.motion_notify_event)
        self.connect('realize',             self.realize)
        self.connect('map',                 self.map)
        self.connect('unmap',               self.unmap)
        self.connect('configure_event',     self.configure_event)
        self.connect('expose_event',        self.expose_event)
        self.connect('destroy',             self.destroy)

        self.glviewer        = None
        self.unit_cell_draw  = None
        self.atom_draw       = None

    def destroy(self, glarea):
        return gtk.TRUE

    def button_press_event(self, glarea, event):
        self.beginx = event.x
        self.beginy = event.y

    def motion_notify_event(self, glarea, event):
        width     = glarea.allocation.width
        height    = glarea.allocation.height

        if (event.state & gtk.gdk.BUTTON1_MASK):
            self.glviewer.roty = self.glviewer.roty + \
                        (((event.x - self.beginx) / float(width)) * 360.0)
            self.glviewer.rotx = self.glviewer.rotx + \
                        (((event.y - self.beginy) / float(height)) * 360.0)

        elif (event.state & gtk.gdk.BUTTON2_MASK):
            self.glviewer.zpos = self.glviewer.zpos + \
                        (((event.y - self.beginy) / float(height)) * 50.0)

            self.glviewer.zpos = max(self.glviewer.zpos, -450.0)

        elif (event.state & gtk.gdk.BUTTON3_MASK):
            self.glviewer.ypos = self.glviewer.ypos - \
                        (((event.y - self.beginy) / float(height)) * 50.0)
            self.glviewer.xpos = self.glviewer.xpos + \
                        (((event.x - self.beginx) / float(width)) * 50.0)

        self.beginx = event.x
        self.beginy = event.y

        self.queue_draw()

    def configure_event(self, glarea, event):
        if not self.glviewer:
            self.glviewer = GLViewer(self.get_gl_context(),
                                     self.get_gl_drawable())

        x, y, width, height = glarea.get_allocation()
        self.glviewer.gl_resize(width, height)

        self.queue_draw()
        return gtk.TRUE


    def map(self, glarea):
        self.add_events(gtk.gdk.BUTTON_PRESS_MASK   |
                        gtk.gdk.BUTTON_RELEASE_MASK |
                        gtk.gdk.BUTTON_MOTION_MASK  |
                        gtk.gdk.POINTER_MOTION_MASK)
        return gtk.TRUE

    def unmap(self, glarea):
        return gtk.TRUE

    def realize(self, glarea):
        self.glviewer.gl_init()
        return gtk.TRUE

    def expose_event(self, glarea, event):
        self.glviewer.gl_draw_lists()
        return gtk.TRUE

    def append_struct(self, struct):
        """Adds a structure to this viewer, and returns the GLStructure
        object so it can be manipulated.
        """
        gl_struct = GLStructure(struct)
        gl_struct.show_axes(True)
        gl_struct.show_unit_cell(True)
        gl_struct.show_atoms()
        
        self.glviewer.append(gl_struct)
        self.queue_draw()

        return gl_struct

    def remove_struct(self, struct):
        """Removes structure from the viewer.
        """
        pass


class StructureDataDialog(gtk.VPaned):
    def __init__(self):
        gtk.VPaned.__init__(self)
        self.set_border_width(3)

        ## make the print box
        sw = gtk.ScrolledWindow()
        sw.set_shadow_type(gtk.SHADOW_ETCHED_IN)
        sw.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        self.add1(sw)

        self.store = gtk.ListStore(gobject.TYPE_STRING, gobject.TYPE_STRING)
     
        treeview = gtk.TreeView(self.store)
        sw.add(treeview)
        treeview.set_rules_hint(gtk.TRUE)
        treeview.set_search_column(0)

        column = gtk.TreeViewColumn("mmLib.Structure Method",
                                    gtk.CellRendererText(), text=0)
        treeview.append_column(column)
    
        column = gtk.TreeViewColumn("Return Value",
                                    gtk.CellRendererText(), text=1)
        treeview.append_column(column)

        ## create the GLViewer
        self.glviewer = GtkGLViewer()
        self.add2(self.glviewer)

    def add_line(self, key, value):
        iter = self.store.append()
        self.store.set(iter, 0, key, 1, str(value))

    def set_structure(self, struct, sel_struct_obj = None):
        self.struct = struct
        self.store.clear()

        if isinstance(sel_struct_obj, Residue):
            self.add_line("Residue.res_name", sel_struct_obj.res_name)
            self.add_line("Residue.get_offset_residue(-1)",
                          sel_struct_obj.get_offset_residue(-1))
            self.add_line("Residue.get_offset_residue(1)",
                          sel_struct_obj.get_offset_residue(1))

        if isinstance(sel_struct_obj, Fragment):
            bonds = ""
            for bond in sel_struct_obj.iter_bonds():
                bonds += str(bond)
            self.add_line("Fragment.iter_bonds()", bonds)

        if isinstance(sel_struct_obj, AminoAcidResidue):
            self.add_line("AminoAcidResidue.calc_mainchain_bond_length()",
                          sel_struct_obj.calc_mainchain_bond_length())
            self.add_line("AminoAcidResidue.calc_mainchain_bond_angle()",
                          sel_struct_obj.calc_mainchain_bond_angle())
            self.add_line("AminoAcidResidue.calc_torsion_psi()",
                          sel_struct_obj.calc_torsion_psi())
            self.add_line("AminoAcidResidue.calc_torsion_phi()",
                          sel_struct_obj.calc_torsion_phi())
            self.add_line("AminoAcidResidue.calc_torsion_omega()",
                          sel_struct_obj.calc_torsion_omega())
            self.add_line("AminoAcidResidue.calc_torsion_chi()",
                          sel_struct_obj.calc_torsion_chi())

        if isinstance(sel_struct_obj, Atom):
            self.add_line("Atom.element", sel_struct_obj.element)
            self.add_line("Atom.name", sel_struct_obj.name)
            self.add_line("Atom.occupancy", sel_struct_obj.occupancy)
            self.add_line("Atom.temp_factor", sel_struct_obj.temp_factor)
            self.add_line("Atom.U", sel_struct_obj.U)
            self.add_line("Atom.position", sel_struct_obj.position)
            self.add_line("len(Atom.bond_list)", len(sel_struct_obj.bond_list))
            self.add_line("Atom.calc_anisotropy()",
                          sel_struct_obj.calc_anisotropy())

        self.glviewer.add_structure(self.struct)


class StructureTreeControl(gtk.TreeView):
    def __init__(self, context):
        self.context     = context
        self.struct_list = []
        
        gtk.TreeView.__init__(self)
        self.get_selection().set_mode(gtk.SELECTION_BROWSE)
        
        self.connect("row-activated", self.row_activated_cb)
        self.connect("button-release-event", self.button_release_event_cb)

        self.model = gtk.TreeStore(
            gobject.TYPE_BOOLEAN,   # 0: viewable
            gobject.TYPE_STRING)    # 1: data name
        self.set_model(self.model)

        cell = gtk.CellRendererText()
        column = gtk.TreeViewColumn("Structure", cell)
        column.add_attribute(cell, "text", 1)
        self.append_column(column)

    def row_activated_cb(self, tree_view, path, column):
        """Retrieve selected node, then call the correct set method for the
        type.
        """
        pass

    
    def button_release_event_cb(self, tree_view, bevent):
        """
        """
        pass
        
    def redraw(self):
        """Clear and refresh the view of the widget according to the
        self.struct_list
        """
        self.model.clear()
        for struct in self.struct_list:
            self.display_struct(struct)

    def display_struct(self, struct):
        iter1 = self.model.append(None)
        self.model.set(iter1, 1, str(struct))

        for chain in struct.iter_chains():
            iter2 = self.model.append(iter1)
            self.model.set(iter2, 1, str(chain))

            for frag in chain.iter_fragments():
                iter3 = self.model.append(iter2)
                self.model.set(iter3, 1, str(frag))

                for atm in frag.iter_all_atoms():
                    iter4 = self.model.append(iter3)
                    self.model.set(iter4, 1, str(atm))
                    
    def append_struct(self, struct):
        self.struct_list.append(struct)
        self.display_struct(struct)

    def remove_struct(self, struct):
        self.struct_list.remove(struct)
        self.redraw()



class StructureGUI:
    def __init__(self):
        self.struct_list = []
    
        ## MAIN WIDGET
        self.hpaned = gtk.HPaned()
        self.hpaned.set_border_width(2)

        ## LEFT HALF: StructureTreeControl
        self.sw1 = gtk.ScrolledWindow()
        self.hpaned.add1(self.sw1)
        self.sw1.set_border_width(3)
        self.sw1.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        
        self.struct_ctrl = StructureTreeControl(self)
        self.sw1.add(self.struct_ctrl)

        ## RIGHT HALF: GtkGLViewer
        self.gtkglviewer = GtkGLViewer()
        self.hpaned.add2(self.gtkglviewer)

    def getWidget(self):
        return self.hpaned

    def load_structure(self, path, update_cb):
        struct = LoadStructure(
            fil              = path,
            update_cb        = update_cb,
            build_properties = ("sequence","bonds"))

        self.struct_ctrl.append_struct(struct)
        self.gtkglviewer.append_struct(struct)


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


        ## Create document
        self.structure_gui = StructureGUI()
        table.attach(self.structure_gui.getWidget(),
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
        self.hbox.pack_start(self.statusbar, gtk.TRUE, gtk.TRUE, 0)
        self.statusbar.set_has_resize_grip(gtk.FALSE)

        self.frame = gtk.Frame()
        self.hbox.pack_start(self.frame, gtk.FALSE, gtk.FALSE, 0)
        self.frame.set_border_width(3)

        self.progress = gtk.ProgressBar()
        self.frame.add(self.progress)
        self.progress.set_size_request(100, 0)

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
        title = "mmView: " + title
        title = title[:50]
        self.window.set_title(title)

    def setStatusBar(self, text):
        self.statusbar.pop(0)
        self.statusbar.push(0, text)

    def load_file(self, path):
        """Loads the structure file specified in the path.
        """
        self.setTitle(path)
        self.setStatusBar("Loading: %s" % (path))
        self.structure_gui.load_structure(path, self.update_cb)
        self.setStatusBar("")
        self.progress.set_fraction(0.0)
        
        return gtk.FALSE
    
    def update_cb(self, percent):
        """Callback for file loading code to inform the GUI of how
        of the file has been read
        """
        self.progress.set_fraction(percent/100.0)
        while gtk.events_pending():
            gtk.main_iteration(gtk.TRUE)


def main(path = None):
    main_window_list = []

    def quit_notify_cb(window, mw):
        main_window_list.remove(mw)
        if len(main_window_list) < 1:
            gtk.main_quit()

    mw = MainWindow(quit_notify_cb)
    main_window_list.append(mw)

    if path:
        gobject.timeout_add(25, mw.load_file, path)

    gtk.main()


if __name__ == "__main__":
    import sys

    try:
        main(sys.argv[1])
    except IndexError:
        main()
