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


class SelectionGLAtomDrawList(GLAtomDrawList):
    def __init__(self):
        GLAtomDrawList.__init__(self)
        self.selected = []
    
    def select_atom_material(self, atm):
        (r, g, b, brightness) = GLAtomDrawList.select_atom_material(self, atm)
        if not atm in self.selected:
            brightness = brightness * 0.2
        return (r, g, b, brightness)
        

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

        self.glviewer = None

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

            self.glviewer.zpos = max(self.glviewer.zpos, -250.0)

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

    def set_unit_cell(self, unit_cell):
        self.glviewer.append(GLUnitCellDrawList(unit_cell))
        self.queue_draw()

    def set_structure(self, struct_obj, sel_struct_obj = None):
        ## remove all previous atom draw lists
        for draw_list in self.glviewer:
            if isinstance(draw_list, GLAtomDrawList):
                self.glviewer.remove(draw_list)
                draw_list.gl_delete_list()

        ## create and add new view group
        draw_list = SelectionGLAtomDrawList()
        self.glviewer.append(draw_list)

        for atm in struct_obj.iter_atoms():
            draw_list.append(atm)

        if sel_struct_obj:
            if isinstance(sel_struct_obj, Atom):
                draw_list.selected=[sel_struct_obj]
            else:
                draw_list.selected=[atm for atm in sel_struct_obj.iter_atoms()]

        self.queue_draw()


class StructurePanel(gtk.VPaned):
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

        self.glviewer.set_unit_cell(self.struct.unit_cell)
        self.glviewer.set_structure(self.struct, sel_struct_obj)


class StructureTreeModel(gtk.GenericTreeModel):
    def __init__(self, structure):
        self.structure_list = [structure]
	gtk.GenericTreeModel.__init__(self)

    def get_iter_root(self):
        return self.structure_list

    def on_get_flags(self):
	'''returns the GtkTreeModelFlags for this particular type of model'''
	return 0

    def on_get_n_columns(self):
	'''returns the number of columns in the model'''
	return 1

    def on_get_column_type(self, index):
	'''returns the type of a column in the model'''
	return gobject.TYPE_STRING

    def on_get_path(self, node):
	'''returns the tree path (a tuple of indices at the various
	levels) for a particular node.
        '''
        if isinstance(node, Structure):
            return ( self.structure_list.index(node), )

        elif isinstance(node, Chain):
            struct_path = self.structure_list.index(node.get_structure())
            chain_path = node.get_structure().index(node)
            return (struct_path, chain_path )

        elif isinstance(node, Fragment):
            struct_path = self.structure_list.index(node.get_structure())
            chain_path  = node.get_structure().index(node.get_chain())
            frag_path   = node.get_chain().index(node)
            return (structure_path, chain_path, frag_path)

        elif isinstance(node, Atom):
            struct_path = self.structure_list.index(node.get_structure())
            chain_path  = node.get_structure().index(node.get_chain())
            frag_path   = node.get_chain().index(node.get_fragment)
            atom_path   = node.get_fragment().index(node)
            return (struct_path, chain_path, frag_path, atom_path)

    def on_get_iter(self, path):
        '''returns the node corresponding to the given path.'''
        struct = self.structure_list[path[0]]
        if len(path) == 1:
            return struct

        chain = struct[path[1]]
        if len(path) == 2:
            return chain

        frag = chain[path[2]]
        if len(path) == 3:
            return frag

        atm = frag[path[3]]
        if len(path) == 4:
            return atm
        
    def on_get_value(self, node, column):
	'''returns the value stored in a particular column for the node'''
	return str(node)

    def on_iter_next(self, node):
	'''returns the next node at this level of the tree'''
        if isinstance(node, Structure):
            i = self.structure_list.index(node) + 1
            try:
                return self.structure_list[i]
            except IndexError:
                return None

        elif isinstance(node, Chain):
            struct = node.get_structure()
            i      = struct.index(node)
            try:
                return struct[i+1]
            except IndexError:
                return None

        elif isinstance(node, Fragment):
            return node.get_offset_fragment(1)
            
        elif isinstance(node, Atom):
            struct = node.get_structure()
            frag   = struct[node.chain_id][node.fragment_id]
            i      = frag.index(node)
            try:
                return frag[i+1]
            except IndexError:
                return None

    def on_iter_children(self, node):
	'''returns the first child of this node'''
        if node == self.structure_list:
            try:
                return self.structure_list[0]
            except IndexError:
                return None

        elif isinstance(node, Structure)  or \
             isinstance(node, Chain)      or \
             isinstance(node, Fragment):
            try:
                return node[0]
            except IndexError:
                return None

        else:
            return None

    def on_iter_has_child(self, node):
	'''returns true if this node has children'''
        if   node == self.structure_list  or \
             isinstance(node, Structure)  or \
             isinstance(node, Chain)      or \
             isinstance(node, Fragment):
            return len(node)
        else:
            return False

    def on_iter_n_children(self, node):
	'''returns the number of children of this node'''
        if   node == self.structure_list  or \
             isinstance(node, Structure)  or \
             isinstance(node, Chain)      or \
             isinstance(node, Fragment):
            print "asdfasdf",len(node)
            return len(node)
        else:
            return 0

    def on_iter_nth_child(self, node, n):
	'''returns the nth child of this node'''
        if   node == self.structure_list  or \
             isinstance(node, Structure)  or \
             isinstance(node, Chain)      or \
             isinstance(node, Fragment):
            return node[n]
        else:
            return None

    def on_iter_parent(self, node):
	'''returns the parent of this node'''
        if node == self.structure_list:
            return None

        elif isinstance(node, Structure):
            return self.structure_list

        elif isinstance(node, Chain):
            return node.get_structure()

        elif isinstance(node, Fragment):
            return node.get_chain()

        elif isinstance(node, Atom):
            return node.get_fragment()


class StructureGUI:
    def __init__(self):
        self.struct = None
    
        ## MAIN WIDGET
        self.hpaned = gtk.HPaned()
        self.hpaned.set_border_width(2)

        ## LEFT HALF
        scrolled_window = gtk.ScrolledWindow()
        scrolled_window.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        self.hpaned.add1(scrolled_window)

        self.struct_tree_view = gtk.TreeView()
        scrolled_window.add(self.struct_tree_view)
        cell = gtk.CellRendererText()
        column = gtk.TreeViewColumn("Structure", cell, text=0)
        self.struct_tree_view.append_column(column)
    
        self.struct_tree_view.connect("row_activated", self.row_activated)

        ## RIGHT HALF
        self.struct_panel = StructurePanel()
        self.hpaned.add2(self.struct_panel)

    def getWidget(self):
        return self.hpaned

    def load_structure(self, path):
        self.struct = LoadStructure(
            fil              = path,
            build_properties = ("sequence","bonds"))

        model = StructureTreeModel(self.struct)
        self.struct_tree_view.set_model(model)
        self.struct_panel.set_structure(self.struct)
        

    def row_activated(self, tree_view, path, column):
        ## find selected node
        model = tree_view.get_model()
        node = model.on_get_iter(path)
        self.struct_panel.set_structure(self.struct, node)


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
        title = "mmView: " + title
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

        self.structure_gui.load_structure(path)

        self.setStatusBar("")
        return gtk.FALSE
    



    
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
