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


class GLViewer(gtk.gtkgl.DrawingArea):
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

        self.rotx        = 0
        self.roty        = 0

        self.xpos        = 0
        self.ypos        = 0
        self.zpos        = -50

        self.draw_list   = []
        
    def destroy(self, glarea):
        return gtk.TRUE

    def button_press_event(self, glarea, event):
        self.beginx = event.x
        self.beginy = event.y

    def motion_notify_event(self, glarea, event):
        width     = glarea.allocation.width
        height    = glarea.allocation.height

        if (event.state & gtk.gdk.BUTTON1_MASK):
            self.roty = self.roty + \
                        (((event.x - self.beginx) / float(width)) * 360.0)
            self.rotx = self.rotx + \
                        (((event.y - self.beginy) / float(height)) * 360.0)

        elif (event.state & gtk.gdk.BUTTON2_MASK):
            self.zpos = self.zpos + \
                        (((event.y - self.beginy) / float(height)) * 50.0)

            self.zpos = max(self.zpos, -250)

        elif (event.state & gtk.gdk.BUTTON3_MASK):
            self.ypos = self.ypos - \
                        (((event.y - self.beginy) / float(height)) * 50.0)
            self.xpos = self.xpos + \
                        (((event.x - self.beginx) / float(width)) * 50.0)

        self.beginx = event.x
        self.beginy = event.y
        self.queue_draw()

    def map(self, glarea):
        self.add_events(gtk.gdk.BUTTON_PRESS_MASK   |
                        gtk.gdk.BUTTON_RELEASE_MASK |
                        gtk.gdk.BUTTON_MOTION_MASK  |
                        gtk.gdk.POINTER_MOTION_MASK)

        return gtk.TRUE

    def unmap(self, glarea):
        return gtk.TRUE

    def set_material(self, r, g, b, bright = 1.0):
        
	glMaterial(GL_FRONT, GL_AMBIENT,
                   [0.2*r*bright, 0.2*g*bright, 0.2*b*bright, 1.0])
	glMaterial(GL_FRONT, GL_DIFFUSE,
                   [0.8*bright, 0.8*bright, 0.8*bright, 1.0])
	glMaterial(GL_FRONT, GL_SPECULAR,
                   [1.0*bright, 1.0*bright, 1.0*bright, 1.0])
        
	glMaterial(GL_FRONT, GL_SHININESS, 50.0*bright)

    def realize(self, glarea):
        # get GLContext and GLDrawable
	glcontext = self.get_gl_context()
	gldrawable = self.get_gl_drawable()
	
	# GL calls
	if not gldrawable.gl_begin(glcontext):
            return gtk.FALSE

        self.set_material(1.0, 1.0, 1.0)
	
	glLight(GL_LIGHT0, GL_AMBIENT,  [1.0, 1.0, 1.0, 1.0])
	glLight(GL_LIGHT0, GL_DIFFUSE,  [1.0, 1.0, 1.0, 1.0])
	glLight(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
	glLight(GL_LIGHT0, GL_POSITION, [1.0, 1.0, 1.0, 0.0])
	
	glLightModel(GL_LIGHT_MODEL_AMBIENT, [0.2, 0.2, 0.2, 1.0])

 	glEnable(GL_LIGHTING)
 	glEnable(GL_LIGHT0)
	
	glDepthFunc(GL_LESS)
	glEnable(GL_DEPTH_TEST)
		
        gldrawable.gl_end()
        return gtk.TRUE

    def configure_event(self, glarea, event):
	# get GLContext and GLDrawable
	glcontext = self.get_gl_context()
	gldrawable = self.get_gl_drawable()
	
	# GL calls
	if not gldrawable.gl_begin(glcontext):
            return gtk.FALSE
	
	x, y, width, height = glarea.get_allocation()
	
	glViewport(0, 0, width, height)
	glMatrixMode(GL_PROJECTION)
	glLoadIdentity()
	if width > height:
            w = float(width) / float(height)
            glFrustum(-w, w, -1.0, 1.0, 3.0, 150.0)
	else:
            h = float(height) / float(width)
            glFrustum(-1.0, 1.0, -h, h, 3.0, 150.0)
	
	glMatrixMode(GL_MODELVIEW)
	
	gldrawable.gl_end()
	
	return gtk.TRUE

    def expose_event(self, glarea, event):
	# get GLContext and GLDrawable
	glcontext = self.get_gl_context()
	gldrawable = self.get_gl_drawable()
	
	# GL calls
	if not gldrawable.gl_begin(glcontext):
            return gtk.FALSE
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
	
	glLoadIdentity()

	glTranslate(self.xpos, self.ypos, self.zpos)
	
	glRotate(self.rotx, 1, 0, 0)
	glRotate(self.roty, 0, 1, 0)

        ## <DRAW>
        for i in self.draw_list: glCallList(i)
        ## </DRAW>
            
	if gldrawable.is_double_buffered():
            gldrawable.swap_buffers()
	else:
            glFlush()
	
	gldrawable.gl_end()
	
	return gtk.TRUE

    def setAtomContainer(self, atom_container, selected = None):
        ## compute the centroid
        if not hasattr(self, "centroid"):
            self.centroid = Vector(0.0, 0.0, 0.0)
            n = 0
            for atm in atom_container.iterAtoms():
                n += 1
                self.centroid += atm.position
            self.centroid = self.centroid / n        

        ## selected map
        select_list = []

        for atm in selected.iterAtoms():
            if atm not in select_list:
                select_list.append(atm)

        ## Temp Factor
        max_tf = 25.0

        ## atom spheres
        visited_bonds = []

        glNewList(1, GL_COMPILE)
        if not self.draw_list: self.draw_list.append(1)

        for atm in atom_container.iterAtoms():
            pos  = atm.position - self.centroid

            tf = math.log10(atm.temp_factor) / max_tf

            ## hilight selections
            try:
                select_list.remove(atm)
            except ValueError:
                br = 0.2
            else:
                br = 1.0

            size = 0.3

            glPushMatrix()
            glTranslatef(pos[0], pos[1], pos[2])

            if   atm.element == "C": self.set_material(1.0, 1.0, 0.0, br)
            elif atm.element == "N": self.set_material(0.0, 0.0, 1.0, br)
            elif atm.element == "O": self.set_material(1.0, 0.0, 0.0, br)
            elif atm.element == "S": self.set_material(0.0, 1.0, 0.0, br)
            else:                    self.set_material(1.0, 1.0, 1.0, 1.0)

            glutSolidSphere(size, 12, 12)
            glPopMatrix()

            ## bonds
            glLineWidth(5.0)

            for bond in atm.iterBonds():
                if bond in visited_bonds: continue
                else:                     visited_bonds.insert(0, bond)

                atm2 = bond.getPartner(atm)
                v1   = atm.position - self.centroid
                v2   = atm2.position - self.centroid

                self.set_material(1.0, 1.0, 1.0, 0.2)

                glBegin(GL_LINES)
                glVertex3f(v1[0], v1[1], v1[2])
                glVertex3f(v2[0], v2[1], v2[2])
                glEnd()
                
        glEndList()

        self.queue_draw()



class AtomContainerPanel(gtk.VPaned):
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
        self.glviewer = GLViewer()
        self.add2(self.glviewer)

    def printBox(self, key, value):
        iter = self.store.append()
        self.store.set(iter, 0, key, 1, str(value))

    def setAtomContainer(self, atom_container, selected = None):
        self.atom_container = atom_container
        self.store.clear()

        if isinstance(atom_container, Residue):
            res = atom_container
            self.printBox("Residue.res_name", res.res_name)
            self.printBox("Residue.getOffsetResidue(-1)",
                          res.getOffsetResidue(-1))
            self.printBox("Residue.getOffsetResidue(1)",
                          res.getOffsetResidue(1))

        if isinstance(atom_container, Fragment):
            frag = atom_container
            bonds = ""
            for bond in frag.iterBonds():
                bonds += str(bond)
            self.printBox("Fragment.iterBonds()", bonds)

        if isinstance(atom_container, AminoAcidResidue):
            aa_res = atom_container
            self.printBox("AminoAcidResidue.calcMainchainBondLength()",
                          aa_res.calcMainchainBondLength())
            self.printBox("AminoAcidResidue.calcMainchainBondAngle()",
                          aa_res.calcMainchainBondAngle())
            self.printBox("AminoAcidResidue.calcTorsionPsi()",
                          aa_res.calcTorsionPsi())
            self.printBox("AminoAcidResidue.calcTorsionPhi()",
                          aa_res.calcTorsionPhi())
            self.printBox("AminoAcidResidue.calcTorsionOmega()",
                          aa_res.calcTorsionOmega())
            self.printBox("AminoAcidResidue.calcTorsionChi()",
                          aa_res.calcTorsionChi())

        if isinstance(atom_container, Atom):
            atm = atom_container
            self.printBox("Atom.element", atm.element)
            self.printBox("Atom.name", atm.name)
            self.printBox("Atom.occupancy", atm.occupancy)
            self.printBox("Atom.temp_factor", atm.temp_factor)
            self.printBox("Atom.U", atm.U)
            self.printBox("Atom.position", atm.position)
            self.printBox("len(Atom.bond_list)", len(atm.bond_list))
            self.printBox("Atom.calcAnisotropy()", atm.calcAnisotropy())

        self.glviewer.setAtomContainer(self.atom_container, selected)



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
	levels) for a particular node.'''

        if isinstance(node, Structure):
            return ( self.structure_list.index(node), )

        elif isinstance(node, Chain):
            struct_path = self.structure_list.index(node.getStructure())
            chain_path = node.getStructure().index(node)
            return (struct_path, chain_path )

        elif isinstance(node, Fragment):
            struct_path = self.structure_list.index(node.getStructure())
            chain_path  = node.getStructure().index(node.getChain())
            frag_path   = node.getChain().index(node)
            return (structure_path, chain_path, frag_path)

        elif isinstance(node, Atom):
            struct_path = self.structure_list.index(node.getStructure())
            chain_path  = node.getStructure().index(node.getChain())
            frag_path   = node.getChain().index(node.getFragment)
            atom_path   = node.getFragment().index(node)
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
            struct = node.getStructure()
            i      = struct.index(node)
            try:
                return struct[i+1]
            except IndexError:
                return None

        elif isinstance(node, Fragment):
            return node.getOffsetFragment(1)
            
        elif isinstance(node, Atom):
            struct = node.getStructure()
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
            return node.getStructure()

        elif isinstance(node, Fragment):
            return node.getChain()

        elif isinstance(node, Atom):
            return node.getFragment()


class StructureGUI:
    def __init__(self):
        self.structure = None
    
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
        self.ac_panel = AtomContainerPanel()
        self.hpaned.add2(self.ac_panel)


    def getWidget(self):
        return self.hpaned

    def loadStructure(self, path):
        self.structure = LoadStructure(
            fil              = path,
            build_properties = ("sequence","bonds"))

        model = StructureTreeModel(self.structure)
        self.struct_tree_view.set_model(model)

    def row_activated(self, tree_view, path, column):
        ## find selected node
        model = tree_view.get_model()
        node = model.on_get_iter(path)
        self.ac_panel.setAtomContainer(self.structure, node)


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
        self.loadFile(path)
        

    def open_cancel_cb(self, cancel_button, file_selector):
        file_selector.destroy()


    def setTitle(self, title):
        title = "mmView: " + title
        title = title[:50]
        self.window.set_title(title)


    def setStatusBar(self, text):
        self.statusbar.pop(0)
        self.statusbar.push(0, text)


    def loadFile(self, path):
        self.setTitle(path)
        self.setStatusBar("Loading %s, please wait..." % (path))

        while gtk.events_pending():
            gtk.main_iteration(gtk.TRUE)

        self.structure_gui.loadStructure(path)

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
        gobject.timeout_add(25, mw.loadFile, path)

    gtk.main()


if __name__ == "__main__":
    import sys

    try:
        main(sys.argv[1])
    except IndexError:
        main()
