#!/usr/bin/env python
## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import pygtk
pygtk.require("2.0")

import gobject
import gtk
import gtk.gtkgl

from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

from Scientific.Geometry import Vector

from mmLib.Structure  import *
from mmLib.FileLoader import LoadStructure, SaveStructure
from mmLib.CCP4       import CCP4MonomerLibrary


try:
    # try double-buffered
    glconfig = gtk.gdkgl.Config(mode = gtk.gdkgl.MODE_RGB    |
                                gtk.gdkgl.MODE_DOUBLE |
                                gtk.gdkgl.MODE_DEPTH)
except gtk.gdkgl.NoMatches:
    # try single-buffered
    glconfig = gtk.gdkgl.Config(mode = gtk.gdkgl.MODE_RGB    |
                                gtk.gdkgl.MODE_DEPTH)



class GLViewer:
    def __init__(self, atom_container):
        self.atom_container = atom_container
        self.calcObjects()
        
        self.rot  = "x"
        self.rotx = 0
        self.roty = 0
        self.is_solid = gtk.TRUE

        ## GL Widget
        self.glarea = gtk.gtkgl.DrawingArea()
        self.glarea.set_size_request(300, 300)

        gtk.gtkgl.widget_set_gl_capability(self.glarea, glconfig)

        self.glarea.connect('realize', self.realize)
        self.glarea.connect('map', self.map)
        self.glarea.connect('unmap', self.unmap)
        self.glarea.connect('configure_event', self.configure_event)
        self.glarea.connect('expose_event', self.expose_event)
        self.glarea.connect('destroy', self.glarea_destroy)
        self.glarea.show()

        ## put the OpenGL widget in a pretty frame
        self.frame = gtk.Frame()
        self.frame.set_shadow_type(gtk.SHADOW_ETCHED_OUT)
        self.frame.add(self.glarea)
        self.frame.connect('destroy', self.destroy)

        
    def __del__(self):
        print "GLViewer Destroyed"


    def getWidget(self):
        return self.frame

    def destroy(self, frame):
        gtk.timeout_remove(self.tn)
        return gtk.TRUE

    def glarea_destroy(self, glarea):
        print "glarea_destroy"
        return gtk.TRUE

    def timeout(self):
        print "timeout"
        if self.rot == "x":
            self.rotx += 1
            if self.rotx > 360:
                self.rotx = 0
                self.rot = "y"
        elif self.rot == "y":
            self.roty += 1
            if self.roty > 360:
                self.roty = 0
                self.rot = "x"
        self.glarea.queue_draw()
        
        return gtk.TRUE

    def map(self, glarea):
        self.tn = gtk.timeout_add(80, self.timeout)
        print "timeout: ", self.tn
        return gtk.TRUE

    def unmap(self, glarea):
        if hasattr(self, "tn"):
            gtk.timeout_remove(self.tn)
            print "removing timeout ",self.tn
        else:
            print "no self.tn"
        return gtk.TRUE

    def realize(self, glarea):
        # get GLContext and GLDrawable
	glcontext = self.glarea.get_gl_context()
	gldrawable = self.glarea.get_gl_drawable()
	
	# GL calls
	if not gldrawable.gl_begin(glcontext):
            return gtk.FALSE
	
	glMaterial(GL_FRONT, GL_AMBIENT,   [0.2, 0.2, 0.2, 1.0])
	glMaterial(GL_FRONT, GL_DIFFUSE,   [0.8, 0.8, 0.8, 1.0])
	glMaterial(GL_FRONT, GL_SPECULAR,  [1.0, 0.0, 1.0, 1.0])
	glMaterial(GL_FRONT, GL_SHININESS, 50.0)
	
	glLight(GL_LIGHT0, GL_AMBIENT,  [0.0, 1.0, 0.0, 1.0])
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
	glcontext = self.glarea.get_gl_context()
	gldrawable = self.glarea.get_gl_drawable()
	
	# GL calls
	if not gldrawable.gl_begin(glcontext):
            return gtk.FALSE
	
	x, y, width, height = glarea.get_allocation()
	
	glViewport(0, 0, width, height)
	glMatrixMode(GL_PROJECTION)
	glLoadIdentity()
	if width > height:
            w = float(width) / float(height)
            glFrustum(-w, w, -1.0, 1.0, 3.0, 100.0)
	else:
            h = float(height) / float(width)
            glFrustum(-1.0, 1.0, -h, h, 3.0, 100.0)
	
	glMatrixMode(GL_MODELVIEW)
	
	gldrawable.gl_end()
	
	return gtk.TRUE

    def expose_event(self, glarea, event):
	# get GLContext and GLDrawable
	glcontext = self.glarea.get_gl_context()
	gldrawable = self.glarea.get_gl_drawable()
	
	# GL calls
	if not gldrawable.gl_begin(glcontext):
            return gtk.FALSE
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
	
	glLoadIdentity()

	glTranslate(0, 0, -10)
	
	glRotate(self.rotx, 1, 0, 0)
	glRotate(self.roty, 0, 1, 0)

        ## draw the atoms
        for (coords, size) in self.sphere_list:
            glPushMatrix()
            apply(glTranslatef, coords)
            glutSolidSphere(size, 16, 16)
            glPopMatrix()

        ## draw the bonds between the atoms
        for (c1, c2, lwidth) in self.line_list:
            glLineWidth(lwidth)
            glBegin(GL_LINES)
            apply(glVertex3f, c1)
            apply(glVertex3f, c2)
            glEnd()
            
	if gldrawable.is_double_buffered():
            gldrawable.swap_buffers()
	else:
            glFlush()
	
	gldrawable.gl_end()
	
	return gtk.TRUE

    def calcObjects(self):
        self.sphere_list = []
        self.line_list = []
        
        ## compute the centroid
        centroid = Vector(0.0, 0.0, 0.0)
        n = 0
        for atm in self.atom_container.atomIterator():
            n += 1
            centroid += atm.getPosition()
        centroid = centroid / n        

        ## atom spheres
        for atm in self.atom_container.atomIterator():
            pos = atm.getPosition() - centroid
            size = atm.el.atomic_number / 10.0
            self.sphere_list.append( ((pos[0], pos[1], pos[2]), size) )

        for bond in self.atom_container.bondIterator():
            ## calculate line width based on the number of standard deviations
            ## from ideal
            (dist, dist_esd) = bond.getDistance()
            if dist and dist_esd:
                sigma = (dist - bond.calcLength()) / dist_esd
                lwidth = 10.0 - (3.0 * sigma)
                if lwidth < 1.0:
                    lwidth = 1.0
            else:
                lwidth = 5.0

            (atm1, atm2) = bond.getAtoms()
            v1 = atm1.getPosition() - centroid
            v2 = atm2.getPosition() - centroid
            self.line_list.append(
               ( (v1[0], v1[1], v1[2]), (v2[0], v2[1], v2[2]), lwidth) )



class DefaultDisplay:
    def __init__(self, chain):
        self.chain = chain

        ## create the GLViewer
        glviewer = GLViewer(chain)

        ## create the frame
        self.frame = gtk.Frame()
        self.frame.set_border_width(3)
        self.frame.add(glviewer.getWidget())
        
    def getWidget(self):
        return self.frame



class AAResidueDisplay:
    def __init__(self, aa_res):
        self.aa_res = aa_res

        ## create the text box
        vbox = gtk.VBox()

        def add_string(s):
            label = gtk.Label(s)
            vbox.add(label)

        add_string("Residue Name: %s" % (self.aa_res.getName()))
        add_string("Sequence ID: %s" % (self.aa_res.getSequenceID()))

        psi = self.aa_res.calcTorsionPsi()
        if psi:
            add_string("PSI Torsion: %f" % (psi))

        phi = self.aa_res.calcTorsionPhi()
        if phi:
            add_string("PHI Torsion: %f" % (phi))


        frame = gtk.Frame()
        frame.set_border_width(3)
        frame.add(vbox)

        ## create the GLViewer
        glviewer = GLViewer(aa_res)

        ## stuff the text box and glviewer in a vbox
        self.vbox = gtk.VBox()
        self.vbox.add(frame)
        self.vbox.add(glviewer.getWidget())
        
    def getWidget(self):
        return self.vbox



class StructureTreeModel(gtk.GenericTreeModel):
    '''This class represents the model of a tree.  The iterators used
    to represent positions are converted to python objects when passed
    to the on_* methods.  This means you can use any python object to
    represent a node in the tree.  The None object represents a NULL
    iterator.

    In this tree, we use simple tuples to represent nodes, which also
    happen to be the tree paths for those nodes.  This model is a tree
    of depth 3 with 5 nodes at each level of the tree.  The values in
    the tree are just the string representations of the nodes.'''

    def __init__(self, structure):
	'''constructor for the model.  Make sure you call
	PyTreeModel.__init__'''
        self.structure = structure
	gtk.GenericTreeModel.__init__(self)

    def get_iter_root(self):
        return self.structure

    # the implementations for TreeModel methods are prefixed with on_
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
	return node.getIndexPath()
    def on_get_iter(self, path):
        '''returns the node corresponding to the given path.'''
        node = self.structure
        for i in path:
            node = node.getChild(i)
        return node
    def on_get_value(self, node, column):
	'''returns the value stored in a particular column for the node'''
	return str(node)
    def on_iter_next(self, node):
	'''returns the next node at this level of the tree'''
        try:
            return node.getSibling(1)
        except IndexError:
            return None
    def on_iter_children(self, node):
	'''returns the first child of this node'''
        return node.getFirstChild()
    def on_iter_has_child(self, node):
	'''returns true if this node has children'''
        if node.getDegree() > 0:
            return 1
        return 0
    def on_iter_n_children(self, node):
	'''returns the number of children of this node'''
        return node.getDegree()
    def on_iter_nth_child(self, node, n):
	'''returns the nth child of this node'''
        if node == None:
            return None
        try:
            return node.getChild(n)
        except IndexError:
            return None
    def on_iter_parent(self, node):
	'''returns the parent of this node'''
        return node.getParent()



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
        self.frame = gtk.Frame()
        self.hpaned.add2(self.frame)
        
        self.data_widget = None
        

    def getWidget(self):
        return self.hpaned


    def loadStructure(self):
        monomer_lib = CCP4MonomerLibrary()
        self.structure = LoadStructure(sys.argv[1])
        for aa_res in self.structure.aminoAcidResidueIterator():
            monomer_lib.setResidueInfo(aa_res)
        
        model = StructureTreeModel(self.structure)
        self.struct_tree_view.set_model(model)


    def row_activated(self, tree_view, path, column):
        ## remove old child
        if self.data_widget:
            self.data_widget.destroy()
            self.data_widget = None
        
        ## find selected node
        node = self.structure
        for i in path:
            node = node.getChild(i)

        ## select display widget
        dis = None
        if isinstance(node, AminoAcidResidue):
            dis = AAResidueDisplay(node)
        else:
            dis = DefaultDisplay(node)
            
        ## add new child, and show
        if dis:
            self.data_widget = dis.getWidget()
            self.frame.add(self.data_widget)
            self.data_widget.show_all()



class MainWindow:
    def __init__(self, quit_notify_cb):
        self.quit_notify_cb = quit_notify_cb

        menu_items = (
            ## FILE Menu
            ('/_File',            None,         None,        0, '<Branch>' ),
            ('/File/_New',        '<control>N',
             self.menuitem_cb, 0, '<StockItem>', gtk.STOCK_NEW),
            ('/File/_Open',       '<control>O',
             self.menuitem_cb, 0, '<StockItem>', gtk.STOCK_OPEN),
            ('/File/_Save',       '<control>S',
             self.menuitem_cb, 0, '<StockItem>', gtk.STOCK_SAVE),
            ('/File/Save _As...', None,
             self.menuitem_cb, 0, '<StockItem>', gtk.STOCK_SAVE),
            ('/File/sep1',        None,
             self.menuitem_cb, 0, '<Separator>'),
            ('/File/_Quit',       '<control>Q',
             self.menuitem_cb, 0, '<StockItem>', gtk.STOCK_QUIT),

            ## HELP Menu
            ('/_Help',       None, None, 0, '<Branch>'),
            ('/Help/_About', None, self.menuitem_cb, 0, ''))


        ## Create the toplevel window

        self.window = gtk.Window()
        self.window.set_title('mmView')
        self.window.set_default_size(200, 200)
        self.window.connect('destroy', self.quit_notify_cb, self)

        table = gtk.Table(1, 4, gtk.FALSE)
        self.window.add(table)

        ## Create the menubar

        accel_group = gtk.AccelGroup()
        self.window.add_accel_group(accel_group)

        item_factory = gtk.ItemFactory(gtk.MenuBar, '<main>', accel_group)

        ## create menu items

        item_factory.create_items(menu_items, self.window)

        table.attach(item_factory.get_widget('<main>'),
                     # X direction           Y direction
                     0, 1,                      0, 1,
                     gtk.EXPAND | gtk.FILL,     0,
                     0,                         0)

        ## Create the toolbar

        toolbar = gtk.Toolbar()
        toolbar.insert_stock(gtk.STOCK_OPEN,
                             "This is a demo button with an 'open' icon",
                             None,
                             self.toolbar_cb,
                             self.window,
                             -1)

        toolbar.insert_stock(gtk.STOCK_CLOSE,
                             "This is a demo button with an 'close' icon",
                             None,
                             self.toolbar_cb,
                             self.window,
                             -1)

        table.attach(toolbar,
                     # X direction           Y direction
                     0, 1,                   1, 2,
                     gtk.EXPAND | gtk.FILL,  0,
                     0,                      0)

        ## Create document

        self.structure_gui = StructureGUI()
        self.structure_gui.loadStructure()
        table.attach(self.structure_gui.getWidget(),
                     # X direction           Y direction
                     0, 1,                   2, 3,
                     gtk.EXPAND | gtk.FILL,  gtk.EXPAND | gtk.FILL,
                     0,                      0)


        ## Create statusbar 

        statusbar = gtk.Statusbar();
        table.attach(statusbar,
                     # X direction           Y direction
                     0, 1,                   3, 4,
                     gtk.EXPAND | gtk.FILL,  0,
                     0,                      0)

##         self.status_buffer = gtk.TextBuffer()
##         self.status_buffer.connect('changed', self.update_statusbar, statusbar)
##         self.status_buffer.connect('mark_set', self.mark_set_callback, statusbar)
##         update_statusbar(self.status_buffer, statusbar)

        self.window.show_all()
 

    def menuitem_cb(self, widget):
        pass


    def toolbar_cb(self, widget):
        pass

    def update_statusbar_cb(self, widget):
        pass

    def mark_set_callback(self, widget):
        pass



def main():
    main_window_list = []

    def quit_notify_cb(window, mw):
        main_window_list.remove(mw)
        if len(main_window_list) < 1:
            gtk.main_quit()
    
    main_window_list.append(MainWindow(quit_notify_cb))
    gtk.main()



if __name__ == "__main__":
    import sys

    try:
        path = sys.argv[1]
    except IndexError:
        print "usage %s <mmCIF file path>" % (sys.argv[0])
        sys.exit(1)

    main()
