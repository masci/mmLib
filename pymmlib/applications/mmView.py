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

from mmLib.Structure      import *
from mmLib.FileLoader     import LoadStructure, SaveStructure
from mmLib.GLViewer       import *
from mmLib.Extensions.TLS import *

try:
    # try double-bufferedq
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

        alist = gl_struct.get_atom_list("a")

        for atm in struct.iter_atoms():
            if atm.name in ["CA", "C", "N"]:
                alist.append(atm)

        self.glviewer.append(gl_struct)

        #gobject.timeout_add(25, self.timeout_cb, alist)
        
        self.queue_draw()

        return gl_struct

    def remove_struct(self, struct):
        """Removes structure from the viewer.
        """
        pass

    def timeout_cb(self, atom_list):
        atom_list.rotx += 1
        self.queue_draw()
        return gtk.TRUE


class StructDetailsDialog(gtk.Dialog):
    def __init__(self, context):
        self.context = context
        
        gtk.Dialog.__init__(self,
                            "Selection Information",
                            self.context.window,
                            gtk.DIALOG_DESTROY_WITH_PARENT)

        self.add_button(gtk.STOCK_CLOSE, gtk.RESPONSE_CLOSE)
        self.set_default_size(400, 400)
        
        self.connect("destroy", self.response_cb)
        self.connect("response", self.response_cb)

        ## make the print box
        sw = gtk.ScrolledWindow()
        self.vbox.add(sw)
        sw.set_border_width(5)
        sw.set_shadow_type(gtk.SHADOW_ETCHED_IN)
        sw.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)


        self.store = gtk.ListStore(gobject.TYPE_STRING, gobject.TYPE_STRING)
     
        treeview = gtk.TreeView(self.store)
        sw.add(treeview)
        treeview.set_rules_hint(gtk.TRUE)
        treeview.set_search_column(0)

        column = gtk.TreeViewColumn(
            "mmLib.Structure Method",
            gtk.CellRendererText(),
            text=0)
        treeview.append_column(column)
    
        column = gtk.TreeViewColumn(
            "Return Value",
            gtk.CellRendererText(),
            text=1)
        treeview.append_column(column)

        self.show_all()

    def response_cb(self, *args):
        self.destroy()
    
    def add_line(self, key, value):
        iter = self.store.append()
        self.store.set(iter, 0, key, 1, str(value))

    def set_struct_obj(self, struct_obj):
        self.store.clear()

        if isinstance(struct_obj, Residue):
            self.add_line("Residue.res_name", struct_obj.res_name)
            self.add_line("Residue.get_offset_residue(-1)",
                          struct_obj.get_offset_residue(-1))
            self.add_line("Residue.get_offset_residue(1)",
                          struct_obj.get_offset_residue(1))

        if isinstance(struct_obj, Fragment):
            bonds = ""
            for bond in struct_obj.iter_bonds():
                bonds += str(bond)
            self.add_line("Fragment.iter_bonds()", bonds)

        if isinstance(struct_obj, AminoAcidResidue):
            self.add_line("AminoAcidResidue.calc_mainchain_bond_length()",
                          struct_obj.calc_mainchain_bond_length())
            self.add_line("AminoAcidResidue.calc_mainchain_bond_angle()",
                          struct_obj.calc_mainchain_bond_angle())
            self.add_line("AminoAcidResidue.calc_torsion_psi()",
                          struct_obj.calc_torsion_psi())
            self.add_line("AminoAcidResidue.calc_torsion_phi()",
                          struct_obj.calc_torsion_phi())
            self.add_line("AminoAcidResidue.calc_torsion_omega()",
                          struct_obj.calc_torsion_omega())
            self.add_line("AminoAcidResidue.calc_torsion_chi()",
                          struct_obj.calc_torsion_chi())

        if isinstance(struct_obj, Atom):
            self.add_line("Atom.element", struct_obj.element)
            self.add_line("Atom.name", struct_obj.name)
            self.add_line("Atom.occupancy", struct_obj.occupancy)
            self.add_line("Atom.temp_factor", struct_obj.temp_factor)
            self.add_line("Atom.U", struct_obj.U)
            self.add_line("Atom.position", struct_obj.position)
            self.add_line("len(Atom.bond_list)", len(struct_obj.bond_list))
            self.add_line("Atom.calc_anisotropy()",
                          struct_obj.calc_anisotropy())


class TLSDialog(gtk.Dialog):
    """Dialog for TLS analysis of a structure.
    """

    def __init__(self, **args):
        self.context        = args["context"]
        self.struct         = args["struct"]
        self.win_title      = args.get("title", str(self.struct))

        self.tls_group_list    = []
        self.gl_tls_group_dict = {}
        
        gtk.Dialog.__init__(self,
                            "TLS Analysis: %s" % (self.win_title),
                            self.context.window,
                            gtk.DIALOG_DESTROY_WITH_PARENT)

        self.add_button(gtk.STOCK_CLOSE, gtk.RESPONSE_CLOSE)
        self.set_default_size(400, 400)
        
        self.connect("destroy", self.response_cb)
        self.connect("response", self.response_cb)

        ## make the print box
        sw = gtk.ScrolledWindow()
        self.vbox.add(sw)
        sw.set_border_width(5)
        sw.set_shadow_type(gtk.SHADOW_ETCHED_IN)
        sw.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)

        self.store = gtk.ListStore(gobject.TYPE_STRING,
                                   gobject.TYPE_STRING,
                                   gobject.TYPE_STRING,
                                   gobject.TYPE_STRING)
     
        treeview = gtk.TreeView(self.store)
        sw.add(treeview)
        treeview.set_rules_hint(gtk.TRUE)
        treeview.set_search_column(0)

        column = gtk.TreeViewColumn(
            "TLS Group",
            gtk.CellRendererText(),
            text=0)
        treeview.append_column(column)
    
        column = gtk.TreeViewColumn(
            "Center of Reaction",
            gtk.CellRendererText(),
            text=1)
        treeview.append_column(column)

        column = gtk.TreeViewColumn(
            "Mean Translation (A^2)",
            gtk.CellRendererText(),
            text=2)
        treeview.append_column(column)

        column = gtk.TreeViewColumn(
            "Mean Libration (deg^2)",
            gtk.CellRendererText(),
            text=3)
        treeview.append_column(column)

        self.show_all()

        self.make_tls_group_segments()
        self.redraw()

    def response_cb(self, *args):
        self.destroy()
    
    def redraw(self):
        """Clear and redisplay TLS groups
        """
        self.store.clear()

        for tls in self.tls_group_list:
            iter = self.store.append()

            calcs = tls.calc_COR()

            trT = trace(tls.T)/3.0
            trL = trace(tls.L * rad2deg2)/3.0
            cor = calcs["COR"]

            self.store.set(iter,
                           0, tls.name,
                           1, "%.3f, %.3f, %.3f" % (cor[0], cor[1], cor[2]),
                           2, "%.3f" % (trT),
                           3, "%.3f" % (trL))

    def segment_list(self, chain, seg_len):
        segment = None
        segment_list = []

        for res in chain.iter_amino_acids():
            if segment == None:
                segment = [res]
                segment_list.append(segment)
            elif len(segment) >= seg_len:
                segment = [res]
                segment_list.append(segment)
            else:
                segment.append(res)

        return segment_list

    def running_segment_list(self, chain, seg_len):
        segment = []
        segment_list = []

        for res in chain.iter_amino_acids():
            segment.append(res)

            if len(segment) < seg_len:
                continue

            segment_list.append(segment)
            segment = segment[1:]

        return segment_list
        

    def make_tls_group_segments(self, seg_len = 5):
        """
        """
        self.tls_group_list = []
        
        for chain in self.struct.iter_chains():    
            for seg in self.running_segment_list(chain, seg_len):

                ## new tls group for segment
                tls = TLSGroup()
                self.tls_group_list.append(tls)
                tls.name = "%s-%s" % (seg[0], seg[-1])

                ## add segment atoms
                for res in seg:
                    for atm in res.iter_atoms():
                        tls.append(atm)

                ## avoid small groups
                if len(tls) < 4:
                    self.tls_group_list.remove(tls)
                    continue

                ## calculate tensors and print
                tls.origin = tls.calc_centroid()
                tls.calc_tls_tensors()

                self.show_tls_group(tls)

    def show_tls_group(self, tls):
        tls.gl_tls = GLTLSGroup(tls_group = tls)
        self.context.struct_gui.gtkglviewer.glviewer.append(tls.gl_tls)
        

class StructureTreeControl(gtk.TreeView):
    def __init__(self, context):
        self.context = context
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
        self.context.set_selected_item(self.resolv_path(path))
    
    def button_release_event_cb(self, tree_view, bevent):
        """
        """
        pass

    def resolv_path(self, path):
        """Get the item in the tree view at the specified path.
        """
        itemx = self.struct_list
        for i in path:
            itemx = itemx[i]
        return itemx
        
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
    def __init__(self, context):
        self.context = context        
    
        ## MAIN WIDGET
        self.hpaned = gtk.HPaned()
        self.hpaned.set_border_width(2)

        ## LEFT HALF: StructureTreeControl
        self.sw1 = gtk.ScrolledWindow()
        self.hpaned.add1(self.sw1)
        self.sw1.set_border_width(3)
        self.sw1.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        
        self.struct_ctrl = StructureTreeControl(self.context)
        self.sw1.add(self.struct_ctrl)

        ## RIGHT HALF: GtkGLViewer
        self.gtkglviewer = GtkGLViewer()
        self.hpaned.add2(self.gtkglviewer)

    def get_widget(self):
        return self.hpaned


class MainWindow:
    def __init__(self, quit_notify_cb):

        self.struct_list         = []
        self.struct_dict         = {}
        
        self.selected_struct_obj = None        
        self.quit_notify_cb      = quit_notify_cb

        ## dialog lists
        self.details_dialog_list = []
        self.tls_dialog_list     = []

        ## Create the toplevel window
        self.window = gtk.Window()
        self.set_title("")
        self.window.set_default_size(500, 400)
        self.window.connect('destroy', self.quit_cb, self)

        table = gtk.Table(1, 5, gtk.FALSE)
        self.window.add(table)

        ### MENUBAR
        MENU_ITEMS = (
            ('/_File',
             None,
             None,
             0,
             '<Branch>' ),

            ('/File/_New Window',
             '',
             self.new_window_cb,
             0,
             '<StockItem>',
             gtk.STOCK_NEW),

            ('/File/_New Tab',
             '',
             self.new_tab_cb,
             0,
             '<StockItem>',
             gtk.STOCK_NEW),
            
            ('/File/_Open',
             '',
             self.open_cb,
             0,
             '<StockItem>',
             gtk.STOCK_OPEN),
            
            ('/File/sep1',
             None,
             None,
             0,
             '<Separator>'),

            ('/File/_Quit',
             '<control>Q',
             self.quit_cb,
             0,
             '<StockItem>',
             gtk.STOCK_QUIT),

            ('/Tools/Selected Item Details...',
             '',
             self.details_dialog_cb,
             0,
             '<StockItem>',
             gtk.STOCK_ADD),

            ('/Tools/TLS Analysis...',
             '',
             self.tls_dialog_cb,
             0,
             '<StockItem>',
             gtk.STOCK_ADD),

            ('/_Help',
             None,
             None,
             0,
             '<Branch>'),
            
            ('/Help/_About',
             None,
             None, #callback
             0,
             ''))

        self.accel_group = gtk.AccelGroup()
        self.window.add_accel_group(self.accel_group)

        self.item_factory = gtk.ItemFactory(
            gtk.MenuBar, '<main>', self.accel_group)
        self.item_factory.create_items(MENU_ITEMS, self)
    
        table.attach(self.item_factory.get_widget('<main>'),
                     # X direction              Y direction
                     0, 1,                      0, 1,
                     gtk.EXPAND | gtk.FILL,     0,
                     0,                         0)

        ### /MENUBAR

        ## Create toolbar
        self.toolbar1 = gtk.Toolbar()
        table.attach(self.toolbar1,
                     # X direction              Y direction
                     0, 1,                      1, 2,
                     gtk.EXPAND | gtk.FILL,     0,
                     0,                         0)


        self.axes_tb = self.toolbar1.append_element(
            gtk.TOOLBAR_CHILD_TOGGLEBUTTON,
            None,
            "Axis",
            "Show/Hode Carteasion Axis",
            "",
            None, #widget/icon
            self.axes_tb_cb,
            None)

        self.unit_cell_tb = self.toolbar1.append_element(
            gtk.TOOLBAR_CHILD_TOGGLEBUTTON,
            None,
            "Unit Cell",
            "Show/Hode Unit Cell",
            "",
            None, #widget/icon
            self.unit_cell_tb_cb,
            None)

        
        ## Second toolbar
        self.hbox1 = gtk.HBox()
        table.attach(self.hbox1,
                     # X direction              Y direction
                     0, 1,                      2, 3,
                     gtk.EXPAND | gtk.FILL,     0,
                     0,                         0)

        self.hbox1.set_border_width(2)

        self.tb_label = gtk.Label("Current Selection:  ")
        self.hbox1.pack_start(self.tb_label, gtk.FALSE, gtk.FALSE, 0)

        self.select_label = gtk.Entry()
        self.select_label.set_editable(gtk.FALSE)
        self.hbox1.pack_start(self.select_label, gtk.TRUE, gtk.TRUE, 0)

        ## Create document
        self.struct_gui = StructureGUI(self)
        table.attach(self.struct_gui.get_widget(),
                     # X direction           Y direction
                     0, 1,                   3, 4,
                     gtk.EXPAND | gtk.FILL,  gtk.EXPAND | gtk.FILL,
                     0,                      0)

        ## Create statusbar 
        self.hbox = gtk.HBox()
        table.attach(self.hbox,
                     # X direction           Y direction
                     0, 1,                   4, 5,
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

    def new_window_cb(self, *args):
        pass

    def new_tab_cb(self, *args):
        pass
        
    def open_cb(self, *args):
        if hasattr(self, "file_selector"):
            self.file_selector.present()
            return

        self.file_selector = gtk.FileSelection("Select file to view");

        ok_button = self.file_selector.ok_button
        ok_button.connect("clicked", self.open_ok_cb, None)

        cancel_button = self.file_selector.cancel_button
        cancel_button.connect("clicked", self.open_cancel_cb, None)
        
        file_selector.present()

    def open_ok_cb(self, *args):
        path = self.file_selector.get_filename()

        self.file_selector.destroy()
        del self.file_selector
        
        self.load_file(path)

    def open_cancel_cb(self, *args):
        file_selector.destroy()

    def quit_cb(self, *args):
        self.quit_notify_cb(self.window, self)

    def error_dialog(self, text):
        print text

    def details_dialog_cb(self, *args):
        if self.selected_struct_obj == None:
            self.error_dialog("No Structure Selected.")
            return
        
        details = StructDetailsDialog(self)
        details.set_struct_obj(self.selected_struct_obj)
        details.present()

    def tls_dialog_cb(self, *args):
        struct = self.get_selected_struct()
        if struct == None:
            return
        tls = TLSDialog(context = self, struct = struct, title = str(struct))
        tls.present()

    def axes_tb_cb(self, *args):
        """Show/Hide axes toggle button.
        """
        struct = self.get_selected_struct()
        if struct == None:
            return
        if self.axes_tb.get_active() == gtk.TRUE:
            self.struct_dict[struct].show_axes(True)
        else:
            self.struct_dict[struct].show_axes(False)
        self.struct_gui.gtkglviewer.queue_draw()

    def unit_cell_tb_cb(self, *args):
        """Show/Hide unit cell.
        """
        struct = self.get_selected_struct()
        if struct == None:
            return
        if self.unit_cell_tb.get_active() == gtk.TRUE:
            self.struct_dict[struct].show_unit_cell(True)
        else:
            self.struct_dict[struct].show_unit_cell(False)
        self.struct_gui.gtkglviewer.queue_draw()

    def set_title(self, title):
        title = "mmView: " + title
        title = title[:50]
        self.window.set_title(title)

    def set_statusbar(self, text):
        self.statusbar.pop(0)
        self.statusbar.push(0, text)

    def load_file(self, path):
        """Loads the structure file specified in the path.
        """
        self.set_title(path)
        self.set_statusbar("Loading: %s" % (path))
        self.load_structure(path)
        self.set_statusbar("")
        self.progress.set_fraction(0.0)
        
        return gtk.FALSE
    
    def load_structure(self, path):
        struct = LoadStructure(
            fil              = path,
            update_cb        = self.update_cb,
            build_properties = ("sequence","bonds"))

        self.struct_list.append(struct)

        ## update the treeview/glviewer
        self.struct_gui.struct_ctrl.append_struct(struct)

        gl_struct = self.struct_gui.gtkglviewer.append_struct(struct)
        self.struct_dict[struct] = gl_struct

    def update_cb(self, percent):
        """Callback for file loading code to inform the GUI of how
        of the file has been read
        """
        self.progress.set_fraction(percent/100.0)
        while gtk.events_pending():
            gtk.main_iteration(gtk.TRUE)

    def set_selected_item(self, struct_obj):
        """Set the currently selected structure item.
        """
        if struct_obj == None:
            self.select_label.set_text("")
            self.axes_cb.set_active(gtk.FALSE)
            self.unit_cell_tb.set_active(gtk.FALSE)
            return

        self.selected_struct_obj = struct_obj
        struct = self.get_selected_struct()
        self.select_label.set_text(str(self.selected_struct_obj))
        gl_struct = self.struct_dict[struct]

        if gl_struct.gl_axes == None:
            self.axes_tb.set_active(gtk.FALSE)
        else:
            self.axes_tb.set_active(gtk.TRUE)

        if gl_struct.gl_unit_cell == None:
            self.unit_cell_tb.set_active(gtk.FALSE)
        else:
            self.unit_cell_tb.set_active(gtk.TRUE)

    def get_selected_struct(self):
        """Returns the structure object of the selected item.
        """
        if self.selected_struct_obj == None:
            self.error_dialog("No Structure Selected.")
            return None
        return self.selected_struct_obj.get_structure()


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
