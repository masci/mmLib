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

from mmLib.PDB            import PDBFile
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



def markup_vector3(vector):
    return "<small>%7.4f %7.4f %7.4f</small>" % (vector[0], vector[1], vector[2])

def markup_matrix3(tensor):
    """Uses pango markup to make the presentation of the tenosr
    look nice.
    """
    return "<small>%7.4f %7.4f %7.4f\n"\
                  "%7.4f %7.4f %7.4f\n"\
                  "%7.4f %7.4f %7.4f</small>" % (
               tensor[0,0], tensor[0,1], tensor[0,2],
               tensor[1,0], tensor[1,1], tensor[1,2],
               tensor[2,0], tensor[2,1], tensor[2,2])




class GtkGLViewer(gtk.gtkgl.DrawingArea, GLViewer):
    def __init__(self):
        gtk.gtkgl.DrawingArea.__init__(self)
        GLViewer.__init__(self)
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


    def get_gl_drawable(self):
        return gtk.gtkgl.DrawingArea.get_gl_drawable(self)

    def get_gl_context(self):
        return gtk.gtkgl.DrawingArea.get_gl_context(self)

    def destroy(self, glarea):
        return gtk.TRUE

    def button_press_event(self, glarea, event):
        self.beginx = event.x
        self.beginy = event.y

    def motion_notify_event(self, glarea, event):
        width     = glarea.allocation.width
        height    = glarea.allocation.height

        if (event.state & gtk.gdk.BUTTON1_MASK):
            self.roty = self.roty + (((event.x - self.beginx) / float(width)) * 360.0)
            self.rotx = self.rotx + (((event.y - self.beginy) / float(height)) * 360.0)

        elif (event.state & gtk.gdk.BUTTON2_MASK):
            self.zpos = self.zpos + (((event.y - self.beginy) / float(height)) * 50.0)
            self.zpos = max(self.zpos, -450.0)

        elif (event.state & gtk.gdk.BUTTON3_MASK):
            self.ypos = self.ypos - (((event.y - self.beginy) / float(height)) * 50.0)
            self.xpos = self.xpos + (((event.x - self.beginx) / float(width)) * 50.0)

        self.beginx = event.x
        self.beginy = event.y

        self.queue_draw()

    def configure_event(self, glarea, event):
        x, y, width, height = glarea.get_allocation()
        self.gl_resize(width, height)
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
        self.gl_init()
        return gtk.TRUE

    def expose_event(self, glarea, event):
        self.gl_render()
        return gtk.TRUE

    def add_struct(self, struct):
        """Adds a structure to this viewer, and returns the GLStructure
        object so it can be manipulated.
        """
        gl_struct = GLStructure(struct=struct)
        self.add_draw_list(gl_struct)
        self.queue_draw()
        return gl_struct

    def remove_struct(self, struct):
        """Removes structure from the viewer.
        """
        pass



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



class GLPropertyEditor(gtk.Frame):
    """Gtk Widget which generates a customized editing widget for a
    GLObject widget supporting GLProperties.
    """
    def __init__(self, gl_object):
        gtk.Frame.__init__(self)
        self.connect("destroy", self.destroy)

        self.gl_object = gl_object
        self.gl_object.properties.add_update_callback(self.properties_update_cb)
        
        self.prop_widget_dict = {}
        num_props = len(self.gl_object.properties.prop_list)

        ## table
        table = gtk.Table(2, num_props, gtk.FALSE)
        self.add(table)
        table.set_border_width(5)
        table.set_row_spacings(5)
        table.set_col_spacings(10)
        table_row = 0

        ## boolean types first since the toggle widgets don't look good mixed
        ## with the entry widgets
        for prop in self.gl_object.properties.prop_list:
            if prop.get("hidden", False)==True:
                continue

            ## only handling boolean right now
            if prop["type"]!="boolean":
                continue

            edit_widget = self.new_property_edit_widget(prop)
            self.prop_widget_dict[prop["name"]] = edit_widget 

            align = gtk.Alignment(0.0, 0.5, 0.0, 0.0)
            align.add(edit_widget)

            table.attach(align, 0, 2, table_row, table_row+1, gtk.EXPAND|gtk.FILL, 0, 0, 0)
            table_row += 1

        ## create and layout the editing widgets
        for prop in self.gl_object.properties.prop_list:
            if prop.get("hidden", False)==True:
                continue
            
            ## boolean types were already handled
            if prop["type"]=="boolean":
                continue

            label_widget = self.new_property_label_widget(prop)
            table.attach(label_widget, 0, 1, table_row,table_row+1, gtk.EXPAND|gtk.FILL, 0, 0, 0)

            edit_widget = self.new_property_edit_widget(prop)
            self.prop_widget_dict[prop["name"]] = edit_widget 

            table.attach(edit_widget, 1, 2, table_row, table_row+1, 0, 0, 0, 0)
            table_row += 1

        self.properties_update_cb()

    def destroy(self, widget):
        """Called after this widget is destroyed.  Make sure to remove the callback we installed
        to monitor property changes.
        """
        self.gl_object.properties.remove_update_callback(self.properties_update_cb)

    def new_property_label_widget(self, prop):
        """Returns the label widget for property editing.
        """
        label = gtk.Label(prop.get("desc", prop["name"]))
        label.set_alignment(0, 1)
        return label

    def new_property_edit_widget(self, prop):
        """Returns the editing widget for a property.
        """
        if prop["type"]=="boolean":
            widget = gtk.CheckButton(prop.get("desc", prop["name"]))
        elif prop["type"]=="integer":
            widget = gtk.Entry()
        elif prop["type"]=="float":
            widget = gtk.Entry()
        elif prop["type"]=="array(3)":
            widget = gtk.Label()
            widget.set_markup(markup_vector3(self.gl_object.properties[prop["name"]]))
        elif prop["type"]=="array(3,3)":
            widget = gtk.Label()
            widget.set_markup(markup_matrix3(self.gl_object.properties[prop["name"]]))
        else:
            text = str(self.gl_object.properties[prop["name"]])
            widget = gtk.Label(text)

        return widget

    def properties_update_cb(self, updates={}, actions=[]):
        """Read the property values and update the widgets to display the values.
        """
        for prop in self.gl_object.properties.prop_list:
            name   = prop["name"]

            try:
                widget = self.prop_widget_dict[name]
            except KeyError:
                continue

            if prop["type"]=="boolean":
                if self.gl_object.properties[name]==True:
                    widget.set_active(gtk.TRUE)
                else:
                    widget.set_active(gtk.FALSE)

            elif prop["type"]=="float":
                text = str(self.gl_object.properties[name])
                widget.set_text(text)

    def update(self):
        """Read values from widgets and apply them to the gl_object
        properties.
        """
        update_dict = {}
        
        for prop in self.gl_object.properties.prop_list:
            name   = prop["name"]

            try:
                widget = self.prop_widget_dict[name]
            except KeyError:
                continue

            if prop["type"]=="boolean":
                if widget.get_active()==gtk.TRUE:
                    update_dict[name] = True
                else:
                    update_dict[name] = False

            elif prop["type"]=="integer":
                try:
                    value = int(widget.get_text())
                except ValueError:
                    pass
                else:
                    update_dict[name] = value
                    
            elif prop["type"]=="float":
                try:
                    value = float(widget.get_text())
                except ValueError:
                    pass
                else:
                    update_dict[name] = value
                    
        self.gl_object.properties.update(**update_dict)


class GLPropertyEditDialog(gtk.Dialog):
    def __init__(self, parent_window, gl_object):
        gtk.Dialog.__init__(self, "Edit GLProperties", parent_window, gtk.DIALOG_DESTROY_WITH_PARENT)
        self.set_resizable(gtk.FALSE)
        self.connect("response", self.response_cb)

        self.gl_prop_editor = GLPropertyEditor(gl_object)
        self.vbox.pack_start(self.gl_prop_editor, gtk.TRUE, gtk.TRUE, 0)

        self.add_button(gtk.STOCK_APPLY, 100)
        self.add_button(gtk.STOCK_CLOSE, gtk.RESPONSE_CLOSE)

        self.show_all()

    def response_cb(self, dialog, response_code):
        if response_code==gtk.RESPONSE_CLOSE:
            self.destroy()
        if response_code==100:
            self.gl_prop_editor.update()
    

class TLSDialog(gtk.Dialog):
    """Dialog for TLS analysis of a structure.
    """    
    def __init__(self, **args):
        self.main_window    = args["main_window"]
        self.struct_context = args["struct_context"]
        self.sel_tls_group  = None
        self.tls_group_list = []

        gtk.Dialog.__init__(
            self,
            "TLS Analysis: %s" % (str(self.struct_context.struct)),
            self.main_window.window,
            gtk.DIALOG_DESTROY_WITH_PARENT)

        self.add_button("Open TLSIN", 100)
        self.add_button("Graphics Properties", 101)
        self.add_button(gtk.STOCK_CLOSE, gtk.RESPONSE_CLOSE)

        self.set_default_size(400, 400)
        
        self.connect("destroy", self.destroy_cb)
        self.connect("response", self.response_cb)

        ## make the print box
        sw = gtk.ScrolledWindow()
        self.vbox.add(sw)
        sw.set_border_width(5)
        sw.set_shadow_type(gtk.SHADOW_ETCHED_IN)
        sw.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)

        self.model = gtk.TreeStore(
            gobject.TYPE_BOOLEAN,   #0
            gobject.TYPE_STRING,    #1
            gobject.TYPE_STRING,    #2
            gobject.TYPE_STRING,    #3
            gobject.TYPE_STRING)    #4
     
        treeview = gtk.TreeView(self.model)
        sw.add(treeview)
        treeview.connect("row-activated", self.row_activated_cb)
        treeview.set_rules_hint(gtk.TRUE)

        cell_rend = gtk.CellRendererToggle()
        column = gtk.TreeViewColumn("Show", cell_rend)
        column.add_attribute(cell_rend, "active", 0)
        treeview.append_column(column)
        cell_rend.connect("toggled", self.view_toggled)

        cell_rend = gtk.CellRendererText()
        column = gtk.TreeViewColumn("TLS Group", cell_rend)
        column.add_attribute(cell_rend, "markup", 1)
        treeview.append_column(column)

        cell_rend = gtk.CellRendererText()
        column = gtk.TreeViewColumn("T(A^2)", cell_rend)
        column.add_attribute(cell_rend, "markup", 2)
        treeview.append_column(column)

        cell_rend = gtk.CellRendererText()
        column = gtk.TreeViewColumn("L(deg^2)", cell_rend)
        column.add_attribute(cell_rend, "markup", 3)
        treeview.append_column(column)

        cell_rend = gtk.CellRendererText()
        column = gtk.TreeViewColumn("S(A*deg)", cell_rend)
        column.add_attribute(cell_rend, "markup", 4)
        treeview.append_column(column)

        self.show_all()

        self.load_PDB(self.struct_context.struct.path)

    def response_cb(self, dialog, response_code):
        if response_code==gtk.RESPONSE_CLOSE:
            self.destroy()
            
        elif response_code==100:
            file_sel = gtk.FileSelection("Select TLSOUT File")
            file_sel.hide_fileop_buttons()
            response = file_sel.run()

            if response==gtk.RESPONSE_OK:
                path = file_sel.get_filename()
                if path!=None and path!="":
                    self.load_TLSOUT(path)
            file_sel.destroy()
            
        elif response_code==101 and self.sel_tls_group!=None:
            prop_dialog = GLPropertyEditDialog(self, self.sel_tls_group.gl_tls)
            prop_dialog.present()

    def destroy_cb(self, *args):
        self.clear_tls_groups()

    def get_tls_group(self, path):
        return self.tls_group_list[int(path)]

    def row_activated_cb(self, tree_view, path, column):
        print "## row activated path = %s" % (str(path))
        row = int(path[0])
        self.sel_tls_group = self.tls_group_list[row]

    def view_toggled(self, cell, path):
        print "## view_toggled path = %s" % (str(path))

        ## why, oh why??
        path = int(path)

        # get toggled iter
        iter = self.model.get_iter((path,))

        show_vis = self.model.get_value(iter, 0)
        tls_group = self.tls_group_list[path]
    
        # do something with the value
        if show_vis==gtk.TRUE:
            tls_group.gl_tls.properties.update(visible=False)
        else:
            tls_group.gl_tls.properties.update(visible=True)

        self.main_window.struct_gui.gtkglviewer.queue_draw()

    def add_tls_group(self, tls_group):
        """Adds the TLS group and creates tls.gl_tls OpenGL
        renderer.
        """
        self.tls_group_list.append(tls_group)
        
        tls_group.gl_tls = GLTLSGroup(tls_group=tls_group)
        tls_group.gl_tls.properties.add_update_callback(self.update_cb)

        gl_viewer = self.main_window.struct_gui.gtkglviewer
        gl_viewer.add_draw_list(tls_group.gl_tls)

        self.redraw_treeview()

    def update_cb(self, updates={}, actions=[]):
        """Property change callback from the GLTLSGroups.
        """
        self.redraw_treeview()

    def clear_tls_groups(self):
        """Remove the current TLS groups, including destroying
        the tls.gl_tls OpenGL renderer
        """
        gl_viewer = self.main_window.struct_gui.gtkglviewer

        for tls_group in self.tls_group_list:
            gl_viewer.remove_draw_list(tls_group.gl_tls)
            tls_group.gl_tls.properties.remove_update_callback(self.update_cb)
            del tls_group.gl_tls
        
        self.tls_group_list = []
        
    def load_PDB(self, path):
        """Load TLS descriptions from PDB REMARK records.
        """
        self.clear_tls_groups()
        
        tls_file = TLSInfoList()
        pdb_file = PDBFile()
        pdb_file.load_file(path)
        pdb_file.record_processor(tls_file)

        for tls_info in tls_file:
            tls_group = tls_info.make_tls_group(self.struct_context.struct)
            tls_group.tls_info = tls_info
            self.add_tls_group(tls_group)
            
    def load_TLSOUT(self, path):
        """Load TLS descriptions from a REMAC/CCP4 TLSOUT file.
        """
        print "## loading = ",path
        
        self.clear_tls_groups()

        tls_file = TLSInfoList()
        tls_file.load_refmac_tlsout_file(open(path, "r"))
    
        for tls_info in tls_file:
            tls_group = tls_info.make_tls_group(self.struct_context.struct)
            tls_group.tls_info = tls_info
            self.add_tls_group(tls_group)

    def markup_tls_name(self, tls_info):
        listx = []
        for (chain_id1, frag_id1, chain_id2, frag_id2, sel) in tls_info.range_list:
            listx.append("%s%s-%s%s %s" % (chain_id1, frag_id1, chain_id2, frag_id2, sel))
        return "<small>"+string.join(listx, "\n")+"</small>"

    def markup_tensor(self, tensor):
        """Uses pango markup to make the presentation of the tenosr
        look nice.
        """
        return "<small>%7.4f %7.4f %7.4f\n"\
                      "%7.4f %7.4f %7.4f\n"\
                      "%7.4f %7.4f %7.4f</small>" % (
                   tensor[0,0], tensor[0,1], tensor[0,2],
                   tensor[1,0], tensor[1,1], tensor[1,2],
                   tensor[2,0], tensor[2,1], tensor[2,2])
           
    def redraw_treeview(self):
        """Clear and redisplay the TLS group treeview list.
        """
        self.model.clear()
        self.model.clear()

        for tls in self.tls_group_list:
            iter = self.model.append(None)

            if tls.gl_tls.properties["visible"]==True:
                self.model.set(iter, 0, gtk.TRUE)
            else:
                self.model.set(iter, 0, gtk.FALSE)

            self.model.set(iter, 1, self.markup_tls_name(tls.tls_info))
            self.model.set(iter, 2, self.markup_tensor(tls.T))
            self.model.set(iter, 3, self.markup_tensor(tls.L*rad2deg2))
            self.model.set(iter, 4, self.markup_tensor(tls.S*rad2deg))
        
    def timeout_cb(self, tls):
        tls.gl_tls.inc_time()
        self.context.struct_gui.gtkglviewer.queue_draw()
        return gtk.TRUE


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
        self.context.set_selected_struct_obj(self.resolv_path(path))
    
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


class StructureGUI(object):
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


class ViewCommands(list):
    def __init__(self):
        view_cmds = [
            { "menu path":       "/View/Carteasion Axes",
              "glstruct property": "axes_visible",
              "action":          100 },

            { "menu path":         "/View/Unit Cell",
              "glstruct property": "unit_cell_visible",
              "action":            101 },

            { "menu path":       "/View/Protein Main Chain",
              "glstruct property": "aa_main_chain_visible",
              "action":          102 },

            { "menu path":       "/View/Protein Side Chain",
              "glstruct property": "aa_side_chain_visible",
              "action":          103 },

            { "menu path":       "/View/DNA Main Chain",
              "glstruct property": "dna_main_chain_visible",
              "action":          104 },

            { "menu path":       "/View/DNA Side Chain",
              "glstruct property": "dna_side_chain_visible",
              "action":          105 },

            { "menu path":       "/View/HET Group",
              "glstruct property": "hetatm_visible",
              "action":          106 },

            { "menu path":       "/View/Water",
              "glstruct property": "water_visible",
              "action":          107},
            
            { "menu path":       "/View/Thermal Axes",
              "glstruct property": "U",
              "action":          108},
            ]

        list.__init__(self, view_cmds)

    def get_action(self, action):
        """Return the cmd dictionary by action number.
        """
        for cmd in self:
            if cmd["action"]==action:
                return cmd
        return None


class ColorCommands(list):
    def __init__(self):
        list.__init__(self, [
            { "menu path":  "/Colors/Color by Atom Type",
              "action":     100 },
            { "menu path":  "/Colors/Color by Chain",
              "action":     101 } ])

    def get_action(self, action):
        """Return the cmd dictionary by action number.
        """
        for cmd in self:
            if cmd["action"]==action:
                return cmd
        return None


class StructureContext(object):
    def __init__(self, struct, gl_struct):
        self.struct    = struct
        self.gl_struct = gl_struct


class MainWindow(object):
    def __init__(self, quit_notify_cb):
        self.struct_context_list = []
        self.sel_struct_context  = None
        self.sel_struct_obj      = None
        self.view_cmds           = ViewCommands()
        self.color_cmds          = ColorCommands()
        self.quit_notify_cb      = quit_notify_cb

        ## dialog lists
        self.details_dialog_list = []
        self.tls_dialog_list     = []

        ## Create the toplevel window
        self.window = gtk.Window()
        self.set_title("")
        self.window.set_default_size(500, 400)
        self.window.connect('destroy', self.quit_cb, self)

        table = gtk.Table(1, 4, gtk.FALSE)
        self.window.add(table)

        ## file menu bar
        file_menu_items = [
            ('/_File',            None,          None,               0,'<Branch>'),
            ('/File/_New Window', None,          self.new_window_cb, 0,'<StockItem>',gtk.STOCK_NEW),
            ('/File/_Open',       None,          self.open_cb,       0,'<StockItem>',gtk.STOCK_OPEN),
            ('/File/sep1',        None,          None,               0,'<Separator>'),
            ('/File/_Quit',      '<control>Q',   self.quit_cb,       0,'<StockItem>',gtk.STOCK_QUIT) ]

        view_menu_items = [
            ('/_View', None, None, 0, '<Branch>') ]
        
        for view_cmd in self.view_cmds:
            view_menu_items.append(
                (view_cmd["menu path"], None, self.view_menu_cb, view_cmd["action"], '<CheckItem>') )

        color_menu_items = [
            ('/_Colors', None, None, 0, '<Branch>') ]
        
        for color_cmd in self.color_cmds:
            color_menu_items.append(
                (color_cmd["menu path"], None, self.color_menu_cb, color_cmd["action"]) )
            
        tools_menu_items = [
            ('/_Tools', None, None, 0, '<Branch>'),
            ('/Tools/Selected Item Details...', None, self.details_dialog_cb, 0, None),
            ('/Tools/TLS Analysis...',          None, self.tls_dialog_cb,     0, None) ]

        help_menu_items = [
            ('/_Help',       None, None, 0, '<Branch>'),
            ('/Help/_About', None, None, 0, None) ]

        menu_items = file_menu_items +\
                     view_menu_items +\
                     color_menu_items +\
                     tools_menu_items +\
                     help_menu_items

        self.accel_group = gtk.AccelGroup()
        self.window.add_accel_group(self.accel_group)
        self.item_factory = gtk.ItemFactory(gtk.MenuBar, '<main>', self.accel_group)
        self.item_factory.create_items(menu_items)
    
        table.attach(self.item_factory.get_widget('<main>'),
                     # X direction              Y direction
                     0, 1,                      0, 1,
                     gtk.EXPAND | gtk.FILL,     0,
                     0,                         0)

        ## Toolbar
        self.hbox1 = gtk.HBox()
        table.attach(self.hbox1,
                     # X direction              Y direction
                     0, 1,                      1, 2,
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
                     0, 1,                   2, 3,
                     gtk.EXPAND | gtk.FILL,  gtk.EXPAND | gtk.FILL,
                     0,                      0)

        ## Create statusbar 
        self.hbox = gtk.HBox()
        table.attach(self.hbox,
                     # X direction           Y direction
                     0, 1,                   3, 4,
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
        self.file_selector.present()

    def destroy_file_selector(self):
        self.file_selector.destroy()
        del self.file_selector

    def open_ok_cb(self, *args):
        path = self.file_selector.get_filename()
        self.destroy_file_selector()
        self.load_file(path)

    def open_cancel_cb(self, *args):
        self.destroy_file_selector()

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
        if self.sel_struct_context==None:
            return
        tls = TLSDialog(
            main_window    = self,
            struct_context = self.sel_struct_context)
        tls.present()

    def view_menu_cb(self, callback_action, widget):
        """Callback for the View menu.
        """
        if self.sel_struct_context==None:
            return
        view_cmd  = self.view_cmds.get_action(callback_action)
        property  = view_cmd["glstruct property"]
        menu_item = self.item_factory.get_item(view_cmd["menu path"])
        if menu_item.get_active() == gtk.TRUE:
            visible = True
        else:
            visible = False
        self.sel_struct_context.gl_struct.properties.update(**{property: visible})            
        self.struct_gui.gtkglviewer.queue_draw()

    def color_menu_cb(self, callback_action, widget):
        """Callback for the Color menu.
        """
        if self.sel_struct_context==None:
            return
        color_cmd = self.color_cmds.get_action(callback_action)

        if callback_action==100:
            self.sel_struct_context.gl_struct.properties.update(color=None)

        elif callback_action==101: 
            gl_struct = self.sel_struct_context.gl_struct

            color_list = [(1.,0.,0.),(0.,1.,0.),(0.,0.,1.),(1.,1.,0.),(0.,1.,1.),(1.,0.,1.),(1.,1.,1.)]
            colori = 0

            chain_ids = gl_struct.gl_chain_dict.keys()
            chain_ids.sort()
            
            for chain_id in chain_ids:
                gl_chain = gl_struct.gl_chain_dict[chain_id]
                
                try:
                    color = color_list[colori]
                except IndexError:
                    color = color_list[-1]
                else:
                    colori += 1
                
                gl_chain.properties.update(color=color)

        self.struct_gui.gtkglviewer.queue_draw()

    def set_title(self, title):
        self.window.set_title("Viewer: %s" % (title[:50]))

    def set_statusbar(self, text):
        self.statusbar.pop(0)
        self.statusbar.push(0, text)

    def load_file(self, path):
        """Loads the structure file specified in the path.
        """
        self.set_title(path)
        self.set_statusbar("Loading: %s" % (path))

        struct = LoadStructure(
            fil              = path,
            update_cb        = self.update_cb,
            build_properties = ("sequence","bonds"))

        struct.path = path
        
        gl_struct = self.struct_gui.gtkglviewer.add_struct(struct)

        struct_context = StructureContext(struct, gl_struct)
        self.struct_context_list.append(struct_context)
        self.struct_gui.struct_ctrl.append_struct(struct_context.struct)
        if self.sel_struct_context==None:
            self.set_selected_struct_obj(struct)

        self.set_statusbar("")
        self.progress.set_fraction(0.0)
        return gtk.FALSE
    
    def update_cb(self, percent):
        """Callback for file loading code to inform the GUI of how
        of the file has been read
        """
        self.progress.set_fraction(percent/100.0)
        while gtk.events_pending():
            gtk.main_iteration(gtk.TRUE)

    def set_selected_struct_obj(self, struct_obj):
        """Set the currently selected structure item.
        """
        self.sel_struct_obj = struct_obj
        
        ## nothing selected
        if struct_obj==None:
            self.sel_struct_context = None
            self.select_label.set_text("")

            for view_cmd in self.view_cmds:
                menu_item = self.item_factory.get_item(view_cmd["menu path"])
                menu_item.set_active(view_cmd["checked"])

        else:
            self.select_label.set_text(str(self.sel_struct_obj))

            struct = struct_obj.get_structure()
            self.sel_struct_context = self.get_struct_context(struct)
            self.set_selected_struct_context(self.sel_struct_context)

    def set_selected_struct_context(self, struct_context):
        """Sets the current struct_context and updates the context
        sensitive GUI widgets.
        """
        for view_cmd in self.view_cmds:
            property  = view_cmd["glstruct property"]
            menu_item = self.item_factory.get_item(view_cmd["menu path"])
            if struct_context.gl_struct.properties[property]==False:
                menu_item.set_active(gtk.FALSE)
            else:
                menu_item.set_active(gtk.TRUE)

    def get_struct_context(self, struct):
        """Returns the struct_context containing struct.
        """
        for struct_context in self.struct_context_list:
            if struct_context.struct==struct:
                return struct_context
        return None


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
