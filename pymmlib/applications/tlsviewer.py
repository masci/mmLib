#!/usr/bin/env python
## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import math
import random

import pygtk
pygtk.require("2.0")

import gobject
import gtk
import gtk.gtkgl

from OpenGL.GL            import *
from OpenGL.GLU           import *
from OpenGL.GLUT          import *

from mmLib.PDB            import PDBFile
from mmLib.Structure      import *
from mmLib.FileLoader     import *
from mmLib.Viewer         import *
from mmLib.R3DDriver      import Raster3DDriver
from mmLib.Extensions.TLS import *


###############################################################################
### Debuggin
###

def print_U(tensor):
    print "%7.4f %7.4f %7.4f\n"\
          "%7.4f %7.4f %7.4f\n"\
          "%7.4f %7.4f %7.4f" % (
        tensor[0,0], tensor[0,1], tensor[0,2],
        tensor[1,0], tensor[1,1], tensor[1,2],
        tensor[2,0], tensor[2,1], tensor[2,2])

    
###############################################################################
### Utility Widgets
###

class DictListTreeView(gtk.TreeView):
    """If you ever thought a multi-column listbox should be as easy as
    creating a list of Python dictionaries, this class is for you!
    """
    def __init__(self,
                 column_list    = ["Empty"],
                 dict_list      = [],
                 selection_mode = gtk.SELECTION_BROWSE):

        ## DistList state
        self.column_list = column_list
        self.dict_list   = []

        ## gtk.TreeView Init
        gtk.TreeView.__init__(self)
        
        self.get_selection().set_mode(selection_mode)

        self.connect("row-activated", self.row_activated_cb)
        self.connect("button-release-event", self.button_release_event_cb)

        model_data_types = []
        for i in range(len(self.column_list)):
            model_data_types.append(gobject.TYPE_STRING)

        self.model = gtk.TreeStore(*model_data_types)
        self.set_model(self.model)

        ## add cell renderers
        for i in range(len(self.column_list)):
            column_desc = self.column_list[i]
            cell        = gtk.CellRendererText()
            column      = gtk.TreeViewColumn(column_desc, cell)

            column.add_attribute(cell, "text", i)
            self.append_column(column)

        ## populate list if needed
        if len(dict_list)>0:
            self.set_dict_list(dict_list)

    def row_activated_cb(self, tree_view, path, column):
        """Retrieve selected node, then call the correct set method for the
        type.
        """
        pass
    
    def button_release_event_cb(self, tree_view, bevent):
        """
        """
        x = int(bevent.x)
        y = int(bevent.y)

        try:
            (path, col, x, y) = self.get_path_at_pos(x, y)
        except TypeError:
            return gtk.FALSE

        self.row_activated_cb(tree_view, path, col)
        return gtk.FALSE

    def set_dict_list(self, dict_list):
        """
        """
        self.dict_list = dict_list
        self.model.clear()

        for dictX in self.dict_list:
            miter = self.model.append(None)

            for i in range(len(self.column_list)):
                column_desc = self.column_list[i]
                self.model.set(miter, i, dictX.get(column_desc, ""))


###############################################################################
### GTK OpenGL Viewer Widget Using GtkGlExt/PyGtkGLExt
###
### see http://gtkglext.sourceforge.net/ for information on GtkGLExt and
### OpenGL binding in GTK 2.0
###

class GtkGLViewer(gtk.gtkgl.DrawingArea, GLViewer):
    def __init__(self):
        self.in_drag = False
        self.beginx  = 0
        self.beginy  = 0

        gtk.gtkgl.DrawingArea.__init__(self)
        GLViewer.__init__(self)

        self.r3d_driver = Raster3DDriver()

        try:
            glconfig = gtk.gdkgl.Config(
                mode = gtk.gdkgl.MODE_RGB|
                gtk.gdkgl.MODE_DOUBLE|
                gtk.gdkgl.MODE_DEPTH)
            
        except gtk.gdkgl.NoMatches:
            glconfig = gtk.gdkgl.Config(
                mode = gtk.gdkgl.MODE_RGB|
                gtk.gdkgl.MODE_DEPTH)

        gtk.gtkgl.widget_set_gl_capability(self, glconfig)
        
        self.set_size_request(300, 300)

        self.connect('button_press_event',   self.gtk_glv_button_press_event)
        self.connect('button_release_event', self.gtk_glv_button_release_event)
        self.connect('motion_notify_event',  self.gtk_glv_motion_notify_event)
        self.connect('realize',              self.gtk_glv_realize)
        self.connect('map',                  self.gtk_glv_map)
        self.connect('configure_event',      self.gtk_glv_configure_event)
        self.connect('expose_event',         self.gtk_glv_expose_event)
        self.connect('destroy',              self.gtk_glv_destroy)

    def gl_begin(self):
        """Sets up the OpenGL drawing context for drawing.  If
        no context is availible (hiddent window, not created yet),
        this method returns False, otherwise it returns True.
        """
        gl_drawable = gtk.gtkgl.DrawingArea.get_gl_drawable(self)
        gl_context  = gtk.gtkgl.DrawingArea.get_gl_context(self)
        if gl_drawable==None or gl_context==None:
            return False
        
	if not gl_drawable.gl_begin(gl_context):
            return False

        return True

    def gl_end(self):
        """Ends a sequence of OpenGL drawing calls, then swaps drawing
        buffers.
        """
        gl_drawable = gtk.gtkgl.DrawingArea.get_gl_drawable(self)
	if gl_drawable.is_double_buffered():
            gl_drawable.swap_buffers()
	else:
            glFlush()
        
        gl_drawable.gl_end()

    def gtk_glv_map(self, glarea):
        """
        """
        self.add_events(gtk.gdk.BUTTON_PRESS_MASK   |
                        gtk.gdk.BUTTON_RELEASE_MASK |
                        gtk.gdk.BUTTON_MOTION_MASK  |
                        gtk.gdk.POINTER_MOTION_MASK)
        return gtk.TRUE

    def gtk_glv_realize(self, glarea):
        """
        """
        if self.gl_begin()==True:
            self.glv_init()
            self.gl_end()
            return gtk.TRUE
        else:
            return gtk.FALSE

    def gtk_glv_configure_event(self, glarea, event):
        """
        """
        x, y, width, height = self.get_allocation()

        if self.gl_begin()==True:
            self.glv_resize(width, height)
            self.gl_end()
            self.queue_draw()
            return gtk.TRUE
        else:
            return gtk.FALSE

    def gtk_glv_expose_event(self, glarea, event):
        """
        """
        if self.gl_begin()==True:
            self.glv_render()
            self.gl_end()
            return gtk.TRUE
        else:
            return gtk.FALSE

    def gtk_glv_button_press_event(self, glarea, event):
        if event.type==gtk.gdk._2BUTTON_PRESS:
            self.glv_render_one(self.r3d_driver)
            return

        self.in_drag = True
        self.beginx  = event.x
        self.beginy  = event.y

    def gtk_glv_button_release_event(self, glarea, event):
        self.in_drag = False
        self.beginx  = 0
        self.beginy  = 0

    def gtk_glv_motion_notify_event(self, glarea, event):
        """
        """
        if not self.in_drag:
            return
        
        if (event.state & gtk.gdk.BUTTON1_MASK):
            self.glv_trackball(self.beginx, self.beginy, event.x, event.y)

        elif (event.state & gtk.gdk.BUTTON2_MASK):
            dx = event.x - self.beginx
            dy = self.beginy - event.y
            self.glv_straif(dx, dy)

        elif (event.state & gtk.gdk.BUTTON3_MASK):
            if (event.state & gtk.gdk.SHIFT_MASK):
                dx = event.x - self.beginx
                dy = self.beginy - event.y
                self.glv_clip(dy, dx)
            else:
                dy = event.y - self.beginy
                self.glv_zoom(dy)

        self.beginx = event.x
        self.beginy = event.y

        self.queue_draw()

    def gtk_glv_destroy(self, glarea):
        ## XXX: delete all opengl draw lists
        return gtk.TRUE

    def glv_redraw(self):
        #self.snapshot()
        self.queue_draw()

    def snapshot(self):
        print [x for x in dir(self) if x.count("pixbuf")]
        
        


###############################################################################
### GLProperty Browser Components for browsing and changing GLViewer.GLObject
### properties
###


def markup_vector3(vector):
    return "<small>%7.4f %7.4f %7.4f</small>" % (
        vector[0], vector[1], vector[2])

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


class EnumeratedStringEntry(gtk.Combo):
    """Allows the entry of a custom string or a choice of a list of
    enumerations.
    """
    def __init__(self, enum_list):
        gtk.Combo.__init__(self)
        self.entry.connect("activate", self.activate, None)

        self.enum_list = enum_list[:]
        self.string    = None

    def activate(self, entry, junk):
        txt = self.entry.get_text()
        self.set_string(txt)

    def set_string(self, string):
        self.string = string
        self.set_popdown_strings(self.enum_list)
        self.entry.set_text(string)

    def get_string(self):
        self.activate(self.entry, None)
        return self.string


class GLPropertyEditor(gtk.Notebook):
    """Gtk Widget which generates a customized editing widget for a
    GLObject widget supporting GLProperties.
    """
    def __init__(self, gl_object):
        gtk.Notebook.__init__(self)
        self.set_scrollable(gtk.TRUE)
        self.connect("destroy", self.destroy)

        self.gl_object = gl_object
        self.gl_object.glo_add_update_callback(self.properties_update_cb)

        ## property name -> widget dictionary
        self.prop_widget_dict = {}

        ## count the number of properties/pages to be displayed
        catagory_dict = self.catagory_dict_sort()
        
        ## sort pages: make sure Show/Hide is the first page
        catagories = catagory_dict.keys()
        if "Show/Hide" in catagories:
            catagories.remove("Show/Hide")
            catagories.insert(0, "Show/Hide")
            
        ## add Notebook pages and tables
        for catagory in catagories:
            table = self.build_catagory_widgets(catagory, catagory_dict)
            self.append_page(table, gtk.Label(catagory))

        self.properties_update_cb(self.gl_object.properties)

    def catagory_dict_sort(self):
        """Returns a tree of dictionaries/lists  catagory_dict[]->[prop list]
        """
        
        ## count the number of properties/pages to be displayed
        catagory_dict = {}
        
        for prop_desc in self.gl_object.glo_iter_property_desc():

            ## skip hidden properties
            if prop_desc.get("hidden", False)==True:
                continue

            ## get the catagory name, default to "Misc"
            catagory = prop_desc.get("catagory", "Misc")
            try:
                catagory_dict[catagory].append(prop_desc)
            except KeyError:
                catagory_dict[catagory] = [prop_desc]

        return catagory_dict

    def build_catagory_widgets(self, catagory, catagory_dict):
        """Returns the GtkTable widget will all the control widgets
        """
        prop_desc_list = catagory_dict[catagory]
        num_properties = len(prop_desc_list)

        ## create table widget
        table = gtk.Table(2, num_properties, gtk.FALSE)
        table.set_border_width(5)
        table.set_row_spacings(5)
        table.set_col_spacings(10)

        ## size group widget to make the edit/display widgets
        ## in the right column all the same size
        size_group = gtk.SizeGroup(gtk.SIZE_GROUP_HORIZONTAL)

        table_row  = 0

        ## boolean types first since the toggle widgets don't look good mixed
        ## with the entry widgets
        for prop_desc in prop_desc_list:
            name = prop_desc["name"]
            
            ## only handling boolean right now
            if prop_desc["type"]!="boolean":
                continue

            ## create the widget for this property
            edit_widget = self.new_property_edit_widget(prop_desc)

            ## add the widget to a class dict so property name->widget is
            ## easy to look up
            self.prop_widget_dict[name] = edit_widget 

            ## use a alignment widget in the table
            align = gtk.Alignment(0.0, 0.5, 1.0, 0.0)
            align.add(edit_widget)

            ## attach to table
            table.attach(align,
                         0, 2, table_row, table_row+1,
                         gtk.FILL|gtk.EXPAND, 0, 0, 0)

            table_row += 1


        ## now create all the widgets for non-boolean types
        for prop_desc in prop_desc_list:
            name = prop_desc["name"]

            ## boolean types were already handled
            if prop_desc["type"]=="boolean":
                continue

            ## create the labell widget for the property and attach
            ## it to the left side of the table
            label_widget = self.new_property_label_widget(prop_desc)
            align = gtk.Alignment(0.0, 0.5, 0.0, 0.0)
            align.add(label_widget)
            table.attach(align,
                         0, 1, table_row, table_row+1,
                         gtk.FILL, 0, 0, 0)

            ## create the edit widget
            edit_widget = self.new_property_edit_widget(prop_desc)
            self.prop_widget_dict[name] = edit_widget 

            ## use alignment widget in table
            align = gtk.Alignment(0.0, 0.5, 1.0, 0.0)
            align.add(edit_widget)

            ## add to size group and attach to table
            size_group.add_widget(edit_widget)
            
            table.attach(align,
                         1, 2, table_row, table_row+1,
                         gtk.EXPAND|gtk.FILL, 0, 0, 0)

            table_row += 1

        table.show_all()
        return table

    def destroy(self, widget):
        """Called after this widget is destroyed.  Make sure to remove the
        callback we installed to monitor property changes.
        """
        self.gl_object.glo_remove_update_callback(self.properties_update_cb)

    def markup_property_label(self, prop_desc, max_len=40):
        """
        """
        listx    = prop_desc.get("desc", prop_desc["name"]).split()
        strx     = ""
        line_len = 0
        
        for word in listx:
            new_line_len = line_len + len(word)
            if new_line_len>max_len:
                strx     += "\n" + word
                line_len = len(word)
            elif line_len==0:
                strx     += word
                line_len += len(word)
            else:
                strx     += " " + word
                line_len += len(word) + 1

        return "<small>%s</small>" % (strx)

    def new_property_label_widget(self, prop):
        """Returns the label widget for property editing.
        """
        label = gtk.Label()
        label.set_markup(self.markup_property_label(prop))
        label.set_alignment(0.0, 0.5)
        return label

    def new_property_edit_widget(self, prop):
        """Returns the editing widget for a property.
        """

        ## BOOLEAN
        if prop["type"]=="boolean":
            widget = gtk.CheckButton(prop.get("desc", prop["name"]))
            widget.get_child().set_markup(
                self.markup_property_label(prop, max_len=50))

        ## INTEGER
        elif prop["type"]=="integer":
            if prop.get("read_only", False)==True:
                widget = gtk.Label()
            else:
                if prop.get("range")!=None:
                    ## range format: min-max,step
                    min_max, step = prop["range"].split(",")
                    min, max      = min_max.split("-")

                    min  = int(min)
                    max  = int(max)
                    step = int(step)

                    widget = gtk.HScale()
                    widget.set_digits(0)
                    widget.set_range(float(min), float(max))
                    widget.set_increments(float(step), float(step))

                elif prop.get("spin")!=None:
                    ## range format: min-max,step
                    min_max, step = prop["spin"].split(",")
                    min, max      = min_max.split("-")
                    
                    min  = int(min)
                    max  = int(max)
                    step = int(step)

                    widget = gtk.SpinButton(climb_rate=float(step), digits=0)
                    widget.set_range(float(min), float(max))
                    widget.set_increments(float(step), float(step*10))
                else:
                    widget = gtk.Entry()

        ## FLOAT
        elif prop["type"]=="float":
            if prop.get("read_only", False)==True:
                widget = gtk.Label()
            else:
                if prop.get("range")!=None:
                    ## range format: min-max,step
                    min_max, step = prop["range"].split(",")
                    min, max      = min_max.split("-")

                    min  = float(min)
                    max  = float(max)
                    step = float(step)

                    widget = gtk.HScale()
                    widget.set_digits(2)
                    widget.set_range(min, max)
                    widget.set_increments(step, step)

                elif prop.get("spin")!=None:
                    ## range format: min-max,step
                    min_max, step = prop["spin"].split(",")
                    min, max      = min_max.split("-")
                    
                    min  = float(min)
                    max  = float(max)
                    step = float(step)

                    widget = gtk.SpinButton(climb_rate=step, digits=2)
                    widget.set_range(min, max)
                    widget.set_increments(step, step*10.0)

                else:
                    widget = gtk.Entry()

        ## ARRAY(3)
        elif prop["type"]=="array(3)":
            widget = gtk.Label()

        ## ARRAY(3,3)
        elif prop["type"]=="array(3,3)":
            widget = gtk.Label()
            
        ## ENUMERATED STRING
        elif prop["type"]=="enum_string":
            widget = EnumeratedStringEntry(prop["enum_list"])
        
        ## WTF?
        else:
            text = str(self.gl_object.properties[prop["name"]])
            widget = gtk.Label(text)

        return widget

    def properties_update_cb(self, updates={}, actions=[]):
        """Read the property values and update the widgets to display the
        values.
        """
        for name in updates:
            try:
                widget = self.prop_widget_dict[name]
            except KeyError:
                continue

            prop_desc = self.gl_object.glo_get_property_desc(name)

            if prop_desc["type"]=="boolean":
                if self.gl_object.properties[name]==True:
                    widget.set_active(gtk.TRUE)
                else:
                    widget.set_active(gtk.FALSE)

            elif prop_desc["type"]=="integer":
                if prop_desc.get("read_only", False)==True:
                    text = "<small>%d</small>" % (
                        self.gl_object.properties[name])
                    widget.set_markup(text)
                else:
                    if prop_desc.get("range")!=None:
                        widget.set_value(
                            float(self.gl_object.properties[name]))
                    elif prop_desc.get("spin")!=None:
                        widget.set_value(
                            float(self.gl_object.properties[name]))
                    else:
                        text = str(self.gl_object.properties[name])
                        widget.set_text(text)

            elif prop_desc["type"]=="float":
                if prop_desc.get("read_only", False)==True:
                    widget.set_markup(
                        "<small>%12.6f</small>" % (
                        self.gl_object.properties[name]))
                else:
                    if prop_desc.get("range")!=None:
                        widget.set_value(self.gl_object.properties[name])
                    elif prop_desc.get("spin")!=None:
                        widget.set_value(self.gl_object.properties[name])
                    else:
                        text = str(self.gl_object.properties[name])
                        widget.set_text(text)

            elif prop_desc["type"]=="array(3)":
                widget.set_markup(
                    markup_vector3(self.gl_object.properties[name]))

            elif prop_desc["type"]=="array(3,3)":
                widget.set_markup(
                    markup_matrix3(self.gl_object.properties[name]))

            elif prop_desc["type"]=="enum_string":
                widget.set_string(self.gl_object.properties[name])

            else:
                widget.set_text(str(self.gl_object.properties[name]))

    def update(self):
        """Read values from widgets and apply them to the gl_object
        properties.
        """
        update_dict = {}
        
        for prop in self.gl_object.glo_iter_property_desc():

            ## skip read_only properties
            if prop.get("read_only", False)==True:
                continue

            ## property name
            name = prop["name"]

            try:
                widget = self.prop_widget_dict[name]
            except KeyError:
                continue

            ## retrieve data based on widget type
            if prop["type"]=="boolean":
                if widget.get_active()==gtk.TRUE:
                    update_dict[name] = True
                else:
                    update_dict[name] = False

            elif prop["type"]=="integer":
                if prop.get("range")!=None:
                    update_dict[name] = int(widget.get_value())
                elif prop.get("spin")!=None:
                    update_dict[name] = int(widget.get_value())
                else:
                    try:
                        update_dict[name] = int(widget.get_text())
                    except ValueError:
                        pass
                    
            elif prop["type"]=="float":
                if prop.get("range")!=None:
                    update_dict[name] = float(widget.get_value())
                elif prop.get("spin")!=None:
                    update_dict[name] = float(widget.get_value())
                else:
                    try:
                        update_dict[name] = float(widget.get_text())
                    except ValueError:
                        pass

            elif prop["type"]=="enum_string":
                update_dict[name] = widget.get_string()
                
        self.gl_object.glo_update_properties(**update_dict)


class GLPropertyTreeControl(gtk.TreeView):
    """Hierarchical tree view of the GLObjects which make up
    the molecular viewer.  The purpose of this widget is to alllow
    easy navigation through the hierarchy.
    """
    
    def __init__(self, gl_object_root, glo_select_cb):
        self.gl_object_root = gl_object_root
        self.glo_select_cb  = glo_select_cb
        self.path_glo_dict  = {}
        
        gtk.TreeView.__init__(self)
        self.get_selection().set_mode(gtk.SELECTION_BROWSE)
        
        self.connect("row-activated", self.row_activated_cb)
        self.connect("button-release-event", self.button_release_event_cb)

        self.model = gtk.TreeStore(gobject.TYPE_STRING)
        self.set_model(self.model)

        cell = gtk.CellRendererText()
        self.column = gtk.TreeViewColumn("Graphics Objects", cell)
        self.column.add_attribute(cell, "text", 0)
        self.append_column(self.column)

        self.rebuild_gl_object_tree()

    def row_activated_cb(self, tree_view, path, column):
        """Retrieve selected node, then call the correct set method for the
        type.
        """
        glo = self.get_gl_object(path)
        self.glo_select_cb(glo)

    def button_release_event_cb(self, tree_view, bevent):
        """
        """
        x = int(bevent.x)
        y = int(bevent.y)

        try:
            (path, col, x, y) = self.get_path_at_pos(x, y)
        except TypeError:
            return gtk.FALSE

        self.row_activated_cb(tree_view, path, col)
        return gtk.FALSE

    def get_gl_object(self, path):
        """Return the GLObject from the treeview path.
        """
        glo = self.path_glo_dict[path]
        return glo

    def rebuild_gl_object_tree(self):
        """Clear and refresh the view of the widget according to the
        self.struct_list
        """

        ## first get a refernence to the current selection
        ## so it can be restored after the reset if it hasn't
        ## been removed
        model, miter = self.get_selection().get_selected()
        if miter!=None:
            selected_glo = self.path_glo_dict[model.get_path(miter)]
        else:
            selected_glo = None

        ## now record how the tree is expanded so the expansion
        ## can be restored after rebuilding
        expansion_dict = {}
        for path, glo in self.path_glo_dict.items():
            expansion_dict[glo] = self.row_expanded(path)

        ## rebuild the tree view using recursion
        self.path_glo_dict = {}
        self.model.clear()
        
        def redraw_recurse(glo, parent_miter):
            miter = self.model.append(parent_miter)
            path = self.model.get_path(miter)

            self.path_glo_dict[path] = glo
            
            self.model.set(miter, 0, glo.glo_name())

            for child in glo.glo_iter_children():
                redraw_recurse(child, miter)
    
        redraw_recurse(self.gl_object_root, None)

        ## restore expansions and selections
        for path, glo in self.path_glo_dict.items():
            try:
                expanded = expansion_dict[glo]
            except KeyError:
                continue
            if expanded==gtk.TRUE:
                for i in range(1, len(path)):
                    self.expand_row(path[:i], gtk.FALSE)

        if selected_glo!=None:
            self.select_gl_object(selected_glo)

    def select_gl_object(self, target_glo):
        """Selects a gl_object which must be a glo_root or a
        decendant of glo_root.
        """
        ## try to find and select target_glo
        for path, glo in self.path_glo_dict.items():
            if id(glo)==id(target_glo):

                for i in range(1, len(path)):
                    self.expand_row(path[:i], gtk.FALSE)

                self.get_selection().select_path(path)
                self.glo_select_cb(glo)
                return

        ## unable to find target_glo, select the root
        self.select_gl_object(self.gl_object_root)
        self.glo_select_cb(self.gl_object_root)


class GLPropertyEditDialog(gtk.Dialog):
    """Open a dialog for editing a single GLObject properties.
    """

    def __init__(self, parent_window, gl_object):
        title = "Edit GLProperties: %s" % (gl_object.glo_name())

        gtk.Dialog.__init__(
            self,
            title,
            parent_window,
            gtk.DIALOG_DESTROY_WITH_PARENT)

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
    

class GLPropertyBrowserDialog(gtk.Dialog):
    """Dialog for manipulating the properties of a GLViewer.
    """
    def __init__(self, **args):

        parent_window  = args["parent_window"]
        gl_object_root = args["glo_root"]
        title          = args.get("title", "")

        gl_object_root.glo_set_name(title)

        title = "Visualization Properties: %s" % (title)

        gtk.Dialog.__init__(
            self,
            title,
            parent_window,
            gtk.DIALOG_DESTROY_WITH_PARENT)

        self.connect("response", self.response_cb)
        self.set_default_size(500, 400)
        self.add_button(gtk.STOCK_APPLY, 100)
        self.add_button(gtk.STOCK_CLOSE, gtk.RESPONSE_CLOSE)

        ## widgets
        self.hpaned = gtk.HPaned()
        self.vbox.pack_start(self.hpaned, gtk.TRUE, gtk.TRUE, 0)
        self.hpaned.set_border_width(2)

        ## property tree control
        self.sw1 = gtk.ScrolledWindow()
        self.hpaned.add1(self.sw1)
        self.sw1.set_size_request(150, -1)
        
        self.gl_tree_ctrl = GLPropertyTreeControl(
            gl_object_root,
            self.gl_tree_ctrl_selected)

        self.sw1.add(self.gl_tree_ctrl)

        ## always start with the root selected
        self.gl_prop_editor = None
        self.select_gl_object(gl_object_root)

        self.show_all()

    def response_cb(self, dialog, response_code):
        if response_code==gtk.RESPONSE_CLOSE:
            self.destroy()
        if response_code==100:
            self.gl_prop_editor.update()    

    def rebuild_gl_object_tree(self):
        """Rebuilds the GLObject tree view.
        """
        self.gl_tree_ctrl.rebuild_gl_object_tree()

    def gl_tree_ctrl_selected(self, glo):
        """Callback invoked by the GLPropertyTreeControl when a
        GLObject is selected
        """
        ## remove the current property editor if there is one
        if self.gl_prop_editor!=None:
            if self.gl_prop_editor.gl_object==glo:
                return

            self.hpaned.remove(self.gl_prop_editor)
            self.gl_prop_editor = None

        ## create new property editor
        self.gl_prop_editor = GLPropertyEditor(glo)
        self.hpaned.add2(self.gl_prop_editor)
        self.gl_prop_editor.show_all()

    def select_gl_object(self, glo):
        """Selects the given GLObject for editing by
        simulating a selection through the GLPropertyTreeControl.
        """
        self.gl_tree_ctrl.select_gl_object(glo)


###############################################################################
### TLS Analysis Dialog 
###

class TLSDialog(gtk.Dialog):
    """Dialog for the visualization and analysis of TLS model parameters
    describing rigid body motion of a structure.  The TLS model parameters
    can either be extracted from PDB REMARK statements, or loaded from a
    CCP4/REFMAC TLSIN file.
    """    
    def __init__(self, **args):
        self.main_window = args["main_window"]
        self.sc = args["struct_context"]
        self.selected_tls   = None
        self.animation_time = 0.0
        self.animation_list = []

        ## map chain_id -> color index so each tls group in a chain
        ## gets its own color
        self.chn_colori     = {}
        self.gl_tls_chain   = {}
        ## master list of tls groups handled by the dialog
        self.tls_list       = []

        gtk.Dialog.__init__(
            self,
            "TLS Analysis: %s" % (str(self.sc.struct)),
            self.main_window.window,
            gtk.DIALOG_DESTROY_WITH_PARENT)

        self.add_button("Open TLSOUT", 100)
        self.add_button("Save TLSOUT", 101)
        self.add_button("Fit TLS Groups", 102)
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
            gobject.TYPE_BOOLEAN,   #1
            gobject.TYPE_STRING,    #2
            gobject.TYPE_STRING,    #3
            gobject.TYPE_STRING,    #4
            gobject.TYPE_STRING)    #5
     
        treeview = gtk.TreeView(self.model)
        sw.add(treeview)
        treeview.connect("row-activated", self.row_activated_cb)
        treeview.connect("button-release-event", self.button_release_event_cb)
        treeview.set_rules_hint(gtk.TRUE)

        cell_rend = gtk.CellRendererToggle()
        column = gtk.TreeViewColumn("Show", cell_rend)
        column.add_attribute(cell_rend, "active", 0)
        treeview.append_column(column)
        cell_rend.connect("toggled", self.view_toggled)

        cell_rend = gtk.CellRendererToggle()
        column = gtk.TreeViewColumn("Animate", cell_rend)
        column.add_attribute(cell_rend, "active", 1)
        treeview.append_column(column)
        cell_rend.connect("toggled", self.animate_toggled)

        cell_rend = gtk.CellRendererText()
        column = gtk.TreeViewColumn("TLS Group", cell_rend)
        column.add_attribute(cell_rend, "markup", 2)
        treeview.append_column(column)

        cell_rend = gtk.CellRendererText()
        column = gtk.TreeViewColumn("", cell_rend)
        column_label = gtk.Label("")
        column_label.set_markup("trace(<b>T</b>) A<sup>2</sup>")
        column_label.show()
        column.set_widget(column_label)
        column.add_attribute(cell_rend, "markup", 3)
        treeview.append_column(column)

        cell_rend = gtk.CellRendererText()
        column = gtk.TreeViewColumn("", cell_rend)
        column.add_attribute(cell_rend, "markup", 4)
        treeview.append_column(column)
        column_label = gtk.Label("")
        column_label.set_markup("trace(<b>L</b>) DEG<sup>2</sup>")
        column_label.show()
        column.set_widget(column_label)

        cell_rend = gtk.CellRendererText()
        column = gtk.TreeViewColumn("Source", cell_rend)
        column.add_attribute(cell_rend, "markup", 5)
        treeview.append_column(column)

        self.show_all()

        gobject.timeout_add(200, self.timeout_cb)

        self.load_PDB(self.sc.struct.path)

    def error_dialog(self, text):
        self.main_window.error_dialog(text)

    def response_cb(self, dialog, response_code):
        """Responses to dialog events.
        """
        if response_code==gtk.RESPONSE_CLOSE:
            self.destroy()
        elif response_code==100:
            self.load_TLSOUT_cb()
        elif response_code==101:
            self.save_TLSOUT_cb()
        elif response_code==102:
            self.load_TLS_fit()

    def load_TLSOUT_cb(self):
        """Loads a TLSOUT File -- callback handler for button.
        """
        file_sel = gtk.FileSelection("Select TLSOUT File To Load")
        file_sel.hide_fileop_buttons()
        response = file_sel.run()

        if response==gtk.RESPONSE_OK:
            path = file_sel.get_filename()
            if path!=None and path!="":
                self.load_TLSOUT(path)

        file_sel.destroy()

    def save_TLSOUT_cb(self):
        """Saves a TLSOUT File -- callback handler for button.
        """
        file_sel  = gtk.FileSelection("Save TLSOUT File")
        response  = file_sel.run()        
        save_path = file_sel.get_filename()
        file_sel.destroy()

        if not save_path or response!=gtk.RESPONSE_OK:
            return

        tls_file = TLSFile()
        tls_file.set_file_format(TLSFileFormatTLSOUT())

        for tls in self.tls_list:
            tls_desc = tls["tls_desc"]
            tls_file.tls_desc_list.append(tls_desc)
            tls_desc.set_name(tls["name"])
            tls_desc.set_tls_group(tls["tls_group"])

        try:
            tls_file.save(open(save_path, "w"))
        except IOError, err:
            self.error_dialog("Error Saving File: %s" % (str(err)))
        
    def destroy_cb(self, *args):
        """Destroy the TLS dialog and everything it has built
        in the GLObject viewer.
        """
        self.clear_tls_groups()

    def row_activated_cb(self, tree_view, path, column):
        """ Use the GLPropertyBrowserDialog associated with the
        application main window to present the properties of the
        gl_tls object for modification.
        """
        index = int(path[0])
        selected_tls = self.tls_list[index]
        self.main_window.autoselect_gl_prop_browser(selected_tls["GLTLSGroup"])

    def button_release_event_cb(self, tree_view, bevent):
        x = int(bevent.x)
        y = int(bevent.y)

        try:
            (path, col, x, y) = tree_view.get_path_at_pos(x, y)
        except TypeError:
            return gtk.FALSE

        self.row_activated_cb(tree_view, path, col)
        return gtk.FALSE

    def view_toggled(self, cell, path):
        """Visible/Hidden TLS representation.
        """
        ## why, oh why?? (is this passed as a string)
        path = int(path)

        ## get toggled iter
        miter = self.model.get_iter((path,))

        show_vis = self.model.get_value(miter, 0)
        tls = self.tls_list[path]
    
        # do something with the value
        if show_vis==gtk.TRUE:
            tls["GLTLSGroup"].glo_update_properties(visible=False)
        else:
            tls["GLTLSGroup"].glo_update_properties(visible=True)

    def animate_toggled(self, cell, path):
        """Start/Stop TLS Animation.
        """
        index = int(path)
        tls   = self.tls_list[index]

        miter   = self.model.get_iter((index,))
        animate = self.model.get_value(miter, 1)

        if animate==gtk.FALSE:
            self.model.set(miter, 1, gtk.TRUE)
            self.animation_list.append(tls)
        elif animate==gtk.TRUE:
            self.model.set(miter, 1, gtk.FALSE)
            self.animation_list.remove(tls)

    def add_tls_group(self, tls):
        """Adds the TLS group and creates tls.gl_tls OpenGL
        renderer.
        """
        self.tls_list.append(tls)
        tls_group = tls["tls_group"]

        ## if no atoms were found in the tls_group,
        ## then do not add it
        try:
            atm0 = tls_group[0]
        except IndexError:
            return

        ## figure out the TLS group chain_id
        chain_id = atm0.chain_id
        tls["chain_id"] = chain_id

        ## if the TLS group is NULL, then perform a LSQ fit of it
        if tls_group.is_null():
            tls_group.calc_TLS_least_squares_fit()
            tls["lsq_fit"] = True
        tls["tls_info"] = tls_group.calc_tls_info()

        ## select color for the TLS group
        if tls.has_key("color_rgbf"):
            tls["color_name"] = "%f,%f,%f" % tls["color_rgbf"]
        else:
            ## set the tls group color starting from color 2, unique
            ## for each chain
            try:
                self.chn_colori[chain_id] += 1
                colori = self.chn_colori[chain_id]
            except KeyError:
                self.chn_colori[chain_id] = colori = 2

            ## calculate importent properties for the TLS group
            tls["color_name"] = COLOR_NAMES_CAPITALIZED[colori]

        ## creat GLTLSGroup for the visualization component
        tls["GLTLSGroup"] = GLTLSGroup(
            tls_group = tls["tls_group"],
            tls_info  = tls["tls_info"],
            tls_name  = tls["name"],
            tls_color = tls["color_name"])

        ## XXX: hack
        if tls.has_key("goodness"):
            trad = tls["GLTLSGroup"].gl_atom_list.properties["trace_radius"]
            trad += tls["goodness"] * 0.1
            tls["GLTLSGroup"].gl_atom_list.properties.update(trace_radius=trad)

        ## get the GLTLSChain to add the GLTLSGroup to
        try:
            gl_tls_chain = self.gl_tls_chain[chain_id]
        except KeyError:
            gl_tls_chain = GLTLSChain(chain_id=chain_id)
            self.sc.gl_struct.glo_add_child(gl_tls_chain)
            self.gl_tls_chain[chain_id] = gl_tls_chain

        tls["GLTLSChain"] = gl_tls_chain
        gl_tls_chain.add_gl_tls_group(tls["GLTLSGroup"])
        tls["GLTLSGroup"].glo_add_update_callback(self.update_cb)

        ## redraw the treeview adds this tls group
        self.redraw_treeview()

        ## rebuild GUI viewers
        tab = self.main_window.get_sc_tab(self.sc)
        if tab.has_key("gl_prop_browser"):
            tab["gl_prop_browser"].rebuild_gl_object_tree()

    def clear_tls_groups(self):
        """Remove the current TLS groups, including destroying
        the tls.gl_tls OpenGL renderer
        """
        gl_viewer = self.sc.gl_struct.glo_get_root()

        ## remove the individual GLTLSGroup visualization objects
        for tls in self.tls_list:
            gl_viewer.glv_remove_draw_list(tls["GLTLSGroup"])
            tls["GLTLSGroup"].glo_remove_update_callback(self.update_cb)
            del tls["GLTLSGroup"]

        ## remove the GLTLSChain objects
        for gl_tls_chain in self.gl_tls_chain.values():
            gl_viewer.glv_remove_draw_list(gl_tls_chain)

        ## re-initalize
        self.gl_tls_chain   = {}
        self.chn_colori     = {}
        self.tls_list       = []
        self.animation_list = []

        ## rebuild GUI viewers
        tab = self.main_window.get_sc_tab(self.sc)
        if tab.has_key("gl_prop_browser"):
            tab["gl_prop_browser"].rebuild_gl_object_tree()
        
        self.redraw_treeview()

    def update_cb(self, updates={}, actions=[]):
        """Property change callback from the GLTLSGroups.
        """
        i = 0
        for tls in self.tls_list:
            miter = self.model.get_iter((i,))

            if tls["GLTLSGroup"].properties["visible"]==True:
                self.model.set(miter, 0, gtk.TRUE)
            else:
                self.model.set(miter, 0, gtk.FALSE)

            i += 1

    def load_PDB(self, path):
        """Load TLS descriptions from PDB REMARK records.
        """
        self.clear_tls_groups()
        
        tls_file = TLSFile()
        tls_file.set_file_format(TLSFileFormatPDB())

        try:
            tls_file.load(open(path, "r"), path)
        except IOError:
            self.error_dialog("File Not Found: %s" % (path))
            return
        except TLSFileFormatError:
            self.error_dialog("File Format Error: %s" % (path))
            return

        if len(tls_file.tls_desc_list)==0:
            self.error_dialog("No TLS Groups Found: %s" % (path))
            return
        
        for tls_desc in tls_file.tls_desc_list:
            tls = {}
            tls["pdb_path"]  = path
            tls["tls_desc"]  = tls_desc
            tls["tls_group"] = tls_desc.generate_tls_group(self.sc.struct)
            tls["name"]      = self.markup_tls_name(tls_desc)
            self.add_tls_group(tls)
            
    def load_TLSOUT(self, path):
        """Load TLS descriptions from a REMAC/CCP4 TLSOUT file.
        """
        self.clear_tls_groups()

        tls_file = TLSFile()
        tls_file.set_file_format(TLSFileFormatTLSOUT())

        try:
            tls_file.load(open(path, "r"), path)
        except IOError:
            self.error_dialog("File Not Found: %s" % (path))
            return
        except TLSFileFormatError:
            self.error_dialog("File Format Error: %s" % (path))
            return

        if len(tls_file.tls_desc_list)==0:
            self.error_dialog("No TLS Groups Found: %s" % (path))
            return
        
        for tls_desc in tls_file.tls_desc_list:
            tls = {}
            tls["tlsout_path"] = path
            tls["tls_desc"]    = tls_desc
            tls["tls_group"]   = tls_desc.generate_tls_group(self.sc.struct)
            tls["name"]        = self.markup_tls_name(tls_desc)
            self.add_tls_group(tls)

    def load_TLS_fit(self):
        """Fit TLS groups to sequence segments
        """
        self.clear_tls_groups()

        dialog = TLSSearchDialog(
            main_window    = self.main_window,
            struct_context = self.sc)

        dialog.run()
        dialog.destroy()

        ## no groups!
        if len(dialog.tls_info_list)==0:
            return

        ## create a tls_list containing the tls description
        ## dictionaries for all fit TLS segments
        min_dp2 = None
        max_dp2 = None

        tls_list = []
        
        for tls_info in dialog.tls_info_list:
            tls = {}
            tls["tls_group"] = tls_info["tls_group"]
            tls["tls_info"]  = tls_info
            tls["tls_desc"]  = tls_info["tls_desc"]
            tls["name"]      = tls_info["name"]
            tls["TLS_fit"]   = True
            tls["lsq_fit"]   = True

            tls_list.append(tls)

            ## find min/max dp2 for coloring
            if min_dp2==None:
                min_dp2 = max_dp2 = tls_info["mean_dp2"]
            else:
                min_dp2 = min(min_dp2, tls_info["mean_dp2"])
                max_dp2 = max(max_dp2, tls_info["mean_dp2"])

        ## color tls groups according to dp2 goodness of fit
        rng = max_dp2 - min_dp2
        for tls in tls_list:
            dp2 = tls["tls_info"]["mean_dp2"]
            gof   = 1.0 - (dp2 - min_dp2)/rng
            tls["goodness"]   = gof
            tls["color_rgbf"] = goodness_color(gof)
            #print x, tls["color_rgbf"]

        for tls in tls_list:
            self.add_tls_group(tls)

    def markup_tls_name(self, tls_desc):
        listx = []
        for (chain_id1,frag_id1,chain_id2,frag_id2,sel) in tls_desc.range_list:
            listx.append("%s%s-%s%s %s" % (
                chain_id1, frag_id1, chain_id2, frag_id2, sel))
        return string.join(listx, "\n")

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

        for tls in self.tls_list:
            miter = self.model.append(None)

            tls_info = tls["tls_info"]
            tr_rT = "%8.4f" % (trace(tls_info["rT'"]))
            tr_L  = "%8.4f" % (trace(tls_info["L'"]*RAD2DEG2))

            if tls["GLTLSGroup"].properties["visible"]==True:
                self.model.set(miter, 0, gtk.TRUE)
            else:
                self.model.set(miter, 0, gtk.FALSE)

            if tls in self.animation_list:
                self.model.set(miter, 1, gtk.TRUE)
            else:
                self.model.set(miter, 1, gtk.FALSE)

            self.model.set(miter, 2, "<small>%s</small>" % (tls["name"]))
            self.model.set(miter, 3, tr_rT)
            self.model.set(miter, 4, tr_L)

            if tls.has_key("pdb_path"):
                source = "PDB File"
                
            elif tls.has_key("tlsout_path"):
                if tls.get("lsq_fit", False):
                    source = "TLSIN/LSQ Fit"
                else:
                    source = "TLSOUT"
                    
            elif tls.has_key("TLS_fit"):
                source = "LSQ Fit"

            self.model.set(miter, 5, source)
        
    def timeout_cb(self):
        """Timer which drives the TLS animation.
        """
        slices = 20
        slice  = 1.0 / slices 

        self.animation_time = (self.animation_time + slice) % 1.0
        
        for tls in self.animation_list:
            tls["GLTLSGroup"].properties.update(time=self.animation_time)
        return gtk.TRUE


class TLSSearchDialog(gtk.Dialog):
    """Dialog interface for the TLS Group Fitting algorithm.
    """    
    def __init__(self, **args):
        self.main_window    = args["main_window"]
        self.sc             = args["struct_context"]
        self.sel_tls_group  = None
        self.tls_info_list  = []
        self.cancel_flag    = False

        gtk.Dialog.__init__(
            self,
            "TLS Search: %s" % (str(self.sc.struct)),
            self.main_window.window,
            gtk.DIALOG_DESTROY_WITH_PARENT)


        self.cancel  = gtk.Button(stock=gtk.STOCK_CANCEL)
        self.back    = gtk.Button(stock=gtk.STOCK_GO_BACK)
        self.forward = gtk.Button(stock=gtk.STOCK_GO_FORWARD)
        self.ok      = gtk.Button(stock=gtk.STOCK_OK)

        self.add_action_widget(self.cancel,  gtk.RESPONSE_CANCEL)
        self.action_area.add(self.back)
        self.action_area.add(self.forward)
        self.add_action_widget(self.ok,      gtk.RESPONSE_OK)

        self.set_default_size(400, 400)

        self.back.connect("clicked", self.back_clicked)
        self.forward.connect("clicked", self.forward_clicked)
        
        self.connect("destroy",  self.destroy_cb)
        self.connect("response", self.response_cb)

        self.build_gui()
        self.switch_page1()
        self.show_all()

    def switch_page1(self):
        self.back.set_sensitive(gtk.FALSE)
        self.forward.set_sensitive(gtk.TRUE)
        self.ok.set_sensitive(gtk.FALSE)
        self.notebook.set_current_page(0)

    def switch_page2(self):
        self.back.set_sensitive(gtk.FALSE)
        self.forward.set_sensitive(gtk.FALSE)
        self.ok.set_sensitive(gtk.FALSE)
        self.notebook.set_current_page(1)

    def switch_page3(self):
        self.back.set_sensitive(gtk.TRUE)
        self.forward.set_sensitive(gtk.FALSE)
        self.ok.set_sensitive(gtk.TRUE)
        self.notebook.set_current_page(2)

    def build_gui(self):
        self.notebook = gtk.Notebook()
        self.vbox.pack_start(self.notebook, gtk.TRUE, gtk.TRUE, 0)
        self.notebook.set_show_tabs(gtk.FALSE)
        self.notebook.set_show_border(gtk.FALSE)

        self.notebook.append_page(self.build_page1(), gtk.Label(""))
        self.notebook.append_page(self.build_page2(), gtk.Label(""))
        self.notebook.append_page(self.build_page3(), gtk.Label(""))

    def build_page1(self):
        ## Selection Table
        frame = gtk.Frame("Rigid Body Search Criteria")
        frame.set_border_width(5)

        ## compute grid size of the table
        table_rows = self.sc.struct.count_chains() + 6
        
        table = gtk.Table(2, table_rows, gtk.FALSE)
        frame.add(table)
        table.set_border_width(5)
        table.set_row_spacings(5)
        table.set_col_spacings(10)

        current_row = 0
        def attach1(_widget, _cr):
            table.attach(_widget, 0, 2, _cr, _cr+1,
                         gtk.FILL|gtk.EXPAND, 0, 0, 0)
        
        ## Chain Checkbuttons
        self.chain_checks = {}

        for chain in self.sc.struct.iter_chains():
            chain_check = gtk.CheckButton(
                "Include Chain %s" % (chain.chain_id))
            self.chain_checks[chain.chain_id] = chain_check

            if chain.count_amino_acids()>10:
                chain_check.set_active(gtk.TRUE)
            else:
                chain_check.set_active(gtk.FALSE)

            align = gtk.Alignment(0.0, 0.5, 1.0, 0.0)
            align.add(chain_check)
            attach1(align, current_row)
            current_row += 1
            
        ##
        self.include_side_chains = gtk.CheckButton("Include Side Chain Atoms")
        align = gtk.Alignment(0.0, 0.5, 1.0, 0.0)
        align.add(self.include_side_chains)
        attach1(align, current_row)
        current_row += 1
        self.include_side_chains.set_active(gtk.TRUE)

        ##
        self.include_disordered = gtk.CheckButton("Include Disordered Atoms")
        align = gtk.Alignment(0.0, 0.5, 1.0, 0.0)
        align.add(self.include_disordered)
        attach1(align, current_row)
        current_row += 1
        self.include_disordered.set_active(gtk.FALSE)

        ##
        self.include_single_bond = gtk.CheckButton(
            "Include Atoms with One Bond")
        align = gtk.Alignment(0.0, 0.5, 1.0, 0.0)
        align.add(self.include_single_bond)
        attach1(align, current_row)
        current_row += 1
        self.include_single_bond.set_active(gtk.TRUE)

        ##
        label = gtk.Label("TLS Segment Residue Width")
        align = gtk.Alignment(0.0, 0.5, 0.0, 0.0)
        align.add(label)
        table.attach(align, 0, 1, current_row, current_row + 1,
                     gtk.FILL, 0, 0, 0)
        
        self.segment_width = gtk.Entry()
        align = gtk.Alignment(0.0, 0.5, 1.0, 0.0)
        align.add(self.segment_width)
        table.attach(align, 1, 2, current_row, current_row + 1,
                     gtk.EXPAND|gtk.FILL, 0, 0, 0)
        current_row += 1
        self.segment_width.set_text("6")

        return frame

    def build_page2(self):
        frame = gtk.Frame("Fitting TLS Groups")
        frame.set_border_width(5)

        self.page2_label = gtk.Label("")
        frame.add(self.page2_label)

        return frame

    def build_page3(self):
        frame = gtk.Frame("TLS Groups Fit")
        frame.set_border_width(5)

        sw = gtk.ScrolledWindow()
        frame.add(sw)
        sw.set_border_width(5)
        sw.set_shadow_type(gtk.SHADOW_ETCHED_IN)
        sw.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
     
        self.treeview = DictListTreeView(
            column_list=["Residue Range", "Atoms", "R", "<DP2>"])
        sw.add(self.treeview)

        return frame

    def response_cb(self, dialog, response_code):
        """Responses to dialog events.
        """
        if response_code==gtk.RESPONSE_CLOSE or \
           response_code==gtk.RESPONSE_CANCEL:
            self.tls_info_list = []
            self.cancel_flag   = True
            
        elif response_code==gtk.RESPONSE_OK:
            pass

    def back_clicked(self, button):
        page_num = self.notebook.get_current_page()
        if page_num==2:
            self.tls_info_list = []
            self.cancel_flag   = False
            self.switch_page1()

    def forward_clicked(self, button):
        page_num = self.notebook.get_current_page()
        if page_num==0:
            self.switch_page2()
            self.run_tls_analysis()
            self.switch_page3()

    def destroy_cb(self, *args):
        """Destroy the TLS dialog and everything it has built
        in the GLObject viewer.
        """
        self.cancel_flag = True
    
    def run_tls_analysis(self):
        """
        """
        ## get the list of chain_ids to fit
        chain_ids = []
        for chain_id, chain_check in self.chain_checks.items():
            if chain_check.get_active()==gtk.TRUE:
                chain_ids.append(chain_id)
        
        if self.include_side_chains.get_active()==gtk.TRUE:
            use_side_chains = True
        else:
            use_side_chains = False

        if self.include_disordered.get_active()==gtk.TRUE:
            include_frac_occupancy = True
        else:
            include_frac_occupancy = False

        if self.include_single_bond.get_active()==gtk.TRUE:
            include_single_bond = True
        else:
            include_single_bond = False
            
        try:
            residue_width = int(self.segment_width.get_text())
        except ValueError:
            residue_width = 6
        
        tls_analysis = TLSStructureAnalysis(self.sc.struct)

        self.tls_info_list = []

        for tls_info in tls_analysis.iter_fit_TLS_segments(
            chain_ids              = chain_ids,
            residue_width          = residue_width,
            use_side_chains        = use_side_chains,
            include_frac_occupancy = include_frac_occupancy,
            include_single_bond    = include_single_bond):

            ## calculations are canceled
            if self.cancel_flag==True:
                break

            self.tls_info_list.append(tls_info)

            ## set some dict values just for the treeview user interface
            ##"Residue Range", "Atoms", "R", "<DP2>", "sig<DP2>"

            tls_info["Residue Range"] = tls_info["name"]
            tls_info["R"]             = "%.3f" % (tls_info["R"])
            tls_info["<DP2>"]         = "%.4f" % (tls_info["mean_dp2"])
            tls_info["Atoms"]         = str(tls_info["num_atoms"])
            
            ## create a TLSGroupDesc object for this tls group
            tls_desc = TLSGroupDesc()
            tls_desc.set_name(tls_info["name"])
            tls_desc.set_tls_group(tls_info["tls_group"])
            tls_desc.add_range(
                tls_info["chain_id"], tls_info["frag_id1"],
                tls_info["chain_id"], tls_info["frag_id2"], "ALL")

            tls_info["tls_desc"] = tls_desc
            
            self.treeview.set_dict_list(self.tls_info_list)

            ## update status 
            ## this takes a long time so process some events
            self.page2_label.set_text("Fit %d Segments" % (
                len(self.tls_info_list)))

            while gtk.events_pending():
                gtk.main_iteration(gtk.TRUE)



###############################################################################
### Hierarchical GTK Treeview control for browsing a list of
### mmLib.Structure objects
###

class StructureTreeControl(gtk.TreeView):
    """
    """
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
        miter1 = self.model.append(None)
        self.model.set(miter1, 1, str(struct))

        for chain in struct.iter_chains():
            miter2 = self.model.append(miter1)
            self.model.set(miter2, 1, str(chain))

            for frag in chain.iter_fragments():
                miter3 = self.model.append(miter2)
                self.model.set(miter3, 1, str(frag))

                for atm in frag.iter_all_atoms():
                    miter4 = self.model.append(miter3)
                    self.model.set(miter4, 1, str(atm))
                    
    def append_struct(self, struct):
        self.struct_list.append(struct)
        self.display_struct(struct)

    def remove_struct(self, struct):
        self.struct_list.remove(struct)
        self.redraw()


class StructDetailsDialog(gtk.Dialog):
    """This was a old dialog used for debugging.
    """
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
        miter = self.store.append()
        self.store.set(miter, 0, key, 1, str(value))

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




###############################################################################
### About Dialog
### 
###


class AboutDialog(gtk.Dialog):
    """About Information...
    """
    def __init__(self, parent):
        gtk.Dialog.__init__(
            self,
            "About mmLib Viewer",
            parent,
            gtk.DIALOG_DESTROY_WITH_PARENT)
        
        self.add_button(gtk.STOCK_CLOSE, gtk.RESPONSE_CLOSE)
        self.connect("response", self.response)
        self.set_resizable(gtk.FALSE)

        frame = gtk.Frame()
        self.vbox.pack_start(frame, gtk.TRUE, gtk.TRUE, 0)
        frame.set_border_width(10)

        label = gtk.Label()
        frame.add(label)
        
        label.set_markup(
            '<span size="large">mmLib Viewer</span>\n'\
            'Ethan Merritt and Jay Painter\n'\
            'http://pymmlib.sourceforge.net/')

        self.show_all()


    def response(self, *args):
        self.destroy()
        

###############################################################################
### Application Window and Multi-Document tabs
### 
###

class StructureContext(object):
    """Keeps track of a loaded Structure's contextual information.
    """
    def __init__(self, struct = None, gl_struct = None):
        self.struct_id = ""
        self.tab_id    = ""
        self.struct    = struct
        self.gl_struct = gl_struct

    def suggest_struct_id(self):
        """Suggests a name for the display of the Structure.
        """
        entry_id = self.struct.cifdb.get_entry_id()
        if entry_id!=None and entry_id!="":
            return entry_id
        if hasattr(self.struct, "path"):
            path = getattr(self.struct, "path")
            dir, basename = os.path.split(path)
            return basename
        return "XXXX"


class MainWindow(object):
    """Main viewer window.
    """
    def __init__(self, quit_notify_cb):
        self.quit_notify_cb      = quit_notify_cb
        
        self.selected_sc         = None
        self.tab_list            = []

        ## dialog lists
        self.details_dialog_list = []
        self.tls_dialog_list     = []

        ## Create the toplevel window
        self.window = gtk.Window()
        self.window.set_title("TLSViewer")
        self.window.set_default_size(500, 400)
        self.window.connect('destroy', self.file_quit, self)

        table = gtk.Table(1, 4, gtk.FALSE)
        self.window.add(table)

        ## file menu bar
        menu_items = [
            ('/_File', None, None, 0,'<Branch>'),
            ('/File/_New Window', None, self.file_new_window,
             0,'<StockItem>', gtk.STOCK_NEW),
            ('/File/_New Tab', None, self.file_new_tab,
             0,'<StockItem>',gtk.STOCK_NEW),
            ('/File/sep1', None, None, 0,'<Separator>'),
            ('/File/_Open', None, self.file_open,
             0,'<StockItem>', gtk.STOCK_OPEN),
            ('/File/sep2', None, None, 0, '<Separator>'),
            ('/File/_Close Structure', None, self.file_close_struct,
             0,'<StockItem>',gtk.STOCK_CLOSE),
            ('/File/_Close Tab', None, self.file_close_tab,
             0,'<StockItem>',gtk.STOCK_CLOSE),
            ('/File/sep3' ,None, None, 0,'<Separator>'),
            ('/File/_Quit', '<control>Q', self.file_quit,
             0,'<StockItem>',gtk.STOCK_QUIT),
            
            ('/_Visualization', None, None, 0, '<Branch>'),
            ('/Visualization/_Properties Browser...', None,
             self.visualization_properties_browser, 0, None),

            ('/_TLS', None, None, 0, '<Branch>'),
            ('/TLS/TLS Analysis...', None, self.tools_tls_analysis, 0, None),

            ('/_Structures', None, None, 0, '<Branch>'),

            ('/_Help',       None, None, 0, '<Branch>'),
            ('/Help/_About', None, self.help_about, 0, None) ]

        self.accel_group = gtk.AccelGroup()
        self.window.add_accel_group(self.accel_group)

        self.item_factory = gtk.ItemFactory(
            gtk.MenuBar, '<main>', self.accel_group)
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

        ## Notebook
        self.notebook = gtk.Notebook()
        table.attach(self.notebook,
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

    def file_new_window(self, *args):
        """File->New Window
        """
        VIEWER_APP.new_window()

    def file_new_tab(self, *args):
        """File->New Tab
        """
        self.add_tab()
    
    def file_open(self, *args):
        """File->Open
        """
        if hasattr(self, "file_selector"):
            self.file_selector.present()
            return

        ## create file selector
        self.file_selector = gtk.FileSelection("Select file to view");

        ok_button     = self.file_selector.ok_button
        cancel_button = self.file_selector.cancel_button

        ok_button.connect("clicked", self.open_ok_cb, None)
        cancel_button.connect("clicked", self.open_cancel_cb, None)        

        self.file_selector.present()
    
    def destroy_file_selector(self):
        """
        """
        self.file_selector.destroy()
        del self.file_selector

    def open_ok_cb(self, *args):
        """
        """
        path = self.file_selector.get_filename()
        self.destroy_file_selector()
        self.load_file(path)

    def open_cancel_cb(self, *args):
        """
        """
        self.destroy_file_selector()

    def file_close_struct(self, *args):
        """File->Close Structure
        Closes the currently selected structure.
        """
        if self.selected_sc!=None:
            self.remove_sc(self.selected_sc)
        
    def file_close_tab(self, *args):
        """File->Close Tab
        Closes the current viewable tab.
        """
        tab = self.get_current_tab()
        if tab==None:
            return

        ## remove all structures 
        for sc in tab["sc_list"]:
            self.remove_sc(sc)

        ## remove the page from the notebook
        page_no = self.notebook.page_num(tab["gl_viewer"])
        self.notebook.remove_page(page_no)

        ## destroy the GLPropertiesDialog for the tab
        if tab.has_key("gl_prop_browser"):
            tab["gl_prop_browser"].destroy()

        self.tab_list.remove(tab)

    def file_quit(self, *args):
        """File->Quit
        """
        self.quit_notify_cb(self.window, self)

    def visualization_properties_browser(self, *args):
        """Edit->Properties (GLPropertyBrowserDialog)
        """
        tab = self.get_current_tab()
        if tab!=None:
            self.present_gl_prop_browser(tab)

    def tools_tls_analysis(self, *args):
        """Tools->TLS Analysis
        Launches the TLS Analysis dialog on the currently selected
        StructureContext.
        """
        if self.selected_sc==None:
            return
        
        tls = TLSDialog(
            main_window    = self,
            struct_context = self.selected_sc)
        
        tls.present()

    def help_about(self, *args):
        """Help->About
        """
        AboutDialog(self.window)

    def set_statusbar(self, text):
        """Sets the text in the status bar of the window.
        """
        self.statusbar.pop(0)
        self.statusbar.push(0, text)

    def update_cb(self, percent):
        """Callback for file loading code to inform the GUI of how
        of the file has been read
        """
        self.progress.set_fraction(percent/100.0)
        while gtk.events_pending():
            gtk.main_iteration(gtk.TRUE)
            
    def error_dialog(self, text):
        """Display modeal error dialog box containing the error text.
        """
        dialog = gtk.MessageDialog(
            self.window,
            gtk.DIALOG_DESTROY_WITH_PARENT,
            gtk.MESSAGE_ERROR,
            gtk.BUTTONS_CLOSE,
            text)
        
        dialog.run()
        dialog.destroy()

    def add_tab(self):
        """Adds a tab page to the tabbed-viewer notebook and places a
        GtkGLViewer widget inside of it.
        """
        ## initalize new tab dictionary
        tab = {}
        tab["sc_list"] = []
        
        ## come up with a new tab name
        tab["num"] = 1
        
        for tabx in self.tab_list:
            if tabx["num"]>=tab["num"]:
                tab["num"] = tabx["num"] + 1

        ## now add the tab dictionary
        self.tab_list.append(tab)

        tab["name"]      = "Tab %d" % (tab["num"])
        tab["gl_viewer"] = GtkGLViewer()

        page_num  = self.notebook.append_page(
            tab["gl_viewer"], gtk.Label(tab["name"]))

        tab["gl_viewer"].show()
        
        return tab

    def get_current_tab(self):
        """Returns the tab dictionary of the current viewable tab.
        """
        page_num = self.notebook.get_current_page()
        if page_num==-1:
            return None

        gl_viewer = self.notebook.get_nth_page(page_num)
        for tab in self.tab_list:
            if tab["gl_viewer"]==gl_viewer:
                return tab

        return None

    def get_sc_tab(self, sc):
        """Return the tab containing the argument StructureContext.
        """
        for tab in self.tab_list:
            for scx in tab["sc_list"]:
                if sc==scx:
                    return tab
        return None

    def remove_sc(self, sc):
        """Completely removes the Structure (referenced by its sc).
        """
        ## if this StructureContext is the selected context, unselect it
        ## before removing
        if self.selected_sc==sc:
            self.set_selected_sc(None)
        
        tab = self.get_sc_tab(sc)
        tab["sc_list"].remove(sc)

        ## remove from GLViewer
        tab["gl_viewer"].glv_remove_draw_list(sc.gl_struct)

        ## remove from /Structure menu
        structs = self.item_factory.get_item("/Structures")
        structs_menu = structs.get_submenu()

        for mi in structs_menu.get_children():
            label = mi.get_child()
            text  = label.get_text()

            if text==sc.struct_id:
                structs_menu.remove(mi)
                mi.destroy()
                break

        ## rebuild the brower for the tab
        if tab.has_key("gl_prop_browser"):
            tab["gl_prop_browser"].rebuild_gl_object_tree()

    def add_struct(self, struct, new_tab=False, new_window=False):
        """Adds the Structure to the window.  Returns the newly created
        StructureContext for the window.
        """
        ## create a StructureContext used by the viewer for
        ## this Structure/GLStructure pair
        sc = StructureContext()

        sc.struct    = struct
        sc.struct_id = sc.suggest_struct_id()
        
        ## get the current notebook tab's GtkGLViewer or create
        ## a new notebook table with a new GtkGLViewer
        if new_tab==True:
            tab = self.add_tab()
        else:
            tab = self.get_current_tab()
            if tab==None:
                tab = self.add_tab()

        ## add the StructureContext to the tab's sc_list
        tab["sc_list"].append(sc)

        ## add the structure to the GLViewer
        sc.gl_struct = tab["gl_viewer"].glv_add_struct(sc.struct)

        ## add the structure to the Structures menu
        structs = self.item_factory.get_item("/Structures")
        structs_menu = structs.get_submenu()

        mi = gtk.CheckMenuItem(sc.struct_id)
        structs_menu.append(mi)
        mi.sc  = sc
        mi.cid = mi.connect("activate", self.struct_item_activate)
        mi.show_all()

        ## rebuild the brower for the tab
        if tab.has_key("gl_prop_browser"):
            tab["gl_prop_browser"].rebuild_gl_object_tree()
            
        return sc

    def load_file(self, path, new_tab=False, new_window=False):
        """Loads the structure file specified in the path.
        """
        self.set_statusbar("Loading: %s" % (path))

        try:
            struct = LoadStructure(
                fil              = path,
                update_cb        = self.update_cb,
                build_properties = ("library_bonds","distance_bonds"))

        except IOError:
            self.set_statusbar("")
            self.progress.set_fraction(0.0)
            self.error_dialog("File Not Found: %s" % (path))
            return

        struct.path = path

        sc = self.add_struct(struct, new_tab, new_window)

        ## autoselect a loaded structure
        self.set_selected_sc(sc)

        ## blank the status bar
        self.set_statusbar("")
        self.progress.set_fraction(0.0)

        return gtk.FALSE

    def struct_item_activate(self, menu_item):
        """Structure->MI
        Called by the Structure menu items when activated.
        """
        self.set_selected_sc(menu_item.sc)

    def set_selected_sc(self, sc):
        """Set the selected StructureContext.
        """
        self.selected_sc = sc

        ## special handling if setting the current StructureContext
        ## to None
        if sc==None:

            ## uncheck all Structure menu items
            structs = self.item_factory.get_item("/Structures")
            structs_menu = structs.get_submenu()
            for mi in structs_menu.get_children():
                mi.disconnect(mi.cid)
                mi.set_active(gtk.FALSE)
                mi.cid = mi.connect("activate", self.struct_item_activate)

            ## blank the label
            self.select_label.set_text("")

            return

        ## set the Structure menu
        structs = self.item_factory.get_item("/Structures")
        structs_menu = structs.get_submenu()

        for mi in structs_menu.get_children():

            label = mi.get_child()
            text  = label.get_text()

            ## XXX: stupid GTK API emits activate signal causing recursion
            mi.disconnect(mi.cid)

            if text==sc.struct_id:
                mi.set_active(gtk.TRUE)
            else:
                mi.set_active(gtk.FALSE)

            ## XXX: reinstall signal handler
            mi.cid = mi.connect("activate", self.struct_item_activate)

        ## set the window bar information
        self.select_label.set_text(sc.struct_id)

        ## switch to the right notebook tab
        for gl_viewer in self.notebook.get_children():
            for glo in gl_viewer.glo_iter_children():
                if glo==sc.gl_struct:
                    page_no = self.notebook.page_num(gl_viewer)
                    self.notebook.set_current_page(page_no)

    def get_selected_sc(self):
        """Returns the currently selected StructureContext
        """
        return self.selected_sc

    def present_gl_prop_browser(self, tab):
        """Creates and/or presents the GLPropertyBrowserDialog for the
        given tab.
        """
        try:
            tab["gl_prop_browser"].present()

        except KeyError:
            tab["gl_prop_browser"] = GLPropertyBrowserDialog(
                parent_window = self.window,
                glo_root      = tab["gl_viewer"],
                title         = tab["name"])

            tab["gl_prop_browser"].connect(
                "destroy", self.gl_prop_browser_destroy, tab)

            tab["gl_prop_browser"].present()

    def autoselect_gl_prop_browser(self, gl_object):
        """Opens/Presents the GLPropertyBrowser Dialog and selects
        gl_object for editing.
        """

        ## the tab containing the gl_object needs to be found
        tab = None
        
        for tabx in self.tab_list:
            ## if the selected gl_object is the gl_viewer
            if tabx["gl_viewer"]==gl_object:
                tab = tabx
                break

            ## iterate the gl_viewer's children looking for the gl_object
            for glo in tabx["gl_viewer"].glo_iter_preorder_traversal():
                if glo==gl_object:
                    tab = tabx
                    break

        self.present_gl_prop_browser(tab)
        tab["gl_prop_browser"].select_gl_object(gl_object)

    def gl_prop_browser_destroy(self, widget, tab):
        """Callback when the GLProeriesBrowser is destroyed.
        """
        del tab["gl_prop_browser"]


class ViewerApp(object):
    """Viewer application manages the windows.  One
    one of these objects is created and used globally.
    """
    def __init__(self, first_path=None):
        self.window_list = []
        self.new_window(first_path)

    def main(self):
        """Runs the GTK mainloop until all windows are closed.
        """
        gtk.main()
        
    def new_window(self, path=None):
        """Creates a new window, and loads the structure at
        the given file path if it is given.
        """
        mw = MainWindow(self.window_quit_notify)
        self.window_list.append(mw)

        if path!=None:
            gobject.timeout_add(25, mw.load_file, path)

    def window_quit_notify(self, window, mw):
        """Callback when a window exits.
        """
        self.window_list.remove(mw)
        if len(self.window_list)==0:
            gtk.main_quit()


### <MAIN>
if __name__=="__main__":
    try:
        first_path = sys.argv[1]
    except IndexError:
        first_path = None
    
    VIEWER_APP = ViewerApp(first_path)

    try:
        VIEWER_APP.main()
    except KeyboardInterrupt:
        sys.exit(1)
### </MAIN>
