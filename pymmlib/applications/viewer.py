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

from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

from mmLib.PDB            import PDBFile
from mmLib.Structure      import *
from mmLib.FileLoader     import LoadStructure, SaveStructure
from mmLib.GLViewer       import *
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
            iter = self.model.append(None)

            for i in range(len(self.column_list)):
                column_desc = self.column_list[i]
                self.model.set(iter, i, dictX.get(column_desc, ""))


###############################################################################
### GTK OpenGL Viewer Widget Using GtkGlExt/PyGtkGLExt
###
### see http://gtkglext.sourceforge.net/ for information on GtkGLExt and
### OpenGL binding in GTK 2.0
###

class GtkGLViewer(gtk.gtkgl.DrawingArea, GLViewer):
    def __init__(self):
        gtk.gtkgl.DrawingArea.__init__(self)
        GLViewer.__init__(self)

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

        self.connect('button_press_event',  self.gtk_glv_button_press_event)
        self.connect('motion_notify_event', self.gtk_glv_motion_notify_event)
        self.connect('realize',             self.gtk_glv_realize)
        self.connect('map',                 self.gtk_glv_map)
        self.connect('configure_event',     self.gtk_glv_configure_event)
        self.connect('expose_event',        self.gtk_glv_expose_event)
        self.connect('destroy',             self.gtk_glv_destroy)

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
        self.beginx = event.x
        self.beginy = event.y

    def gtk_glv_motion_notify_event(self, glarea, event):
        width     = glarea.allocation.width
        height    = glarea.allocation.height

        x = 0.0
        y = 0.0
        z = 0.0

        rotx = 0.0
        roty = 0.0
        rotz = 0.0
        
        if (event.state & gtk.gdk.BUTTON1_MASK):
            roty += 360.0 * ((event.x - self.beginx) / float(width)) 
            rotx += 360.0 * ((event.y - self.beginy) / float(height))

        elif (event.state & gtk.gdk.BUTTON2_MASK):
            z = 50.0 * ((event.y - self.beginy) / float(height))
            rotz += 360.0 * ((event.x - self.beginx) / float(width)) 

        elif (event.state & gtk.gdk.BUTTON3_MASK):
            x = 50.0 * ((event.x - self.beginx) / float(height))
            y = 50.0 * (-(event.y - self.beginy) / float(height))

        self.glv_translate(x, y, z)
        self.glv_rotate(rotx, roty, rotz)

        self.beginx = event.x
        self.beginy = event.y

        self.queue_draw()

    def gtk_glv_destroy(self, glarea):
        ## XXX: delete all opengl draw lists
        return gtk.TRUE

    def glv_redraw(self):
        self.queue_draw()


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


class ColorSelection(gtk.Combo):
    def __init__(self):
        gtk.Combo.__init__(self)

        self.color = None

        self.color_dict = {
            "white":   (1.0, 1.0, 1.0),
            "red":     (1.0, 0.0, 0.0),
            "green":   (0.0, 1.0, 0.0),
            "blue":    (0.0, 0.0, 1.0) }
        
        self.set_popdown_strings(
            ["Random",
             "Color By Atom",
             "White",
             "Red",
             "Green",
             "Blue" ])

        self.entry.connect("activate", self.activate, None)

    def activate(self, entry, junk):
        txt = self.entry.get_text()

        ## color format: r, g, b
        try:
            r, g, b = txt.split(",")
            color = (float(r), float(g), float(b))
        except ValueError:
            pass
        else:
            self.set_color(color)
            return

        ## color format: Random
        if txt=="Random":
            color = (random.random(),
                     random.random(),
                     random.random())
            self.set_color(color)
            return

        ## color format: keywords
        if txt=="Color By Atom":
            self.set_color(None)
            return

        ## color fomat: match color name
        try:
            color = self.color_dict[txt.lower()]
        except KeyError:
            pass
        else:
            self.set_color(color)
            return

    def match_color_name(self, color):
        if color==None:
            return "Color By Atom"

        for color_name, color_val in self.color_dict.items():
            if color==color_val:
                return color_name

        return "%3.2f,%3.2f,%3.2f" % (color[0], color[1], color[2])
    
    def set_color(self, color):
        self.color = color
        name       = self.match_color_name(color)
        self.entry.set_text(name)

    def get_color(self):
        self.activate(self.entry, None)
        return self.color


class GLAtomColorSelection(gtk.Combo):

    color_dict = {
        "white":   (1.0, 1.0, 1.0),
        "red":     (1.0, 0.0, 0.0),
        "green":   (0.0, 1.0, 0.0),
        "blue":    (0.0, 0.0, 1.0) }
        
    def __init__(self):
        gtk.Combo.__init__(self)
        self.color = None
        self.atm_color = None
        self.entry.connect("activate", self.activate, None)

    def activate(self, entry, junk):
        txt = self.entry.get_text()

        ## color format: r, g, b
        try:
            r, g, b = txt.split(",")
            color = (float(r), float(g), float(b))
        except ValueError:
            pass
        else:
            self.set_color_val(color)
            return

        ## color format: Random
        if txt=="random":
            color = (random.random(),
                     random.random(),
                     random.random())
            self.set_color_val(color)
            return

        ## color fomat: match color name
        try:
            color = self.color_dict[txt.lower()]
        except KeyError:
            pass
        else:
            self.set_color_val(color)
            return

        ## keyword color settings
        self.set_color_val(txt)

    def match_color_name(self, color):
        if type(color)==StringType:
            return color
        elif type(color)==TupleType:
            for color_name, color_val in self.color_dict.items():
                if color==color_val:
                    return color_name

            return "%3.2f,%3.2f,%3.2f" % (color[0], color[1], color[2])

        return "WTF??"
    
    def set_color_val(self, color):
        self.color     = color
        self.atm_color = GLAtomColor(color)

        ## go through the GLAtomColor enumeration list and add the
        ## color selection options
        self.set_popdown_strings(
            ["random"] + self.color_dict.keys() +\
            self.atm_color.settings.keys())
        
        name = self.match_color_name(color)
        self.entry.set_text(name)

    def set_gl_atom_color(self, gl_atom_color):
        self.set_color_val(gl_atom_color.setting)

    def get_gl_atom_color(self):
        self.activate(self.entry, None)
        return self.atm_color


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

    def markup_property_label(self, prop_desc, max_len = 20):
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
                self.markup_property_label(prop, max_len=40))

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

        ## COLOR
        elif prop["type"]=="color":
            widget = ColorSelection()
            
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

            elif prop_desc["type"]=="color":
                widget.set_color(self.gl_object.properties[name])

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

            elif prop["type"]=="color":
                update_dict[name] = widget.get_color()

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
        self.column = gtk.TreeViewColumn("OpenGL Renderer", cell)
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
        model, iter = self.get_selection().get_selected()
        if iter!=None:
            selected_glo = self.path_glo_dict[model.get_path(iter)]
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
        
        def redraw_recurse(glo, parent_iter):
            iter = self.model.append(parent_iter)
            path = self.model.get_path(iter)

            self.path_glo_dict[path] = glo
            
            self.model.set(iter, 0, glo.glo_name())

            for child in glo.glo_iter_children():
                redraw_recurse(child, iter)
    
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
        self.sel_tls_group = None

        self.animation_time = 0.0
        self.animation_list = []
        self.tls_group_list = []

        gtk.Dialog.__init__(
            self,
            "TLS Analysis: %s" % (str(self.sc.struct)),
            self.main_window.window,
            gtk.DIALOG_DESTROY_WITH_PARENT)

        self.add_button("Open TLSIN", 100)
        self.add_button("Graphics Properties", 101)
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
        column = gtk.TreeViewColumn("T(A^2)", cell_rend)
        column.add_attribute(cell_rend, "markup", 3)
        treeview.append_column(column)

        cell_rend = gtk.CellRendererText()
        column = gtk.TreeViewColumn("L(deg^2)", cell_rend)
        column.add_attribute(cell_rend, "markup", 4)
        treeview.append_column(column)

        cell_rend = gtk.CellRendererText()
        column = gtk.TreeViewColumn("S(A*deg)", cell_rend)
        column.add_attribute(cell_rend, "markup", 5)
        treeview.append_column(column)

        self.show_all()

        gobject.timeout_add(200, self.timeout_cb)

        self.load_PDB(self.sc.struct.path)

    def response_cb(self, dialog, response_code):
        """Responses to dialog events.
        """
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
            ## use the GLPropertyBrowserDialog associated with
            ## the application main window to present the
            ## properties of the gl_tls object for modification
            
            self.main_window.autoselect_gl_prop_browser(
                self.sel_tls_group.gl_tls)

        elif response_code==102:
            self.load_TLS_fit()

    def destroy_cb(self, *args):
        """Destroy the TLS dialog and everything it has built
        in the GLObject viewer.
        """
        self.clear_tls_groups()

    def get_tls_group(self, path):
        return self.tls_group_list[int(path)]

    def row_activated_cb(self, tree_view, path, column):
        row = int(path[0])
        self.sel_tls_group = self.tls_group_list[row]

    def view_toggled(self, cell, path):
        """Visible/Hidden TLS representation.
        """
        ## why, oh why??
        path = int(path)

        # get toggled iter
        iter = self.model.get_iter((path,))

        show_vis = self.model.get_value(iter, 0)
        tls_group = self.tls_group_list[path]
    
        # do something with the value
        if show_vis==gtk.TRUE:
            tls_group.gl_tls.glo_update_properties(visible=False)
        else:
            tls_group.gl_tls.glo_update_properties(visible=True)

    def animate_toggled(self, cell, path):
        """Start/Stop TLS Animation.
        """
        path      = int(path)
        tls_group = self.tls_group_list[path]

        iter      = self.model.get_iter((path,))
        animate   = self.model.get_value(iter, 1)

        if animate==gtk.FALSE:
            self.model.set(iter, 1, gtk.TRUE)
            self.animation_list.append(tls_group)
        elif animate==gtk.TRUE:
            self.model.set(iter, 1, gtk.FALSE)
            self.animation_list.remove(tls_group)

    def add_tls_group(self, tls_group):
        """Adds the TLS group and creates tls.gl_tls OpenGL
        renderer.
        """
        self.tls_group_list.append(tls_group)
        
        tls_group.gl_tls = GLTLSGroup(tls_group=tls_group)
        self.sc.gl_struct.glo_add_child(tls_group.gl_tls)
        tls_group.gl_tls.glo_add_update_callback(self.update_cb)

        ## rebuild GUI viewers
        tab = self.main_window.get_sc_tab(self.sc)

        if tab.has_key("gl_prop_browser"):
            tab["gl_prop_browser"].rebuild_gl_object_tree()

        self.redraw_treeview()

    def clear_tls_groups(self):
        """Remove the current TLS groups, including destroying
        the tls.gl_tls OpenGL renderer
        """
        gl_viewer = self.sc.gl_struct.glo_get_root()

        for tls_group in self.tls_group_list:
            gl_viewer.glv_remove_draw_list(tls_group.gl_tls)
            tls_group.gl_tls.glo_remove_update_callback(self.update_cb)
            del tls_group.gl_tls
        
        self.tls_group_list = []
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
        for tls in self.tls_group_list:
            iter = self.model.get_iter((i,))

            if tls.gl_tls.properties["visible"]==True:
                self.model.set(iter, 0, gtk.TRUE)
            else:
                self.model.set(iter, 0, gtk.FALSE)

            i += 1

    def load_PDB(self, path):
        """Load TLS descriptions from PDB REMARK records.
        """
        self.clear_tls_groups()
        
        tls_file = TLSInfoList()
        pdb_file = PDBFile()
        pdb_file.load_file(path)
        pdb_file.record_processor(tls_file)

        for tls_info in tls_file:
            tls_group = tls_info.make_tls_group(self.sc.struct)
            tls_group.tls_info = tls_info
            self.add_tls_group(tls_group)
            
    def load_TLSOUT(self, path):
        """Load TLS descriptions from a REMAC/CCP4 TLSOUT file.
        """
        self.clear_tls_groups()
        try:
            fil = open(path, "r")
        except IOError:
            return

        tls_file = TLSInfoList()
        tls_file.load_refmac_tlsout_file(fil)
        
        for tls_info in tls_file:
            tls_group = tls_info.make_tls_group(self.sc.struct)
            tls_group.tls_info = tls_info
            self.add_tls_group(tls_group)

    def load_TLS_fit(self):
        """Fit TLS groups to sequence segments
        """
        self.clear_tls_groups()

        dialog = TLSSearchDialog(
            main_window    = self.main_window,
            struct_context = self.sc)

        dialog.run()
        dialog.destroy()

        for stats in dialog.tls_stats_list:
            self.add_tls_group(stats["tls"])
        
    def markup_tls_name(self, tls_info):
        listx = []
        for (chain_id1,frag_id1,chain_id2,frag_id2,sel) in tls_info.range_list:
            listx.append("%s%s-%s%s %s" % (
                chain_id1, frag_id1, chain_id2, frag_id2, sel))
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

        for tls in self.tls_group_list:
            iter = self.model.append(None)

            if tls.gl_tls.properties["visible"]==True:
                self.model.set(iter, 0, gtk.TRUE)
            else:
                self.model.set(iter, 0, gtk.FALSE)

            if tls in self.animation_list:
                self.model.set(iter, 1, gtk.TRUE)
            else:
                self.model.set(iter, 1, gtk.FALSE)


            if hasattr(tls, "tls_info"):
                self.model.set(iter, 2, self.markup_tls_name(tls.tls_info))
            elif hasattr(tls, "name"):
                self.model.set(iter, 2, tls.name)
            else:
                self.model.set(iter, 2, "name here")

            self.model.set(iter, 3, self.markup_tensor(tls.T))
            self.model.set(iter, 4, self.markup_tensor(tls.L*rad2deg2))
            self.model.set(iter, 5, self.markup_tensor(tls.S*rad2deg))
        
    def timeout_cb(self):
        """Timer which drives the TLS animation.
        """
        self.animation_time += 0.005
        for tls_group in self.animation_list:
            tls_group.gl_tls.properties.update(time=self.animation_time)
        return gtk.TRUE


class TLSSearchDialog(gtk.Dialog):
    """Dialog interface for the TLS Group Fitting algorithm.
    """    
    def __init__(self, **args):
        self.main_window    = args["main_window"]
        self.sc = args["struct_context"]
        self.sel_tls_group  = None

        self.tls_stats_list = []

        gtk.Dialog.__init__(
            self,
            "TLS Search: %s" % (str(self.sc.struct)),
            self.main_window.window,
            gtk.DIALOG_DESTROY_WITH_PARENT)

        self.add_button("Accept", gtk.RESPONSE_OK)
        self.add_button("Cancel", gtk.RESPONSE_CANCEL)

        self.set_default_size(400, 400)
        
        self.connect("destroy", self.destroy_cb)
        self.connect("response", self.response_cb)

        ## Selection Table
        frame = gtk.Frame("Rigid Body Search Criteria")
        self.vbox.pack_start(frame, gtk.FALSE, gtk.FALSE, 0)
        frame.set_border_width(5)
        
        table = gtk.Table(2, 7, gtk.FALSE)
        frame.add(table)
        table.set_border_width(5)
        table.set_row_spacings(5)
        table.set_col_spacings(10)
        
        self.include_side_chains = gtk.CheckButton("Include Side Chain Atoms")
        align = gtk.Alignment(0.0, 0.5, 1.0, 0.0)
        align.add(self.include_side_chains)
        table.attach(align, 0, 2, 0, 1, gtk.FILL|gtk.EXPAND, 0, 0, 0)
        self.include_side_chains.set_active(gtk.TRUE)

        self.include_disordered = gtk.CheckButton("Include Disordered Atoms")
        align = gtk.Alignment(0.0, 0.5, 1.0, 0.0)
        align.add(self.include_disordered)
        table.attach(align, 0, 2, 1, 2, gtk.FILL|gtk.EXPAND, 0, 0, 0)
        self.include_disordered.set_active(gtk.FALSE)

        self.include_single_bond = gtk.CheckButton(
            "Include Atoms with One Bond")
        align = gtk.Alignment(0.0, 0.5, 1.0, 0.0)
        align.add(self.include_single_bond)
        table.attach(align, 0, 2, 2, 3, gtk.FILL|gtk.EXPAND, 0, 0, 0)
        self.include_single_bond.set_active(gtk.TRUE)
        
        label = gtk.Label("TLS Segment Residue Width")
        align = gtk.Alignment(0.0, 0.5, 0.0, 0.0)
        align.add(label)
        table.attach(align, 0, 1, 3, 4, gtk.FILL, 0, 0, 0)
        self.segment_width = gtk.Entry()
        align = gtk.Alignment(0.0, 0.5, 1.0, 0.0)
        align.add(self.segment_width)
        table.attach(align, 1, 2, 3, 4, gtk.EXPAND|gtk.FILL, 0, 0, 0)
        self.segment_width.set_text("6")
        
        label = gtk.Label("Maximum dP2 Value")
        align = gtk.Alignment(0.0, 0.5, 0.0, 0.0)
        align.add(label)
        table.attach(align, 0, 1, 4, 5, gtk.FILL, 0, 0, 0)
        self.dp2_entry = gtk.Entry()
        align = gtk.Alignment(0.0, 0.5, 1.0, 0.0)
        align.add(self.dp2_entry)
        table.attach(align, 1, 2, 4, 5, gtk.EXPAND|gtk.FILL, 0, 0, 0)
        self.dp2_entry.set_text("0.03")

        self.run_button = gtk.Button("Run Analysis")
        align = gtk.Alignment(1.0, 0.5, 0.0, 0.0) 
        align.add(self.run_button)
        table.attach(align, 0, 2, 5, 6, gtk.FILL|gtk.EXPAND, 0, 0, 0)
        self.run_button.connect("clicked", self.run_tls_analysis, None)
        
        ## make the print box
        sw = gtk.ScrolledWindow()
        self.vbox.pack_end(sw, gtk.TRUE, gtk.TRUE, 0)
        sw.set_border_width(5)
        sw.set_shadow_type(gtk.SHADOW_ETCHED_IN)
        sw.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
     
        self.treeview = DictListTreeView(
            column_list=["Residue Range", "Atoms", "R", "<DP2>", "sig<DP2>"])
        sw.add(self.treeview)
       
        self.show_all()

    def response_cb(self, dialog, response_code):
        """Responses to dialog events.
        """
        if response_code==gtk.RESPONSE_CLOSE or \
           response_code==gtk.RESPONSE_CANCEL:
            pass
        elif response_code==gtk.RESPONSE_OK:
            pass

    def destroy_cb(self, *args):
        """Destroy the TLS dialog and everything it has built
        in the GLObject viewer.
        """
        pass
    
    def run_tls_analysis(self, *args):
        """
        """
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

        try:
            max_DP2 = float(self.dp2_entry.get_text())
        except ValueError:
            max_DP2 = 1.0
        
        tls_analysis = TLSStructureAnalysis(self.sc.struct)
##         stats_list   = tls_analysis.fit_TLS_segments(
##             residue_width          = residue_width,
##             use_side_chains        = use_side_chains,
##             include_frac_occupancy = include_frac_occupancy,
##             include_single_bond    = include_single_bond)

        stats_list = tls_analysis.fit_common_TLS(
            residue_width          = residue_width )

        self.tls_stats_list = []
        for stats in stats_list:

            if stats["mean_DP2"]>max_DP2:
                continue

            stats["Residue Range"] = stats["name"]
            stats["R"]             = "%.3f" % (stats["R"])
            stats["<DP2>"]         = "%.4f" % (stats["mean_DP2"])
            stats["Atoms"]         = str(stats["num_atoms"])

            self.tls_stats_list.append(stats)

        self.treeview.set_dict_list(self.tls_stats_list)


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
        self.window.set_title("mmLib Viewer")
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

            ('/_Tools', None, None, 0, '<Branch>'),
            ('/Tools/TLS Analysis...', None, self.tools_tls_analysis, 0, None),

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

        tab["name"]      = "Page %d" % (tab["num"])
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
                build_properties = ("sequence","bonds"))

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

try:
    first_path = sys.argv[1]
except IndexError:
    first_path = None
    
VIEWER_APP = ViewerApp(first_path)
VIEWER_APP.main()

### </MAIN>
