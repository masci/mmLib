## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""OpenGL rendering classes.
"""
from __future__  import generators

from OpenGL.GL   import *
from OpenGL.GLU  import *
from OpenGL.GLUT import *
from mmTypes     import *
from Structure   import *
from Extensions.TLS import *


## <hack>

GL_MATERIALS = [
    { "name":         "emerald",
      "GL_AMBIENT":   (0.0215, 0.1745, 0.0215),
      "GL_DIFFUSE":   (0.07568, 0.61424, 0.07568),
      "GL_SPECULAR":  (0.633, 0.727811, 0.633),
      "GL_SHININESS": 0.6 },

    { "name":         "jade",
      "GL_AMBIENT":   (0.135, 0.2225, 0.1575),
      "GL_DIFFUSE":   (0.54, 0.89, 0.63),
      "GL_SPECULAR":  (0.316228, 0.316228, 0.316228),
      "GL_SHININESS": 0.1 },

    { "name":         "obsidian",
      "GL_AMBIENT":   (0.05375, 0.05, 0.06625),
      "GL_DIFFUSE":   (0.18275, 0.17, 0.22525),
      "GL_SPECULAR":  (0.332741, 0.328634, 0.346435),
      "GL_SHININESS":  0.3 },

    { "name":         "pearl",
      "GL_AMBIENT":   (0.25, 0.20725, 0.20725),
      "GL_DIFFUSE":   (1.0, 0.829, 0.829),
      "GL_SPECULAR":  (0.296648, 0.296648, 0.296648),
      "GL_SHININESS": 0.088 },

    { "name":         "ruby",
      "GL_AMBIENT":   (0.1745, 0.01175, 0.01175),
      "GL_DIFFUSE":   (0.61424, 0.04136, 0.04136),
      "GL_SPECULAR":  (0.727811, 0.626959, 0.626959),
      "GL_SHININESS": 0.6 },

    { "name":         "turquoise",
      "GL_AMBIENT":   (0.1, 0.18725, 0.1745),
      "GL_DIFFUSE":   (0.396, 0.74151, 0.69102),
      "GL_SPECULAR":  (0.297254, 0.30829, 0.306678),
      "GL_SHININESS": 0.1 },

    { "name":         "brass",
      "GL_AMBIENT":   (0.329412, 0.223529, 0.027451),
      "GL_DIFFUSE":   (0.780392, 0.568627, 0.113725),
      "GL_SPECULAR":  (0.992157, 0.941176, 0.807843),
      "GL_SHININESS": 0.21794872 },

    { "name":         "bronze",
      "GL_AMBIENT":   (0.2125, 0.1275, 0.054),
      "GL_DIFFUSE":   (0.714, 0.4284, 0.18144),
      "GL_SPECULAR":  (0.393548, 0.271906, 0.166721),
      "GL_SHININESS": 0.2 },

    { "name":         "chrome",
      "GL_AMBIENT":   (0.25, 0.25, 0.25,),
      "GL_DIFFUSE":   (0.4, 0.4, 0.4),
      "GL_SPECULAR":  (0.774597, 0.774597, 0.774597),
      "GL_SHININESS": 0.6 },

    { "name":         "copper",
      "GL_AMBIENT":   (0.19125, 0.0735, 0.0225),
      "GL_DIFFUSE":   (0.7038, 0.27048, 0.0828),
      "GL_SPECULAR":  (0.256777, 0.137622, 0.086014),
      "GL_SHININESS": 0.1 },

    { "name":         "gold",
      "GL_AMBIENT":   (0.24725, 0.1995, 0.0745),
      "GL_DIFFUSE":   (0.75164, 0.60648, 0.22648),
      "GL_SPECULAR":  (0.628281, 0.555802, 0.366065),
      "GL_SHININESS": 0.4 },

    { "name":         "silver",
      "GL_AMBIENT":   (0.19225, 0.19225, 0.19225),
      "GL_DIFFUSE":   (0.50754, 0.50754, 0.50754),
      "GL_SPECULAR":  (0.508273, 0.508273, 0.508273),
      "GL_SHININESS": 0.4 },

    { "name":         "black",
      "GL_AMBIENT":   (0.02, 0.02, 0.02),
      "GL_DIFFUSE":   (0.01, 0.01, 0.01),
      "GL_SPECULAR":  (0.4, 0.4, 0.4),
      "GL_SHININESS": 0.078125 },

    { "name":         "cyan",
      "GL_AMBIENT":   (0.0, 0.1, 0.06),
      "GL_DIFFUSE":   (0.0, 0.50980392, 0.50980392),
      "GL_SPECULAR":  (0.50196078, 0.50196078, 0.50196078),
      "GL_SHININESS": 0.25 },

    { "name":         "green",
      "GL_AMBIENT":   (0.0, 0.0, 0.0),
      "GL_DIFFUSE":   (0.1, 0.35, 0.1),
      "GL_SPECULAR":  (0.45, 0.55, 0.45),
      "GL_SHININESS": 0.25 },

    { "name":         "red",
      "GL_AMBIENT":   (0.0, 0.0, 0.0),
      "GL_DIFFUSE":   (0.5, 0.0, 0.0),
      "GL_SPECULAR":  (0.7, 0.6, 0.6,),
      "GL_SHININESS": 0.25 },

    { "name":         "white",
      "GL_AMBIENT":   (0.0, 0.0, 0.0),
      "GL_DIFFUSE":   (0.55, 0.55, 0.55),
      "GL_SPECULAR":  (0.70, 0.70, 0.70),
      "GL_SHININESS": 0.25 },

    { "name":         "yellow",
      "GL_AMBIENT":   (0.0, 0.0, 0.0),
      "GL_DIFFUSE":   (0.5, 0.5, 0.0),
      "GL_SPECULAR":  (0.60, 0.60, 0.50),
      "GL_SHININESS": 0.25 },
    ]


GL_MATERIALS_DICT = {}
for gl_material in GL_MATERIALS:
    GL_MATERIALS_DICT[gl_material["name"]] = gl_material
    
## </hack>


class GLPropertyDict(dict):
    """Property cache/routing dictionary
    """
    def __init__(self, gl_object):
        self.gl_object = gl_object

    def update(self, **args):
        self.gl_object.glo_update_properties(**args)


class GLObject(object):
    """Base class for all OpenGL rendering objects.  It combines a
    composite-style tree structure with a system for setting properties.
    The properties are used for the specific OpenGL drawing objects
    to control color, position, line width, etc...  Implementing properties
    requres the GLProperties object which is the access object for the
    properties.
    """
    def __init__(self, **args):
        object.__init__(self)

        self.__globject_parent               = None
        self.__globject_children             = []

        self.properties                      = GLPropertyDict(self)
        self.__globject_properties_id        = None
        self.__globject_properties           = []
        self.__globject_properties_callbacks = []

        self.glo_install_properties()

    def glo_add_child(self, child):
        assert isinstance(child, GLObject)
        assert child.__globject_parent==None
        child.__globject_parent = self
        self.__globject_children.append(child)
        
    def glo_remove_child(self, child):
        assert isinstance(child, GLObject)
        assert child.__globject_parent==self
        child.__globject_parent = None
        self.__globject_children.remove(child)

    def glo_iter_children(self):
        """Iterate immediate children.
        """
        return iter(self.__globject_children)

    def glo_iter_preorder_traversal(self):
        """Preorder Traversal for GLObject composite.
        """
        for child1 in self.__globject_children:
            yield child1
            for child2 in child1.glo_iter_preorder_traversal():
                yield child2

    def glo_get_depth(self):
        """Returns the depth, the root composite is depth 0.
        """
        depth = 0
        ancestor = self.__globject_parent
        while ancestor:
            depth += 1
            ancestor = ancestor.parent
        return depth

    def glo_get_degree(self):
        """Returns the number of children (degree).
        """
        return len(self.__globject_children)

    def glo_count_descendants(self):
        """Counts all decendant GLObjects.
        """
        n = self.getDegree()
        for child in self.__globject_children:
            n += child.glo_count_descendants()
        return n

    def glo_get_root(self):
        """Returns the root GLObject.
        """
        gl_object = self
        while gl_object.__globject_parent:
            gl_object = gl_object.__globject_parent
        return gl_object

    def glo_get_parent(self):
        """Returns the parent GLObject.
        """
        return self.__globject_parent

    def glo_get_path(self):
        """Returns the tree-path to the composite as a list of its
        parent composites.
        """
        path_list = [self]
        parent = self.__globject_parent
        while parent:
            path_list.insert(0, parent)
            parent = parent.__parent
        return path_list

    def glo_get_index_path(self):
        """Returns the tree-path to the GLObject as a list of its
        integer indexes.
        """
        ipath_list = []
        child = self
        parent = child.__parent
        while parent:
            ipath_list.insert(0, parent.__child_list.index(child))
            child = parent
            parent = parent.__parent
        return ipath_list

    def glo_get_parent_list(self):
        """Returns a list of the parent GLObjects back to the root.
        """
        list = []
        composite = self
        while composite.__parent:
            composite = composite.__parent
            list.append(composite)
        return list

    def glo_get_lowest_common_ancestor(self, gl_object):
        """Returns the lowest common ancesotry of self and argument
        composite.
        """
        assert isinstance(gl_object, GLObject)

        pl1 = self.getParentList()
        pl2 = gl_object.getParentList()
        pl1.reverse()
        pl2.reverse()

        ancestor = None
        for i in range(min(len(pl1), len(pl2))):
            if pl1[i] == pl2[i]:
                ancestor = pl1[i]
            else:
                break

        return ancestor

    def glo_is_descendant_of(self, gl_object):
        """Returns true if self composite is a decent of argument GLObject.
        """
        assert isinstance(gl_object, GLObject)
        
        ancestor = self.__globject_parent
        while ancestor:
            if ancestor==gl_object:
                return True
            ancestor = ancestor.__globject_parent
        return False

    def glo_add_child(self, child):
        assert isinstance(child, GLObject)
        assert child!=self
        assert self.glo_is_descendant_of(child)==False
        assert child.__globject_parent==None
        child.__globject_parent = self
        self.__globject_children.append(child)

    def glo_prepend_child(self, child):
        """Adds a child GLObject to the beginning of the GLObject's
        child list.
        """
        assert isinstance(child, GLObject)
        assert child!=self
        assert self.glo_is_descendant_of(child)==False
        assert child.__globject_parent==None
        child.__globject_parent = self
        self.__globject_children.insert(0, child)
            
    def glo_append_child(self, child):
        """Adds a child GLObject to the end of the GLObject's child list.
        """
        assert isinstance(child, GLObject)
        assert child!=self
        assert self.glo_is_descendant_of(child)==False
        assert child.__globject_parent==None
        child.__globject_parent = self
        self.__globject_children.append(child)
        
    def glo_remove_child(self, child):
        assert isinstance(child, GLObject)
        assert child.__globject_parent==self
        child.__globject_parent = None
        self.__globject_children.remove(child)

    def glo_remove(self):
        """Removes GLObject.
        """
        parent = self.__globject_parent
        if parent==None:
            return
        parent.glo_remove_child(self)

    def glo_set_properties_id(self, gl_object_id):
        """Set the property name for this GLObject.
        """
        self.__globject_properties_id = gl_object_id

    def glo_get_properties_id(self):
        return self.__globject_properties_id

    def glo_install_properties(self):
        """Called by GLObject.__init__ to install properties.
        """
        pass

    def glo_add_property(self, prop_desc):
        """Adds a new property to the GLObject.
        """
        assert prop_desc["name"] not in self.properties
        self.__globject_properties.append(prop_desc)
        self.properties[prop_desc["name"]] = prop_desc["default"]

    def glo_iter_property_desc(self):
        """Iterates over all property descriptions.
        """
        return iter(self.__globject_properties)
        
    def glo_get_property_desc(self, name):
        """Return the property description dictionary for the given
        property name.
        """
        for prop_desc in self.__globject_properties:
            if prop_desc["name"]==name:
                return prop_desc
        return None

    def glo_link_child_property(self, name, child_gl_object_id, child_name):
        """Link the value of the GLObject's property to the value of
        a child property.
        """
        prop_desc = self.glo_get_property_desc(name)

        link_dict = {"gl_object": child_gl_object_id,
                     "name":      child_name}

        try:
            prop_desc["link"].append(link_dict)
        except KeyError:
            prop_desc["link"] = [link_dict]

    def glo_get_child(self, gl_object_id):
        """Returns the child GLObject matching the given gl_object_id.
        """
        for gl_object in self.glo_iter_children():
            if gl_object.__globject_properties_id==gl_object_id:
                return gl_object
        return None

    def glo_init_properties(self, **args):
        """This is a special form of update which propagates all linked
        values, not just the changed ones.
        """
        for prop_desc in self.__globject_properties:
            name = prop_desc["name"]

            try:
                self.properties[name] = args[name]
            except KeyError:
                pass

            try:
                linked_props = prop_desc["link"]
            except KeyError:
                pass
            else:
                for linked_prop in linked_props:
                    gl_object_id = linked_prop["gl_object"]
                    child_name   = linked_prop["name"]
                    gl_object    = self.glo_get_child(gl_object_id)
                    child_args   = {child_name: self.properties[name]}
                    gl_object.glo_init_properties(**child_args)

            if prop_desc.get("cascade", False)==True:
                child_args = {name: self.properties[name]}
                for gl_object in self.glo_iter_preorder_traversal():
                    gl_object.glo_init_properties(**child_args)

    def glo_update_properties(self, **args):
        """Update property values and trigger update callbacks.
        """        
        updates = {}
        actions = []

        ## update properties
        for prop_desc in self.__globject_properties:
            name = prop_desc["name"]
            old_value = self.properties[name]

            try:
                self.properties[name] = args[name]
            except KeyError:
                continue

            if self.properties[name]!=old_value:
                updates[name] = self.properties[name]
                if prop_desc["action"] not in actions:
                    actions.append(prop_desc["action"])

            ## linked properties
            try:
                linked_props = prop_desc["link"]
            except KeyError:
                pass
            else:
                for linked_prop in linked_props:
                    gl_object_id = linked_prop["gl_object"]
                    child_name   = linked_prop["name"]
                    gl_object    = self.glo_get_child(gl_object_id)
                    child_args   = {child_name: self.properties[name]}
                    gl_object.glo_update_properties(**child_args)

            ## cascadeing properties
            if prop_desc.get("cascade", False)==True:
                child_args = {name: self.properties[name]}
                for gl_object in self.glo_iter_preorder_traversal():
                    gl_object.glo_init_properties(**child_args)

        if len(updates):
            for func in self.__globject_properties_callbacks:
                func(updates, actions)

    def glo_add_update_callback(self, func):
        """Adds a function which is called whenever property values change.
        The function is called with two arguments: a updates dictionary
        containing all updated properties and the values they were changed
        to, and a actions list which contains a unique list of action
        key words forme self.prop_list = []
        self.callback_list = [] from the action fields of the updated
        properties.
        """
        self.__globject_properties_callbacks.append(func)

    def glo_remove_update_callback(self, func):
        """Removes the update callback.
        """
        self.__globject_properties_callbacks.remove(func)

    def glo_redraw(self):
        """Triggers a redraw of the GLViewer
        """
        gl_viewer = self.glo_get_root()
        if isinstance(gl_viewer, GLViewer):
            gl_viewer.gl_redraw()

    def glo_get_glstructure(self):
        """Returns the parent GLStructure object, or None if the GLObject
        is not a child of a GLStructure.
        """
        gl_object = self
        while gl_object!=None and isinstance(gl_object, GLStructure)==False:
            gl_object = gl_object.__globject_parent
        return gl_object


class GLDrawList(GLObject):
    """Fundamental OpenGL rigid entity.
    """
    def __init__(self, **args):
        GLObject.__init__(self, **args)
        self.gl_name = None
        self.glo_add_update_callback(self.update_cb)
        self.glo_init_properties(**args)

    def update_cb(self, updates, actions):
        if "recompile" in actions:
            self.gl_delete_list()
            self.glo_redraw()
        if "redraw" in actions:
            self.glo_redraw()

    def glo_install_properties(self):
        self.glo_add_property(
            { "name" :      "visible",
              "desc":       "Visible",
              "type":       "boolean",
              "default":    True,
              "action":     "redraw" })
        self.glo_add_property(
            { "name" :      "origin",
              "desc":       "Origin",
              "type":       "array(3)",
              "hidden":     True,
              "default":    zeros(3, Float),
              "action":     "redraw" })
        self.glo_add_property(
            { "name" :      "axes",
              "desc":       "Rotation Axes",
              "type":       "array(3,3)",
              "hidden":     True,
              "default":    identity(3),
              "action":     "redraw" })
        self.glo_add_property(
            { "name" :      "rot_x",
              "desc":       "Degrees Rotation About X Axis",
              "type":       "float",
              "hidden":     True,
              "default":    0.0,
              "action":     "redraw" })
        self.glo_add_property(
            { "name" :      "rot_y",
              "desc":       "Degrees Rotation About Y Axis",
              "type":       "float",
              "hidden":     True,
              "default":    0.0,
              "action":     "redraw" })
        self.glo_add_property(
            { "name" :      "rot_z",
              "desc":       "Degrees Rotation About Z Axis",
              "type":       "float",
              "hidden":     True,
              "default":    0.0,
              "action":     "redraw" })
        
    def gl_push_matrix(self):
        """Rotate and translate to the correct position for drawing.
        """
        glPushMatrix()

        if not allclose(self.properties["origin"], zeros(3, Float)):
            glTranslatef(*self.properties["origin"])

        axes = self.properties["axes"]
        if not allclose(axes, identity(3, Float)):
            glRotatef(self.properties["rot_x"], *axes[0])
            glRotatef(self.properties["rot_y"], *axes[1])
            glRotatef(self.properties["rot_z"], *axes[2])

    def gl_pop_matrix(self):
        """Pop the roatated/translated position.
        """
        glPopMatrix()

    def gl_render(self):
        """Compile or force a recompile of this object's gl_draw list, and
        render the scene.  Rendering the scene can be bypassed if
        this method is called with render = False.
        """
        if self.properties["visible"]==False:
            return

        if self.gl_name == None:
            self.gl_compile_list()

        self.gl_push_matrix()
        self.gl_call_list()
        self.gl_pop_matrix()

    def gl_call_list(self):
        """Calls glCallList on the GLDrawList's compiled list.
        """
        glCallList(self.gl_name)

    def gl_compile_list(self):
        """Compile a OpenGL draw list for this object.
        """
        if self.gl_name != None:
            self.gl_delete_list()
    
        self.gl_name = glGenLists(1)
        glNewList(self.gl_name, GL_COMPILE)
        self.gl_draw()
        glEndList()

    def gl_delete_list(self):
        """Delete the currently compiled glList for this object.  This forces
        a recompile next time gl_call_list_render() is invoked.  It also
        should be called before a GLDrawList is dereferenced so that unused
        compiled glLists are not left in the OpenGL libraries.
        """
        if self.gl_name != None:
            glDeleteLists(self.gl_name, 1)
            self.gl_name = None

    def gl_draw(self):
        """Implement in subclass to draw somthing.
        """
        pass


class GLDrawListContainer(GLDrawList):
    """Base class for draw lists which contain other draw lists.
    """
    def gl_render(self):
        if self.properties["visible"]==False:
            return
        
        self.gl_push_matrix()

        for draw_list in self.glo_iter_children():
            draw_list.gl_render()

        if self.gl_name==None:
            self.gl_compile_list()

        self.gl_call_list()

        self.gl_pop_matrix()

    def gl_delete_list(self):
        GLDrawList.gl_delete_list(self)
        for draw_list in self.glo_iter_children():
            draw_list.gl_delete_list()


class GLAxes(GLDrawList):
    """Draw orthogonal axes in red = x, green = y, blue = z.
    """
    def __init__(self, **args):
        GLDrawList.__init__(self, **args)
        self.glo_init_properties(**args)

    def glo_install_properties(self):
        GLDrawList.glo_install_properties(self)

        self.glo_add_property(
            { "name":       "line_length",
              "type":       "float",
              "default":    200.0,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "line_width",
              "type":       "float",
              "default":    10.0,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "color_x",
              "type":       "color",
              "default":    (1.0, 0.0, 0.0),
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "color_y",
              "type":       "color",
              "default":    (0.0, 1.0, 0.0),
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "color_z",
              "type":       "color",
              "default":    (0.0, 0.0, 1.0),
              "action":     "recompile" })

    def gl_draw(self):
        glDisable(GL_LIGHTING)

        glLineWidth(self.properties["line_width"])
        
        def axis_line(u):
            glBegin(GL_LINES)
            glVertex3f(0.0, 0.0, 0.0)
            glVertex3f(*u)
            glEnd()
                    
        glColor3f(*self.properties["color_x"])
        axis_line(array([self.properties["line_length"], 0.0, 0.0]))

        glColor3f(*self.properties["color_y"])
        axis_line(array([0.0, self.properties["line_length"], 0.0]))

        glColor3f(*self.properties["color_z"])
        axis_line(array([0.0, 0.0, self.properties["line_length"]]))

        glEnable(GL_LIGHTING)


class GLUnitCell(GLDrawList):
    """Draw unit cell.
    """
    def __init__(self, **args):
        GLDrawList.__init__(self, **args)
        self.unit_cell = args["unit_cell"]
        self.glo_init_properties(**args)

    def glo_install_properties(self):
        GLDrawList.glo_install_properties(self)
        
        self.glo_add_property(
            { "name":       "line_width",
              "type":       "float",
              "default":    2.0,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "color",
              "type":       "color",
              "default":    (1.0, 1.0, 1.0),
              "action":     "recompile" })

    def draw_cell(self, x1, y1, z1, x2, y2, z2):
        """Draw the unit cell lines in a rectangle starting at fractional
        integer coordinates x1, y1, z1, ending at x2, y2, z2.  The first set of
        coordinates must be <= the second set.
        """
        assert x1 <= x2 and y1 <= y2 and z1 <= z2

        a = self.unit_cell.calc_frac_to_orth(array([1.0, 0.0, 0.0]))
        b = self.unit_cell.calc_frac_to_orth(array([0.0, 1.0, 0.0]))
        c = self.unit_cell.calc_frac_to_orth(array([0.0, 0.0, 1.0]))

        glDisable(GL_LIGHTING)
        glLineWidth(self.properties["line_width"])
        glColor3f(*self.properties["color"])

        for k in range(z1, z2+2):
            for j in range(y1, y2+2):
                glBegin(GL_LINES)
                glVertex3f(*     x1*a + j*b + k*c)
                glVertex3f(* (x2+1)*a + j*b + k*c)
                glEnd()

        for k in range(z1, z2+2):
            for i in range(x1, x2+2):
                glBegin(GL_LINES)
                glVertex3f(* i*a +     y1*b + k*c)
                glVertex3f(* i*a + (y2+1)*b + k*c)
                glEnd()

        for j in range(y1, y2+2):
            for i in range(x1, x2+2):
                glBegin(GL_LINES)
                glVertex3f(* i*a + j*b +     z1*c)
                glVertex3f(* i*a + j*b + (z2+1)*c)
                glEnd()

        glEnable(GL_LIGHTING)

    def draw_sphere(self):
        radius = max(self.unit_cell.a,
                     max(self.unit_cell.b, self.unit_cell.c))

        glEnable(GL_LIGHTING)

        gl_material = GL_MATERIALS_DICT["pearl"]
        alpha = 0.25
        
        r, g, b = gl_material["GL_AMBIENT"]
        glMaterialfv(GL_FRONT, GL_AMBIENT, (r, g, b, alpha))
        r, g, b = gl_material["GL_DIFFUSE"]
        glMaterialfv(GL_FRONT, GL_DIFFUSE, (r, g, b, alpha))
        r, g, b = gl_material["GL_SPECULAR"]
        glMaterialfv(GL_FRONT, GL_SPECULAR, (r, g, b, alpha))
        glMaterialf(GL_FRONT, GL_SHININESS, gl_material["GL_SHININESS"]*128.0)

        glutSolidSphere(radius, 128, 128)

    def gl_draw(self):
        self.draw_cell(-1, -1, -1, 0, 0, 0)
        #self.draw_sphere()


class GLAtomList(GLDrawList):
    """OpenGL renderer for a list of atoms.  Optional arguments iare:
    color, U, U_color.
    """
    def __init__(self, **args):
        GLDrawList.__init__(self, **args)

        self.atom_list      = AtomList()
        self.el_color_cache = {}
        self.symop          = None

        self.glo_init_properties(**args)

    def glo_install_properties(self):
        GLDrawList.glo_install_properties(self)
        
        self.glo_add_property(
            { "name":       "sphere_quality",
              "desc":       "CPK Sphere Quality",
              "type":       "integer",
              "default":    12,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "line_width",
              "desc":       "Atom Line Drawing Width",
              "type":       "float",
              "default":    3.0,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "atom_origin",
              "type":       "array(3)",
              "hidden":     True,
              "default":    None,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "lines",
              "desc":       "Draw Atom Bond Lines",
              "type":       "boolean",
              "default":    True,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "trace",
              "desc":       "Draw Backbone Trace",
              "type":       "boolean",
              "default":    False,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "cpk",
              "desc":       "Draw CPK Spheres",
              "type":       "boolean",
              "default":    False,
              "action":    "recompile" })
        self.glo_add_property(
            { "name":       "cpk_opacity",
              "desc":       "CPK Sphere Opacity",
              "type":       "float",
              "default":    0.10,
              "action":    "recompile" })
        self.glo_add_property(
            { "name":      "U",
              "desc":      "Draw ADP Axes",
              "type":      "boolean",
              "default":   False,
              "action":    "recompile" })
        self.glo_add_property(
            { "name":      "atm_U_attr",
              "type":      "string",
              "hidden":    True,
              "default":   "U",
              "action":    "recompile" })
        self.glo_add_property(
            { "name":      "U_color",
              "desc":      "Thermal Axes Color",
              "type":      "color",
              "default":   (1.0, 1.0, 1.0),
              "action":    "recompile" })
        self.glo_add_property(
            { "name":      "color",
              "desc":      "Solid Color",
              "type":      "color",
              "default":   None,
              "action":    "recompile" })
        self.glo_add_property(
            { "name":      "color_func",
              "type":      "function",
              "hidden":    True,
              "default":   None,
              "action":    "recompile" })
        self.glo_add_property(
            { "name":      "symmetry",
              "type":      "boolean",
              "default":   True,
              "action":    "recompile" })

    def gl_call_list(self):
        """Specialized draw list invokation to recycle the draw list for
        symmetry related copies.  Cartesian versions of the symmetry rotation
        and translation operators are generated by GLStructure/UnitCell
        classes.
        """
        gl_struct = self.glo_get_glstructure()
        if gl_struct==None or self.properties["symmetry"]==False:
            GLDrawList.gl_call_list(self)
            return

        for symop in gl_struct.iter_orth_symops():

##             import random
##             color = (random.random(), random.random(), random.random())
##             if allclose(symop.R, identity(3)):
##                 continue

##             save = {}
##             for atm in self.atom_list:
##                 save[atm] = atm.position
##                 atm.position = symop(atm.position)

##             self.draw_trace(color)

##             for atm in self.atom_list:
##                 atm.position = save[atm]
            
            glPushMatrix()

            if self.properties["atom_origin"]!=None:
                glTranslatef(* -self.properties["atom_origin"])

            glMultMatrixf(
                (symop.R[0,0], symop.R[1,0], symop.R[2,0], 0.0,
                 symop.R[0,1], symop.R[1,1], symop.R[2,1], 0.0,
                 symop.R[0,2], symop.R[1,2], symop.R[2,2], 0.0,
                 symop.t[0],   symop.t[1],   symop.t[2],   1.0) )

            if self.properties["atom_origin"]!=None:
                glTranslatef(*self.properties["atom_origin"])

            #self.draw_trace(color, 2.0)
            GLDrawList.gl_call_list(self)

            glPopMatrix()

    def gl_draw_old(self):
        """Perform the OpenGL drawing operations to render this atom list
        with the current settings.
        """
        if self.properties["symmetry"]==True:
            gl_struct = self.glo_get_glstructure()

            if gl_struct==None:
                self.draw_all()
            else:
                for symop in gl_struct.iter_orth_symops():
                    self.symop = symop
                    self.draw_all()

        else:
            self.draw_all()

    def gl_draw(self):
        """Draw all selected representations.
        """
        ## per-atom draw functions
        draw_funcs = []
        
        if self.properties["lines"]==True:
            draw_funcs.append(self.draw_lines)
        if self.properties["cpk"]==True:
            draw_funcs.append(self.draw_cpk)
        if self.properties["U"]==True:
            draw_funcs.append(self.draw_U_axes)

        for atm in self.atom_list:
            for draw_func in draw_funcs:
                draw_func(atm)

        ## drawing functions requiring neighboring atoms
        if self.properties["trace"]==True:
            self.draw_trace()

    def calc_position(self, pos):
        """
        """
        if self.properties["atom_origin"]!=None:
            pos = pos - self.properties["atom_origin"]
        #if self.symop!=None:
        #    pos = self.symop(pos)
        return pos
                
    def get_color(self, atm):
        """Sets the open-gl color for the atom.
        """
        ## case 1: set color by color callback function
        if self.properties["color_func"]!=None:
            color_func = self.properties["color_func"]
            color = color_func(atm)
            return color

        ## case 2: set color by solid color
        elif self.properties["color"]!=None:
            color = self.properties["color"]
            return color

        ## case 3: set color by atom type
        else:
            try:
                return self.el_color_cache[atm.element]
            except KeyError:
                elem = atm.get_structure().library.get_element(atm.element)
                if elem!=None:
                    color = elem.color
                else:
                    color = (1.0, 1.0, 1.0)

                self.el_color_cache[atm.element] = color
                return color

    def draw_cpk(self, atm):
        """Draw a atom as a CPK sphere.
        """
        elem = atm.get_structure().library.get_element(atm.element)
        if elem:
            radius = elem.van_der_waals_radius
        else:
            radius = 2.0

        r, g, b = self.get_color(atm)
        a       = self.properties["cpk_opacity"]

        ambient  = array([r, g, b, a])
        diffuse  = ambient
        specular = ambient * 1.25
        
        sphere_quality = self.properties["sphere_quality"]

        glEnable(GL_LIGHTING)
        
	glMaterial(GL_FRONT, GL_AMBIENT,   ambient)
	glMaterial(GL_FRONT, GL_DIFFUSE,   diffuse)
	glMaterial(GL_FRONT, GL_SPECULAR,  specular)
	glMaterial(GL_FRONT, GL_SHININESS, 50.0)
        
        glPushMatrix()
        glTranslatef(*self.calc_position(atm.position))
        glutSolidSphere(radius, sphere_quality, sphere_quality)
        glPopMatrix()

    def draw_lines(self, atm):
        """Draw a atom using bond lines only.
        """
        pos = self.calc_position(atm.position)

        glDisable(GL_LIGHTING)
        glColor3f(*self.get_color(atm))
        glLineWidth(self.properties["line_width"])
        
        ## if there are bonds, then draw the lines 1/2 way to the
        ## bonded atoms
        if len(atm.bond_list) > 0:
            glBegin(GL_LINES)
            
            for bond in atm.iter_bonds():
                atm2 = bond.get_partner(atm)
                pos2 = self.calc_position(atm2.position)

                start = pos
                end   = start + ((pos2 - pos) / 2)

                glVertex3f(*start)
                glVertex3f(*end)

            glEnd()

        ## if there are no bonds, draw a small cross-point 
        else:
            self.draw_cross(atm)
            self.draw_cpk(atm)

    def draw_cross(self, atm):
        """Draws atom with a cross of lines.
        """
        position = self.calc_position(atm.position)

        glDisable(GL_LIGHTING)
        glColor3f(*self.get_color(atm))        
        glLineWidth(self.properties["line_width"])

        vx = array([0.25, 0.0,  0.0])
        vy = array([0.0,  0.25, 0.0])
        vz = array([0.0,  0.0,  0.25])

        glBegin(GL_LINES)
        
        start = position - vx
        end   = position + vx
        glVertex3f(*start)
        glVertex3f(*end)
        
        start = position - vy
        end   = position + vy
        glVertex3f(*start)
        glVertex3f(*end)
        
        start = position - vz
        end   = position + vz
        glVertex3f(*start)
        glVertex3f(*end)

        glEnd()

    def draw_trace(self, color=(1.0, 1.0, 1.0), line_width=5.0):
        """Draws trace over all CA atoms.
        """
        glDisable(GL_LIGHTING)
        glColor3f(*color)
        glLineWidth(self.properties["line_width"])
        glLineWidth(line_width)
        
        glBegin(GL_LINE_STRIP)

        for atm in self.atom_list:
            if atm.name in ["N", "CA", "C"]:
                glVertex3f(*self.calc_position(atm.position))

        glEnd()

    def draw_U_axes(self, atm):
        """Draw the anisotropic axies of the atom with R.M.S.*sqrt(2)
        magnitude.
        """
        U = getattr(atm, self.properties["atm_U_attr"], None)
        if U==None:
            return

        eigen_values, eigen_vectors = eigenvectors(U)

        v0_peak = 1.414 * math.sqrt(abs(eigen_values[0]))
        v1_peak = 1.414 * math.sqrt(abs(eigen_values[1]))
        v2_peak = 1.414 * math.sqrt(abs(eigen_values[2]))
        
        v0 = eigen_vectors[0] * v0_peak
        v1 = eigen_vectors[1] * v1_peak
        v2 = eigen_vectors[2] * v2_peak

        position = self.calc_position(atm.position)


        glDisable(GL_LIGHTING)
        glColor3f(*self.properties["U_color"])
        glLineWidth(1.0)
        
        glBegin(GL_LINES)

        glVertex3f(*position - v0)
        glVertex3f(*position + v0)
        glVertex3f(*position - v1)
        glVertex3f(*position + v1)
        glVertex3f(*position - v2)
        glVertex3f(*position + v2)

        glEnd()



class GLSymmetry(GLDrawListContainer):
    """
    """
    def gl_call_list(self):
        gl_struct = self.glo_get_glstructure()
        if gl_struct==None or self.properties["symmetry"]==False:
            glPushMatrix()
    

class GLTLSAtomList(GLAtomList):
    """
    """
    def glo_install_properties(self):
        GLAtomList.glo_install_properties(self)

        self.glo_add_property(
            { "name":        "L_eigen_vec",
              "desc":        "L Eigen Vectors", 
              "type":        "array(3,3)",
              "default":     identity(3),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1",
              "desc":        "L1 Rotation", 
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L2",
              "desc":        "L2 Rotation", 
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L3",
              "desc":        "L3 Rotation", 
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })

    def update_cb(self, updates, actions):
        GLAtomList.update_cb(self, updates, actions)
        if "trace" in updates:
            self.properties.update(lines=False)
 
    def gl_call_list(self):
        gl_struct = self.glo_get_glstructure()
        if gl_struct==None or self.properties["symmetry"]==False:
            glPushMatrix()
            
            axes = self.properties["L_eigen_vec"]
            glRotatef(self.properties["L1"], *axes[0])
            glRotatef(self.properties["L2"], *axes[1])
            glRotatef(self.properties["L3"], *axes[2])
            
            GLDrawList.gl_call_list(self)
            glPopMatrix()
            return

        for symop in gl_struct.iter_orth_symops():
            glPushMatrix()

            if self.properties["atom_origin"]!=None:
                glTranslatef(* -self.properties["atom_origin"])

            glMultMatrixf(
                (symop.R[0,0], symop.R[1,0], symop.R[2,0], 0.0,
                 symop.R[0,1], symop.R[1,1], symop.R[2,1], 0.0,
                 symop.R[0,2], symop.R[1,2], symop.R[2,2], 0.0,
                 symop.t[0],   symop.t[1],   symop.t[2],   1.0) )

            if self.properties["atom_origin"]!=None:
                glTranslatef(*self.properties["atom_origin"])

            axes = self.properties["L_eigen_vec"]
            glRotatef(self.properties["L1"], *axes[0])
            glRotatef(self.properties["L2"], *axes[1])
            glRotatef(self.properties["L3"], *axes[2])
            
            GLDrawList.gl_call_list(self)

            glPopMatrix()
                
        
class GLTLSGroup(GLDrawListContainer):
    """Draws TLS group
    """
    def __init__(self, **args):
        GLDrawListContainer.__init__(self)

        ## TLS calculations

        ## step 1: copy the TLS group
        orig_tls_group = args["tls_group"]

        self.tls_group = TLSGroup()

        for atm in orig_tls_group:
            self.tls_group.append(atm)
        
        self.tls_group.origin = orig_tls_group.origin.copy()
        self.tls_group.T = orig_tls_group.T.copy()
        self.tls_group.L = orig_tls_group.L.copy()
        self.tls_group.S = orig_tls_group.S.copy()

        ## step 2: transform T,L,S to center of reaction
        self.calcs = self.tls_group.calc_COR()
        self.tls_group.origin = self.calcs["COR"].copy()
        self.tls_group.T = self.calcs["T'"].copy()
        self.tls_group.L = self.calcs["L'"].copy()
        self.tls_group.S = self.calcs["S'"].copy()

        (eigen_values, eigen_vectors) = eigenvectors(self.tls_group.L)
        self.evalL = eigen_values

        ## add child GLAtomList
        self.gl_atom_list = GLTLSAtomList(
            atm_U_attr  = "Utls",
            atom_origin = self.tls_group.origin)

        for atm, Utls in self.tls_group.iter_atm_Ucalc():
            atm.Utls = Utls
            self.gl_atom_list.atom_list.append(atm)

        self.gl_atom_list.glo_set_properties_id("gl_atom_list")
        self.glo_add_child(self.gl_atom_list)

        self.glo_link_child_property(
            "L_eigen_vec", "gl_atom_list", "L_eigen_vec")
        self.glo_link_child_property(
            "L1", "gl_atom_list", "L1")
        self.glo_link_child_property(
            "L2", "gl_atom_list", "L2")
        self.glo_link_child_property(
            "L3", "gl_atom_list", "L3")
        self.glo_link_child_property(
            "atom_color", "gl_atom_list", "color")
        self.glo_link_child_property(
            "atom_line_width", "gl_atom_list", "line_width")
        self.glo_link_child_property(
            "U", "gl_atom_list", "U")
        self.glo_link_child_property(
            "U_color", "gl_atom_list", "U_color")
        self.glo_link_child_property(
            "symmetry", "gl_atom_list", "symmetry")
        self.glo_link_child_property(
            "trace", "gl_atom_list", "trace")
 
        ## initalize properties
        self.glo_add_update_callback(self.tls_update_cb)
        self.glo_init_properties(
            origin      = self.tls_group.origin,
            L_eigen_vec = eigen_vectors,
            **args)

    def glo_install_properties(self):
        GLDrawListContainer.glo_install_properties(self)

        self.glo_add_property(
            { "name":        "L_eigen_vec",
              "desc":        "L Eigen Vectors", 
              "type":        "array(3,3)",
              "default":     identity(3),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1",
              "desc":        "L1 Rotation", 
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L2",
              "desc":        "L2 Rotation", 
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L3",
              "desc":        "L3 Rotation", 
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })        
        self.glo_add_property(
            { "name":        "time",
              "desc":        "Simulation Time", 
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":       "tensor_line_width",
              "desc":       "Line Width of Tensors",
              "type":       "float",
              "default":    2.0,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":        "TLS_visible",
              "desc":        "Show TLS Tensors",
              "type":        "boolean",
              "default":     True,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "scale",
              "desc":        "Tensor Display Scaling Factor",
              "type":        "float",
              "default":     1.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "T_color",
              "desc":        "T Tensor Color",
              "type":        "color",
              "default":     (0.0, 1.0, 0.0),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":       "L_color",
              "desc":       "L Tensor Color",
              "type":       "color",
              "default":    (1.0, 0.0, 0.0),
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "S_color",
              "desc":       "S Tensor Color",
              "type":       "color",
              "default":    (0.0, 0.0, 1.0),
              "action":     "recompile" })
        self.glo_add_property(
            { "name":        "CA_line_visible",
              "desc":        "Show Lines to C-Alpha Atoms",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "fan_visible",
              "desc":        "Show Fans to Backbone",
              "type":        "boolean",
              "default":     True,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "fan_opacity",
              "desc":        "Set Backbone Fan Transparancy",
              "type":        "float",
              "default":     0.60,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "fan_material",
              "desc":        "Set Backbone Fan Material",
              "type":        "material",
              "default":     "emerald",
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "atom_color",
              "desc":        "Atom Color",
              "type":        "color",
              "default":     (1.0, 1.0, 1.0),
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "atom_line_width",
              "desc":        "Atom Line Width",
              "type":        "float",
              "default":     2.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "U",
              "desc":        "Show ADP Utls Axes",
              "type":        "boolean",
              "default":     False,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "symmetry",
              "desc":        "Show Symmetry Related Molecules",
              "type":        "boolean",
              "default":     True,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "U_color",
              "desc":        "ADP Axes Color",
              "type":        "color",
              "default":     (0.0, 1.0, 0.0),
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "trace",
              "desc":        "CA Backbone Trace", 
              "type":        "boolean",
              "default":     False,
              "action":      "recompile" })
         
    def tls_update_cb(self, updates, actions):
        if "time" in updates:
            self.update_time()

    def gl_call_list(self):
        gl_struct = self.glo_get_glstructure()
        if gl_struct==None or self.properties["symmetry"]==False:
            glPushMatrix()

            axes = self.properties["L_eigen_vec"]
            glRotatef(self.properties["L1"], *axes[0])
            glRotatef(self.properties["L2"], *axes[1])
            glRotatef(self.properties["L3"], *axes[2])
            
            GLDrawList.gl_call_list(self)

            glPopMatrix()
            return

        for symop in gl_struct.iter_orth_symops():
            glPushMatrix()

            glTranslatef(*-self.properties["origin"])

            glMultMatrixf(
                (symop.R[0,0], symop.R[1,0], symop.R[2,0], 0.0,
                 symop.R[0,1], symop.R[1,1], symop.R[2,1], 0.0,
                 symop.R[0,2], symop.R[1,2], symop.R[2,2], 0.0,
                 symop.t[0],   symop.t[1],   symop.t[2],   1.0) )
            
            glTranslatef(*self.properties["origin"])

            axes = self.properties["L_eigen_vec"]
            glRotatef(self.properties["L1"], *axes[0])
            glRotatef(self.properties["L2"], *axes[1])
            glRotatef(self.properties["L3"], *axes[2])
            
            GLDrawList.gl_call_list(self)

            glPopMatrix()

    def gl_draw(self):
        ## draw TLS axes
        if self.properties["TLS_visible"]==True:
            self.draw_tensors()

        ## draw a line from the COR (center of reaction)
        ## to all CA atoms in the TLS group
        if self.properties["CA_line_visible"]==True:
            self.draw_CA_lines()

        ## draw a transparent fan from the COR to all backbone atoms
        ## in the TLS group
        if self.properties["fan_visible"]==True:
            self.draw_fan()

    def draw_CA_lines(self):
        glDisable(GL_LIGHTING)
        glColor3f(0.5, 0.5, 0.5)
        glLineWidth(1.0)
        
        for atm in self.tls_group:
            if atm.name in ["CA"]:
                glBegin(GL_LINES)
                glVertex3f(0.0, 0.0, 0.0)
                glVertex3f(*atm.position - self.properties["origin"])
                glEnd()

        glEnable(GL_LIGHTING)

    def draw_fan(self):
        glEnable(GL_LIGHTING)

        gl_material = GL_MATERIALS_DICT[self.properties["fan_material"]]
        alpha = self.properties["fan_opacity"]

        r, g, b = gl_material["GL_AMBIENT"]
        glMaterialfv(GL_FRONT, GL_AMBIENT, (r, g, b, alpha))
        r, g, b = gl_material["GL_DIFFUSE"]
        glMaterialfv(GL_FRONT, GL_DIFFUSE, (r, g, b, alpha))
        r, g, b = gl_material["GL_SPECULAR"]
        glMaterialfv(GL_FRONT, GL_SPECULAR, (r, g, b, alpha))
        glMaterialf(GL_FRONT, GL_SHININESS, gl_material["GL_SHININESS"]*128.0)
        
        glBegin(GL_TRIANGLE_FAN)
        glVertex3f(0.0, 0.0, 0.0)
        
        for atm in self.tls_group:
            if atm.name in ["N", "CA", "C"]:
                glVertex3f(*atm.position - self.properties["origin"])

        glEnd()

    def draw_tensors(self):
        """Draw tensor axis.
        """
        scale = self.properties["scale"]

        glDisable(GL_LIGHTING)
        glLineWidth(self.properties["tensor_line_width"])

        ## T: units (A^2)
        glColor3f(*self.properties["T_color"])
        (eigen_values, eigen_vectors) = eigenvectors(self.tls_group.T)
        
        for i in range(3):
            amplitude = 1.414 * math.sqrt(abs(eigen_values[i]))
            scaled_amplitude = scale * amplitude
            v = scaled_amplitude * array(eigen_vectors[i])

            glBegin(GL_LINES)
            glVertex3f(*-v)
            glVertex3f(*v)
            glEnd()

        ## L: units (RAD^2)
        glColor3f(*self.properties["L_color"])
        (eigen_values, eigen_vectors) = eigenvectors(self.tls_group.L)
        
        for i in range(3):
            amplitude = 1.414 * rad2deg * math.sqrt(abs(eigen_values[i]))
            scaled_amplitude = scale * amplitude
            v = scaled_amplitude * array(eigen_vectors[i])

            glBegin(GL_LINES)
            glVertex3f(*-v)
            glVertex3f(*v)
            glEnd()

        ## S: units (A*RAD)
        glColor3f(*self.properties["S_color"])
        (eigen_values, eigen_vectors) = eigenvectors(self.tls_group.S)
        
        for i in range(3):
            amplitude = 1.414 * rad2deg * abs(eigen_values[i])
            scaled_amplitude = scale * amplitude
            v = scaled_amplitude * array(eigen_vectors[i])

            glBegin(GL_LINES)
            glVertex3f(*-v)
            glVertex3f(*v)
            glEnd()

    def update_time(self):
        """Changes the time of the TLS group simulating harmonic motion.
        """
        L1_peak = 1.414 * math.sqrt(abs(self.evalL[0]*rad2deg2))
        L2_peak = 1.414 * math.sqrt(abs(self.evalL[1]*rad2deg2))
        L3_peak = 1.414 * math.sqrt(abs(self.evalL[2]*rad2deg2))

        sin_tm = math.sin(3.0 * self.properties["time"] * 2 * math.pi)

        L1 = L1_peak * sin_tm 
        L2 = L2_peak * sin_tm
        L3 = L3_peak * sin_tm

##         Spc = self.calcs["S'^"] * rad2deg

##         dSp = array(
##             [ (Lx * Spc[0,0]) + (Ly * Spc[1,0]) + (Lz * Spc[2,0]),
##               (Lx * Spc[0,1]) + (Ly * Spc[1,1]) + (Lz * Spc[2,1]),
##               (Lx * Spc[0,2]) + (Ly * Spc[1,2]) + (Lz * Spc[2,2]) ])

##         dS = matrixmultiply(transpose(self.axes), dSp)

        dS = zeros(3)
        origin = array(self.tls_group.origin + dS)
        self.glo_update_properties(L1=L1, L2=L2, L3=L3)


class GLFragment(GLDrawListContainer):
    def __init__(self, **args):
        self.fragment = args["fragment"]

    def glo_install_properties(self):
        pass

    def is_amino_acid(self):
        return isinstance(self.fragment, AminoAcidResidue)

    def is_nucleic_acid(self):
        return isinstance(self.fragment, NucleicAcidResidue)

    def is_water(self):
        return self.fragment.is_water()


class GLChain(GLDrawListContainer):
    def __init__(self, **args):
        GLDrawListContainer.__init__(self, **args)
        self.chain          = args["chain"]

        self.aa_main_chain  = GLAtomList()
        self.aa_side_chain  = GLAtomList()
        self.dna_main_chain = GLAtomList()
        self.dna_side_chain = GLAtomList()
        self.water          = GLAtomList()
        self.hetatm         = GLAtomList()

        ## separate atoms in the chain into groups
        for frag in self.chain.iter_fragments():
            if isinstance(frag, AminoAcidResidue):
                for atm in frag.iter_atoms():
                    if atm.name in ["C", "O", "CA", "N"]:
                        self.aa_main_chain.atom_list.append(atm)
                    else:
                        self.aa_side_chain.atom_list.append(atm)

            elif isinstance(frag, NucleicAcidResidue):
                for atm in frag.iter_atoms():
                    self.dna_main_chain.atom_list.append(atm)

            elif frag.is_water():
                for atm in frag.iter_atoms(): 
                    self.water.atom_list.append(atm)

            else:
                for atm in frag.iter_atoms():
                    self.hetatm.atom_list.append(atm)

        if len(self.aa_main_chain.atom_list)>0:
            self.aa_main_chain.glo_set_properties_id("aa_main_chain")
            self.glo_add_child(self.aa_main_chain)
            self.glo_link_child_property(
                "aa_main_chain_visible", "aa_main_chain", "visible")
        else:
            self.aa_main_chain = None

        if len(self.aa_side_chain.atom_list)>0:
            self.aa_side_chain.glo_set_properties_id("aa_side_chain")
            self.glo_add_child(self.aa_side_chain)
            self.glo_link_child_property(
                "aa_side_chain_visible", "aa_side_chain", "visible")
        else:
            self.aa_side_chain = None

        if len(self.dna_main_chain.atom_list)>0:
            self.dna_main_chain.glo_set_properties_id("dna_main_chain")
            self.glo_add_child(self.dna_main_chain)
            self.glo_link_child_property(
                "dna_main_chain_visible", "dna_main_chain", "visible")
        else:
            self.dna_main_chain = None

        if len(self.dna_side_chain.atom_list)>0:
            self.dna_side_chain.glo_set_properties_id("dna_side_chain")
            self.glo_add_child(self.dna_side_chain)
            self.glo_link_child_property(
                "dna_side_chain_visible", "dna_side_chain", "visible")
        else:
            self.dna_side_chain = None

        if len(self.hetatm.atom_list)>0:
            self.hetatm.glo_set_properties_id("hetatm")
            self.glo_add_child(self.hetatm)
            self.glo_link_child_property(
                "hetatm_visible", "hetatm", "visible")
        else:
            self.hetatm = None

        if len(self.water.atom_list)>0:
            self.water.glo_set_properties_id("water")
            self.glo_add_child(self.water)
            self.glo_link_child_property(
                "water_visible", "water", "visible")
        else:
            self.water = None

        for gl_object in self.glo_iter_children():
            gl_object_id = gl_object.glo_get_properties_id()
            self.glo_link_child_property("color", gl_object_id, "color")
            self.glo_link_child_property("U", gl_object_id, "U")

        self.glo_init_properties(**args)

    def glo_install_properties(self):
        GLDrawListContainer.glo_install_properties(self)

        self.glo_add_property(
            { "name":      "aa_main_chain_visible",
              "type":      "boolean",
              "default":   True,
              "action":    "redraw" })

        self.glo_add_property(
            { "name":      "aa_side_chain_visible",
              "type":      "boolean",
              "default":   True,
              "action":    "redraw" })

        self.glo_add_property(
            { "name":      "dna_main_chain_visible",
              "type":      "boolean",
              "default":   True,
              "action":    "redraw" })

        self.glo_add_property(
            { "name":      "dna_side_chain_visible",
              "type":      "boolean",
              "default":   True,
              "action":    "redraw" })
        
        self.glo_add_property(
            { "name":      "hetatm_visible",
              "type":      "boolean",
              "default":   True,
              "action":    "redraw" })
        
        self.glo_add_property(
            { "name":      "water_visible",
              "type":      "boolean",
              "default":   True,
              "action":    "redraw" })
        
        self.glo_add_property(
             { "name":      "color",
               "desc":      "Atom Color",
               "type":      "color",
               "default":   None,
               "action":    "redraw" })

        self.glo_add_property(
            { "name":      "U",
              "desc":      "Show U Axes",
              "type":      "matrix(3,3)",
              "default":   False,
              "action":    "redraw" })


class GLStructure(GLDrawListContainer):
    def __init__(self, **args):
        GLDrawListContainer.__init__(self, **args)
        self.struct = args["struct"]
        
        ## add GLObject children

        ## structure axes
        self.gl_axes = GLAxes()
        self.gl_axes.glo_set_properties_id("gl_axes")
        self.glo_add_child(self.gl_axes)
        self.glo_link_child_property(
            "axes_visible", "gl_axes", "visible")

        ## unit cell
        self.gl_unit_cell = GLUnitCell(unit_cell=self.struct.unit_cell)
        self.gl_unit_cell.glo_set_properties_id("gl_unit_cell")
        self.glo_add_child(self.gl_unit_cell)
        self.glo_link_child_property(
            "unit_cell_visible", "gl_unit_cell", "visible")

        ## GLChains 
        for chain in self.struct.iter_chains():
            gl_chain = GLChain(chain=chain)
            gl_chain.glo_set_properties_id(str(chain))
            self.glo_add_child(gl_chain)

            self.glo_link_child_property(
                "aa_main_chain_visible", str(chain), "aa_main_chain_visible")
            self.glo_link_child_property(
                "aa_side_chain_visible", str(chain), "aa_side_chain_visible")
            self.glo_link_child_property(
                "dna_main_chain_visible", str(chain), "dna_main_chain_visible")
            self.glo_link_child_property(
                "dna_side_chain_visible", str(chain), "dna_side_chain_visible")
            self.glo_link_child_property(
                "hetatm_visible", str(chain), "hetatm_visible")
            self.glo_link_child_property(
                "water_visible", str(chain), "water_visible")
            self.glo_link_child_property(
                "color", str(chain), "color")

        ## init properties
        self.glo_init_properties(**args)

    def glo_install_properties(self):
        GLDrawListContainer.glo_install_properties(self)

        self.glo_add_property(
            { "name":     "axes_visible",
              "desc":     "Show Cartesian Axes",
              "type":     "boolean",
              "default":  True,
              "action":  "redraw" })

        self.glo_add_property(
            { "name":     "unit_cell_visible",
              "desc":     "Show Unit Cell",
              "type":     "boolean",
              "default":  True,
              "action":   "redraw" })
        
        self.glo_add_property(
            { "name":     "aa_main_chain_visible",
              "desc":     "Show Amino Acid Main Chain Atoms",
              "default":  True,
              "action":   "redraw" })

        self.glo_add_property(
            { "name":     "aa_side_chain_visible",
              "desc":     "Show Amino Acid Side Chain Atoms",
              "default":  True,
              "action":   "redraw" })

        self.glo_add_property(
            { "name":              "dna_main_chain_visible",
              "default":           True,
              "action":            "redraw" })

        self.glo_add_property(
            { "name":              "dna_side_chain_visible",
              "default":           True,
              "action":            "redraw" })
        
        self.glo_add_property(
            { "name":              "hetatm_visible",
              "default":           True,
              "action":            "redraw" })
        
        self.glo_add_property(
            { "name":              "water_visible",
              "default":           True,
              "action":            "redraw" })
        
        self.glo_add_property(
            { "name":              "color",
              "default":           None,
              "action":            "redraw" })
        
        self.glo_add_property(
            { "name":              "U",
              "default":           False,
              "action":            "redraw" })

    def color_by_residue_chem_type(self, atm):
        """GLAtomList color callback for coloring a structure by
        residue chical type: aliphatic, aromatic, sulfer-containing,
        alchols, acids, bases, amides.
        """
        mon = self.struct.library.get_monomer(atm.res_name)
        if mon.is_amino_acid()==False or mon.chem_type=="":
            return (1.0, 1.0, 1.0)

        chem_type_color_dict = {
            "aliphatic":         (0.50, 0.50, 0.50),
            "aromatic":          (0.75, 0.75, 0.75),
            "sulfer-containing": (0.00, 1.00, 0.00),
            "alchols":           (0.75, 0.75, 1.00),
            "acids":             (1.00, 0.00, 0.00),
            "bases":             (0.00, 0.00, 1.00),
            "amides":            (1.00, 0.75, 0.75)}
        
        return chem_type_color_dict[mon.chem_type]

    def iter_orth_symops(self):
        """Iterate orthogonal-space symmetry operations useful for
        displaying symmetry-equivelant molecules without having to
        calculate new draw lists.
        """
        if hasattr(self, "orth_symop_cache"):
            for symop in self.orth_symop_cache:
                yield symop
        else:
            self.orth_symop_cache = []
            uc = self.struct.unit_cell

            for symop in uc.iter_struct_orth_symops(self.struct):
                self.orth_symop_cache.append(symop)
                yield symop
        

class GLViewer(GLObject):
    """This class renders a list of GLDrawList (or subclasses of) onto
    the given glcontext and gldrawable objects.  The glcontext and gldrawable
    must be created by the underling GUI toolkit, or perhaps the GLUT
    libraries.  This class is completely platform and tookit independent
    once it is passed the glcontext and gldrawable.  The design of this
    class and the associated GLDrawList classes incorporates some basic
    OpenGL drawing optimizations.  The GLDrawList objects are drawn and
    compiled into OpenGL draw lists, and have their own
    transformation/rotation operators WRT the GLViewer origin, allowing
    each GLDrawList to be redrawn very quickly as long as it moves as a
    rigid body.
    """
    def __init__(self):
        GLObject.__init__(self)
        
        ## position and rotation of viewer window
        self.xpos = 0.0
        self.ypos = 0.0
        self.zpos = -50.0
        self.rotx = 0.0
        self.roty = 0.0
        self.rotz = 0.0

    def add_draw_list(self, draw_list):
        """Append a GLDrawList.
        """
        assert isinstance(draw_list, GLDrawList)
        self.glo_add_child(draw_list)
        self.gl_redraw()

    def remove_draw_list(self, draw_list):
        """Remove a GLDrawList.
        """
        assert isinstance(draw_list, GLDrawList)

        ## delete the compiled GL draw lists
        draw_list.gl_delete_list()

        ## remove and trigger redraw
        draw_list.glo_remove()
        self.gl_redraw()

    def get_gl_context(self):
        """Implement in subclass to return the gl_context.
        """
        pass

    def get_gl_drawable(self):
        """Implement in subclass to reutrn the gl_drawable.
        """
        pass

    def gl_redraw(self):
        """Redraw.
        """
        pass

    def gl_init(self):
        """Called once to initalize the GL scene before drawing.
        """
        gl_drawable = self.get_gl_drawable()
        gl_context  = self.get_gl_context()
        if gl_drawable==None or gl_context==None:
            return
        
        if not gl_drawable.gl_begin(gl_context):
            return        

        ambient        = [0.0, 0.0, 0.0, 1.0]
        diffuse        = [1.0, 1.0, 1.0, 1.0]
        specular       = [1.0, 1.0, 1.0, 1.0]
        position       = [0.0, 3.0, 3.0, 0.0]
        lmodel_ambient = [0.5, 0.5, 0.5, 1.0]
        local_view     = [0.0]
        
        glLightfv(GL_LIGHT0, GL_AMBIENT, ambient)
        glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse)
        glLightfv(GL_LIGHT0, GL_POSITION, position)
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient)
        glLightModelfv(GL_LIGHT_MODEL_LOCAL_VIEWER, local_view)
        
        glFrontFace(GL_CW)
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        glEnable(GL_AUTO_NORMAL)
        glEnable(GL_NORMALIZE)
        glEnable(GL_DEPTH_TEST)	


	glDepthFunc(GL_LESS)
	glEnable(GL_DEPTH_TEST)


        ## FOG
        #glEnable(GL_FOG)
        #glFogf(GL_FOG_MODE, GL_EXP)
        #glFogf(GL_FOG_DENSITY, 0.01)
        #glFogf(GL_FOG_START, 20.0)

        ## ANTI-ALIASING
        glEnable(GL_LINE_SMOOTH)
        glEnable(GL_POINT_SMOOTH)
        glEnable(GL_POLYGON_SMOOTH)

        ## ALPHA BLENDING
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

        gl_drawable.gl_end()
    
    def gl_resize(self, width, height):
        """Called to set the size of the OpenGL window this class is
        drawing on.
        """
        gl_drawable = self.get_gl_drawable()
        gl_context  = self.get_gl_context()
        if gl_drawable==None or gl_context==None:
            return
        
	if not gl_drawable.gl_begin(gl_context):
            return
	
	glViewport(0, 0, width, height)
	glMatrixMode(GL_PROJECTION)
	glLoadIdentity()

	if width > height:
            w = float(width) / float(height)
            glFrustum(-w, w, -1.0, 1.0, 3.0, 5000.0)
	else:
            h = float(height) / float(width)
            glFrustum(-1.0, 1.0, -h, h, 3.0, 5000.0)
	
	glMatrixMode(GL_MODELVIEW)
	gl_drawable.gl_end()


    def ROT(self):
        alpha = - self.rotx * deg2rad
        beta  = - self.roty * deg2rad
        gamma = - self.rotz * deg2rad

        R1 = ( math.cos(alpha) * identity(3) +

               math.sin(alpha) * array([[0.0,  0.0,  0.0],
                                        [0.0,  0.0, -1.0],
                                        [0.0,  1.0,  0.0]]) +

             (1.0 - math.cos(alpha)) * array([[1.0, 0.0, 0.0],
                                              [0.0, 0.0, 0.0],
                                              [0.0, 0.0, 0.0]]) )
        
        R2 = ( math.cos(beta) * identity(3) +
               
               math.sin(beta) * array([[ 0.0,  0.0, 1.0],
                                       [ 0.0,  0.0, 0.0],
                                       [-1.0,  0.0, 0.0]]) +

             (1.0 - math.cos(beta)) * array([[0.0, 0.0, 0.0],
                                             [0.0, 1.0, 0.0],
                                             [0.0, 0.0, 0.0]]) )
        
               
        R3 = ( math.cos(gamma) * identity(3) +
               
               math.sin(gamma) * array([[ 0.0, -1.0, 0.0],
                                        [ 1.0,  0.0, 0.0],
                                        [ 0.0,  0.0, 0.0]]) +

             (1.0 - math.cos(gamma)) * array([[0.0, 0.0, 0.0],
                                              [0.0, 0.0, 0.0],
                                              [0.0, 0.0, 1.0]]) )

        return matrixmultiply(R3, matrixmultiply(R2, R1))

    def gl_translate(self, x, y, z):
        """
        """
        delta = matrixmultiply(self.ROT(), array([x,y,z]))

        self.xpos += delta[0]
        self.ypos += delta[1]
        self.zpos += delta[2]

        
    def gl_render(self):
        """Draw all GLDrawList objects onto the given glcontext/gldrawable.
        If the GLDrawList objects are not yet compiled into OpenGL draw
        lists, they will be compiled while they are drawn, since this is
        a useful optimization.
        """
        gl_drawable = self.get_gl_drawable()
        gl_context  = self.get_gl_context()
        if gl_drawable==None or gl_context==None:
            return
        
	if not gl_drawable.gl_begin(gl_context):
            return

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
	glLoadIdentity()

        glTranslatef(0.0, 0.0, -10.0)

	glRotatef(self.rotx, 1.0, 0.0, 0.0)
	glRotatef(self.roty, 0.0, 1.0, 0.0)
        glRotatef(self.rotz, 0.0, 0.0, 1.0)

        glTranslatef(self.xpos, self.ypos, self.zpos + 10.0)

        for draw_list in self.glo_iter_children():
            draw_list.gl_render()
            
	if gl_drawable.is_double_buffered():
            gl_drawable.swap_buffers()
	else:
            glFlush()
        
        gl_drawable.gl_end()
