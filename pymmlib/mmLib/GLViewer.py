## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""OpenGL rendering classes.
"""
from __future__  import generators

import random

from OpenGL.GL      import *
from OpenGL.GLU     import *
from OpenGL.GLUT    import *
from mmTypes        import *
from Structure      import *
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

    def glo_name(self):
        if self.__globject_properties_id!=None:
            return "%s(%s)" % (self.__class__.__name__,
                               self.__globject_properties_id)
        else:
            return self.__class__.__name__

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
            gl_viewer.glv_redraw()

    def glo_get_glviewer(self):
        """Returns the root GLViewer object.
        """
        return self.glo_get_root()

    def glo_get_glstructure(self):
        """Returns the parent GLStructure object, or None if the GLObject
        is not a child of a GLStructure.
        """
        gl_object = self
        while gl_object!=None and isinstance(gl_object, GLStructure)==False:
            gl_object = gl_object.__globject_parent
        return gl_object


class OpenGLRenderMethods(object):
    """OpenGL renderer methods.  Eventually, all OpenGL rendering will
    be done through methods in this object.
    """
    def glr_set_material_name(self, material_name, opacity):
        """Sets the given named material with the given opacity
        """
        gl_material = GL_MATERIALS_DICT[material_name]
        alpha       = opacity

        r, g, b = gl_material["GL_AMBIENT"]
        glMaterialfv(GL_FRONT, GL_AMBIENT, (r, g, b, alpha))
        r, g, b = gl_material["GL_DIFFUSE"]
        glMaterialfv(GL_FRONT, GL_DIFFUSE, (r, g, b, alpha))
        r, g, b = gl_material["GL_SPECULAR"]
        glMaterialfv(GL_FRONT, GL_SPECULAR, (r, g, b, alpha))
        glMaterialf(GL_FRONT, GL_SHININESS, gl_material["GL_SHININESS"]*128.0)

    def glr_set_material_rgb(self, r, g, b, opacity):
        """Creates a stock rendering material colored according to the given
        RGB values.
        """
        ambient  = array([r, g, b, opacity])
        diffuse  = ambient
        specular = ambient * 1.25
                
	glMaterial(GL_FRONT, GL_AMBIENT,   ambient)
	glMaterial(GL_FRONT, GL_DIFFUSE,   diffuse)
	glMaterial(GL_FRONT, GL_SPECULAR,  specular)
	glMaterial(GL_FRONT, GL_SHININESS, 50.0)

    def glr_set_material_atom(self, atm, opacity):
        """Sets a material for rendering a atom which uses the atom's
        element color for the material selection.
        """
        elem = atm.get_structure().library.get_element(atm.element)
        if elem!=None:
            r, g, b = elem.color
        else:
            r, g, b = 1.0, 1.0, 1.0

        self.glr_set_material_rgb(r, g, b, opacity)
    
    def glr_axis(self, position, axis, radius):
        """Draw a vector axis using the current set material at position
        with the given radius.
        """
        axis_cap_length     = 0.01 * length(axis)
        axis_length         = length(axis) - 2*axis_cap_length
        axis_cap_tip_radius = 0.75 * radius

        glEnable(GL_LIGHTING)
        quad = gluNewQuadric()
	gluQuadricNormals(quad, GLU_SMOOTH)

        glPushMatrix()

        ## rotation matrix to align the current OpenGL coordinate
        ## system such that the vector axis is along the z axis
        R = inverse(rmatrixz(axis))
        glMultMatrixf(
            (R[0,0],           R[1,0],      R[2,0], 0.0,
             R[0,1],           R[1,1],      R[2,1], 0.0,
             R[0,2],           R[1,2],      R[2,2], 0.0,
             position[0], position[1], position[2], 1.0) )

        gluDisk(quad, 0.0, axis_cap_tip_radius, 10, 5)
        gluCylinder(quad, axis_cap_tip_radius, radius, axis_cap_length, 10, 5)
        glTranslatef(0.0, 0.0, axis_cap_length)
	gluCylinder(quad, radius, radius, axis_length , 10, 1) 
        glTranslatef(0.0, 0.0, axis_length)
        gluCylinder(quad, radius, axis_cap_tip_radius, axis_cap_length, 10, 5)
        glTranslatef(0.0, 0.0, axis_cap_length)
        gluDisk(quad, 0.0, axis_cap_tip_radius, 10, 5)

        glPopMatrix()

    def glr_cpk(self, position, radius, quality):
        """Draw a atom as a CPK sphere.
        """
        glEnable(GL_LIGHTING)
        glPushMatrix()
        glTranslatef(*position)
        glutSolidSphere(radius, quality, quality)
        glPopMatrix()

    def glr_cross(self, position, color, line_width):
        """Draws atom with a cross of lines.
        """
        glDisable(GL_LIGHTING)
        glColor3f(*color)
        glLineWidth(line_width)

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

    def glr_U_axes(self, position, U, color, line_width):
        """Draw the anisotropic axies of the atom with R.M.S.
        magnitude.
        """
        eigen_values, eigen_vectors = eigenvectors(U)

        v0_peak = math.sqrt(abs(eigen_values[0]))
        v1_peak = math.sqrt(abs(eigen_values[1]))
        v2_peak = math.sqrt(abs(eigen_values[2]))
        
        v0 = eigen_vectors[0] * v0_peak
        v1 = eigen_vectors[1] * v1_peak
        v2 = eigen_vectors[2] * v2_peak

        glDisable(GL_LIGHTING)
        glColor3f(*color)
        glLineWidth(line_width)
        
        glBegin(GL_LINES)

        glVertex3f(*position - v0)
        glVertex3f(*position + v0)
        glVertex3f(*position - v1)
        glVertex3f(*position + v1)
        glVertex3f(*position - v2)
        glVertex3f(*position + v2)

        glEnd()


class GLDrawList(GLObject, OpenGLRenderMethods):
    """Fundamental OpenGL rigid entity.
    """
    def __init__(self, **args):
        GLObject.__init__(self, **args)
        self.gl_name   = None
        self.glo_add_update_callback(self.gldl_update_cb)
        self.glo_init_properties(**args)

    def gldl_update_cb(self, updates, actions):
        if "recompile" in actions:
            self.gldl_delete_list()
            self.glo_redraw()
        if "redraw" in actions:
            self.glo_redraw()

    def glo_install_properties(self):
        self.glo_add_property(
            { "name" :      "visible",
              "desc":       "Visible",
              "catagory":   "Show/Hide",
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
                
    def gldl_push_matrix(self):
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

    def gldl_pop_matrix(self):
        """Pop the roatated/translated position.
        """
        glPopMatrix()

    def gldl_render(self):
        """Compile or force a recompile of this object's gl_draw list, and
        render the scene.  Rendering the scene can be bypassed if
        this method is called with render = False.
        """
        if self.properties["visible"]==False:
            return
        
        ## pop the matrix pushed by properties["origin"], properties["rot_x"]
        ## properties["rot_y"], and properties["rot_z"]
        self.gldl_push_matrix()

        ## render this GLDrawList, compiling a new OpenGL draw list if needed
        if self.gl_name==None:
            self.gl_name = glGenLists(1)
            glNewList(self.gl_name, GL_COMPILE)
            self.gldl_draw()
            glEndList()

        ## support multiple rendering images by implementing class
        ## iterators gldl_iter_multidraw_all() for multiple
        ## rendering iterations of the GLDrawList and all its children,
        ## or gldl_iter_multidraw_self() for multiple images of just
        ## this GLDrawList, rendering the children just once
        for draw_flag_multi in self.gldl_iter_multidraw_all():

            for draw_flag_self in self.gldl_iter_multidraw_self():
                glCallList(self.gl_name)

            ## render first-level children of this GLDrawList
            ## which, in turn, will render their children
            for draw_list in self.glo_iter_children():
                draw_list.gldl_render()

        self.gldl_pop_matrix()

    def gldl_iter_multidraw_all(self):
        """When implemented as a iterator in a subclass, each time yield
        is invoked the GLDrawList and all its decendants will be rendered
        from whatever OpenGL coordinate system is set in the iterator.
        """
        yield True

    def gldl_iter_multidraw_self(self):
        """Similar to gldl_iter_multidraw_all, but only this GLDrawList is
        rendered.  The decendant GLDrawLists are rendered normally.
        """
        yield True

    def gldl_delete_list(self):
        """Delete the currently compiled glList for this object.  This forces
        a recompile of the OpenGL draw list next time gl_call_list_render()
        is invoked.  It also should be called before a GLDrawList is
        dereferenced so that unused compiled glLists are not left in the
        OpenGL libraries.  This method also deletes the draw lists for all
        children of the GLDrawList.
        """
        if self.gl_name != None:
            glDeleteLists(self.gl_name, 1)
            self.gl_name = None
            
        for draw_list in self.glo_iter_children():
            draw_list.gldl_delete_list()

    def gldl_draw(self):
        """Implement in subclass to draw somthing.
        """
        pass


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
              "desc":       "Axis Length",
              "type":       "float",
              "default":    20.0,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "line_width",
              "desc":       "Axis Radius",
              "type":       "float",
              "default":    0.1,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "color_x",
              "desc":       "X Axis Color",
              "catagory":   "Colors",
              "type":       "color",
              "default":    (1.0, 0.0, 0.0),
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "color_y",
              "desc":       "Y Axis Color",
              "catagory":   "Colors",
              "type":       "color",
              "default":    (0.0, 1.0, 0.0),
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "color_z",
              "desc":       "Z Axis Color",
              "catagory":   "Colors",
              "type":       "color",
              "default":    (0.0, 0.0, 1.0),
              "action":     "recompile" })

    def gldl_draw(self):
        line_length = self.properties["line_length"]
        line_width  = self.properties["line_width"]
        
        (r, g, b) = self.properties["color_x"]
        self.glr_set_material_rgb(r, g, b, 1.0)
        self.glr_axis(zeros(3), array([line_length, 0.0, 0.0]), line_width) 

        (r, g, b) = self.properties["color_y"]
        self.glr_set_material_rgb(r, g, b, 1.0)
        self.glr_axis(zeros(3), array([0.0, line_length, 0.0]), line_width) 

        (r, g, b) = self.properties["color_z"]
        self.glr_set_material_rgb(r, g, b, 1.0)
        self.glr_axis(zeros(3), array([0.0, 0.0, line_length]), line_width) 


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
              "desc":       "Line Width",
              "type":       "float",
              "default":    1.0,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "color",
              "desc":       "Line Color",
              "type":       "color",
              "default":    (1.0, 1.0, 1.0),
              "action":     "recompile" })

    def gldl_draw(self):
        self.draw_cell(-1, -1, -1, 0, 0, 0)
        
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


class GLAtomList(GLDrawList):
    """OpenGL renderer for a list of atoms.  Optional arguments iare:
    color, U, U_color.
    """
    def __init__(self, **args):
        GLDrawList.__init__(self, **args)

        self.atom_list      = AtomList()
        self.el_color_cache = {}

        self.glo_init_properties(**args)

    def glo_install_properties(self):
        GLDrawList.glo_install_properties(self)


        ## global
        self.glo_add_property(
            { "name":      "symmetry",
              "desc":      "Show Symmetry Equivelants",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   False,
              "action":    "redraw" })
        self.glo_add_property(
            { "name":       "atom_origin",
              "desc":       "Atom Calculation Origin",
              "type":       "array(3)",
              "hidden":     True,
              "default":    None,
              "action":     "recompile" })
        
        self.glo_add_property(
            { "name":      "color",
              "desc":      "Atom Color",
              "catagory":  "Colors",
              "type":      "color",
              "default":   None,
              "action":    "recompile" })
        self.glo_add_property(
            { "name":      "color_func",
              "type":      "function",
              "hidden":    True,
              "default":   None,
              "action":    "recompile" })

        ## lines
        self.glo_add_property(
            { "name":       "lines",
              "desc":       "Draw Atom Bond Lines",
              "catagory":   "Show/Hide",
              "type":       "boolean",
              "default":    True,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "line_width",
              "desc":       "Bond Line Drawing Width",
              "catagory":   "Bond Lines",
              "type":       "float",
              "default":    1.0,
              "action":     "recompile" })

        ## cpk
        self.glo_add_property(
            { "name":       "cpk",
              "desc":       "Draw CPK Spheres",
              "catagory":   "Show/Hide",
              "type":       "boolean",
              "default":    False,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "cpk_scale_radius",
              "desc":       "Scale CPK Radius",
              "catagory":   "CPK",
              "type":       "float",
              "default":    1.0,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "cpk_opacity",
              "desc":       "CPK Sphere Opacity",
              "catagory":   "CPK",
              "type":       "float",
              "default":    0.50,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "sphere_quality",
              "desc":       "CPK Sphere Quality",
              "catagory":   "CPK",
              "type":       "integer",
              "default":    12,
              "action":     "recompile" })

        ## trace           
        self.glo_add_property(
            { "name":       "trace",
              "desc":       "Draw Backbone Trace",
              "catagory":   "Show/Hide",
              "type":       "boolean",
              "default":    False,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "trace_line_width",
              "desc":       "Backbone Trace Line Width",
              "catagory":   "Trace",
              "type":       "float",
              "default":    1.0,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "trace_color",
              "desc":       "Backbone Trace Color",
              "catagory":   "Trace",
              "type":       "color",
              "default":    (1.0, 1.0, 1.0),
              "action":     "recompile" })

        ## ADPs
        self.glo_add_property(
            { "name":      "atm_U_attr",
              "type":      "string",
              "hidden":    True,
              "default":   "U",
              "action":    "recompile" })

        self.glo_add_property(
            { "name":      "U",
              "desc":      "Draw ADP Axes",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   False,
              "action":    "recompile" })
        self.glo_add_property(
            { "name":      "U_color",
              "desc":      "Thermal Axes Color",
              "catagory":  "ADP",
              "type":      "color",
              "default":   (1.0, 1.0, 1.0),
              "action":    "recompile" })

    def gldl_iter_multidraw_self(self):
        """Specialized draw list invokation to recycle the draw list for
        symmetry related copies.  Cartesian versions of the symmetry rotation
        and translation operators are generated by GLStructure/UnitCell
        classes.
        """
        if self.properties["symmetry"]==False:
            yield True
            
        else:

            gl_struct = self.glo_get_glstructure()
            if gl_struct==None:
                yield True

            else:
                for symop in gl_struct.iter_orth_symops():
                    glPushMatrix()

                    if self.properties["atom_origin"]!=None:
                        glTranslate(*-self.properties["atom_origin"])

                    glMultMatrixf(
                        (symop.R[0,0], symop.R[1,0], symop.R[2,0], 0.0,
                         symop.R[0,1], symop.R[1,1], symop.R[2,1], 0.0,
                         symop.R[0,2], symop.R[1,2], symop.R[2,2], 0.0,
                         symop.t[0],   symop.t[1],   symop.t[2],   1.0) )

                    if self.properties["atom_origin"]!=None:
                        glTranslate(*self.properties["atom_origin"])
                    
                    yield True
                    glPopMatrix()

    def gldl_draw(self):
        """Draw all selected representations.
        """
        #self.draw_symmetry_debug()
        
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
        """Calculate a position vector with respect to the
        proeprty: atom_origin.
        """
        if self.properties["atom_origin"]!=None:
            pos = pos - self.properties["atom_origin"]
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
        self.glr_set_material_rgb(r, g, b, self.properties["cpk_opacity"])

        self.glr_cpk(
            self.calc_position(atm.position),
            self.properties["cpk_scale_radius"] * radius,
            self.properties["sphere_quality"])

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

    def draw_cross(self, atm):
        """Draws atom with a cross of lines.
        """
        self.glr_cross(
            self.calc_position(atm.position),
            self.get_color(atm),
            self.properties["line_width"])

    def draw_trace(self, color=(1.0, 1.0, 1.0)):
        """Draws trace over all protein backbone atoms.
        """
        glDisable(GL_LIGHTING)
        glColor3f(*self.properties["trace_color"])
        glLineWidth(self.properties["trace_line_width"])
        
        glBegin(GL_LINE_STRIP)

        for atm in self.atom_list:
            if atm.name in ["N", "CA", "C"]:
                glVertex3f(*self.calc_position(atm.position))

        glEnd()

    def draw_U_axes(self, atm):
        """Draw the anisotropic axies of the atom with R.M.S.
        magnitude.
        """
        U = getattr(atm, self.properties["atm_U_attr"], None)
        if U==None:
            return
        
        self.glr_U_axes(
            self.calc_position(atm.position),
            U,
            self.properties["U_color"],
            1.0)

    def draw_symmetry_debug(self):
        """Draws crosses where all the symmetry equivalent atoms should be.
        To make sure the real-space symmetry operations are working.
        """
        gl_struct = self.glo_get_glstructure()
        unit_cell = gl_struct.struct.unit_cell

        for atm in self.atom_list:
            frac = unit_cell.calc_orth_to_frac(atm.position)

            for xyz in unit_cell.space_group.iter_equivalent_positions(frac):

                xyz  = array([xyz[0]%1.0, xyz[1]%1.0, xyz[2]%1.0])
                orth = unit_cell.calc_frac_to_orth(xyz)
                
                self.glr_cross(orth, (1.0, 1.0, 1.0), 1.0) 


class GLTLSScrewRotation(GLDrawList):
    """Handles visualization and animation of a single TLS screw rotation
    axis.  Three of these make up a full TLS group animation.
    """
    def __init__(self, **args):
        GLDrawList.__init__(self, **args)
        self.atom_list = args["atom_list"]

        ## create a unique list of bonds which will be used to
        ## render the TLS surface; this list may be passed in a argument
        ## to avoid multiple calculations for each screw-rotation axis
        try:
            self.bond_list = args["bond_list"]
        except KeyError:
            self.bond_list = []
            in_dict        = {}

            for atm in self.atom_list:
                in_dict[atm] = True

            for atm in self.atom_list:
                for bond in atm.iter_bonds():
                    if in_dict.has_key(bond.get_partner(atm)):
                        self.bond_list.append(bond)

        self.glo_init_properties(**args)

    def glo_install_properties(self):
        GLDrawList.glo_install_properties(self)

        self.glo_add_property(
            { "name":        "symmetry",
              "desc":        "Show Symmetry Equivelants",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "scale_rot",
              "desc":        "Scale Rotation",
              "catagory":    "Show/Hide",
              "type":        "float",
              "default":     1.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "COR",
              "desc":        "TLS Center of Reaction", 
              "catagory":    "TLS Analysis",
              "type":        "array(3)",
              "default":     zeros(3),
              "action":      "recompile" })        
        self.glo_add_property(
            { "name":        "Lx_eigen_val",
              "desc":        "Libration Axis Eigenvalue", 
              "catagory":    "TLS Analysis",
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "Lx_eigen_vec",
              "desc":        "Libration Axis Normalized Eigenvector",  
              "catagory":    "TLS Analysis",
              "type":        "array(3)",
              "default":     zeros(3),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "Lx_rho",
              "desc":        "Libration Axis Translation from COR",  
              "catagory":    "TLS Analysis",
              "type":        "array(3)",
              "default":     zeros(3),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "Lx_pitch",
              "desc":        "Libration screw pitch (DEG/A)",  
              "catagory":    "TLS Analysis",
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "Lx_rot",
              "desc":        "Current Setting for Libration Rotation",
              "catagory":    "Simulation",
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })

    def gldl_iter_multidraw_self(self):
        """Specialized draw list invokation to recycle the draw list for
        symmetry related copies.  Cartesian versions of the symmetry rotation
        and translation operators are generated by GLStructure/UnitCell
        classes.
        """
        if self.properties["symmetry"]==False:
            yield True
            
        else:

            gl_struct = self.glo_get_glstructure()
            if gl_struct==None:
                yield True

            else:
                for symop in gl_struct.iter_orth_symops():
                    glPushMatrix()

                    glMultMatrixf(
                        (symop.R[0,0], symop.R[1,0], symop.R[2,0], 0.0,
                         symop.R[0,1], symop.R[1,1], symop.R[2,1], 0.0,
                         symop.R[0,2], symop.R[1,2], symop.R[2,2], 0.0,
                         symop.t[0],   symop.t[1],   symop.t[2],   1.0) )

                    yield True
                    glPopMatrix()
                    
    def gldl_draw(self):
        self.draw_TLS_surface(
            self.properties["Lx_eigen_vec"],
            self.properties["Lx_eigen_val"],
            self.properties["Lx_rho"],
            self.properties["Lx_pitch"])
            
    def draw_TLS_surface(self, Lx_eigen_vec, Lx_eigen_val, Lx_rho, Lx_pitch):
        """Draws the TLS probability surface for a single non-intersecting
        screw axis.  Lx_eigen_val is the vaiance (mean square deviation MSD)
        of the rotation about the Lx_eigen_vec axis.
        """
        Lx_s          = math.sqrt(abs(Lx_eigen_val * deg2rad2))
        Lx_s         *= self.properties["scale_rot"]
        COR           = self.properties["COR"]
        steps         = 10
        rot_step      = Lx_s / steps
        sq2pi         = math.sqrt(2.0 * math.pi)
        opacity_scale = None

        glEnable(GL_LIGHTING)
        glBegin(GL_QUADS)

        p = 1.0
        
        for step in range(steps):

            step1 = rot_step * step
            step2 = rot_step * (1 + step)

            ## caclulate a opacity proportional to the
            ## gaussian PDF of the rotation
            Lx_s2 = Lx_s * Lx_s 
            p1 = (1.0/(sq2pi*Lx_s))* math.exp(-(step1*step1)/(2.0*Lx_s2))
            p2 = (1.0/(sq2pi*Lx_s))* math.exp(-(step2*step2)/(2.0*Lx_s2))

            dstep = step2 - step1
            p = (p2*dstep) + 0.5*(p1-p2)*dstep

            if opacity_scale==None:
                opacity_scale = 0.6 / p

            opacity = p * opacity_scale

#            print "## step = %f Lx = %f opacity = %f" % (
#                step2, Lx_s, opacity)

            for sign in (-1, 1):

                rot1   = step1 * sign
                rot2   = step2 * sign

                Rstep1 = rmatrixu(Lx_eigen_vec, rot1)
                Rstep2 = rmatrixu(Lx_eigen_vec, rot2)

                Dstep1 = dmatrixu(Lx_eigen_vec, rot1)
                Dstep2 = dmatrixu(Lx_eigen_vec, rot2)

                screw1 = Lx_eigen_vec * rot1 * Lx_pitch
                screw2 = Lx_eigen_vec * rot2 * Lx_pitch

                rho1 = matrixmultiply(Dstep1, -Lx_rho)
                rho2 = matrixmultiply(Dstep2, -Lx_rho)

                rho_screw1 = rho1 + screw1
                rho_screw2 = rho2 + screw2

                for bond in self.bond_list:

                    pos1 = bond.atom1.position - COR
                    pos2 = bond.atom2.position - COR

                    v1 = matrixmultiply(Rstep1, pos1) + rho_screw1
                    v2 = matrixmultiply(Rstep2, pos1) + rho_screw2
                    
                    v3 = matrixmultiply(Rstep2, pos2) + rho_screw2
                    v4 = matrixmultiply(Rstep1, pos2) + rho_screw1
                    
                    v14 = v1 + (0.5 * (v4 - v1))
                    v23 = v2 + (0.5 * (v3 - v2))

                    self.glr_set_material_atom(bond.atom1, opacity)

                    glVertex3f(*v1  + COR)
                    glVertex3f(*v2  + COR)
                    glVertex3f(*v23 + COR)
                    glVertex3f(*v14 + COR)
                    
                    self.glr_set_material_atom(bond.atom2, opacity)

                    glVertex3f(*v3  + COR)
                    glVertex3f(*v4  + COR)
                    glVertex3f(*v14 + COR)
                    glVertex3f(*v23 + COR)

        glEnd()


class GLTLSAtomList(GLAtomList):
    """OpenGL visualizations of TLS group atoms.
    """
    def __init__(self, **args):
        GLAtomList.__init__(self, **args)
        self.atom_list = args["atom_list"]        
        self.glo_init_properties(**args)
    
    def glo_install_properties(self):
        GLAtomList.glo_install_properties(self)

        self.glo_add_property(
            { "name":        "COR",
              "desc":        "TLS Center of Reaction", 
              "catagory":    "TLS Analysis",
              "type":        "array(3)",
              "default":     zeros(3),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_animation_visible",
              "desc":        "Show L1 Animation",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     True,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L2_animation_visible",
              "desc":        "Show L2 Animation",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     True,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L3_animation_visible",
              "desc":        "Show L3 Animation",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     True,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L1_eigen_vec",
              "desc":        "L1 Eigen Vector", 
              "type":        "array(3)",
              "default":     zeros(3),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_eigen_vec",
              "desc":        "L2 Eigen Vector", 
              "catagory":    "TLS Analysis",
              "type":        "array(3)",
              "default":     zeros(3),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_eigen_vec",
              "desc":        "L3 Eigen Vector", 
              "catagory":    "TLS Analysis",
              "type":        "array(3)",
              "default":     zeros(3),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_eigen_val",
              "desc":        "L1 Eigen Value", 
              "catagory":    "TLS Analysis",
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_eigen_val",
              "desc":        "L2 Eigen Value", 
              "catagory":    "TLS Analysis",
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_eigen_val",
              "desc":        "L3 Eigen Value", 
              "catagory":    "TLS Analysis",
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_rho",
              "desc":        "L1 translation from COR", 
              "catagory":    "TLS Analysis",
              "type":        "array(3)",
              "default":     zeros(3),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_rho",
              "desc":        "L2 translation from COR", 
              "catagory":    "TLS Analysis",
              "type":        "array(3)",
              "default":     zeros(3),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_rho",
              "desc":        "L3 translation from COR", 
              "catagory":    "TLS Analysis",
              "type":        "array(3)",
              "default":     zeros(3),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_pitch",
              "desc":        "L1 screw pitch (A/DEG)", 
              "catagory":    "TLS Analysis",
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_pitch",
              "desc":        "L2 screw pitch (A/DEG)", 
              "catagory":    "TLS Analysis",
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_pitch",
              "desc":        "L3 screw pitch (A/DEG)", 
              "catagory":    "TLS Analysis",
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_rot",
              "desc":        "L1 Rotation", 
              "catagory":    "Simulation State",
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L2_rot",
              "desc":        "L2 Rotation", 
              "catagory":    "Simulation State",
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L3_rot",
              "desc":        "L3 Rotation", 
              "catagory":    "Simulation State",
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })

    def gldl_update_cb(self, updates, actions):
        GLAtomList.gldl_update_cb(self, updates, actions)
        if "trace" in updates:
            self.properties.update(lines=False)

    def gldl_iter_multidraw_self(self):
        """
        """
        for Lx_axis, Lx_rho, Lx_pitch, Lx_rot in (
            ("L1_eigen_vec", "L1_rho", "L1_pitch", "L1_rot"),
            ("L2_eigen_vec", "L2_rho", "L2_pitch", "L2_rot"),
            ("L3_eigen_vec", "L3_rho", "L3_pitch", "L3_rot") ):

            if Lx_axis=="L1_eigen_vec" and \
               self.properties["L1_animation_visible"]==False:
                continue
            if Lx_axis=="L2_eigen_vec" and \
               self.properties["L2_animation_visible"]==False:
                continue
            if Lx_axis=="L3_eigen_vec" and \
               self.properties["L3_animation_visible"]==False:
                continue

            Lx_axis  = self.properties[Lx_axis]
            Lx_rho   = self.properties[Lx_rho]
            Lx_pitch = self.properties[Lx_pitch]
            Lx_rot   = self.properties[Lx_rot]

            screw = Lx_axis * Lx_rot * Lx_pitch

            ## positive rotation
            glPushMatrix()
            glTranslatef(*Lx_rho + screw)
            glRotatef(Lx_rot, *Lx_axis)
            glTranslatef(*-Lx_rho)

            for draw_flag in GLAtomList.gldl_iter_multidraw_self(self):
                yield True

            glPopMatrix()

            ## negitive rotation
            glPushMatrix()
            glTranslatef(*Lx_rho - screw)
            glRotatef(-Lx_rot, *Lx_axis)
            glTranslatef(*-Lx_rho)

            for draw_flag in GLAtomList.gldl_iter_multidraw_self(self):
                yield True

            glPopMatrix()

    def gldl_iter_multidraw_self_old_slow(self):
        """
        """        
        for Lx_axis, Lx_rho, Lx_pitch, Lx_rot in (
            ("L1_eigen_vec", "L1_rho", "L1_pitch", "L1_rot"),
            ("L2_eigen_vec", "L2_rho", "L2_pitch", "L2_rot"),
            ("L3_eigen_vec", "L3_rho", "L3_pitch", "L3_rot") ):

            if Lx_axis=="L1_eigen_vec" and \
               self.properties["L1_animation_visible"]==False:
                continue
            if Lx_axis=="L2_eigen_vec" and \
               self.properties["L2_animation_visible"]==False:
                continue
            if Lx_axis=="L3_eigen_vec" and \
               self.properties["L3_animation_visible"]==False:
                continue
            
            Lx_axis  = self.properties[Lx_axis]
            Lx_rho   = self.properties[Lx_rho]
            Lx_pitch = self.properties[Lx_pitch]
            Lx_rot   = self.properties[Lx_rot]

            rot = 0.0
            
            while 1:
                screw = Lx_axis * rot * Lx_pitch

                ## positive rotation
                glPushMatrix()
                glTranslatef(*Lx_rho + screw)
                glRotatef(rot, *Lx_axis)
                glTranslatef(*-Lx_rho)

                for draw_flag in GLAtomList.gldl_iter_multidraw_self(self):
                    yield True

                glPopMatrix()

                ## negitive rotation
                glPushMatrix()
                glTranslatef(*Lx_rho - screw)
                glRotatef(-rot, *Lx_axis)
                glTranslatef(*-Lx_rho)

                for draw_flag in GLAtomList.gldl_iter_multidraw_self(self):
                    yield True

                glPopMatrix()

                if rot >= Lx_rot:
                    break

                rot += 0.50
                if rot > Lx_rot:
                    rot = Lx_rot


class GLTLSGroup(GLDrawList):
    """Top level visualization object for a TLS group.
    """
    def __init__(self, **args):
        GLDrawList.__init__(self)

        ## TLS calculations

        ## step 1: copy the TLS group
        orig_tls_group        = args["tls_group"]
        calcs                 = orig_tls_group.calc_COR()

        ## create this object's TLS group at the center of reaction
        self.tls_group        = TLSGroup(orig_tls_group)      
        self.tls_group.origin = calcs["COR"].copy()
        self.tls_group.T      = calcs["T'"].copy()
        self.tls_group.L      = calcs["L'"].copy()
        self.tls_group.S      = calcs["S'"].copy()

        L_eigen_vals, L_eigen_vecs = eigenvectors(self.tls_group.L)

        ## step 3: add visualization objects for each TLS screw rotation axis
        self.gl_screw_rot1 = GLTLSScrewRotation(atom_list = self.tls_group)
        self.gl_screw_rot2 = GLTLSScrewRotation(atom_list = self.tls_group)
        self.gl_screw_rot3 = GLTLSScrewRotation(atom_list = self.tls_group)
        
        self.gl_screw_rot1.glo_set_properties_id("gl_screw_rot1")
        self.gl_screw_rot2.glo_set_properties_id("gl_screw_rot2")
        self.gl_screw_rot3.glo_set_properties_id("gl_screw_rot3")
        
        self.glo_add_child(self.gl_screw_rot1)
        self.glo_add_child(self.gl_screw_rot2)
        self.glo_add_child(self.gl_screw_rot3)

        self.glo_link_child_property(
            "screw1_visible", "gl_screw_rot1", "visible")
        self.glo_link_child_property(
            "screw2_visible", "gl_screw_rot2", "visible")
        self.glo_link_child_property(
            "screw3_visible", "gl_screw_rot3", "visible")

        self.glo_link_child_property(
            "symmetry", "gl_screw_rot1", "symmetry")
        self.glo_link_child_property(
            "symmetry", "gl_screw_rot2", "symmetry")
        self.glo_link_child_property(
            "symmetry", "gl_screw_rot3", "symmetry")

        self.glo_link_child_property(
            "COR", "gl_screw_rot1", "COR")
        self.glo_link_child_property(
            "COR", "gl_screw_rot2", "COR")
        self.glo_link_child_property(
            "COR", "gl_screw_rot3", "COR")

        self.glo_link_child_property(
            "L1_eigen_vec", "gl_screw_rot1", "Lx_eigen_vec")
        self.glo_link_child_property(
            "L2_eigen_vec", "gl_screw_rot2", "Lx_eigen_vec")
        self.glo_link_child_property(
            "L3_eigen_vec", "gl_screw_rot3", "Lx_eigen_vec")
        self.glo_link_child_property(
            "L1_eigen_val", "gl_screw_rot1", "Lx_eigen_val")
        self.glo_link_child_property(
            "L2_eigen_val", "gl_screw_rot2", "Lx_eigen_val")
        self.glo_link_child_property(
            "L3_eigen_val", "gl_screw_rot3", "Lx_eigen_val")
        self.glo_link_child_property(
            "L1_rho", "gl_screw_rot1", "Lx_rho")
        self.glo_link_child_property(
            "L2_rho", "gl_screw_rot2", "Lx_rho")
        self.glo_link_child_property(
            "L3_rho", "gl_screw_rot3", "Lx_rho")
        self.glo_link_child_property(
            "L1_pitch", "gl_screw_rot1", "Lx_pitch")
        self.glo_link_child_property(
            "L2_pitch", "gl_screw_rot2", "Lx_pitch")
        self.glo_link_child_property(
            "L3_pitch", "gl_screw_rot3", "Lx_pitch")
        self.glo_link_child_property(
            "L1_rot", "gl_screw_rot1", "Lx_rot")
        self.glo_link_child_property(
            "L2_rot", "gl_screw_rot2", "Lx_rot")
        self.glo_link_child_property(
            "L3_rot", "gl_screw_rot3", "Lx_rot")

        ## step 4: add child GLTLSAtomList 
        self.gl_atom_list = GLTLSAtomList(
            atom_list   = self.tls_group,
            origin      = calcs["COR"],
            atom_origin = calcs["COR"])

        self.gl_atom_list.glo_set_properties_id("gl_atom_list")
        self.glo_add_child(self.gl_atom_list)
        
        self.glo_link_child_property(
            "animation_visible", "gl_atom_list", "visible")
        self.glo_link_child_property(
            "COR", "gl_atom_list", "atom_list")
        self.glo_link_child_property(
            "COR", "gl_atom_list", "COR")
        self.glo_link_child_property(
            "L1_eigen_vec", "gl_atom_list", "L1_eigen_vec")
        self.glo_link_child_property(
            "L2_eigen_vec", "gl_atom_list", "L2_eigen_vec")
        self.glo_link_child_property(
            "L3_eigen_vec", "gl_atom_list", "L3_eigen_vec")
        self.glo_link_child_property(
            "L1_eigen_val", "gl_atom_list", "L1_eigen_val")
        self.glo_link_child_property(
            "L2_eigen_val", "gl_atom_list", "L2_eigen_val")
        self.glo_link_child_property(
            "L3_eigen_val", "gl_atom_list", "L3_eigen_val")
        self.glo_link_child_property(
            "L1_rho", "gl_atom_list", "L1_rho")
        self.glo_link_child_property(
            "L2_rho", "gl_atom_list", "L2_rho")
        self.glo_link_child_property(
            "L3_rho", "gl_atom_list", "L3_rho")
        self.glo_link_child_property(
            "L1_pitch", "gl_atom_list", "L1_pitch")
        self.glo_link_child_property(
            "L2_pitch", "gl_atom_list", "L2_pitch")
        self.glo_link_child_property(
            "L3_pitch", "gl_atom_list", "L3_pitch")
        self.glo_link_child_property(
            "L1_rot", "gl_atom_list", "L1_rot")
        self.glo_link_child_property(
            "L2_rot", "gl_atom_list", "L2_rot")
        self.glo_link_child_property(
            "L3_rot", "gl_atom_list", "L3_rot")
        self.glo_link_child_property(
            "atom_color", "gl_atom_list", "color")
        self.glo_link_child_property(
            "atom_line_width", "gl_atom_list", "line_width")
        self.glo_link_child_property(
            "symmetry", "gl_atom_list", "symmetry")
        self.glo_link_child_property(
            "trace", "gl_atom_list", "trace")
 
        ## initalize properties
        self.glo_add_update_callback(self.tls_update_cb)
        self.glo_init_properties(
            COR          = calcs["COR"],
            T            = self.tls_group.T,
            L            = self.tls_group.L * rad2deg2,
            S            = self.tls_group.S * rad2deg,
            L1_eigen_vec = L_eigen_vecs[0],
            L2_eigen_vec = L_eigen_vecs[1],
            L3_eigen_vec = L_eigen_vecs[2],
            L1_eigen_val = L_eigen_vals[0] * rad2deg2,
            L2_eigen_val = L_eigen_vals[1] * rad2deg2,
            L3_eigen_val = L_eigen_vals[2] * rad2deg2,
            L1_rho       = calcs["L1_rho"],
            L2_rho       = calcs["L2_rho"],
            L3_rho       = calcs["L3_rho"],
            L1_pitch     = calcs["L1_pitch"] * (1.0/rad2deg),
            L2_pitch     = calcs["L2_pitch"] * (1.0/rad2deg),
            L3_pitch     = calcs["L3_pitch"] * (1.0/rad2deg),
            **args)

    def glo_install_properties(self):
        GLDrawList.glo_install_properties(self)

        self.glo_add_property(
            { "name":        "symmetry",
              "desc":        "Show Symmetry Equivelant",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "COR",
              "desc":        "TLS Center of Reaction",
              "catagory":    "TLS Analysis",
              "type":        "array(3)",
              "default":     zeros(3),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "T",
              "desc":        "Translation Tensor (A*A)",
              "catagory":    "TLS Analysis",
              "type":        "array(3,3)",
              "default":     zeros((3,3)),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L",
              "desc":        "Libration Tensor (DEG*DEG)",
              "catagory":    "TLS Analysis",
              "type":        "array(3,3)",
              "default":     zeros((3,3)),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "S",
              "desc":        "Skew Tensor (A*DEG)",
              "catagory":    "TLS Analysis",
              "type":        "array(3,3)",
              "default":     zeros((3,3)),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_eigen_vec",
              "desc":        "L1 Eigen Vector",
              "catagory":    "TLS Analysis",
              "type":        "array(3)",
              "default":     zeros(3),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_eigen_vec",
              "desc":        "L2 Eigen Vector",
              "catagory":    "TLS Analysis",
              "type":        "array(3)",
              "default":     zeros(3),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_eigen_vec",
              "desc":        "L3 Eigen Vector", 
              "catagory":    "TLS Analysis",
              "type":        "array(3)",
              "default":     zeros(3),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_eigen_val",
              "desc":        "L1 Eigen Value", 
              "catagory":    "TLS Analysis",
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_eigen_val",
              "desc":        "L2 Eigen Value", 
              "catagory":    "TLS Analysis",
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_eigen_val",
              "desc":        "L3 Eigen Value", 
              "catagory":    "TLS Analysis",
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_rho",
              "desc":        "L1 translation from COR", 
              "catagory":    "TLS Analysis",
              "type":        "array(3)",
              "default":     zeros(3),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_rho",
              "desc":        "L2 translation from COR", 
              "catagory":    "TLS Analysis",
              "type":        "array(3)",
              "default":     zeros(3),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_rho",
              "desc":        "L3 translation from COR", 
              "catagory":    "TLS Analysis",
              "type":        "array(3)",
              "default":     zeros(3),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_pitch",
              "desc":        "L1 screw pitch (A/DEG)", 
              "catagory":    "TLS Analysis",
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_pitch",
              "desc":        "L2 screw pitch (A/DEG)", 
              "catagory":    "TLS Analysis",
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_pitch",
              "desc":        "L3 screw pitch (A/DEG)", 
              "catagory":    "TLS Analysis",
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_rot",
              "desc":        "L1 Rotation",
              "catagory":    "Simulation State",
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L2_rot",
              "desc":        "L2 Rotation", 
              "catagory":    "Simulation State",
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L3_rot",
              "desc":        "L3 Rotation",
              "catagory":    "Simulation State",
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })        
        self.glo_add_property(
            { "name":        "time",
              "desc":        "Simulation Time",
              "catagory":    "Simulation State",
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })

        ## visualization properties
        self.glo_add_property(
            { "name":        "TLS_visible",
              "desc":        "Show TLS Tensors",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     True,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":       "tensor_line_width",
              "desc":       "Line Width of Tensors",
              "catagory":   "Tensors",
              "type":       "float",
              "default":    1.0,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "tensor_axis_radius",
              "desc":       "Radius of Tensor Axes",
              "catagory":   "Tensors",
              "type":       "float",
              "default":    0.05,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":        "scale",
              "desc":        "Tensor Display Scaling Factor",
              "catagory":    "Tensors",
              "type":        "float",
              "default":     1.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "T_color",
              "desc":        "T Tensor Color",
              "catagory":    "Colors",
              "type":        "color",
              "default":     (0.0, 1.0, 0.0),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":       "L_color",
              "desc":       "L Tensor Color",
              "catagory":   "Colors",
              "type":       "color",
              "default":    (1.0, 0.0, 0.0),
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "S_color",
              "desc":       "S Tensor Color",
              "catagory":   "Colors",
              "type":       "color",
              "default":    (0.0, 0.0, 1.0),
              "action":     "recompile" })
        self.glo_add_property(
            { "name":        "CA_line_visible",
              "desc":        "Show Lines to C-Alpha Atoms",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "fan_visible",
              "desc":        "Show Fans to Backbone",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "fan_opacity",
              "desc":        "Set Backbone Fan Transparancy",
              "catagroy":    "Fans",
              "type":        "float",
              "default":     0.60,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "fan_material",
              "desc":        "Set Backbone Fan Material",
              "catagory":    "Colors",
              "type":        "material",
              "default":     "emerald",
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "atom_color",
              "desc":        "Atom Color",
              "catagory":    "Colors",
              "type":        "color",
              "default":     (1.0, 1.0, 1.0),
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "atom_line_width",
              "desc":        "Atom Line Width",
              "type":        "float",
              "default":     1.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "U",
              "desc":        "Show ADP Utls Axes",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "U_color",
              "desc":        "ADP Axes Color",
              "catagory":    "Colors",
              "type":        "color",
              "default":     (0.0, 1.0, 0.0),
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "trace",
              "desc":        "Protein Backbone Trace Only", 
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "animation_visible",
              "desc":        "Show Animated Atoms", 
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "screw1_visible",
              "desc":        "Show Screw Axis 1 Gaussian Surface", 
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "screw2_visible",
              "desc":        "Show Screw Axis 2 Gaussian Surface", 
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "screw3_visible",
              "desc":        "Show Screw Axis 3 Gaussian Surface",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "redraw" })
         
    def tls_update_cb(self, updates, actions):
        if "time" in updates:
            self.update_time()

    def update_time(self):
        """Changes the time of the TLS group simulating harmonic motion.
        """
        sin_tm = math.sin(3.0 * self.properties["time"] * 2 * math.pi)
        sin_tm = abs(sin_tm)

        ## L Tensor
        L1_peak = math.sqrt(abs(self.properties["L1_eigen_val"]))
        L2_peak = math.sqrt(abs(self.properties["L2_eigen_val"]))
        L3_peak = math.sqrt(abs(self.properties["L3_eigen_val"]))

        L1_rot  = abs(L1_peak * sin_tm)
        L2_rot  = abs(L2_peak * sin_tm)
        L3_rot  = abs(L3_peak * sin_tm)

        self.glo_update_properties(L1_rot=L1_rot, L2_rot=L2_rot, L3_rot=L3_rot)

    def gldl_iter_multidraw_self(self):
        """Specialized draw list invokation to recycle the draw list for
        symmetry related copies.  Cartesian versions of the symmetry rotation
        and translation operators are generated by GLStructure/UnitCell
        classes.
        """
        if self.properties["symmetry"]==False:
            yield True
            
        else:

            gl_struct = self.glo_get_glstructure()
            if gl_struct==None:
                yield True

            else:
                for symop in gl_struct.iter_orth_symops():
                    glPushMatrix()

                    glMultMatrixf(
                        (symop.R[0,0], symop.R[1,0], symop.R[2,0], 0.0,
                         symop.R[0,1], symop.R[1,1], symop.R[2,1], 0.0,
                         symop.R[0,2], symop.R[1,2], symop.R[2,2], 0.0,
                         symop.t[0],   symop.t[1],   symop.t[2],   1.0) )

                    yield True
                    glPopMatrix()

    def gldl_draw(self):
        ## everything is drawn from the center of reaction
        glPushMatrix()
        glTranslatef(*self.properties["COR"])

        ## draw TLS axes
        if self.properties["TLS_visible"]==True:
            self.draw_tensors()

        ## draw Utls thermal axes
        if self.properties["U"]==True:
            self.draw_Utls()

        ## draw a line from the COR (center of reaction)
        ## to all CA atoms in the TLS group
        if self.properties["CA_line_visible"]==True:
            self.draw_CA_lines()

        ## draw a transparent fan from the COR to all backbone atoms
        ## in the TLS group
        if self.properties["fan_visible"]==True:
            self.draw_fan()

        glPopMatrix()

    def draw_Utls(self):
        """Render the anisotropic thremal axes calculated from the TLS
        model.
        """
        COR   = self.properties["COR"]
        color = self.properties["U_color"]

        for atm, Utls in self.tls_group.iter_atm_Utls():
            self.glr_U_axes(
                atm.position - COR,
                Utls,
                color,
                1.0)
        
    def draw_CA_lines(self):
        glDisable(GL_LIGHTING)
        glColor3f(0.5, 0.5, 0.5)
        glLineWidth(1.0)
        
        for atm in self.tls_group:
            if atm.name in ["CA"]:
                glBegin(GL_LINES)
                glVertex3f(0.0, 0.0, 0.0)
                glVertex3f(*atm.position - self.properties["COR"])
                glEnd()

        glEnable(GL_LIGHTING)

    def draw_fan(self):
        glEnable(GL_LIGHTING)

        self.glr_set_material_name(
            self.properties["fan_material"],
            self.properties["fan_opacity"])
        
        glBegin(GL_TRIANGLE_FAN)
        glVertex3f(0.0, 0.0, 0.0)
        for atm in self.tls_group:
            if atm.name in ["N", "CA", "C"]:
                glVertex3f(*atm.position - self.properties["COR"])
        glEnd()

    def draw_tensors(self):
        """Draw tensor axis.
        """
        scale = self.properties["scale"]

        glDisable(GL_LIGHTING)
        glLineWidth(self.properties["tensor_line_width"])

        ## T: units (A^2)
        glColor3f(*self.properties["T_color"])
        (eigen_values, eigen_vectors) = eigenvectors(self.properties["T"])
        
        for i in range(3):
            v  = array(eigen_vectors[i])
            v *= math.sqrt(abs(eigen_values[i]))
            v *= self.properties["scale"]
            v *= 0.5
            
            glBegin(GL_LINES)
            glVertex3f(*-v)
            glVertex3f(*v)
            glEnd()

            glEnable(GL_LIGHTING)
            r, g, b = self.properties["T_color"]
            self.glr_set_material_rgb(r, g, b, 1.0)
            self.glr_axis(-v, 2*v, self.properties["tensor_axis_radius"])
            glDisable(GL_LIGHTING)

        ## L: units (RAD^2)
        glColor3f(*self.properties["L_color"])
        (eigen_values, eigen_vectors) = eigenvectors(self.tls_group.L)
        
        for Lx_eigen_val, Lx_eigen_vec, Lx_rho in [
            ("L1_eigen_val", "L1_eigen_vec", "L1_rho"),
            ("L2_eigen_val", "L2_eigen_vec", "L2_rho"),
            ("L3_eigen_val", "L3_eigen_vec", "L3_rho")]:

            v  = self.properties[Lx_eigen_vec]
            v *= math.sqrt(abs(self.properties[Lx_eigen_val]))
            v *= self.properties["scale"]
            v *= 0.5

            rho = self.properties[Lx_rho]

            glBegin(GL_LINES)

            glVertex3f(0.0, 0.0, 0.0)
            glVertex3f(*rho)
            
            glVertex3f(*-v + rho)
            glVertex3f(* v + rho)

            glEnd()

            glEnable(GL_LIGHTING)
            r, g, b = self.properties["L_color"]
            self.glr_set_material_rgb(r, g, b, 1.0)
            self.glr_axis(rho-v, 2*v, self.properties["tensor_axis_radius"])
            glDisable(GL_LIGHTING)


class GLFragment(GLDrawList):
    """Not implemented yet.
    """
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


class GLChain(GLDrawList):
    """Visualization object for mmLib.Structure.Chain.
    """
    def __init__(self, **args):
        GLDrawList.__init__(self, **args)
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
        GLDrawList.glo_install_properties(self)

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
               "catagory":  "Colors",
               "type":      "color",
               "default":   None,
               "action":    "redraw" })
        self.glo_add_property(
            { "name":      "U",
              "desc":      "Show U Axes",
              "type":      "boolean",
              "default":   False,
              "action":    "redraw" })

    def glo_name(self):
        return "Chain %s" % (self.chain.chain_id)


class GLStructure(GLDrawList):
    """Visualization object for a mmLib.Structure.Structure.
    """
    def __init__(self, **args):
        GLDrawList.__init__(self, **args)
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
                "U", str(chain), "U")
            self.glo_link_child_property(
                "color", str(chain), "color")

        ## init properties
        self.glo_init_properties(**args)

    def glo_install_properties(self):
        GLDrawList.glo_install_properties(self)

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
              "desc":     "Show Amino Acid Main Chains",
              "type":     "boolean",
              "default":  True,
              "action":   "redraw" })
        self.glo_add_property(
            { "name":     "aa_side_chain_visible",
              "desc":     "Show Amino Acid Side Chains",
              "type":     "boolean",
              "default":  True,
              "action":   "redraw" })
        self.glo_add_property(
            { "name":              "dna_main_chain_visible",
              "desc":              "Show DNA Main Chains",
              "type":              "boolean",
              "default":           True,
              "action":            "redraw" })
        self.glo_add_property(
            { "name":              "dna_side_chain_visible",
              "desc":              "Show DNA Side Chains",
              "type":              "boolean",
              "default":           True,
              "action":            "redraw" })
        self.glo_add_property(
            { "name":              "hetatm_visible",
              "desc":              "Show HET (Non-Standard) Groups",
              "type":              "boolean",
              "default":           True,
              "action":            "redraw" })
        self.glo_add_property(
            { "name":              "water_visible",
              "desc":              "Show Waters",
              "type":              "boolean",
              "default":           True,
              "action":            "redraw" })
        self.glo_add_property(
            { "name":              "color",
              "desc":              "Set Solid Color for Structure",
              "catagory":          "Colors",
              "type":              "color",
              "default":           None,
              "action":            "redraw" })
        self.glo_add_property(
            { "name":              "U",
              "desc":              "Show ADP Axes",
              "type":              "boolean",
              "default":           False,
              "action":            "redraw" })

    def glo_name(self):
        return "%s" % (self.struct.cifdb.get_entry_id())

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
        self.glo_add_update_callback(self.glv_update_cb)
        self.glo_init_properties()

    def glo_install_properties(self):
        GLObject.glo_install_properties(self)

        self.glo_add_property(
            { "name":     "R",
              "desc":     "View Window Rotation Matrix",
              "type":     "array(3,3)",
              "default":  identity(3),
              "action":   "redraw" })
        self.glo_add_property(
            { "name":     "t",
              "desc":     "View Window Translation Vector",
              "type":     "array(3)",
              "default":  array([0.0, 0.0, -50]),
              "action":   "redraw" })

    def glo_name(self):
        return "GLViewer"

    def glv_update_cb(self, updates, actions):
        if "redraw" in actions:
            self.glv_redraw()
            
    def glv_add_draw_list(self, draw_list):
        """Append a GLDrawList.
        """
        assert isinstance(draw_list, GLDrawList)
        self.glo_add_child(draw_list)
        self.glv_redraw()

    def glv_remove_draw_list(self, draw_list):
        """Remove a GLDrawList.
        """
        assert isinstance(draw_list, GLDrawList)

        ## delete the compiled GL draw lists
        draw_list.gldl_delete_list()

        ## remove and trigger redraw
        draw_list.glo_remove()
        self.glv_redraw()

    def glv_add_struct(self, struct):
        """Adds the visualization for a mmLib.Structure.Structure object
        to the GLViewer.  It returns the GLStructure object created to
        visualize the Structure object.
        """
        assert isinstance(struct, Structure)
        
        gl_struct = GLStructure(struct=struct)
        self.glv_add_draw_list(gl_struct)
        return gl_struct

    def glv_redraw(self):
        """This method is called by GLViewer children to trigger a redraw
        in the toolkit embedding the GLViewer object.  It needs to be
        re-implemented when subclassed to call the tookit's widget redraw
        method.
        """
        pass

    def glv_init(self):
        """Called once to initalize the GL scene before drawing.
        """
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
        #glEnable(GL_POLYGON_SMOOTH)

        ## ALPHA BLENDING
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

    def glv_resize(self, width, height):
        """Called to set the size of the OpenGL window this class is
        drawing on.
        """
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

    def glv_translate(self, x, y, z):
        """Translate the current view by vector (x,y,z).
        """
        R     = self.properties["R"]
        t     = self.properties["t"]
        delta = array([x, y, z])
        
        upt = matrixmultiply(inverse(R), delta) + t
        self.properties.update(t=upt)

    def glv_rotate(self, alpha, beta, gamma):
        """Change the viewing position of the structure.  Changes to the
        current viewport are given relative to the current view vector.
        """
        R = self.properties["R"]

        Rx = rmatrixu(
            matrixmultiply(inverse(R), array([1.0, 0.0, 0.0])),
            math.radians(alpha))

        Ry = rmatrixu(
            matrixmultiply(inverse(R), array([0.0, 1.0, 0.0])),
            math.radians(beta))

        Rz = rmatrixu(
            matrixmultiply(inverse(R), array([0.0, 0.0, 1.0])),
            math.radians(gamma))

        Rxyz = matrixmultiply(Rz, matrixmultiply(Ry, Rx))
        upR  = matrixmultiply(R, Rxyz)

        self.properties.update(R=upR)
        
    def glv_render(self):
        """Draw all GLDrawList objects onto the given glcontext/gldrawable.
        If the GLDrawList objects are not yet compiled into OpenGL draw
        lists, they will be compiled while they are drawn, since this is
        a useful optimization.
        """
        R = self.properties["R"]
        t = self.properties["t"]

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
	glLoadIdentity()

        glTranslatef(0.0, 0.0, -10.0)

        glMultMatrixf(
            (R[0,0],           R[1,0],      R[2,0], 0.0,
             R[0,1],           R[1,1],      R[2,1], 0.0,
             R[0,2],           R[1,2],      R[2,2], 0.0,
                0.0,              0.0,         0.0, 1.0) )

        glTranslatef(t[0], t[1], t[2] + 10.0)

        for draw_list in self.glo_iter_children():
            draw_list.gldl_render()
