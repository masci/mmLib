## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""OpenGL rendering classes.
"""
from __future__  import generators

import copy

from OpenGL.GL      import *
from OpenGL.GLU     import *
from OpenGL.GLUT    import *
from mmTypes        import *
from Structure      import *
from Extensions.TLS import *
from GeometryDict   import *

try:
    import glaccel
except ImportError:
    GLACCEL_EXISTS = False
else:
    GLACCEL_EXISTS = True


## MISC Constents
PROP_OPACITY_RANGE    = "0.0-1.0,0.1"
PROP_LINE_RANGE       = "1.0-10.0,1.0"
PROP_PROBABILTY_RANGE = "1-99,1"


## trivariate gaussian critical values
GAUSS3C = {
    1:  0.3389,
    2:  0.4299,
    3:  0.4951,
    4:  0.5479,
    5:  0.5932,

    6:  0.6334,
    7:  0.6699,
    8:  0.7035,
    9:  0.7349,
    10: 0.7644,

    11: 0.7924,
    12: 0.8192,
    13: 0.8447,
    14: 0.8694,
    15: 0.8932,

    16: 0.9162,
    17: 0.9386,
    18: 0.9605,
    19: 0.9818,
    20: 1.0026,

    21: 1.0230,
    22: 1.0430,
    23: 1.0627,
    24: 1.0821,
    25: 1.1012,

    26: 1.1200,
    27: 1.1386,
    28: 1.1570,
    29: 1.1751,
    30: 1.1932,

    31: 1.2110,
    32: 1.2288,
    33: 1.2464,
    34: 1.2638,
    35: 1.2812,

    36: 1.2985,
    37: 1.3158,
    38: 1.3330,
    39: 1.3501,
    40: 1.3627,

    41: 1.3842,
    42: 1.4013,
    43: 1.4183,
    44: 1.4354,
    45: 1.4524,

    46: 1.4695,
    47: 1.4866,
    48: 1.5037,
    49: 1.5209,
    50: 1.8724,

    51: 1.5555,
    52: 1.5729,
    53: 1.5904,
    54: 1.6080,
    55: 1.6257,

    56: 1.6436,
    57: 1.6616,
    58: 1.6797,
    59: 1.6980,
    60: 1.7164,

    61: 1.7351,
    62: 1.7540,
    63: 1.7730,
    64: 1.7924,
    65: 1.8119,

    66: 1.8318,
    67: 1.8519,
    68: 1.8724,
    69: 1.8932,
    70: 1.9144,

    71: 1.9360,
    72: 1.9580,
    73: 1.9840,
    74: 2.0034,
    75: 2.0269,

    76: 2.0510,
    77: 2.0757,
    78: 2.1012,
    79: 2.1274,
    80: 2.1544,

    81: 2.1824,
    82: 2.2114,
    83: 2.2416,
    84: 2.2730,
    85: 2.3059,

    86: 2.3404,
    87: 2.3767,
    88: 2.4153,
    89: 2.4563,
    90: 2.5003,

    91: 2.5478,
    92: 2.5997,
    93: 2.6571,
    94: 2.7216,
    95: 2.7955,

    96: 2.8829,
    97: 2.9912,
    98: 3.1365,
    99: 3.3682 }


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
        self.__globject_properties_name      = None
        self.__globject_properties           = []
        self.__globject_properties_callbacks = []

        self.glo_install_properties()

    def glo_name(self):
        """Returns the GLObject name.
        """
        if self.__globject_properties_name!=None:
            return self.__globject_properties_name
        elif self.__globject_properties_id!=None:
            return "%s(%s)" % (self.__class__.__name__,
                               self.__globject_properties_id)
        else:
            return self.__class__.__name__

    def glo_set_name(self, name):
        """Sets the GLObject name.
        """
        self.__globject_properties_name = name

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
        """Removes the child GLObject.
        """
        assert isinstance(child, GLObject)
        assert child.__globject_parent==self
        child.__globject_parent = None
        self.__globject_children.remove(child)

    def glo_remove(self):
        """The GLObject removes itself from its parent.
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
        n = self.glo_get_degree()
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
        parent_list = []
        composite = self
        while composite.__parent:
            composite = composite.__parent
            parent_list.append(composite)
        return parent_list

    def glo_get_lowest_common_ancestor(self, gl_object):
        """Returns the lowest common ancesotry of self and argument
        composite.
        """
        assert isinstance(gl_object, GLObject)

        pl1 = self.glo_get_parent_list()
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
        """Adds a new property to the GLObject.  The prop_desc is a dictionary
        with attributes describing the property.  See comments in source code
        for a description of the key values for property descriptions.
        """
        assert prop_desc["name"] not in self.properties

        ## is the proprty marked read-only, this is only a hint to
        ## the user interface, not a actual read-only property
        prop_desc["read_only"] = prop_desc.get(
            "read_only", False)

        ## the property triggers update callbacks when the
        ## value is changed
        prop_desc["update_on_changed"] = prop_desc.get(
            "update_on_changed", True)

        ## the property triggers update callbacks when the
        ## property is set with comparison to the old value
        prop_desc["update_on_set"] = prop_desc.get(
            "update_on_set", False)

        ## the property triggers update callbacks on initalization
        prop_desc["update_on_init"] = prop_desc.get(
            "update_on_init", True)

        self.__globject_properties.append(prop_desc)

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
        if prop_desc==None:
            raise ValueError,\
                  "GLObject.glo_link_child_property(x, y, z) "\
                  "parent has no property: %s" % (name)

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
        updates = {}
        actions = []

        for prop_desc in self.__globject_properties:
            name = prop_desc["name"]

            ## set the property value
            try:
                self.properties[name] = args[name]
            except KeyError:
                self.properties[name] = prop_desc["default"]

            ## if the update callbacks are to be triggered on initialization
            if prop_desc["update_on_init"]==True:
                 updates[name] = self.properties[name]

                 ## add changes to actions list
                 ## case 1: action is a string
                 if type(prop_desc["action"])==StringType:
                     if prop_desc["action"] not in actions:
                         actions.append(prop_desc["action"])
                 ## case 2: action is a list of strings
                 elif type(prop_desc["action"])==ListType:
                     for prop_action in prop_desc["action"]:
                         if prop_action not in actions:
                             actions.append(prop_action)

            ## propagate linked values
            try:
                linked_props = prop_desc["link"]
            except KeyError:
                pass
            else:
                for linked_prop in linked_props:
                    child = self.glo_get_child(linked_prop["gl_object"])
                    child_name = linked_prop["name"]
                    child.glo_update_properties(
                        **{ child_name: self.properties[name] })
                    
        if len(updates)>0:
            for func in self.__globject_properties_callbacks:
                func(updates, actions)
                
    def glo_update_properties(self, **args):
        """Update property values and trigger update callbacks.
        """        
        updates = {}
        actions = []

        ## update properties
        for prop_desc in self.__globject_properties:
            name = prop_desc["name"]

            ## continue if this property is not being updated
            if args.has_key(name)==False:
                continue

            ## update_on_set:
            ##     If True, always trigger update callbacks when
            ##     the value is set.
            ## update_on_changed:
            ##     If true, trigger update callbacks when a value
            ##     is changed.
            if ( (prop_desc["update_on_set"]==True) or
                 (prop_desc["update_on_changed"]==True and
                  self.properties[name]!=args[name]) ):

                self.properties[name] = updates[name] = args[name]

                ## now update the actions taken when a property changes
                ## case 1: action is a string
                if type(prop_desc["action"])==StringType:
                    if prop_desc["action"] not in actions:
                        actions.append(prop_desc["action"])
                ## case 2: action is a list of strings
                elif type(prop_desc["action"])==ListType:
                    for prop_action in prop_desc["action"]:
                        if prop_action not in actions:
                            actions.append(prop_action)

            ## propagate updates for linked properties
            try:
                linked_props = prop_desc["link"]
            except KeyError:
                pass
            else:
                for linked_prop in linked_props:
                    child = self.glo_get_child(linked_prop["gl_object"])
                    child_name = linked_prop["name"]                    
                    child.glo_update_properties(
                        **{ child_name: self.properties[name] })

        if len(updates)>0:
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
    def glr_translate(self, t):
        """Translates the scene by vector t.
        """
        glTranslatef(*t)
    
    def glr_mult_matrix(self, R, t=zeros(3, Float)):
        """Multiplies the current matrix by rotation matrix R and translates
        by t
        """
        ## OpenGL wants the matrix in column-major form
        glMultMatrixf(
            (R[0,0], R[1,0], R[2,0], 0.0,
             R[0,1], R[1,1], R[2,1], 0.0,
             R[0,2], R[1,2], R[2,2], 0.0,
             t[0],   t[1],   t[2],   1.0) )
    
    def glr_set_material_rgb(self, r, g, b, a):
        """Creates a stock rendering material colored according to the given
        RGB values.
        """
        glColor3f(r, g, b)

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, (r, g, b, a))
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, (1.0, 1.0, 1.0, 1.0))
        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, (0.0, 0.0, 0.0, 1.0))

        if a<1.0:
            glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, 128.0)
        else:
            glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, 100.0)

    def glr_set_material_atom(self, atm, a):
        """Sets a material for rendering a atom which uses the atom's
        element color for the material selection.
        """
        elem = atm.get_structure().library.get_element(atm.element)
        if elem!=None:
            r, g, b = elem.color
        else:
            r, g, b = 1.0, 1.0, 1.0

        self.glr_set_material_rgb(r, g, b, a)

    def glr_text(self, text):
        """Renders a text string.
        """
        glDisable(GL_LIGHTING)
        glColor3f(1.0, 1.0, 1.0)
        glLineWidth(2.0)

        glPushMatrix()
        
        s = 1.0 / 1000.0
        glScalef(s, s, s)

        for c in text:
            glutStrokeCharacter(GLUT_STROKE_ROMAN, ord(c))

        glPopMatrix()
            
    def glr_axis(self, position, axis, radius):
        """Draw a vector axis using the current set material at position
        with the given radius.
        """
        ## don't bother redering small axes -- they look like junk
        if allclose(length(axis), 0.0):
            return
        
        end = position + axis

        try:
            glaccel.rod(
                position[0], position[1], position[2],
                end[0], end[1], end[2],
                radius)

        except NameError:
            glDisable(GL_LIGHTING)
            glLineWidth(1.0)
            glBegin(GL_LINES)
            glVertex3f(*position)
            glVertex3f(*end)
            glEnd()

    def glr_tube(self, pos1, pos2, radius):
        """Draws a hollow tube beginning at pos1, and ending at pos2.
        """
        try:
            glaccel.tube(
                pos1[0], pos1[1], pos1[2],
                pos2[0], pos2[1], pos2[2],
                radius)

        except NameError:
            glDisable(GL_LIGHTING)
            glLineWidth(1.0)
            glBegin(GL_LINES)
            glVertex3f(*pos1)
            glVertex3f(*pos2)
            glEnd()

    def glr_sphere(self, position, radius, quality):
        """Draw a atom as a CPK sphere.
        """
        try:
            glaccel.sphere(position[0], position[1], position[2],
                           radius, quality)

        except NameError:
            pass

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

    def glr_Uaxes(self, position, U, prob, color, line_width):
        """Draw the anisotropic axies of the atom at the given probability.
        """
        C = GAUSS3C[prob]        
        eval_U, evec_U = eigenvectors(U)

        try:
            v0_peak = C * math.sqrt(eval_U[0])
        except ValueError:
            v0_peak = 0.0
        try:
            v1_peak = C * math.sqrt(eval_U[1])
        except ValueError:
            v1_peak = 0.0
        try:
            v2_peak = C * math.sqrt(eval_U[2])
        except ValueError:
            v2_peak = 0.0
        
        v0 = evec_U[0] * v0_peak
        v1 = evec_U[1] * v1_peak
        v2 = evec_U[2] * v2_peak

        glDisable(GL_LIGHTING)
        glColor3f(*color)
        glLineWidth(line_width)

        glPushMatrix()
        glTranslatef(*position)

        glBegin(GL_LINES)
        glVertex3f(*-v0)
        glVertex3f(* v0)
        glVertex3f(*-v1)
        glVertex3f(* v1)
        glVertex3f(*-v2)
        glVertex3f(* v2)
        glEnd()
        
        glPopMatrix()

    def glr_Uellipse(self, position, U, prob):
        """Renders the ellipsoid enclosing the given fractional probability
        given the gaussian variance-covariance matrix U at the given position.
        C=1.8724 = 68%
        """
        try:
            glaccel.Uellipse(
                position[0], position[1], position[2],
                U[0,0], U[1,1], U[2,2], U[0,1], U[0,2], U[1,2],
                GAUSS3C[prob], 3)

        except NameError:
            pass

    def glr_Urms(self, position, U):
        """Renders the root mean square (one standard deviation) surface of the
        gaussian variance-covariance matrix U at the given position.  This
        is a peanut-shaped surface. (Note: reference the peanut paper!)
        """
        try:
            glaccel.Upeanut(
                position[0], position[1], position[2],
                U[0,0], U[1,1], U[2,2], U[0,1], U[0,2], U[1,2],
                3)

        except NameError:
            pass


class GLDrawList(GLObject, OpenGLRenderMethods):
    """Fundamental OpenGL rigid entity.
    """
    def __init__(self, **args):
        self.__draw_method_list = []
        
        GLObject.__init__(self, **args)
        self.glo_add_update_callback(self.gldl_update_cb)
        self.glo_init_properties(**args)

        self.gldl_install_draw_methods()

    def glo_remove_child(self, child):
        """Override GLObject's remove to also delete the compiled OpenGL
        draw lists.
        """
        self.gldl_delete_list()
        GLObject.glo_remove_child(self, child)

    def glo_remove(self):
        """Override GLObject's remove to also delete the compiled OpenGL
        draw lists.
        """
        self.gldl_delete_list()
        GLObject.glo_remove(self)

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
              "default":    identity(3, Float),
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
        self.glo_add_property(
            { "name" :      "visible_distance",
              "desc":       "",
              "type":       "float",
              "hidden":     True,
              "default":    0.0,
              "action":     "redraw" })

    def gldl_redraw(self):
        """Triggers a redraw of the GLViewer
        """
        gl_viewer = self.glo_get_root()
        if isinstance(gl_viewer, GLViewer):
            gl_viewer.glv_redraw()

    def gldl_get_glviewer(self):
        """Returns the root GLViewer object.
        """
        return self.glo_get_root()

    def gldl_install_draw_methods(self):
        """Override in children to install draw methods for a GLDrawList.
        """
        pass

    def gldl_update_cb(self, updates, actions):
        if "recompile" in actions:
            self.gldl_delete_list()
            self.gldl_redraw()

        elif "redraw" in actions:
            self.gldl_redraw()

        for draw_method in self.__draw_method_list:
            if draw_method["recompile_action"] in actions:
                self.gldl_draw_method_gl_draw_list_delete(draw_method)

                op = draw_method["opacity_property"]
                if op!=None:
                    draw_method["transparent"] = self.properties[op]<1.0

    def gldl_draw_method_install(self, draw_method):
        """Installs a draw method to compile and render a OpenGL draw listlist.
        keys:
           name:       text description of the method
           func:       the method to invoke to render the draw list
           tranparent: True if the draw list is drawing transparent

        private values:
           gl_draw_list_id: OpenGL Drawlist ID
        """
        assert draw_method.has_key("name")
        assert draw_method.has_key("func")

        draw_method["transparent"] = draw_method.get(
            "transparent", False)

        draw_method["no_gl_compile"] = draw_method.get(
            "no_gl_compile", False)

        draw_method["gl_draw_list_id"] = None

        draw_method["visible_property"] =  draw_method.get(
            "visible_property", None)

        draw_method["recompile_action"] = draw_method.get(
            "recompile_action", None)

        draw_method["opacity_property"] = draw_method.get(
            "opacity_property", None)

        draw_method["multidraw_iter"]  = draw_method.get(
            "multipdraw_iter", None)

        draw_method["multidraw_all_iter"] = draw_method.get(
            "multipdraw_all_iter", None)

        self.__draw_method_list.append(draw_method)

    def gldl_draw_method_get(self, draw_method_name):
        """Returns the draw metod of the given name or None if not found.
        """
        for draw_method in self.__draw_method_list:
            if draw_method["name"]==draw_method_name:
                return draw_method
        return None

    def gldl_draw_method_remove(self, draw_method_name):
        """Removes a draw method by name.
        """
        draw_method = self.gldl_draw_method_get(draw_method_name)
        if draw_method==None:
            raise ValueError, "GLDrawList.gl_draw_method_remove(x) x "\
                  "not in GLDrawList"

        self.__draw_method_list.remove(draw_method)
        if draw_method["gl_draw_list_id"]!=None:
            glDeleteLists(draw_method["gl_draw_list_id"], 1)
        del draw_method["gl_draw_list_id"]

    def gldl_draw_method_recompile(self, draw_method_name):
        """Deletes the current compiloed OpenGL draw list for the draw method
        """
        draw_method = self.gldl_draw_method_get(draw_method_name)
        if draw_method==None:
            raise ValueError,\
                  "GLDrawList.gl_draw_method_redraw(x) x not in GLDrawList"

        self.gldl_draw_method_gl_draw_list_delete(draw_method)

    def gldl_draw_method_recompile_all(self):
        """Deletes all currently compiled OpenGL compiled draw lists.
        """
        for draw_method in self.__draw_method_list:
            self.gldl_draw_method_gl_draw_list_delete(draw_method)

    def gldl_draw_method_set_transparent(self, draw_method_name, transparent):
        """Sets the transparent flag for a draw_method.
        """
        draw_method = self.gldl_draw_method_get(draw_method_name)
        if draw_method==None:
            raise ValueError,\
                  "GLDrawList.gl_draw_method_set_transparent(x, y) x "\
                  "not in GLDrawList"

        if draw_method["transparent"]!=transparent:
            self.gldl_draw_method_gl_draw_list_delete(draw_method)
            draw_method["transparent"] = transparent

    def gldl_draw_method_check_visible(self, draw_method):
        """Returns True if the draw method is visible, False if not.
        """
        vis_name = draw_method["visible_property"]
        if vis_name!=None:
            return self.properties[vis_name]
        return True

    def gldl_draw_method_gl_draw_list_compile(self, draw_method):
        """Performs a OpenGL Compile of the draw method if needed.
        """
        assert draw_method["gl_draw_list_id"]==None
        
        draw_method["gl_draw_list_id"] = glGenLists(1)
        glNewList(draw_method["gl_draw_list_id"], GL_COMPILE)
        draw_method["func"]()
        glEndList()

        self.gldl_redraw()

    def gldl_draw_method_gl_draw_list_delete(self, draw_method):
        """Deletes the compiled OpenGL draw list for the draw method.
        """
        if draw_method["gl_draw_list_id"]!=None:
            glDeleteLists(draw_method["gl_draw_list_id"], 1)
            draw_method["gl_draw_list_id"] = None

        self.gldl_redraw()

    def gldl_render_draw_methods(self, transparent):
        """Render all draw methods.
        """
        for draw_method in self.__draw_method_list:
            if self.gldl_draw_method_check_visible(draw_method)==False:
                continue

            ## transparent methods are only drawn when during the second
            ## rednering pass
            if draw_method["transparent"]!=transparent:
                continue

            ## some draw lists may be not be compiled into a OpenGL draw
            ## list, these have to be redrawn every time
            if draw_method["no_gl_compile"]==True:
                draw_method["func"]()

            else:
                if draw_method["gl_draw_list_id"]==None:
                    self.gldl_draw_method_gl_draw_list_compile(draw_method)
                glCallList(draw_method["gl_draw_list_id"])

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

    def gldl_render(self, transparent=False):
        """Compile or force a recompile of this object's gl_draw list, and
        render the scene.  Rendering the scene can be bypassed if
        this method is called with render = False.
        """
        if self.properties["visible"]==False:
            return

        self.gldl_push_matrix()

        ## support multiple rendering images by implementing class
        ## iterators gldl_iter_multidraw_all() for multiple
        ## rendering iterations of the GLDrawList and all its children,
        ## or gldl_iter_multidraw_self() for multiple images of just
        ## this GLDrawList, rendering the children just once
        for draw_flag_multi in self.gldl_iter_multidraw_all():

            for draw_flag_self in self.gldl_iter_multidraw_self():
                self.gldl_render_draw_methods(transparent)

            ## render first-level children of this GLDrawList
            ## which, in turn, will render their children
            for draw_list in self.glo_iter_children():
                if isinstance(draw_list, GLDrawList):
                    draw_list.gldl_render(transparent)

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
        ## this deletes all the compiled draw lists
        self.gldl_draw_method_recompile_all()
            
        for draw_list in self.glo_iter_children():
            if isinstance(draw_list, GLDrawList):
                draw_list.gldl_delete_list()

    def gldl_draw(self):
        """Implement in subclass to draw somthing.
        """
        pass

    def gldl_draw_transparent(self):
        """Implement in subclass to draw transparent objects.
        """
        pass


class GLAxes(GLDrawList):
    """Draw orthogonal axes in red = x, green = y, blue = z.
    """
    def __init__(self, **args):
        GLDrawList.__init__(self, **args)
        self.glo_set_name("Cartesian Axes")
        self.glo_init_properties(**args)

    def glo_install_properties(self):
        GLDrawList.glo_install_properties(self)

        self.glo_add_property(
            { "name":       "line_length",
              "desc":       "Axis Length",
              "type":       "float",
              "spin":       "1.0-500.0,10.0",
              "default":    20.0,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "line_width",
              "desc":       "Axis Radius",
              "type":       "float",
              "spin":      "0.0-5.0,0.1",
              "default":    0.1,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "color_x",
              "desc":       "X Axis Color",
              "catagory":   "Colors",
              "type":       "color",
              "default":    (1.0, 0.4, 0.4),
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "color_y",
              "desc":       "Y Axis Color",
              "catagory":   "Colors",
              "type":       "color",
              "default":    (0.4, 1.0, 0.4),
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "color_z",
              "desc":       "Z Axis Color",
              "catagory":   "Colors",
              "type":       "color",
              "default":    (0.4, 0.4, 1.0),
              "action":     "recompile" })

    def gldl_install_draw_methods(self):
        self.gldl_draw_method_install(
            { "name":        "axes",
              "func":        self.draw_axes,
              "transparent": False })

    def draw_axes(self):
        line_length = self.properties["line_length"]
        line_width  = self.properties["line_width"]
        
        (r, g, b) = self.properties["color_x"]
        self.glr_set_material_rgb(r, g, b, 1.0)
        self.glr_axis(zeros(3, Float),
                      array([line_length, 0.0, 0.0]), line_width) 

        (r, g, b) = self.properties["color_y"]
        self.glr_set_material_rgb(r, g, b, 1.0)
        self.glr_axis(zeros(3, Float),
                      array([0.0, line_length, 0.0]), line_width) 

        (r, g, b) = self.properties["color_z"]
        self.glr_set_material_rgb(r, g, b, 1.0)
        self.glr_axis(zeros(3, Float),
                      array([0.0, 0.0, line_length]), line_width) 


class GLUnitCell(GLDrawList):
    """Draw unit cell.
    """
    def __init__(self, **args):
        self.unit_cell = args["unit_cell"]

        GLDrawList.__init__(self, **args)
        self.glo_set_name("Unit Cell")        
        self.glo_init_properties(**args)

    def glo_install_properties(self):
        GLDrawList.glo_install_properties(self)
        
        self.glo_add_property(
            { "name":       "line_width",
              "desc":       "Line Width",
              "type":       "float",
              "spin":       PROP_LINE_RANGE,
              "default":    1.0,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "color",
              "desc":       "Line Color",
              "type":       "color",
              "default":    (1.0, 1.0, 1.0),
              "action":     "recompile" })

    def gldl_install_draw_methods(self):
        self.gldl_draw_method_install(
            { "name":        "unit_cell",
              "func":        self.draw_unit_cell,
              "transparent": False })

    def draw_unit_cell(self):
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


class GLAtomColor(GLObject):
    """This is going to replace the 3-tuple RGB color...
    """
    ## if any color functions fail, default to this
    ## color
    fail_color = (1.0, 1.0, 1.0)

    ## XXX: move into monomers.cif library!
    aa_type_color_dict = {
        "aliphatic":         (0.50, 0.50, 0.50),
        "aromatic":          (0.75, 0.75, 0.75),
        "sulfer-containing": (0.20, 1.00, 0.20),
        "alchols":           (1.00, 0.60, 0.60),
        "acids":             (1.00, 0.25, 0.25),
        "bases":             (0.25, 0.25, 1.00),
        "amides":            (0.60, 0.60, 1.00)}
        
    def __init__(self):
        GLObject.__init__(self)
        
        self.value      = None
        self.value_desc = None

        self.enum_list = [
            "Random"]

        self.enum_dict = {}

        self.install_color_method(
            "Color By Element", self.color_by_element)
        self.install_color_method(
            "Color By Residue Type", self.color_by_residue)
        self.install_color_method(
            "Color By Temp Factor", self.color_by_temp)
        self.install_color_method(
            "Color By Anisotropy", self.color_by_anisotropy)

    def install_color_method(self, name, func):
        """Adds a color function to the 
        """
        assert name not in self.enum_list
        self.enum_list.append(name)
        self.enum_dict[name.lower()] = func

    def get_enum_list(self):
        """Return a list of the values as a enumeration.
        """
        return self.enum_list

    def get_value_desc(self):
        """Returns the description string of the current enumerated value.
        """
        return self.value_desc

    def get_value(self):
        return self.value

    def set_value(self, value):
        """
        """
        if type(value)==StringType:
            value = value.lower()

            ## select a random color
            if value=="random":
                self.value = (random.random(),random.random(),random.random())
                self.value_desc = "Random RGB: %4.2f %4.2f %4.2f" % self.value
                    
            ## selects one of the named color functions
            elif value in self.enum_dict:
                self.value = value
                for upper in self.enum_list:
                    if upper.lower()==self.value:
                        self.value_desc = upper
                        break

            else:
                ## maybe the color is in R,G,B format
                try:
                    r, g, b = value.split(",")
                except ValueError:
                    pass
                else:
                    self.value = (float(r), float(g), float(b))
                    self.value_desc = "RGB: %4.2f %4.2f %4.2f" % self.value

        elif type(value)==TupleType:
            self.value = value
            self.value_desc = "RGB: %4.2f %4.2f %4.2f" % self.value

    def calc_atom_color(self, atom):
        """This is what this class is for!  Returns the r,g,b 3-tuple for
        the atom.
        """
        if type(self.value)==TupleType:
            return self.value
        else:
            return self.enum_dict[self.value](atom)

    def color_by_element(self, atom):
        """Returns the color based on the element type.  The color definitions
        for all the elements are in the mmLib/Data/elements.cif file.
        """
        element = atom.get_structure().library.get_element(atom.element)
        if element!=None:
            return element.color
        else:
            return self.fail_color

    def color_by_residue(self, atom):
        """Returns the color based on the residue type.
        """
        library = atom.get_structure().library
        monomer = library.get_monomer(atom.res_name)
        
        if monomer.is_amino_acid()==False:
            return self.fail_color

        try:
            return self.aa_type_color_dict[monomer.chem_type]
        except KeyError:
            return self.fail_color

    def color_by_temp(self, atom):
        """Returns the color based on the tempature factor.
        """
        if atom.temp_factor!=None:
            return self.color_by_range(10.0, 25.0, atom.temp_factor)
        else:
            return self.fail_color

    def color_by_anisotropy(self, atom):
        """Color by anisotropy.
        """
        aniso = atom.calc_anisotropy()
        return self.color_by_range(0.0, 1.0, 1.0 - aniso)

    def color_by_range(self, min, max, value):
        """Blue->Red
        """
        if value<min:
            b = 1.0
        if value>max:
            b = 0.0
        else:
            b = 1.0 - (value - min) / (max - min)

        return (1.0-b, 0.0, b)


class GLAtomList(GLDrawList):
    """OpenGL renderer for a list of atoms.  Optional arguments iare:
    color, U, U_color.
    """
    def __init__(self, **args):
        ## setup atom coloring object
        self.glal_atom_color         = GLAtomColor()
        self.glal_hidden_atoms_dict  = None
        self.glal_visible_atoms_dict = None
        self.glal_xyzdict            = None
        
        GLDrawList.__init__(self, **args)

        ## install atom coloring object
        self.glal_atom_color.glo_set_name("Atom Color Properties")
        self.glal_atom_color.glo_set_properties_id("glal_atom_color")
        self.glo_add_child(self.glal_atom_color)

        self.glo_add_update_callback(self.glal_update_properties)
        self.glo_init_properties(**args)

    def glo_install_properties(self):
        GLDrawList.glo_install_properties(self)

        ## Show/Hide
        self.glo_add_property(
            { "name":       "atom_origin",
              "desc":       "Atom Calculation Origin",
              "type":       "array(3)",
              "hidden":     True,
              "default":    None,
              "action":     ["recompile", "recalc_positions"] })

        self.glo_add_property(
            { "name":      "color_setting",
              "desc":      "Atom Color Setting",
              "catagory":  "Colors",
              "type":      "string",
              "read_only": True,
              "default":   "",
              "action":    "" })
        self.glo_add_property(
            { "name":      "color",
              "desc":      "Atom Color",
              "catagory":  "Colors",
              "type":      "enum_string",
              "default":   "Color By Element",
              "enum_list": self.glal_atom_color.get_enum_list(),
              "action":    "recompile" })

        self.glo_add_property(
            { "name":      "symmetry",
              "desc":      "Show Symmetry Equivelants",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   False,
              "action":    "redraw" })
        self.glo_add_property(
            { "name":      "main_chain_visible",
              "desc":      "Show Main Chain Atoms",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   True,
              "action":    ["recompile", "recalc_positions"] })
        self.glo_add_property(
            { "name":      "side_chain_visible",
              "desc":      "Show Side Chain Atoms",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   True,
              "action":    ["recompile", "recalc_positions"] })
        self.glo_add_property(
            { "name":      "hetatm_visible",
              "desc":      "Show Hetrogen Atoms",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   True,
              "action":    ["recompile", "recalc_positions"] })
        self.glo_add_property(
            { "name":      "water_visible",
              "desc":      "Show Waters",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   False,
              "action":    ["recompile", "recalc_positions"] })
        self.glo_add_property(
            { "name":      "hydrogen_visible",
              "desc":      "Show Hydrogens",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   False,
              "action":    ["recompile", "recalc_positions"] })

        ## labels
        self.glo_add_property(
            { "name":      "labels",
              "desc":      "Show Atom Lables",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   False,
              "action":    "recompile_labels" })

        self.glo_add_property(
            { "name":       "auto_label_distance",
              "desc":       "Label Auto-Show Distance",
              "catagory":   "Labels",
              "type":       "float",
              "default":    10.0,
              "action":     "redraw" })

        ## lines
        self.glo_add_property(
            { "name":       "lines",
              "desc":       "Draw Atom Bond Lines",
              "catagory":   "Show/Hide",
              "type":       "boolean",
              "default":    True,
              "action":     "recompile_lines" })
        self.glo_add_property(
            { "name":       "line_width",
              "desc":       "Bond Line Drawing Width",
              "catagory":   "Bond Lines",
              "type":       "float",
              "spin":       PROP_LINE_RANGE,
              "default":    1.0,
              "action":     "recompile_lines" })

        ## Ball/Stick
        self.glo_add_property(
            { "name":       "ball_stick",
              "desc":       "Draw Ball/Sticks",
              "catagory":   "Show/Hide",
              "type":       "boolean",
              "default":    False,
              "action":     "recompile_ball_stick" })
        self.glo_add_property(
            { "name":       "ball_radius",
              "desc":       "Atom (Ball) Radius",
              "catagory":   "Ball/Stick",
              "type":       "float",
              "default":    0.1,
              "action":     "recompile_ball_stick" })
        self.glo_add_property(
            { "name":       "stick_radius",
              "desc":       "Bond (Stick) Radius",
              "catagory":   "Ball/Stick",
              "type":       "float",
              "default":    0.1,
              "action":     "recompile_ball_stick" })

        ## cpk
        self.glo_add_property(
            { "name":       "cpk",
              "desc":       "Draw CPK Spheres",
              "catagory":   "Show/Hide",
              "type":       "boolean",
              "default":    False,
              "action":     "recompile_cpk" })
        self.glo_add_property(
            { "name":       "cpk_opacity_occupancy",
              "desc":       "Set Opacity by Atom Occupancy",
              "catagory":   "CPK",
              "type":       "boolean",
              "default":    False,
              "action":     "recompile_cpk" })
        self.glo_add_property(
            { "name":       "cpk_scale_radius",
              "desc":       "Scale CPK Radius",
              "catagory":   "CPK",
              "type":       "float",
              "spin":       "0.0-5.0,0.1",
              "default":    1.0,
              "action":     "recompile_cpk" })
        self.glo_add_property(
            { "name":       "cpk_opacity",
              "desc":       "CPK Sphere Opacity",
              "catagory":   "CPK",
              "type":       "float",
              "range":      PROP_OPACITY_RANGE,
              "default":    1.00,
              "action":     "recompile_cpk" })
        self.glo_add_property(
            { "name":       "sphere_quality",
              "desc":       "CPK Sphere Quality",
              "catagory":   "CPK",
              "type":       "integer",
              "range":      "5-36,5",
              "default":    10,
              "action":     "recompile_cpk" })

        ## trace           
        self.glo_add_property(
            { "name":       "trace",
              "desc":       "Draw Backbone Trace",
              "catagory":   "Show/Hide",
              "type":       "boolean",
              "default":    False,
              "action":     "recompile_trace" })
        self.glo_add_property(
            { "name":       "trace_radius",
              "desc":       "Trace Stick Radius",
              "catagory":   "Trace",
              "type":       "float",
              "default":    0.2,
              "action":     "recompile_trace" })
        self.glo_add_property(
            { "name":       "trace_color",
              "desc":       "Backbone Trace Color",
              "catagory":   "Trace",
              "type":       "color",
              "default":    (1.0, 1.0, 1.0),
              "action":     "recompile_trace" })

        ## ADPs
        self.glo_add_property(
            { "name":       "adp_prob",
              "desc":       "Contour Probability",
              "catagory":   "ADP",
              "type":       "integer",
              "range":      PROP_PROBABILTY_RANGE,
              "default":    50,
              "action":     ["recompile_Uaxes", "recompile_Uellipse"] })
        self.glo_add_property(
            { "name":      "U",
              "desc":      "Show Thermal Axes",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   False,
              "action":    "recompile_Uaxes" })
        self.glo_add_property(
            { "name":      "U_color",
              "desc":      "Thermal Axes Color",
              "catagory":  "ADP",
              "type":      "color",
              "default":   (1.0, 1.0, 1.0),
              "action":    "recompile_Uaxes" })
        self.glo_add_property(
            { "name":       "ellipse",
              "desc":       "Show Thermal Ellipseoids",
              "catagory":   "Show/Hide",
              "type":       "boolean",
              "default":    False,
              "action":     "recompile_Uellipse" })
        self.glo_add_property(
            { "name":       "ellipse_opacity",
              "desc":       "Thermal Ellipseoid Opacity",
              "catagory":   "ADP",
              "type":       "float",
              "range":      PROP_OPACITY_RANGE,
              "default":    1.0,
              "action":     "recompile_Uellipse" })
        self.glo_add_property(
            { "name":       "rms",
              "desc":       "Show Thermal Peanuts (RMS Deviation Surface)",
              "catagory":   "Show/Hide",
              "type":       "boolean",
              "default":    False,
              "action":     "recompile_Urms" })
        self.glo_add_property(
            { "name":       "rms_opacity",
              "desc":       "Peanut Surface Opacity",
              "catagory":   "ADP",
              "type":       "float",
              "range":      PROP_OPACITY_RANGE,
              "default":    1.0,
              "action":     "recompile_Urms" })

    def gldl_install_draw_methods(self):
        self.gldl_draw_method_install(
            { "name":                "labels",
              "func":                self.draw_labels,
              "no_gl_compile":          True,
              "transparent":         False,
              "visible_property":    "labels",
              "recompile_action":    "recompile_labels" })
        self.gldl_draw_method_install(
            { "name":                "lines",
              "func":                self.draw_lines,
              "transparent":         False,
              "visible_property":    "lines",
              "recompile_action":    "recompile_lines" })
        self.gldl_draw_method_install(
            { "name":                "trace",
              "func":                self.draw_trace,
              "transparent":         False,
              "visible_property":    "trace",
              "recompile_action":    "recompile_trace" })
        self.gldl_draw_method_install(
            { "name":                "ball_stick",
              "func":                self.draw_ball_stick,
              "transparent":         False,
              "visible_property":    "ball_stick",
              "recompile_action":    "recompile_ball_stick" })
        self.gldl_draw_method_install(
            { "name":                "cpk",
              "func":                self.draw_cpk,
              "visible_property":    "cpk",
              "opacity_property":    "cpk_opacity",
              "recompile_action":    "recompile_cpk" })
        self.gldl_draw_method_install(
            { "name":                "Uaxes",
              "func":                self.draw_Uaxes,
              "transparent":         False,
              "visible_property":    "U",
              "recompile_action":    "recompile_Uaxes" })
        self.gldl_draw_method_install(
            { "name":                "Uellipse",
              "func":                self.draw_Uellipse,
              "visible_property":    "ellipse",
              "opacity_property":    "ellipse_opacity",
              "recompile_action":    "recompile_Uellipse" })
        self.gldl_draw_method_install(
            { "name":                "Urms",
              "func":                self.draw_Urms,
              "visible_property":    "rms",
              "opacity_property":    "rms_opacity",
              "recompile_action":    "recompile_Urms" })
        
    def glal_update_properties(self, updates, actions):
        ## rebuild visible/hidden dictionaries
        if "recalc_positions" in actions:
             self.glal_hidden_atoms_dict  = None
             self.glal_visible_atoms_dict = None
             self.glal_xyzdict           = None

        ## update color enumeration
        if "color" in updates:
            self.glal_atom_color.set_value(updates["color"])

            if updates["color"].lower()=="random":
                new_color_str = "%4.2f,%4.2f,%4.2f" % (
                    self.glal_atom_color.value[0],
                    self.glal_atom_color.value[1],
                    self.glal_atom_color.value[2])

                self.properties.update(color=new_color_str)
                return
                
            self.glo_update_properties(
                color_setting = str(self.glal_atom_color.get_value_desc()))

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

    def glal_iter_atoms(self):
        """Implement in a subclass to iterate over all atoms which
        need to be drawn.
        """
        pass

    def glal_iter_fragments(self):
        """Implement in a subclass to iterate one Fragment object
        at a time, in order.  This implementation works with any
        implementation of glal_iter_atoms, but is very inefficent.
        """
        struct = Structure()
        memo   = {}
        
        for atm in self.glal_iter_atoms():
            struct.add_atom(copy.deepcopy(atm, memo))

        for frag in struct.iter_fragments():
            yield frag

    def glal_iter_chains(self):
        """Implement in a subclass to iterate one Chain object
        at a time, in order.  This implementation works with any
        implementation of glal_iter_atoms, but is very inefficent.
        """
        struct = Structure()
        memo   = {}
        
        for atm in self.glal_iter_atoms():
            struct.add_atom(copy.deepcopy(atm, memo))

        for chain in struct.iter_chains():
            yield chain

    def glal_iter_models(self):
        """Implement in a subclass to iterate one Model object
        at a time, in order.  This implementation works with any
        implementation of glal_iter_atoms, but is very inefficent.
        """
        struct = Structure()
        memo   = {}
        
        for atm in self.glal_iter_atoms():
            struct.add_atom(copy.deepcopy(atm, memo))

        for model in struct.iter_models():
            yield model

    def glal_iter_atoms_filtered(self):
        """Iterate all atoms and yield the tuble (atom, visible_flag).
        """
        aa_bb_atoms = ("N", "CA", "C", "O")
        na_bb_atoms = ("P", "O5*", "C5*", "C4*", "C3*", "O3*")

        main_chain_visible = self.properties["main_chain_visible"]
        side_chain_visible = self.properties["side_chain_visible"]
        hetatm_visible     = self.properties["hetatm_visible"]
        water_visible      = self.properties["water_visible"]
        hydrogen_visible   = self.properties["hydrogen_visible"]

        for atm in self.glal_iter_atoms():
            if hydrogen_visible==False:
                if atm.element=="H":
                    yield atm, False
                    continue

            frag = atm.get_fragment()
            
            if frag.is_amino_acid()==True:
                if main_chain_visible==True and side_chain_visible==True:
                    yield atm, True
                elif main_chain_visible==True and side_chain_visible==False:
                    if atm.name in aa_bb_atoms:
                        yield atm, True
                elif main_chain_visible==False and side_chain_visible==True:
                    if atm.name not in aa_bb_atoms:
                        yield atm, True
                yield atm, False

            elif frag.is_nucleic_acid()==True:
                if main_chain_visible==True and side_chain_visible==True:
                    yield atm, True
                elif main_chain_visible==True and side_chain_visible==False:
                    if atm.name in na_bb_atoms:
                        yield atm, True
                elif main_chain_visible==False and side_chain_visible==True:
                    if atm.name not in na_bb_atoms:
                        yield atm, True
                yield atm, False

            elif frag.is_water()==True:
                if water_visible==True:
                    yield atm, True
                yield atm, False

            elif hetatm_visible==True:
                yield atm, True

            else:
                yield atm, False

    def glal_rebuild_atom_dicts(self):
        """When a atom selection setting or origin changes, the atom
        dictionaries need to be rebuilt.
        """
        self.glal_hidden_atoms_dict  = {}
        self.glal_visible_atoms_dict = {}
        self.glal_xyzdict           = XYZDict(5.0)

        for atm, visible in self.glal_iter_atoms_filtered():
            pos = self.glal_calc_position(atm.position)

            if visible==True:
                self.glal_visible_atoms_dict[atm] = pos
                self.glal_xyzdict.add(pos, atm)
            else:
                self.glal_hidden_atoms_dict[atm] = pos

    def glal_iter_visible_atoms(self):
        """Iterate over all visible atoms yielding the 2-tuple (atm, position).
        """
        if self.glal_visible_atoms_dict==None:
            self.glal_rebuild_atom_dicts()

        for atm, pos in self.glal_visible_atoms_dict.items():
            yield atm, pos
                    
    def glal_calc_position(self, position):
        """Calculate a position vector with respect to the
        proeprty: atom_origin.
        """
        if self.properties["atom_origin"]!=None:
             return position - self.properties["atom_origin"]
        return position
    
    def glal_calc_color(self, atom):
        """Sets the open-gl color for the atom.
        """
        return self.glal_atom_color.calc_atom_color(atom)

    def glal_calc_U(self, atom):
        """Return the ADP U tensor for the atom
        """
        return atom.get_U()

    def draw_labels(self):
        """Draws atom lables.
        """
        auto_label_dist = self.properties["auto_label_distance"]

        viewer = self.gldl_get_glviewer()
        R = viewer.properties["R"]
        t = viewer.properties["t"]

        Ri = transpose(R)
        tp = matrixmultiply(R, t)

        ## this vecor is needed to adjust for origin/atom_origin
        ## changes, but it does not account for object rotation yet
        cv = self.properties["origin"]
        
        glPushMatrix()

        ## shift back to the view window coordinate system so the
        ## labels can be drawn in the plane perpendicular to the viewer's
        ## z axis
        self.glr_translate(-t)
        self.glr_mult_matrix(Ri)

        for atm, pos in self.glal_iter_visible_atoms():
            relative_pos = matrixmultiply(R, pos + cv) + tp
            if length(relative_pos)>auto_label_dist:
                continue

            if atm.alt_loc=="":
                text = "%s %s %s %s" % (
                    atm.name, atm.res_name, atm.fragment_id, atm.chain_id)
            else:
                text = "%s(%s) %s %s %s" % (
                    atm.name, atm.alt_loc, atm.res_name, atm.fragment_id,
                    atm.chain_id)

            glPushMatrix()
            glTranslatef(*relative_pos + array([0.0, 0.0, 0.5]))
            self.glr_text(text)
            glPopMatrix()

        glPopMatrix()

    def draw_cpk(self):
        """Draw a atom as a CPK sphere.
        """
        for atm, pos in self.glal_iter_visible_atoms():
            elem = atm.get_structure().library.get_element(atm.element)
            if elem:
                radius = elem.van_der_waals_radius
            else:
                radius = 2.0

            r, g, b = self.glal_calc_color(atm)

            a = self.properties["cpk_opacity"]
            self.glr_set_material_rgb(r, g, b, a)

            self.glr_sphere(
                pos,
                self.properties["cpk_scale_radius"] * radius,
                self.properties["sphere_quality"])

    def draw_Uaxes(self):
        """Draw thermal axes at the given ADP probability level.
        """
        for atm, pos in self.glal_iter_visible_atoms():
            U = self.glal_calc_U(atm)
            if U==None:
                continue

            self.glr_Uaxes(
                pos, U,
                self.properties["adp_prob"],
                self.properties["U_color"],
                1.0)

    def draw_Uellipse(self):
        """Draw the ADP determined probability ellipseoid.
        """
        for atm, pos in self.glal_iter_visible_atoms():
            U = self.glal_calc_U(atm)
            if U==None:
                continue

            r, g, b = self.glal_calc_color(atm) 
            self.glr_set_material_rgb(
                r, g, b, self.properties["ellipse_opacity"])
            self.glr_Uellipse(pos, U, self.properties["adp_prob"])

    def draw_Urms(self):
        """Draw the ADP determined RMS displacement surface.
        """
        for atm, pos in self.glal_iter_visible_atoms():
            U = self.glal_calc_U(atm)
            if U==None:
                continue

            r, g, b = self.glal_calc_color(atm)        
            self.glr_set_material_rgb(r, g, b, self.properties["rms_opacity"])
            self.glr_Urms(pos, U)

    def draw_lines(self):
        """Draw a atom using bond lines only.
        """
        glDisable(GL_LIGHTING)
        glLineWidth(self.properties["line_width"])

        for atm1, pos1 in self.glal_iter_visible_atoms():
            glColor3f(*self.glal_calc_color(atm1))

            if len(atm1.bond_list)>0:
                ## if there are bonds, then draw the lines 1/2 way to the
                ## bonded atoms
                for bond in atm1.iter_bonds():
                    atm2 = bond.get_partner(atm1)

                    try:
                        pos2 = self.glal_visible_atoms_dict[atm2]
                    except KeyError:
                        if self.glal_hidden_atoms_dict.has_key(atm2):
                            continue
                        else:
                            pos2 = self.glal_calc_position(atm2.position)
                    
                    end = pos1 + ((pos2 - pos1) / 2)

                    glBegin(GL_LINES)
                    glVertex3f(*pos1)
                    glVertex3f(*end)
                    glEnd()
                    
            else:
                self.glr_cross(
                    pos1,
                    self.glal_calc_color(atm1),
                    self.properties["line_width"])

    def draw_ball_stick(self):
        """Draw atom with ball/stick model.
        """
        ball_radius  = self.properties["ball_radius"]
        stick_radius = self.properties["stick_radius"]

        for atm1, pos1 in self.glal_iter_visible_atoms():
            r, g, b = self.glal_calc_color(atm1)
            self.glr_set_material_rgb(r, g, b, 1.0)

            ## if there are bonds, then draw the lines 1/2 way to the
            ## bonded atoms
            for bond in atm1.iter_bonds():
                atm2 = bond.get_partner(atm1)

                try:
                    pos2 = self.glal_visible_atoms_dict[atm2]
                except KeyError:
                    if self.glal_hidden_atoms_dict.has_key(atm2):
                        continue
                    else:
                        pos2 = self.glal_calc_position(atm2.position)

                end = pos1 + ((pos2 - pos1) / 2)
                self.glr_tube(pos1, end, stick_radius)

            ## draw ball
            self.glr_sphere(pos1, ball_radius, 10)

    def draw_Udiff(self, atoms):
        ball_radius  = self.properties["ball_radius"]
        stick_radius = self.properties["stick_radius"]

        visited = {}

        for atm, pos in atoms.items():
            ## if there are bonds, then draw the lines 1/2 way to the
            ## bonded atoms
            for bond in atm.iter_bonds():

                if visited.has_key(bond):
                    continue

                atm2 = bond.get_partner(atm)

                try:
                    pos2 = atoms[atm2]
                except KeyError:
                    continue

                v    = pos2 - pos
                n    = normalize(v)
                end  = pos + v

                U  = atm.get_U()
                U2 = atm2.get_U()

                rms  = math.sqrt(matrixmultiply(n, matrixmultiply(U, n)))
                rms2 = math.sqrt(matrixmultiply(n, matrixmultiply(U2, n)))

                rms_diff = abs(rms - rms2)

                rms_diff = rms_diff  * 50.0
                self.glr_set_material_rgb(
                    1.0 - rms_diff, 1.0 - rms_diff, 1.0, 1.0)

                self.glr_tube(pos, end, stick_radius)

            self.glr_set_material_rgb(1.0, 1.0, 1.0, 1.0)
            self.glr_sphere(pos, ball_radius, 10)

    def draw_Udiff2(self, atoms):
        ball_radius  = self.properties["ball_radius"] + 0.01
        stick_radius = self.properties["stick_radius"] + 0.01

        visited = {}
        atom_list = atoms.items()

        for atm1, pos1 in atom_list:
            visited[atm1] = True
            
            for atm2, pos2 in atom_list:
                if visited.has_key(atm2):
                    continue

                if atm1.alt_loc!=atm2.alt_loc:
                    if atm1.alt_loc!="" and atm2.alt_loc!="":
                        continue

                v = pos2 - pos1

                if length(v) > 5.0:
                    continue

                if atm1.get_bond(atm2)==None:
                    tube_radius = 0.25 * stick_radius
                else:
                    tube_radius = stick_radius
                
                n    = normalize(v)
                end  = pos1 + v

                U1 = atm1.get_U()
                U2 = atm2.get_U()

                rms1 = math.sqrt(matrixmultiply(n, matrixmultiply(U1, n)))
                rms2 = math.sqrt(matrixmultiply(n, matrixmultiply(U2, n)))

                rms_diff = abs(rms1 - rms2)

                if rms_diff > 0.005:
                    continue

                rms_diff = rms_diff  * 20.0
                self.glr_set_material_rgb(
                    1.0 - rms_diff, 1.0 - rms_diff, 1.0, 1.0)

                self.glr_tube(pos1, end, tube_radius)

            self.glr_set_material_rgb(1.0, 1.0, 1.0, 1.0)
            self.glr_sphere(pos1, ball_radius, 10)

    def draw_cross(self, atm, pos):
        """Draws atom with a cross of lines.
        """
        self.glr_cross(
            pos,
            self.glal_calc_color(atm),
            self.properties["line_width"])

    def draw_trace(self):
        """Draws trace over all polymer backbone atoms.
        """
        trace_radius = self.properties["trace_radius"]
        r, g, b      = self.properties["trace_color"]
        
        self.glr_set_material_rgb(r, g, b, 1.0)

        for chain in self.glal_iter_chains():

            last_atm = None

            for frag in chain.iter_fragments():

                if frag.is_amino_acid()==True:
                    backbone_atoms = ("N", "CA", "C")
                elif frag.is_nucleic_acid()==True:
                    backbone_atoms = ("P", "O5*", "C5*", "C4*", "C3*", "O3*")
                else:
                    last_atm = None
                    continue

                for name in backbone_atoms:
                    try:
                        atm = frag[name]
                    except KeyError:
                        last_atm = None
                        continue

                    if last_atm==None:
                        last_atm = atm
                        continue

                    ## if there are alternate conformations, make sure to
                    ## trace them all in the backbone trace
                    
                    if last_atm.alt_loc=="" and atm.alt_loc=="":
                        lpos = self.glal_calc_position(
                            last_atm.position)
                        pos = self.glal_calc_position(
                            atm.position)

                        self.glr_sphere(lpos, trace_radius, 12)
                        self.glr_tube(lpos, pos, trace_radius)

                    elif last_atm.alt_loc=="" and atm.alt_loc!="":
                        lpos = self.glal_calc_position(
                            last_atm.position)

                        for aa in atm.iter_alt_loc():
                            pos = self.glal_calc_position(
                                aa.position)

                            self.glr_sphere(lpos, trace_radius, 12)
                            self.glr_tube(lpos, pos, trace_radius)

                    elif last_atm.alt_loc!="" and atm.alt_loc=="":
                        pos = self.glal_calc_position(
                            atm.position)

                        for laa in last_atm.iter_alt_loc():
                            lpos = self.glal_calc_position(
                                laa.position)

                            self.glr_sphere(lpos, trace_radius, 12)
                            self.glr_tube(lpos, pos, trace_radius)

                    elif last_atm.alt_loc!="" and atm.alt_loc!="":
                        for aa in atm.iter_alt_loc():
                            for laa in last_atm.iter_alt_loc():

                                if aa.alt_loc!=laa.alt_loc:
                                    continue

                                lpos = self.glal_calc_position(
                                    laa.position)
                                pos = self.glal_calc_position(
                                    aa.position)

                                self.glr_sphere(lpos, trace_radius, 12)
                                self.glr_tube(lpos, pos, trace_radius)

                    last_atm = atm

            if last_atm!=None:
                for laa in last_atm.iter_alt_loc():
                    lpos = self.glal_calc_position(laa.position)
                    self.glr_sphere(lpos, trace_radius, 10)

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

        
class GLTLSAtomList(GLAtomList):
    """OpenGL visualizations of TLS group atoms.
    """
    def __init__(self, **args):
        self.tls_group = args["tls_group"]
        
        GLAtomList.__init__(self, **args)
        self.glo_init_properties(**args)
    
    def glo_install_properties(self):
        GLAtomList.glo_install_properties(self)

        ## Show/Hide
        self.glo_add_property(
            { "name":        "fan_visible",
              "desc":        "Show Fans from COR to Backbone",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     True,
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

        ## TLS
        self.glo_add_property(
            { "name":        "fan_color",
              "desc":        "COR/Backbone Color",
              "catagory":    "TLS",
              "type":        "color",
              "default":     (0.5, 1.0, 0.5),
              "action":      "recompile_fan" })
        self.glo_add_property(
            { "name":        "fan_opacity",
              "desc":        "COR/Backbone Fan Opacity",
              "catagory":    "TLS",
              "type":        "float",
              "range":       PROP_OPACITY_RANGE,
              "default":     0.25,
              "action":      "recompile_fan" })

        self.glo_add_property(
            { "name":        "L1_scale",
              "desc":        "Scale L1 Rotation", 
              "catagory":    "TLS",
              "type":        "float",
              "default":     1.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L2_scale",
              "desc":        "Scale L2 Rotation", 
              "catagory":    "TLS",
              "type":        "float",
              "default":     1.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L3_scale",
              "desc":        "Scale L3 Rotation", 
              "catagory":    "TLS",
              "type":        "float",
              "default":     1.0,
              "action":      "redraw" })

        ## TLS Analysis
        self.glo_add_property(
            { "name":        "COR",
              "desc":        "TLS Center of Reaction", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "T",
              "desc":        "Translation Tensor (A*A)",
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3,3)",
              "default":     zeros((3,3), Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L",
              "desc":        "Libration Tensor (DEG*DEG)",
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3,3)",
              "default":     zeros((3,3), Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "S",
              "desc":        "Skew Tensor (A*DEG)",
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3,3)",
              "default":     zeros((3,3), Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_eigen_vec",
              "desc":        "L1 Eigen Vector", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_eigen_vec",
              "desc":        "L2 Eigen Vector", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_eigen_vec",
              "desc":        "L3 Eigen Vector", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_eigen_val",
              "desc":        "L1 Eigen Value", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_eigen_val",
              "desc":        "L2 Eigen Value", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_eigen_val",
              "desc":        "L3 Eigen Value", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_rho",
              "desc":        "L1 translation from COR", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_rho",
              "desc":        "L2 translation from COR", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_rho",
              "desc":        "L3 translation from COR", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_pitch",
              "desc":        "L1 screw pitch (A/DEG)", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_pitch",
              "desc":        "L2 screw pitch (A/DEG)", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_pitch",
              "desc":        "L3 screw pitch (A/DEG)", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })

        ## Simulation State
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

    def gldl_install_draw_methods(self):
        GLAtomList.gldl_install_draw_methods(self)
        
        self.gldl_draw_method_install(
            { "name":                "fan",
              "func":                self.draw_fan,
              "visible_property":    "fan_visible",
              "opacity_property":    "fan_opacity",
              "recompile_action":    "recompile_fan" })

    def gldl_iter_multidraw_self(self):
        for draw_flag in GLAtomList.gldl_iter_multidraw_self(self):
            for draw_flag2 in self.gldl_iter_multidraw_animate():
                yield True
            
    def gldl_iter_multidraw_animate(self):
        """
        """
        for Lx_axis, Lx_rho, Lx_pitch, Lx_rot, Lx_scale in (
            ("L1_eigen_vec", "L1_rho", "L1_pitch", "L1_rot", "L1_scale"),
            ("L2_eigen_vec", "L2_rho", "L2_pitch", "L2_rot", "L2_scale"),
            ("L3_eigen_vec", "L3_rho", "L3_pitch", "L3_rot", "L3_scale") ):

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
            Lx_rot   = self.properties[Lx_rot] * self.properties[Lx_scale]

            screw = Lx_axis * Lx_rot * Lx_pitch

            glPushMatrix()
            glTranslatef(*Lx_rho + screw)
            glRotatef(Lx_rot, *Lx_axis)
            glTranslatef(*-Lx_rho)
            yield True
            glPopMatrix()

    def glal_iter_atoms(self):
        for atm in self.tls_group:
            yield atm

    def glal_calc_U(self, atom):
        """Always return the reduced T tensor.
        """
        return self.properties["T"]

    def draw_fan(self):
        """Draws a fan from the TLS group center of reaction to the
        TLS group backbone atoms.
        """
        glEnable(GL_LIGHTING)

        COR     = self.properties["COR"]
        r, g, b = self.properties["fan_color"]
        a       = self.properties["fan_opacity"]

        self.glr_set_material_rgb(r, g, b, a)
        
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE)
        glBegin(GL_TRIANGLE_FAN)

        v1 = None
        v2 = None

        for atm in self.tls_group:
            if atm.name not in ("N", "CA", "C"):
                continue

            if v1==None:
                v1 = atm.position - COR
                continue
            elif v2==None:
                v2 = atm.position - COR
                glNormal3f(*cross(v1, v2))
                glVertex3f(0.0, 0.0, 0.0)                
            else:
                v1 = v2
                v2 = atm.position - COR

            glNormal3f(*cross(v1, v2))    
            glVertex3f(*v1)
            glVertex3f(*v2)

        glEnd()

        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE)


class GLTLSGroup(GLDrawList):
    """Top level visualization object for a TLS group.
    """
    def __init__(self, **args):
        GLDrawList.__init__(self)

        orig_tls_group = args["tls_group"]
        self.orig_tls_group = orig_tls_group

        if orig_tls_group.name!="":
            self.glo_set_name("TLS Group: %s" % (orig_tls_group.name))
        
        ## TLS calculations        
        ## step 1: copy the TLS group
        calcs                 = orig_tls_group.calc_COR()

        ## create this object's TLS group at the center of reaction
        self.tls_group        = TLSGroup(orig_tls_group)
        self.tls_group.origin = calcs["COR"].copy()
        self.tls_group.T      = calcs["T'"].copy()
        self.tls_group.L      = calcs["L'"].copy()
        self.tls_group.S      = calcs["S'"].copy()

        L_eigen_vals, L_eigen_vecs = eigenvectors(self.tls_group.L)

        ## step 4: add child GLTLSAtomList 
        self.gl_atom_list = GLTLSAtomList(
            tls_group        = self.tls_group,
            origin           = calcs["COR"].copy(),
            atom_origin      = calcs["COR"].copy(),
            trace            = True,
            trace_line_width = 10.0,
            trace_color      = (0.5, 1.0, 0.5),
            lines            = False,
            fan_visible      = True)

        self.gl_atom_list.glo_set_name("TLS Atom Animation")
        self.gl_atom_list.glo_set_properties_id("gl_atom_list")
        self.glo_add_child(self.gl_atom_list)
        
        self.glo_link_child_property(
            "animation_visible", "gl_atom_list", "visible")

        self.glo_link_child_property(
            "COR", "gl_atom_list", "COR")
        self.glo_link_child_property(
            "T_reduced", "gl_atom_list", "T")
        self.glo_link_child_property(
            "L", "gl_atom_list", "L")
        self.glo_link_child_property(
            "S", "gl_atom_list", "S")

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
            "symmetry", "gl_atom_list", "symmetry")
 
        ## initalize properties
        self.glo_add_update_callback(self.tls_update_cb)
        self.glo_init_properties(
            COR          = calcs["COR"].copy(),
            T            = self.tls_group.T.copy(),
            T_reduced    = calcs["rT'"].copy(),
            L            = self.tls_group.L.copy() * rad2deg2,
            S            = self.tls_group.S.copy() * rad2deg,
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

        ## TLS Analysis
        self.glo_add_property(
            { "name":        "COR",
              "desc":        "TLS Center of Reaction",
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "T",
              "desc":        "Translation Tensor (A*A)",
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3,3)",
              "default":     zeros((3,3), Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "T_reduced",
              "desc":        "Reduced Translation Tensor (A*A)",
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3,3)",
              "default":     zeros((3,3), Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L",
              "desc":        "Libration Tensor (DEG*DEG)",
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3,3)",
              "default":     zeros((3,3), Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "S",
              "desc":        "Skew Tensor (A*DEG)",
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3,3)",
              "default":     zeros((3,3), Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_eigen_vec",
              "desc":        "L1 Eigen Vector",
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_eigen_vec",
              "desc":        "L2 Eigen Vector",
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_eigen_vec",
              "desc":        "L3 Eigen Vector", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_eigen_val",
              "desc":        "L1 Eigen Value", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_eigen_val",
              "desc":        "L2 Eigen Value", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_eigen_val",
              "desc":        "L3 Eigen Value", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_rho",
              "desc":        "L1 Translation Vector from COR", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_rho",
              "desc":        "L2 Translation Vector from COR", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_rho",
              "desc":        "L3 Translation Vector from COR", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_pitch",
              "desc":        "L1 screw pitch (A/DEG)", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L2_pitch",
              "desc":        "L2 screw pitch (A/DEG)", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L3_pitch",
              "desc":        "L3 screw pitch (A/DEG)", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "L1_rot",
              "desc":        "L1 Axis Viewing Libration (DEG)",
              "catagory":    "Simulation State",
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L2_rot",
              "desc":        "L2 Axis Viewing Libration (DEG)", 
              "catagory":    "Simulation State",
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L3_rot",
              "desc":        "L2 Axis Viewing Libration (DEG)",
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

        ## Show/Hide
        self.glo_add_property(
            { "name":        "symmetry",
              "desc":        "Show Symmetry Equivelant",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "TLS_visible",
              "desc":        "Show TLS Tensors",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     True,
              "action":      "recompile_tensors" })
        self.glo_add_property(
            { "name":        "main_chain_only",
              "desc":        "Show Protein Mainchain Atoms Only",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     True,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "Utls_U_diff",
              "desc":        "Show Utls vs. Experimental U Diff",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile_Utls_ellipse" })
        self.glo_add_property(
            { "name":        "U",
              "desc":        "Show TLS Thermal Axes",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile_Utls_axes" })
        self.glo_add_property(
            { "name":        "ellipse",
              "desc":        "Show TLS Thermal Ellipsoids",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile_Utls_ellipse" })
        self.glo_add_property(
            { "name":        "rms",
              "desc":        "Show TLS Thermal Peanuts",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile_Utls_rms" })
        self.glo_add_property(
            { "name":        "animation_visible",
              "desc":        "Show TLS Animated Atoms", 
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     True,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "L1_visible",
              "desc":        "Show Screw L1 Displacement Surface", 
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile_surface" })
        self.glo_add_property(
            { "name":        "L2_visible",
              "desc":        "Show Screw L1 Displacement Surface", 
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile_surface" })
        self.glo_add_property(
            { "name":        "L3_visible",
              "desc":        "Show Screw L1 Displacement Surface",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile_surface" })
        
        ## TLS
        self.glo_add_property(
            { "name":        "add_biso",
              "desc":        "Add Atom Biso to Utls (Required for REFMAC)",
              "catagory":    "TLS",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":       "adp_prob",
              "desc":       "Contour Probability",
              "catagory":   "TLS",
              "type":       "integer",
              "range":      PROP_PROBABILTY_RANGE,
              "default":    50,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "T_line_width",
              "desc":       "T Axes Line Width",
              "catagory":   "TLS",
              "type":       "float",
              "default":    1.0,
              "action":     "recompile_tensors" })
        self.glo_add_property(
            { "name":        "T_color",
              "desc":        "T Axes Color",
              "catagory":    "TLS",
              "type":        "color",
              "default":     (0.5, 0.5, 1.0),
              "action":      "recompile_tensors" })
        self.glo_add_property(
            { "name":       "L_axis_radius",
              "desc":       "Screw Axes Radius",
              "catagory":   "TLS",
              "type":       "float",
              "default":    0.05,
              "action":     "recompile_tensors" })
        self.glo_add_property(
            { "name":       "L_color",
              "desc":       "Screw Axes Color",
              "catagory":   "TLS",
              "type":       "color",
              "default":    (1.0, 0.5, 0.5),
              "action":     "recompile_tensors" })
        self.glo_add_property(
            { "name":        "Utls_color",
              "desc":        "Thermal Utls Color",
              "catagory":    "TLS",
              "type":        "color",
              "default":     (0.10, 0.75, 0.35),
              "action":      ["recompile_Utls_axes", "recompile_Utls_ellipse",
                              "recompile_Utls_rms"] })
        self.glo_add_property(
            { "name":        "U_color",
              "desc":        "Thermal Atom U Color for Diff Mode",
              "catagory":    "TLS",
              "type":        "color",
              "default":     (1.0, 1.0, 1.0),
              "action":      "recompile_Utls_ellipse" })
        self.glo_add_property(
            { "name":        "ellipse_opacity",
              "desc":        "TLS Thermal Ellipseoid Opacity",
              "catagory":    "TLS",
              "type":        "float",
              "range":       PROP_OPACITY_RANGE,
              "default":     1.0,
              "action":      "recompile_Utls_ellipse" })
        self.glo_add_property(
            { "name":        "rms_opacity",
              "desc":        "TLS Thermal Peanut Opacity",
              "catagory":    "TLS",
              "type":        "float",
              "range":       PROP_OPACITY_RANGE,
              "default":     1.0,
              "action":      "recompile_Utls_rms" })
        self.glo_add_property(
            { "name":        "surface_opacity",
              "desc":        "Screw Surface Opacity",
              "catagory":    "TLS",
              "type":        "float",
              "range":      PROP_OPACITY_RANGE,
              "default":     1.0,
              "action":      "recompile_surface" })

    def gldl_install_draw_methods(self):
        self.gldl_draw_method_install(
            { "name":                "tls_tensors",
              "func":                self.draw_tensors,
              "transparent":         False,
              "visible_property":    "TLS_visible",
              "recompile_action":    "recompile_tensors" })
        self.gldl_draw_method_install(
            { "name":                "Utls_axes",
              "func":                self.draw_Utls_axes,
              "transparent":         False,
              "visible_property":    "U",
              "recompile_action":    "recompile_Utls_axes" })
        self.gldl_draw_method_install(
            { "name":                "Utls_ellipse",
              "func":                self.draw_Utls_ellipse,
              "visible_property":    "ellipse",
              "opacity_property":    "ellipse_opacity",
              "recompile_action":    "recompile_Utls_ellipse" })
        self.gldl_draw_method_install(
            { "name":                "Utls_rms",
              "func":                self.draw_Utls_rms,
              "visible_property":    "rms",
              "opacity_property":    "rms_opacity",
              "recompile_action":    "recompile_Utls_rms" })
        self.gldl_draw_method_install(
            { "name":                "L1_surface",
              "func":                self.draw_L1_surface,
              "visible_property":    "L1_visible",
              "opacity_property":    "surface_opacity",
              "recompile_action":    "recompile_surface" })
        self.gldl_draw_method_install(
            { "name":                "L2_surface",
              "func":                self.draw_L2_surface,
              "visible_property":    "L2_visible",
              "opacity_property":    "surface_opacity",
              "recompile_action":    "recompile_surface" })
        self.gldl_draw_method_install(
            { "name":                "L3_surface",
              "func":                self.draw_L3_surface,
              "visible_property":    "L3_visible",
              "opacity_property":    "surface_opacity",
              "recompile_action":    "recompile_surface" })
         
    def tls_update_cb(self, updates, actions):
        if "time" in updates:
            self.update_time()

    def update_time(self):
        """Changes the time of the TLS group simulating harmonic motion.
        """
        sin_tm = math.sin(3.0 * self.properties["time"] * 2 * math.pi)
        sin_tm = sin_tm

        ## L Tensor
        C = GAUSS3C[self.properties["adp_prob"]]

        try:
            L1_peak = C * math.sqrt(self.properties["L1_eigen_val"])
        except ValueError:
            L1_peak = 0.0

        try:
            L2_peak = C * math.sqrt(self.properties["L2_eigen_val"])
        except ValueError:
            L2_peak = 0.0

        try:
            L3_peak = C * math.sqrt(self.properties["L3_eigen_val"])
        except ValueError:
            L3_peak = 0.0

        L1_rot  = L1_peak * sin_tm
        L2_rot  = L2_peak * sin_tm
        L3_rot  = L3_peak * sin_tm

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

    def gltls_iter_atoms(self):
        """Special atom iterator for the TLS drawing functions yields:
        atm, Utls
        """
        for atm, Utls in self.tls_group.iter_atm_Utls():

            if self.properties["main_chain_only"]==True and \
               atm.name not in ("C", "CA", "N"):
                continue
            
            yield atm, Utls
    
    def draw_tensors(self):
        """Draw tensor axis.
        """
        glPushMatrix()
        glTranslatef(*self.properties["COR"])

        ## put a small sphere at the center of reaction
        r, g, b = self.properties["Utls_color"]        
        self.glr_set_material_rgb(r, g, b, 1.0)
        self.glr_sphere(array([0.0, 0.0, 0.0]), 0.05, 12)
        
        ## T: units (A^2)
        r, g, b = self.properties["T_color"]
        self.glr_set_material_rgb(r, g, b, 1.0)
        self.glr_Uellipse(
            (0.0, 0.0, 0.0),
            self.properties["T_reduced"],
            self.properties["adp_prob"])

        ## L: units (DEG^2)
        for Lx_eigen_val, Lx_eigen_vec, Lx_rho, Lx_pitch in [
            ("L1_eigen_val", "L1_eigen_vec", "L1_rho", "L1_pitch"),
            ("L2_eigen_val", "L2_eigen_vec", "L2_rho", "L2_pitch"),
            ("L3_eigen_val", "L3_eigen_vec", "L3_rho", "L3_pitch")]:


            L_eigen_vec = self.properties[Lx_eigen_vec]
            L_eigen_val = self.properties[Lx_eigen_val]
            L_rho       = self.properties[Lx_rho]
            L_pitch     = self.properties[Lx_pitch]
            
            if L_eigen_val<=0.0:
                continue

            C = GAUSS3C[self.properties["adp_prob"]]
            
            L_rot = C * math.sqrt(L_eigen_val)
            L_v   = L_eigen_vec * L_rot

            ## line from COR to center of screw/rotation axis
            glColor3f(*self.properties["L_color"])

            ## draw lines from COR to the axis
            glDisable(GL_LIGHTING)
            glBegin(GL_LINES)
            glVertex3f(0.0, 0.0, 0.0)
            glVertex3f(*L_rho)
            glEnd()

            ## draw axis
            r, g, b = self.properties["L_color"]
            self.glr_set_material_rgb(r, g, b, 1.0)
            self.glr_axis(L_rho - (0.5 * L_v), L_v,
                          self.properties["L_axis_radius"])

            ## draw disks with translational displacement
            L_screw_dis = L_eigen_vec * L_rot * L_pitch

            self.glr_axis(
                L_rho - (0.5 * L_screw_dis),
                L_screw_dis,
                1.5 * self.properties["L_axis_radius"])

        glPopMatrix()

    def draw_Utls_axes(self):
        """Render the anisotropic thremal axes calculated from the TLS
        model.
        """
        prob  = self.properties["adp_prob"]
        color = self.properties["Utls_color"]

        for atm, Utls in self.gltls_iter_atoms():
            self.glr_Uaxes(atm.position, Utls, prob, color, 1.0)

    def draw_Utls_ellipse(self):
        """Render the anisotropic thremal ellipseoids at the given probability
        contour calculated from the TLS model.
        """
        prob    = self.properties["adp_prob"]
        r, g, b = self.properties["Utls_color"]
        a       = self.properties["ellipse_opacity"]

        #dp2_mean, dp2_sig = self.tls_group.calc_mean_DP2()

        for atm, Utls in self.gltls_iter_atoms():

            if self.properties["add_biso"]==True:
                Uiso = identity(3, Float)*(atm.temp_factor / (8.0*math.pi**2))
                Utls = Utls + Uiso

##            dp2 = calc_DP2uij(atm.get_U(), Utls)

##             if dp2<=dp2_mean:
##                 dp2_scale = 1.0
##             elif dp2<=dp2_mean+dp2_sig:
##                 dp2_scale = 0.5
##             else:
##                 dp2_scale = 0.1
            
##             self.glr_set_material_rgb(
##                 r * dp2_scale,  g * dp2_scale,  b * dp2_scale, a)

            self.glr_set_material_rgb(r, g, b, a)
            self.glr_Uellipse(atm.position, Utls, prob)

            if self.properties["Utls_U_diff"]:
                r, g, b = self.properties["U_color"]
                self.glr_set_material_rgb(r, g, b, a)
                self.glr_Uellipse(
                    atm.position,
                    atm.get_U(),
                    self.properties["adp_prob"])

    def draw_Utls_rms(self):
        """Render the anisotropic thremal peanuts calculated from the TLS
        model.
        """
        r, g, b = self.properties["Utls_color"]
        a       = self.properties["rms_opacity"]
        self.glr_set_material_rgb(r, g, b, a)

        for atm, Utls in self.gltls_iter_atoms():
            self.glr_Urms(atm.position, Utls)

    def draw_L1_surface(self):
         self.draw_tls_surface(
             self.properties["L1_eigen_vec"],
             self.properties["L1_eigen_val"],
             self.properties["L1_rho"],
             self.properties["L1_pitch"])

    def draw_L2_surface(self):
         self.draw_tls_surface(
             self.properties["L2_eigen_vec"],
             self.properties["L2_eigen_val"],
             self.properties["L2_rho"],
             self.properties["L2_pitch"])

    def draw_L3_surface(self):
         self.draw_tls_surface(
             self.properties["L3_eigen_vec"],
             self.properties["L3_eigen_val"],
             self.properties["L3_rho"],
             self.properties["L3_pitch"])

    def draw_tls_surface(self, Lx_eigen_vec, Lx_eigen_val, Lx_rho, Lx_pitch):
        """Draws the TLS probability surface for a single non-intersecting
        screw axis.  Lx_eigen_val is the vaiance (mean square deviation MSD)
        of the rotation about the Lx_eigen_vec axis.
        """
        ## create a unique list of bonds which will be used to
        ## render the TLS surface; this list may be passed in a argument
        ## to avoid multiple calculations for each screw-rotation axis
        bond_list = []
        in_dict   = {}

        for atm, Utls in self.gltls_iter_atoms():
            in_dict[atm] = True

        for atm, Utls in self.gltls_iter_atoms():
            for bond in atm.iter_bonds():
                if in_dict.has_key(bond.get_partner(atm)):
                    bond_list.append(bond)
        
        ## this just won't work...
        if Lx_eigen_val==0.0:
            return

        C = GAUSS3C[self.properties["adp_prob"]]
        try:
            Lx_s = C * math.sqrt(Lx_eigen_val * deg2rad2)
        except ValueError:
            return

        COR           = self.properties["COR"]
        steps         = 10
        rot_step      = Lx_s / steps
        sq2pi         = math.sqrt(2.0 * math.pi)

        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE)
        glEnable(GL_LIGHTING)

        glBegin(GL_QUADS)

        for step in range(steps):
            step1 = rot_step * step
            step2 = rot_step * (1 + step)

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

                for bond in bond_list:

                    pos1 = bond.atom1.position - COR
                    pos2 = bond.atom2.position - COR

                    v1 = matrixmultiply(Rstep1, pos1) + rho_screw1
                    v2 = matrixmultiply(Rstep2, pos1) + rho_screw2
                    
                    v3 = matrixmultiply(Rstep2, pos2) + rho_screw2
                    v4 = matrixmultiply(Rstep1, pos2) + rho_screw1
                    
                    v14 = v1 + (0.5 * (v4 - v1))
                    v23 = v2 + (0.5 * (v3 - v2))

                    self.glr_set_material_atom(
                        bond.atom1, self.properties["surface_opacity"])

                    glNormal3f(*cross(v2-v1,v3-v1))

                    glVertex3f(*v1  + COR)
                    glVertex3f(*v2  + COR)
                    glVertex3f(*v23 + COR)
                    glVertex3f(*v14 + COR)
                    
                    self.glr_set_material_atom(
                        bond.atom2, self.properties["surface_opacity"])

                    glVertex3f(*v3  + COR)
                    glVertex3f(*v4  + COR)
                    glVertex3f(*v14 + COR)
                    glVertex3f(*v23 + COR)

        glEnd()

        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE)

            
class GLAtom(GLObject):
    def __init__(self, **args):
        GLObject.__init__(self, **args)
        self.atom = args["atom"]

        self.glo_init_properties(
            name     = self.atom.name,
            position = self.atom.position,
            **args)

    def glo_install_properties(self):
        GLObject.glo_install_properties(self)

        self.glo_add_property(
            { "name":        "name",
              "desc":        "Atom Name",
              "catagory":    "Atom",
              "type":        "string",
              "default":     None,

              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "position",
              "desc":        "Position",
              "catagory":    "Atom",
              "type":        "array(3)",
              "default":     None,
              "action":      "redraw" })

    def glo_name(self):
        return self.atom.name


class GLFragment(GLObject):
    def __init__(self, **args):
        GLObject.__init__(self, **args)
        self.fragment = args["fragment"]

        self.glo_init_properties(
            fragment_id = self.fragment.fragment_id,
            **args)

    def glo_install_properties(self):
        GLObject.glo_install_properties(self)

        self.glo_add_property(
            { "name":        "fragment_id",
              "desc":        "Fragment ID",
              "catagory":    "Fragment",
              "type":        "string",
              "default":     "",
              "action":      "redraw" })

    def glo_name(self):
        return "%5s %3s" % (self.fragment.fragment_id,
                            self.fragment.res_name)


class GLChain(GLAtomList):
    """Visualization object for mmLib.Structure.Chain.
    """
    def __init__(self, **args):
        GLAtomList.__init__(self, **args)
        self.chain = args["chain"]
        self.glo_init_properties(**args)

    def glo_install_properties(self):
        GLAtomList.glo_install_properties(self)

    def glo_name(self):
        return "Chain %s" % (self.chain.chain_id)

    def glal_iter_atoms(self):
        for frag in self.chain.iter_fragments():
            for atm in frag.iter_all_atoms():
                yield atm

    def glal_iter_fragments(self):
        """The GLAtomList implementation of this is slow.
        """
        for frag in self.chain.iter_fragments():
            yield frag

    def glal_iter_chains(self):
        """The GLAtomList implementation of this is slow.
        """
        yield self.chain


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

        ## init properties
        self.glo_init_properties(**args)

    def glo_install_properties(self):
        GLDrawList.glo_install_properties(self)

        self.glo_add_property(
            { "name":     "axes_visible",
              "desc":     "Show Cartesian Axes",
              "catagory": "Show/Hide",
              "type":     "boolean",
              "default":  True,
              "action":  "redraw" })
        self.glo_add_property(
            { "name":     "unit_cell_visible",
              "desc":     "Show Unit Cell",
              "catagory": "Show/Hide",
              "type":     "boolean",
              "default":  True,
              "action":   "redraw" })        

    def glo_name(self):
        return "%s" % (self.struct.cifdb.get_entry_id())

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


class GLViewer(GLObject, OpenGLRenderMethods):
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
        self.quat = array((0.0, 0.0, 0.0, 1.0), Float)
        
        GLObject.__init__(self)
        self.glo_add_update_callback(self.glv_update_cb)
        self.glo_init_properties()

    def glo_install_properties(self):
        GLObject.glo_install_properties(self)

        ## View
        self.glo_add_property(
            { "name":      "R",
              "desc":      "View Window Rotation Matrix",
              "catagory":  "View",
              "read_only": True,
              "type":      "array(3,3)",
              "default":   identity(3, Float),
              "action":    "redraw" })
        self.glo_add_property(
            { "name":      "cor",
              "desc":      "Center of Rotation",
              "catagory":  "View",
              "read_only": True,
              "type":      "array(3)",
              "default":   zeros(3, Float),
              "action":    "redraw" })
        self.glo_add_property(
            { "name":      "width",
              "desc":      "Window Width",
              "catagory":  "View",
              "read_only": True,
              "type":      "integer",
              "default":   1,
              "action":    "redraw" })
        self.glo_add_property(
            { "name":      "height",
              "desc":      "Window Height",
              "catagory":  "View",
              "read_only": True,
              "type":      "integer",
              "default":   1,
              "action":    "redraw" })
        self.glo_add_property(
            { "name":      "near",
              "desc":      "Near Clipping Plane",
              "catagory":  "View",
              "type":      "float",
              "default":   -20.0,
              "action":    "redraw" })
        self.glo_add_property(
            { "name":      "far",
              "desc":      "Far Clipping Plane",
              "catagory":  "View",
              "type":      "float",
              "default":   1000.0,
              "action":    "redraw" })
        self.glo_add_property(
            { "name":      "zoom",
              "desc":      "Zoom",
              "catagory":  "View",
              "type":      "float",
              "default":   20.0,
              "action":    "redraw" })
                
        ## OpenGL Lighting
        self.glo_add_property(
            { "name":      "GL_AMBIENT",
              "desc":      "Ambient Light",
              "catagory":  "OpenGL Lighting",
              "type":      "float",
              "range":     "0.0-1.0,0.1",
              "default":   0.2,
              "action":    "redraw" })
        self.glo_add_property(
            { "name":      "GL_SPECULAR",
              "desc":      "Specular Light",
              "catagory":  "OpenGL Lighting",
              "type":      "float",
              "range":     "0.0-1.0,0.1",
              "default":   1.0,
              "action":    "redraw" })
        self.glo_add_property(
            { "name":      "GL_DIFFUSE",
              "desc":      "Diffuse Light",
              "catagory":  "OpenGL Lighting",
              "type":      "float",
              "range":     "0.0-1.0,0.1",
              "default":   1.0,
              "action":    "redraw" })

        ## High-Performance OpenGL Features
        self.glo_add_property(
            { "name":      "GL_LINE_SMOOTH",
              "desc":      "Smooth Lines",
              "catagory":  "OpenGL Performance",
              "type":      "boolean",
              "default":   False,
              "action":    "redraw" })
        self.glo_add_property(
            { "name":      "GL_POINT_SMOOTH",
              "desc":      "Smooth Points",
              "catagory":  "OpenGL Performance",
              "type":      "boolean",
              "default":   False,
              "action":    "redraw" })
        self.glo_add_property(
            { "name":      "GL_POLYGON_SMOOTH",
              "desc":      "Smooth Polygons",
              "catagory":  "OpenGL Performance",
              "type":      "boolean",
              "default":   False,
              "action":    "redraw" })
        self.glo_add_property(
            { "name":      "GL_BLEND",
              "desc":      "Alpha-Blending (required for Fog)",
              "catagory":  "OpenGL Performance",
              "type":      "boolean",
              "default":   True,
              "action":    "redraw" })
        self.glo_add_property(
            { "name":      "GL_FOG",
              "desc":      "Enable Fog",
              "catagory":  "OpenGL Performance",
              "type":      "boolean",
              "default":   False,
              "action":    "redraw" })

        ## Fog Properties
        self.glo_add_property(
            { "name":      "GL_FOG_START",
              "desc":      "Fog Start Parameter",
              "catagory":  "Fog",
              "type":      "float",
              "default":   20.0,
              "action":    "redraw" })
        self.glo_add_property(
            { "name":      "GL_FOG_END",
              "desc":      "Fog End Parameter",
              "catagory":  "Fog",
              "type":      "float",
              "default":   50.0,
              "action":    "redraw" })

        
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
        draw_list.glo_remove()
        self.glv_redraw()

    def glv_add_struct(self, struct):
        """Adds the visualization for a mmLib.Structure.Structure object
        to the GLViewer.  It returns the GLStructure object created to
        visualize the Structure object.
        """
        assert isinstance(struct, Structure)

        ## add the structure
        gl_struct = GLStructure(struct=struct)
        self.glv_add_draw_list(gl_struct)
        
        ## select the center of the new structure as the rotation center
        n        = 0
        centroid = zeros(3, Float)

        max_x = 0.0
        max_y = 0.0
        max_z = 0.0

        min_x = 0.0
        min_y = 0.0
        min_z = 0.0

        for atm in struct.iter_atoms():
            n        += 1
            centroid += atm.position
            
            max_x = max(max_x, atm.position[0])
            max_y = max(max_y, atm.position[1])
            max_z = max(max_z, atm.position[2])

            min_x = min(min_x, atm.position[0])
            min_y = min(min_y, atm.position[1])
            min_z = min(min_z, atm.position[2])

        centroid = centroid / float(n)

        m1   = max(max_x, max(max_y, max_z))
        m2   = min(min_x, min(min_y, min_z))
        zoom = m1 - m2

        near = - zoom / 2.0
        far  =   zoom / 2.0

        self.properties.update(
            cor  = centroid,
            zoom = zoom,
            near = near,
            far  = far)
        
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
        pass

    def glv_resize(self, width, height):
        """Called to set the size of the OpenGL window this class is
        drawing on.
        """
        self.properties.update(width=width, height=height)

    def glv_clip(self, near, far):
        """Adjust near/far clipping planes.
        """
        width  = self.properties["width"]
        zoom   = self.properties["zoom"]

        angstrom_per_pixel = zoom / float(width)

        nearA = angstrom_per_pixel * float(near)
        farA  = angstrom_per_pixel * float(far)

        n = self.properties["near"] + nearA
        f = self.properties["far"]  + farA
        self.properties.update(near=n, far=f)

    def glv_zoom(self, z):
        """Adjust zoom levels.
        """
        width  = self.properties["width"]
        zoom   = self.properties["zoom"]

        angstrom_per_pixel = zoom / float(width)
        
        zoom = self.properties["zoom"]
        zoom += angstrom_per_pixel * float(z)
        
        if zoom<1.0:
            zoom = 1.0
        
        self.properties.update(zoom=zoom)

    def glv_straif(self, x, y):
        """Move in XY plane.
        """
        ## figure out A/pixel, multipy straif by pixes to get the
        ## the translation

        width  = self.properties["width"]
        zoom   = self.properties["zoom"]

        angstrom_per_pixel = zoom / float(width)

        xA = angstrom_per_pixel * float(x)
        yA = angstrom_per_pixel * float(y)

        ## XY translational shift
        dt = array((xA, yA, 0.0), Float)

        ## change the center of rotation
        R = self.properties["R"]

        ## shift in the XY plane by chainging the position of the
        ## center of rotation
        cor = self.properties["cor"] - matrixmultiply(transpose(R), dt)

        self.properties.update(cor=cor)

    def glv_trackball(self, x1, y1, x2, y2):
        """Virtual Trackball
        """
        def project_to_sphere(r, x, y):
            d = math.sqrt(x*x + y*y)
            if d<(r*0.707):
                return math.sqrt(r*r - d*d)
            else:
                return (r/1.414)**2 / d

        width  = self.properties["width"]
        height = self.properties["height"]

        ## determine the circle where the trackball is vs. the edges
        ## which are z-rotation
        square = min(width, height)
        radius = int(square * 0.7)

        x1c = x1 - width/2
        y1c = y1 - width/2
        d = math.sqrt(x1c*x1c + y1c*y1c)

        ## Z rotation
        if d>=radius:
            x2c = x2 - width/2
            y2c = y2 - width/2

            p1 = normalize(array((x1c, y1c, 0.0), Float))
            p2 = normalize(array((x2c, y2c, 0.0), Float))

            c = cross(p2, p1)

            a = normalize(c)
            theta = length(c) * math.pi/2.0

        ## XY trackball rotation
        else:

            x1 = (2.0*x1 - width)  / width
            y1 = (height - 2.0*y1) / height
            x2 = (2.0*x2 - width)  / width
            y2 = (height - 2.0*y2) / height

            tb_size = 1.0

            ## check for zero rotation
            if x1==x2 and y1==y2:
                return

            p1 = array((x1, y1, project_to_sphere(tb_size, x1, y1)), Float)
            p2 = array((x2, y2, project_to_sphere(tb_size, x2, y2)), Float)

            a = cross(p1, p2)
            d = p1 - p2
            t = length(d) / (2.0 * tb_size)

            if t>1.0:
                t - 1.0
            if t<-1.0:
                t = -1.0

            theta = 2.0 * math.asin(t)

        ## convert rotation axis a and rotation theta to a quaternion
        R = self.properties["R"]
        a = matrixmultiply(transpose(R), a)
        q = rquaternionu(a, theta)

        self.quat = addquaternion(q, self.quat)
        R         = rmatrixquaternion(self.quat)
        
        self.properties.update(R=R)
        
    def glv_pre_render(self):
        """Sets up lighting and OpenGL options before scene rendering.
        """
        ## OpenGL Features
        glEnable(GL_NORMALIZE)
        glEnable(GL_DEPTH_TEST)
        glDepthFunc(GL_LESS)
        
        ## Enables Smooth Color Shading
	#glShadeModel(GL_SMOOTH)
        
        ## inital orientation
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        
        ## lighting
        ambient  = (self.properties["GL_AMBIENT"],
                    self.properties["GL_AMBIENT"],
                    self.properties["GL_AMBIENT"],
                    1.0)
        diffuse  = (self.properties["GL_DIFFUSE"],
                    self.properties["GL_DIFFUSE"],
                    self.properties["GL_DIFFUSE"],
                    1.0)
        specular = (self.properties["GL_SPECULAR"],
                    self.properties["GL_SPECULAR"],
                    self.properties["GL_SPECULAR"],
                    1.0)

        ## use model abient light
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient)

        ## light 0
        glEnable(GL_LIGHT0)
        glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse)
        glLightfv(GL_LIGHT0, GL_SPECULAR, specular)
        glLightfv(GL_LIGHT0, GL_POSITION, (0.0, 0.0, 10.0, 0.0))
        
        ## light 1
        glDisable(GL_LIGHT1)
        glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse)
        glLightfv(GL_LIGHT1, GL_SPECULAR, specular)
        glLightfv(GL_LIGHT1, GL_POSITION, (0.0, 0.0, 100.0, 0.0))


        ## ANTI-ALIASING
        if self.properties["GL_LINE_SMOOTH"]==True:
            glEnable(GL_LINE_SMOOTH)
        else:
            glDisable(GL_LINE_SMOOTH)

        if self.properties["GL_POINT_SMOOTH"]==True:
            glEnable(GL_POINT_SMOOTH)
        else:
            glDisable(GL_POINT_SMOOTH)

        if self.properties["GL_POLYGON_SMOOTH"]==True:
            glEnable(GL_POLYGON_SMOOTH)
        else:
            glDisable(GL_POLYGON_SMOOTH)


        ## ALPHA BLENDING
        if self.properties["GL_BLEND"]==True:
            glEnable(GL_BLEND)
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
            
            ## FOG
            if self.properties["GL_FOG"]==True:
                glEnable(GL_FOG)
                glFogf(GL_FOG_MODE,    GL_LINEAR)
                glFogf(GL_FOG_START,   self.properties["GL_FOG_START"])
                glFogf(GL_FOG_END,     self.properties["GL_FOG_END"])
            else:
                glDisable(GL_FOG)
            
        else:
            glDisable(GL_BLEND)

    def glv_render(self):
        """Draw all GLDrawList objects onto the given glcontext/gldrawable.
        If the GLDrawList objects are not yet compiled into OpenGL draw
        lists, they will be compiled while they are drawn, since this is
        a useful optimization.
        """
        ## setup vieweport
        width  = self.properties["width"]
        height = self.properties["height"]
        glViewport(0, 0, width, height)

        ## setup perspective matrix
	glMatrixMode(GL_PROJECTION)
 	glLoadIdentity()

        zoom  = self.properties["zoom"]/2.0
        ratio = float(height)/float(width)

        glOrtho(-zoom,       zoom,
                -ratio*zoom, ratio*zoom,
                self.properties["near"],
                self.properties["far"])
        
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()

        self.glv_pre_render()
        
        ## inital orientation
        R   = self.properties["R"]
        cor = self.properties["cor"]

        ## sphere at 0,0,0
        self.glr_sphere((0.0, 0.0, 0.0), 1.0, 20)
        
        glMultMatrixf(
            (R[0,0],           R[1,0],      R[2,0], 0.0,
             R[0,1],           R[1,1],      R[2,1], 0.0,
             R[0,2],           R[1,2],      R[2,2], 0.0,
                0.0,              0.0,         0.0, 1.0) )

        glTranslatef(*-cor)

        ## render solid objects
        for draw_list in self.glo_iter_children():
            draw_list.gldl_render()

        ## render transparent objects
        for draw_list in self.glo_iter_children():
            draw_list.gldl_render(True) 

