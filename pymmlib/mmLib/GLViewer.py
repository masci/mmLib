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
        if prop_desc==None:
            raise "glo_link_child_property", \
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


class GLColor(object):
    """This is going to replace the 3-tuple RGB color...
    """
    def __init__(self, **args):
        self.rgb = args.get("rgb")


class OpenGLRenderMethods(object):
    """OpenGL renderer methods.  Eventually, all OpenGL rendering will
    be done through methods in this object.
    """
    def glr_set_material_rgb(self, r, g, b, a):
        """Creates a stock rendering material colored according to the given
        RGB values.
        """
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, (r, g, b, a))
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, (1.0, 1.0, 1.0, 1.0))
        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, (0.0, 0.0, 0.0, 1.0))
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
    
    def glr_axis(self, position, axis, radius):
        """Draw a vector axis using the current set material at position
        with the given radius.
        """
        if allclose(length(axis), 0.0):
            return

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

    def glr_Uaxes(self, position, U, prob, color, line_width):
        """Draw the anisotropic axies of the atom at the given probability.
        """
        C = GAUSS3C[prob]        
        eval, evec = eigenvectors(U)

        try:
            v0_peak = C * math.sqrt(eval[0])
        except ValueError:
            v0_peak = 0.0
        try:
            v1_peak = C * math.sqrt(eval[1])
        except ValueError:
            v1_peak = 0.0
        try:
            v2_peak = C * math.sqrt(eval[2])
        except ValueError:
            v2_peak = 0.0
        
        v0 = evec[0] * v0_peak
        v1 = evec[1] * v1_peak
        v2 = evec[2] * v2_peak

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
        C2 = GAUSS3C[prob]**2
        Ui = inverse(U)
        
        def ellipse_surface_func(x):
            d = matrixmultiply(x, matrixmultiply(Ui, x))
            try:
                vec = math.sqrt(C2 / d) * x
                return vec

            except ValueError:
                print "ellipse_surface_func: U not positive definate"
                return zeros(3, Float)

        self.glr_radial_surface(position, ellipse_surface_func, 3)

    def glr_Urms(self, position, U):
        """Renders the root mean square (one standard deviation) surface of the
        gaussian variance-covariance matrix U at the given position.  This
        is a peanut-shaped surface. (Note: reference the peanut paper!)
        """
        def rms_surface_func(v):
            return v * math.sqrt(abs(matrixmultiply(v, matrixmultiply(U, v))))
        self.glr_radial_surface(position, rms_surface_func, 3)

    def glr_radial_surface(self, position, func, depth):
        """
        """
        delta = 1.0e-7
        
        glPushMatrix()
        glTranslatef(*position)
        
        glEnable(GL_LIGHTING)
        glEnable(GL_NORMALIZE) ## let OpenGL normalize the normals
        glBegin(GL_TRIANGLES)
        
        for v1, v2, v3 in self.glr_iter_triangle_tesselation_cache(depth):
            u1 = func(v1)
            u2 = func(v2)
            u3 = func(v3)

            ## normals are computed numericly
            v21 = delta * (v2 - v1)
            v31 = delta * (v3 - v1)
            v32 = delta * (v3 - v2)

            d21 = func(normalize(v1 + v21))
            d31 = func(normalize(v1 + v31))
            n1  = cross(d21-u1, d31-u1)

            d32 = func(normalize(v2 + v32))
            d12 = func(normalize(v2 - v21))
            n2  = cross(d32-u2, d12-u2)

            d13 = func(normalize(v3 - v31))
            d23 = func(normalize(v3 - v32))
            n3  = cross(d13-u3, d23-u3)
                
            glNormal3f(*n1)
            glVertex3f(*u1)

            glNormal3f(*n2)
            glVertex3f(*u2)

            glNormal3f(*n3)
            glVertex3f(*u3)

        glEnd()
        glDisable(GL_NORMALIZE)
        glPopMatrix()

    def glr_iter_triangle_tesselation_cache(self, depth):
        """Caches expensive tesselation vericies.
        """
        if hasattr(self, "glr_tess_cache")==False:
            self.glr_tess_cache = {}

        try:

            return iter(self.glr_tess_cache[depth])
        except KeyError:
            print "createing cache"
            self.glr_tess_cache[depth] = list(
                self.glr_iter_triangle_tesselation(depth))
            return self.glr_tess_cache[depth]

    def glr_iter_triangle_tesselation(self, depth):
        """Iterates 3-tuples of unit vector length triangle vericies.  The
        triangle verticies are all yielded in counter-clockwise order, and
        the depth specifies the recursion depth of the algorithm.  
        """

        def tesselate(rdepth, v1, v2, v3):
            ## step 1: bisect the triangle edges
            v12 = (0.5 * (v2 - v1)) + v1
            v23 = (0.5 * (v3 - v2)) + v2
            v31 = (0.5 * (v1 - v3)) + v3

            ## step 2: form 4 triangles and either yield or recuse

            ## if we are at the given depth, then yield the triangle
            ## vertex
            rdepth -= 1

            if rdepth==0:
                yield normalize(v1),  normalize(v12), normalize(v31)
                yield normalize(v2),  normalize(v23), normalize(v12)
                yield normalize(v3),  normalize(v31), normalize(v23)
                yield normalize(v12), normalize(v23), normalize(v31)
            else:
                for tri1, tri2, tri3 in tesselate(rdepth, v1, v12, v31):
                    yield tri1, tri2, tri3
                for tri1, tri2, tri3 in tesselate(rdepth, v2, v23, v12):
                    yield tri1, tri2, tri3
                for tri1, tri2, tri3 in tesselate(rdepth, v3, v31, v23):
                    yield tri1, tri2, tri3       
                for tri1, tri2, tri3 in tesselate(rdepth, v12, v23, v31):
                    yield tri1, tri2, tri3
        
        x = array([1.0, 0.0, 0.0])
        y = array([0.0, 1.0, 0.0])
        z = array([0.0, 0.0, 1.0])
        
        for tri1, tri2, tri3 in tesselate(depth, x, y, z):
            yield tri1, tri2, tri3
        for tri1, tri2, tri3 in tesselate(depth, x, -z, y):
            yield tri1, tri2, tri3
        for tri1, tri2, tri3 in tesselate(depth, x, z, -y):
            yield tri1, tri2, tri3
        for tri1, tri2, tri3 in tesselate(depth, x, -y, -z):
            yield tri1, tri2, tri3
        for tri1, tri2, tri3 in tesselate(depth, -x, z, y):
            yield tri1, tri2, tri3
        for tri1, tri2, tri3 in tesselate(depth, -x, y, -z):
            yield tri1, tri2, tri3
        for tri1, tri2, tri3 in tesselate(depth, -x, -y, z):
            yield tri1, tri2, tri3
        for tri1, tri2, tri3 in tesselate(depth, -x, -z, -y):
            yield tri1, tri2, tri3


class GLDrawList(GLObject, OpenGLRenderMethods):
    """Fundamental OpenGL rigid entity.
    """
    def __init__(self, **args):
        GLObject.__init__(self, **args)
        self.__draw_list_id = None
        self.__transparent_draw_list_id = None
        
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
        
        ## render this GLDrawList, compiling a new OpenGL draw list if needed
        if self.__draw_list_id==None:
            self.__draw_list_id = glGenLists(1)
            glNewList(self.__draw_list_id, GL_COMPILE)
            self.gldl_draw()
            glEndList()

            self.__transparent_draw_list_id = glGenLists(1)
            glNewList(self.__transparent_draw_list_id, GL_COMPILE)
            self.gldl_draw_transparent()
            glEndList()

        self.gldl_push_matrix()

        ## support multiple rendering images by implementing class
        ## iterators gldl_iter_multidraw_all() for multiple
        ## rendering iterations of the GLDrawList and all its children,
        ## or gldl_iter_multidraw_self() for multiple images of just
        ## this GLDrawList, rendering the children just once
        for draw_flag_multi in self.gldl_iter_multidraw_all():

            for draw_flag_self in self.gldl_iter_multidraw_self():
                if transparent==True:
                    glCallList(self.__transparent_draw_list_id)
                else:
                    glCallList(self.__draw_list_id)

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
        if self.__draw_list_id!=None:
            glDeleteLists(self.__draw_list_id, 1)
            glDeleteLists(self.__transparent_draw_list_id, 1)
            
            self.__draw_list_id = None
            self.__transparent_draw_list_id = None
            
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

    def gldl_draw(self):
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
        GLDrawList.__init__(self, **args)
        self.unit_cell = args["unit_cell"]
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
        self.el_color_cache = {}
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
              "action":    "recompile" })
        self.glo_add_property(
            { "name":      "side_chain_visible",
              "desc":      "Show Side Chain Atoms",
              "type":      "boolean",
              "default":   True,
              "action":    "recompile" })
        self.glo_add_property(
            { "name":      "hetatm_visible",
              "desc":      "Show Hetrogen Atoms",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   True,
              "action":    "recompile" })
        self.glo_add_property(
            { "name":      "water_visible",
              "desc":      "Show Waters",
              "catagory":  "Show/Hide",
              "type":      "boolean",
              "default":   False,
              "action":    "redraw" })

        ## lines
        self.glo_add_property(
            { "name":       "lines",
              "desc":       "Draw Atom Bond Lines",
              "catagory":   "Show/Hide",
              "type":       "boolean",
              "default":    False,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "line_width",
              "desc":       "Bond Line Drawing Width",
              "catagory":   "Bond Lines",
              "type":       "float",
              "spin":       PROP_LINE_RANGE,
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
            { "name":       "cpk_opacity_occupancy",
              "desc":       "Set Opacity by Atom Occupancy",
              "catagory":   "CPK",
              "type":       "boolean",
              "default":    False,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "cpk_scale_radius",
              "desc":       "Scale CPK Radius",
              "catagory":   "CPK",
              "type":       "float",
              "spin":       "0.0-5.0,0.1",
              "default":    1.0,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "cpk_opacity",
              "desc":       "CPK Sphere Opacity",
              "catagory":   "CPK",
              "type":       "float",
              "range":      PROP_OPACITY_RANGE,
              "default":    1.00,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "sphere_quality",
              "desc":       "CPK Sphere Quality",
              "catagory":   "CPK",
              "type":       "integer",
              "range":      "5-36,5",
              "default":    10,
              "action":     "recompile" })

        ## trace           
        self.glo_add_property(
            { "name":       "trace",
              "desc":       "Draw Backbone Trace",
              "catagory":   "Show/Hide",
              "type":       "boolean",
              "default":    True,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "trace_line_width",
              "desc":       "Backbone Trace Line Width",
              "catagory":   "Trace",
              "type":       "float",
              "spin":       PROP_LINE_RANGE,
              "default":    4.0,
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
            { "name":       "adp_prob",
              "desc":       "Contour Probability",
              "catagory":   "ADP",
              "type":       "integer",
              "range":      PROP_PROBABILTY_RANGE,
              "default":    50,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":      "U",
              "desc":      "Show Thermal Axes",
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
        self.glo_add_property(
            { "name":       "ellipse",
              "desc":       "Show Thermal Ellipseoids",
              "catagory":   "Show/Hide",
              "type":       "boolean",
              "default":    False,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "ellipse_opacity",
              "desc":       "Thermal Ellipseoid Opacity",
              "catagory":   "ADP",
              "type":       "float",
              "range":      PROP_OPACITY_RANGE,
              "default":    1.0,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "rms",
              "desc":       "Show Thermal Peanuts (RMS Deviation Surface)",
              "catagory":   "Show/Hide",
              "type":       "boolean",
              "default":    False,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "rms_opacity",
              "desc":       "Peanut Surface Opacity",
              "catagory":   "ADP",
              "type":       "float",
              "range":      PROP_OPACITY_RANGE,
              "default":    1.0,
              "action":     "recompile" })

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
            yield chain

    def gldl_draw(self):
        """Draw all selected representations.
        """
        #self.draw_symmetry_debug()
        
        ## drawing functions requiring neighboring atoms
        if self.properties["trace"]==True:
            self.draw_trace()
            
        ## per-atom draw functions
        draw_funcs = []
        
        if self.properties["lines"]==True:
            draw_funcs.append(self.draw_lines)

        if self.properties["U"]==True:
            draw_funcs.append(self.draw_U_axes)

        if self.properties["cpk"]==True and \
           self.properties["cpk_opacity"]>=1.0:
            draw_funcs.append(self.draw_cpk)

        if self.properties["ellipse"]==True and \
           self.properties["ellipse_opacity"]>=1.0:
            draw_funcs.append(self.draw_ellipse)

        if self.properties["rms"]==True and \
           self.properties["rms_opacity"]>=1.0:
            draw_funcs.append(self.draw_rms)

        ## execute draw functions for each atom
        for atm in self.glal_iter_atoms():
            for draw_func in draw_funcs:
                draw_func(atm)

    def gldl_draw_transparent(self):
        """
        """
        draw_funcs = []

        if self.properties["cpk"]==True and \
           self.properties["cpk_opacity"]<1.0:
            draw_funcs.append(self.draw_cpk)
            
        if self.properties["ellipse"]==True and \
           self.properties["ellipse_opacity"]<1.0:
            draw_funcs.append(self.draw_ellipse)

        if self.properties["rms"]==True and \
           self.properties["rms_opacity"]<1.0:
            draw_funcs.append(self.draw_rms)

        for atm in self.glal_iter_atoms():
            for draw_func in draw_funcs:
                draw_func(atm)
    
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

        if self.properties["cpk_opacity_occupancy"]==True:
            a = atm.occupancy
            self.glr_set_material_rgb(r, g, b, a)

            for atm_alt in atm.iter_alt_loc():
                self.glr_cpk(
                    self.calc_position(atm_alt.position),
                    self.properties["cpk_scale_radius"] * radius,
                    self.properties["sphere_quality"])
        else:
            a = self.properties["cpk_opacity"]
            self.glr_set_material_rgb(r, g, b, a)

            self.glr_cpk(
                self.calc_position(atm.position),
                self.properties["cpk_scale_radius"] * radius,
                self.properties["sphere_quality"])

    def draw_ellipse(self, atm):
        """Draw the ADP determined probability ellipseoid.
        """
        r, g, b = self.get_color(atm)        
        self.glr_set_material_rgb(r, g, b, self.properties["ellipse_opacity"])
        self.glr_Uellipse(
            self.calc_position(atm.position),
            atm.get_U(),
            self.properties["adp_prob"])

    def draw_rms(self, atm):
        """Draw the ADP determined RMS displacement surface.
        """
        r, g, b = self.get_color(atm)        
        self.glr_set_material_rgb(r, g, b, self.properties["rms_opacity"])
        self.glr_Urms(self.calc_position(atm.position), atm.get_U())

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

        for chain in self.glal_iter_chains():
            glBegin(GL_LINE_STRIP)

            for frag in chain.iter_fragments():
                ## PROTEIN
                if frag.is_amino_acid()==True:
                    for name in ("N", "CA", "C"):
                        try:
                            atm = frag[name]
                        except KeyError:
                            pass
                        else:
                            glVertex3f(*self.calc_position(atm.position))
                ## DNA/RNA
                elif frag.is_nucleic_acid()==True:
                    pass

            glEnd()

    def draw_U_axes(self, atm):
        """Draw the anisotropic axies of the atom with R.M.S.
        magnitude.
        """
        U = getattr(atm, self.properties["atm_U_attr"], None)
        if U==None:
            return
        
        self.glr_Uaxes(
            self.calc_position(atm.position),
            U,
            self.properties["adp_prob"],
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
            { "name":        "Lx_prob",
              "desc":        "Surface Probability Range",
              "catagory":    "TLS",
              "type":        "integer",
              "range":       PROP_PROBABILTY_RANGE,
              "default":     50,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "surface_opacity",
              "desc":        "Screw Surface Opacity",
              "catagory":    "TLS",
              "type":        "float",
              "range":      PROP_OPACITY_RANGE,
              "default":     1.0,
              "action":      "recompile" })

        self.glo_add_property(
            { "name":        "COR",
              "desc":        "TLS Center of Reaction", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })        
        self.glo_add_property(
            { "name":        "Lx_eigen_val",
              "desc":        "Libration Axis Eigenvalue", 
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "Lx_eigen_vec",
              "desc":        "Libration Axis Normalized Eigenvector",  
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "Lx_rho",
              "desc":        "Libration Axis Translation from COR",  
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "array(3)",
              "default":     zeros(3, Float),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "Lx_pitch",
              "desc":        "Libration screw pitch (DEG/A)",  
              "catagory":    "TLS Analysis",
              "read_only":   True,
              "type":        "float",
              "default":     0.0,
              "action":      "recompile" })

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
        if self.properties["surface_opacity"]==1.0:
            self.draw_TLS_surface(
                self.properties["Lx_eigen_vec"],
                self.properties["Lx_eigen_val"],
                self.properties["Lx_rho"],
                self.properties["Lx_pitch"])
    
    def gldl_draw_transparent(self):
        if self.properties["surface_opacity"]<1.0:
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

        C = GAUSS3C[self.properties["Lx_prob"]]
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
        glEnable(GL_NORMALIZE)

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

                for bond in self.bond_list:

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

        glDisable(GL_NORMALIZE)
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE)

        
class GLTLSAtomList(GLAtomList):
    """OpenGL visualizations of TLS group atoms.
    """
    def __init__(self, **args):
        GLAtomList.__init__(self, **args)
        self.tls_group = args["tls_group"]        
        self.glo_init_properties(**args)
    
    def glo_install_properties(self):
        GLAtomList.glo_install_properties(self)

        ## Show/Hide
        self.glo_add_property(
            { "name":        "CA_line_visible",
              "desc":        "Show Lines from COR to Backbone",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "fan_visible",
              "desc":        "Show Fans from COR to Backbone",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     True,
              "action":      "recompile" })

        self.glo_add_property(
            { "name":        "Ut_axes",
              "desc":        "Show T Thermal Axes",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "Ut_ellipse",
              "desc":        "Show T Thermal Ellipsoids",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "Ut_rms",
              "desc":        "Show T Thermal Peanuts",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
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
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "fan_opacity",
              "desc":        "COR/Backbone Fan Opacity",
              "catagory":    "TLS",
              "type":        "float",
              "range":       PROP_OPACITY_RANGE,
              "default":     0.25,
              "action":      "recompile" })

        self.glo_add_property(
            { "name":        "Ut_color",
              "desc":        "T Color",
              "catagory":    "TLS",
              "type":        "color",
              "default":     (0.5, 0.5, 1.0),
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "Ut_ellipse_opacity",
              "desc":        "T Thermal Ellipse Opacity",
              "catagory":    "TLS",
              "type":        "float",
              "range":       PROP_OPACITY_RANGE,
              "default":     1.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "Ut_rms_opacity",
              "desc":        "T Thermal RMS Peanut Opacity",
              "catagory":    "TLS",
              "type":        "float",
              "range":       PROP_OPACITY_RANGE,
              "default":     1.0,
              "action":      "recompile" })
        
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

    def gldl_update_cb(self, updates, actions):
        GLAtomList.gldl_update_cb(self, updates, actions)
        if "trace" in updates:
            self.properties.update(lines=False)

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

    def gldl_draw(self):
        GLAtomList.gldl_draw(self)

        if self.properties["Ut_axes"]:
            self.draw_Ut_axes()
            
        ## draw a line from the COR (center of reaction)
        ## to all CA atoms in the TLS group
        if self.properties["CA_line_visible"]==True:
            self.draw_CA_lines()

        ## draw a transparent fan from the COR to all backbone atoms
        ## in the TLS group
        if self.properties["fan_visible"]==True and \
           self.properties["fan_opacity"]>=1.0:
            self.draw_fan()

        if self.properties["Ut_ellipse"]==True and \
           self.properties["Ut_ellipse_opacity"]>=1.0:
            self.draw_Ut_ellipse()

        if self.properties["Ut_rms"]==True and \
           self.properties["Ut_rms_opacity"]>=1.0:
            self.draw_Ut_rms()

    def gldl_draw_transparent(self):
        GLAtomList.gldl_draw_transparent(self)

        if self.properties["fan_visible"]==True and \
           self.properties["fan_opacity"]<1.0:
            self.draw_fan()

        if self.properties["Ut_ellipse"]==True and \
           self.properties["Ut_ellipse_opacity"]<1.0:
            self.draw_Ut_ellipse()

        if self.properties["Ut_rms"]==True and \
           self.properties["Ut_rms_opacity"]<1.0:
            self.draw_Ut_rms()

    def draw_Ut_axes(self):
        """Draw the anisotropic axies of the T tensor only.
        """
        T     = self.properties["T"]
        color = self.properties["Ut_color"]

        for atm in self.tls_group:
            self.glr_Uaxes(
                self.calc_position(atm.position),
                T,
                self.properties["Ut_prob"],
                color,
                1.0)

    def draw_Ut_ellipse(self):
        """Draw the ADP determined probability ellipseoid of the T tensor only.
        """
        T       = self.properties["T"]
        r, g, b = self.properties["Ut_color"]
        a       = self.properties["Ut_ellipse_opacity"]
        
        self.glr_set_material_rgb(r, g, b, a)

        for atm in self.tls_group:
            self.glr_Uellipse(
                self.calc_position(atm.position),
                T,
                self.properties["Ut_prob"])

    def draw_Ut_rms(self, atm):
        """Draw the ADP determined RMS displacement surface of the T tensor
        only.
        """
        T       = self.properties["T"]
        r, g, b = self.properties["Ut_color"]
        a       = self.properties["Ut_rms_opacity"]
        
        self.glr_set_material_rgb(r, g, b, a)

        for atm in self.atom_list:
            self.glr_Urms(self.calc_position(atm.position), T)
   
    def draw_CA_lines(self):
        glDisable(GL_LIGHTING)        
        COR = self.properties["COR"]
        glColor3f(*self.properties["fan_color"])
        glLineWidth(1.0)
        
        for atm in self.tls_group:
            if atm.name=="CA":
                glBegin(GL_LINES)
                glVertex3f(0.0, 0.0, 0.0)
                glVertex3f(*atm.position - COR)
                glEnd()

    def draw_fan(self):
        glEnable(GL_LIGHTING)

        COR     = self.properties["COR"]
        r, g, b = self.properties["fan_color"]
        a       = self.properties["fan_opacity"]

        self.glr_set_material_rgb(r, g, b, a)
        
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE)
        glEnable(GL_NORMALIZE) ## let OpenGL normalize the normals

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
        glDisable(GL_NORMALIZE)


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

        ## step 3: add visualization objects for each TLS screw rotation axis
        self.gl_screw_rot1 = GLTLSScrewRotation(atom_list = self.tls_group)
        self.gl_screw_rot2 = GLTLSScrewRotation(atom_list = self.tls_group)
        self.gl_screw_rot3 = GLTLSScrewRotation(atom_list = self.tls_group)

        self.gl_screw_rot1.glo_set_name("L1 Screw Rotation")
        self.gl_screw_rot2.glo_set_name("L2 Screw Rotation")
        self.gl_screw_rot3.glo_set_name("L3 Screw Rotation")
        
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

        self.glo_link_child_property(
            "adp_prob", "gl_screw_rot1", "Lx_prob")
        self.glo_link_child_property(
            "adp_prob", "gl_screw_rot2", "Lx_prob")
        self.glo_link_child_property(
            "adp_prob", "gl_screw_rot3", "Lx_prob")

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
              "action":      "recompile" })
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
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "U",
              "desc":        "Show TLS Thermal Axes",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "ellipse",
              "desc":        "Show TLS Thermal Ellipsoids",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "rms",
              "desc":        "Show TLS Thermal Peanuts",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "animation_visible",
              "desc":        "Show TLS Animated Atoms", 
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     True,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "screw1_visible",
              "desc":        "Show Screw Axis 1 Traced Surface", 
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "screw2_visible",
              "desc":        "Show Screw Axis 2 Traced Surface", 
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "screw3_visible",
              "desc":        "Show Screw Axis 3 Traced Surface",
              "catagory":    "Show/Hide",
              "type":        "boolean",
              "default":     False,
              "action":      "redraw" })
        
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
              "action":     "recompile" })
        self.glo_add_property(
            { "name":        "T_color",
              "desc":        "T Axes Color",
              "catagory":    "TLS",
              "type":        "color",
              "default":     (0.5, 0.5, 1.0),
              "action":      "recompile" })

        self.glo_add_property(
            { "name":       "L_axis_radius",
              "desc":       "Screw Axes Radius",
              "catagory":   "TLS",
              "type":       "float",
              "default":    0.05,
              "action":     "recompile" })
        self.glo_add_property(
            { "name":       "L_color",
              "desc":       "Screw Axes Color",
              "catagory":   "TLS",
              "type":       "color",
              "default":    (1.0, 0.5, 0.5),
              "action":     "recompile" })

        self.glo_add_property(
            { "name":        "Utls_color",
              "desc":        "Thermal Utls Color",
              "catagory":    "TLS",
              "type":        "color",
              "default":     (0.10, 0.75, 0.35),
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "U_color",
              "desc":        "Thermal Atom U Color for Diff Mode",
              "catagory":    "TLS",
              "type":        "color",
              "default":     (1.0, 1.0, 1.0),
              "action":      "redraw" })

        self.glo_add_property(
            { "name":        "ellipse_opacity",
              "desc":        "TLS Thermal Ellipseoid Opacity",
              "catagory":    "TLS",
              "type":        "float",
              "range":       PROP_OPACITY_RANGE,
              "default":     1.0,
              "action":      "recompile" })
        self.glo_add_property(
            { "name":        "rms_opacity",
              "desc":        "TLS Thermal Peanut Opacity",
              "catagory":    "TLS",
              "type":        "float",
              "range":       PROP_OPACITY_RANGE,
              "default":     1.0,
              "action":      "recompile" })

         
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

    def gldl_draw(self):
        draw_func_list = []

        if self.properties["U"]==True:
            draw_func_list.append(self.draw_Utls_axes)
        
        if self.properties["ellipse"]==True and \
           self.properties["ellipse_opacity"]>=1.0:
            draw_func_list.append(self.draw_Utls_ellipse)

        if self.properties["rms"]==True and \
           self.properties["rms_opacity"]>=1.0:
            draw_func_list.append(self.draw_Utls_rms)

        ## everything is drawn from the center of reaction
        COR = self.properties["COR"]

        glPushMatrix()
        glTranslatef(*COR)

        ## draw TLS axes
        if self.properties["TLS_visible"]==True:
            self.draw_tensors()

        if len(draw_func_list)>0:
            for atm, Utls in self.tls_group.iter_atm_Utls():
                if self.properties["main_chain_only"]==True and \
                    atm.name not in ("C", "CA", "N"):
                    continue
                for draw_func in draw_func_list:
                    draw_func(atm, Utls)

        glPopMatrix()

    def gldl_draw_transparent(self):
        draw_func_list = []
        
        if self.properties["ellipse_opacity"]<1.0:
            draw_func_list.append(self.draw_Utls_ellipse)

        if self.properties["rms_opacity"]<1.0:
            draw_func_list.append(self.draw_Utls_rms)

        if len(draw_func_list)==0:
            return

        ## everything is drawn from the center of reaction
        COR = self.properties["COR"]

        glPushMatrix()
        glTranslatef(*COR)

        for atm, Utls in self.tls_group.iter_atm_Utls():
            if self.properties["main_chain_only"]==True and \
               atm.name not in ("C", "CA", "N"):
                continue
            for draw_func in draw_func_list:
                draw_func(atm, Utls)
            
        glPopMatrix()
    
    def draw_Utls_axes(self, atm, Utls):
        """Render the anisotropic thremal axes calculated from the TLS
        model.
        """
        color = self.properties["Utls_color"]
        self.glr_Uaxes(atm.position - self.properties["COR"],
                       Utls,
                       self.properties["adp_prob"],
                       color,
                       1.0)

    def draw_Utls_ellipse(self, atm, Utls):
        """Render the anisotropic thremal ellipseoids/RMS surface from
        TLS calculated U values.
        """
        r, g, b = self.properties["Utls_color"]
        a       = self.properties["ellipse_opacity"]
        self.glr_set_material_rgb(r, g, b, a)

        if self.properties["add_biso"]==True:
            Uiso = identity(3, Float) * (atm.temp_factor / (24.0*math.pi**2))
            Utls = Utls + Uiso

        self.glr_Uellipse(
            atm.position - self.properties["COR"],
            Utls,
            self.properties["adp_prob"])

        if self.properties["Utls_U_diff"]:
            r, g, b = self.properties["U_color"]
            self.glr_set_material_rgb(r, g, b, a)
            self.glr_Uellipse(
                atm.position - self.properties["COR"],
                atm.get_U(),
                self.properties["adp_prob"])

    def draw_Utls_rms(self, atm, Utls):
        """Render the anisotropic thremal ellipseoids/RMS surface from
        TLS calculated U values.
        """
        r, g, b = self.properties["Utls_color"]
        a       = self.properties["rms_opacity"]
        self.glr_set_material_rgb(r, g, b, a)

        self.glr_Urms(atm.position - self.properties["COR"], Utls)

    def draw_tensors(self):
        """Draw tensor axis.
        """
        ## put a small sphere at the center of reaction
        r, g, b = self.properties["Utls_color"]        
        self.glr_set_material_rgb(r, g, b, 1.0)
        self.glr_cpk(array([0.0, 0.0, 0.0]), 0.05, 12)
        
        ## T: units (A^2)
        self.glr_Uaxes(
            (0.0, 0.0, 0.0),
            self.properties["T_reduced"],
            self.properties["adp_prob"],
            self.properties["T_color"],
            self.properties["T_line_width"])

        ## L: units (DEG^2)
        for Lx_eigen_val, Lx_eigen_vec, Lx_rho in [
            ("L1_eigen_val", "L1_eigen_vec", "L1_rho"),
            ("L2_eigen_val", "L2_eigen_vec", "L2_rho"),
            ("L3_eigen_val", "L3_eigen_vec", "L3_rho")]:

            Lx_eigen_val = self.properties[Lx_eigen_val]

            if Lx_eigen_val<=0.0:
                continue

            v  = self.properties[Lx_eigen_vec].copy() *\
                 GAUSS3C[self.properties["adp_prob"]] *\
                 math.sqrt(Lx_eigen_val)

            rho = self.properties[Lx_rho]

            ## line from COR to center of screw/rotation axis
            glColor3f(*self.properties["L_color"])

            glDisable(GL_LIGHTING)
            glBegin(GL_LINES)
            glVertex3f(0.0, 0.0, 0.0)
            glVertex3f(*rho)
            glEnd()

            r, g, b = self.properties["L_color"]

            ## for negitive librations, make the L axis really dark
            if Lx_eigen_val<0.0:
                self.glr_set_material_rgb(r*0.25, g*0.25, b*0.25, 1.0)
            else:
                self.glr_set_material_rgb(r, g, b, 1.0)
            
            self.glr_axis(rho-v, 2*v, self.properties["L_axis_radius"])


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

        for frag in self.chain.iter_fragments():
            gl_frag = GLFragment(fragment=frag)
            self.glo_add_child(gl_frag)

        self.glo_init_properties(**args)

    def glo_install_properties(self):
        GLAtomList.glo_install_properties(self)

    def glo_name(self):
        return "Chain %s" % (self.chain.chain_id)

    def glal_iter_atoms(self):
        for frag in self.chain.iter_fragments():
            for atm in frag.iter_all_atoms():
                yield atm

##     def glal_iter_fragments(self):
##         for frag in self.chain.iter_fragments():
##             yield frag


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
              "type":     "boolean",
              "default":  True,
              "action":  "redraw" })
        self.glo_add_property(
            { "name":     "unit_cell_visible",
              "desc":     "Show Unit Cell",
              "type":     "boolean",
              "default":  True,
              "action":   "redraw" })        

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
              "read_only":True,
              "type":     "array(3,3)",
              "default":  identity(3, Float),
              "action":   "redraw" })
        self.glo_add_property(
            { "name":     "t",
              "desc":     "View Window Translation Vector",
              "read_only":True,
              "type":     "array(3)",
              "default":  array([0.0, 0.0, -50]),
              "action":   "redraw" })

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
              "default":   True,
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
        pass

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

        ## OpenGL Features
        #glEnable(GL_AUTO_NORMAL)
        #glEnable(GL_NORMALIZE)
        glEnable(GL_DEPTH_TEST)
        glDepthFunc(GL_LESS)

        ## inital orientation
        R = self.properties["R"]
        t = self.properties["t"]
        
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
	glLoadIdentity()
        
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
        glLightfv(GL_LIGHT0, GL_POSITION, (0.0, 0.0, 1.0, 0.0))

        ## orientation
        glTranslatef(0.0, 0.0, -10.0)

        glMultMatrixf(
            (R[0,0],           R[1,0],      R[2,0], 0.0,
             R[0,1],           R[1,1],      R[2,1], 0.0,
             R[0,2],           R[1,2],      R[2,2], 0.0,
                0.0,              0.0,         0.0, 1.0) )

        glTranslatef(t[0], t[1], t[2] + 10.0)

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

        ## render solid objects
        for draw_list in self.glo_iter_children():
            draw_list.gldl_render()

        ## render transparent objects
        for draw_list in self.glo_iter_children():
            draw_list.gldl_render(True) 
