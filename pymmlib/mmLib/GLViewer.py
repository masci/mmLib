## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""OpenGL rendering classes.
"""
from __future__  import generators
import weakref

from OpenGL.GL   import *
from OpenGL.GLU  import *
from OpenGL.GLUT import *
from mmTypes     import *
from Structure   import *


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
            for child2 in child1.iter_breath_first():
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
            ancestor = ancestor.__parent
        return False

    def glo_add_child(self, child):
        assert isinstance(child, GLObject)
        assert child.__globject_parent==None
        child.__globject_parent = self
        self.__globject_children.append(child)
        
    def glo_remove_child(self, child):
        """Removes child GLObject.
        """
        assert isinstance(child, GLObject)
        assert child.__globject_parent==self
        child.__globject_parent = None
        self.__globject_children.remove(child)

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
            else:
                print "## glo_init_properties(%s) %s=%s" % (
                    str(self.__globject_properties_id), name, str(self.properties[name]))

            try:
                linked_props = prop_desc["link"]
            except KeyError:
                continue

            for linked_prop in linked_props:
                gl_object_id = linked_prop["gl_object"]
                child_name   = linked_prop["name"]
                gl_object    = self.glo_get_child(gl_object_id)
                child_args   = {child_name: self.properties[name]}
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
                print "## glo_update_properties(%s) %s=%s" % (
                    str(self.__globject_properties_id), name, str(self.properties[name]))
                
                updates[name] = self.properties[name]
                if prop_desc["action"] not in actions:
                    actions.append(prop_desc["action"])

            try:
                linked_props = prop_desc["link"]
            except KeyError:
                continue

            for linked_prop in linked_props:
                gl_object_id = linked_prop["gl_object"]
                child_name   = linked_prop["name"]
                gl_object    = self.glo_get_child(gl_object_id)
                child_args   = {child_name: self.properties[name]}
                gl_object.glo_update_properties(**child_args)

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
        glTranslatef(*self.properties["origin"])
        axes = self.properties["axes"]
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
        glCallList(self.gl_name)
        self.gl_pop_matrix()

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

        glCallList(self.gl_name)

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

    def gl_draw(self):
        self.draw_cell(-1, -1, -1, 0, 0, 0)


class GLAtomList(GLDrawList, AtomList):
    """OpenGL renderer for a list of atoms.  Optional arguments iare:
    color, U, U_color.
    """
    def __init__(self, **args):
        GLDrawList.__init__(self, **args)
        AtomList.__init__(self)
        self.el_color_cache = {}
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
            { "name":       "cpk",
              "desc":       "Draw CPK Spheres",
              "type":       "boolean",
              "default":    False,
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

    def gl_draw(self):
        """Perform the OpenGL drawing operations to render this atom list
        with the current settings.
        """
        if self.properties["lines"]==True:
            for atm in self:
                self.draw_lines(atm)

        if self.properties["cpk"]==True:
            for atm in self:
                self.draw_cpk(atm)

        if self.properties["U"]==True:
            for atm in self:
                self.draw_U_axes(atm)

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

        color = self.get_color(atm)
        self.set_material(r=color[0], g=color[1], b=color[2], br=0.75)

        glPushMatrix()
        if self.properties["atom_origin"]!=None:
            glTranslatef(*atm.position - self.properties["atom_origin"])
        else:
            glTranslatef(*atm.position)

        glEnable(GL_LIGHTING)

        sphere_quality = self.properties["sphere_quality"]
        glutSolidSphere(radius, sphere_quality, sphere_quality)
        glPopMatrix()

    def set_material(self, r=1.0, g=1.0, b=1.0, br=1.0, alpha=1.0):
        """Utility function used to set material for rendering CPK spheres.
        """
        ambient  = 0.5
        diffuse  = 0.5
        specular = 0.5

        amb_r = ambient * r * br
        amb_g = ambient * g * br
        amb_b = ambient * b * br

        dif_r = diffuse * br
        spe_r = specular * br

        shine = 50.0 * br
        
	glMaterial(GL_FRONT, GL_AMBIENT, [amb_r, amb_g, amb_b, alpha])
	glMaterial(GL_FRONT, GL_DIFFUSE, [dif_r, dif_r, dif_r, alpha])
	glMaterial(GL_FRONT, GL_SPECULAR,[spe_r, spe_r, spe_r, alpha])
	glMaterial(GL_FRONT, GL_SHININESS, shine)

    def draw_lines(self, atm):
        """Draw a atom using bond lines only.
        """
        glColor3f(*self.get_color(atm))

        if self.properties["atom_origin"]!=None:
            position = atm.position - self.properties["atom_origin"]
        else:
            position = atm.position
        
        glLineWidth(self.properties["line_width"])
        glDisable(GL_LIGHTING)
        
        ## if there are bonds, then draw the lines 1/2 way to the
        ## bonded atoms
        if len(atm.bond_list) > 0:
            for bond in atm.iter_bonds():
                atm2 = bond.get_partner(atm)

                start = position
                end   = start + ((atm2.position - atm.position) / 2)

                glBegin(GL_LINES)
                glVertex3f(*start)
                glVertex3f(*end)
                glEnd()

        ## if there are no bonds, draw a small cross-point 
        else:
            vx = array([0.25, 0.0,  0.0])
            vy = array([0.0,  0.25, 0.0])
            vz = array([0.0,  0.0,  0.25])

            start = position - vx
            end   = position + vx
            glBegin(GL_LINES)
            glVertex3f(*start)
            glVertex3f(*end)
            glEnd()

            start = position - vy
            end   = position + vy
            glBegin(GL_LINES)
            glVertex3f(*start)
            glVertex3f(*end)
            glEnd()

            start = position - vz
            end   = position + vz
            glBegin(GL_LINES)
            glVertex3f(*start)
            glVertex3f(*end)
            glEnd()

        glEnable(GL_LIGHTING)

    def draw_U_axes(self, atm):
        """Draw the anisotropic axies of the atom with root mean square
        """
        if atm.U == None:
            return
        
        if self.properties["U"] == True:
            U = atm.U
        else:
            U = getattr(atm, self.properties["atm_U_attr"])

        eval, evec = eigenvectors(U)

        v0_peak = 1.414 * math.sqrt(eval[0])
        v1_peak = 1.414 * math.sqrt(eval[1])
        v2_peak = 1.414 * math.sqrt(eval[2])
        
        v0 = evec[0] * v0_peak
        v1 = evec[1] * v1_peak
        v2 = evec[2] * v2_peak

        if self.properties["atom_origin"]:
            position = atm.position - self.properties["atom_origin"]
        else:
            position = atm.position

        glColor3f(*self.properties["U_color"])

        glLineWidth(1.0)
        glDisable(GL_LIGHTING)
        
        glBegin(GL_LINES)
        glVertex3f(*position - v0)
        glVertex3f(*position + v0)
        glEnd()

        glBegin(GL_LINES)
        glVertex3f(*position - v1)
        glVertex3f(*position + v1)
        glEnd()

        glBegin(GL_LINES)
        glVertex3f(*position - v2)
        glVertex3f(*position + v2)
        glEnd()

        glEnable(GL_LIGHTING)

    def draw_atom_old(self, atm, symm = True):
        ## draw symmetry equivelent positions
        if symm:
            
            uc = atm.get_structure().unit_cell
            orig_pos = atm.position
            frac_pos = uc.calc_orth_to_frac(atm.position)

            for fpos in uc.space_group.iter_equivalent_positions(frac_pos):
                (x, y, z) = fpos
                
                if x < 0.0: x += 1.0
                elif x >= 1.0: x -= 1.0

                if y < 0.0: y += 1.0
                elif y >= 1.0: y -= 1.0

                if z < 0.0: z += 1.0
                elif z >= 1.0: z -= 1.0

                for (tx, ty, tz) in [(0,0,0), (0,-1,0), (-1,0,0),(-1,-1,0)]:
                    v = array([x+tx, y+ty, z+tz])
                    atm.position = uc.calc_frac_to_orth(v)
                    self.draw_atom(atm, False)

            atm.position = orig_pos
            return
        

class GLTLSGroup(GLDrawListContainer):
    """Draws TLS group
    """
    def __init__(self, **args):
        GLDrawListContainer.__init__(self)

        ## TLS calculations 
        self.tls_group = args["tls_group"]
        self.calcs     = self.tls_group.calc_COR()        
        (eigen_values, eigen_vectors) = eigenvectors(self.tls_group.L)
        self.evalL     = eigen_values

        ## add child GLAtomList
        self.gl_atom_list = GLAtomList(
            U_atm_arg   = "Ucalc",
            atom_origin = self.calcs["COR"])

        for atm, Ucalc in self.tls_group.iter_atm_Ucalc():
            atm.Ucalc = Ucalc
            self.gl_atom_list.append(atm)

        self.gl_atom_list.glo_set_properties_id("gl_atom_list")
        self.glo_add_child(self.gl_atom_list)

        self.glo_link_child_property(
            "atom_color", "gl_atom_list", "color")
        self.glo_link_child_property(
            "atom_line_width", "gl_atom_list", "line_width")
        self.glo_link_child_property(
            "U", "gl_atom_list", "U")
        self.glo_link_child_property(
            "U_color", "gl_atom_list", "U_color")

        ## initalize properties
        self.glo_add_update_callback(self.tls_update_cb)
        self.glo_init_properties(
            origin = self.calcs["COR"],
            axes   = eigen_vectors,
            **args)

    def glo_install_properties(self):
        GLDrawListContainer.glo_install_properties(self)

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
            { "name":        "atom_color",
              "desc":        "Atom Color",
              "type":        "color",
              "default":     (1.0, 1.0, 1.0),
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "atom_line_width",
              "desc":        "Atom Line Width",
              "type":        "float",
              "default":     10.0,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "U",
              "desc":        "Show ADP Utls Axes",
              "type":        "boolean",
              "default":     False,
              "action":      "redraw" })
        self.glo_add_property(
            { "name":        "U_color",
              "desc":        "ADP Axes Color",
              "type":        "color",
              "default":     (0.0, 1.0, 0.0),
              "action":      "redraw" })
        
    def tls_update_cb(self, updates, actions):
        print "## tls update cb = %s" % (str(updates))
        if "time" in updates:
            print "## updating time"
            self.update_time()

    def iter_draw_lists(self):
        yield self.gl_atom_list

    def gl_draw(self):
        ## draw TLS axes
        if self.properties["TLS_visible"]==True:
            
            glColor3f(*self.properties["T_color"])
            self.draw_tensor(self.calcs["T'"])

            glColor3f(*self.properties["L_color"])
            self.draw_tensor(self.calcs["L'"], 0.05 * rad2deg2)

            glColor3f(*self.properties["S_color"])
            self.draw_tensor(self.calcs["S'"], rad2deg)

        ## draw a line from the center of reaction to all CA atoms in the
        ## TLS group
        if self.properties["CA_line_visible"]==True:
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

    def draw_tensor(self, ten, scale = 1.0):
        """Draw tensor axis.
        """
        (eval, evec) = eigenvectors(ten)

        glDisable(GL_LIGHTING)
        glLineWidth(self.properties["tensor_line_width"])
        
        for i in range(3):
            v = scale * eval[i] * array([evec[i,0],evec[i,1],evec[i,2]])
            glBegin(GL_LINES)
            glVertex3f(*-v)
            glVertex3f(*v)
            glEnd()

        glEnable(GL_LIGHTING)

    def update_time(self):
        """Changes the time of the TLS group simulating harmonic motion.
        """
        Lx_peak = 1.414 * math.sqrt(self.evalL[0]*rad2deg2)
        Ly_peak = 1.414 * math.sqrt(self.evalL[1]*rad2deg2)
        Lz_peak = 1.414 * math.sqrt(self.evalL[2]*rad2deg2)

        sin_tm = math.sin(self.properties["time"] * 2 * math.pi)

        Lx = Lx_peak * sin_tm 
        Ly = Ly_peak * sin_tm
        Lz = Lz_peak * sin_tm

##         Spc = self.calcs["S'^"] * rad2deg

##         dSp = array(
##             [ (Lx * Spc[0,0]) + (Ly * Spc[1,0]) + (Lz * Spc[2,0]),
##               (Lx * Spc[0,1]) + (Ly * Spc[1,1]) + (Lz * Spc[2,1]),
##               (Lx * Spc[0,2]) + (Ly * Spc[1,2]) + (Lz * Spc[2,2]) ])

##         dS = matrixmultiply(transpose(self.axes), dSp)

        dS = zeros(3)
        origin = array(self.calcs["COR"] + dS)

        self.glo_update_properties(origin=origin, rot_x=Lx, rot_y=Ly, rot_z=Lz)



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
                        self.aa_main_chain.append(atm)
                    else:
                        self.aa_side_chain.append(atm)

            elif isinstance(frag, NucleicAcidResidue):
                for atm in frag.iter_atoms():
                    self.dna_main_chain.append(atm)

            elif frag.is_water():
                for atm in frag.iter_atoms(): 
                    self.water.append(atm)

            else:
                for atm in frag.iter_atoms():
                    self.hetatm.append(atm)

        if len(self.aa_main_chain)>0:
            self.aa_main_chain.glo_set_properties_id("aa_main_chain")
            self.glo_add_child(self.aa_main_chain)
            self.glo_link_child_property(
                "aa_main_chain_visible", "aa_main_chain", "visible")
        else:
            self.aa_main_chain = None

        if len(self.aa_side_chain)>0:
            self.aa_side_chain.glo_set_properties_id("aa_side_chain")
            self.glo_add_child(self.aa_side_chain)
            self.glo_link_child_property(
                "aa_side_chain_visible", "aa_side_chain", "visible")
        else:
            self.aa_side_chain = None

        if len(self.dna_main_chain)>0:
            self.dna_main_chain.glo_set_properties_id("dna_main_chain")
            self.glo_add_child(self.dna_main_chain)
            self.glo_link_child_property(
                "dna_main_chain_visible", "dna_main_chain", "visible")
        else:
            self.dna_main_chain = None

        if len(self.dna_side_chain)>0:
            self.dna_side_chain.glo_set_properties_id("dna_side_chain")
            self.glo_add_child(self.dna_side_chain)
            self.glo_link_child_property(
                "dna_side_chain_visible", "dna_side_chain", "visible")
        else:
            self.dna_side_chain = None

        if len(self.hetatm)>0:
            self.hetatm.glo_set_properties_id("hetatm")
            self.glo_add_child(self.hetatm)
            self.glo_link_child_property(
                "hetatm_visible", "hetatm", "visible")
        else:
            self.hetatm = None

        if len(self.water)>0:
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

    def iter_draw_lists(self):
        """Iterate over all GL Lists.
        """
        if self.aa_main_chain!=None:  yield self.aa_main_chain
        if self.aa_side_chain!=None:  yield self.aa_side_chain
        if self.dna_main_chain!=None: yield self.dna_main_chain
        if self.dna_side_chain!=None: yield self.dna_side_chain
        if self.hetatm!=None:         yield self.hetatm
        if self.water!=None:          yield self.water
        

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
                "U", str(chain), "U")

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
              "cascade":           True,
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
        self.glo_remove_child(draw_list)
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

	glLight(GL_LIGHT0, GL_AMBIENT,  [1.0, 1.0, 1.0, 1.0])
	glLight(GL_LIGHT0, GL_DIFFUSE,  [1.0, 1.0, 1.0, 1.0])
	glLight(GL_LIGHT0, GL_SPECULAR, [0.2, 0.2, 0.2, 1.0])
	glLight(GL_LIGHT0, GL_POSITION, [1.0, 1.0, 1.0, 0.0])
	
	glLightModel(GL_LIGHT_MODEL_AMBIENT, [0.2, 0.2, 0.2, 1.0])

 	glEnable(GL_LIGHTING)
 	glEnable(GL_LIGHT0)
	
	glDepthFunc(GL_LESS)
	glEnable(GL_DEPTH_TEST)
		
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

	glTranslatef(self.xpos, self.ypos, self.zpos)

	glRotatef(self.rotx, 1.0, 0.0, 0.0)
	glRotatef(self.roty, 0.0, 1.0, 0.0)
        glRotatef(self.rotz, 0.0, 0.0, 1.0)

        for draw_list in self.glo_iter_children():
            draw_list.gl_render()
            
	if gl_drawable.is_double_buffered():
            gl_drawable.swap_buffers()
	else:
            glFlush()

        gl_drawable.gl_end()
