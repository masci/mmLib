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


class GLProperties(dict):
    """A very specialized dictionary object for handling the many selectable
    properties of the drawing objects.
    """
    def __init__(self):
        dict.__init__(self)
        self.prop_list     = []
        self.callback_list = []

    def get_property(self, name):
        """Returns a property dictionary by property name.
        """
        for prop in self.prop_list:
            if prop["name"]==name:
                return prop
        return None

    def add(self, prop):
        """Add a property dictionary.  Do not allow properties with the same
        name.
        """
        assert self.get_property(prop["name"])==None
        
        self.prop_list.append(prop)
        self[prop["name"]] =  prop["default"]

    def add_update_callback(self, func):
        """Adds a function which is called whenever property values change.
        The function is called with two arguments: a updates dictionary
        containing all updated properties and the values they were changed
        to, and a actions list which contains a unique list of action
        key words formed from the action fields of the updated properties.
        """
        self.callback_list.append(func)

    def remove_update_callback(self, func):
        """Removes the update callback.
        """
        self.callback_list.remove(func)

    def init(self, **args):
        """This is a special form of update which propagates all linked
        values, not just the changed ones.
        """
        for prop in self.prop_list:
            name       = prop["name"]
            self[name] = args.get(name, prop["default"])

            try:
                linked_props = prop["link"]
            except KeyError:
                pass
            else:
                for linked_prop in linked_props:
                    properties = linked_prop["gl_properties"]
                    args = {linked_prop["name"]: self[name]}
                    properties.init(**args)

    def update(self, **args):
        """Update the properties, and any linked values
        """
        print "## update(%s)" % (str(args))
        
        updates = {}
        actions = []

        for prop in self.prop_list:
            name      = prop["name"]
            old_value = self[name]

            try:
                self[name] = args[name]
            except KeyError:
                pass
            else:
                if self[name]!=old_value or prop.has_key("link"):
                    updates[name] = self[name]

                    if prop["action"] not in actions:
                        actions.append(prop["action"])

                try:
                    linked_props = prop["link"]
                except KeyError:
                    pass
                else:
                    for linked_prop in linked_props:
                        properties = linked_prop["gl_properties"]
                        args = {linked_prop["name"]: self[name]}
                        properties.update(**args)

        for func in self.callback_list:
            func(updates, actions)



class GLDrawListProperties(GLProperties):
    def __init__(self):
        GLProperties.__init__(self)
        
        self.add(
            { "name" :      "visible",
              "desc":       "Visible",
              "type":       "boolean",
              "default":    True,
              "action":     "redraw" })
        self.add(
            { "name" :      "origin",
              "desc":       "Origin",
              "type":       "array(3)",
              "hidden":     True,
              "default":    zeros(3, Float),
              "action":     "redraw" })
        self.add(
            { "name" :      "axes",
              "desc":       "Rotation Axes",
              "type":       "array(3,3)",
              "hidden":     True,
              "default":    identity(3),
              "action":     "redraw" })
        self.add(
            { "name" :      "rot_x",
              "desc":       "Degrees Rotation About X Axis",
              "type":       "float",
              "hidden":     True,
              "default":    0.0,
              "action":     "redraw" })
        self.add(
            { "name" :      "rot_y",
              "desc":       "Degrees Rotation About Y Axis",
              "type":       "float",
              "hidden":     True,
              "default":    0.0,
              "action":     "redraw" })
        self.add(
            { "name" :      "rot_z",
              "desc":       "Degrees Rotation About Z Axis",
              "type":       "float",
              "hidden":     True,
              "default":    0.0,
              "action":     "redraw" })
        

class GLDrawList(object):
    """Fundamental OpenGL rigid entity.
    """
    def __init__(self, **args):
        self.gl_name = None

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
    def iter_draw_lists(self):
        """Iterates all draw lists.
        """
        pass
    
    def gl_render(self):
        if self.properties["visible"]==False:
            return
        
        self.gl_push_matrix()

        for draw_list in self.iter_draw_lists():
            draw_list.gl_render()

        if self.gl_name==None:
            self.gl_compile_list()

        glCallList(self.gl_name)

        self.gl_pop_matrix()

    def gl_delete_list(self):
        for draw_list in self.iter_draw_lists():
            draw_list.gl_delete_list()
        GLDrawList.gl_delete_list(self)


class GLAxesProperties(GLDrawListProperties):
    def __init__(self):
        GLDrawListProperties.__init__(self)

        self.add(
            { "name":       "line_length",
              "type":       "float",
              "default":    200.0,
              "action":     "recompile" })
        self.add(
            { "name":       "line_width",
              "type":       "float",
              "default":    10.0,
              "action":     "recompile" })
        self.add(
            { "name":       "color_x",
              "type":       "color",
              "default":    (1.0, 0.0, 0.0),
              "action":     "recompile" })
        self.add(
            { "name":       "color_y",
              "type":       "color",
              "default":    (0.0, 1.0, 0.0),
              "action":     "recompile" })
        self.add(
            { "name":       "color_z",
              "type":       "color",
              "default":    (0.0, 0.0, 1.0),
              "action":     "recompile" })


class GLAxes(GLDrawList):
    """Draw orthogonal axes in red = x, green = y, blue = z.
    """
    def __init__(self, **args):
        GLDrawList.__init__(self, **args)

        self.properties = args.get("properties", GLAxesProperties())
        self.properties.add_update_callback(self.update_cb)
        self.properties.init(**args)

    def update_cb(self, updates, actions):
        if "recompile" in actions:
            self.gl_delete_list()
        
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



class GLUnitCellProperties(GLDrawListProperties):
    def __init__(self):
        GLDrawListProperties.__init__(self)
        
        self.add(
            { "name":       "line_width",
              "type":       "float",
              "default":    2.0,
              "action":     "recompile" })
        self.add(
            { "name":       "color",
              "type":       "color",
              "default":    (1.0, 1.0, 1.0),
              "action":     "recompile" })


class GLUnitCell(GLDrawList):
    """Draw unit cell.
    """
    def __init__(self, **args):
        GLDrawList.__init__(self, **args)

        self.properties = args.get("properties", GLUnitCellProperties())
        self.properties.add_update_callback(self.update_cb)
        self.properties.init(**args)

        self.unit_cell = args["unit_cell"]

    def update_cb(self, updates, actions):
        if "recompile" in actions:
            self.gl_delete_list()

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


class GLAtomListProperties(GLDrawListProperties):
    def __init__(self):
        GLDrawListProperties.__init__(self)
        
        self.add(
            { "name":       "sphere_quality",
              "desc":       "CPK Sphere Quality",
              "type":       "integer",
              "default":    12,
              "action":     "recompile" })
        self.add(
            { "name":       "line_width",
              "desc":       "Atom Line Drawing Width",
              "type":       "float",
              "default":    3.0,
              "action":     "recompile" })
        self.add(
            { "name":       "atom_origin",
              "type":       "array(3)",
              "hidden":     True,
              "default":    None,
              "action":     "recompile" })
        self.add(
            { "name":       "lines",
              "desc":       "Draw Atom Bond Lines",
              "type":       "boolean",
              "default":    True,
              "action":     "recompile" })
        self.add(
            { "name":       "cpk",
              "desc":       "Draw CPK Spheres",
              "type":       "boolean",
              "default":    False,
              "action":    "recompile" })
        self.add(
            { "name":      "U",
              "desc":      "Draw ADP Axes",
              "type":      "boolean",
              "default":   False,
              "action":    "recompile" })
        self.add(
            { "name":      "atm_U_attr",
              "type":      "string",
              "hidden":    True,
              "default":   "U",
              "action":    "recompile" })
        self.add(
            { "name":      "U_color",
              "desc":      "Thermal Axes Color",
              "type":      "color",
              "default":   (1.0, 1.0, 1.0),
              "action":    "recompile" })
        self.add(
            { "name":      "color",
              "desc":      "Solid Color",
              "type":      "color",
              "default":   None,
              "action":    "recompile" })
        self.add(
            { "name":      "color_func",
              "type":      "function",
              "hidden":    True,
              "default":   None,
              "action":    "recompile" })
        

class GLAtomList(GLDrawList, AtomList):
    """OpenGL renderer for a list of atoms.  Optional arguments iare:
    color, U, U_color.
    """
    def __init__(self, **args):
        GLDrawList.__init__(self, **args)
        AtomList.__init__(self)

        self.properties = args.get("properties", GLAtomListProperties())
        self.properties.add_update_callback(self.update_cb)
        self.properties.init(**args)

        self.el_color_cache = {}

    def update_cb(self, updates, actions):
        print "## GLAtomList.update(%s)" % (str(updates))
        if "recompile" in actions:
            self.gl_delete_list()

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
        

class GLTLSGroupProperties(GLDrawListProperties):
    def __init__(self):
        GLDrawListProperties.__init__(self)
        
        self.add(
            { "name":        "time",
              "desc":        "Simulation Time", 
              "type":        "float",
              "default":     0.0,
              "action":      "redraw" })
        self.add(
            { "name":       "tensor_line_width",
              "desc":       "Line Width of Tensors",
              "type":       "float",
              "default":    2.0,
              "action":     "recompile" })
        self.add(
            { "name":        "TLS_visible",
              "desc":        "Show TLS Tensors",
              "type":        "boolean",
              "default":     True,
              "action":      "recompile" })
        self.add(
            { "name":        "T_color",
              "desc":        "T Tensor Color",
              "type":        "color",
              "default":     (0.0, 1.0, 0.0),
              "action":      "recompile" })
        self.add(
            { "name":       "L_color",
              "desc":       "L Tensor Color",
              "type":       "color",
              "default":    (1.0, 0.0, 0.0),
              "action":     "recompile" })
        self.add(
            { "name":       "S_color",
              "desc":       "S Tensor Color",
              "type":       "color",
              "default":    (0.0, 0.0, 1.0),
              "action":     "recompile" })
        self.add(
            { "name":        "CA_line_visible",
              "desc":        "Show Lines to C-Alpha Atoms",
              "type":        "boolean",
              "default":     False,
              "action":      "recompile" })


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

        ## properties 
        self.properties = args.get("properties", GLTLSGroupProperties())
        self.properties.add_update_callback(self.update_cb)
        
        self.properties.init(
            origin = self.calcs["COR"],
            axes   = eigen_vectors,
            **args)

        ## set up a GLAtomList to display the group
        self.gl_atom_list = GLAtomList(
            U           = True,
            U_atm_arg   = "Ucalc",
            U_color     = (0.,1.,0.),
            color       = (1.,1.,1.),
            line_width  = 10.0,
            atom_origin = self.calcs["COR"])

        for atm, Ucalc in self.tls_group.iter_atm_Ucalc():
            atm.Ucalc = Ucalc
            self.gl_atom_list.append(atm)

    def update_cb(self, updates, actions):
        if "recompile" in actions:
            self.gl_delete_list()
        if "time" in updates:
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

        sin_tm = math.sin(self.properties["time"])

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

        self.properties.update(origin=origin, rot_x=Lx, rot_y=Ly, rot_z=Lz)



class GLChain(GLDrawListContainer):
    def __init__(self, **args):
        GLDrawListContainer.__init__(self, **args)

        self.properties = args.get("properties", GLDrawListProperties())
        self.properties.add_update_callback(self.update_cb)

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
            self.properties.add(
                { "name":     "aa_main_chain_visible",
                  "type":     "boolean",
                  "default":  True,
                  "action":   "redraw",
                  "link":     [{"gl_properties": self.aa_main_chain.properties,
                                "name":          "visible" }] })
        else:
            self.aa_main_chain = None

        if len(self.aa_side_chain)>0:
            self.properties.add(
                { "name":     "aa_side_chain_visible",
                  "type":     "boolean",
                  "default":  True,
                  "action":   "redraw",
                  "link":     [{"gl_properties": self.aa_side_chain.properties,
                                "name":          "visible" }] })
        else:
            self.aa_side_chain = None

        if len(self.dna_main_chain)>0:
            self.properties.add(
                { "name":      "dna_main_chain_visible",
                  "type":      "boolean",
                  "default":   True,
                  "action":    "redraw",
                  "link":      [{"gl_properties": self.dna_main_chain.properties,
                                 "name":          "visible" }] })
        else:
            self.dna_main_chain = None

        if len(self.dna_side_chain)>0:
            self.properties.add(
                { "name":      "dna_side_chain_visible",
                  "type":      "boolean",
                  "default":   True,
                  "action":    "redraw",
                  "link":      [{"gl_properties": self.dna_side_chain.properties,
                                 "name":          "visible" }] })
        else:
            self.dna_side_chain = None

        if len(self.hetatm)>0:
            self.properties.add(
                { "name":      "hetatm_visible",
                  "type":      "boolean",
                  "default":   True,
                  "action":    "redraw",
                  "link":      [{"gl_properties": self.hetatm.properties,
                                 "name":          "visible" }] })
        else:
            self.hetatm = None

        if len(self.water)>0:
            self.properties.add(
                { "name":      "water_visible",
                  "type":      "boolean",
                  "default":   True,
                  "action":    "redraw",
                  "link":      [{"gl_properties": self.water.properties,
                                 "name":          "visible" }] })
        else:
            self.water = None

        ## link GLChain's color property to the atom lists
        def linked_values(property):
            linked = []
            for gl_atom_list in self.iter_draw_lists():
                linked.append({"gl_properties": gl_atom_list.properties,
                               "name": property})
            return linked

        self.properties.add(
            { "name":      "color",
              "type":      "color",
              "default":   None,
              "action":    "redraw",
              "link":      linked_values("color") })

        self.properties.add(
            { "name":      "U",
              "type":      "matrix(3,3)",
              "default":   False,
              "action":    "redraw",
              "link":      linked_values("U") })

        self.properties.init(**args)

    def update_cb(self, updates, actions):
        print "GLChain.update(%s)" % (str(updates))
        if "recompile" in actions:
            self.gl_delete_list()
        
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

        self.properties = args.get("properties", GLDrawListProperties())
        self.properties.add_update_callback(self.update_cb)
        
        self.struct         = args["struct"]
        self.gl_axes        = GLAxes()
        self.gl_unit_cell   = GLUnitCell(unit_cell=self.struct.unit_cell)
        self.gl_chain_dict  = {}

        for chain in self.struct.iter_chains():
            self.gl_chain_dict[chain.chain_id] = GLChain(chain=chain)

        self.properties.add(
            { "name":     "axes_visible",
              "type":     "boolean",
              "default":  True,
              "action":  "redraw",
              "link":    [{ "gl_properties": self.gl_axes.properties,
                            "name":          "visible" }] })

        self.properties.add(
            { "name":     "unit_cell_visible",
              "type":     "boolean",
              "default":  True,
              "action":   "redraw",
              "link":     [{ "gl_properties": self.gl_unit_cell.properties,
                             "name":          "visible" }] })

        def gl_chain_linked(property):
            linked = []
            for gl_chain in self.gl_chain_dict.values():
                linked.append({ "gl_properties": gl_chain.properties,
                                "name":          property })
            return linked

        self.properties.add(
            { "name":              "aa_main_chain_visible",
              "default":           True,
              "action":            "redraw",
              "link": gl_chain_linked("aa_main_chain_visible") })
        self.properties.add(
            { "name":              "aa_side_chain_visible",
              "default":           True,
              "action":            "redraw" ,
              "link": gl_chain_linked("aa_side_chain_visible") })
        self.properties.add(
            { "name":              "dna_main_chain_visible",
              "default":           True,
              "action":            "redraw" ,
              "link": gl_chain_linked("dna_main_chain_visible") })
        self.properties.add(
            { "name":              "dna_side_chain_visible",
              "default":           True,
              "action":            "redraw",
              "link": gl_chain_linked("dna_side_chain_visible") })
        self.properties.add(
            { "name":              "hetatm_visible",
              "default":           True,
              "action":            "redraw",
              "link": gl_chain_linked("hetatm_visible") })
        self.properties.add(
            { "name":              "water_visible",
              "default":           True,
              "action":            "redraw",
              "link": gl_chain_linked("water_visible") })
        self.properties.add(
            { "name":              "color",
              "default":           None,
              "action":            "redraw",
              "link": gl_chain_linked("color") })
        self.properties.add(
            { "name":              "U",
              "default":           False,
              "action":            "redraw",
              "link": gl_chain_linked("U") })

        self.properties.init(**args)

    def update_cb(self, updates, actions):
        print "## GLStructure.update(%s)" % (str(updates))
        
    def iter_draw_lists(self):
        """Iterate over all GL Lists.
        """
        yield self.gl_axes
        yield self.gl_unit_cell
        for gl_chain in self.gl_chain_dict.values():
            yield gl_chain

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


class GLViewer(object):
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
        self.gl_draw_lists = []

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
        self.gl_draw_lists.append(draw_list)

    def remove_draw_list(self, draw_list):
        """Remove a GLDrawList.
        """
        assert isinstance(draw_list, GLDrawList)
        self.gl_draw_lists.remove(draw_list)

    def get_gl_context(self):
        """Implement in subclass to return the gl_context.
        """
        pass

    def get_gl_drawable(self):
        """Implement in subclass to reutrn the gl_drawable.
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

        for draw_list in self.gl_draw_lists:
            draw_list.gl_render()
            
	if gl_drawable.is_double_buffered():
            gl_drawable.swap_buffers()
	else:
            glFlush()

        gl_drawable.gl_end()
