## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""OpenGL rendering classes.
"""
from __future__ import generators
from OpenGL.GL   import *
from OpenGL.GLU  import *
from OpenGL.GLUT import *
from mmTypes import *
from Structure import *


class GLProperties(dict):
    """A very specialized dictionary object for handling the many selectable
    properties of the drawing objects.
    """
    def __init__(self):
        dict.__init__(self)
        self.prop_list = []

    def add(self, prop):
        """Add a property.
        """
        self.prop_list.append(prop)
        self[prop["name"]] =  prop["default"]

    def init(self, **args):
        """This is a special form of update which propagates all linked
        values, not just the changed ones.
        """
        for prop in self.prop_list:
            name       = prop["name"]
            self[name] = args.get(name, prop["default"])

            try:
                linked_props = prop["linked properties"]
            except KeyError:
                pass
            else:
                for linked_prop in linked_props:
                    object = linked_prop["object"]
                    args   =  {linked_prop["name"]: self[name]}
                    object.properties.init(**args)

    def update(self, **args):
        """Update the properties, and any linked values
        """
        action_list = []

        for prop in self.prop_list:
            name      = prop["name"]
            old_value = self[name]

            try:
                self[name] = args[name]
            except KeyError:
                pass
            else:
                if self[name]!=old_value:
                    if prop["action"] not in action_list:
                        action_list.append(prop["action"])

                try:
                    linked_props = prop["linked properties"]
                except KeyError:
                    pass
                else:
                    for linked_prop in linked_props:
                        object = linked_prop["object"]
                        args   =  {linked_prop["name"]: self[name]}
                        object.update(**args)

        return action_list


class GLDrawList(object):
    """Fundamental OpenGL rigid entity.
    """
    def __init__(self, **args):
        self.gl_name = None

        ## reflection properties of material
        self.ambient  = 1.0
        self.diffuse  = 1.0
        self.specular = 0.2
        self.material = (1.0, 1.0, 1.0, 1.0)

        self.properties = GLProperties()
        self.properties.add(
            { "name" :      "visible",
              "default":    True,
              "action":     "redraw" })
        self.properties.add(
            { "name" :      "origin",
              "default":    zeros(3, Float),
              "action":     "redraw" })
        self.properties.add(
            { "name" :      "axes",
              "default":    identity(3),
              "action":     "redraw" })
        self.properties.add(
            { "name" :      "rot_x",
              "default":    0.0,
              "action":     "redraw" })
        self.properties.add(
            { "name" :      "rot_y",
              "default":    0.0,
              "action":     "redraw" })
        self.properties.add(
            { "name" :      "rot_z",
              "default":    0.0,
              "action":     "redraw" })

        self.properties.init(**args)
        
    def update(self, **args):
        self.properties.update(**args)

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

    def set_material(self, r=1.0, g=1.0, b=1.0, br=1.0, alpha=1.0):
        """Utility function used in subclasses to set material properties
        before drawing.
        """
        amb_r = self.ambient * r * br
        amb_g = self.ambient * g * br
        amb_b = self.ambient * b * br

        dif_r = self.diffuse * br
        spe_r = self.specular * br

        shine = 50.0 * br
        
	glMaterial(GL_FRONT, GL_AMBIENT, [amb_r, amb_g, amb_b, alpha])
	glMaterial(GL_FRONT, GL_DIFFUSE, [dif_r, dif_r, dif_r, alpha])
	glMaterial(GL_FRONT, GL_SPECULAR,[spe_r, spe_r, spe_r, alpha])
	glMaterial(GL_FRONT, GL_SHININESS, shine)


class GLAxes(GLDrawList):
    """Draw orthogonal axes in red = x, green = y, blue = z.
    """
    def __init__(self, **args):
        GLDrawList.__init__(self)

        self.properties.add(
            { "name":       "line_length",
              "default":    200.0,
              "action":     "recompile" })
        self.properties.add(
            { "name":       "line_width",
              "default":    10.0,
              "action":     "recompile" })
        self.properties.add(
            { "name":       "color_x",
              "default":    (1.0, 0.0, 0.0),
              "action":     "recompile" })
        self.properties.add(
            { "name":       "color_y",
              "default":    (0.0, 1.0, 0.0),
              "action":     "recompile" })
        self.properties.add(
            { "name":       "color_z",
              "default":    (0.0, 0.0, 1.0),
              "action":     "recompile" })

        self.properties.update(**args)
        
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
    def __init__(self, unit_cell, **args):
        GLDrawList.__init__(self)
        self.unit_cell = unit_cell

        self.properties.add(
            { "name":       "line_width",
              "default":    2.0,
              "action":     "recompile" })
        self.properties.add(
            { "name":       "color",
              "default":    (1.0, 1.0, 1.0),
              "action":     "recompile" })

        self.properties.init(**args)

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
        GLDrawList.__init__(self)
        AtomList.__init__(self)

        self.properties.add(
            { "name":       "sphere_quality",
              "default":    12,
              "action":     "recompile" })
        self.properties.add(
            { "name":       "line_width",
              "default":    3.0,
              "action":     "recompile" })
        self.properties.add(
            { "name":       "atom_origin",
              "default":    None,
              "action":     "recompile" })
        self.properties.add(
            { "name":       "lines",
              "default":    True,
              "action":     "recompile" })
        self.properties.add(
            { "name":       "cpk",
              "default":    False,
              "action":    "recompile" })
        self.properties.add(
            { "name":      "U",
              "default":   False,
              "action":    "recompile" })
        self.properties.add(
            { "name":      "atm_U_attr",
              "default":   "U",
              "action":    "recompile" })
        self.properties.add(
            { "name":      "U_color",
              "default":   (1.0, 1.0, 1.0),
              "action":    "recompile" })
        self.properties.add(
            { "name":      "color",
              "default":   None,
              "action":    "recompile" })
        self.properties.add(
            { "name":      "color_func",
              "default":   None,
              "action":    "recompile" })
        
        self.properties.init(**args)
        self.el_color_cache = {}

    def update(self, **args):
        #print "## GLAtomList.update(%s)" % (str(args))
        actions = self.properties.update(**args)
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
        

class GLTLSGroup(GLDrawList):
    """Draws TLS group
    """
    def __init__(self, **args):
        GLDrawList.__init__(self)

        self.properties.add(
            { "name":       "time",
              "default":    0.0,
              "action":     "redraw" })
        self.properties.add(
            { "name":       "line_width",
              "default":    2.0,
              "action":     "recompile" })
        self.properties.add(
            { "name":       "TLS_visible",
              "default":    True,
              "action":     "recompile" })
        self.properties.add(
            { "name":       "T_color",
              "default":    (0.0, 1.0, 0.0),
              "action":     "recompile" })
        self.properties.add(
            { "name":       "L_color",
              "default":    (1.0, 0.0, 0.0),
              "action":     "recompile" })
        self.properties.add(
            { "name":       "S_color",
              "default":    (0.0, 0.0, 1.0),
              "action":     "recompile" })
        self.properties.add(
            { "name":       "CA_line_visible",
              "default":    False,
              "action":     "recompile" })

        ## TLS calculations 
        self.tls_group = args["tls_group"] 
        self.calcs     = self.tls_group.calc_COR()        
        (eigen_values, eigen_vectors) = eigenvectors(self.tls_group.L)
        self.evalL   = eval

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
            origin      = self.calcs["COR"], 
            atom_origin = self.calcs["COR"])

        for atm, Ucalc in self.tls_group.iter_atm_Ucalc():
            atm.Ucalc = Ucalc
            self.gl_atom_list.append(atm)

    def update(self, **args):
        self.gl_atom_list.update(**args)
        actions = self.properties.update(**args)
        if "recompile" in actions:
            self.gl_delete_list()
        if "update_time" in actions:
            self.update_time()

    def gl_render(self):
        GLDrawList.gl_render(self)
        self.gl_atom_list.gl_render()

    def gl_delete_list(self):
        GLDrawList.gl_delete_list(self)
        self.gl_atom_list.gl_delete_list()
        
    def gl_draw(self):
        ## draw atoms
        self.gl_atom_list.gl_draw()
        
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
        glLineWidth(self.properties["line_width"])
        
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

        self.update(origin=origin, rot_x=Lx, rot_y=Ly, rot_z=Lz)


class GLChain(GLDrawList):
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
                { "name":              "aa_main_chain_visible",
                  "default":           True,
                  "action":            "redraw",
                  "linked properties": [{"object": self.aa_main_chain,
                                         "name":   "visible" }] })
        else:
            self.aa_main_chain = None

        if len(self.aa_side_chain)>0:
            self.properties.add(
                { "name":              "aa_side_chain_visible",
                  "default":           True,
                  "action":            "redraw",
                  "linked properties": [{"object": self.aa_side_chain,
                                        "name":   "visible" }] })
        else:
            self.aa_side_chain = None

        if len(self.dna_main_chain)>0:
            self.properties.add(
                { "name":              "dna_main_chain_visible",
                  "default":           True,
                  "action":            "redraw",
                  "linked properties": [{"object": self.dna_main_chain,
                                        "name":   "visible" }] })
        else:
            self.dna_main_chain = None

        if len(self.dna_side_chain)>0:
            self.properties.add(
                { "name":              "dna_side_chain_visible",
                  "default":           True,
                  "action":            "redraw",
                  "linked properties": [{"object": self.dna_side_chain,
                                        "name":   "visible" }] })
        else:
            self.dna_side_chain = None

        if len(self.hetatm)>0:
            self.properties.add(
                { "name":              "hetatm_visible",
                  "default":           True,
                  "action":            "redraw",
                  "linked properties": [{"object": self.hetatm,
                                        "name":   "visible" }] })
        else:
            self.hetatm = None

        if len(self.water)>0:
            self.properties.add(
                { "name":              "water_visible",
                  "default":           True,
                  "action":            "redraw",
                  "linked properties": [{"object": self.water,
                                        "name":   "visible" }] })
        else:
            self.water = None

        ## link GLChain's color property to the atom lists
        def linked_values(property):
            linked = []
            for gl_atom_list in self.iter_draw_lists():
                linked.append({"object": gl_atom_list, "name": property})
            return linked

        self.properties.add(
            { "name":              "color",
              "default":           None,
              "action":            "redraw",
              "linked properties": linked_values("color") })

        self.properties.add(
            { "name":              "U",
              "default":           False,
              "action":            "redraw",
              "linked properties": linked_values("U") })

        self.properties.init(**args)

    def update(self, **args):
        self.properties.update(**args)
        
    def iter_draw_lists(self):
        """Iterate over all GL Lists.
        """
        if self.aa_main_chain!=None:  yield self.aa_main_chain
        if self.aa_side_chain!=None:  yield self.aa_side_chain
        if self.dna_main_chain!=None: yield self.dna_main_chain
        if self.dna_side_chain!=None: yield self.dna_side_chain
        if self.hetatm!=None:         yield self.hetatm
        if self.water!=None:          yield self.water
        
    def gl_render(self):
        self.gl_push_matrix()
        for draw_list in self.iter_draw_lists():
            draw_list.gl_render()
        self.gl_pop_matrix()

    def gl_draw(self):
        for draw_list in self.iter_draw_lists():
            draw_list.gl_draw()

    def gl_delete_list(self):
        """Delete all OpenGL draw lists.
        """
        for draw_list in self.iter_draw_lists():
            draw_list.gl_delete_list()


class GLStructure(GLDrawList):
    def __init__(self, **args):
        GLDrawList.__init__(self)
        self.struct         = args["struct"]
        self.gl_axes        = GLAxes()
        self.gl_unit_cell   = GLUnitCell(self.struct.unit_cell)
        self.gl_chain_dict  = {}

        for chain in self.struct.iter_chains():
            self.gl_chain_dict[chain.chain_id] = GLChain(chain=chain)

        self.properties.add(
            { "name":     "axes_visible",
              "default":  True,
              "action":  "redraw",
              "linked properties": [{ "object": self.gl_axes,
                                      "name":   "visible" }] })
        self.properties.add(
            { "name":       "unit_cell_visible",
              "default":    True,
              "action":     "redraw",
              "linked properties": [{ "object": self.gl_unit_cell,
                                      "name":   "visible" }] })

        def gl_chain_linked(property):
            linked = []
            for gl_chain in self.gl_chain_dict.values():
                linked.append({ "object": gl_chain,
                                "name":   property })
            return linked

        self.properties.add(
            { "name":              "aa_main_chain_visible",
              "default":           True,
              "action":            "redraw",
              "linked properties": gl_chain_linked("aa_main_chain_visible") })
        self.properties.add(
            { "name":              "aa_side_chain_visible",
              "default":           True,
              "action":            "redraw" ,
              "linked properties": gl_chain_linked("aa_side_chain_visible") })
        self.properties.add(
            { "name":              "dna_main_chain_visible",
              "default":           True,
              "action":            "redraw" ,
              "linked properties": gl_chain_linked("dna_main_chain_visible") })
        self.properties.add(
            { "name":              "dna_side_chain_visible",
              "default":           True,
              "action":            "redraw",
              "linked properties": gl_chain_linked("dna_side_chain_visible") })
        self.properties.add(
            { "name":              "hetatm_visible",
              "default":           True,
              "action":            "redraw",
              "linked properties": gl_chain_linked("hetatm_visible") })
        self.properties.add(
            { "name":              "water_visible",
              "default":           True,
              "action":            "redraw",
              "linked properties": gl_chain_linked("water_visible") })
        self.properties.add(
            { "name":              "color",
              "default":           None,
              "action":            "redraw",
              "linked properties": gl_chain_linked("color") })
        self.properties.add(
            { "name":              "U",
              "default":           False,
              "action":            "redraw",
              "linked properties": gl_chain_linked("U") })

        self.properties.init(**args)

    def update(self, **args):
        #print "## GLStructure.update(%s)" % (str(args))
        self.properties.update(**args)
        
    def iter_draw_lists(self):
        """Iterate over all GL Lists.
        """
        yield self.gl_axes
        yield self.gl_unit_cell
        for gl_chain in self.gl_chain_dict.values():
            yield gl_chain

    def gl_render(self):
        self.gl_push_matrix()
        for draw_list in self.iter_draw_lists():
            draw_list.gl_render()
        self.gl_pop_matrix()

    def gl_draw(self):
        for draw_list in self.iter_draw_lists():
            draw_list.gl_draw()

    def gl_delete_list(self):
        """Delete all OpenGL draw lists.
        """
        for draw_list in self.iter_draw_lists():
            draw_list.gl_delete_list()

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


class GLViewer(list):
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
    def __init__(self, glcontext, gldrawable):
        list.__init__(self)
        
        ## gldrawable and glcontext to draw on
        self.glcontext  = glcontext
        self.gldrawable = gldrawable

        ## position and rotation of viewer window
        self.xpos = 0.0
        self.ypos = 0.0
        self.zpos = -50.0
        self.rotx = 0.0
        self.roty = 0.0
        self.rotz = 0.0

    def append(self, draw_list):
        """Append a GLDrawList.
        """
        assert isinstance(draw_list, GLDrawList)
        list.append(self, draw_list)

    def remove(self, draw_list):
        """Remove a GLDrawList.
        """
        assert isinstance(draw_list, GLDrawList)
        list.remove(self, draw_list)

    def gl_init(self):
        """Called once to initalize the GL scene before drawing.
        """
        if not self.gldrawable.gl_begin(self.glcontext):
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
		
        self.gldrawable.gl_end()
    
    def gl_resize(self, width, height):
        """Called to set the size of the OpenGL window this class is
        drawing on.
        """
	if not self.gldrawable.gl_begin(self.glcontext):
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
	self.gldrawable.gl_end()

    def gl_draw_lists(self):
        """Draw all GLDrawList objects onto the given glcontext/gldrawable.
        If the GLDrawList objects are not yet compiled into OpenGL draw
        lists, they will be compiled while they are drawn, since this is
        a useful optimization.
        """
	if not self.gldrawable.gl_begin(self.glcontext):
            return

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
	glLoadIdentity()

	glTranslatef(self.xpos, self.ypos, self.zpos)

	glRotatef(self.rotx, 1.0, 0.0, 0.0)
	glRotatef(self.roty, 0.0, 1.0, 0.0)
        glRotatef(self.rotz, 0.0, 0.0, 1.0)

        for draw_list in self:
            draw_list.gl_render()
            
	if self.gldrawable.is_double_buffered():
            self.gldrawable.swap_buffers()
	else:
            glFlush()

        self.gldrawable.gl_end()
