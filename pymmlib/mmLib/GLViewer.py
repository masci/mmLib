## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""OpenGL rendering classes.
"""
from OpenGL.GL   import *
from OpenGL.GLU  import *
from OpenGL.GLUT import *
from mmTypes import *
from Structure import *


class GLDrawList:
    def __init__(self):
        self.name   = None
        self.origin = Vector(0.0, 0.0, 0.0)
        self.axes   = identity(3)
        self.rotx   = 0.0
        self.roty   = 0.0
        self.rotz   = 0.0

    def gl_compile_list(self, execute = False):
        if self.name != None:
            self.gl_delete_list()
        self.name = glGenLists(1)
        debug("gl_compile_list name="+str(self.name))
        if execute:
            glNewList(self.name, GL_COMPILE_AND_EXECUTE)
        else:
            glNewList(self.name, GL_COMPILE)
        self.gl_draw()
        glEndList()

    def gl_delete_list(self):
        if self.name != None:
            glDeleteLists(self.name, 1)
            self.name = None

    def gl_call_list(self):
        glPushMatrix()
        glTranslatef(*self.origin)
	glRotatef(self.rotx, self.axes[0,0], self.axes[0,1], self.axes[0,2])
	glRotatef(self.roty, self.axes[1,0], self.axes[1,1], self.axes[1,2])
        glRotatef(self.rotz, self.axes[2,0], self.axes[2,1], self.axes[2,2])
        glCallList(self.name)
        glPopMatrix()
        
    def gl_draw(self):
        pass


class GLUnitCellDrawList(GLDrawList):
    def __init__(self, unit_cell):
        GLDrawList.__init__(self)
        self.unit_cell = unit_cell

        ## draw constents
        self.cell_line_width = 2.0
        self.axis_line_width = 10.0

        ## reflection properties of material
        self.ambient  = 0.4
        self.diffuse  = 0.8
        self.specular = 1.0

    def set_material(self, r, g, b, br):
        alpha = 1.0

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

    def draw_cell(self, x1, y1, z1, x2, y2, z2):
        """Draw the unit cell lines in a rectangle starting at fractional
        integer coordinates x1, y1, z1, ending at x2, y2, z2.  The first set of
        coordinates must be <= the second set.
        """
        assert x1 <= x2 and y1 <= y2 and z1 <= z2

        a = self.unit_cell.calc_frac_to_orth(Vector(1.0, 0.0, 0.0))
        b = self.unit_cell.calc_frac_to_orth(Vector(0.0, 1.0, 0.0))
        c = self.unit_cell.calc_frac_to_orth(Vector(0.0, 0.0, 1.0))

        glLineWidth(self.cell_line_width)

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

    def draw_axis(self):
        def axis_line(v1, v2):
            glLineWidth(self.axis_line_width)
            glBegin(GL_LINES)
            glVertex3f(*v1)
            glVertex3f(*v2)
            glEnd()
                    
        self.set_material(1.0, 0.0, 0.0, 1.0)
        axis_line(Vector(0.0, 0.0, 0.0), Vector(200.0, 0.0, 0.0))

        self.set_material(0.0, 1.0, 0.0, 1.0)
        axis_line(Vector(0.0, 0.0, 0.0), Vector(1.0, 200.0, 0.0))

        self.set_material(0.0, 0.0, 1.0, 1.0)
        axis_line(Vector(0.0, 0.0, 0.0), Vector(0.0, 0.0, 200.0))

    def gl_draw(self):
        self.set_material(1.0, 1.0, 1.0, 1.0)
        self.draw_cell(-1, -1, -1, 0, 0, 0)
        self.draw_axis()


class GLAtomDrawList(GLDrawList, AtomList):
    """OpenGL renderer for a list of atoms.
    """
    def __init__(self):
        GLDrawList.__init__(self)
        AtomList.__init__(self)

        ## draw settings
        self.settings = {}

        ## reflection properties of material
        self.ambient  = 0.4
        self.diffuse  = 0.8
        self.specular = 1.0

        ## defaults
        self.sphere_quality = 12
        self.line_width     = 3.0
        self.atom_draw_func = self.draw_bond
        self.atom_origin    = None

    def set_material(self, r, g, b, br):
        alpha = 1.0

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

    def draw_cpk(self, atm, symop = None):
        """Draw a atom as a CPK sphere.
        """
        el = atm.get_structure().library.get_element(atm.element)
        if el:
            (r, g, b) = el.color
            self.set_material(r, g, b, 1.0)
            radius = el.van_der_waals_radius
        else:
            self.set_material(1.0, 1.0, 1.0, 1.0)
            radius = 2.0

        glPushMatrix()
        if self.atom_origin:
            glTranslatef(*atm.position - self.atom_origin)
        else:
            glTranslatef(*atm.position)
        glutSolidSphere(radius, self.sphere_quality, self.sphere_quality)
        glPopMatrix()

    def draw_bond(self, atm, symop = None):
        """Draw a atom using bond lines only.
        """
        el = atm.get_structure().library.get_element(atm.element)
        if el:
            (r, g, b) = el.color
            self.set_material(r, g, b, 1.0)
        else:
            self.set_material(1.0, 1.0, 1.0, 1.0)
        
        glLineWidth(self.line_width)
        if atm.bond_list:
            for bond in atm.iter_bonds():
                atm2 = bond.get_partner(atm)

                if bond.bond_type != None:
                    glLineWidth(10.0)
                else:
                    glLineWidth(self.line_width)

                if self.atom_origin:
                    start = atm.position - self.atom_origin
                else:
                    start = atm.position

                end   = start + ((atm2.position - atm.position) / 2)

                glBegin(GL_LINES)
                glVertex3f(*start)
                glVertex3f(*end)
                glEnd()

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
                    v = Vector(x+tx, y+ty, z+tz)
                    atm.position = uc.calc_frac_to_orth(v)
                    self.draw_atom(atm, False)

            atm.position = orig_pos
            return
        
    def gl_draw(self):
        """Perform the OpenGL drawing operations to render this atom list
        with the current settings.
        """
        if self.settings.get("main_chain"):
            for atm in self:
                frag = atm.get_fragment()
                if isinstance(frag, AminoAcidResidue):
                    if atm.name in ["C", "CA", "O", "N"]:
                        self.atom_draw_func(atm)
                else:
                    self.atom_draw_func(atm)
        else:
            for atm in self:            
                self.atom_draw_func(atm)

    def set_setting(self, key, setting):
        """Set settings for the drawing engine.
        """
        success = True

        ## set the origin of the atoms in the atom list
        if key == "atom_origin":
            assert len(setting) == 3
            self.atom_origin = setting

        ## set the draw style of the atom list
        elif key == "draw_style":
            assert setting == "cpk" or setting == "bond"
            if setting == "cpk":
                self.atom_draw_func = self.draw_cpk
            elif setting == "bond":
                self.atom_draw_func = self.draw_bond

        ## draw main chain atoms only
        elif key == "main_chain":
            self.settings["main_chain"] = setting

        ## draw symmetry equivelant atoms into unit cells listed as
        ## 3-tuples of neighboring cell, the origin unit cell is (0,0,0)
        elif key == "draw_symmetry_equivalents":
            self.settings["draw_symmetry_equivalents"] = setting

        else:
            success = False

        ## if a setting was reset, then delete the current draw list
        ## so it will be regenerated
        if success:
            self.gl_delete_list()

        return success

    def unset_setting(self, key):
        """Unset drawing engine settings.
        """
        success = True

        if key == "atom_origin":
            self.atom_origin = None

        elif key == "draw_symmetry_equivalents":
            del self.setting["draw_symmetry_equivalents"]

        elif key == "main_chain":
            try: del self.settings["main_chain"]
            except KeyError: pass
            
        if success:
            self.gl_delete_list()

        return success


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
	glLight(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
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
            if draw_list.name == None:
                draw_list.gl_compile_list(execute = 1)
            else:
                draw_list.gl_call_list()
            
	if self.gldrawable.is_double_buffered():
            self.gldrawable.swap_buffers()
	else:
            glFlush()

        self.gldrawable.gl_end()
