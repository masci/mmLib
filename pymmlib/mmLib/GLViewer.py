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


class GLDrawList(object):
    """Fundamental OpenGL rigid entity.
    """
    def __init__(self):
        self.gl_name = None
        self.show    = True
        self.origin  = Vector(0.0, 0.0, 0.0)
        self.axes    = identity(3)
        self.rotx    = 0.0
        self.roty    = 0.0
        self.rotz    = 0.0


        ## reflection properties of material
        self.ambient  = 1.0
        self.diffuse  = 1.0
        self.specular = 0.2
        self.material = (1.0, 1.0, 1.0, 1.0)

    def gl_push_matrix(self):
        """Rotate and translate to the correct position for drawing.
        """
        glPushMatrix()
        glTranslatef(*self.origin)
        glRotatef(
            self.rotx, self.axes[0,0], self.axes[0,1], self.axes[0,2])
        glRotatef(
            self.roty, self.axes[1,0], self.axes[1,1], self.axes[1,2])
        glRotatef(
            self.rotz, self.axes[2,0], self.axes[2,1], self.axes[2,2])

    def gl_pop_matrix(self):
        """Pop the roatated/translated position.
        """
        glPopMatrix()

    def gl_call_list_render(self, render = True):
        """Compile or force a recompile of this object's gl_draw list, and
        render the scene.  Rendering the scene can be bypassed if
        this method is called with render = False.
        """
        if self.show == False:
            return
        
        if self.gl_name != None:
            self.gl_push_matrix()
            glCallList(self.gl_name)
            self.gl_pop_matrix()
            return

        self.gl_name = glGenLists(1)
        debug("gl_compile_list name = " + str(self.gl_name))

        if render:
            self.gl_push_matrix()
            glNewList(self.gl_name, GL_COMPILE_AND_EXECUTE)
            self.gl_draw()
            glEndList()
            self.gl_pop_matrix()
        else:
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
    def __init__(self):
        GLDrawList.__init__(self)
        self.axis_line_width = 10.0

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
        self.set_material(*self.material)
        self.draw_axis()


class GLUnitCell(GLDrawList):
    """Draw unit cell.
    """
    def __init__(self, unit_cell):
        GLDrawList.__init__(self)
        self.unit_cell = unit_cell

        ## draw constents
        self.cell_line_width = 2.0

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

    def gl_draw(self):
        self.set_material(*self.material)
        self.draw_cell(-1, -1, -1, 0, 0, 0)


class GLTLSGroup(GLDrawList):
    """Draws TLS group
    """
    def __init__(self, **args):
        GLDrawList.__init__(self)

        self.tls_group = args["tls_group"] 

    def gl_draw(self):
        ## calcuate center of reaction and related tensors

        print self.tls_group.name

        try:
            calcs = self.tls_group.calc_COR()
        except:
            return

        cor = Vector(calcs["COR"][0], calcs["COR"][1], calcs["COR"][2])
        
##         ## L axes
##         (evalL, evecL) = eigenvectors(self.tls_group.L)
##         self.set_material(0.0, 0.5, 0.2, 1.0)
##         for v in evecL:
##             evec = Vector(v[0], v[1], v[2])
##             print evalL
            
##             glBegin(GL_LINES)
##             glVertex3f(* cor - evec)
##             glVertex3f(* cor + evec)
##             glEnd()

        self.set_material(0.0, 0.5, 0.2, 1.0)
        glPushMatrix()
        glTranslatef(*cor)
        glutSolidSphere(trace(self.tls_group.L)*rad2deg2/24.0, 16, 16)
        glPopMatrix()

        ## draw a line from the center of reaction to each atom
        self.set_material(0.5, 0.5, 0.5, 1.0)
        for atm in self.tls_group:
            if atm.name in ["CA"]:
                glBegin(GL_LINES)
                glVertex3f(*cor)
                glVertex3f(*atm.position)
                glEnd()


class GLAtomList(GLDrawList, AtomList):
    """OpenGL renderer for a list of atoms.
    """
    def __init__(self, **args):
        GLDrawList.__init__(self)
        AtomList.__init__(self)

        ## draw settings
        if args.has_key("atom_style"):
            if args["atom_style"] == "bonds":
                self.atom_draw_func = self.draw_bond
            elif args["atom_style"] == "cpk":
                self.atom_draw_func = self.draw_cpk
            else:
                self.atom_draw_func = self.draw_bond
        else:
            self.atom_draw_func = self.draw_bond

        ## reflection properties of material
        self.ambient  = 0.4
        self.diffuse  = 0.8
        self.specular = 1.0
        self.material = (1.0, 1.0, 1.0, 1.0)

        ## defaults
        self.sphere_quality = 12
        self.line_width     = 3.0
        self.atom_draw_func = self.draw_bond
        self.atom_origin    = None

    def gl_draw(self):
        """Perform the OpenGL drawing operations to render this atom list
        with the current settings.
        """
        for atm in self:            
            self.atom_draw_func(atm)

    def draw_cpk(self, atm, symop = None):
        """Draw a atom as a CPK sphere.
        """
        el = atm.get_structure().library.get_element(atm.element)
        if el:
            (r, g, b) = el.color
            self.set_material(r, g, b, 1.0)
            radius = el.van_der_waals_radius
        else:
            self.set_material(*self.material)
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

        ## if there are bonds, then draw the lines 1/2 way to the
        ## bonded atoms
        if len(atm.bond_list) > 0:
            for bond in atm.iter_bonds():
                atm2 = bond.get_partner(atm)

                if self.atom_origin:
                    start = atm.position - self.atom_origin
                else:
                    start = atm.position

                end   = start + ((atm2.position - atm.position) / 2)

                glBegin(GL_LINES)
                glVertex3f(*start)
                glVertex3f(*end)
                glEnd()

        ## if there are no bonds, draw a small cross-point 
        else:
            start = atm.position - Vector(0.25, 0.0, 0.0)
            end   = atm.position + Vector(0.25, 0.0, 0.0)
            glBegin(GL_LINES)
            glVertex3f(*start)
            glVertex3f(*end)
            glEnd()
            
            start = atm.position - Vector(0.0, 0.25, 0.0)
            end   = atm.position + Vector(0.0, 0.25, 0.0)
            glBegin(GL_LINES)
            glVertex3f(*start)
            glVertex3f(*end)
            glEnd()

            start = atm.position - Vector(0.0, 0.0, 0.25)
            end   = atm.position + Vector(0.0, 0.0, 0.25)
            glBegin(GL_LINES)
            glVertex3f(*start)
            glVertex3f(*end)
            glEnd()

    def draw_U_axes(self, atm, symop = None):
        """Draw principal thermal axes for atom.
        """
        if atm.U == None:
            return

        evec, eval = eigenvectors(atm.U)

        glLineWidth(self.line_width)
        
        glBegin(GL_LINES)
        glVertex3f(*atm.position - evec[0])
        glVertex3f(*atm.position + evec[0])
        glEnd()

        glBegin(GL_LINES)
        glVertex3f(*atm.position - evec[1])
        glVertex3f(*atm.position + evec[1])
        glEnd()

        glBegin(GL_LINES)
        glVertex3f(*atm.position - evec[2])
        glVertex3f(*atm.position + evec[2])
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
        

class GLStructure(GLDrawList):
    def __init__(self, struct):
        GLDrawList.__init__(self)

        self.struct         = struct
        self.gl_axes        = None
        self.gl_unit_cell   = None
        self.draw_lists     = []
        self.atom_list_dict = {}

    def gl_call_list_render(self, render=True):
        self.gl_push_matrix()
        for draw_list in self.draw_lists:
            draw_list.gl_call_list_render(render=render)
        self.gl_pop_matrix()

    def gl_draw(self):
        for draw_list in self.draw_lists:
            draw_list.gl_draw()

    def gl_delete_list(self):
        """Delete all OpenGL draw lists.
        """
        for draw_list in self.draw_lists:
            draw_list.gl_delete_list()

        self.draw_lists   = []
        self.gl_axes      = None
        self.gl_unit_cell = None

    def show_axes(self, show):
        """True/False to show the coordinate axes.
        """
        if show == True and self.gl_axes == None:
            self.gl_axes = GLAxes()
            self.draw_lists.append(self.gl_axes)
            
        elif show == False and self.gl_axes != None:
            self.draw_lists.remove(self.gl_axes)
            self.gl_axes.gl_delete_list()
            self.gl_axes = None

    def show_unit_cell(self, show):
        """True/False to show the unit cell.
        """
        if show == True and self.gl_unit_cell == None:
            self.gl_unit_cell = GLUnitCell(self.struct.unit_cell)
            self.draw_lists.append(self.gl_unit_cell)
            
        elif show == False and self.gl_unit_cell != None:
            self.draw_lists.remove(self.gl_unit_cell)
            self.gl_unit_cell.gl_delete_list()
            self.gl_unit_cell = None

    def get_atom_list(self, list_id):
        """Get or create a new GLAtomList and return it.
        """
        try:
            return self.atom_list_dict[list_id]
        except KeyError:
            pass
        
        gl_atom_list = GLAtomList()
        self.atom_list_dict[list_id] = gl_atom_list
        self.draw_lists.append(gl_atom_list)

        return gl_atom_list

    



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
            draw_list.gl_call_list_render(render=True)
            
	if self.gldrawable.is_double_buffered():
            self.gldrawable.swap_buffers()
	else:
            glFlush()

        self.gldrawable.gl_end()
