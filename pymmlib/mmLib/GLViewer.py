## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

from OpenGL.GL   import *
from OpenGL.GLU  import *
from OpenGL.GLUT import *

from mmTypes import *


class GLDrawList:
    def __init__(self):
        self.name = None
        self.origin = Vector(0.0, 0.0, 0.0)
        self.axes = identity(3)
        self.rotx = 0.0
        self.roty = 0.0
        self.rotz = 0.0

    def gl_compile_list(self, execute = 0):
        if self.name != None:
            gl_delete_lists(self.name, 1)
        
        self.name = glGenLists(1)

        print "gl_compile_list name = ", self.name

        if execute: glNewList(self.name, GL_COMPILE_AND_EXECUTE)
        else:       glNewList(self.name, GL_COMPILE)
        
        self.gl_draw()
        glEndList()

    def gl_delete_list(self):
        if self.name != None:
            gl_delete_lists(self.name, 1)
            self.name = None

    def gl_call_list(self):
        #print "gl_call_list = ",self.name

        glPushMatrix()

        apply(glTranslatef, self.origin)

	glRotatef(self.rotx, self.axes[0,0], self.axes[0,1], self.axes[0,2])
	glRotatef(self.roty, self.axes[1,0], self.axes[1,1], self.axes[1,2])
        glRotatef(self.rotz, self.axes[2,0], self.axes[2,1], self.axes[2,2])
        
        gl_call_list(self.name)

        glPopMatrix()
        
    def gl_draw(self):
        pass

    def set_origin(self, origin):
        """Reset the origin of the draw list.
        """
        self.origin = origin




class GLUnitCellDrawList(GLDrawList):
    def __init__(self, unit_cell):
        GLDrawList.__init__(self)
        
        self.unit_cell = unit_cell

    def gl_draw(self):
        self.set_material(1.0, 1.0, 1.0, 1.0)

        def cell_line(v1, v2):
            glLineWidth(2.0)
            glBegin(GL_LINES)
            apply(glVertex3f, v1)
            apply(glVertex3f, v2)
            glEnd()

        m = self.unit_cell.calc_cartisian_unit_cell_axes()

        def draw_cell(i, j, k):
            o = i*m[0] + j*m[1] + k*m[2]
            cell_line(o, o+m[0])
            cell_line(o, o+m[1])
            cell_line(o, o+m[2])

        rng = range(-1, 2)

        for i in rng:
            for j in rng:
                for k in rng:
                    draw_cell(i, j, k)

        ## draw cartesian axes
        def axis_line(v1, v2):
            glLineWidth(5.0)
            glBegin(GL_LINES)
            apply(glVertex3f, v1)
            apply(glVertex3f, v2)
            glEnd()
                    
        self.set_material(1.0, 0.0, 0.0, 1.0)
        axis_line(Vector(0.0, 0.0, 0.0), Vector(100.0, 0.0, 0.0))

        self.set_material(0.0, 1.0, 0.0, 1.0)
        axis_line(Vector(0.0, 0.0, 0.0), Vector(1.0, 100.0, 0.0))

        self.set_material(0.0, 0.0, 1.0, 1.0)
        axis_line(Vector(0.0, 0.0, 0.0), Vector(0.0, 0.0, 100.0))

        
    def set_material(self, r, g, b, brightness):
	glMaterial(GL_FRONT, GL_AMBIENT,
                   [0.2 * r * brightness,
                    0.2 * g * brightness,
                    0.2 * b * brightness,
                    1.0])
        
	glMaterial(GL_FRONT, GL_DIFFUSE,
                   [0.8 * brightness,
                    0.8 * brightness,
                    0.8 * brightness,
                    1.0])
        
	glMaterial(GL_FRONT, GL_SPECULAR,
                   [1.0 * brightness,
                    1.0 * brightness,
                    1.0 * brightness,
                    1.0])
        
	glMaterial(GL_FRONT, GL_SHININESS, 50.0 * brightness)


class GLAtomDrawList(GLDrawList, list):
    def __init__(self):
        GLDrawList.__init__(self)
        list.__init__(self)
        self.atom_origin = Vector(0.0, 0.0, 0.0)
        self.draw_u = 0

    def set_material(self, r, g, b, brightness):
	glMaterial(GL_FRONT, GL_AMBIENT,
                   [0.2 * r * brightness,
                    0.2 * g * brightness,
                    0.2 * b * brightness,
                    1.0])
        
	glMaterial(GL_FRONT, GL_DIFFUSE,
                   [0.8 * brightness,
                    0.8 * brightness,
                    0.8 * brightness,
                    1.0])
        
	glMaterial(GL_FRONT, GL_SPECULAR,
                   [1.0 * brightness,
                    1.0 * brightness,
                    1.0 * brightness,
                    1.0])
        
	glMaterial(GL_FRONT, GL_SHININESS, 50.0 * brightness)

    def select_atom_material(self, atm):
        r = 1.0
        g = 1.0
        b = 1.0
        brightness = 1.0
        
        if   atm.element == "C":
            b = 1.0
            g = 1.0
            b = 0.0

        elif atm.element == "N":
            r = 0.0
            g = 0.0
            b = 1.0

        elif atm.element == "O":
            r = 1.0
            g = 0.0
            b = 0.0

        elif atm.element == "S":
            r = 0.0
            g = 1.0
            b = 0.0

        return (r, g, b, brightness)

    def draw_atom(self, atm):
        glPushMatrix()

        apply(glTranslatef, atm.position - self.atom_origin)

        (r, g, b, brightness) = self.select_atom_material(atm)
        self.set_material(r, g, b, brightness)
        glutSolidSphere(0.1, 12, 12)

        ## U axes
        m = array([[ atm.U[0], atm.U[3], atm.U[4] ],
                   [ atm.U[3], atm.U[1], atm.U[5] ],
                   [ atm.U[4], atm.U[5], atm.U[2] ]])
        
        (eval, evec) = eigenvectors(m)

        if self.draw_u:
            self.set_material(1.0, 1.0, 1.0, 1.0)
            for i in range(3):
                glLineWidth(1.0)
                glBegin(GL_LINES)
                apply(glVertex3f,  evec[i])
                apply(glVertex3f, -evec[i])
                glEnd()

        glPopMatrix()
        
    def draw_bond(self, atm1, atm2):
        glLineWidth(5.0)
        glBegin(GL_LINES)
        apply(glVertex3f, atm1.position - self.atom_origin)
        apply(glVertex3f, atm2.position - self.atom_origin)
        glEnd()

    def gl_draw(self):
        ## draw origin
        self.set_material(1.0, 1.0, 1.0, 1.0)
        glutSolidSphere(1.0, 32, 32)

        ## draw atoms and bonds
        visited_bonds = []

        for atm in self:
            ## draw atom
            self.draw_atom(atm)

            ## draw bonds
            for bond in atm.iter_bonds():
                if bond in visited_bonds: continue
                else:                     visited_bonds.insert(0, bond)

                atm2 = bond.get_partner(atm)
                self.draw_bond(atm, atm2)

    def set_atom_origin(self, atom_origin):
        self.atom_origin = atom_origin
        self.gl_delete_list()

        

class GLViewer(list):
    """Inherits from Python's list object.  To draw """
    
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
        assert isinstance(draw_list, GLDrawList)
        list.append(self, draw_list)

    def remove(self, draw_list):
        assert isinstance(draw_list, GLDrawList)
        list.remove(self, draw_list)

    def gl_init(self):
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
	if not self.gldrawable.gl_begin(self.glcontext):
            return
	
	glViewport(0, 0, width, height)
	glMatrixMode(GL_PROJECTION)
	glLoadIdentity()

	if width > height:
            w = float(width) / float(height)
            glFrustum(-w, w, -1.0, 1.0, 3.0, 300.0)
	else:
            h = float(height) / float(width)
            glFrustum(-1.0, 1.0, -h, h, 3.0, 300.0)
	
	glMatrixMode(GL_MODELVIEW)
	self.gldrawable.gl_end()

    def gl_drawLists(self):
        #print "begin"
        
	if not self.gldrawable.gl_begin(self.glcontext):
            return

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
	glLoadIdentity()

	glTranslatef(self.xpos, self.ypos, self.zpos)

	glRotatef(self.rotx, 1.0, 0.0, 0.0)
	glRotatef(self.roty, 0.0, 1.0, 0.0)
        glRotatef(self.rotz, 0.0, 0.0, 1.0)

        #glutSolidSphere(3.0, 12, 12)

        for draw_list in self:
            if draw_list.name == None: draw_list.gl_compile_list(execute = 1)
            else:                      draw_list.gl_call_list()
            
	if self.gldrawable.is_double_buffered(): self.gldrawable.swap_buffers()
	else:                                    glFlush()

        self.gldrawable.gl_end()

        #print "end"
