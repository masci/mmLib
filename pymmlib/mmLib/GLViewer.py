## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

from OpenGL.GL   import *
from OpenGL.GLU  import *
from OpenGL.GLUT import *


class GLDrawList:
    def __init__(self):
        self.name = None
        self.xpos = 0.0
        self.ypos = 0.0
        self.zpos = 0.0
        self.rotx = 0.0
        self.roty = 0.0
        self.rotz = 0.0

    def glCompileList(self, execute = 0):
        if self.name != None:
            glDeleteLists(self.name, 1)
        
        self.name = glGenLists(1)

        print "glCompileList name = ", self.name

        if execute: glNewList(self.name, GL_COMPILE_AND_EXECUTE)
        else:       glNewList(self.name, GL_COMPILE)

        glPushMatrix()
	glTranslatef(self.xpos, self.ypos, self.zpos)
	glRotatef(self.rotx, 1.0, 0.0, 0.0)
	glRotatef(self.roty, 0.0, 1.0, 0.0)
        glRotatef(self.rotz, 0.0, 0.0, 1.0)
        
        self.glDraw()

        glPopMatrix()
        glEndList()

    def glDeleteList(self):
        glDeleteLists(self.name, 1)
        self.name = None

    def glDraw(self):
        pass


class GLUnitCellDrawList(GLDrawList):
    def __init__(self, unit_cell):
        GLDrawList.__init__(self)
        
        self.unit_cell = unit_cell

    def glDraw(self):
        lw = 2.0
        self.set_material(1.0, 1.0, 1.0, 1.0)

        def cell_line(v1, v2):
            glLineWidth(lw)
            glBegin(GL_LINES)
            glVertex3f(v1[0], v1[1], v1[2])
            glVertex3f(v2[0], v2[1], v2[2])
            glEnd()

        m = self.unit_cell.calcCartisianUnitCellAxes()

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

        glTranslatef(atm.position[0], atm.position[1], atm.position[2])
        (r, g, b, brightness) = self.select_atom_material(atm)
        self.set_material(r, g, b, brightness)
        glutSolidSphere(0.5, 12, 12)

        glPopMatrix()
        
    def draw_bond(self, atm1, atm2):
        glLineWidth(5.0)
        glBegin(GL_LINES)
        glVertex3f(atm1.position[0], atm1.position[1], atm1.position[2])
        glVertex3f(atm2.position[0], atm2.position[1], atm2.position[2])
        glEnd()

    def glDraw(self):
        visited_bonds = []
        
        for atm in self:
            ## draw atom
            self.draw_atom(atm)

            ## draw bonds
            for bond in atm.iterBonds():
                if bond in visited_bonds: continue
                else:                     visited_bonds.insert(0, bond)

                atm2 = bond.getPartner(atm)
                self.draw_bond(atm, atm2)


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

    def glInit(self):
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
    
    def glResize(self, width, height):
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

    def glDrawLists(self):
        print "begin"
        
	if not self.gldrawable.gl_begin(self.glcontext):
            return

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
	glLoadIdentity()

	glTranslatef(self.xpos, self.ypos, self.zpos)

	glRotatef(self.rotx, 1.0, 0.0, 0.0)
	glRotatef(self.roty, 0.0, 1.0, 0.0)
        glRotatef(self.rotz, 0.0, 0.0, 1.0)

        glutSolidSphere(3.0, 12, 12)

        for draw_list in self:
            print "glDrawLists name = ", draw_list.name
            
            if draw_list.name == None: draw_list.glCompileList(execute = 1)
            else:                      glCallList(draw_list.name)
            
	if self.gldrawable.is_double_buffered(): self.gldrawable.swap_buffers()
	else:                                    glFlush()

        self.gldrawable.gl_end()

        print "end"
