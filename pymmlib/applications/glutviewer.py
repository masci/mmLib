#!/usr/bin/env python
## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys
import math
import random

from OpenGL.GL            import *
from OpenGL.GLU           import *
from OpenGL.GLUT          import *

from mmLib.PDB            import PDBFile
from mmLib.Structure      import *
from mmLib.FileLoader     import *
from mmLib.Viewer         import *
from mmLib.R3DDriver      import Raster3DDriver
from mmLib.OpenGLDriver   import OpenGLDriver
from mmLib.Extensions.TLS import *

##
## constants
##
WELCOME = """\
mmLIB TLSViewer: GLUT Version

Structure Navigation:

Terminal Commands:
    [ESC]         press escape to show/hide terminal
    load <path>   load a PDB or mmCIF structure file

"""


GL_DOUBLE_BUFFER = True

CHAR_ASCENT = 119.05
CHAR_DECENT = 33.3
CHAR_HEIGHT = CHAR_ASCENT + CHAR_DECENT


##
## console output
##
def error(text):
    sys.stderr.write("[GV ERROR] %s\n" % (text))
    sys.stderr.flush()
    
def info(text):
    sys.stderr.write("[GV INFO]  %s\n" % (text))
    sys.stderr.flush()


##
## OpenGL Terminal Hack
##
class Terminal(object):
    """Terminal window for controlling mmLib.Viewer options.
    """
    def __init__(self):
        self.visible = True
        self.width   = 0
        self.height  = 0
        self.zplane  = 5000.0

        self.term_alpha = 0.5

        self.char_width = 80
        
        self.wind_border = 5.0
        self.term_border = 1.0

        self.prompt = "> "
        self.lines = []

    def keypress(self, key):
        ascii = ord(key)
        #info("keypress: %d" % (ascii))

        ## enter
        if ascii==13:
            self.lines.insert(0, self.prompt)
        ## backspace
        elif ascii==8 or ascii==127:
            ln = self.lines[0]
            if len(ln)>len(self.prompt):
                self.lines[0] = ln[:-1]
        ## anything else
        else:
            self.lines[0] += key

    def write(self, text):
        """Writes text to the terminal.
        """
        for ln in text.split("\n"):
            self.lines.insert(0, ln)
        self.lines.insert(0, self.prompt)

    def opengl_render(self):
        ## setup viewport
        glViewport(0, 0, self.width, self.height)

        ## setup perspective matrix
	glMatrixMode(GL_PROJECTION)
 	glLoadIdentity()

        zplane = self.zplane
        near   = self.zplane + 1.0
        far    = self.zplane - 1.0

        ratio = float(self.height) / float(self.width)
        glwidth  = 80.0
        glheight = ratio * glwidth
        glOrtho(0.0, glwidth, 0.0, glheight, -near, -far)

        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()

        glClear(GL_DEPTH_BUFFER_BIT)

        ## OpenGL Features
        glEnable(GL_NORMALIZE)
        glShadeModel(GL_SMOOTH)
        glEnable(GL_DEPTH_TEST)
        glDepthFunc(GL_LESS)

        ## alpha blending func
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

        ## light 0
        glEnable(GL_LIGHT0)
        glLightfv(GL_LIGHT0, GL_AMBIENT, (0.0, 0.0, 0.0, 1.0))
        glLightfv(GL_LIGHT0, GL_DIFFUSE, (1.0, 1.0, 1.0, 1.0))
        glLightfv(GL_LIGHT0, GL_SPECULAR, (1.0, 1.0, 1.0, 1.0))
        glLightfv(GL_LIGHT0, GL_POSITION, (0.0, 0.0, zplane + 10.0, 0.0))

        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, (0.2, 0.2, 0.2, 1.0))

        ## light 1 disable
        glDisable(GL_LIGHT1)

        ##
        ## draw background
        ##
        glEnable(GL_LIGHTING)
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE)

        glMaterialfv(
            GL_FRONT,
            GL_AMBIENT,
            (0.0, 0.0, 0.0, self.term_alpha))
        
        glMaterialfv(
            GL_FRONT,
            GL_DIFFUSE,
            (0.1, 0.1, 0.1, self.term_alpha))
        
	glMaterialfv(
            GL_FRONT,
            GL_SPECULAR,
            (0.0, 0.1, 0.0, self.term_alpha))

        glMaterialfv(
            GL_FRONT,
            GL_EMISSION,
            (0.0, 0.1, 0.0, self.term_alpha))

        glMaterialfv(GL_FRONT, GL_SHININESS, 128.0)
        
        glBegin(GL_QUADS)

        glNormal3f(0.0, 0.0, 1.0)

        #glNormal3f(-0.5, -0.5, 1.0)
        glVertex3f(
            self.wind_border,
            self.wind_border,
            zplane)

        #glNormal3f(0.5, -0.5, 1.0)
        glVertex3f(
            glwidth - self.wind_border,
            self.wind_border,
            zplane)

        #glNormal3f(0.5, 0.5, 1.0)
        glVertex3f(
            glwidth - self.wind_border,
            glheight - self.wind_border,
            zplane)

        #glNormal3f(-0.5, 0.5, 1.0)
        glVertex3f(
            self.wind_border,
            glheight - self.wind_border,
            zplane)

        glEnd()

        ## outline the screen
        glDisable(GL_LIGHTING)
        glColor3f(0.0, 1.0, 0.0)
        glLineWidth(1.0)

        glBegin(GL_LINE_LOOP)

        glVertex3f(
            self.wind_border,
            self.wind_border,
            zplane + 0.1)

        glVertex3f(
            glwidth - self.wind_border,
            self.wind_border,
            zplane + 0.1)

        glVertex3f(
            glwidth - self.wind_border,
            glheight - self.wind_border,
            zplane + 0.1)

        glVertex3f(
            self.wind_border,
            glheight - self.wind_border,
            zplane + 0.1)

        glEnd()

        ## draw text lines
        glDisable(GL_LIGHTING)
        glColor3f(0.0, 1.0, 0.0)
        glLineWidth(1.0)

        ## compute a charactor scale which makes the height 1.0
        char1_scale = 1.0 / CHAR_HEIGHT

        i = 0
        for ln in self.lines:
            ypos = self.wind_border + self.term_border + 1.0*i
            if ypos>(glheight - self.wind_border - self.term_border):
                break

            glPushMatrix()
            glTranslatef(
                self.wind_border + self.term_border,
                self.wind_border + self.term_border + 1.0*i,
                zplane + 0.1)

            ## scale the chara
            glScalef(char1_scale, char1_scale, char1_scale)
            
            for c in ln:
                glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, ord(c))

            glPopMatrix()

            i += 1


class GLUT_Viewer(GLViewer):
    """The main OpenGL Viewer using GLUT.
    """
    def __init__(self):
        self.glut_init_done = False
        self.struct_desc_list = []
        self.opengl_driver = OpenGLDriver()

        ## pixel width and height of window
        self.width  = 640
        self.height = 480

        ## mouse navigation state
        self.in_drag         = False
        self.navigation_mode = None
        self.beginx          = 0
        self.beginy          = 0

        ## terminal
        self.term = Terminal()
        self.term.write(WELCOME)

        GLViewer.__init__(self)
        #self.properties.update(bg_color="White")
        
    def load_struct(self, path):
        """Loads the requested structure.
        """
        info("loading: %s" % (path))
        
        try:
            struct = LoadStructure(
                fil              = path,
                build_properties = ("library_bonds","distance_bonds"))
        except IOError:
            error("file not found: %s" % (path))
            return

        struct_desc = {}
        struct_desc["struct"] = struct
        struct_desc["path"] = path

        glstruct = self.glv_add_struct(struct)
        for glchain in glstruct.glo_iter_children():
            glchain.properties.update(
                ball_stick = True,
                ellipse    = True)

    def glv_render(self):
        """
        """
        self.glv_render_one(self.opengl_driver)

        if self.term.visible:
            self.term.opengl_render()
            
        if GL_DOUBLE_BUFFER:
            glutSwapBuffers()

        glFlush()

    def glv_redraw(self):
        """
        """
        if self.glut_init_done:
            glutPostRedisplay()

    def glut_init(self):
        """
        """
        glutInit(sys.argv)

        ## initalize OpenGL display mode
        if GL_DOUBLE_BUFFER:
            glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB)
        else:
            glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB)

        glutInitWindowSize(self.width, self.height); 
        glutInitWindowPosition(100, 100)
        glutCreateWindow("GLUT TLSViewer")

        glutDisplayFunc(self.glut_display)
        glutReshapeFunc(self.glut_reshape)
        glutMouseFunc(self.glut_mouse)
        glutMotionFunc(self.glut_motion)
        glutKeyboardFunc(self.glut_keyboard)

        self.glut_init_done = True

    def glut_main(self):
        """Run the GLUT main event loop.
        """
        glutMainLoop()

    def glut_display(self):
        """Render the scene.
        """
        self.glv_render()

    def glut_reshape(self, width, height):
        """Reshape the viewering window.
        """
        self.width  = width
        self.height = height

        self.term.width  = width
        self.term.height = height

        self.glv_resize(self.width, self.height)

    def glut_mouse(self, button, state, x, y):
        """Mouse button press callback.
        """
        if state==GLUT_UP:
            self.navigation_mode = None
            self.in_drag = False
            return

        self.in_drag = True
        
        if button==GLUT_LEFT_BUTTON:
            self.navigation_mode = "trackball"
        elif button==GLUT_MIDDLE_BUTTON:
            self.navigation_mode = "straif"
        elif button==GLUT_RIGHT_BUTTON:
            self.navigation_mode = "zoom"

        self.beginx = x
        self.beginy = y

    def glut_motion(self, x, y):
        """Mouse motion callback when one of the mouse buttons is pressed.
        """
        if not self.in_drag:
            return
        
        if self.navigation_mode=="trackball":
            self.glv_trackball(self.beginx, self.beginy, x, y)

        elif self.navigation_mode=="straif":
            dx = x - self.beginx
            dy = self.beginy - y
            self.glv_straif(dx, dy)

        elif self.navigation_mode=="zoom":
            if False:
                dx = event.x - self.beginx
                dy = self.beginy - event.y
                self.glv_clip(dy, dx)
            else:
                dy = y - self.beginy
                self.glv_zoom(dy)

        self.beginx = x
        self.beginy = y

    def glut_keyboard(self, key, x, y):
        """Keyboard press events.
        """
        ## toggle terminal/command mode
        if ord(key)==27:
            self.term.visible = not self.term.visible
            glutPostRedisplay()

        ## terminal mode -- route keystrokes to the terminal
        elif self.term.visible:
            self.term.keypress(key)
            glutPostRedisplay()

        ## command mode
        else:
            key = key.lower()
        
            ## quit
            if key=="q":
                sys.exit(0)
    

def main():
    gv = GLUT_Viewer()
    gv.glut_init()


    for path in sys.argv[1:]:
        gv.load_struct(path)
    
    gv.glut_main()


if __name__=="__main__":
    main()
