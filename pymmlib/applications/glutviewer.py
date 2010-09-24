#!/usr/bin/env python
## Copyright 2002-2010 by PyMMLib Development Group (see AUTHORS file)
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


from mmLib import FileIO, Structure, Viewer, OpenGLDriver


## tweakers (tweakable constants)
##
BG_COLOR         = "0.3,0.3,0.3"
GL_DOUBLE_BUFFER = True


##
## constants
##
WELCOME = """\
WARNING!!! This application is incomplete!

===============================================================================
                        mmLIB TLSViewer: GLUT Version
-------------------------------------------------------------------------------
Terminal Commands:
    [ESC]                press escape to show/hide terminal
    load <path>          load a PDB or mmCIF structure file
    bg_color [color]     show/set the background color
    exit/quit            exits the viewer

mmLib Virtual Filesystem:
    ls                   list contents of current virtual file system
                         directory
    cd [directory]       change to the given virtual file system directory
    
Mouse Navigation:
    Right Mouse Button   "trackball" style rotation of structure
    Middle Mouse Button  x,y translation of structure
    Left Mouse Button    zoom in/zoom out
===============================================================================

"""

CHAR_ASCENT = 119.05
CHAR_DECENT = 33.3
CHAR_HEIGHT = CHAR_ASCENT + CHAR_DECENT


##
## console output
##
def error(text):
    sys.stderr.write("[VIEWER:ERROR] %s\n" % (text))
    sys.stderr.flush()
    
def info(text):
    sys.stderr.write("[VIEWER:INFO]  %s\n" % (text))
    sys.stderr.flush()

def draw_panel(z, x1, y1, x2, y2, r, g, b, a):
    """Draw transparent glass panel with outline.
    """
    ## transparent smoky glass
    glEnable(GL_LIGHTING)
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE)

    glMaterialfv(GL_FRONT, GL_AMBIENT,  (0.0, 0.0, 0.0, a))
    glMaterialfv(GL_FRONT, GL_DIFFUSE,  (0.1, 0.1, 0.1, a))
    glMaterialfv(GL_FRONT, GL_SPECULAR, (0.1*r, 0.1*g, 0.1*b, a))
    glMaterialfv(GL_FRONT, GL_EMISSION, (0.05*r, 0.05*g, 0.05*b, a))
    glMaterialfv(GL_FRONT, GL_SHININESS, 128.0)

    glNormal3f(0.0, 0.0, 1.0)

    glBegin(GL_QUADS)
    glVertex3f(x1, y1, z - 0.0001)
    glVertex3f(x2, y1, z - 0.0001)
    glVertex3f(x2, y2, z - 0.0001)
    glVertex3f(x1, y2, z - 0.0001)
    glEnd()

    ## outline the screen
    glDisable(GL_LIGHTING)
    glColor3f(r, g, b)
    glLineWidth(1.0)

    glBegin(GL_LINE_LOOP)
    glVertex3f(x1, y1, z)
    glVertex3f(x2, y1, z)
    glVertex3f(x2, y2, z)
    glVertex3f(x1, y2, z)
    glEnd()
    

##
## OpenGL Terminal Hack
##
class Terminal(object):
    """Terminal window for controlling mmLib.Viewer options.
    """
    def __init__(self, cobj):
        self.visible = True
        self.width   = 0
        self.height  = 0
        self.zplane  = 0.0

        self.term_alpha = 0.90

        self.char_width = 80
        
        self.wind_border = 5.0
        self.term_border = 1.0

        self.cobj = cobj
        self.prompt = ""
        self.construct_prompt()
        
        self.lines = []

    def construct_prompt(self):
        self.prompt = "%s# " % (self.cobj.pwd)

    def keypress(self, key):
        ascii = ord(key)
        #info("keypress: %d" % (ascii))

        ## enter
        if ascii==13:
            ln = self.lines[0][len(self.prompt):]
            output = self.cobj.command(ln)
            self.write(output)
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
        if text!=None and len(text)>0:
            for ln in text.split("\n"):
                self.lines.insert(0, ln)
        self.lines.insert(0, self.prompt)

    def opengl_render(self):
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
        glEnable(GL_CULL_FACE)

        ## alpha blending func
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

        ## light 0 disable
        glDisable(GL_LIGHT0)

        ## light 1
        glEnable(GL_LIGHT1)
        glLightfv(GL_LIGHT1, GL_AMBIENT, (0.0, 0.0, 0.0, 1.0))
        glLightfv(GL_LIGHT1, GL_DIFFUSE, (1.0, 1.0, 1.0, 1.0))
        glLightfv(GL_LIGHT1, GL_SPECULAR, (1.0, 1.0, 1.0, 1.0))

        glLightfv(
            GL_LIGHT1,
            GL_POSITION,
            (glwidth/2.0, glheight/2.0, zplane + 20.0, 0.0))

        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, (0.2, 0.2, 0.2, 1.0))

        ##
        ## draw background
        ##
        draw_panel(zplane,
                   self.wind_border, self.wind_border,
                   glwidth - self.wind_border, glheight - self.wind_border,
                   0.0, 1.0, 0.0, 0.8)

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

            ## draw a cursor on the first line
            if i==0:
                draw_panel(0.0,
                           10.0, -CHAR_DECENT, 70.0, CHAR_ASCENT,
                           1.0, 1.0, 1.0, 0.8)
                
                glColor3f(0.0, 1.0, 0.0)

            glPopMatrix()
            i += 1


class GLUT_Viewer(Viewer.GLViewer):
    """The main OpenGL Viewer using GLUT.
    """
    def __init__(self):
        self.glut_init_done = False
        self.struct_desc_list = []
        self.opengl_driver = OpenGLDriver.OpenGLDriver()

        ## pixel width and height of window
        self.width  = 640
        self.height = 480

        ## mouse navigation state
        self.in_drag         = False
        self.navigation_mode = None
        self.beginx          = 0
        self.beginy          = 0

        ## current virutal directory
        self.pwd = "/"

        ## terminal
        self.term = Terminal(self)
        self.term.write(WELCOME)

        Viewer.GLViewer.__init__(self)
        self.properties.update(bg_color=BG_COLOR)

    def command(self, cmd):
        cmd = cmd.strip()
        
        if cmd.startswith("load"):
            return self.command_load_struct(cmd)

        elif cmd.startswith("ls"):
            return self.command_ls(cmd)

        elif cmd.startswith("cd"):
            return self.command_cd(cmd)

        elif cmd.startswith("bg_color"):
            return self.command_bg_color(cmd)

        elif cmd.startswith("exit") or cmd.startswith("quit"):
            return self.command_exit(cmd)

        if len(cmd.strip())==0:
            return ""
        
        return "unknown command: %s" % (cmd)

    def command_bg_color(self, cmd):
        cmd = cmd[8:].strip()

        if len(cmd)==0:
            return self.properties["bg_color"]

        self.properties.update(bg_color=cmd)
        return "background color: %s" % (self.properties["bg_color"])

    def command_load_struct(self, cmd):
        path = cmd[5:].strip()
        struct = self.load_struct(path)
        if struct==None:
            return "file not found: %s" % (path)

        return "loaded %s:%s" % (path, struct.structure_id)

    def command_ls(self, cmd):
        if self.pwd=="/":
            x = ""
            for glstruct in self.glo_iter_children():
                x += "%s\n" % (glstruct.struct.structure_id)
            return x

    def command_exit(self, cmd):
        sys.exit(0)

    def command_cd(self, cmd):
        path = cmd[3:].strip()
        for dir in path.split("/"):
            self.command_cd_one(dir)        
        self.term.construct_prompt()

    def command_cd_one(self, dir):
        dir = dir.strip()

        if len(dir)==0:
            return (True, "")

        if self.pwd=="/":
            if dir=="..":
                return (False, "invalid path")

            for glstruct in self.glo_iter_children():
                if glstruct.struct.structure_id==dir:
                    return (True, dir)

    def load_struct(self, path):
        """Loads the requested structure.
        """
        info("loading: %s" % (path))
        
        try:
            struct = FileIO.LoadStructure(
                fil = path,
                distance_bonds = True)
        except IOError:
            error("file not found: %s" % (path))
            return None

        struct_desc = {}
        struct_desc["struct"] = struct
        struct_desc["path"] = path

        glstruct = self.glv_add_struct(struct)
        return struct

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
            glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
        else:
            glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH)

        glutInitWindowSize(self.width, self.height); 
        glutInitWindowPosition(100, 100)
        glutCreateWindow("GLUT TLSViewer")

        glutDisplayFunc(self.glut_display)
        glutReshapeFunc(self.glut_reshape)
        glutMouseFunc(self.glut_mouse)
        glutMotionFunc(self.glut_motion)
        glutKeyboardFunc(self.glut_keyboard)

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
        ## GLUT initalization isn't really done until the
        ## window has been shaped
        self.glut_init_done = True
        
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
