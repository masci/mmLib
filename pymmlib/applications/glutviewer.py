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
GL_DOUBLE_BUFFER = False


class Terminal(object):
    """Terminal window for controlling mmLib.Viewer options.
    """
    def __init__(self, glut_viewer):
        self.x      = 0
        self.y      = 0
        self.width  = 0
        self.height = 0


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

        GLViewer.__init__(self)
        
    def error(self, text):
        sys.stderr.write("[GV ERROR] %s\n" % (text))
        sys.stderr.flush()

    def info(self, text):
        sys.stderr.write("[GV INFO]  %s\n" % (text))
        sys.stderr.flush()

    def load_struct(self, path):
        """Loads the requested structure.
        """
        self.info("loading: %s" % (path))
        
        try:
            struct = LoadStructure(
                fil              = path,
                build_properties = ("library_bonds","distance_bonds"))
        except IOError:
            self.error("file not found: %s" % (path))
            return

        struct_desc = {}
        struct_desc["struct"] = struct
        struct_desc["path"] = path

        self.glv_add_struct(struct)

    def glv_render(self):
        self.glv_render_one(self.opengl_driver)
        
        if GL_DOUBLE_BUFFER:
            glutSwapBuffers()
        else:
            glFlush()

    def glv_redraw(self):
        if self.glut_init_done:
            glutPostRedisplay()

    def glut_init(self):
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
        #glutSpecialFunc(self.glut_special) 

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
        key = key.lower()
        
        ## quit
        if key=="q":
            sys.exit(0)

        ## toggle fullscreen mode:
        elif key=="o":
            pass

    def glut_special(self, key, x, y):
        pass
    

def main():
    gv = GLUT_Viewer()
    gv.glut_init()


    for path in sys.argv[1:]:
        gv.load_struct(path)
    
    gv.glut_main()


if __name__=="__main__":
    main()
