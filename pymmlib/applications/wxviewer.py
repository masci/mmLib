#!/usr/bin/env python

## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys

import wx
import wx.xrc
import wx.glcanvas

from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

from mmLib.PDB            import PDBFile
from mmLib.Structure      import *
from mmLib.FileLoader     import LoadStructure, SaveStructure
from mmLib.GLViewer       import *


class wxGLViewer(wx.glcanvas.GLCanvas):
    def __init__(self, parent):
        self.init = False
        self.mmlib_glviewer = GLViewer()
        self.x = 0
        self.y = 0
        
        wx.glcanvas.GLCanvas.__init__(self, parent, -1)
        
        wx.EVT_ERASE_BACKGROUND(self, self.on_erase_background)
        wx.EVT_SIZE(self, self.on_size)
        wx.EVT_PAINT(self, self.on_paint)
        wx.EVT_LEFT_DOWN(self, self.on_mouse_down)
        wx.EVT_MIDDLE_DOWN(self, self.on_mouse_down)
        wx.EVT_RIGHT_DOWN(self, self.on_mouse_down)
        wx.EVT_LEFT_UP(self, self.on_mouse_up)
        wx.EVT_MIDDLE_UP(self, self.on_mouse_up)
        wx.EVT_RIGHT_UP(self, self.on_mouse_up)
        wx.EVT_MOTION(self, self.on_mouse_motion)

    def on_erase_background(self, event):
        """Prevents blinking on Microsoft Windows.
        """
        pass

    def on_size(self, event):
        width, height = self.GetClientSize()

        if not self.GetContext():
            return
        
        self.SetCurrent()
        self.mmlib_glviewer.gl_resize(width, height)
        self.SwapBuffers()

    def on_paint(self, event):
        if not self.GetContext():
            return
        
        self.SetCurrent()

        if self.init==False:
            self.mmlib_glviewer.gl_init()
            self.init = True

            w, h = self.GetClientSize()
            self.mmlib_glviewer.gl_resize(w, h)
        
        self.mmlib_glviewer.gl_render()
        self.SwapBuffers()

    def on_mouse_down(self, event):
        self.CaptureMouse()
        self.beginx, self.beginy = event.GetPosition()

    def on_mouse_up(self, event):
        self.ReleaseMouse()

    def on_mouse_motion(self, event):
        if event.Dragging()==False:
            return

        x = 0.0
        y = 0.0
        z = 0.0

        evx,   evy    = event.GetPosition()
        width, height = self.GetClientSize()
            
        if event.LeftIsDown():
            self.mmlib_glviewer.roty += 360.0 * ((evx - self.beginx) / float(width)) 
            self.mmlib_glviewer.rotx += 360.0 * ((evy - self.beginy) / float(height))

        elif event.MiddleIsDown():
            z = 50.0 * ((evy - self.beginy) / float(height))

        elif event.RightIsDown():
            x = 50.0 * ( (evx - self.beginx) / float(height))
            y = 50.0 * (-(evy - self.beginy) / float(height))

        self.mmlib_glviewer.gl_translate(x, y, z)
        self.beginx, self.beginy = evx, evy

        self.Refresh(False)


class AppFrame(wx.Frame):
    """Frame class.
    """
    def __init__(self, parent=None, id=-1, title='Title', pos=wx.DefaultPosition, size=(400, 200)):
        """Create a Frame instance.
        """
        wx.Frame.__init__(self, parent, id, title, pos, size)

        ## <menu bar>
        menu_bar = wx.MenuBar()

        file_menu = wx.Menu()
        file_menu.Append(101, "&Quit")
        menu_bar.Append(file_menu, "&File")
        wx.EVT_MENU(self, 101, self.evt_file_quit)

        self.SetMenuBar(menu_bar)
        ## </menu bar>
        
        self.glviewer = wxGLViewer(self)
        self.load_structure(sys.argv[1])

    def evt_file_quit(self, event):
        self.Close()
        
    def load_structure(self, path):
        struct = LoadStructure(
            fil              = path,
            update_cb        = self.update_cb,
            build_properties = ("sequence","bonds"))
        struct.path = path
        self.add_struct(struct)
        return struct

    def update_cb(self, percent):
        pass

    def add_struct(self, struct):
        """Adds a structure to this viewer, and returns the GLStructure
        object so it can be manipulated.
        """
        gl_struct = GLStructure(struct=struct)
        self.glviewer.mmlib_glviewer.add_draw_list(gl_struct)


class App(wx.App):
    """Application class.
    """
    def OnInit(self):
        self.frame = AppFrame()
        self.frame.Show()
        self.SetTopWindow(self.frame)
        return True

def main():
    app = App()
    app.MainLoop()

if __name__ == '__main__':
    main()


