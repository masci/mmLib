## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""OpenGL rendering classes.
"""
from __future__  import generators

import copy

from OpenGL.GL      import *
from OpenGL.GLU     import *
from OpenGL.GLUT    import *
from Gaussian       import *
from Structure      import *


try:
    import glaccel
except ImportError:
    GLACCEL_EXISTS = False
else:
    GLACCEL_EXISTS = True


class OpenGLDriver(object):
    """OpenGL render driver for Viewer.py
    """
    def __init__(self):
        self.gl_draw_list_id            = {}
        self.gl_draw_list_expiration_id = {}

    def glr_compile_supported(self):
        """Returns True if draw compiling is supported by the driver.
        """
        return True

    def glr_compile_start(self, draw_method):
        """
        """
        mid = id(draw_method)
        expiration_id = draw_method["expiration_id"]
        assert mid not in self.gl_draw_list_id

        draw_list_id = glGenLists(1)
        self.gl_draw_list_id[mid] = draw_list_id
        self.gl_draw_list_expiration_id[mid] = expiration_id
        glNewList(draw_list_id, GL_COMPILE)

        return mid

    def glr_compile_end(self):
        """
        """
        glEndList()
        
    def glr_compile_delete(self, draw_method):
        """
        """
        mid = id(draw_method)
        assert mid in self.gl_draw_list_id

        draw_list_id = self.gl_draw_list_id[mid] 
        glDeleteLists(draw_list_id, 1)
        del self.gl_draw_list_id[mid]
        del self.gl_draw_list_expiration_id[mid]

    def glr_compile_exists(self, draw_method):
        """
        """
        mid = id(draw_method)
        return mid in self.gl_draw_list_id

    def glr_compile_current(self, draw_method):
        """
        """
        mid = id(draw_method)
        expiration_id = draw_method["expiration_id"]
        assert mid in self.gl_draw_list_id
        return self.gl_draw_list_expiration_id[mid]==expiration_id
    
    def glr_compile_render(self, draw_method):
        """
        """
        mid = id(draw_method)
        assert mid in self.gl_draw_list_id

        draw_list_id = self.gl_draw_list_id[mid] 
        glCallList(draw_list_id)
        
    def glr_render_begin(
        self,
        width              = 200,
        height             = 100,
        zoom               = 50,
        near               = 0,
        far                = 0,
        bg_color_rgbf      = (0.0, 0.0, 0.0),
        ambient_light      = 1.0,
        diffuse_light      = 1.0,
        specular_light     = 1.0,
        gl_line_smooth     = False,
        gl_point_smooth    = False,
        gl_polygon_smooth  = False,
        gl_blend           = True,
        gl_fog             = False,
        gl_fog_start       = 0.0,
        gl_fog_end         = 0.0,
        **args):
        """Sets up lighting and OpenGL options before scene rendering.
        """
        
        ## setup vieweport
        glViewport(0, 0, width, height)

        ## setup perspective matrix
	glMatrixMode(GL_PROJECTION)
 	glLoadIdentity()

        zoom  = zoom / 2.0
        ratio = float(height) / float(width)
        glOrtho(-zoom, zoom, -ratio*zoom, ratio*zoom, -near, -far)
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        
        ## OpenGL Features
        glEnable(GL_NORMALIZE)
        glShadeModel(GL_SMOOTH)
        glEnable(GL_DEPTH_TEST)
        glDepthFunc(GL_LESS)
        
        ## background color
        glClearColor(bg_color_rgbf[0],
                     bg_color_rgbf[1],
                     bg_color_rgbf[2],
                     0.0)
        
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT)
            
        ## lighting
        ambient  = (ambient_light, ambient_light, ambient_light, 1.0)
        diffuse  = (diffuse_light, diffuse_light, diffuse_light, 1.0)
        specular = (specular_light, specular_light, specular_light, 1.0)

        ## use model abient light
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient)

        ## light 0
        glEnable(GL_LIGHT0)
        glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse)
        glLightfv(GL_LIGHT0, GL_SPECULAR, specular)
        glLightfv(GL_LIGHT0, GL_POSITION, (0.0, 0.0, 10.0, 0.0))
        
        ## light 1
        glDisable(GL_LIGHT1)
        glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse)
        glLightfv(GL_LIGHT1, GL_SPECULAR, specular)
        glLightfv(GL_LIGHT1, GL_POSITION, (0.0, 0.0, 100.0, 0.0))

        ## ANTI-ALIASING
        if gl_line_smooth:
            glEnable(GL_LINE_SMOOTH)
        else:
            glDisable(GL_LINE_SMOOTH)

        if gl_point_smooth:
            glEnable(GL_POINT_SMOOTH)
        else:
            glDisable(GL_POINT_SMOOTH)

        if gl_polygon_smooth:
            glEnable(GL_POLYGON_SMOOTH)
        else:
            glDisable(GL_POLYGON_SMOOTH)

        ## ALPHA BLENDING
        if gl_blend:
            glEnable(GL_BLEND)
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
            
            ## FOG
            if gl_fog:
                glEnable(GL_FOG)
                glFogf(GL_FOG_MODE,    GL_LINEAR)
                glFogf(GL_FOG_START,   gl_fog_start)
                glFogf(GL_FOG_END,     gl_fog_end)
            else:
                glDisable(GL_FOG)
        else:
            glDisable(GL_BLEND)

    def glr_render_end(self):
        """
        """
        pass

    def glr_push_matrix(self):
        """
        """
        glPushMatrix()

    def glr_pop_matrix(self):
        """
        """
        glPopMatrix()
        
    def glr_translate(self, t):
        """Translates the scene by vector t.
        """
        glTranslatef(*t)

    def glr_translate3(self, x, y, z):
        """
        """
        glTranslatef(x, y, z)
    
    def glr_mult_matrix_Rt(self, R, t):
        """Return the current matrix as a 3x3 rotation matrix R and 3x1
        translation vector t.
        """
        ## OpenGL wants the matrix in column-major form
        glMultMatrixf(
            (R[0,0], R[1,0], R[2,0], 0.0,
             R[0,1], R[1,1], R[2,1], 0.0,
             R[0,2], R[1,2], R[2,2], 0.0,
             t[0],   t[1],   t[2],   1.0) )
        
    def glr_mult_matrix_R(self, R):
        """Multiplies the current matrix by rotation matrix R and translates
        by t
        """
        ## OpenGL wants the matrix in column-major form
        glMultMatrixf(
            (R[0,0], R[1,0], R[2,0], 0.0,
             R[0,1], R[1,1], R[2,1], 0.0,
             R[0,2], R[1,2], R[2,2], 0.0,
             0.0,    0.0,    0.0,    1.0) )

    def glr_rotate_axis(self, deg, axis):
        """
        """
        glRotatef(deg, *axis)

    def glr_lighting_enable(self):
        """
        """
        glEnable(GL_LIGHTING)

    def glr_lighting_disable(self):
        """
        """
        glDisable(GL_LIGHTING)

    def glr_set_line_width(self, width):
        """
        """
        glLineWidth(width)

    def glr_begin_lines(self):
        """
        """
        glBegin(GL_LINES)

    def glr_begin_triangles(self):
        """
        """
        self.vertex_draw_type = "triangles"

    def glr_begin_triangle_fan(self):
        """
        """
        glBegin(GL_TRIANGLE_FAN)
        
    def glr_begin_quads(self):
        """
        """
        glBegin(GL_QUADS)
        
    def glr_end(self):
        """
        """
        glEnd()
    
    def glr_vertex(self, position):
        """
        """
        glVertex3f(*position)

    def glr_vertex3(self, x, y, z):
        """
        """
        glVertex3f(x, y, z)

    def glr_normalize_enable(self):
        """
        """
        glEnable(GL_NORMALIZE)

    def glr_normalize_disable(self):
        """
        """
        glDisable(GL_NORMALIZE)

    def glr_normal(self, n):
        """
        """
        glNormal3f(*n)

    def glr_normal3(self, x, y, z):
        """
        """
        glNormal3f(x, y, z)

    def glr_light_two_sides_enable(self):
        """
        """
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE)

    def glr_light_two_sides_disable(self):
        """
        """
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE)
    
    def glr_set_material_rgb(self, r, g, b):
        """Creates a stock rendering material colored according to the given
        RGB values.
        """
        glColor3f(r, g, b)
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,(r, g, b, 1.0))
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, (1.0, 1.0, 1.0, 1.0))
        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, (0.0, 0.0, 0.0, 1.0))
        glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, 100.0)

    def glr_set_material_rgba(self, r, g, b, a):
        """Creates a stock rendering material colored according to the given
        RGB values.
        """
        glColor3f(r, g, b)

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, (r, g, b, a))
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, (1.0, 1.0, 1.0, 1.0))
        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, (0.0, 0.0, 0.0, 1.0))

        if a<1.0:
            glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, 128.0)
        else:
            glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, 100.0)

    def glr_begin_lines(self):
        """
        """
        glBegin(GL_LINES)

    def glr_end(self):
        """
        """
        glEnd()

    def glr_vertex(self, position):
        """
        """
        glVertex3f(*position)

    def glr_vertex3(self, x, y, z):
        """
        """
        glVertex3f(x, y, z)

    def glr_line(self, position1, position2):
        """Draws a single line.
        """
        glBegin(GL_LINES)
        glVertex3f(*position1)
        glVertex3f(*position2)
        glEnd()

    def glr_text(self, text, scale):
        """Renders a text string.
        """
        glDisable(GL_LIGHTING)
        glLineWidth(2.0)

        glPushMatrix()
        
        s = scale / 1000.0
        glScalef(s, s, s)

        for c in text:
            glutStrokeCharacter(GLUT_STROKE_ROMAN, ord(c))

        glPopMatrix()
            
    def glr_axis(self, position, axis, radius):
        """Draw a vector axis using the current set material at position
        with the given radius.
        """
        ## don't bother redering small axes -- they look like junk
        if allclose(length(axis), 0.0):
            return
        
        end = position + axis

        glaccel.rod(
            position[0], position[1], position[2],
            end[0], end[1], end[2],
            radius)

    def glr_tube(self, pos1, pos2, radius):
        """Draws a hollow tube beginning at pos1, and ending at pos2.
        """
        try:
            glaccel.tube(
                pos1[0], pos1[1], pos1[2],
                pos2[0], pos2[1], pos2[2],
                radius)

        except NameError:
            glDisable(GL_LIGHTING)
            glLineWidth(1.0)
            glBegin(GL_LINES)
            glVertex3f(*pos1)
            glVertex3f(*pos2)
            glEnd()

    def glr_sphere(self, position, radius, quality):
        """Draw a atom as a CPK sphere.
        """
        try:
            glaccel.sphere(position[0], position[1], position[2],
                           radius, quality)

        except NameError:
            pass

    def glr_cross(self, position, color, line_width):
        """Draws atom with a cross of lines.
        """
        glDisable(GL_LIGHTING)
        glColor3f(*color)
        glLineWidth(line_width)

        vx = array([0.25, 0.0,  0.0])
        vy = array([0.0,  0.25, 0.0])
        vz = array([0.0,  0.0,  0.25])

        glBegin(GL_LINES)
        
        start = position - vx
        end   = position + vx
        glVertex3f(*start)
        glVertex3f(*end)
        
        start = position - vy
        end   = position + vy
        glVertex3f(*start)
        glVertex3f(*end)
        
        start = position - vz
        end   = position + vz
        glVertex3f(*start)
        glVertex3f(*end)

        glEnd()

    def glr_Uaxes(self, position, U, prob, color, line_width):
        """Draw the anisotropic axies of the atom at the given probability.
        """
        C = GAUSS3C[prob]        
        eval_U, evec_U = eigenvectors(U)

        try:
            v0_peak = C * math.sqrt(eval_U[0])
        except ValueError:
            v0_peak = 0.0
        try:
            v1_peak = C * math.sqrt(eval_U[1])
        except ValueError:
            v1_peak = 0.0
        try:
            v2_peak = C * math.sqrt(eval_U[2])
        except ValueError:
            v2_peak = 0.0
        
        v0 = evec_U[0] * v0_peak
        v1 = evec_U[1] * v1_peak
        v2 = evec_U[2] * v2_peak

        glDisable(GL_LIGHTING)
        glColor3f(*color)
        glLineWidth(line_width)

        glPushMatrix()
        glTranslatef(*position)

        glBegin(GL_LINES)
        glVertex3f(*-v0)
        glVertex3f(* v0)
        glVertex3f(*-v1)
        glVertex3f(* v1)
        glVertex3f(*-v2)
        glVertex3f(* v2)
        glEnd()
        
        glPopMatrix()

    def glr_Uellipse(self, position, U, prob):
        """Renders the ellipsoid enclosing the given fractional probability
        given the gaussian variance-covariance matrix U at the given position.
        C=1.8724 = 68%
        """
        glaccel.Uellipse(
            position[0], position[1], position[2],
            U[0,0], U[1,1], U[2,2], U[0,1], U[0,2], U[1,2],
            GAUSS3C[prob], 3)
        
    def glr_Urms(self, position, U):
        """Renders the root mean square (one standard deviation) surface of
        the gaussian variance-covariance matrix U at the given position.  This
        is a peanut-shaped surface. (Note: reference the peanut paper!)
        """
        glaccel.Upeanut(
            position[0], position[1], position[2],
            U[0,0], U[1,1], U[2,2], U[0,1], U[0,2], U[1,2],
            3)
        
