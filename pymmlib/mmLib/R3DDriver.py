## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""Viewer.py graphics driver for producing a output file
for the Raster3D ray tracer.
"""
from __future__  import generators
import copy
import random

from Gaussian       import *
from Structure      import *


## constants
MARGIN = 1.20
BASE_LINE_WIDTH = 0.05


def matrixmultiply43(M, x):
    """M is a 4x4 rotation-translation matrix, x is a 3x1 vector.  Returns
    the 3x1 vector of x transformed by M.
    """
    x4    = array((x[0],x[1],x[2],1.0), Float)
    y     = matrixmultiply(M, x4)
    return array((y[0],y[1],y[2]), Float)


class Raster3DDriver(object):
    """Viewer.py graphics driver for producing a output file
    for the Raster3D ray tracer.
    """
    def __init__(self):
        self.glr_init_state()

    def glr_init_state(self):
        """Re-initalizes driver state variables.
        """
        self.matrix           = identity(4, Float)
        self.matrix_stack     = []

        self.line_width       = 1.0 * BASE_LINE_WIDTH

        self.normal           = None
        self.normalize        = False

        self.light_two_sides  = False
        
        self.width            = 400
        self.height           = 400
        self.bg_color_rgbf    = (0.0, 0.0, 0.0)
        self.ambient_light    = 0.2
        self.specular_light   = 1.0
        
        self.material_color_r = 1.0
        self.material_color_g = 1.0
        self.material_color_b = 1.0
        self.material_alpha   = 1.0

    def glr_clear_objects(self):
        """Clears out the current object list.
        """
        self.object_list = []

    def glr_compile_supported(self):
        """Returns True if draw compiling is supported by the driver.
        """
        return False
        
    def glr_render_begin(
        self,
        clear_object_list  = True,
        width              = 200,
        height             = 100,
        zoom               = 50,
        near               = 0,
        far                = 0,
        bg_color_rgbf      = (0.0, 0.0, 0.0),
        ambient_light      = 0.2,
        diffuse_light      = 1.0,
        specular_light     = 1.0,
        **args):
        """Sets up lighting and OpenGL options before scene rendering.
        """
        self.glr_init_state()

        if clear_object_list:
            self.glr_clear_objects()
        
        self.width            = width
        self.height           = height
        self.zoom             = zoom
        self.bg_color_rgbf    = bg_color_rgbf

        ## the lighting model for Raster3D is not quite the same as
        ## OpenGL; this conversion gets it close
        total_light = ambient_light + diffuse_light + specular_light
        
        self.ambient          = ambient_light  / total_light
        self.specular         = specular_light / total_light 
        self.phong            = 5

    def glr_construct_header(self):
        """Creates the header for the render program.
        """
        tsz_width   = 3
        tsz_height  = 3
        
        tile_width   = int(round(self.width  / float(tsz_width)))
        tile_height  = int(round(self.height / float(tsz_height)))

        pixel_width  = tile_width  * 3
        pixel_height = tile_height * 3

        ## self.zoom is the horizontal number of Angstroms shown in the
        ## image, this must be converted to the Raster3D zoom parameter
        ## which is the number of Angstroms of the shortest dimention
        if pixel_width>pixel_height:
            r = float(pixel_height) / float(pixel_width)
            z = self.zoom * r
        else:
            z = self.zoom

        self.header_list = [
            "mmLib Generated Raster3D Output", 
            "%d %d     tiles in x,y" % (tile_width, tile_height), 
            "%d %d     pixels (x,y) per tile" % (tsz_width, tsz_height),
            "4         anti-aliasing level 4; 3x3->2x2",
            "%4.2f %4.2f %4.2f background(rgb)" % self.bg_color_rgbf,
            "F 	       no shadows cast", 
            "%2d       Phong power" % (self.phong), 
            "0.20      secondary light contribution",
            "%4.2f     ambient light contribution" % (self.ambient),
            "%4.2f     specular reflection component" % (self.specular), 
            "0.0       eye position(no perspective)", 
            "1 1 1     main light source position", 
            "1 0 0 0   4x4 view matrix",
            "0 1 0 0", 
            "0 0 1 0",
            "0 0 0 %f" % (z),
            "3 	       mixed objects", 
            "* 	      (free format triangle and plane descriptors)",
            "* 	      (free format sphere descriptors",
            "* 	      (free format cylinder descriptors)",
            ]

    def glr_render_end(self):
        """Write out the input file for the render program.
        """
        ## open r3d file, write header
        fil = open("raster.r3d", "w")

        ## add required hader for the render program
        self.glr_construct_header()
        
        for ln in self.header_list:
            fil.write(ln + "\n")

        ## write objects
        for gob in self.object_list:

            gob_type = gob["type"]

            if gob_type=="cylinder":
                p1 = gob["position1"]
                p2 = gob["position2"]
                r  = gob["radius"]
                
                fil.write(
                    "5\n"\
                    "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f "\
                    "%4.2f %4.2f %4.2f\n" % (
                    p1[0], p1[1], p1[2],
                    r,
                    p2[0], p2[1], p2[2],
                    r,
                    gob["r"], gob["g"], gob["b"]))

            elif gob_type=="sphere":
                p = gob["position"]

                fil.write(
                    "2\n"\
                    "%8.3f %8.3f %8.3f %8.3f %4.2f %4.2f %4.2f\n" % (
                    p[0], p[1], p[2],
                    gob["radius"],
                    gob["r"], gob["g"], gob["b"]))
                
            elif gob_type=="triangle":
                v1 = gob["vertex1"]
                v2 = gob["vertex2"]
                v3 = gob["vertex3"]

                fil.write(
                    "1\n"\
                    "%8.3f %8.3f %8.3f "\
                    "%8.3f %8.3f %8.3f "\
                    "%8.3f %8.3f %8.3f "\
                    "%4.2f %4.2f %4.2f\n" % (
                     v1[0], v1[1], v1[2],
                     v2[0], v2[1], v2[2],
                     v3[0], v3[1], v3[2],
                     gob["r"], gob["g"], gob["b"]))

            elif gob_type=="normal":
                n1 = gob["normal1"]
                n2 = gob["normal2"]
                n3 = gob["normal3"]
                
                fil.write(
                    "7\n"\
                    "%8.3f %8.3f %8.3f "\
                    "%8.3f %8.3f %8.3f "\
                    "%8.3f %8.3f %8.3f\n" % (
                     n1[0], n1[1], n1[2],
                     n2[0], n2[1], n2[2],
                     n3[0], n3[1], n3[2]))

            elif gob_type=="ellipse":
                p = gob["position"]
                q = gob["quadric"]
                
                fil.write(
                    "14\n"\
                    "%8.3f %8.3f %8.3f "\
                    "%8.3f "\
                    "%4.2f %4.2f %4.2f "\
                    "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f 0 0 0 %8.3f\n" % (
                     p[0], p[1], p[2],
                     gob["limit_radius"],
                     gob["r"], gob["g"], gob["b"],
                     q[0,0], q[1,1], q[2,2], q[0,1], q[0,2], q[1,2],
                     gob["prob"]))

            elif gob_type=="material_properties":
                if gob["two_sided"]==True:
                    two_sided_flag = 2
                else:
                    two_sided_flag = 0

                fil.write(
                    "8\n"\
                    "-1 -1  -1 -1 -1  %4.2f  %1d 0 0 0\n" % (
                    gob["clrity"], two_sided_flag))

                
        fil.close()

    def glr_push_matrix(self):
        """
        """
        assert len(self.matrix_stack)<=25
        self.matrix_stack.append(self.matrix.copy())

    def glr_pop_matrix(self):
        """
        """
        self.matrix       = self.matrix_stack[-1]
        self.matrix_stack = self.matrix_stack[:-1]
        
    def glr_translate(self, t):
        """Translates the scene by vector t.
        """
        M = array( [[   1.0,    0.0,    0.0, t[0]],
                    [   0.0,    1.0,    0.0, t[1]],
                    [   0.0,    0.0,    1.0, t[2]],
                    [   0.0,    0.0,    0.0,  1.0]], Float)

        self.matrix = matrixmultiply(self.matrix, M)

    def glr_translate3(self, x, y, z):
        """
        """
        M = array( [[   1.0,    0.0,    0.0,    x],
                    [   0.0,    1.0,    0.0,    y],
                    [   0.0,    0.0,    1.0,    z],
                    [   0.0,    0.0,    0.0,  1.0]], Float)

        self.matrix = matrixmultiply(self.matrix, M)

    def glr_mult_matrix_Rt(self, R, t):
        """Return the current matrix as a 3x3 rotation matrix R and 3x1
        translation vector t.
        """
        M = array( [[R[0,0], R[0,1], R[0,2], t[0]],
                    [R[1,0], R[1,1], R[1,2], t[1]],
                    [R[2,0], R[2,1], R[2,2], t[2]],
                    [   0.0,    0.0,    0.0,  1.0]], Float)

        self.matrix = matrixmultiply(self.matrix, M)
        
    def glr_mult_matrix_R(self, R):
        """Multiplies the current matrix by rotation matrix R and translates
        by t
        """
        M = array( [[R[0,0], R[0,1], R[0,2], 0.0],
                    [R[1,0], R[1,1], R[1,2], 0.0],
                    [R[2,0], R[2,1], R[2,2], 0.0],
                    [   0.0,    0.0,    0.0, 1.0]], Float)

        self.matrix = matrixmultiply(self.matrix, M)
        
    def glr_rotate_axis(self, deg, axis):
        """
        """
        R = rmatrixu(axis, deg*DEG2RAD)
        self.glr_mult_matrix_R(R)

    def glr_lighting_enable(self):
        """
        """
        pass

    def glr_lighting_disable(self):
        """
        """
        pass

    def glr_set_line_width(self, width):
        """
        """
        self.line_width = width * BASE_LINE_WIDTH
    
    def glr_set_material_rgb(self, r, g, b):
        """Creates a stock rendering material colored according to the given
        RGB values.
        """
        self.material_color_r = r
        self.material_color_g = g
        self.material_color_b = b
        
        if self.material_alpha<1.0:
            self.material_alpha = 1.0

            self.object_list.append(
                {"type":       "material_properties",
                 "clrity":     0.0,
                 "two_sided":  self.light_two_sides })

    def glr_set_material_rgba(self, r, g, b, a):
        """Creates a stock rendering material colored according to the given
        RGB values.
        """
        self.material_color_r = r
        self.material_color_g = g
        self.material_color_b = b

        if self.material_alpha!=a:
            self.material_alpha = a

            self.object_list.append(
                {"type":       "material_properties",
                 "clrity":     1.0 - self.material_alpha,
                 "two_sided":  self.light_two_sides })

    def glr_vertex(self, vertex):
        """
        """
        self.glr_vertex_func(vertex)
        
    def glr_vertex3(self, x, y, z):
        """
        """
        self.glr_vertex_func((x,y,z))

    def glr_begin_lines(self):
        """
        """
        raise FinishMe()

    def glr_begin_triangles(self):
        """
        """
        raise FinishMe()

    def glr_begin_quads(self):
        self.glr_vertex_func  = self.glr_vertex_quads_1
        self.glr_end          = self.glr_end_quads
        self.vertex_1         = None
        self.normal_1         = None
        self.vertex_2         = None
        self.normal_2         = None
        self.vertex_3         = None
        self.normal_3         = None

    def glr_end_quads(self):
        del self.glr_vertex_func
        del self.glr_end
        del self.vertex_1
        del self.normal_1
        del self.vertex_2
        del self.normal_2
        del self.vertex_3
        del self.normal_3

    def glr_vertex_quads_1(self, vertex):
        self.glr_vertex_func = self.glr_vertex_quads_2
        self.normal_1        = self.normal
        self.vertex_1        = matrixmultiply43(self.matrix, vertex)

    def glr_vertex_quads_2(self, vertex):
        self.glr_vertex_func = self.glr_vertex_quads_3
        self.normal_2        = self.normal
        self.vertex_2        = matrixmultiply43(self.matrix, vertex)
        
    def glr_vertex_quads_3(self, vertex):
        self.glr_vertex_func = self.glr_vertex_quads_4
        self.normal_3        = self.normal
        self.vertex_3        = matrixmultiply43(self.matrix, vertex)

    def glr_vertex_quads_4_test(self, vertex):
        self.glr_vertex_func = self.glr_vertex_quads_1

        normal_4 = self.normal
        vertex_4 = matrixmultiply43(self.matrix, vertex)

        r = 0.05

        self.object_list.append(
            {"type":       "cylinder",
             "position1":  self.vertex_1,
             "position2":  self.vertex_2,
             "radius":     r,
             "r":          1.0,
             "g":          0.0,
             "b":          0.0 })
        self.object_list.append(
            {"type":       "cylinder",
             "position1":  self.vertex_2,
             "position2":  self.vertex_3,
             "radius":     r,
             "r":          0.0,
             "g":          1.0,
             "b":          0.0 })
        self.object_list.append(
            {"type":       "cylinder",
             "position1":  self.vertex_3,
             "position2":  vertex_4,
             "radius":     r,
             "r":          0.0,
             "g":          0.0,
             "b":          1.0 })
        self.object_list.append(
            {"type":       "cylinder",
             "position1":  vertex_4,
             "position2":  self.vertex_1,
             "radius":     r,
             "r":          1.0,
             "g":          1.0,
             "b":          1.0 })
        
        self.object_list.append(
            {"type":       "cylinder",
             "position1":  self.vertex_1,
             "position2":  self.vertex_1 + self.normal_1,
             "radius":     r,
             "r":          1.0,
             "g":          0.0,
             "b":          0.0 })
        self.object_list.append(
            {"type":       "cylinder",
             "position1":  self.vertex_2,
             "position2":  self.vertex_2 + self.normal_2,
             "radius":     r,
             "r":          0.0,
             "g":          1.0,
             "b":          0.0 })
        self.object_list.append(
            {"type":       "cylinder",
             "position1":  self.vertex_3,
             "position2":  self.vertex_3 + self.normal_3,
             "radius":     r,
             "r":          0.0,
             "g":          0.0,
             "b":          1.0 })
        self.object_list.append(
            {"type":       "cylinder",
             "position1":  vertex_4,
             "position2":  vertex_4 + normal_4,
             "radius":     r,
             "r":          1.0,
             "g":          1.0,
             "b":          1.0 })

    def glr_vertex_quads_4(self, vertex):
        self.glr_vertex_func = self.glr_vertex_quads_1

        normal_4 = self.normal
        vertex_4 = matrixmultiply43(self.matrix, vertex)
 
        self.object_list.append(
            {"type":      "triangle",
             "vertex1":   self.vertex_1,
             "vertex2":   self.vertex_2,
             "vertex3":   self.vertex_3,
             "r":         self.material_color_r,
             "g":         self.material_color_g,
             "b":         self.material_color_b })

        self.object_list.append(
            {"type":      "normal",
             "normal1":   self.normal_1,
             "normal2":   self.normal_2,
             "normal3":   self.normal_3 })

        self.object_list.append(
            {"type":      "triangle",
             "vertex1":   self.vertex_1,
             "vertex2":   self.vertex_3,
             "vertex3":   vertex_4,
             "r":         self.material_color_r,
             "g":         self.material_color_g,
             "b":         self.material_color_b })

        self.object_list.append(
            {"type":      "normal",
             "normal1":   self.normal_1,
             "normal2":   self.normal_3,
             "normal3":   normal_4 })

    def glr_begin_triangle_fan(self):
        """
        """
        self.glr_vertex_func  = self.glr_vertex_triangle_fan_1
        self.glr_end          = self.glr_end_triangle_fan
        self.vertex_1         = None
        self.vertex_2         = None
        self.normal_1         = None
        self.normal_2         = None
    
    def glr_end_triangle_fan(self):
        """
        """
        del self.glr_vertex_func
        del self.glr_end
        del self.vertex_1
        del self.vertex_2
        del self.normal_1
        del self.normal_2
        
    def glr_vertex_triangle_fan_1(self, vertex):
        """Get (first) common fan vertex.
        """
        self.glr_vertex_func = self.glr_vertex_triangle_fan_2
        
        self.vertex_1 = matrixmultiply43(self.matrix, vertex)
        self.normal_1 = self.normal

    def glr_vertex_triangle_fan_2(self, vertex):
        """Get second vertex.
        """
        self.glr_vertex_func = self.glr_vertex_triangle_fan_3

        self.vertex_2 = matrixmultiply43(self.matrix, vertex)
        self.normal_2 = self.normal

    def glr_vertex_triangle_fan_3(self, vertex):
        """Get third vertex and beyond: construct triangles.
        """
        vertex_3 = matrixmultiply43(self.matrix, vertex)
        normal_3 = self.normal

        self.object_list.append(
            {"type":      "triangle",
             "vertex1":   self.vertex_1,
             "vertex2":   self.vertex_2,
             "vertex3":   vertex_3,
             "r":         self.material_color_r,
             "g":         self.material_color_g,
             "b":         self.material_color_b })

        self.object_list.append(
            {"type":      "normal",
             "normal1":   self.normal_1,
             "normal2":   self.normal_2,
             "normal3":   normal_3 })

        self.vertex_2 = vertex_3
        self.normal_2 = normal_3

    def glr_normalize_enable(self):
        """
        """
        self.normalize = True

    def glr_normalize_disable(self):
        """
        """
        self.normalize = False

    def glr_normal(self, n):
        """
        """
        ## just rotate the normal
        R  = self.matrix[:3,:3]
        nr = matrixmultiply(R, n)

        if self.normalize==True:
            self.normal = normalize(nr)
        else:
            self.normal = nr

    def glr_normal3(self, x, y, z):
        """
        """
        ## just rotate the normal
        R  = self.matrix[:3,:3]
        nr = matrixmultiply(R, n)

        if self.normalize==True:
            self.normal = normalize(nr)
        else:
            self.normal = nr
        
    def glr_light_two_sides_enable(self):
        """
        """
        self.light_two_sides = True

        self.object_list.append(
            {"type":       "material_properties",
             "clrity":     1.0 - self.material_alpha,
             "two_sided":  True })
        
    def glr_light_two_sides_disable(self):
        """
        """
        self.light_two_sides = False

        self.object_list.append(
            {"type":       "material_properties",
             "clrity":     1.0 - self.material_alpha,
             "two_sided":  False })
        
    def glr_line(self, position1, position2):
        """Draws a single line.
        """
        self.object_list.append(
            {"type":       "cylinder",
             "position1":  matrixmultiply43(self.matrix, position1),
             "position2":  matrixmultiply43(self.matrix, position2),
             "radius":     self.line_width,
             "r":          self.material_color_r,
             "g":          self.material_color_g,
             "b":          self.material_color_b })
            
    def glr_text(self, text, scale):
        """Renders a text string.
        """
        pass
            
    def glr_axis(self, position, axis, radius):
        """Draw a vector axis using the current set material at position
        with the given radius.
        """
        ## don't bother redering small axes -- they look like junk
        if allclose(length(axis), 0.0):
            return

        self.object_list.append(
            {"type":       "cylinder",
             "position1":  matrixmultiply43(self.matrix, position),
             "position2":  matrixmultiply43(self.matrix, position + axis),
             "radius":     radius,
             "r":          self.material_color_r,
             "g":          self.material_color_g,
             "b":          self.material_color_b })

    def glr_tube(self, position1, position2, radius):
        """Draws a hollow tube beginning at pos1, and ending at pos2.
        """
        self.object_list.append(
            {"type":       "cylinder",
             "position1":  matrixmultiply43(self.matrix, position1),
             "position2":  matrixmultiply43(self.matrix, position2),
             "radius":     radius,
             "r":          self.material_color_r,
             "g":          self.material_color_g,
             "b":          self.material_color_b })

    def glr_sphere(self, position, radius, quality):
        """Draw a atom as a CPK sphere.
        """
        self.object_list.append(
            {"type":       "sphere",
             "position":   matrixmultiply43(self.matrix, position),
             "radius":     radius,
             "r":          self.material_color_r,
             "g":          self.material_color_g,
             "b":          self.material_color_b })

    def glr_cross(self, position, color, line_width):
        """Draws atom with a cross of lines.
        """
        pass

    def glr_Uaxes(self, position, U, prob, color, line_width):
        """Draw the anisotropic axies of the atom at the given probability.
        """
        pass
        
    def glr_Uellipse(self, position, U, prob):
        """Renders the ellipsoid enclosing the given fractional probability
        given the gaussian variance-covariance matrix U at the given position.
        C=1.8724 = 68%
        """
        ## rotate U
        R  = self.matrix[:3,:3]
        Ur = matrixmultiply(matrixmultiply(R, U), transpose(R))
        
        max_eigenvalue = max(eigenvalues(Ur))
        limit_radius = GAUSS3C[prob] * MARGIN * math.sqrt(max_eigenvalue)

        try:
            quadric = inverse(Ur)
        except LinAlgError:
            return
        
        self.object_list.append(
            {"type":         "ellipse",
             "position":     matrixmultiply43(self.matrix, position),
             "limit_radius": limit_radius,
             "quadric":      quadric,
             "prob":         -GAUSS3C[prob]**2,
             "r":            self.material_color_r,
             "g":            self.material_color_g,
             "b":            self.material_color_b })
    
    def glr_Urms(self, position, U):
        """Renders the root mean square (one standard deviation) surface of
        the
        gaussian variance-covariance matrix U at the given position.  This
        is a peanut-shaped surface. (Note: reference the peanut paper!)
        """
        pass

