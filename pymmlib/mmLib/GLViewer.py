## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""OpenGL rendering classes.
"""
from __future__ import generators
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

        glRotatef(self.rotx,
                  self.axes[0,0],
                  self.axes[0,1],
                  self.axes[0,2])

        glRotatef(self.roty,
                  self.axes[1,0],
                  self.axes[1,1],
                  self.axes[1,2])

        glRotatef(self.rotz,
                  self.axes[2,0],
                  self.axes[2,1],
                  self.axes[2,2])

    def gl_pop_matrix(self):
        """Pop the roatated/translated position.
        """
        glPopMatrix()

    def gl_render(self):
        """Compile or force a recompile of this object's gl_draw list, and
        render the scene.  Rendering the scene can be bypassed if
        this method is called with render = False.
        """
        if self.show == False:
            return        

        if self.gl_name == None:
            self.gl_compile_list()

        self.gl_push_matrix()
        glCallList(self.gl_name)
        self.gl_pop_matrix()

    def gl_compile_list(self):
        """Compile a OpenGL draw list for this object.
        """
        if self.gl_name != None:
            self.gl_delete_list()
    
        self.gl_name = glGenLists(1)
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

    def gl_draw(self):
        glDisable(GL_LIGHTING)
        
        def axis_line(v1, v2):
            glLineWidth(self.axis_line_width)
            glBegin(GL_LINES)
            glVertex3f(*v1)
            glVertex3f(*v2)
            glEnd()
                    
        glColor3f(1.0, 0.0, 0.0)
        axis_line(Vector(0.0, 0.0, 0.0), Vector(200.0, 0.0, 0.0))

        glColor3f(0.0, 1.0, 0.0)
        axis_line(Vector(0.0, 0.0, 0.0), Vector(1.0, 200.0, 0.0))

        glColor3f(0.0, 0.0, 1.0)
        axis_line(Vector(0.0, 0.0, 0.0), Vector(0.0, 0.0, 200.0))

        glEnable(GL_LIGHTING)


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

        glDisable(GL_LIGHTING)
        glLineWidth(self.cell_line_width)
        glColor3f(1.0, 1.0, 1.0)

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

        glEnable(GL_LIGHTING)

    def gl_draw(self):
        self.draw_cell(-1, -1, -1, 0, 0, 0)


class GLTLSGroup(GLDrawList):
    """Draws TLS group
    """
    def __init__(self, **args):
        GLDrawList.__init__(self)

        self.tls_group = args["tls_group"] 
        self.calcs     = self.tls_group.calc_COR()
        self.origin    = Vector(self.calcs["COR"])

        
        (eval, evec) = eigenvectors(self.tls_group.L)
        self.axes    = evec
        self.evalL   = eval
        self.tm      = 0.0

        ## set up the atom list too
        self.gl_atom_list             = GLAtomList(U       = "Ucalc",
                                                   U_color = (0.,1.,0.),
                                                   color   = (1.,1.,1.),
                                                   line_width = 10.)
        self.gl_atom_list.origin      = self.origin
        self.axes                     = evec
        self.gl_atom_list.atom_origin = self.origin

        for atm, Ucalc in self.tls_group.iter_atm_Ucalc():
            atm.Ucalc = Ucalc
            self.gl_atom_list.append(atm)

    def gl_render(self):
        GLDrawList.gl_render(self)
        self.gl_atom_list.gl_render()

    def gl_delete_list(self):
        GLDrawList.gl_delete_list(self)
        self.gl_atom_list.gl_delete_list()
        
    def gl_draw(self):
        self.gl_atom_list.gl_draw()
        
        ## T axes
        glColor3f(0.0, 1.0, 0.0)
        self.draw_tensor(self.calcs["T'"])

        ## L axes
        glColor3f(1.0, 0.0, 0.0)
        self.draw_tensor(self.calcs["L'"], 0.05 * rad2deg2)

        ## S axes
        glColor3f(0.0, 0.0, 1.0)
        self.draw_tensor(self.calcs["S'"], rad2deg)

        return

        ## draw a line from the center of reaction to each atom
        self.set_material(1.0, 1.0, 1.0, 1.0)
        for atm in self.tls_group:
            if atm.name in ["CA"]:
                glBegin(GL_LINES)
                glVertex3f(0.0, 0.0, 0.0)
                glVertex3f(*atm.position - self.origin)
                glEnd()

    def draw_tensor(self, ten, scale = 1.0):
        """Draw tensor axis.
        """
        (eval, evec) = eigenvectors(ten)

        glDisable(GL_LIGHTING)
        glLineWidth(1.0)
        
        for i in range(3):
            v = scale * eval[i] * Vector(evec[i,0],evec[i,1],evec[i,2])
            glBegin(GL_LINES)
            glVertex3f(*-v)
            glVertex3f(*v)
            glEnd()

        glEnable(GL_LIGHTING)

    def inc_time(self):
        """val cycles every 2pi
        """
        return
        
        if self.gl_name == None:
            return
        
        self.tm += math.pi / 50.0

        try:
            Lx = self.evalL[0] * rad2deg2 * math.sin(self.tm) 
            Ly = self.evalL[1] * rad2deg2 * math.sin(self.tm)
            Lz = self.evalL[2] * rad2deg2 * math.sin(self.tm)
        except ValueError, err:
            print "inc_time: ",err
            return

##         Spc = self.calcs["S'^"] * rad2deg

##         dSp = array(
##             [ (Lx * Spc[0,0]) + (Ly * Spc[1,0]) + (Lz * Spc[2,0]),
##               (Lx * Spc[0,1]) + (Ly * Spc[1,1]) + (Lz * Spc[2,1]),
##               (Lx * Spc[0,2]) + (Ly * Spc[1,2]) + (Lz * Spc[2,2]) ])

##         dS = matrixmultiply(transpose(self.axes), dSp)

##         self.origin = Vector(self.calcs["COR"] + dS)

        self.rotx = Lx
        self.roty = Ly
        self.rotz = Lz

        self.gl_atom_list.origin = self.origin
        
        self.gl_atom_list.rotx = self.rotx
        self.gl_atom_list.roty = self.roty
        self.gl_atom_list.rotz = self.rotz


class GLAtomList(GLDrawList, AtomList):
    """OpenGL renderer for a list of atoms.  Optional arguments iare:
    color, U, U_color.
    """
    def __init__(
        self,
        sphere_quality = 12,
        line_width     = 3.0,
        atom_origin    = None,

        **args):

        GLDrawList.__init__(self)
        AtomList.__init__(self)

        ## draw settings
        self.args = args
        if not self.args.has_key("lines"):
            self.args["lines"] = True

        ## defaults
        self.sphere_quality = sphere_quality
        self.line_width     = line_width
        self.atom_origin    = atom_origin
        self.el_color_cache = {}

    def gl_draw(self):
        """Perform the OpenGL drawing operations to render this atom list
        with the current settings.
        """
        if self.args.get("lines", False):
            for atm in self:
                self.draw_lines(atm)

        if self.args.get("cpk", False):
            for atm in self:
                self.draw_cpk(atm)

        if self.args.get("U", False):
            for atm in self:
                self.draw_U_axes(atm)

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

    def draw_lines(self, atm, symop = None):
        """Draw a atom using bond lines only.
        """
        try:
            glColor3f(*self.el_color_cache[atm.element])
        except KeyError:
            if self.args.has_key("color"):
                cx = self.args["color"]
                self.el_color_cache[atm.element] = cx
                glColor3f(*cx)
            else:
                el = atm.get_structure().library.get_element(atm.element)
                if el != None:
                    self.el_color_cache[atm.element] = el.color
                    glColor3f(*el.color)
                else:
                    self.el_color_cache[atm.element] = (1.0, 1.0, 1.0)
                    glColor3f(1.0, 1.0, 1.0)

        if self.atom_origin:
            position = atm.position - self.atom_origin
        else:
            position = atm.position
        
        glLineWidth(self.line_width)
        glDisable(GL_LIGHTING)
        
        ## if there are bonds, then draw the lines 1/2 way to the
        ## bonded atoms
        if len(atm.bond_list) > 0:
            for bond in atm.iter_bonds():
                atm2 = bond.get_partner(atm)

                start = position
                end   = start + ((atm2.position - atm.position) / 2)

                glBegin(GL_LINES)
                glVertex3f(*start)
                glVertex3f(*end)
                glEnd()

        ## if there are no bonds, draw a small cross-point 
        else:
            vx = Vector(0.25, 0.0,  0.0)
            vy = Vector(0.0,  0.25, 0.0)
            vz = Vector(0.0,  0.0,  0.25)

            start = position - vx
            end   = position + vx
            glBegin(GL_LINES)
            glVertex3f(*start)
            glVertex3f(*end)
            glEnd()

            start = position - vy
            end   = position + vy
            glBegin(GL_LINES)
            glVertex3f(*start)
            glVertex3f(*end)
            glEnd()

            start = position - vz
            end   = position + vz
            glBegin(GL_LINES)
            glVertex3f(*start)
            glVertex3f(*end)
            glEnd()

        glEnable(GL_LIGHTING)

    def draw_U_axes(self, atm, symop = None):
        """Draw principal thermal axes for atom.
        """
        if atm.U == None:
            return
        
        if self.args["U"] == True:
            U = atm.U
        else:
            U = getattr(atm, self.args["U"])

        evec, eval = eigenvectors(U)

        v0 = Vector(evec[0] * eval[0])
        v1 = Vector(evec[1] * eval[1])
        v2 = Vector(evec[2] * eval[2])

        if self.atom_origin:
            position = atm.position - self.atom_origin
        else:
            position = atm.position

        try:
            glColor3f(*self.args["U_color"])
        except KeyError:
            glColor3f(1.0, 1.0, 1.0)

        glLineWidth(1.0)
        glDisable(GL_LIGHTING)
        
        glBegin(GL_LINES)
        glVertex3f(*position - v0)
        glVertex3f(*position + v0)
        glEnd()

        glBegin(GL_LINES)
        glVertex3f(*position - v1)
        glVertex3f(*position + v1)
        glEnd()

        glBegin(GL_LINES)
        glVertex3f(*position - v2)
        glVertex3f(*position + v2)
        glEnd()

        glEnable(GL_LIGHTING)

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
        self.gl_axes        = GLAxes()
        self.gl_unit_cell   = GLUnitCell(self.struct.unit_cell)
        self.aa_main_chain  = {}
        self.aa_side_chain  = {}
        self.dna_main_chain = {}
        self.dna_side_chain = {}
        self.water          = {}
        self.hetatm         = {}

        for chain in self.struct.iter_chains():
            aa_main_chain  = GLAtomList(U = True)
            aa_side_chain  = GLAtomList(U = True)
            dna_main_chain = GLAtomList()
            dna_side_chain = GLAtomList()
            water          = GLAtomList()
            hetatm         = GLAtomList()
            
            for frag in chain.iter_fragments():
                if isinstance(frag, AminoAcidResidue):
                    for atm in frag.iter_atoms():
                        if atm.name in ["C", "O", "CA", "N"]:
                            aa_main_chain.append(atm)
                        else:
                            aa_side_chain.append(atm)

                elif isinstance(frag, NucleicAcidResidue):
                    for atm in frag.iter_atoms():
                        dna_main_chain.append(atm)

                elif frag.is_water():
                    for atm in frag.iter_atoms(): 
                        water.append(atm)

                else:
                    for atm in frag.iter_atoms():
                        hetatm.append(atm)

            ## delete empty draw lists
            if len(aa_main_chain) > 0:
                self.aa_main_chain[chain.chain_id] = aa_main_chain
            if len(aa_side_chain) > 0:
                self.aa_side_chain[chain.chain_id] = aa_side_chain
            if len(dna_main_chain) > 0:
                self.dna_main_chain[chain.chain_id] = dna_main_chain
            if len(dna_side_chain) > 0:
                self.dna_side_chain[chain.chain_id] = dna_side_chain
            if len(water) > 0:
                self.water[chain.chain_id] = water
            if len(hetatm) > 0:
                self.hetatm[chain.chain_id] = hetatm

    def iter_draw_lists(self):
        """Iterate over all GL Lists.
        """
        yield self.gl_axes
        yield self.gl_unit_cell

        for gl_draw_list in self.aa_main_chain.values():
            yield gl_draw_list
        for gl_draw_list in self.aa_side_chain.values():
            yield gl_draw_list
        for gl_draw_list in self.dna_main_chain.values():
            yield gl_draw_list
        for gl_draw_list in self.dna_side_chain.values():
            yield gl_draw_list
        for gl_draw_list in self.water.values():
            yield gl_draw_list
        for gl_draw_list in self.hetatm.values():
            yield gl_draw_list

    def gl_render(self):
        self.gl_push_matrix()
        for draw_list in self.iter_draw_lists():
            draw_list.gl_render()
        self.gl_pop_matrix()

    def gl_draw(self):
        for draw_list in self.iter_draw_lists():
            draw_list.gl_draw()

    def gl_delete_list(self):
        """Delete all OpenGL draw lists.
        """
        for draw_list in self.iter_draw_lists():
            draw_list.gl_delete_list()

        self.gl_axes        = None
        self.gl_unit_cell   = None
        self.aa_main_chain  = {}
        self.aa_side_chain  = {}
        self.dna_main_chain = {}
        self.dna_side_chain = {}
        self.water          = {}
        self.hetatm         = {}

    def show_axes(self, show = None):
        """True/False to show the coordinate axes.
        """
        if show == None:
            return self.gl_axes.show
        else:
            self.gl_axes.show = show

    def show_unit_cell(self, show = None):
        """True/False to show the unit cell.
        """
        if show == None:
            return self.gl_unit_cell.show
        else:
            self.gl_unit_cell.show = show

    def show_chain_dict(self, show_dict, show, chain_id):
        if show == None:
            if chain_id == None:
                for gl_draw_list in show_dict.values():
                    if gl_draw_list.show == False:
                        return False
                return True
            else:
                return show_dict[chain_id].show
        else:
            if chain_id == None:
                for gl_draw_list in show_dict.values():
                    gl_draw_list.show = show
            else:
                show_dict[chain_id].show = show

    def show_aa_main_chain(self, show = None, chain_id = None):
        """Show/Hide amino acid main chain atoms, with optional specification
        of chain_id.
        """
        return self.show_chain_dict(self.aa_main_chain, show, chain_id)

    def show_aa_side_chain(self, show = None, chain_id = None):
        """Show/Hide amino acid side chain atoms, with optional specification
        of chain_id.
        """
        return self.show_chain_dict(self.aa_side_chain, show, chain_id)

    def show_dna_main_chain(self, show = None, chain_id = None):
        """Show/Hide nucleic acid main chain atoms, with optional specification
        of chain_id.
        """
        return self.show_chain_dict(self.dna_main_chain, show, chain_id)

    def show_dna_side_chain(self, show = None, chain_id = None):
        """Show/Hide nucleic acid side chain atoms, with optional specification
        of chain_id.
        """
        return self.show_chain_dict(self.dna_side_chain, show, chain_id)

    def show_water(self, show = None, chain_id = None):
        """Show/Hide water, with optional specification
        of chain_id.
        """
        return self.show_chain_dict(self.water, show, chain_id)

    def show_hetatm(self, show = None, chain_id = None):
        """Show/Hide het-group atoms, with optional specification
        of chain_id.
        """
        return self.show_chain_dict(self.hetatm, show, chain_id)


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
            draw_list.gl_render()
            
	if self.gldrawable.is_double_buffered():
            self.gldrawable.swap_buffers()
	else:
            glFlush()

        self.gldrawable.gl_end()
