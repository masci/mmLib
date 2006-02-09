## TLS Minimized Domains (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

###############################################################################
## Report Directory Generation
##

import os
import popen2
import time
import numpy

## Python Imaging Library imports
import Image
import ImageDraw

## mmLib
from mmLib.mmTypes import *
from mmLib import Colors, Viewer, R3DDriver, Structure, Gaussian, FileLoader
from mmLib.Extensions import TLS

## tlsmdlib
import misc
import const
import conf
import gnuplots
import sequence_plot

from captions             import *
from tls_animate          import TLSAnimate, TLSAnimateFailure


## JMol URL path
JMOL_DIR = "../../../jmol"

## the pixel width of the TLS visualization rendered ray traces
VIS_WIDTH = 550

## the JMol viewer is a square window, generated with this pixel size
JMOL_SIZE = 600


def calc_inertia_tensor(atom_iter):
    """Calculate moment of inertia tensor at the centroid
    of the atoms.
    """
    al              = Structure.AtomList(atom_iter)
    centroid        = al.calc_centroid()

    I = numpy.zeros((3,3), float)
    for atm in al:
        x = atm.position - centroid

        I[0,0] += x[1]**2 + x[2]**2
        I[1,1] += x[0]**2 + x[2]**2
        I[2,2] += x[0]**2 + x[1]**2

        I[0,1] += - x[0]*x[1]
        I[1,0] += - x[0]*x[1]
        
        I[0,2] += - x[0]*x[2]
        I[2,0] += - x[0]*x[2]

        I[1,2] += - x[1]*x[2]
        I[2,1] += - x[1]*x[2]

    evals, evecs = numpy.linalg.eigenvectors(I)

    elist = [(evals[0], evecs[0]),
             (evals[1], evecs[1]),
             (evals[2], evecs[2])]

    elist.sort()

    R = numpy.array((elist[0][1], elist[1][1], elist[2][1]), float)

    ## make sure the tensor uses a right-handed coordinate system
    if numpy.allclose(numpy.linalg.determinant(R), -1.0):
        I = numpy.identity(3, float)
        I[0,0] = -1.0
        R = numpy.matrixmultiply(I, R)
    assert numpy.allclose(numpy.linalg.determinant(R), 1.0)

    return centroid, R


def calc_orientation(struct, chain):
    """Orient the structure based on a moment-of-inertia like tensor
    centered at the centroid of the structure.
    """
    ori = {}

    def iter_protein_atoms(sobjx):
        for fragx in sobjx.iter_amino_acids():
            for atmx in fragx.iter_atoms():
                if atmx.name in ["N","CA","C","O","CB"]: yield atmx

    centroids, Rs = calc_inertia_tensor(iter_protein_atoms(struct))
    centroidc, Rc = calc_inertia_tensor(iter_protein_atoms(chain))

    R = Rs
    centroid = centroidc

    ## now calculate a rectangular box
    first_atm = True

    min_x = 0.0
    max_x = 0.0
    min_y = 0.0
    max_y = 0.0
    min_z = 0.0
    max_z = 0.0

    for atm in iter_protein_atoms(chain):
        x  = numpy.matrixmultiply(R, atm.position - centroid)

        if first_atm==True:
            first_atm = False

            min_x = max_x = x[0]
            min_y = max_y = x[1]
            min_z = max_z = x[2]
        else:
            min_x = min(min_x, x[0])
            max_x = max(max_x, x[0])
            min_y = min(min_y, x[1])
            max_y = max(max_y, x[1])
            min_z = min(min_z, x[2])
            max_z = max(max_z, x[2])

    ## border around the structure
    border = 2.0
    min_x -= border
    max_x += border
    min_y -= border
    max_y += border
    min_z -= border
    max_z += border

    ## determine the center of the image as a xy offset from the centroid
    width  = max_x - min_x
    height = max_y - min_y
    zoom   = width

    xcenter = min_x + (width / 2.0)
    ycenter = min_y + (height / 2.0)

    xydelta  = numpy.array((xcenter, ycenter, 0.0), float)

    pwidth = VIS_WIDTH
    pheight = pwidth * (height / width)

    ori["R"]        = R
    ori["centroid"] = centroid + numpy.matrixmultiply(numpy.transpose(R), xydelta)
    ori["pwidth"]   = pwidth
    ori["pheight"]  = pheight 
    ori["hzoom"]    = zoom

    ## calculate near, far clipping plane
    ori["near"] = max_z
    ori["far"]  = min_z
    
    return ori


def html_tls_group_table(chainopt, tlsopt, ntls, report_root=None):
    """Generate HTML for a table containing the details of the ntls-group partitioning of the given chain.
    """
    ## inspect the first tls group dictionary to determine TLS model type
    try:
        tls = tlsopt.tls_list[0]
    except IndexError:
        return ""

    tls_model = tls["tls_group"].model
    
    if tls_model=="ISOT":
        t_head = 'T<sup>r</sup> <var>B</var>'

    elif tls_model=="ANISO":
        t_head = 'eval(T<sup>r</sup>) <var>B</var>'

    else:
        return ""

    l = ['<table width="100%" border=0 style="background-color:#eeeeee; font-size:x-small">',
         '<tr>',
         '<th align="center" colspan="12">Analysis of TLS Group Chain Segments</th>',
         '</tr>',

         '<tr>',
         '<th colspan="6" style="background-color:#aaaaaa">Input Structure</th>',
         '<th colspan="5" style="background-color:#bbbbbb">TLS Predictions</th>',
         '</tr>',
         
         '<tr style="background-color:#bbbbbb">',
         '<th>Color</th>',
         '<th>Segment</th>',
         '<th>Residues</th>',
         '<th>Atoms</th>',
         '<th>&#60;B&#62;</th>',
         '<th>&#60;Aniso&#62;</th>',
         '<th>RMSD B</th>',
         '<th>%s</th>' % (t_head),
         '<th>eval(L) <var>DEG<sup>2</sup></var></th>',
         '<th>&#60;B&#62;</th>',
         '<th>&#60;Aniso&#62;</th>',
         '</tr>' ]

    bgcolor_flag = True

    for tls in tlsopt.tls_list:
        
        tls_group = tls["tls_group"]
        tls_info  = tls["tls_info"]
        mtls_info = tls["model_tls_info"]

        L1 = mtls_info["L1_eigen_val"] * RAD2DEG2
        L2 = mtls_info["L2_eigen_val"] * RAD2DEG2
        L3 = mtls_info["L3_eigen_val"] * RAD2DEG2

        if tls_model=="ISOT":
            t_data = "%5.1f" % (mtls_info["Tr1_eigen_val"] * U2B)

        else:
            Tr1 = mtls_info["Tr1_eigen_val"] * U2B
            Tr2 = mtls_info["Tr2_eigen_val"] * U2B
            Tr3 = mtls_info["Tr3_eigen_val"] * U2B
            t_data = '%5.1f, %5.1f, %5.1f' % (Tr1, Tr2, Tr3),

        ## alternate row background color
        if bgcolor_flag:
            l.append('<tr style="background-color:#dddddd">')
        else:
            l.append('<tr>')
        bgcolor_flag = not bgcolor_flag

        ## path to color thumbnail
        if report_root:
            cpath = os.path.join(report_root, tls["color"]["thumbnail_path"])
        else:
            cpath = tls["color"]["thumbnail_path"]

        l += ['<td align="center" valign="middle"><img src="%s" alt="%s"></td>' % (cpath, tls["color"]["name"]),
              '<td>%s-%s</td>' % (tls["frag_id1"], tls["frag_id2"]),
              '<td>%d</td>'    % (len(tls["segment"])),
              '<td>%d</td>'    % (len(tls_group)),
              '<td>%5.1f</td>' % (tls_info["exp_mean_temp_factor"]),
              '<td>%4.2f</td>' % (tls_info["exp_mean_anisotropy"]),
              '<td>%5.2f</td>' % (tls["rmsd_b"]),
              '<td>%s</td>'    % (t_data),
              '<td>%5.2f, %5.2f, %5.2f</td>' % (L1, L2, L3),
              '<td>%5.1f</td>' % (tls_info["tls_mean_temp_factor"]),
              '<td>%4.2f</td>' % (tls_info["tls_mean_anisotropy"]),
              '</tr>']

    l.append('</table>')
    return "".join(l)


class Report(object):
    """Base class of HTML Report generating objects.
    """
    def html_head(self, title):
        """Header for all HTML pages.
        """
        l = ['<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">',
             '<html>',
             '<head>',
             '<title>%s</title>' % (title),
             '<style type="text/css" media="screen">',
             '<!-- ',
             'BODY {background-color:white; margin-left:5%; margin-right:5%; border-left:5%; border-right:5%; margin-top:2%; border-top:2%;}',
             '-->',
             '</style>',
             '</head>',
             '<body>']

        return "".join(l)

    def html_title(self, title):
        timestr = time.strftime("%d %b %Y", time.localtime(conf.globalconf.start_time))
        
        l  = ['<table border="0" width="100%" style="background-color:#eeeeee"><tr>',
              '<td align="left" valign="top"><font size="-5">%s</font></td>' % (timestr),
              '<td align="right" valign="top"><font size="-5">TLSMD Version %s</font></td>' % (const.VERSION),
              '</tr></table>',
              '<center><font size="+2">%s</font></center><br>' % (title)]
        
        return "".join(l)

    def html_foot(self):
        """Footer for all HTML pages.
        """
        timestr = time.strftime("%d %b %Y", time.localtime(conf.globalconf.start_time))
        
        l  = ['<table border="0" width="100%" style="background-color:#eeeeee"><tr>',
              '<td align="left"><font size="-5">%s</font></td>' % (timestr),
              '<td align="center"><font size="-5">Released %s by <i>%s</i></font></td>' % (const.RELEASE_DATE, const.EMAIL),
              '<td align="right"><font size="-5">TLSMD Version %s</font></td>' % (const.VERSION),
              '</tr></table>',
              '</body></html>']
        
        return "".join(l)


class HTMLReport(Report):
    """Create a through HTML report it its own subdirectory.
    """
    def __init__(self, struct_tls_analysis):
        Report.__init__(self)
        
        self.struct          = struct_tls_analysis.struct
        self.struct_path     = struct_tls_analysis.struct_path
        self.struct_id       = struct_tls_analysis.struct_id
        self.chains          = struct_tls_analysis.chains

        self.page_multi_chain_alignment  = None
        self.pages_chain_motion_analysis = []
        self.page_refinement_prep        = None
        
    def write(self, report_dir):
        """Write out the TLSMD report to the given directory.
        """
        ## create new directory and move into it
        old_dir = os.getcwd()
        if not os.path.isdir(report_dir):
            os.mkdir(report_dir)
        os.chdir(report_dir)

        self.write_cwd()
        
        ## change back to original directory
        os.chdir(old_dir)

    def init_chain_optimization(self, chain):
        """Returns a dictionary with all the calculations on the tls_graph
        needed for the HTML rendering so that no calculations should be
        performed inside the html generating methods.
        """
        chainopt                 = {}
        chainopt["struct"]       = self.struct
        chainopt["struct_path"]  = self.struct_path
        chainopt["chain"]        = chain
        chainopt["struct_id"]    = self.struct_id
        chainopt["chain_id"]     = chain.chain_id
        chainopt["max_ntls"]     = conf.globalconf.nparts
        chainopt["ntls_list"]    = []
        chainopt["tlsopt"]       = {}

        ## calculate the maximum interesting ntls
        minimizer = chain.tls_chain_minimizer

        ## generate the minimized, segmentd TLS groups for 1 TLS
        ## group up to max_ntls and store it in chainopt["ntls_list"]
        num_valid_configurations = 0
        
        for ntls_constraint in range(1, chainopt["max_ntls"]+1):
            tlsopt = minimizer.calc_tls_optimization(ntls_constraint)

            if tlsopt==None:
                continue
            if not tlsopt.is_valid():
                continue

            num_valid_configurations += 1
            
            chainopt["ntls_list"].append((ntls_constraint, tlsopt))
            chainopt["tlsopt"][ntls_constraint] = tlsopt
            
            ## assign a unique color to each tls group in a
            ## chain spanning set of tls groupos
            tlsi = 0
            for tls in tlsopt.tls_list:
                if tls["method"]=="TLS":
                    tls["color"] = self.get_tls_color(tlsi)
                    tlsi += 1
                else:
                    tls["color"] = self.colors[0]

        chainopt["num_valid_configurations"] = num_valid_configurations

        return chainopt

    def init_colors(self):
        """Generated the self.colors dictionary of colors for the report,
        and also writes thumbnail .png images of all the colors.
        """
        thumbnail_dir   = "colors"
        thumbnail_size  = (25, 25)

        ## make the thumbnail subdirectory
        if not os.path.isdir(thumbnail_dir):
            os.mkdir(thumbnail_dir)

        ## clear the colors list
        self.colors = []

        ## skip the first two colors which are black/white
        for i in range(len(Colors.COLORS)):
            cname, rgbf = Colors.COLORS[i]
            
            color = {}
            self.colors.append(color)

            color["index"] = i
            color["name"]  = cname
            color["rgbf"]  = rgbf
            color["rgbi"]  = misc.rgb_f2i(rgbf)

            rgbs = "#%2x%2x%2x" % misc.rgb_f2i(rgbf)
            rgbs = rgbs.replace(" ", "0")
            color["rgbs"]  = rgbs

            ## generate thumbnail image
            color["thumbnail_path"] = os.path.join(thumbnail_dir, "%s.png" % (color["name"]))

            img = Image.new("RGBA", thumbnail_size, color["rgbi"])
            img.save(color["thumbnail_path"], "png")

    def get_tls_color(self, tls_index):
        """Returns the color dict description for a TLS segment of the
        given index, starting at 0.
        """
        ## skip the first two colors; they are black and white
        i = tls_index + 2
        return self.colors[i]

    def write_cwd(self):
        """Write out all the files in the report.
        """
        ## write a local copy of the Structure
        self.struct_path = "%s.pdb" % (self.struct_id)
        FileLoader.SaveStructure(fil=self.struct_path, struct=self.struct)

        ## generate small .png images so  they can be placed in the
        ## TLS group tables to identify the TLS group tabular data
        ## with the generated visualization
        self.init_colors()

        ## all TLSGraph objects get their calculations out of the
        ## way before writing HTML
        chainopt_list = []
        for chain in self.chains:
            misc.begin_chain_timing(chain.chain_id)
            chainopt = self.init_chain_optimization(chain)
	    misc.end_chain_timing(chain.chain_id)
            chainopt_list.append(chainopt)

        ## a report page comparing the tls group segments of all
        ## chains against each other
        self.write_multi_chain_alignment(chainopt_list)

        ## write out all TLSGraph reports
        for chainopt in chainopt_list:
            self.write_tls_chain_optimization(chainopt)

        self.write_refinement_prep(chainopt_list)
 
        ## write out index page
        self.write_index()

    def write_index(self):
        """Writes the main index.html file of the report.
        """
        fil = open("index.html","w")
        fil.write(self.html_index())
        fil.close()

    def html_index(self):
        """Generate and returns the HTML string for the main index.html
        file of the report.
        """
        title = "TLSMD Thermal Parameter Analysis of Structure %s" % (self.struct_id)

        l = [self.html_head(title),
             self.html_title(title),

             ## OPTIMIZATION PARAMETERS
             self.html_globals(),
             '<br>',
        
             ## MOTION ANALYSIS
             '<center>',
             '<h3>TLS Partitions and Motion Analysis of Individual Chains</h3>',
             '</center>',
             '<p style="font-size:small">%s</p>' % (MOTION_ANALYSIS_TEXT)]

        for xdict in self.pages_chain_motion_analysis:
            l.append('<p><a href="%s">%s</a></p>' % (xdict["href"], xdict["title"]))

        l +=['<br>',
             ## MULTI CHAIN ALIGNMENT
             '<center><h3>Multi-Chain TLS Group Alignment</h3></center>',
             '<p style="font-size:small">%s</p>' % (MULTI_CHAIN_ALIGNMENT_TEXT)]

        if self.page_multi_chain_alignment!=None:
            l.append('<p><a href="%s">%s</a></p>' % (self.page_multi_chain_alignment["href"], self.page_multi_chain_alignment["title"]))
        else:
            l.append('<p><u>Only one chain was analyized in this structure, so the multi-chain alignment analysis was not performed.</u></p>')

        l +=['<br>',
             ## REFINEMENT PREP
             '<center><h3>Use Optimal TLS Groups with Refmac5 TLS Refinement</h3></center>',
             '<p style="font-size:small">%s</p>' % (REFINEMENT_PREP_TEXT),
             '<p><a href="%s">%s</a></p>' % (self.page_refinement_prep["href"], self.page_refinement_prep["title"]),
             self.html_foot()]
        
        return "".join(l)

    def html_globals(self):
        """
        """
        if conf.globalconf.tls_model in ["ISOT", "NLISOT"]:
            tls_model = "Isotropic"
        elif conf.globalconf.tls_model in ["ANISO", "NLANISO"]:
            tls_model = "Anisotropic"

        if conf.globalconf.weight_model == "UNIT":
            weight = 'Unit Weights (All Weights 1.0)'
        elif conf.globalconf.weight_model == "IUISO":
            weight = 'Input Structure Atoms Weighted by <var>1.0/B<sub>iso</sub></var>'

        if conf.globalconf.include_atoms == "MAINCHAIN":
            include_atoms = "Main Chain Protein Atoms (N, CA, C, O, CB)"
        else:
            include_atoms = "All Protein Atoms"

        l = ['<center>',
             '<table border="0" cellpadding="3" width="75%" style="background-color:#eeeeee; font-size:small">',
             '<tr style="background-color:#cccccc"><th>Program Option</th><th>Setting</th></tr>',
             '<tr style="background-color:#dddddd"><td>Temperature Factors</td><td><b>%s</b></td></tr>' % (tls_model),
             '<tr><td>Minimum TLS Segment Length</td><td><b>%s Residues</b></td></tr>' % (conf.globalconf.min_subsegment_size),
             '<tr style="background-color:#dddddd"><td>Atoms Analyzed</td><td><b>%s</b></td></tr>' % (conf.globalconf.include_atoms),
             '</table>',
             '</center>']

        return "".join(l)
    
    def write_tls_chain_optimization(self, chainopt):
        """Writes the HTML report analysis of a single TLS graphed chain.
        """
        misc.begin_chain_timing(chainopt["chain_id"])
            
        path  = "%s_CHAIN%s_ANALYSIS.html" % (self.struct_id, chainopt["chain_id"])
        title = "Chain %s TLS Analysis" % (chainopt["chain_id"])

        self.pages_chain_motion_analysis.append(
            {"title": title,
             "href":  path })

        fil = open(path, "w")
        fil.write(self.html_tls_chain_optimization(chainopt))
        fil.close()

        misc.end_chain_timing(chainopt["chain_id"])
        
    def html_tls_chain_optimization(self, chainopt):
        """Generates and returns the HTML string report analysis of a
        single TLS graphed chain.
        """
        title = "Chain %s TLS Analysis of %s" % (chainopt["chain_id"], self.struct_id)
        
        x  = self.html_head(title)
        x += self.html_title(title)
        x += '<center><a href="index.html">Back to Index</a></center>'
        x += '<br>'

        ## if there were no valid chain configurations found
        ## then write out a useful error message
        if chainopt["num_valid_configurations"]==0:
            x += '<p>%s</p>' % (NO_VALID_CONFIGURATIONS)
            x += self.html_foot()
            return x

        ## TLS Segments vs. Residual
        x += self.html_chain_lsq_residual_plot(chainopt)
        
        ## generate a plot comparing all segmentations
        x += self.html_chain_alignment_plot(chainopt)

        ## add tables for all TLS group selections using 1 TLS group
        ## up to max_ntls
        for ntls_constraint in range(1, chainopt["max_ntls"]+1):
            tmp = self.html_tls_graph_path(chainopt, ntls_constraint)
            if tmp!=None:
                x += tmp

            ## maybe this will help with the memory problems...
            import gc
            gc.collect()

        x += self.html_foot()
        return x

    def html_chain_lsq_residual_plot(self, chainopt):
        """Generates the Gnuplot/PNG image plot, and returns the HTML
        fragment for its display in a web page.
        """
        gp = gnuplots.LSQR_vs_TLS_Segments_Plot(chainopt)

        title = 'Chain %s Optimization Residual' % (chainopt["chain_id"])

        l = ['<center>',
             gp.html_markup(title, LSQR_CAPTION),
             '</center>']
        
        return "".join(l)

    def html_chain_alignment_plot(self, chainopt):
        """generate a plot comparing all segmentations
        """
        plot = sequence_plot.TLSSegmentAlignmentPlot()
        
        for ntls, tlsopt in chainopt["ntls_list"]:
            plot.add_tls_segmentation(chainopt, ntls)

        ## create filename for plot PNG image file
        plot_path = "%s_CHAIN%s_ALIGN.png" % (self.struct_id, chainopt["chain_id"])
        
        plot.plot(plot_path)

        l = ['<center><h3>TLS Partition Segment Alignment of Chain %s</h3></center>' % (chainopt["chain_id"]),
             '<center>',
             '<table border="0" style="background-color:#dddddd">',
             '<tr><th># of TLS<br>Groups</th>',
             '<th>Segment/Sequence Alignment</th></tr>',
             '<tr>',
        
             '<td align="right" valign="top">',
             '<p style="font-size:xx-small; margin-top:%dpx; line-height:18px">' % (plot.border_width)]

        for ntls, tlsopt in chainopt["ntls_list"]:
            l.append('<a href="#NTLS%d">%d</a><br>' % (ntls, ntls))

        l +=['</td>',
             '<td><img src="%s" alt="Sequence Alignment Plot"></td>' % (plot_path),
             '</tr>',
             '</table>',
             '</center>',
             '<p>%s</p>' % (SEG_ALIGN_CAPTION)]
        
        return "".join(l)

    def html_tls_graph_path(self, chainopt, ntls):
        """Generates the HTML table describing the path (set of tls groups)
        for the given number of segments(h, or ntls)
        """
        ## select the correct TLSChainDescription() for the number of ntls
        if not chainopt["tlsopt"].has_key(ntls):
            return None

        tlsopt = chainopt["tlsopt"][ntls]

        ## write out PDB file
        self.write_tls_pdb_file(chainopt, tlsopt, ntls)

        ## Raster3D Image
        pml_path, png_path = self.raster3d_render_tls_graph_path(chainopt, tlsopt, ntls)

        ## JMol Viewer Page
        jmol_path = self.jmol_html(chainopt, tlsopt, ntls)
        jmol_animate_path = self.jmol_animate_html(chainopt, tlsopt)

        ## tlsout file
        tlsout_path = self.write_tlsout_file(chainopt, tlsopt, ntls)

        ## detailed analysis of all TLS groups
        ntls_analysis = self.chain_ntls_analysis(chainopt, tlsopt)

        ## BMean Plot
        ntls_analysis.bmean_plot.width = 640
        ntls_analysis.bmean_plot.height = 250
        ntls_analysis.bmean_plot.tls_group_titles = False
        ntls_analysis.bmean_plot.output_png()

        
        l = ['<hr>',
             
             '<center style="page-break-before:always; font-size:small">',
             '<a name="NTLS%d" style="font-size:large">' % (ntls),
             'Summary of Optimal TLS Group Partition using %d Groups' % (ntls),
             '</a>',
             '<br>',

             ## navigation links
             '<a href="%s">Full TLS Partition Analysis</a>' % (ntls_analysis.url),

             '&nbsp;&nbsp;&nbsp;&nbsp;',
         
             '<a href="." onClick="window.open(&quot;%s&quot;,&quot;&quot;,&quot;' % (jmol_path),
             'width=%d,height=%d,screenX=10,screenY=10,left=10,top=10&quot;);' % (JMOL_SIZE, JMOL_SIZE),
             'return false;">View with JMol</a>',

             '&nbsp;&nbsp;&nbsp;&nbsp;',

             '<a href="." onClick="'\
             'window.open('\
             '&quot;%s&quot;,'\
             '&quot;&quot;,'\
             '&quot;width=%d,height=%d,screenX=10,'\
             'screenY=10,left=10,top=10&quot;);'\
             'return false;">Animate Screw Displacement with JMol</a>' % (jmol_animate_path, JMOL_SIZE, JMOL_SIZE),

             '<br>',
             '<a href="%s">Download TLSOUT File for TLSView</a>' % (tlsout_path),
             '&nbsp;&nbsp;&nbsp;&nbsp;',
             '<a href="%s">Generate PDBIN/TLSIN Files for REFMAC5</a>' % ("%s_REFINEMENT_PREP.html" % (self.struct_id)),
             
             '</center>',

             ## raytraced image
             '<table style="background-color:white" width="100%" border=0>',
             '<tr><th>',
             '<center><img src="%s" alt="structimage"></center><br>' % (png_path),
             '</th></tr>',

             '<tr><th><center>',
             ntls_analysis.bmean_plot.html_link(),
             '</center></th></tr>',
             '</table>',

             ## now the table
             html_tls_group_table(chainopt, tlsopt, ntls),
             
             '<br clear="all">']
        
        return "".join(l)

    def raster3d_render_tls_graph_path(self, chainopt, tlsopt, ntls):
        """Render TLS visualizations using Raster3D.
        """
        basename = "%s_CHAIN%s_NTLS%d" % (self.struct_id, chainopt["chain_id"], ntls)
        png_path = "%s.png"   % (basename)

        misc.start_timing()
        print "Raster3D: rendering %s..." % (basename)

        struct_id = self.struct_id
        chain     = chainopt["chain"]
        chain_id  = chainopt["chain_id"]

        driver = R3DDriver.Raster3DDriver()

        ## XXX: Size hack: some structures have too many chains,
        ## or are just too large
        show_chain = {}
        for chx in self.struct.iter_chains():
            if chx.chain_id==chain_id:
                show_chain[chx.chain_id] = True
                continue
            
            if chx.count_fragments()>=20:
                show_chain[chx.chain_id] = False
                continue
            
            show_chain[chx.chain_id] = True
        ## end size hack

        viewer = Viewer.GLViewer()
        gl_struct = viewer.glv_add_struct(self.struct)

        ## orient the structure with the super-spiffy orientation algorithm
        ## which hilights the chain we are examining
        ori = calc_orientation(self.struct, chain)
        viewer.glo_update_properties(
            R         = ori["R"],
            cor       = ori["centroid"],
            zoom      = ori["hzoom"],
            near      = ori["near"],
            far       = ori["far"],
            width     = ori["pwidth"],
            height    = ori["pheight"],
            bg_color  = "White")

        ## turn off axes and unit cell visualization
        gl_struct.glo_update_properties_path("gl_axes/visible", False)
        gl_struct.glo_update_properties_path("gl_unit_cell/visible", False)

        ## setup base structural visualization
        for gl_chain in gl_struct.glo_iter_children():
            if not isinstance(gl_chain, Viewer.GLChain):
                continue

            ## chain is hidden
            if show_chain.get(gl_chain.chain.chain_id, False)==False:
                gl_chain.properties.update(visible=False)
                continue

            if gl_chain.chain.chain_id==chain_id:
                gl_chain.properties.update(
                    oatm_visible       = False,
                    side_chain_visible = False,
                    hetatm_visible     = True,
                    color              = "0.20,0.20,0.20",
                    lines              = False,
                    ball_stick         = False,
                    trace              = True,
                    trace_radius       = 0.20,
                    trace_color         = "0.20,0.20,0.20", )
            else:
                gl_chain.properties.update(
                    visible       = True,
                    ball_stick    = False,
                    cpk           = True)
                
        ## add the TLS group visualizations
        for tls in tlsopt.tls_list:
            if tls["method"]!="TLS":
                continue
            
            tls_name = "TLS_%s_%s" % (tls["frag_id1"], tls["frag_id2"])
            
            gl_tls_group = TLS.GLTLSGroup(
                oatm_visible       = False,
                side_chain_visible = False,
                hetatm_visible     = True,
                adp_prob           = conf.ADP_PROB,
                L_axis_scale       = 2.0,
                L_axis_radius      = 0.25,
		both_phases        = True,
                tls_group          = tls["tls_group"],
                tls_info           = tls["model_tls_info"],
                tls_name           = tls_name,
                tls_color          = tls["color"]["name"])
            gl_struct.glo_add_child(gl_tls_group)

            ## set width of trace according to the group's translational tensor trace
            mtls_info = tls["model_tls_info"]
            tiso = (mtls_info["Tr1_eigen_val"] + mtls_info["Tr2_eigen_val"] + mtls_info["Tr3_eigen_val"]) / 3.0

            ## too big usually for good visualization -- cheat and scale it down
            radius = 0.5 * Gaussian.GAUSS3C[conf.ADP_PROB] * TLS.calc_rmsd(tiso)
            radius = max(radius, 0.30)
            radius = 0.50
            
            gl_tls_group.gl_atom_list.properties.update(trace_radius = radius)
            gl_tls_group.glo_update_properties(time = 0.25)
            
        driver.glr_set_render_png_path(png_path)
        viewer.glv_render_one(driver)
        print misc.end_timing()

        return "", png_path

    def write_tls_pdb_file(self, chainopt, tlsopt, ntls):
        """Write out a PDB file with the TLS predicted anisotropic ADPs for
        this segmentation.
        """
        basename = "%s_CHAIN%s_NTLS%d_UTLS"  % (self.struct_id, chainopt["chain_id"], ntls)
        pdb_path = "%s.pdb" % (basename)

        ## temporarily set the atom temp_factor and U tensor to the Utls value
        old_temp_factor = {}
        old_U = {}
        for tls in tlsopt.tls_list:
            tls_group = tls["tls_group"]
            
            for atm, Utls in tls_group.iter_atm_Utls():
                old_temp_factor[atm] = atm.temp_factor
                old_U[atm] = atm.U

                atm.temp_factor = U2B * (numpy.trace(Utls)/3.0)
                atm.U = Utls

        FileLoader.SaveStructure(fil=pdb_path, struct=self.struct)

        ## restore atom temp_factor and U
        for atm, temp_factor in old_temp_factor.iteritems():
            atm.temp_factor = temp_factor
            atm.U = old_U[atm]

    def write_tlsout_file(self, chainopt, tlsopt, ntls):
        """Writes the TLSOUT file for the segmentation.
        """
        basename = "%s_CHAIN%s_NTLS%d" % (self.struct_id, chainopt["chain_id"], ntls)
        tlsout_path = "%s.tlsout" % (basename)

        struct_id = self.struct_id
        chain_id  = chainopt["chain_id"]

        tls_file = TLS.TLSFile()
        tls_file.set_file_format(TLS.TLSFileFormatTLSOUT())

        for tls in tlsopt.tls_list:
            ## don't write out bypass edges
            if tls["method"]!="TLS":
                continue
            
            tls_desc = TLS.TLSGroupDesc()
            tls_file.tls_desc_list.append(tls_desc)
            
            tls_desc.set_tls_group(tls["tls_group"])
            tls_desc.add_range(chain_id, tls["frag_id1"], chain_id, tls["frag_id2"], "ALL")

        tls_file.save(open(tlsout_path, "w"))

        return tlsout_path

    def chain_ntls_analysis(self, chainopt, tlsopt):
        """Generate ntls optimization constraint report and free memory.
        """
        report = ChainNTLSAnalysisReport(chainopt, tlsopt, tlsopt.ntls)
        return report

    def jmol_html(self, chainopt, tlsopt, ntls):
        """Writes out the HTML page which will display the
        structure using the JMol Applet.
        """
        jmol_path = "%s_CHAIN%s_NTLS%d_JMOL.html"  % (self.struct_id, chainopt["chain_id"], ntls)

        ## create the JMol script using cartoons and consistant
        ## coloring to represent the TLS groups
        js = ['load %s;' % (self.struct_path),
              'select *;',
              'cpk off;',
              'wireframe off;',
              'select protein;',
              'cartoon on;']

        ## loop over TLS groups and color
        for tls in tlsopt.tls_list:
            js.append('select %s-%s:%s;' % (tls["frag_id1"], tls["frag_id2"], tls["chain_id"]))
            js.append('color [%d,%d,%d];' % (tls["color"]["rgbi"]))

        ## select non-protein non-solvent and display
        js +=['select not protein and not solvent;',
              'color CPK;',
              'wireframe on; wireframe 0.5;',
              'spacefill 80%;',
              'spacefill on;']
        
        ## write the HTML page to render the script in
        l = ['<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">',
             '<html>',
             '<head>',
             '<title>Chain %s using %d TLS Groups</title>' % (chainopt["chain_id"], ntls),
             '<script type="text/javascript" src="%s/Jmol.js">' % (JMOL_DIR),
             '</script>',
             '</head>',
             '<body>',
             '<script type="text/javascript">',
             'jmolInitialize("%s");' % (JMOL_DIR),
             'jmolSetAppletColor("white");',
             'jmolApplet(%d, "%s");' % (JMOL_SIZE, "".join(js)),
             '</script>',
             '</body>',
             '</html>']

        open(jmol_path, "w").write("".join(l))
        return jmol_path

    def jmol_animate_html(self, chainopt, tlsopt):
        """Writes out the HTML page which will display the
        structure using the JMol Applet.
        """
        basename = "%s_CHAIN%s_NTLS%d_ANIMATE" % (self.struct_id, chainopt["chain_id"], tlsopt.ntls)

        html_path = "%s.html" % (basename)
        pdb_path  = "%s.pdb" % (basename)

        ## generate animation PDB file

        try:
            print "TLSAnimate: creating animation PDB file..."
            misc.start_timing()
            tlsa = TLSAnimate(self.struct, chainopt, tlsopt)
            tlsa.construct_animation(pdb_path)
            print misc.end_timing()
        except TLSAnimateFailure:
            pass
        
        ## create the JMol script using cartoons and consistant
        ## coloring to represent the TLS groups
        js = ['load %s;' % (pdb_path),
              'select *;',
              'cpk off;',
              'wireframe off;',
              'select protein;',
              'trace on;']

        ## loop over TLS groups and color
        for tls in tlsopt.tls_list:
            chain_ids = [tlsa.L1_chain.chain_id,
                         tlsa.L2_chain.chain_id,
                         tlsa.L3_chain.chain_id]

            for chain_id in chain_ids:
                js.append('select %s-%s:%s;' % (tls["frag_id1"], tls["frag_id2"], chain_id))
                js.append('color [%d,%d,%d];' % (tls["color"]["rgbi"]))

        ## select non-protein non-solvent and display
        js +=['select not protein and not solvent;',
              'color CPK;',
              'wireframe on;',
              'wireframe 0.5;',
              'spacefill 80%;',
              'spacefill on;',
              'anim fps 2;',
              'anim mode loop 0 0;',
              'anim on;']

        ## write the HTML page to render the script in
        l = ['<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">',
             '<html>',
             '<head>',
             '<title>Chain %s using %d TLS Groups</title>' % (chainopt["chain_id"], tlsopt.ntls),
             '<script type="text/javascript" src="%s/Jmol.js">' % (JMOL_DIR),
             '</script>',
             '</head>',
             '<body>',
             '<script type="text/javascript">',
             'jmolInitialize("%s");' % (JMOL_DIR),
             'jmolSetAppletColor("white");',
             'jmolApplet(%d, "%s");' % (JMOL_SIZE, "".join(js)),
             '</script>',
             '</body>',
             '</html>']

        ## manually free memory
        tlsa = None
        import gc
        gc.collect()

        open(html_path, "w").write("".join(l))
        return html_path

    def write_multi_chain_alignment(self, chainopt_list):
        """Write out the chain residue alignment page.
        """
        ## only write out the comparison page if there is more than one
        ## chain analyzed in the structure
        if len(chainopt_list)<2:
            return
        
        path  = "%s_CHAIN_COMP.html" % (self.struct_id)
        title = "Multi-Chain Alignment Analysis of TLS Groups"

        self.page_multi_chain_alignment = {
            "title": title,
            "href":  path }

        fil = open(path, "w")
        fil.write(self.html_multi_chain_alignment(chainopt_list))
        fil.close()
    
    def html_multi_chain_alignment(self, chainopt_list):
        """Write out all HTML/PDB/TLSIN files which compare
        chains in the structure.
        """
        title = self.page_multi_chain_alignment["title"]

        x  = self.html_head(title)
        x += self.html_title(title)

        x += '<center>'
        x += '<a href="index.html">Back to Index</a>'
        x += '</center>'
        x += '<br>'

        ## figure out the maximum number of ntls in all chains
        max_ntls = 0
        for chainopt in chainopt_list:
            max_ntls = max(max_ntls, chainopt["max_ntls"])

        ## generate ntls number of plots and add them to the
        ## HTML document
        for ntls in range(1, max_ntls+1):

            ## create a 2-tuple list of (chain_id, chainopt) for
            ## each chain which a a TLSMD segmentation of h groups
            seg_list = []
            for chainopt in chainopt_list:
                if chainopt["tlsopt"].has_key(ntls):
                    seg_list.append((chainopt["chain_id"], 
                                     chainopt["tlsopt"][ntls]))

            ## generate PDB and TLSIN files containing the TLS
            ## predicted anisotropic ADPs for all chains for the
            ## given number of tls segments
            basename    = "%s_NTLS%d"  % (self.struct_id, ntls)
            tlsout_path = "%s.tlsout" % (basename)
            pdb_path    = "%s.pdb" % (basename)

            old_temp_factor = {}
            old_U = {}

            for chain_id, tlsopt in seg_list:
                for tls in tlsopt.tls_list:
                    tls_group = tls["tls_group"]

                    for atm, Utls in tls_group.iter_atm_Utls():
                        old_temp_factor[atm] = atm.temp_factor
                        old_U[atm] = atm.U
                        atm.temp_factor = U2B * (Utls[0,0] + Utls[1,1] + Utls[2,2]) / 3.0
                        atm.U = Utls

            ## save the structure file
            FileLoader.SaveStructure(fil=pdb_path, struct=self.struct)

            ## restore atom temp_factor and U
            for atm, temp_factor in old_temp_factor.iteritems():
                atm.temp_factor = temp_factor
                atm.U = old_U[atm]

            ## generate the TLS segmentation alignment plot for all chains
            plot = segment_plot.TLSSegmentAlignmentPlot()
            chain_id_list = []

            for chainopt in chainopt_list:
                if not chainopt["tlsopt"].has_key(ntls):
                    continue
                
                chain_id_list.append(chainopt["chain_id"])
                plot.add_tls_segmentation(chainopt, ntls)

            plot_path = "%s_ALIGN.png" % (basename)
            plot.plot(plot_path)

            ## write HTML
            x += '<h3>Chains Alignment using %d TLS Groups</h3>' % (ntls)

            ## plot table
            x += '<table border="1">'
            x += '<tr><th>Chain</th><th>Chain Alignment</th></tr>'
            x += '<tr>'
        
            x += '<td align="center">'
            x += '<table border="0" cellspacing="0" cellpadding="0">'

            for chain_id in chain_id_list:
                x += '<tr><td align="right" valign="middle" height="20"><font size="-20">%s</font></td></tr>' % (chain_id)

            x += '</table>'
            x += '</td>'

            x += '<td><img src="%s" alt="Segmentation Plot"></td>' % (plot_path)

            x += '</tr>'
            x += '</table>'

        x += self.html_foot()
        return x

    def write_refinement_prep(self, chainopt_list):
        """Generate form to allow users to select the number of TLS groups
        to use per chain.
        """
        path  = "%s_REFINEMENT_PREP.html" % (self.struct_id)
        title = "Use Optimal TLS Groups with Refmac5 TLS Refinement"
 
        self.page_refinement_prep = {
            "title": title,
            "href" : path }

        fil = open(path, "w")
        fil.write(self.html_refinement_prep(chainopt_list))
        fil.close()

    def html_refinement_prep(self, chainopt_list):
        title = self.page_refinement_prep["title"]

        x  = self.html_head(title)
        x += self.html_title(title)

        x += '<center><h3>'
        x += 'Step 1: Select the number of TLS groups for each chain'
        x += '</h3></center>'

        x += '<center>'
        x += '<a href="index.html">Back to Index</a>'
        x += '</center>'
        x += '<br>'

        x += '<form enctype="multipart/form-data" action="%s" method="post">' % (conf.REFINEPREP_URL)
        x += '<input type="hidden" name="job_id" value="%s">' % (conf.globalconf.job_id)

        x += '<p>%s</p>' % (REFINEMENT_PREP_INFO)
        
        x += '<center><table><tr><td>'
        
        plot = gnuplots.LSQR_vs_TLS_Segments_All_Chains_Plot(chainopt_list)
        x += plot.html_link()

        x += '</td></tr><tr><td>'

        x += '<table width="100%" border="1">'
        x += '<tr><th>'
        x += '<p>Select the Number of TLS Groups per Chain</p>'
        x += '</th></tr>'

        x += '<tr><td align="center">'

        x += '<table cellspacing="5">'
        for chainopt in chainopt_list:
            chain_id = chainopt["chain_id"]
            
            x += '<tr><td>'
            x += 'Number of TLS Groups for Chain %s' % (chain_id)
            x += '</td><td>'
        
            x += '<select name="NTLS_CHAIN%s">' % (chain_id)
            for ntls, tlsopt in chainopt["ntls_list"]:
                x += '<option value="%d">%d</option>' % (ntls, ntls)
            x += '</select>'

            x += '</td></tr>'
        x += '</table>'
        
        x += '</td></tr>'

        x += '<tr><td align="right">'
        x += '<input type="submit" value="OK">'
        x += '</td></tr></table>'
        
        x += '</td></tr></table></center>'

        x += '</form>'
        
        x += self.html_foot()
        return x


class ChainNTLSAnalysisReport(Report):
    """Writes a HTML report detailing one given TLS segmentation of a chain.
    """
    def __init__(self, chainopt, tlsopt, ntls):
        Report.__init__(self)

        self.struct      = chainopt["struct"]
        self.struct_id   = chainopt["struct_id"]
        self.struct_path = chainopt["struct_path"]
        self.chain_id    = chainopt["chain_id"]

        
        self.chainopt = chainopt
        self.tlsopt   = tlsopt
        self.ntls     = ntls

        self.root  = ".."
        self.dir   = "%s_CHAIN%s_NTLS%d"  % (self.struct_id, self.chain_id, self.ntls)
        self.index = "%s.html" % (self.dir)
        self.url   = "%s/%s" % (self.dir, self.index)

        self.write_report()

    def write_report(self):
        if not os.path.isdir(self.dir):
            os.mkdir(self.dir)
        os.chdir(self.dir)

        self.write_all_files()

        os.chdir(self.root)

    def write_all_files(self):
        """Writes analysis details of each TLS group.
        """
        title = "Chain %s Partitioned by %d TLS Groups" % (self.chain_id, self.ntls)
        path = "%s_CHAIN%s_ANALYSIS.html" % (self.struct_id, self.chain_id)


        l = [self.html_head(title),
             self.html_title(title),

             '<center>',
             '<a href="../index.html">Back to Index</a>',
             '&nbsp;&nbsp;&nbsp;&nbsp;',
             '<a href="../%s">Back to Chain %s Analysis</a>' % (path, self.chain_id),
             '</center>',
             
             '<br>',

             self.html_tls_group_table(),'<br>',
             self.html_bmean(),'<br>',
             self.html_translation_analysis(),'<br>',
             self.html_libration_analysis(),'<br>',
             self.html_ca_differance(),'<br>',
             self.html_rmsd_plot()]
        
        for tls in self.tlsopt.tls_list:
            ## don't write out bypass edges
            if tls["method"]!="TLS": continue

            l.append(self.html_tls_fit_histogram(tls))

        l.append(self.html_foot())
        
        open(self.index, "w").write("".join(l))

    def html_tls_group_table(self):
        return html_tls_group_table(self.chainopt, self.tlsopt, self.ntls, "..")

    def html_translation_analysis(self):
        """Perform a translation analysis of the protein chain as
        spanned by the tlsopt TLS groups.
        """
        tanalysis = gnuplots.TranslationAnalysis(self.chainopt, self.tlsopt)

        l = ['<center>',
             tanalysis.html_markup("Translation Analysis of T<sup>r</sup>", TRANSLATION_GRAPH_CAPTION),
             '</center>']
        
        return "".join(l)

    def html_libration_analysis(self):
        """Perform a libration analysis of the protein chain as
        spanned by the tlsopt TLS groups.
        """
        libration_analysis = gnuplots.LibrationAnalysis(self.chainopt, self.tlsopt)
        
        l = ['<center>',
             libration_analysis.html_markup("Screw Displacement Analysis", LIBRATION_GRAPH_CAPTION),
             '</center>']
        
        return "".join(l)

    def html_ca_differance(self):
        """Perform a fit analysis of the protein chain as
        spanned by the tlsopt TLS groups.
        """
        plot = gnuplots.CA_TLS_Differance_Plot(self.chainopt, self.tlsopt)
        
        l = ['<center>',
             plot.html_markup("Deviation of Observed CA Atom B-Factors From TLS Model", FIT_GRAPH_CAPTION),
             '</center>']

        return "".join(l)

    def html_bmean(self):
        """Mean B-Factor per residue.
        """
        self.bmean_plot = gnuplots.BMeanPlot(self.chainopt, self.tlsopt)

        l = ['<center>',
             self.bmean_plot.html_markup("Mean BFactor Analysis", "Comparison of TLS predicted B factors with experimental (input) B factors."),
             '</center>']

        return "".join(l)

    def html_rmsd_plot(self):
        rmsd_plot = gnuplots.RMSDPlot(self.chainopt, self.tlsopt)

        l = ['<center>',
             rmsd_plot.html_markup("RMSD Deviation of Observed vs. TLS Predicted B Factors", ""),
             '</center>']

        return "".join(l)

    def html_tls_fit_histogram(self, tls):
        """histogram of atomic U_ISO - U_TLS_ISO
        """
        
        his = UIso_vs_UtlsIso_Histogram(self.chainopt, self.tlsopt, tls)

        title = 'Distribution Histogram of TLS Group %s%s-%s%s' % (self.chain_id, tls["frag_id1"], self.chain_id, tls["frag_id2"])
        
        l = ['<center>',
             his.html_markup(title, ""),
             '</center>']

        return "".join(l)



