## TLS Motion Determination (TLSMD)
## Copyright 2002-2009 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

###############################################################################
## Report Directory Generation
##

## Python
import os
import time
import numpy
import shutil    ## for copying files around (e.g., jmol)
import sys       ## for try/except
import itertools ## used in selecting backbone atoms, 2009-01-13
import re        ## used in summary_file_update(), 2009-05-07

## Python Imaging Library imports
import Image
import ImageDraw

## Pymmlib
from mmLib import Constants, Colors, Viewer, R3DDriver, Structure, Gaussian, FileIO, TLS

## TLSMD
import misc
import const
import conf
import console
import gnuplots
import sequence_plot
import table
import captions
from tls_animate import TLSAnimate, TLSAnimateFailure


class MyError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def calc_inertia_tensor(atom_iter):
    """Calculate moment of inertia tensor at the centroid
    of the atoms.
    """
    al       = Structure.AtomList(atom_iter)
    centroid = al.calc_centroid()

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

    evals, evecs = numpy.linalg.eig(I)

    ## FIXME: Fix "Warning: invalid value encountered in divide", 2009-05-25
    elist = [(evals[0], evecs[0]),
             (evals[1], evecs[1]),
             (evals[2], evecs[2])]

    elist.sort()

    R = numpy.array((elist[0][1], elist[1][1], elist[2][1]), float)

    ## make sure the tensor uses a right-handed coordinate system
    if numpy.allclose(numpy.linalg.det(R), -1.0):
        I = numpy.identity(3, float)
        I[0,0] = -1.0
        R = numpy.dot(I, R)
    assert numpy.allclose(numpy.linalg.det(R), 1.0)

    return centroid, R


def calc_orientation(struct, chain):
    """Orient the structure based on a moment-of-inertia like tensor
    centered at the centroid of the structure.
    """
    ori = {}

    def iter_atoms(sobjx):
        for fragx in sobjx.iter_standard_residues():
            for atmx in fragx.iter_atoms():
                yield atmx

    centroids, Rs = calc_inertia_tensor(iter_atoms(struct))
    centroidc, Rc = calc_inertia_tensor(iter_atoms(chain))

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

    for atm in iter_atoms(chain):
        x  = numpy.dot(R, atm.position - centroid)

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

    pheight = conf.VIS_WIDTH
    pwidth  = pheight * (width/height)

    ori["R"]        = R
    ori["centroid"] = centroid + numpy.dot(numpy.transpose(R), xydelta)
    ori["pwidth"]   = pwidth
    ori["pheight"]  = pheight 
    ori["hzoom"]    = zoom

    ## calculate near, far clipping plane
    ori["near"] = max_z
    ori["far"]  = min_z

    return ori


class ColorInfo(object):
    def __init__(self, index, name, rgbf, thumbnail_dir, thumbnail_size = (25, 25)):
        self.index = index
        self.name = name
        self.rgbf = rgbf
        self.rgbi = misc.rgb_f2i(rgbf)
        self.rgbs = misc.rgb_f2s(rgbf)
        self.thumbnail_dir = thumbnail_dir
        self.thumbnail_size = thumbnail_size

        self.thumbnail_path = os.path.join(thumbnail_dir, "%s.png" % (name))
        img = Image.new("RGBA", thumbnail_size, self.rgbi)
        img.save(self.thumbnail_path, "png")

def html_tls_group_table(ntls, chain, cpartition, report_root = None, detail = None):
    """Generate HTML for a table containing the details of the ntls-group
    partitioning of the given chain.
    """
    ## inspect the first tls group dictionary to determine TLS model type
    try:
        tls = cpartition.tls_list[0]
    except IndexError:
        return ""

    tls_model = tls.tls_group.model

    if tls_model == "ISOT":
        t_head = 'T<sup>r</sup> <var>B</var>'

    elif tls_model == "ANISO":
        t_head = 'eval(T<sup>r</sup>) <var>B</var>'

    else:
        return ""

    ## XXX
    console.stdoutln("RMSD_B-%s: %.2f" % (ntls, cpartition.rmsd_b()))
    console.stdoutln("RESIDUAL-%s: %.2f" % (ntls, cpartition.residual()))

    l = ['<table class="tls_segments">',
         '<tr>',
         '<th align="center" colspan="12">Analysis of TLS Group %s Chain Segments ' % ntls,
         '(overall rmsd_b=%.2f and residual=%.2f)</th>' % (
             cpartition.rmsd_b(), cpartition.residual()),
         '</tr>',

         '<tr>',
         '<th colspan="7" style="background-color:#aaaaaa">Input Structure</th>',
         '<th colspan="5" style="background-color:#bbbbbb">TLS Predictions</th>',
         '</tr>',

         '<tr style="background-color:#bbbbbb">',
         '<th>Color</th>',
         '<th>Segment</th>',
         '<th>Residues</th>',
         '<th>Atoms</th>',
         '<th>&#60;B&#62;</th>',
         '<th>B<sub>rmsd</sub></th>',
         '<th>&#60;Aniso&#62;</th>',
         '<th>RMSD B</th>',
         '<th>%s</th>' % (t_head),
         '<th>eval(L) <var>DEG<sup>2</sup></var></th>',
         '<th>&#60;B&#62;</th>',
         '<th>&#60;Aniso&#62;</th>',
         '</tr>' ]

    bgcolor_flag = True

    i = 1 ## which segment within the partition
    for tls in cpartition.iter_tls_segments():

        ## Calculate the stddev for all temperature factors in a given segment
        tmp_temp_factor = []
        for atm, Utls in tls.tls_group.iter_atm_Utls():
            tmp_temp_factor.append(atm.temp_factor)
        stddev = numpy.std(tmp_temp_factor)

        tls_group = tls.tls_group
        mtls_info = tls.model_tls_info
        ## EAM DEBUG - I think this results from a previous exception in 
        ## html_tls_graph_path()
        if mtls_info == None:
            l += ['<tr style="background-color:#ffeeee">',
                  '<td colspan="12" align-text="center">Error</td></tr>']
            continue

        L1 = mtls_info["L1_eigen_val"] * Constants.RAD2DEG2
        L2 = mtls_info["L2_eigen_val"] * Constants.RAD2DEG2
        L3 = mtls_info["L3_eigen_val"] * Constants.RAD2DEG2

        if tls_model=="ISOT":
            t_data = "%5.1f" % (mtls_info["Tr1_eigen_val"] * Constants.U2B)

        else:
            Tr1 = mtls_info["Tr1_eigen_val"] * Constants.U2B
            Tr2 = mtls_info["Tr2_eigen_val"] * Constants.U2B
            Tr3 = mtls_info["Tr3_eigen_val"] * Constants.U2B
            t_data = '%5.1f, %5.1f, %5.1f' % (Tr1, Tr2, Tr3),

        ## alternate row background color
        if bgcolor_flag:
            l.append('<tr style="background-color:#dddddd">')
        else:
            l.append('<tr>')
        bgcolor_flag = not bgcolor_flag

        ## path to color thumbnail
        if report_root:
            cpath = os.path.join(report_root, tls.color.thumbnail_path)
        else:
            cpath = tls.color.thumbnail_path

        ##======================================================================
        ##<FLATFILE>
        ## NOTE: This is in a standalone def "html_tls_group_table"
        if detail:
            flatfile_name = "%s/%s.dat" % (os.getcwd(), conf.globalconf.job_id)
            flatfile = open(flatfile_name, "a+")

            chain_ntls = "%s,%s.%s" % (chain.chain_id, ntls, i)
            flatfile.write("\nTIME %s [%s] html_tls_group_table" % (
                chain_ntls, misc.timestamp()))
            flatfile.write("\nCCCC %s Analysis of TLS Group %s Chain Segments" % (chain_ntls, ntls))
            flatfile.write("\nDATA %s SEGMENT: %s" % (chain_ntls, tls.display_label()))
            flatfile.write("\nDATA %s RESIDUES: %d" % (chain_ntls, tls.num_residues()))
            flatfile.write("\nDATA %s ATOMS: %d" % (chain_ntls, tls.num_atoms()))
            flatfile.write("\nDATA %s MEAN_B: %.1f" % (chain_ntls, tls.mean_b()))
            flatfile.write("\nDATA %s Brmsd: %.2f" % (chain_ntls, stddev))
            flatfile.write("\nDATA %s MEAN_ANISO: %.2f" % (chain_ntls, tls.mean_anisotropy()))
            flatfile.write("\nDATA %s RMSD_B: %.2f" % (chain_ntls, tls.rmsd_b))
            flatfile.write("\nDATA %s T^rB: %s" % (chain_ntls, t_data))
            flatfile.write("\nDATA %s EVAL_L: %.2f,%.2f,%.2f" % (chain_ntls, L1, L2, L3))
            flatfile.write("\nDATA %s TLS_MEAN_B: %.1f" % (chain_ntls, tls.tls_mean_b()))
            flatfile.write("\nDATA %s TLS_MEAN_ANISO: %.2f" % (chain_ntls, tls.tls_mean_anisotropy()))

            ## XXX: Are the following redundant?
            flatfile.write("\nDATA %s CHECK_RMSD_B: %.2f" % (chain_ntls, cpartition.rmsd_b()))
            flatfile.write("\nDATA %s CHECK_RESIDUAL: %.2f" % (chain_ntls, cpartition.residual()))
            flatfile.write("\nTLST %s" % t_data)

            flatfile.close()
            i += 1 ## increment segment within the partition
        ##</FLATFILE>
        ##======================================================================

        ## "Analysis of TLS Group n Chain Segments" table
        ## FIXME: Why are the 'tls.rmsd_b'/"RMSD B" values different from the
        ## off-diagonal matrix values? The "RMSD B Values of Combined TLS Groups"
        ## values.
        l += ['<td align="center" valign="middle"><img src="%s" alt="%s"/></td>' % (
             cpath, tls.color.name),
             ## Input Structure ================================================
             '<td>%s</td>'    % (tls.display_label()),       ## "Segment"
             '<td>%d</td>'    % (tls.num_residues()),        ## "Residues"
             '<td>%d</td>'    % (tls.num_atoms()),           ## "Atoms"
             '<td>%5.1f</td>' % (tls.mean_b()),              ## "<B>"
             '<td>%5.2f</td>' % (stddev),                    ## "Brmsd", 2008-04-15
             '<td>%4.2f</td>' % (tls.mean_anisotropy()),     ## "<Aniso>"
             ## TLS Predictions ================================================
             '<td>%5.2f</td>' % (tls.rmsd_b),                ## "RMSD B"
             '<td>%s</td>'    % (t_data),                    ## "T^rB"
             '<td>%5.2f, %5.2f, %5.2f</td>' % (L1, L2, L3),  ## "eval(L) DEG^2"
             '<td>%5.1f</td>' % (tls.tls_mean_b()),          ## "<B>"
             '<td>%4.2f</td>' % (tls.tls_mean_anisotropy()), ## "<Aniso>"
             '</tr>']

    l.append('</table>')
    return "".join(l)


_REPORT_CSS_STYLES = """\
BODY {background-color:white; margin:2% 5% 0 5%; border:2% 5% 0 5%;}
a.structimage {text-decoration:none;}
img {border:0;padding:4px 0 0 0;}
span.small {font-size:0.8em;font-weight:normal;text-align:center;}
table.title {border:0;width:100%;background-color:#eee;}
td.title {font-size:0.7em;}
table.report_globals {
    padding:3px; width:75%; background-color:#eee; 
    font-size:small; border:0 
    }
tr.report_globals { background-color:#ddd; }
table.tls_segments {
    width:100%; border:0;
    background-color:#eeeeee; font-size:x-small;
    }
table.matrix {
    background-color:#bbbbff;
    border:thin solid black;
    }
td.matrix {
    padding:5px;font-size:x-small;
    }
td.l {text-align:left;}
td.c {text-align:center;}
td.r {text-align:right;}
td.l, td.c, td.r { width:33%; }
h2 {text-align:center;}
ul { list-style-type:none; margin:0px; padding:0px; }
li { line-height:20px }

p.captions {
    font-size:small; text-align:left;
    }
p.gnuplot_captions {
    padding:2%; background-color:#eeeeee;
    border:thin dashed black;
    }
p.notes {
    width:70%; padding:10px;
    text-align:center;
    font-size:large bold;
    border:thin solid #f00;
    background-color:#faa;
    }
.prog-border {
    margin:0;padding:0;height:100%;width:150px;
    background:#fff;border:1px solid black;text-align:left;
    }
.prog-bar {
    margin:0px;padding:0;height:10px;
    background:#8FBC8F;
    }
"""

_REPORT_CSS_PRINT_STYLES = """\
body { margin:0; font-size:9pt; }
a.imageview { list-style:none; border:0; }
img.structimage { vertical-align:bottom; border:0; }
table.tls_segments { margin:0; width:100%; font-size:9pt; }
div.links, a.links, span.print { display:none; }
div.clear { clear:both; }
"""

class Report(object):
    """Base class of HTML Report generating objects.
    """
    def html_head(self, title):
        """Header for all HTML pages.
        """
        l = ['<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">\n',
             '<html>',
             '<head>',
             '<title>%s</title>\n' % (title),
             '<style type="text/css" media="screen">',
             '<!-- ',
             '%s' % _REPORT_CSS_STYLES,
             '-->',
             '<style type="text/css" media="print">',
             '<!-- ',
             '%s' % _REPORT_CSS_PRINT_STYLES,
             '-->',
             '</style>\n',
             '</head>',
             '<body>\n']

        return "".join(l)

    def html_title(self, title):
        """Title for all HTML pages.
        """
        l  = ['<table class="title"><tr>\n',
              '<td class="l title">%s</td>' % (misc.start_time()),
              '<td class="c title">JobID: %s</td>' % (
                  conf.globalconf.job_id),
              '<td class="r title">TLSMD Version %s</td>\n' % (
                  const.VERSION),
              '</tr></table>\n',
              '<h2>%s</h2><br/>\n' % (title)]

        return "".join(l)

    def html_foot(self):
        """Footer for all HTML pages.
        """
        l  = ['<table class="title"><tr>',
              '<td class="l title">%s</td>' % (misc.start_time()),
              '<td class="c title">JobID: %s</td>' % (
                  conf.globalconf.job_id),
              '<td class="r title">TLSMD Version %s Released %s</td>' % (
                  const.VERSION, const.RELEASE_DATE),
              '</tr></table>',
              '</body></html>\n']

        return "".join(l)


class HTMLSummaryReport(Report):
    """Create a summary HTML report.
    """
    def __init__(self, tlsmd_analysis):
        Report.__init__(self)

        self.tlsmd_analysis = tlsmd_analysis
        self.struct = tlsmd_analysis.struct
        self.struct_id = tlsmd_analysis.struct_id
        self.struct_path = "%s.pdb" % (self.struct_id)
        self.job_id = conf.globalconf.job_id

        self.flatfile_name = "%s.dat" % self.job_id

        self.page_multi_chain_alignment  = None
        self.pages_chain_motion_analysis = []
        self.page_refinement_prep        = None

    def write_summary(self, report_dir):
        """Write out the TLSMD report to the given directory.
        """
        ## create new directory and move into it
        old_dir = os.getcwd()
        if not os.path.isdir(report_dir):
            os.mkdir(report_dir)
        os.chdir(report_dir)

        analysis_dir = os.getcwd()

        self.flatfile_globals()

        ## These are the Jmol Java files needed for the viewer and animator
        shutil.copy(conf.JMOL_PATH + "/JmolApplet.jar", analysis_dir)
        shutil.copy(conf.JMOL_PATH + "/Jmol.jar", analysis_dir)
        shutil.copy(conf.JMOL_PATH + "/Jmol.js", analysis_dir)

        ## This is a script that allows the user to animate a given partition
        ## of a given chain into the 8 phases of its associated libration.
        try:
            shutil.copy(conf.PDB_ANIMATE_SCRIPT, analysis_dir)
        except:
            console.stdoutln("NOTE: Could not find %s" % (
                conf.PDB_ANIMATE_SCRIPT))

        ## Create preliminary summary.png plot
        min_residuals = []
        max_residuals = []
        max_ntls = []
        list_stddev = []
        residual_log = open(conf.RESIDUALS_LOG_FILE, "a+")
        for chain in self.tlsmd_analysis.iter_chains():
            #self.pre_tls_chain_optimization(chain)

            residual_log.write("%s %s " % (self.job_id, chain.chain_id))

            ## Collect data for logfile and Berkeley DB
            segs = 0
            tmp_min = 100.0
            tmp_max = 0.0
            for ntls, cpartition in chain.partition_collection.iter_ntls_chain_partitions():
                ##fields: cpartition.rmsd_b(), cpartition.residual())
                segs += 1 ## for max seg reached per chain

                ## log all residuals (i.e., the residual for each partition)
                residual_log.write("%.2f " % cpartition.rmsd_b())

                ## Roundabout way to find min/max values
                if float(cpartition.rmsd_b()) >= tmp_max:
                    tmp_max = cpartition.rmsd_b()
                if float(cpartition.rmsd_b()) <= tmp_min:
                    tmp_min = cpartition.rmsd_b()

                ## Calculate the stddev for all temperature factors in a given
                ## chain (for the first partition only)
                if int(ntls) == 1:
                    for tls in cpartition.iter_tls_segments():
                        tmp_temp_factor = []
                        for atm, Utls in tls.tls_group.iter_atm_Utls():
                            tmp_temp_factor.append(atm.temp_factor)
                        list_stddev.append("%s:%.2f" % (
                            chain.chain_id, numpy.std(tmp_temp_factor)))

            residual_log.write("\n")
            min_residuals.append("%s:%.2f" % (chain.chain_id, float(tmp_min)))
            max_residuals.append("%s:%.2f" % (chain.chain_id, float(tmp_max)))
            max_ntls.append("%s:%s" % (chain.chain_id, segs))

            ## add tables for all TLS group selections using 1 TLS group
            ## up to max_ntls
            for ntls in chain.partition_collection.iter_ntls():
                gp = gnuplots.LSQR_vs_TLS_Segments_Plot(chain)
                ## maybe this will help with the memory problems...
                import gc
                gc.collect()
        plot = gnuplots.LSQR_vs_TLS_Segments_All_Chains_Plot(self.tlsmd_analysis)

        residual_log.close()

        ## This will store the initial + final residual for each chain in the
        ## flatfile, as well as the stddev(Bfact) for each chain.
        initial_residuals = ";".join(max_residuals)
        final_residuals   = ";".join(min_residuals)
        stddev_bfact      = ";".join(list_stddev)
        chain_max_segs    = ";".join(max_ntls)
        console.stdoutln("RESIDUALS: INITIAL = %s" % initial_residuals)
        console.stdoutln("RESIDUALS: FINAL = %s" % final_residuals)
        console.stdoutln("STDDEV_BFACT: %s" % stddev_bfact)
        console.stdoutln("MAX_SEGS: %s" % chain_max_segs)

        ##======================================================================
        ##<FLATFILE>
        flatfile = open(self.flatfile_name, "a+")
        flatfile.write("\nGENR INITIAL_RESIDUALS: %s" % initial_residuals)
        flatfile.write("\nGENR FINAL_RESIDUALS: %s" % final_residuals)
        flatfile.write("\nGENR STDDEV_BFACT: %s" % stddev_bfact)
        flatfile.write("\nGENR MAX_SEGS: %s" % chain_max_segs)
        flatfile.close()
        ## </FLATFILE>
        ##======================================================================

        self.write_summary_index()

        ## change back to original directory
        os.chdir(old_dir)

    def write_summary_index(self):
        """Writes the summary index.html file of the report.
        """
        fil = open("index.html", "w")
        fil.write(self.html_summary_index())
        fil.close()
        console.stdoutln("HTML: Saving summary index.html")

    def html_summary_index(self):
        """Generate and returns the HTML string for the summary index.html
        file of the report.
        """
        title = "TLSMD Thermal Parameter Analysis of Structure %s" % (
            self.struct_id)

        l = [self.html_head(title),
             self.html_title(title),

             ## link back to job summary page, 2009-05-26
             '<center><a href="%s?page=explore&amp;job_id=%s">' % (
                 conf.WEBTLSMD_URL, self.job_id),
             'Back to job summary page</a></center>',

             ## OPTIMIZATION PARAMETERS
             self.html_globals(),
             '<br/>\n',

             ## MOTION ANALYSIS
             '<center>',
             '<h3>TLS Partitions and Motion Analysis of Individual Chains</h3>',
             '</center>\n',
             '<table><tr><td valign=top>',
             '<p class="captions">%s</p>\n' % (captions.MOTION_ANALYSIS_TEXT)]

        l += ['<center><p class="notes">',
             'Your job is still running and the analysis is incomplete. ',
             'This is a summary page.</p>\n']

        ## progress bar
        try:
            prog_file = open("../progress", 'r')
            progress = int(float(prog_file.read().strip())*100)
            prog_file.close()
        except:
            progress = 20
        l += ['<p><div class="prog-border">',
             '<div class="prog-bar" style="width: %s%%;">' % progress,
             '</div></div>%s%% Complete</p>\n' % progress,
             '</center>']

        l += ['</td><td valign=top><img src="summary.png" /></td></tr></table>\n']

        """
        l +=['<br/>',
             ## MULTI CHAIN ALIGNMENT
             '<center><h3>Multi-Chain TLS Group Alignment</h3></center>',
             '<p class="captions">%s</p>' % (captions.MULTI_CHAIN_ALIGNMENT_TEXT)]

        if self.page_multi_chain_alignment != None:
            l.append('<p><a href="%s">%s</a></p>' % (
                self.page_multi_chain_alignment["href"], 
                self.page_multi_chain_alignment["title"]))
        else:
            l.append('<p><u>Only one chain was analyized in this structure, ')
            l.append('so the multi-chain alignment analysis was not performed.</u></p>')

        if self.page_refinement_prep is not None:
            l +=['<br/>',
                ## REFINEMENT PREP
                '<center><h3>Generate input files for multigroup TLS Refinement</h3></center>',
                '<p class="captions">%s</p>' % (
                    captions.REFINEMENT_PREP_TEXT),
                '<p><a href="%s">%s</a></p>' % (
                self.page_refinement_prep["href"], 
                self.page_refinement_prep["title"])]
        """

        l += [self.html_foot()]

        return "".join(l)

    def html_globals(self):
        """Output a HTML table displaying global TLSMD settings.
        """
        ## NOTE: class HTMLSummaryReport()
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

        l = ['<center><table class="report_globals">\n']

        if self.struct.title:
            l += ['<tr class="report_globals">',
                  '<td>Title</td><td><b>%s</b></td></tr>\n' % (
                  self.struct.title)]

        if self.struct.header:
            l += ['<tr class="report_globals">',
                  '<td>Heading Summary</td><td><b>%s</b></td></tr>\n' % (
                  self.struct.header)]

        if self.struct.experimental_method:
            l += ['<tr class="report_globals">',
                  '<td>Experimental Method</td><td><b>%s</b></td></tr>\n' % (
                  self.struct.experimental_method)]

        l +=['<tr class="report_globals">',
             '<td>Temperature Factors</td><td><b>%s</b></td></tr>\n' % (
             tls_model),
             '<tr class="report_globals">',
             '<td>Minimum TLS Segment Length</td>',
             '<td><b>%s Residues</b></td></tr>\n' % (
             conf.globalconf.min_subsegment_size),
             '<tr class="report_globals">',
             '<td>Atoms Analyzed</td>',
             '<td><b>%s</b></td></tr>\n' % (
             conf.globalconf.include_atoms),
             '</table>',
             '</center>\n']

        return "".join(l)

    def flatfile_globals(self):
        """Store globals in flatfile
        """
        ## NOTE: class HTMLSummaryReport()
        ##======================================================================
        ##<FLATFILE>
        flatfile = open(self.flatfile_name, "a+")

        flatfile.write("\nGLOB JOB_ID: %s" % conf.globalconf.job_id)
        flatfile.write("\nGLOB STRUCT_ID: %s" % conf.globalconf.struct_id)
        flatfile.write("\nGLOB TLS_MODEL: %s" % conf.globalconf.tls_model)
        flatfile.write("\nGLOB WEIGHT_MODEL: %s" % conf.globalconf.weight_model)
        flatfile.write("\nGLOB INCLUDE_ATOMS: %s" % conf.globalconf.include_atoms)
        flatfile.write("\nGLOB MIN_SUBSEGMENT_SIZE: %s" % conf.globalconf.min_subsegment_size)
        flatfile.write("\nGLOB ADP_PROB: %s" % conf.globalconf.adp_prob)
        flatfile.write("\nGLOB ADP_SMOOTHING: %s" % conf.globalconf.adp_smoothing)
        flatfile.write("\nGLOB NPARTS: %s" % conf.globalconf.nparts)
        flatfile.write("\nGLOB USE_SVG: %s" % conf.globalconf.use_svg)
        flatfile.write("\nGLOB SKIP_HTML: %s" % conf.globalconf.skip_html)
        flatfile.write("\nGLOB SKIP_JMOL: %s" % conf.globalconf.skip_jmol)
        flatfile.write("\nGLOB GENERATE_JMOL_VIEW: %s" % conf.globalconf.generate_jmol_view)
        flatfile.write("\nGLOB GENERATE_JMOL_ANIMATE: %s" % conf.globalconf.generate_jmol_animate)
        flatfile.write("\nGLOB GENERATE_HISTOGRAM: %s" % conf.globalconf.generate_histogram)
        flatfile.write("\nGLOB CROSS_CHAIN_ANALYSIS: %s" % conf.globalconf.cross_chain_analysis)
        flatfile.write("\nGLOB RECOMBINATION: %s" % conf.globalconf.recombination)
        flatfile.write("\nGLOB TARGET_STRUCT_PATH: %s" % conf.globalconf.target_struct_path)
        flatfile.write("\nGLOB TARGET_STRUCT_CHAIN_ID: %s" % conf.globalconf.target_struct_chain_id)

        flatfile.close()
        ##</FLATFILE>
        ##======================================================================


class HTMLReport(Report):
    """Create a thorough HTML report it its own subdirectory.
    """
    def __init__(self, tlsmd_analysis):
        Report.__init__(self)

        self.tlsmd_analysis = tlsmd_analysis
        self.struct = tlsmd_analysis.struct
        self.struct_id = tlsmd_analysis.struct_id
        self.struct_path = "%s.pdb" % (self.struct_id)
        self.job_id = conf.globalconf.job_id

        self.page_multi_chain_alignment  = None
        self.pages_chain_motion_analysis = []
        self.page_refinement_prep        = None

        self.flatfile_name = "%s/%s/%s.dat" % (
            conf.TLSMD_WORK_DIR, conf.globalconf.job_id,
            conf.globalconf.job_id)

        self.orient = {}
        self.r3d_header_file = None

    def write(self, report_dir):
        """Write out the TLSMD report to the given directory.
        """
        ## class HTMLReport()
        ## create new directory and move into it
        console.stdoutln("REPORT DIR: %s" % report_dir)
        old_dir = os.getcwd()
        if not os.path.isdir(report_dir):
            os.mkdir(report_dir)
        os.chdir(report_dir)

        self.write_cwd() ## NOTE: This is the very last step of TLSMD

        ## change back to original directory
        os.chdir(old_dir)

    def init_colors(self):
        """Generated the self.colors dictionary of colors for the report,
        and also writes thumbnail .png images of all the colors.
        """
        ## class HTMLReport()
        thumbnail_dir = "colors"

        ## make the thumbnail subdirectory
        if not os.path.isdir(thumbnail_dir):
            os.mkdir(thumbnail_dir)

        ## clear the colors list
        self.colors = []

        ## skip the first two colors, which are black/white
        ## NOTE: Generates only 52 colors. If there are more TLS partitions to
        ## segment, the code will reuse the colors by looping around.
        for i in xrange(len(Colors.COLORS)):
            cname, rgbf = Colors.COLORS[i]
            color = ColorInfo(i, cname, rgbf, thumbnail_dir)
            self.colors.append(color)

        ## assign a unique color to each tls group in a
        ## chain spanning set of tls groups
        for chain in self.tlsmd_analysis.iter_chains():
            for cpartition in chain.partition_collection.iter_chain_partitions():

                for tlsi, tls in enumerate(cpartition.iter_tls_segments()):
                    if tls.method == "TLS":
                        tls.set_color(self.get_tls_color(tlsi))
                    else:
                        tls.set_color(self.colors[0])

    def get_tls_color(self, tls_index):
        """Returns the color dict description for a TLS segment of the
        given index, starting at 0.
        """
        ## class HTMLReport()
        ## skip the first two colors; they are black and white
        i = tls_index + 2
        return self.colors[i % 50] ## switched to mod(50) for wrap-around

    def write_cwd(self):
        """Write out all the files in the report.
        """
        ## class HTMLReport()

        ## write a local copy of the Structure (adds "HEADER" + "CRYST1" data
        ## even if nothing existed in the original, user-contributed PDB file)
        FileIO.SaveStructure(fil = self.struct_path, struct = self.struct)

        ## generate small .png images so they can be placed in the
        ## TLS group tables to identify the TLS group tabular data
        ## with the generated visualization
        self.init_colors()

        ## a report page comparing the tls group segments of all
        ## chains against each other
        try:
            console.stdoutln("Writing multi-chain alignment")
            self.write_multi_chain_alignment()
        except:
            console.stdoutln("        Error: Couldn't deal with multi-chain alignment")
            print console.formatExceptionInfo()
            pass

        ## Progress tracking:
        ##    - assume this portion of the run occupies 0.5 -> 1.0 of the
        ##      total time
        progress = 0.5

        ## write out all TLSGraph reports
        for chain in self.tlsmd_analysis.iter_chains():
            self.write_tls_chain_optimization(chain)

            ## Track progress
            progress += 0.5/self.tlsmd_analysis.num_chains()
            progress_report = open("../progress","w+")
            print >> progress_report, progress
            progress_report.close()

        self.write_refinement_prep()

        ## write out index page
        self.write_index()

        ## update summary index.html to show completed status
        #self.summary_file_update()

    def write_index(self):
        """Writes the main index.html file of the report.
        """
        ## class HTMLReport()
        fil = open("index.html","w")
        fil.write(self.html_index())
        fil.close()
        console.stdoutln("HTML: Saving main index.html")

    def summary_file_update(self):
        """Updates the summary index.html to reflect status.
        """
        ## class HTMLReport()

        link = '<a href="%s/jobs/%s/ANALYSIS/index.html">' % (
            conf.TLSMD_PUBLIC_URL, conf.globalconf.job_id)
        data = ""
        lines = open("index.html", 'r').readlines()
        for line in lines:
            if re.match(r'(.*job is) still running and the analysis is incomplete(.*)', line):
                data += re.sub(r'(.*job is) still running and the analysis is incomplete(.*)', '\\1 %scompleted</a>\\2' % link, line)
            elif re.match(r'.*#8FBC8F;width: 50%;"></div></div>50% (Complete.)', line):
                data += re.sub(r'(.*#8FBC8F;width:) 50%(;"\>\<\/div\>\<\/div)\>50% (.*)', '\\1 100%\\2>100% \\3', line)
            else:
                data += line

        outfile = open("index.html", 'w')
        outfile.write(data)
        outfile.close()

    def html_index(self):
        """Generate and returns the HTML string for the main index.html
        file of the report.
        """
        ## class HTMLReport()
        title = "TLSMD Thermal Parameter Analysis of Structure %s" % (
            self.struct_id)

        l = [self.html_head(title),
             self.html_title(title),

             ## link back to job summary page, 2009-05-26
             '<center><a href="%s?page=explore&amp;job_id=%s">Back to job summary page</a></center>' % (
                 conf.WEBTLSMD_URL, self.job_id),

             ## OPTIMIZATION PARAMETERS
             self.html_globals(),
             '<br/>',

             ## MOTION ANALYSIS
             '<center>',
             '<h3>TLS Partitions and Motion Analysis of Individual Chains</h3>',
             '</center>',
             '<table><tr><td valign=top>',
             '<p class="captions">%s</p>' % (captions.MOTION_ANALYSIS_TEXT)]

        for xdict in self.pages_chain_motion_analysis:
            l.append('<p><a href="%s">%s</a></p>' % (xdict["href"], xdict["title"]))

        l +=['</td><td valign=top><img src="summary.png"></td></tr></table>']

        l +=['<br/>',
             ## MULTI CHAIN ALIGNMENT
             '<center><h3>Multi-Chain TLS Group Alignment</h3></center>',
             '<p class="captions">%s</p>' % (
                 captions.MULTI_CHAIN_ALIGNMENT_TEXT)]

        if self.page_multi_chain_alignment!=None:
            l.append('<p><a href="%s">%s</a></p>' % (
                self.page_multi_chain_alignment["href"], 
                self.page_multi_chain_alignment["title"]))
        else:
            l.append('<p><u>Only one chain was analyized in this structure, ')
            l.append('so the multi-chain alignment analysis was not performed.</u></p>')

        l +=['<br/>',
             ## REFINEMENT PREP
             '<center><h3>Generate input files for multigroup TLS Refinement</h3></center>',
             '<p class="captions">%s</p>' % (captions.REFINEMENT_PREP_TEXT),
             '<p><a href="%s">%s</a></p>' % (
                 self.page_refinement_prep["href"], 
                 self.page_refinement_prep["title"]),
             self.html_foot()]

        return "".join(l)

    def html_globals(self):
        """Output a HTML table displaying global TLSMD settings.
        """
        ## class HTMLReport()
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

        l = ['<center><table class="report_globals">']

        if self.struct.title:
            l += ['<tr class="report_globals"><td>Title</td>',
                  '<td><b>%s</b></td></tr>' % (self.struct.title)]

        if self.struct.header:
            l += ['<tr class="report_globals"><td>Heading Summary</td>',
                  '<td><b>%s</b></td></tr>' % (self.struct.header)]

        if self.struct.experimental_method:
            l += ['<tr class="report_globals"><td>Experimental Method</td>',
                  '<td><b>%s</b></td></tr>' % (
                  self.struct.experimental_method)]

        l +=['<tr class="report_globals"><td>Temperature Factors</td>',
             '<td><b>%s</b></td></tr>' % (tls_model),
             '<tr class="report_globals">',
             '<td>Minimum TLS Segment Length</td>',
             '<td><b>%s Residues</b></td></tr>' % (
                 conf.globalconf.min_subsegment_size),
             '<tr class="report_globals"><td>Atoms Analyzed</td>',
             '<td><b>%s</b></td></tr>' % (conf.globalconf.include_atoms),
             '</table>',
             '</center>']

        return "".join(l)

    def write_tls_chain_optimization(self, chain):
        """Writes the HTML report analysis of a single TLS graphed chain.
        """
        ## class HTMLReport()
        path  = "%s_CHAIN%s_ANALYSIS.html" % (self.struct_id, chain.chain_id)
        title = "Chain %s TLS Analysis" % (chain.chain_id)

        self.pages_chain_motion_analysis.append(
            {"title": title,
             "href":  path })

        fil = open(path, "w")
        fil.write(self.html_tls_chain_optimization(chain))
        fil.close()
        console.stdoutln("HTML: Saving %s" % path)

    def html_tls_chain_optimization(self, chain):
        """Generates and returns the HTML string report analysis of a
        single TLS graphed chain.
        """
        ## class HTMLReport()
        title = "Chain %s TLS Analysis of %s" % (chain.chain_id, self.struct_id)

        l  = [self.html_head(title),
              self.html_title(title),
              '<center><a href="index.html">Back to Index</a></center>',
              '<br/>' ]

        ## if there were no valid chain configurations found
        ## then write out a useful error message
        if chain.partition_collection.num_chain_partitions() == 0:
            l += ['<p>%s</p>' % (NO_VALID_CONFIGURATIONS),
                  self.html_foot()]
            return "".join(l)

        ## TLS Segments vs. Residual
        l += [self.html_chain_lsq_residual_plot(chain),
              self.html_chain_alignment_plot(chain)]

        ## orient the structure with the super-spiffy orientation algorithm
        ## which highlights the chain we are examining
        try:
            self.orient = calc_orientation(self.struct, chain)
            console.debug_stdoutln("[%s] Raster3D: calculating orientation" % (chain.chain_id))
            #    R         = ori["R"],
            #    cor       = ori["centroid"],
            #    zoom      = ori["hzoom"],
            #    near      = ori["near"],
            #    far       = ori["far"],
            #    width     = ori["pwidth"],
            #    height    = ori["pheight"],
            #    bg_color  = "White")
        except:
            console.stdoutln("     Warning: failed to find orientation for graphics output")
            print console.formatExceptionInfo()
            pass

        ## Generate Raster3D header for a given chain
        if conf.RENDER_SKIP:
            self.r3d_header_file = None
        else:
            try:
                self.r3d_header_file = self.generate_raster3d_header(chain)
            except:
                console.stdoutln("     Warning: failed to generate Raster3D header")
                print console.formatExceptionInfo()

        ## add tables for all TLS group selections using 1 TLS group
        ## up to max_ntls
        for ntls in chain.partition_collection.iter_ntls():
            console.stdoutln("=" * 80) ## Entering new segment
            try:
                tmp = self.html_tls_graph_path(chain, ntls)
                if tmp != None:
                    l.append(tmp)
            except:
                print console.formatExceptionInfo()

            ## maybe this will help with the memory problems...
            import gc
            gc.collect()

        l.append(self.html_foot())
        return "".join(l)

    def html_chain_lsq_residual_plot(self, chain):
        """Generates the Gnuplot/PNG image plot, and returns the HTML
        fragment for its display in a web page.
        """
        ## class HTMLReport()
        ## Least SQuare Residual (LSQR) vs. Number of TLS Segments
        gp = gnuplots.LSQR_vs_TLS_Segments_Plot(chain)

        title = 'Chain %s Optimization Residual' % (chain.chain_id)

        l = ['<center>',
             gp.html_markup(title, captions.LSQR_CAPTION),
             '</center>']

        return "".join(l)

    def html_chain_alignment_plot(self, chain):
        """generate a plot comparing all segmentations
        """
        ## class HTMLReport()
        plot = sequence_plot.TLSSegmentAlignmentPlot()

        for ntls, cpartition in chain.partition_collection.iter_ntls_chain_partitions():
            plot.add_tls_segmentation(cpartition)

        ## create filename for plot PNG image file
        plot_file = "%s_CHAIN%s_ALIGN.png" % (self.struct_id, chain.chain_id)

        plot.plot(plot_file)

        l = ['<center><h3>TLS Partition Segment Alignment of Chain %s</h3></center>' % (
                 chain.chain_id),
             '<center>',
             '<table border="0" style="background-color:#dddddd">',
             '<tr><th># of TLS<br/>Groups</th>',
             '<th>Segment/Sequence Alignment</th></tr>',
             '<tr>',

             '<td align="right" valign="top">',
             '<p style="font-size:xx-small; margin-top:%dpx">' % (
                 plot.border_width)]

        l.append('<ul>')
        for ntls, cpartition in chain.partition_collection.iter_ntls_chain_partitions():
            l.append('<li><a href="#NTLS%d">%d</a><br/>' % (ntls, ntls))
        l.append('</ul>')

        l +=['</td>',
             '<td><img src="%s" alt="Sequence Alignment Plot"></td>' % (
                 plot_file),
             '</tr>',
             '</table>',
             '</center>',
             '<p>%s</p>' % (captions.SEG_ALIGN_CAPTION)]

        return "".join(l)

    def html_tls_graph_path(self, chain, ntls):
        """Generates the HTML table describing the path (set of tls groups)
        for the given number of segments(h, or ntls)
        """
        ## class HTMLReport()
        cpartition = chain.partition_collection.get_chain_partition(ntls)
        if cpartition == None:
            return None

        #self.write_tls_pdb_file(chain, cpartition) ## write out PDB file

        ## Find generate Jmol viewer/animate in globals
        try:
            job_id = conf.globalconf.job_id
            jmol_view_toggle = conf.globalconf.generate_jmol_view
            jmol_animate_toggle = conf.globalconf.generate_jmol_animate
        except:
            jmol_view_toggle = ""
            jmol_animate_toggle = ""
            console.stdoutln("     Warning: couldn't find Jmol-viewer settings")
            print console.formatExceptionInfo()
            pass

        ## Jmol Viewer Script
        if conf.globalconf.generate_jmol_view == False:
            jmol_file = ""
            console.stdoutln("NOTE: Skipping Jmol-viewer section")
        else:
            try:
                jmol_file = self.jmol_html(chain, cpartition)
            except:
                ## EAM FIXME:  But really something must have gone wrong 
                ## before this point.
                console.stdoutln("     Warning: Jmol setup failed in jmol_html")
                print console.formatExceptionInfo()
                jmol_file = ""
                pass

        ## Jmol Animation Script
        if conf.globalconf.generate_jmol_animate == False:
            console.stdoutln("NOTE: Skipping Jmol-animation section")
            jmol_animate_file = ""
            raw_r3d_file, r3d_body_file =\
                self.generate_raw_backbone_file(chain, cpartition)
        else:
            try:
                jmol_animate_file, raw_r3d_file, r3d_body_file =\
                    self.jmol_animate_html(chain, cpartition)
            except:
                ## EAM FIXME:  But really something must have gone wrong 
                ## before this point.
                jmol_animate_file = ""
                console.stdoutln("     Warning: Jmol setup failed in jmol_animate_html")
                print console.formatExceptionInfo()
                pass

        ## New Raster3D section (bypassing pymmlib)
        if conf.RENDER_SKIP:
            png_file = ""
            pml_file = ""
        else:
            try:
                basename = "%s_CHAIN%s_NTLS%d" % (
                    self.struct_id, chain.chain_id, 
                    cpartition.num_tls_segments())
                png_file = "%s.png" % (basename)
                pml_file = "" ## never used

                if os.path.isfile("../struct.r3d"):
                    ## only run if *.r3d files exist
                    gen_r3d_body_cmd = "%s < %s > %s 2> /dev/null" % (
                        conf.TLSANIM2R3D, raw_r3d_file, r3d_body_file)
                    os.system(gen_r3d_body_cmd)

                    if os.path.isfile("../bases.r3d"):
                        ## there are nucleic acids, so
                        ## cat header.r3d static.r3d animate.r3d grey.r3d \
                        ## bases.r3d sugars.r3d
                        render_cmd = "cat %s ../struct.r3d %s %s ../bases.r3d ../sugars.r3d | %s > %s 2> /dev/null" % (
                            self.r3d_header_file, r3d_body_file, 
                            conf.GREY_R3D_FILE, conf.RENDER, png_file)
                    else:
                        ## there are _not_ any nucleic acids, so
                        render_cmd = "cat %s ../struct.r3d %s | %s > %s 2> /dev/null" % (
                            self.r3d_header_file, r3d_body_file,
                            conf.RENDER, png_file)

                    start = time.time()
                    os.system(render_cmd)
                    elapsed = (time.time() - start)
                    console.stdoutln("[%s,%s] RENDERING TIME for %s: %.2f s" % (
                        chain.chain_id, cpartition.num_tls_segments(), 
                        png_file, elapsed))
            except:
                raw_r3d_file = ""
                r3d_body_file = ""
                png_file = ""
                console.stdoutln("     Warning: failure in rendering PNG file")
                print console.formatExceptionInfo()
                pass

        ## Refmac/Phenix files
        if conf.REFMAC_SKIP:
            tlsout_file = ""
            phenixout_file = ""
            console.stdoutln("NOTE: Skipping Refmac/Phenix section")
        else:
            try:
                ## tlsout (for refmac) and phenix files
                tlsout_file = self.write_tlsout_file(chain, cpartition)
                phenixout_file = self.write_phenixout_file(chain, cpartition)
            except:
                tlsout_file = ""
                phenixout_file = ""
                console.stdoutln("     Warning: failed to create Refmac/Phenix files")
                print console.formatExceptionInfo()
                pass

        ## detailed analysis of all TLS groups
        try:
            ntls_analysis = self.chain_ntls_analysis(chain, cpartition)
        except:
            console.stderr("ERROR: Analysis of this partition failed")
            print console.formatExceptionInfo()
            pass

        ## BMean Plot
        ntls_analysis.bmean_plot.width = conf.BMEAN_PLOT_WIDTH
        ntls_analysis.bmean_plot.height = conf.BMEAN_PLOT_HEIGHT
        ntls_analysis.bmean_plot.tls_group_titles = conf.BMEAN_PLOT_GROUP_TITLES
        ntls_analysis.bmean_plot.output_png()

        l = ['<hr/>',

             '<center style="page-break-before:always; font-size:small">',
             '<a name="NTLS%d" style="font-size:large">' % (ntls),
             'Summary of Optimal TLS Group Partition using %d Groups' % (ntls),
             '</a>',
             '<br/>',

             ## navigation links
             '<a href="%s">Full TLS Partition Analysis</a>' % (
                 ntls_analysis.url),

             '&nbsp;&nbsp;&nbsp;&nbsp;',

             '<a href="." onClick="window.open(&quot;%s&quot;,&quot;&quot;,&quot;' % (
                 jmol_file),
             'width=%d,height=%d,screenX=10,screenY=10,left=10,top=10&quot;);' % (
                 conf.JMOL_SIZE, conf.JMOL_SIZE),
             'return false;">View with Jmol</a>',

             '&nbsp;&nbsp;&nbsp;&nbsp;',

             '<a href="." onClick="'\
             'window.open('\
             '&quot;%s&quot;,'\
             '&quot;&quot;,'\
             '&quot;width=%d,height=%d,screenX=10,'\
             'screenY=10,left=10,top=10&quot;);'\
             'return false;">Animate Screw Displacement with Jmol</a>' % (
                 jmol_animate_file, conf.JMOL_SIZE, conf.JMOL_SIZE),

             '<br/>',
             '<a href="%s">Download TLSOUT File for TLSView</a>' % (tlsout_file),
             '&nbsp;&nbsp;&nbsp;&nbsp;',
             '<a href="%s">Group description for PHENIX</a>' % (phenixout_file),
             '&nbsp;&nbsp;&nbsp;&nbsp;',
             '<a href="%s">Generate PDBIN/TLSIN Files for REFMAC5/PHENIX</a>' % (
             "%s_REFINEMENT_PREP.html" % (self.struct_id)),

             '</center>',

             ## raytraced image
             '<table style="background-color:white" width="100%" border=0>',
             '<tr><th>',
             '<center><a href="%s" class="imageview">' % (png_file),
             '<img src="%s" height="320" alt="structimage"/></a>' % (png_file),
             '<span class="small print">Click on image for expanded view</span>',
             '</center></th></tr>',

             '<tr><th><center>',
             ntls_analysis.bmean_plot.html_link(),
             '</center></th></tr>',
             '</table>',

             ## now the table. Pass "ntls" as well, 2008-04-05
             ## NOTE: This is for the main "ANALYSIS" page.
             html_tls_group_table(ntls, chain, cpartition, detail = True),

             '<br clear="all">']

        return "".join(l)

    #def raster3d_render_tls_graph_path(self, chain, cpartition):
    def generate_raster3d_header(self, chain):
        """Render TLS visualizations using Raster3D.
        """
        ## class HTMLReport()
        r3d_header_file = "%s_CHAIN%s_header.r3d" % (
            self.struct_id, chain.chain_id)
        console.stdoutln("Raster3D: generating header %s..." % r3d_header_file)

        ori = self.orient
        struct_id = self.struct_id

        ## XXX: Size hack: some structures have too many chains,
        ## or are just too large
        show_chain = {}
        for chx in self.struct.iter_chains():
            if chx.chain_id == chain.chain_id:
                show_chain[chx.chain_id] = True
                continue

            if chx.count_amino_acids() >= 100:
                show_chain[chx.chain_id] = False
                continue

            show_chain[chx.chain_id] = True
        ## end size hack

        ## got target chain?
        if self.tlsmd_analysis.target_chain != None:
            for atm in self.tlsmd_analysis.target_chain.iter_all_atoms():
                atm.orig_position = atm.position
                atm.position = chain.target_chain_sresult.transform(atm.position)

        ##<Raster3D Header>=====================================================
        object_list    = []
        width          = 200
        height         = 100
        zoom           = 50
        near           = 0
        far            = 0
        bg_color_rgbf  = (0.0, 0.0, 0.0)
        ambient_light  = 0.2
        diffuse_light  = 1.0
        specular_light = 1.0

        front_clip = near
        back_clip  = far

        ## Raster3D assumes r,g,b triplits are squared
        r, g, b = bg_color_rgbf
        r = r*r
        g = g*g
        b = b*b
        bg_color_rgbf = (r,g,b)

        ## the lighting model for Raster3D is not quite the same as
        ## OpenGL; this conversion gets it close
        total_light = ambient_light + diffuse_light + specular_light

        ambient  = ambient_light  / total_light
        specular = specular_light / total_light
        phong    = 3

        ## initial material state
        object_list.append((8, 0.0, 0, front_clip, back_clip))

        ## now create the header for the render program
        tsz_width   = 16
        tsz_height  = 16

        xtiles = int(round(width  / float(tsz_width)))
        ytiles = int(round(height / float(tsz_height)))

        pixel_width  = xtiles * tsz_width
        pixel_height = ytiles * tsz_height

        ## zoom is the horizontal number of Angstroms shown in the
        ## image, this must be converted to the Raster3D zoom parameter
        ## which is the number of Angstroms of the shortest dimention
        if pixel_width > pixel_height:
            r = float(pixel_height) / float(pixel_width)
            z = zoom * r
        else:
            z = zoom

        R = ori["R"] ## Rotation matrix
        cor = ori["centroid"]

        header_list = [
            "mmLib Generated Raster3D Output",
            "%d %d     tiles in x,y" % (2*ori["pwidth"], 2*ori["pheight"]),
            "0 0       pixels (x,y) per tile",
            "4         anti-aliasing level 4; 3x3->2x2",
            "1 1 1     background",
            "F         no shadows cast",
            "%2d       Phong power" % (phong),
            "0.20      secondary light contribution",
            "%4.2f     ambient light contribution" % (ambient),
            "%4.2f     specular reflection component" % (specular),
            "0.0       eye position(no perspective)",
            "1 1 1     main light source position",
            "1 0 0 0",
            "0 1 0 0",
            "0 0 1 0",
            "%s %s %s %f" % (-cor[0], -cor[1], -cor[2], ori["hzoom"]),
            "3         mixed objects",
            "*        (free format triangle and plane descriptors)",
            "*        (free format sphere descriptors",
            "*        (free format cylinder descriptors)",
            "# Auto-orientation matrix",
            "16",
            "ROTATION",
            "%s %s %s " % (R[0,0], R[0,1], R[0,2]),
            "%s %s %s " % (R[1,0], R[1,1], R[1,2]),
            "%s %s %s " % (R[2,0], R[2,1], R[2,2]),
            "# End of orientation matrix",
            ""
            ]

        r3d_file = open(r3d_header_file, "w")
        r3d_file.write("\n".join(header_list))
        r3d_file.close()
        console.stdoutln("Raster3D: Finish creating r3d_header %s..." % r3d_header_file)

        return r3d_header_file
        ##=====================================================================
        ## NOTE: None of the following is used anymore

        ## turn off axes and unit cell visualization
        gl_struct.glo_update_properties_path("gl_axes/visible", False)
        gl_struct.glo_update_properties_path("gl_unit_cell/visible", False)

        ## setup base structural visualization
        console.debug_stdoutln(">html.py->Raster3D: setup base structural visualization")
        for gl_chain in gl_struct.glo_iter_children():
            if not isinstance(gl_chain, Viewer.GLChain):
                continue

            ## chain is hidden
            if show_chain.get(gl_chain.chain.chain_id, False) == False:
                gl_chain.properties.update(visible = False)
                continue

            if gl_chain.chain.chain_id == chain.chain_id:

                if gl_chain.chain.has_amino_acids():
                    gl_chain.properties.update(
                        lines              = False,
                        trace              = True,
                        trace_radius       = 0.25,
                        trace_color        = "0.80,0.80,0.80")
                elif gl_chain.chain.has_nucleic_acids():
                    gl_chain.properties.update(
                        lines              = False,
                        ball_stick         = True,
                        ball_stick_radius  = 0.25,
                        color              = "0.80,0.80,0.80")
            else:
                if gl_chain.chain.has_amino_acids():
                    gl_chain.properties.update(
                        lines              = False,
                        trace              = True,
                        trace_radius       = 0.25,
                        trace_color        = "0.80,0.80,0.80")
                elif gl_chain.chain.has_nucleic_acids():
                    gl_chain.properties.update(
                        hetatm_visible     = True,
                        trace              = True,
                        trace_radius       = 0.35,
                        trace_color        = "0.80,0.80,0.80",
                        ball_stick         = True,
                        ball_stick_radius  = 0.25,
                        color              = "0.60,0.60,0.70")
                else:
                    gl_chain.properties.update(
                        visible           = True,
                        ball_stick        = True,
                        ball_stick_radius = 0.25,
                        cpk               = False)

        ## add the TLS group visualizations
        has_amino_acids = cpartition.chain.has_amino_acids()
        has_nucleic_acids = cpartition.chain.has_nucleic_acids()

        console.debug_stdoutln(">html.py->Raster3D: add the TLS group visualizations")
        for tls in cpartition.iter_tls_segments():
            if tls.method != "TLS":
                continue

            if self.tlsmd_analysis.target_chain is not None:
                if tls.rmsd_pre_alignment <= 0.8:
                    continue
                if (tls.rmsd_pre_alignment - tls.sresult.rmsd) < 0.5:
                    continue

            tls_name = "TLS_%s" % (tls.filename_label())
            gl_tls_group = TLS.GLTLSGroup(
                oatm_visible       = False,
                side_chain_visible = False,
                hetatm_visible     = True,
                adp_prob           = conf.ADP_PROB,
                L_axis_scale       = 2.0,
                L_axis_radius      = 0.20,
                both_phases        = True,
                tls_group          = tls.tls_group,
                tls_info           = tls.model_tls_info,
                tls_name           = tls_name,
                tls_color          = tls.color.name)

            gl_struct.glo_add_child(gl_tls_group)

            if tls.superposition_vscrew != None:
                gl_tls_group.properties.update(COR_vector = tls.superposition_vscrew)

            ## set width of trace according to the group's translationral tensor trace
            mtls_info = tls.model_tls_info
            tiso = (mtls_info["Tr1_eigen_val"] +\
                    mtls_info["Tr2_eigen_val"] +\
                    mtls_info["Tr3_eigen_val"]) / 3.0

            ## too big usually for good visualization; cheat and scale it down
            radius = 0.30

            if has_amino_acids:
                gl_tls_group.gl_atom_list.properties.update(trace_radius = radius)

            elif has_nucleic_acids:
                gl_tls_group.gl_atom_list.properties.update(
                    oatm_visible = True,
                    side_chain_visible = True,
                    trace = True,
                    trace_radius = 0.25,
                    ball_stick = True,
                    ball_stick_radius = radius)

            gl_tls_group.glo_update_properties(time = 0.25)

        ## got target chain?
        console.debug_stdoutln(">html.py->Raster3D: Viewer.GLChain()") ## DEBUG
        if self.tlsmd_analysis.target_chain is not None:
            gl_chain = Viewer.GLChain(chain = self.tlsmd_analysis.target_chain)
            gl_chain.properties.update(
                    oatm_visible       = False,
                    side_chain_visible = False,
                    hetatm_visible     = True,
                    lines              = False,
                    ball_stick         = False,
                    trace              = True,
                    trace_radius       = 0.25,
                    trace_color        = "0.40,0.40,0.40")
            gl_struct.glo_add_child(gl_chain)

        driver.glr_set_render_png_path(png_file)
        viewer.glv_render_one(driver)

        ## got target chain?
        if self.tlsmd_analysis.target_chain != None:
            for atm in self.tlsmd_analysis.target_chain.iter_all_atoms():
                atm.position = atm.orig_position
                del atm.orig_position

        console.stdoutln("Raster3D: Finish rendering %s..." % png_file)
        return "", png_file

    def write_tls_pdb_file(self, chain, cpartition):
        """Write out a PDB file with the TLS predicted anisotropic ADPs for
        this segmentation.
        """
        ## class HTMLReport()
        basename = "%s_CHAIN%s_NTLS%d_UTLS"  % (
                   self.struct_id, chain.chain_id,
                   cpartition.num_tls_segments())
        pdb_file = "%s.pdb" % (basename)

        ## temporarily set the atom temp_factor and U tensor to the Utls value
        old_temp_factor = {}
        old_U = {}
        for tls in cpartition.iter_tls_segments():
            for atm, Utls in tls.tls_group.iter_atm_Utls():
                old_temp_factor[atm] = atm.temp_factor
                old_U[atm] = atm.U

                atm.temp_factor = Constants.U2B * (numpy.trace(Utls)/3.0)
                atm.U = Utls

        FileIO.SaveStructure(fil = pdb_file, struct = self.struct)

        console.stdoutln("[%s,%s] PDB: Saving %s" % (
            chain.chain_id, cpartition.num_tls_segments(), pdb_file))

        ## restore atom temp_factor and U
        for atm, temp_factor in old_temp_factor.iteritems():
            atm.temp_factor = temp_factor
            atm.U = old_U[atm]

    def write_tlsout_file(self, chain, cpartition):
        """Writes the TLSOUT file for the segmentation.
        """
        ## class HTMLReport()
        basename = "%s_CHAIN%s_NTLS%d" % (
                   self.struct_id, chain.chain_id,
                   cpartition.num_tls_segments())
        tlsout_file = "%s.tlsout" % (basename)

        struct_id = self.struct_id
        chain_id  = chain.chain_id

        tls_file = TLS.TLSFile()
        tls_file.set_file_format(TLS.TLSFileFormatTLSOUT())

        ##======================================================================
        ##<FLATFILE>
        chain_ntls = "%s,%s" % (chain_id, cpartition.num_tls_segments())
        if os.path.basename(os.getcwd()) != "ANALYSIS":
            flatfile = open("../%s.dat" % conf.globalconf.job_id, "a+")
        else:
            flatfile = open("%s.dat" % conf.globalconf.job_id, "a+")
        flatfile.write("\nTIME %s.0 [%s] write_tlsout_file" % (
            chain_ntls, misc.timestamp()))

        for tls in cpartition.iter_tls_segments():
            ## don't write out bypass edges
            if tls.method != "TLS":
                continue

            tls_desc = TLS.TLSGroupDesc()
            tls_file.tls_desc_list.append(tls_desc)
            tls_desc.set_tls_group(tls.tls_group)
            for frag_id1, frag_id2 in tls.iter_segment_ranges():
                tls_desc.add_range(chain_id, frag_id1, chain_id, frag_id2, "ALL")
                flatfile.write("\nCCCC %s.0 TLSOUT %s:%s %f" % (
                    chain_ntls, frag_id1, frag_id2, tls.tls_group.itls_T))

        flatfile.close()
        tls_file.save(open(tlsout_file, "w"))

        console.stdoutln("[%s,%s] REFMAC: Saving %s" % (
            chain_id, cpartition.num_tls_segments(), tlsout_file))

        return tlsout_file

    def write_phenixout_file(self, chain, cpartition):
        """Writes the TLS-PHENIX-OUT file for the segmentation.
        """
        ## class HTMLReport()
        basename = "%s_CHAIN%s_NTLS%d" % (
                   self.struct_id, chain.chain_id,
                   cpartition.num_tls_segments())
        phenixout_file = "%s.phenixout" % (basename)

        struct_id = self.struct_id
        chain_id  = chain.chain_id

        phenix_file = TLS.TLSFile()
        phenix_file.set_file_format(TLS.TLSFileFormatPHENIXOUT())

        for tls in cpartition.iter_tls_segments():
            ## don't write out bypass edges
            if tls.method != "TLS":
                continue

            tls_desc = TLS.TLSGroupDesc()
            phenix_file.tls_desc_list.append(tls_desc)
            tls_desc.set_tls_group(tls.tls_group)
            for frag_id1, frag_id2 in tls.iter_segment_ranges():
                tls_desc.add_range(chain_id, frag_id1, chain_id, frag_id2, "ALL")

        phenix_file.save(open(phenixout_file, "w"))

        console.stdoutln("[%s,%s] PHENIX: Saving %s" % (
            chain_id, cpartition.num_tls_segments(), phenixout_file))

        return phenixout_file

    def chain_ntls_analysis(self, chain, cpartition):
        """Generate ntls optimization constraint report and free memory.
        """
        ## class HTMLReport()
        report = ChainNTLSAnalysisReport(chain, cpartition)
        return report

    def jmol_html(self, chain, cpartition):
        """Writes out the HTML page which will display the
        structure using the Jmol Applet.
        """
        ## class HTMLReport()
        jmol_file = "%s_CHAIN%s_NTLS%d_JMOL.html"  % (
                    self.struct_id, chain.chain_id,
                    cpartition.num_tls_segments())

        ## SEE:
        ##    - http://wiki.jmol.org:81/index.php/AtomSets
        ##    - http://chemapps.stolaf.edu/jmol/docs/index.htm

        ## create the Jmol script using cartoons and consistant
        ## coloring to represent the TLS groups
        js = ['load %s;' % (self.struct_path),
              'select *;',
              'cpk off;',
              'wireframe off;',
              'select protein;',
              'cartoon on;']

        ## loop over TLS groups and color
        for tls in cpartition.iter_tls_segments():
            js.append('select %s;' % (tls.jmol_select()))
            js.append('color [%d,%d,%d];' % (tls.color.rgbi))

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
             '<title>Chain %s using %d TLS Groups</title>' % (
             chain.chain_id, cpartition.num_tls_segments()),
             '<script type="text/javascript" src="Jmol.js">',
             '</script>',
             '</head>',
             '<body>',
             '<script type="text/javascript">',
             'jmolInitialize("./");',
             'jmolSetAppletColor("white");',
             'jmolApplet(%d, "%s");' % (conf.JMOL_SIZE, "".join(js)),
             '</script>',
             '</body>',
             '</html>']

        open(jmol_file, "w").write("".join(l))

        console.stdoutln("[%s,%s] JMOL: Saving %s" % (
            chain.chain_id, cpartition.num_tls_segments(), jmol_file))

        return jmol_file

    def iter_filter_atoms(self, atom_iter):
        """Filters out any non-backbone atoms.
        """
        ## class HTMLReport()
        filter = lambda atm: conf.DISPLACE_ATOM_NAME_DICT.has_key(atm.name)
        return itertools.ifilter(filter, atom_iter)

    def generate_raw_backbone_file(self, chain, cpartition):
        """Writes out the 'raw' backbone data for the tlsanim2r3d program.
        NOTE: This should only be run if JMOL_SKIP == True.
        """
        ## class HTMLReport()
        basename = "%s_CHAIN%s_NTLS%d" % (
            self.struct_id, chain.chain_id, cpartition.num_tls_segments())
        raw_r3d_file  = "%s.raw" % (basename)
        r3d_body_file = "%s.r3d" % (basename)

        this_seg = ""
        which_ntls = 0
        raw_r3d_filename = open(raw_r3d_file, "w")
        for tls in cpartition.iter_tls_segments():
            if str(this_seg) != str(tls):
                which_ntls += 1
            this_seg = str(tls)

            for frag_id1, frag_id2 in tls.iter_segment_ranges():
                for frag in Structure.iter_fragments(chain.iter_fragments(), frag_id1, frag_id2):
                    for atm in self.iter_filter_atoms(frag.iter_atoms()):
                        ## Raw input file for tlsanim2r3d->Raster3D.
                        ## Only save backbones atoms.
                        ## on/off model chain segment libration x y z
                        ## E.g., "1 0 A 0 0 7.069 -24.991 -2.991"
                        raw_r3d_filename.write("1 1 %s %s 1 %.3f %.3f %.3f\n" % (
                            chain.chain_id, which_ntls,
                            atm.position[0], atm.position[1], atm.position[2]))
        raw_r3d_filename.close()

        return raw_r3d_file, r3d_body_file

    def jmol_animate_html(self, chain, cpartition):
        """Writes out the HTML page which will display the
        structure using the Jmol Applet.
        """
        ## class HTMLReport()
        basename = "%s_CHAIN%s_NTLS%d_ANIMATE" % (
            self.struct_id, chain.chain_id, cpartition.num_tls_segments())

        ## TODO: Include http://www.wwpdb.org/documentation/format32/sect8.html
        html_file     = "%s.html" % (basename)
        pdb_file      = "%s.pdb" % (basename)
        raw_r3d_file  = "%s.raw" % (basename)
        r3d_body_file = "%s.r3d" % (basename)

        ## generate animation PDB file
        try:
            console.stdoutln("[%s,%s] TLSAnimate: creating animation PDB file..." % (
                chain.chain_id, cpartition.num_tls_segments()))
            tlsa = TLSAnimate(chain, cpartition)
            tlsa.construct_animation(pdb_file, raw_r3d_file)
        except TLSAnimateFailure:
            console.stdoutln("     Warning: failed to create animation PDB file.")
            print console.formatExceptionInfo()
            pass

        ## create the Jmol script using cartoons and consistant
        ## coloring to represent the TLS groups
        js = ['load %s;' % (pdb_file),
              'select *;',
              'cpk off;',
              'wireframe off;',
              'select protein;',
              'trace on;']

        ## figure out which libration eigen value is the largest and
        ## use that value in the animation
        n = 0 ## counter for which max_libration to use
        max_libration = []
        for tls in cpartition.iter_tls_segments():
            tls_info = tls.model_tls_info
            L1_val = float(tls_info["L1_eigen_val"]) * Constants.RAD2DEG2
            L2_val = float(tls_info["L2_eigen_val"]) * Constants.RAD2DEG2
            L3_val = float(tls_info["L3_eigen_val"]) * Constants.RAD2DEG2
            max = 0.00
            for val in L1_val, L2_val, L3_val:
                if val >= max:
                    ## TODO: Store max libration chain_id to ANIMATE.txt, 2009-08-05
                    max = val
            max_libration.append(float(max))
            n += 1

        use_chain = ""
        n = 0
        ## loop over TLS groups and color _only_ the max_libration
        for tls in cpartition.iter_tls_segments():
            tls_info = tls.model_tls_info
            L1_val = float(tls_info["L1_eigen_val"]) * Constants.RAD2DEG2
            L2_val = float(tls_info["L2_eigen_val"]) * Constants.RAD2DEG2
            L3_val = float(tls_info["L3_eigen_val"]) * Constants.RAD2DEG2

            if L1_val == max_libration[n]:
                use_chain = tlsa.L1_chain.chain_id
            elif L2_val == max_libration[n]:
                use_chain = tlsa.L2_chain.chain_id
            elif L3_val == max_libration[n]:
                use_chain = tlsa.L3_chain.chain_id
            n += 1

            l = []
            for frag_id1, frag_id2 in tls.iter_segment_ranges():
                l.append("%s-%s:%s" % (frag_id1, frag_id2, use_chain))
            select = ",".join(l)
            use_chain = ""

            js.append('select %s;' % select)
            js.append('color [%d,%d,%d];' % (tls.color.rgbi))

        ##=====================================================================
        ## NOTE: Old code to loop over TLS groups and color
        #for tls in cpartition.iter_tls_segments():
        #    chain_ids = [tlsa.L1_chain.chain_id,
        #                 tlsa.L2_chain.chain_id,
        #                 tlsa.L3_chain.chain_id]

        #    for chain_id in chain_ids:
        #        l = []
        #        for frag_id1, frag_id2 in tls.iter_segment_ranges():
        #            l.append("%s-%s:%s" % (frag_id1, frag_id2, chain_id))
        #        select = ",".join(l)

        #        ## switched to selecting all three libration states
        #        #js.append('select %s;' % (tls.jmol_select()))
        #        js.append('select %s;' % select)
        #        #console.stdoutln("JMOL-select %s vs. %s vs. %s" % (
        #            chain_id, tls.jmol_select(), select))

        #        js.append('color [%d,%d,%d];' % (tls.color.rgbi))
        ##=====================================================================

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
             '<title>Chain %s using %d TLS Groups</title>' % (
                 chain.chain_id, cpartition.num_tls_segments()),
             '<script type="text/javascript" src="%s/Jmol.js">' % (
                 conf.JMOL_DIR),
             '</script>',
             '</head>',
             '<body>',
             '<script type="text/javascript">',
             'jmolInitialize("%s");' % (conf.JMOL_DIR),
             'jmolSetAppletColor("white");',
             'jmolApplet(%d, "%s");' % (conf.JMOL_SIZE, "".join(js)),
             '</script>',
             '</body>',
             '</html>']

        ## manually free memory
        import gc
        gc.collect()

        open(html_file, "w").write("".join(l))

        console.stdoutln("[%s,%s] HTML: Saving %s" % (
            chain.chain_id, cpartition.num_tls_segments(), html_file))

        return html_file, raw_r3d_file, r3d_body_file

    def write_multi_chain_alignment(self):
        """Write out the chain residue alignment page.
        """
        ## class HTMLReport()
        ## only write out the comparison page if there is more than one
        ## chain analyzed in the structure
        if self.tlsmd_analysis.num_chains() < 2:
            return

        path  = "%s_CHAIN_COMP.html" % (self.struct_id)
        title = "Multi-Chain Alignment Analysis of TLS Groups"

        self.page_multi_chain_alignment = {
            "title": title,
            "href":  path }

        fil = open(path, "w")
        fil.write(self.html_multi_chain_alignment())
        fil.close()

    def html_multi_chain_alignment(self):
        """Write out all HTML/PDB/TLSIN files which compare
        chains in the structure.
        """
        ## class HTMLReport()
        title = self.page_multi_chain_alignment["title"]

        l = [ self.html_head(title),
              self.html_title(title),

              '<center>',
              '<a href="index.html">Back to Index</a>',
              '</center>',
              '<br/>' ]

        ## figure out the maximum number of ntls in all chains
        max_ntls = 0
        for chain in self.tlsmd_analysis.iter_chains():
            max_ntls = max(chain.partition_collection.max_ntls(), max_ntls)

        ## generate ntls number of plots and add them to the
        ## HTML document
        for ntls in xrange(1, max_ntls + 1):

            ## create a 2-tuple list of (chain_id, cpartiton) for
            ## each chain which a a TLSMD segmentation of h groups
            seg_list = []
            for chain in self.tlsmd_analysis.iter_chains():
                cpartition = chain.partition_collection.get_chain_partition(ntls)
                if cpartition != None:
                    seg_list.append((chain.chain_id, cpartition))

            ## generate the TLS segmentation alignment plot for all chains
            plot = sequence_plot.TLSSegmentAlignmentPlot()
            chain_id_list = []

            for chain in self.tlsmd_analysis.iter_chains():
                cpartition = chain.partition_collection.get_chain_partition(ntls)
                if cpartition != None:
                    chain_id_list.append(chain.chain_id)
                    plot.add_tls_segmentation(cpartition)

            basename  = "%s_NTLS%d"  % (self.struct_id, ntls)
            plot_file = "%s_ALIGN.png" % (basename)
            plot.plot(plot_file)

            ## write HTML
            l.append('<h3>Chains Alignment using %d TLS Groups</h3>' % (ntls))

            ## plot table
            l += ['<table border="1">',
                  '<tr><th>Chain</th><th>Chain Alignment</th></tr>',
                  '<tr>',
                  '<td align="center">',
                  '<table border="0" cellspacing="0" cellpadding="0">' ]

            for chain_id in chain_id_list:
                l.append('<tr><td align="right" valign="middle" height="20">')
                l.append('<font size="-20">%s</font></td></tr>' % (chain_id))

            l += ['</table>'
                  '</td>',
                  '<td><img src="%s" alt="Segmentation Plot" /></td>' % (plot_file),
                  '</tr>',
                  '</table>']

        l.append(self.html_foot())

        return "".join(l)

    def write_refinement_prep(self):
        """Generate form to allow users to select the number of TLS groups
        to use per chain.
        """
        ## class HTMLReport()
        path  = "%s_REFINEMENT_PREP.html" % (self.struct_id)
        title = "Generate input files for multigroup TLS Refinement"

        self.page_refinement_prep = {
            "title": title,
            "href" : path }

        fil = open(path, "w")
        fil.write(self.html_refinement_prep())
        fil.close()

        console.stdoutln("HTML: Saving %s" % path)

    def html_refinement_prep(self):
        ## class HTMLReport()
        title = self.page_refinement_prep["title"]
        plot = gnuplots.LSQR_vs_TLS_Segments_All_Chains_Plot(self.tlsmd_analysis)

        l = [self.html_head(title),
             self.html_title(title),
             '<center><h3>',
             'Step 1: Select the number of TLS groups for each chain',
             '</h3></center>',
             '<center>',
             '<a href="index.html">Back to Index</a>',
             '</center>',
             '<br/>',
             '<form enctype="multipart/form-data" action="%s" method="post">' % (
                 conf.REFINEPREP_URL),
             '<input type="hidden" name="job_id" value="%s">' % (
                 conf.globalconf.job_id),
             '<p>%s</p>' % (captions.REFINEMENT_PREP_INFO),
             '<center><table><tr><td>',
             plot.html_link(),
             '<br/>',
             '</td></tr><tr><td>',
             '<table width="100%" border="1">',
             '<tr><th>',
             '<p>Select the Number of TLS Groups per Chain</p>',
             '</th></tr>',
             '<tr><td align="center">',
             '<table cellspacing="5">' ]

        for chain in self.tlsmd_analysis.iter_chains():
            l += ['<tr><td>',
                  'Number of TLS Groups for Chain %s' % (chain.chain_id),
                  '</td><td>',
                  '<select name="NTLS_CHAIN%s">' % (chain.chain_id) ]
            for ntls in chain.partition_collection.iter_ntls():
                l.append('<option value="%d">%d</option>' % (ntls, ntls))
            l += ['</select>',
                  '</td></tr>' ]

        l += ['</table>',
              '</td></tr>',
              '<tr><td align="right">',
              '<input type="submit" value="OK">',
              '</td></tr></table>',
              '</td></tr></table></center>',
              '</form>' ]

        l.append(self.html_foot())
        return "".join(l)


class ChainNTLSAnalysisReport(Report):
    """Writes a HTML report detailing one given TLS segmentation of a chain.
    """
    def __init__(self, chain, cpartition):
        Report.__init__(self)

        self.chain = chain
        self.cpartition = cpartition

        self.struct = chain.struct
        self.struct_id = chain.struct.structure_id
        self.chain_id = chain.chain_id
        self.ntls = cpartition.num_tls_segments()

        self.root  = ".."
        self.dir   = "%s_CHAIN%s_NTLS%d"  % (
            self.struct_id, self.chain_id, self.ntls)
        self.index = "%s.html" % (self.dir)
        self.url   = "%s/%s" % (self.dir, self.index)

        self.write_report()

    def write_report(self):
        if not os.path.isdir(self.dir):
            os.mkdir(self.dir)
        os.chdir(self.dir)

        try:
            console.stdoutln("HTML: Writing %s" % self.index)
            self.write_all_files()
        except:
            print console.formatExceptionInfo()
        finally:
            os.chdir(self.root)

    def write_all_files(self):
        """Writes analysis details of each TLS group.
        """
        title = "Chain %s Partitioned by %d TLS Groups" % (
            self.chain_id, self.ntls)
        path = "%s_CHAIN%s_ANALYSIS.html" % (self.struct_id, self.chain_id)

        l = [self.html_head(title),
             self.html_title(title),

             '<center>',
             '<a href="../index.html">Back to Index</a>',
             '&nbsp;&nbsp;&nbsp;&nbsp;',
             '<a href="../%s">Back to Chain %s Analysis</a>' % (
                 path, self.chain_id),
             '</center>',

             '<br/>']

        ## NOTE: The following table is for the "Full TLS Partition Analysis"
        ## pages
        l += [self.html_tls_group_table(detail = False),'<br/><hr>']
        l += [self.html_bmean(),'<br/><hr>'] ## REQUIRED!
        l += [self.tls_segment_recombination(),'<br/><hr>']
        l += [self.html_translation_analysis(),'<br/><hr>']
        l += [self.html_libration_analysis(),'<br/><hr>']
        l += [self.html_ca_differance(),'<br/><hr>']
        l += [self.html_rmsd_plot()] ## REQUIRED!

        ## EAM April 2008 - This sure looks redundant to me.
        ## Let's try doing without it
        #self.tls_segment_recombination()

        ## Create histogram plots
        job_id = conf.globalconf.job_id
        if conf.globalconf.generate_histogram == False:
            console.stdoutln("NOTE: Skipping Histogram section")
        else:
            try:
                for tls in self.cpartition.iter_tls_segments():
                    ## don't write out bypass edges
                    if tls.method != "TLS":
                        continue
                    l.append('<hr>')
                    l.append(self.html_tls_fit_histogram(tls))
            except:
                console.stdoutln("     Warning: failure in html_tls_fit_histogram()")
                print console.formatExceptionInfo()
                pass

        l.append(self.html_foot())

        open(self.index, "w").write("".join(l))
        console.stdoutln("HTML: Saving %s" % path)

    def html_tls_group_table(self, detail):
        ## Pass "ntls" as well, 2008-04-05
        return html_tls_group_table(self.ntls, self.chain, self.cpartition, "..", detail)

    def html_translation_analysis(self):
        """Perform a translation analysis of the protein chain as
        spanned by the tlsopt TLS groups.
        """
        tanalysis = gnuplots.TranslationAnalysis(self.chain, self.cpartition)
        l = ['<center>',
             tanalysis.html_markup("Translation Analysis of T<sup>r</sup>",
             captions.TRANSLATION_GRAPH_CAPTION),
             '</center>']
        return "".join(l)

    def html_libration_analysis(self):
        """Perform a libration analysis of the protein chain as
        spanned by the tlsopt TLS groups.
        """
        libration_analysis = gnuplots.LibrationAnalysis(self.chain, self.cpartition)
        l = ['<center>',
             libration_analysis.html_markup("Screw Displacement Analysis",
             captions.LIBRATION_GRAPH_CAPTION),
             '</center>']
        return "".join(l)

    def html_ca_differance(self):
        """Perform a fit analysis of the protein chain as
        spanned by the tlsopt TLS groups.
        """
        plot = gnuplots.CA_TLS_Differance_Plot(self.chain, self.cpartition)
        l = ['<center>',
             plot.html_markup("Deviation of Observed CA Atom B-Factors From TLS Model",
             captions.FIT_GRAPH_CAPTION),
             '</center>']
        return "".join(l)

    def html_bmean(self):
        """Mean B-Factor per residue.
        """
        self.bmean_plot = gnuplots.BMeanPlot(self.chain, self.cpartition)
        l = ['<center>',
             self.bmean_plot.html_markup("Mean BFactor Analysis",
                 "Comparison of TLS predicted B factors with experimental (input) B factors."),
             '</center>']
        return "".join(l)

    def html_rmsd_plot(self):
        rmsd_plot = gnuplots.RMSDPlot(self.chain, self.cpartition)
        l = ['<center>',
             rmsd_plot.html_markup("RMSD Deviation of Observed vs. TLS Predicted B Factors",
                                   captions.RMSD_DEVIATION_GRAPH_CAPTION),
             '</center>']
        return "".join(l)

    def html_tls_fit_histogram(self, tls):
        """histogram of atomic U_ISO - U_TLS_ISO
        """
        his = gnuplots.UIso_vs_UtlsIso_Histogram(self.chain, self.cpartition, tls)
        title = 'Distribution Histogram of TLS Group %s' % (
            tls.display_label())
        l = ['<center>',
             his.html_markup(title, ""),
             '</center>']
        return "".join(l)

    def tls_segment_recombination(self):
        """The TLSMD optimization algorithm models TLS groups as sequential 
        segments of a protein or DNA/RNA chain. This matrix shows the 
        RMSD B values of the individual groups on the diagonal, and the 
        RMSD B values of combined groups as off-diagonal elements. This 
        helps identify non-contiguous protein segments which may be 
        combined into a single TLS group.
        """
        if not hasattr(self.cpartition, "rmsd_b_mtx"):
            return ""

        rmatrix = self.cpartition.rmsd_b_mtx
        m, n = rmatrix.shape
        chain_ntls = "%s,%s" % (
            self.chain_id, self.cpartition.num_tls_segments())

        tbl = table.StringTableFromMatrix(rmatrix)
        ## E.g., for three segments, the file would look something like the
        ## following:
        ## 7.11514881827     9.86985699559     10.1031000732
        ## 9.86985699559     9.94840587829     8.74779028204
        ## 10.1031000732     8.74779028204     9.61743006191
        filename = "%s_CHAIN%s_NTLS%s_RECOMBINATION.txt" % (
            self.struct_id, self.chain_id, self.cpartition.num_tls_segments())

        open(filename, "w").write(str(tbl))
        console.stdoutln("[%s] RECOMBINATION: Saving %s" % (
            chain_ntls, filename))

        ##======================================================================
        ##<FLATFILE>
        if os.path.basename(os.getcwd()) != "ANALYSIS":
            flatfile = open("../%s.dat" % conf.globalconf.job_id, "a+")
        else:
            flatfile = open("%s.dat" % conf.globalconf.job_id, "a+")

        flatfile.write("\nTIME %s.0 [%s] tls_segment_recombination" % (
            chain_ntls, misc.timestamp()))
        flatfile.write("\nCCCC %s.0 RECOMBINATION" % chain_ntls)
        ##</FLATFILE>
        ##======================================================================

        tbl = table.StringTable(m + 1, n + 1, '<td></td>')

        for i, tls in enumerate(self.cpartition.iter_tls_segments()):
            ##self.cpartition = "(B:1-10)(B:11-15)(B:16-21)(B:22-27)(B:28-39)"

            ## label top and left table margins
            tbl[i + 1, 0] = '<td style="background-color:%s">%s</td>' % (
                tls.color.rgbs, tls.display_label())
            tbl[0, i + 1] = tbl[i + 1, 0]

            ##<FLATFILE>
            flatfile.write("\nRECB %s.%s RMSD_B: %.2f" % (
                chain_ntls, i+1, tls.rmsd_b))
            flatfile.write("\nRECB %s.%s RESIDUAL_RMSD_B: %.2f" % (
                chain_ntls, i+1, rmatrix[i, i]))
            flatfile.write("\nRECB %s.%s T = %.8f" % (
                chain_ntls, i+1, tls.tls_group.itls_T))
            ## L = | 0,0 0,1 0,2 |
            ##     | 1,0 1,1 1,2 |
            ##     | 2,0 2,1 2,2 |
            ## L = 0,0 1,1 2,2 1,0 2,1 2,0
            ## L = [0][0] [1][1] [2][2] [1][0] [2][1] [2][0]
            flatfile.write("\nRECB %s.%s L = %.8f,%.8f,%.8f,%.8f,%.8f,%.8f" % (
                chain_ntls, i+1,
                tls.tls_group.itls_L[0][0],
                tls.tls_group.itls_L[1][1],
                tls.tls_group.itls_L[2][2],
                tls.tls_group.itls_L[1][0],
                tls.tls_group.itls_L[2][1],
                tls.tls_group.itls_L[2][0]))
            flatfile.write("\nRECB %s.%s S = %.8f,%.8f,%.8f" % (
                chain_ntls, i+1,
                tls.tls_group.itls_S[0],
                tls.tls_group.itls_S[1],
                tls.tls_group.itls_S[2]))
            flatfile.write("\nRECB %s.%s ORIGIN = %.3f,%.3f,%.3f" % (
                chain_ntls, i+1,
                tls.tls_group.origin[0],
                tls.tls_group.origin[1],
                tls.tls_group.origin[2]))
            ##</FLATFILE>

        ## recombination values
        console.debug_stdoutln(">MATRIX: Creating RMSD B values table...")

        ## determine min/max values for coloring
        min_rmsd_b = rmatrix[0, 0]
        max_rmsd_b = rmatrix[0, 0]
        for i in xrange(m):
            flatfile.write('\nMATX %s.%s ' % (chain_ntls, i+1)) ## FLATFILE
            for j in xrange(n):
                flatfile.write('%.2f,' % rmatrix[i, j]) ## FLATFILE
                min_rmsd_b = min(min_rmsd_b, rmatrix[i, j])
                max_rmsd_b = max(max_rmsd_b, rmatrix[i, j])
        b_range = max_rmsd_b - min_rmsd_b

        flatfile.write("\n")
        flatfile.close()

        for i in xrange(m):
            for j in xrange(n):
                scale = (rmatrix[i, j] - min_rmsd_b) / b_range
                bright = 1.0 - (0.75 * scale)
                rgbs = misc.rgb_f2s((bright, bright, bright))
                tbl[i + 1, j + 1] = '<td class="matrix" style="background-color:%s;">%6.2f</td>' % (
                    rgbs, rmatrix[i, j])

        l = ['<table class="matrix">\n']
        for i in xrange(m + 1):
            l.append('<tr>')
            for j in xrange(n + 1):
                l.append(tbl[i, j])
            l.append('</tr>\n')
        l += ['</table>\n']

        l2 = ['<center>',
              gnuplots.FormatFigureHTML("RMSD B Values of Combined TLS Groups",
                                        captions.TLS_GROUP_RECOMBINATION,
                                        "".join(l)),
              '</center>']

        return "".join(l2)
