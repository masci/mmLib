## TLS Minimized Domains (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

###############################################################################
## Report Directory Generation
##

import popen2

## Python Imaging Library imports
import Image
import ImageDraw

## mmLib
from mmLib.Colors         import *
from mmLib.Viewer         import *
from mmLib.R3DDriver      import Raster3DDriver

## tlsmdlib
from misc                 import *
from captions             import *
from tls_animate          import TLSAnimate, TLSAnimateFailure

## program paths
GNUPLOT_PATH = "gnuplot"
GNUPLOT_FONT = "/home/jpaint/tlsmd/fonts/LucidaSansOblique.ttf"
GNUPLOT_FONT_SIZE = "10"

JMOL_DIR     = "../../../jmol"

## constants

## the pixel width of the TLS visualization rendered ray traces
VIS_WIDTH = 800

## pixel size of the gnuplot generated images
GNUPLOT_WIDTH = 600

## target pixel width, height, and spacing of sequence
## alignment plots
ALIGN_TARGET_WIDTH = 500
ALIGN_HEIGHT       = 15
ALIGN_SPACING      = 5

## the JMol viewer is a square window, generated with
## this pixel size
JMOL_SIZE = 600

## the isoprobability contour level for all
## visualizations
ADP_PROB = 85


def rgb_f2i(rgb):
    """Transforms the float 0.0-1.0 RGB color values to
    integer 0-255 RGB values.
    """
    r, g, b = rgb
    ri = int(255.0 * r)
    gi = int(255.0 * g)
    bi = int(255.0 * b)
    return (ri, gi, bi)


def calc_inertia_tensor(atom_iter):
    """Calculate moment of inertia tensor at the centroid
    of the atoms.
    """
    al              = AtomList(atom_iter)
    centroid        = al.calc_centroid()

    I = zeros((3,3), Float)
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

    evals, evecs = eigenvectors(I)

    elist = [(evals[0], evecs[0]),
             (evals[1], evecs[1]),
             (evals[2], evecs[2])]

    elist.sort()

    R = array((elist[0][1], elist[1][1], elist[2][1]), Float)

    ## make sure the tensor uses a right-handed coordinate system
    if allclose(determinant(R), -1.0):
        I = identity(3, Float)
        I[0,0] = -1.0
        R = matrixmultiply(I, R)
    assert allclose(determinant(R), 1.0)

    return centroid, R


def calc_orientation(struct, chain):
    """Orient the structure based on a moment-of-intertia like tensor
    centered at the centroid of the structure.
    """
    ori = {}

    def iter_protein_atoms(sobjx):
        for fragx in sobjx.iter_amino_acids():
            for atmx in fragx.iter_atoms():
                yield atmx
                
    str_centroid, str_R = calc_inertia_tensor(iter_protein_atoms(struct))
    chn_centroid, chn_R = calc_inertia_tensor(iter_protein_atoms(chain))

    ## now calculate a rectangular box
    min_x = 0.0
    max_x = 0.0
    min_y = 0.0
    max_y = 0.0
    min_z = 0.0
    max_z = 0.0

    for atm in chain.iter_all_atoms():
        x     = matrixmultiply(str_R, atm.position - chn_centroid)
        min_x = min(min_x, x[0])
        max_x = max(max_x, x[0])
        min_y = min(min_y, x[1])
        max_y = max(max_y, x[1])
        min_z = min(min_z, x[2])
        max_z = max(max_z, x[2])

    ## add slop splace around the edges
    slop   = 2.0

    min_x -= slop
    max_x += slop
    min_y -= slop
    max_y += slop
    min_z -= slop
    max_z += slop

    ## calculate the zoom based on a target width
    target_pwidth = VIS_WIDTH

    hwidth  = max(abs(min_x),abs(max_x))
    hheight = max(abs(min_y),abs(max_y))
    pheight = target_pwidth * (hheight / hwidth)
    hzoom   = 2.0 * hwidth

    ori["R"]        = str_R
    ori["centroid"] = chn_centroid
    ori["pwidth"]   = target_pwidth
    ori["pheight"]  = pheight 
    ori["hzoom"]    = hzoom

    ## calculate near, far clipping blane
    ori["near"] = max_z
    ori["far"]  = min_z
    
    return ori


class FragmentID(object):
    """A fragment ID class acts a lot like a string, but separates the
    res_seq and icode internally.
    """
    def __init__(self, frag_id):
        self.res_seq = 1
        self.icode = ""
        try:
            self.res_seq = int(frag_id)
        except ValueError:
            try:
                self.res_seq = int(frag_id[:-1])
            except ValueError:
                pass
            else:
                self.icode = frag_id[-1]
    def __str__(self):
        return str(self.res_seq) + self.icode
    def __lt__(self, other):
        return (self.res_seq, self.icode) < (other.res_seq, other.icode)
    def __le__(self, other):
        return (self.res_seq, self.icode) <= (other.res_seq, other.icode)
    def __eq__(self, other):
        return (self.res_seq, self.icode) == (other.res_seq, other.icode)
    def __ne__(self, other):
        return (self.res_seq, self.icode) != (other.res_seq, other.icode)
    def __gt__(self, other):
        return (self.res_seq, self.icode) > (other.res_seq, other.icode)
    def __ge__(self, other):
        return (self.res_seq, self.icode) >= (other.res_seq, other.icode)



class GNUPlot(object):
    """Provides useful methods for subclasses which need to run gnuplot.
    """
    def gnuplot_mkscript(self, template, replace_dict):
        """Replaces key strings in the template with the values in
        the replac_dict.  Returns the modificed template.
        """
        for key, val in replace_dict.items():
            template = template.replace(key, val)
        return template
    
    def gnuplot_run(self, script, basename=None):
        """Runs gnuplot
        """
        ## if a basename is given, then write the GnuPlot script
        ## as a file
        if basename!=None:
            open("%s.plot" % (basename), "w").write(script)
        
        ## run gnuplot
        stdout, stdin, stderr = popen2.popen3((GNUPLOT_PATH, ), 32768)
        stdin.write(script)
        stdin.close()

        #stdout.read()
        stdout.close()
        #stderr.read()
        stderr.close()
        

_LSQR_VS_TLS_SEGMENTS_TEMPLATE = """\
set xlabel "Number of TLS Segments"
set xrange [1:20]
set ylabel "Residual"
set format y "%5.2f"
set style line 1 lw 3
set term png enhanced font "<font>" <fontsize>
set output "<pngfile>"
set title "<title>"
plot "<txtfile>" using 1:2 title "Minimization (Weighted) Residual" ls 1 with linespoints
"""
        
class LSQR_vs_TLS_Segments_Plot(GNUPlot):
    def __init__(self, chainopt):
        ## generate data and png paths
        basename = "%s_CHAIN%s_RESID" % (
            chainopt["struct_id"] , chainopt["chain_id"])

        self.txt_path = "%s.txt" % (basename)
        self.png_path = "%s.png" % (basename)

        fil = open(self.txt_path, "w")
        for h, tlsopt in chainopt["ntls_list"]:
            fil.write("%10d %f\n" % (h, tlsopt.residual))
        fil.close()

        ## modify script template
        script = _LSQR_VS_TLS_SEGMENTS_TEMPLATE
        script = script.replace("<font>", GNUPLOT_FONT)
        script = script.replace("<fontsize>", GNUPLOT_FONT_SIZE)
        script = script.replace("<txtfile>", self.txt_path)
        script = script.replace("<pngfile>", self.png_path)
        script = script.replace(
            "<title>", "Least Squares Residual vs. Number of TLS "\
            "Segments for %s Chain %s " % (
            chainopt["struct_id"], chainopt["chain_id"]))

        self.gnuplot_run(script, basename)



_LSQR_VS_TLS_SEGMENTS_ALL_CHAINS_TEMPLATE = """\
set xlabel "Number of TLS Segments"
set xrange [1:20]
set ylabel "Minimization (Weighted) LSQR Residual"
set format y "%5.2f"
set term png enhanced font "<font>" <fontsize>
set output "<pngfile>"
set title "<title>"
"""

class LSQR_vs_TLS_Segments_All_Chains_Plot(GNUPlot):
    def __init__(self, chainopt_list):
        struct_id = chainopt_list[0]["struct_id"]
        
        ## generate data and png paths
        basename = "%s_RESID" % (struct_id)
        self.png_path = "%s.png" % (basename)

        ## prepare gnuplot script
        script = _LSQR_VS_TLS_SEGMENTS_ALL_CHAINS_TEMPLATE
        script = script.replace("<font>", GNUPLOT_FONT)
        script = script.replace("<fontsize>", GNUPLOT_FONT_SIZE)
        script = script.replace("<pngfile>", self.png_path)
        script = script.replace(
            "<title>", "Least Squares Residual vs. Number of TLS "\
            "Segments of %s" % (struct_id))

        ## re-use the data files of LSQRvsNTLS from the individual
        ## graphs; to do this the filenames have to be re-constructed
        plist = []
        for chainopt in chainopt_list:
            chain_id = chainopt["chain_id"]
            filename = "%s_CHAIN%s_RESID.txt" % (struct_id, chain_id)
            x = '"%s" using 1:2 title "Chain %s" lw 3 with linespoints' % (
                filename, chain_id)
            plist.append(x)
        script += "plot " + string.join(plist, ",\\\n\t") + "\n"

        self.gnuplot_run(script, basename)


_TRANSLATION_ANALYSIS_TEMPLATE = """\
set xlabel "Residue"
set xrange [<xrng1>:<xrng2>]
set ylabel "Angstroms Displacement"
set format y "%4.2f"
set term png enhanced font "<font>" <fontsize>
set output "<pngfile>"
set title "<title>"
"""

class TranslationAnalysis(GNUPlot):
    def __init__(self, chainopt, tlsopt):
        basename = "%s_CHAIN%s_NTLS%s_TRANSLATION" % (
            chainopt["struct_id"], chainopt["chain_id"], tlsopt.ntls)

        self.png_path = "%s.png" % (basename)

        data_file_list = []
        for tls in tlsopt.tls_list:
            filename = self.write_data_file(chainopt, tls)
            data_file_list.append(filename)

        script = _TRANSLATION_ANALYSIS_TEMPLATE
        script = script.replace("<font>", GNUPLOT_FONT)
        script = script.replace("<fontsize>", GNUPLOT_FONT_SIZE)
        script = script.replace("<fontsize>", GNUPLOT_FONT_SIZE)

        script = script.replace("<xrng1>", tlsopt.tls_list[0]["frag_id1"])
        script = script.replace("<xrng2>", tlsopt.tls_list[-1]["frag_id2"])
        
        script = script.replace("<pngfile>", self.png_path)
        script = script.replace(
            "<title>",
            "Translation Displacement Analysis of Atoms for "\
            "%d TLS Groups" % (tlsopt.ntls))

        ## line style
        ls = 0
        for tls in tlsopt.tls_list:
            ls += 1
            script += 'set style line %d lc rgb "%s" lw 3\n' % (
                ls, tls["color"]["rgbs"])

        ## plot list
        plist = []
        ls = 0
        for filename in data_file_list:
            ls += 1
            for n in (2,3,4):
                x = '"%s" using 1:%d notitle ls %d with lines' % (
                    filename, n, ls)
                plist.append(x)

        script += "plot " + string.join(plist, ",\\\n\t") + "\n"
           
        self.gnuplot_run(script, basename)

    def write_data_file(self, chainopt, tls):
        """Generate the data file and return the filename.
        """
        tls_group = tls["tls_group"]
        tls_info  = tls["tls_info"]

        ## generate a sorted list of fragment IDs from the TLS group atoms
        fid_list = []
        for atm in tls_group:
            fid = FragmentID(atm.fragment_id)
            if fid not in fid_list:
                fid_list.append(fid)
        fid_list.sort()

        ## determine Tr translational eigenvalues
        evals = eigenvalues(tls_info["rT'"])
        t1    = GAUSS3C[ADP_PROB] * math.sqrt(evals[0])
        t2    = GAUSS3C[ADP_PROB] * math.sqrt(evals[1])
        t3    = GAUSS3C[ADP_PROB] * math.sqrt(evals[2])
        
        ## write data file
        filename  = "%s_CHAIN%s_TLS%s_%s_TRANSLATION.txt" % (
            chainopt["struct_id"], chainopt["chain_id"],
            tls["frag_id1"], tls["frag_id2"])

        fil = open(filename, "w")
        for fid in fid_list:
            fil.write("%s %f %f %f\n" % (fid, t1, t2, t3))
        fil.close()

        return filename



_LIBRATION_ANALYSIS_TEMPLATE = """\
set xlabel "Residue"
set xrange [<xrng1>:<xrng2>]
set ylabel "Angstroms Displacement"
set format y "%4.2f"
set term png enhanced font "<font>" <fontsize>
set output "<pngfile>"
set title "<title>"
"""

class LibrationAnalysis(GNUPlot):
    def __init__(self, chainopt, tlsopt):      
        
        basename = "%s_CHAIN%s_NTLS%s_LIBRATION" % (
            chainopt["struct_id"], chainopt["chain_id"], tlsopt.ntls)

        self.png_path = "%s.png" % (basename)

        data_file_list = []
        for tls in tlsopt.tls_list:
            filename = self.write_data_file(chainopt, tls)
            data_file_list.append(filename)

        script = _LIBRATION_ANALYSIS_TEMPLATE
        script = script.replace("<font>", GNUPLOT_FONT)
        script = script.replace("<fontsize>", GNUPLOT_FONT_SIZE)

        script = script.replace("<xrng1>", tlsopt.tls_list[0]["frag_id1"])
        script = script.replace("<xrng2>", tlsopt.tls_list[-1]["frag_id2"])
        
        script = script.replace("<pngfile>", self.png_path)
        script = script.replace(
            "<title>",
            "Screw Displacment Analysis of backbone Atoms using "\
            "%d TLS Groups" % (tlsopt.ntls))

        ## line style
        ls = 0
        for tls in tlsopt.tls_list:
            ls += 1
            script += 'set style line %d lc rgb "%s" lw 3\n' % (
                ls, tls["color"]["rgbs"])

        ## plot list
        plist = []
        ls = 0
        for filename in data_file_list:
            ls += 1
            for n in (4,5,6):
                x = '"%s" using 3:%d smooth bezier '\
                    'notitle ls %d with lines' % (
                    filename, n, ls)
                plist.append(x)

        script += "plot " + string.join(plist, ",\\\n\t") + "\n"
           
        self.gnuplot_run(script, basename)

    def write_data_file(self, chainopt, tls):
        """Generate the data file and return the filename.
        """
        tls_group = tls["tls_group"]
        tls_info  = tls["tls_info"]
        cor       = tls_info["COR"]

        frag_dict = {}
        for atm in tls_group:
            if atm.name in ["N","CA","C"]:
                frag_dict[atm] = [FragmentID(atm.fragment_id), 0.0, 0.0, 0.0] 
                
        for n, Lx_val, Lx_vec, Lx_rho, Lx_pitch in [
            (1, "L1_eigen_val", "L1_eigen_vec", "L1_rho", "L1_pitch"),
            (2, "L2_eigen_val", "L2_eigen_vec", "L2_rho", "L2_pitch"),
            (3, "L3_eigen_val", "L3_eigen_vec", "L3_rho", "L3_pitch") ]:

            Lval   = tls_info[Lx_val]
            Lvec   = tls_info[Lx_vec]
            Lrho   = tls_info[Lx_rho]
            Lpitch = tls_info[Lx_pitch]

            for atm, frag_rec in frag_dict.items():
                d = calc_LS_displacement(
                    cor, Lval, Lvec, Lrho, Lpitch, atm.position, ADP_PROB)
                frag_rec[n] = length(d)

        ## write data file
        filename  = "%s_CHAIN%s_TLS%s_%s_LIBRATION.txt" % (
            chainopt["struct_id"], chainopt["chain_id"],
            tls["frag_id1"], tls["frag_id2"])

        fil = open(filename, "w")

        listx = []
        for atm, frag_rec in frag_dict.items():
            if   atm.name=="N":  i = 1
            elif atm.name=="CA": i = 2
            elif atm.name=="C":  i = 3
            listx.append((frag_rec[0], i, frag_rec[1],frag_rec[2],frag_rec[3]))
        listx.sort()

        for frag_rec in listx:
            fid, i, d1, d2, d3 = frag_rec

            try:
                fidf = float(str(fid))
            except ValueError:
                continue

            if i==1: fidf -= 0.33
            if i==3: fidf += 0.33
            
            fil.write("%s %1d %f %f %f %f\n" % (fid, i, fidf, d1, d2, d3))

        fil.close()

        return filename

     

_FIT_ANALYSIS_TEMPLATE = """\
set xlabel "Residue"
set xrange [<xrng1>:<xrng2>]
set ylabel "B_{obs} - B_{calc}"
set format y "%5.2f"
set term png enhanced font "<font>" <fontsize>
set output "<pngfile>"
set title "<title>"
"""

class FitAnalysis(GNUPlot):
    def __init__(self, chainopt, tlsopt):
        
        basename = "%s_CHAIN%s_NTLS%s_FIT" % (
            chainopt["struct_id"], chainopt["chain_id"], tlsopt.ntls)

        self.png_path = "%s.png" % (basename)

        data_file_list = []
        for tls in tlsopt.tls_list:
            filename = self.write_data_file(chainopt, tls)
            data_file_list.append(filename)

        script = _FIT_ANALYSIS_TEMPLATE
        script = script.replace("<font>", GNUPLOT_FONT)
        script = script.replace("<fontsize>", GNUPLOT_FONT_SIZE)

        script = script.replace("<xrng1>", tlsopt.tls_list[0]["frag_id1"])
        script = script.replace("<xrng2>", tlsopt.tls_list[-1]["frag_id2"])
        
        script = script.replace("<pngfile>", self.png_path)
        script = script.replace(
            "<title>",
            "TLS Model Fit Analysis of Backbone Atoms for "\
            "%d TLS Groups" % (tlsopt.ntls))

        ## line style
        ls = 0
        for tls in tlsopt.tls_list:
            ls += 1
            script += 'set style line %d lc rgb "%s" lw 3\n' % (
                ls, tls["color"]["rgbs"])

        ## plot list
        plist = []
        ls = 0
        for filename in data_file_list:
            ls += 1
            x = '"%s" using 1:2 smooth bezier '\
                'notitle ls %d with lines' % (
                filename, ls)

            plist.append(x)

        script += "plot " + string.join(plist, ",\\\n\t") + "\n"
           
        self.gnuplot_run(script, basename)

    def write_data_file(self, chainopt, tls):
        """Generate the data file and return the filename.
        """
        tls_group = tls["tls_group"]
        tls_info  = tls["tls_info"]

        filename  = "%s_CHAIN%s_TLS%s_%s_FIT.txt" % (
            chainopt["struct_id"], chainopt["chain_id"],
            tls["frag_id1"], tls["frag_id2"])

        fil = open(filename, "w")

        for atm, U in tls_group.iter_atm_Utls():
            if atm.name not in ["N", "CA", "C"]:
                continue

            utls_temp_factor = U2B * trace(U)/3.0
            bdiff = atm.temp_factor - utls_temp_factor
            
            try:
                fidf = float(str(atm.fragment_id))
            except ValueError:
                continue

            if atm.name=="N": fidf -= 0.33
            if atm.name=="C": fidf += 0.33

            fil.write("%f %f\n" % (fidf, bdiff))

        fil.close()

        return filename



_UISO_VS_UTLSISO_HISTOGRAM_TEMPLATE = """\
set xlabel "B_{obs} - B_{calc}"
set ylabel "Number of Atoms"
set style line 1 lc rgb "<rgb>" lw 3
set term png enhanced font "<font>" <fontsize>
set output "<pngfile>"
set title "<title>"
plot "<txtfile>" using 1:2 ls 1 notitle with histeps
"""

class UIso_vs_UtlsIso_Hisotgram(GNUPlot):
    def __init__(self, chainopt, tlsopt, tls):
        ## generate data and png paths
        basename  = "%s_CHAIN%s_TLS%s_%s_BoBc" % (
            chainopt["struct_id"],
            chainopt["chain_id"],
            tls["frag_id1"],
            tls["frag_id2"])
        
        self.txt_path = "%s.txt" % (basename)
        self.png_path = "%s.png" % (basename)

        ## write out the data file
        tls_group = tls["tls_group"]

        ## create a histogram of (Uiso - Utls_iso)
        bdiff_min = 0.0
        bdiff_max = 0.0

        for atm, Utls in tls_group.iter_atm_Utls():
            u_tls_iso = (trace(Utls) / 3.0) * U2B
            bdiff = atm.temp_factor - u_tls_iso

            bdiff_min = min(bdiff_min, bdiff)
            bdiff_max = max(bdiff_max, bdiff)

        ## compute the bin width and range to bin over
        brange    = (bdiff_max - bdiff_min) + 2.0
        num_bins  = int(brange)
        bin_width = brange / float(num_bins)
        bins      = [0 for n in range(num_bins)]

        ## name the bins with their mean value
        bin_names = []
        for n in range(num_bins):
            bin_mean = bdiff_min + (float(n) * bin_width) + (bin_width / 2.0)
            bin_names.append(bin_mean)

        ## count the bins
        for atm, Utls in tls_group.iter_atm_Utls():
            u_tls_iso = (trace(Utls) / 3.0) * U2B
            bdiff = atm.temp_factor - u_tls_iso
            bin = int((bdiff - bdiff_min)/ bin_width)
            bins[bin] += 1

        ## write out the gnuplot input file
        fil = open(self.txt_path, "w")
        fil.write("## Histogram of atoms in the TLS group binned by\n")
        fil.write("## the difference of their isotropic tempature factors\n")
        fil.write("## from the isotropic values predicted from the TLS model.\n")
        fil.write("##\n")
        fil.write("## Structure ----------------: %s\n" % (
            chainopt["struct_id"]))
        fil.write("## Chain --------------------: %s\n" % (
            chainopt["chain_id"]))
        fil.write("## Number of TLS Groups -----: %d\n" % (tlsopt.ntls))
        fil.write("## TLS Group Residue Range --: %s-%s\n" % (
            tls["frag_id1"], tls["frag_id2"]))

        for i in range(len(bins)):
            fil.write("%f %d\n" % (bin_names[i], bins[i]))

        fil.close()

        ## modify script template
        script = _UISO_VS_UTLSISO_HISTOGRAM_TEMPLATE
        script = script.replace("<font>", GNUPLOT_FONT)
        script = script.replace("<fontsize>", GNUPLOT_FONT_SIZE)
        script = script.replace("<txtfile>", self.txt_path)
        script = script.replace("<pngfile>", self.png_path)

        title = "Histogram of Observed B_{iso} - Calculated TLS B_{iso} "\
                "for TLS Group %s%s-%s%s" % (
            tls["chain_id"], tls["frag_id1"], tls["chain_id"], tls["frag_id2"])
        script = script.replace("<title>", title)

        script = script.replace("<rgb>", tls["color"]["rgbs"])

        self.gnuplot_run(script, basename)


class TLSSegmentAlignmentPlot(object):
    """Step 1: add all chains, generate unique list of ordered fragment ids,
               and hash all chain+frag_id->tls color
       Step 2: generate graphs
    """
    def __init__(self):
        ## border pixels
        self.border_width = 3
        ## bars are 15 pixels heigh
        self.pheight    = ALIGN_HEIGHT
        ## spacing pixels between stacked bars
        self.spacing    = ALIGN_SPACING
        ## background color
        self.bg_color   = rgb_f2i((1.0, 1.0, 1.0))
        
        self.frag_list     = []
        self.configurations = []

    def add_tls_segmentation(self, chainopt, ntls):
        """Add a TLS optimization to the alignment plot.
        """
        tlsopt = chainopt["tlsopt"][ntls]
        
        ## get the list of TLS segments for the specified number of
        ## segments (ntls)
        tls_seg_desc = {}
        self.configurations.append(tls_seg_desc)
        tls_seg_desc["chainopt"] = chainopt
        tls_seg_desc["ntls"]     = ntls
        tls_seg_desc["tlsopt"]   = tlsopt
        
        ## update the master fragment_list
        self.__update_frag_list(chainopt["chain"], tlsopt)

    def __update_frag_list(self, chain, tlsopt):
        """Add any fragment_ids found in the tls segments to the master
        self.frag_list and sort it.
        """
        for frag in chain.iter_fragments():
            fid = FragmentID(frag.fragment_id)
            if fid not in self.frag_list:
                self.frag_list.append(fid)

        self.frag_list.sort()
       
    def plot(self, path):
        """Plot and write the png plot image to the specified path.
        """
        if len(self.frag_list)==0 or len(self.configurations)==0:
            return False
        
        nfrag = len(self.frag_list)
        target_width = 500
        fw = int(round(float(ALIGN_TARGET_WIDTH) / nfrag))
        one_frag_width = max(1, fw)

        ## calculate with pixel width/fragment
        ## adjust the width of the graph as necessary
        pheight = self.pheight
        img_width = (2 * self.border_width) + \
                    (one_frag_width * len(self.frag_list))
            
        ## calculate the totoal height of the image
        num_plots = len(self.configurations)
        img_height = (2 * self.border_width) + \
                     (pheight * num_plots) + \
                     (self.spacing * (num_plots-1)) 

        ## create new image and 2D drawing object
        assert img_width>0 and img_height>0
        image = Image.new("RGBA", (img_width, img_height), self.bg_color)
        idraw = ImageDraw.Draw(image)
        idraw.setfill(True)

        ## draw plots
        for i in range(len(self.configurations)):
            tls_seg_desc = self.configurations[i]
            
            xo = self.border_width
            yo = self.border_width + (i * pheight) + (i * self.spacing)

            self.__plot_segmentation(
                idraw,
                img_width - 2*self.border_width,
                one_frag_width,
                (xo, yo),
                tls_seg_desc)

        image.save(path, "png")
        return True

    def __plot_segmentation(self, idraw, pwidth, fwidth, offset, tls_seg_desc):
        pheight = self.pheight
        nfrag   = len(self.frag_list)

        ## x/y offsets
        xo, yo = offset
        
        ## iterate over tls segments, draw and color
        tlsopt = tls_seg_desc["tlsopt"]

        ## draw a gray background for the lengt of the chain;
        ## the colored TLS segments will be drawn on top of this
        chain = tls_seg_desc["chainopt"]["chain"]
        fid1 = FragmentID(chain[0].fragment_id)
        fid2 = FragmentID(chain[-1].fragment_id)

        i1 = self.frag_list.index(fid1)
        i2 = self.frag_list.index(fid2)

        x1 = i1       * fwidth
        x2 = (i2 + 1) * fwidth

        idraw.setink((128, 128, 128))
        outline_width = 1
        idraw.rectangle((x1 + xo - outline_width,
                         yo - outline_width,
                         x2 + xo + outline_width,
                         pheight + yo + outline_width))

        ## draw colored TLS segments
        for tls in tlsopt.tls_list:
            fid1 = FragmentID(tls["frag_id1"])
            fid2 = FragmentID(tls["frag_id2"])

            i1 = self.frag_list.index(fid1)
            i2 = self.frag_list.index(fid2)

            x1 = i1       * fwidth
            x2 = (i2 + 1) * fwidth

            idraw.setink(tls["color"]["rgbi"])
            idraw.rectangle((x1+xo, yo, x2+xo, pheight+yo))



class Report(object):
    """Base class of HTML Report generating objects.
    """
    def html_head(self, title):
        """Header for all HTML pages.
        """
        x  = '<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" '
        x += '"http://www.w3.org/TR/html4/loose.dtd">\n\n'
        
        x += '<html>\n'
        x += '<head>\n'
        x += '  <title>%s</title>\n' % (title)
        x += '  <style type="text/css" media=screen>\n'
        x += '  <!-- \n'
        x += '  BODY {background-color:white;'
        x += '        margin-left:5%;margin-right:5%;'
        x += '        border-left:5%;border-right:5%;'
        x += '        margin-top:2%;border-top:2%;}\n'
        x += '  -->\n'
        x += '  </style>\n'

        x += '</head>\n'
        x += '<body>\n'
        return x

    def html_title(self, title):
        x  = ''
        x += '<center>'
        x += '<h1>%s</h1>' % (title)
        x += '</center>'
        return x

    def html_foot(self):
        """Footer for all HTML pages.
        """
        x  = ''
        x += '<br>'
        x += '<center><small>'
        x += 'TLSMD Version v%s Released %s ' % (
            GLOBALS["VERSION"], GLOBALS["RELEASE_DATE"])
        x += 'by %s <i>%s</i>' % (GLOBALS["AUTHOR"], GLOBALS["EMAIL"])
        x += '</small></center>'
        x += '</body></html>\n'
        return x


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
        chainopt["max_ntls"]     = 20
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
        for i in range(len(COLORS)):
            cname, rgbf = COLORS[i]
            
            color = {}
            self.colors.append(color)

            color["index"] = i
            color["name"]  = cname
            color["rgbf"]  = rgbf
            color["rgbi"]  = rgb_f2i(rgbf)

            rgbs = "#%2x%2x%2x" % rgb_f2i(rgbf)
            rgbs = rgbs.replace(" ", "0")
            color["rgbs"]  = rgbs

            ## generate thumbnail image
            color["thumbnail_path"] = os.path.join(
                thumbnail_dir, "%s.png" % (color["name"]))

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
        SaveStructure(fil=self.struct_path, struct=self.struct)

        ## generate small .png images so  they can be placed in the
        ## TLS group tables to identify the TLS group tabular data
        ## with the generated visualization
        self.init_colors()

        ## all TLSGraph objects get their calculations out of the
        ## way before writing HTML
        chainopt_list = []
        for chain in self.chains:
            begin_chain_timing(chain.chain_id)
            chainopt = self.init_chain_optimization(chain)
	    end_chain_timing(chain.chain_id)
            chainopt_list.append(chainopt)

        ## write out all TLSGraph reports
        for chainopt in chainopt_list:
            self.write_tls_chain_optimization(chainopt)

        ## a report page comparing the tls group segments of all
        ## chains aginst eachother
        self.write_multi_chain_alignment(chainopt_list)
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
        title = "TLSMD Rigid Body Analysis of %s" % (self.struct_id)

        x  = ''
        x += self.html_head(title)
        x += self.html_title(title)

        ## OPTIMIZATION PARAMETERS
        x += self.html_globals()
        x += '<br>'
        
        ## MOTION ANALYSIS
        x += '<center>'
        x += '<h3>TLS Partitions and Motion Analysis of Individual Chains</h3>'
        x += '</center>'
        x += '<p>%s</p>' % (MOTION_ANALYSIS_TEXT)

        for xdict in self.pages_chain_motion_analysis:
            x += '<p><a href="%s">%s</a></p>\n' % (
                xdict["href"], xdict["title"])

        x += '<br>'

        ## MULTI CHAIN ALIGNMENT
        x += '<center>'
        x += '<h3>Multi-Chain TLS Group Alignment</h3>'
        x += '</center>'
        x += '<p>%s</p>' % (MULTI_CHAIN_ALIGNMENT_TEXT)

        if self.page_multi_chain_alignment!=None:
            x += '<p><a href="%s">%s</a></p>\n' % (
                self.page_multi_chain_alignment["href"],
                self.page_multi_chain_alignment["title"])
        else:
            x += '<p><u>Only one chain was analyized in this '
            x += 'structure, so the multi-chain alignment analyisis '
            x += 'was not performed.'
            x += '</u></p>'

        x += '<br>'

        ## REFINEMENT PREP
        x += '<center>'
        x += '<h3>Use Optimal TLS Groups with Refmac5 TLS Refinement</h3>'
        x += '</center>'
        x += '<p>%s</p>' % (REFINEMENT_PREP_TEXT)

        x += '<p><a href="%s">%s</a></p>\n' % (
            self.page_refinement_prep["href"],
            self.page_refinement_prep["title"])
            
        x += self.html_foot()
        return x

    def html_globals(self):
        x  = ''
        
        x += '<center>'
        x += '<h3>Optimization Parameters</h3>'
        x += '</center>'

        x += '<p>%s</p>' % (OPTIMIZATION_PARAMS_TEXT)

        x += '<table border="1" cellpadding="3">'
        x += '<tr>'
        x += '<td>Form of TLS Model</td>'
        if GLOBALS["TLS_MODEL"]=="HYBRID":
            tls_model = 'For Input Structure with Isotropic ADPs'
        elif GLOBALS["TLS_MODEL"]=="ANISO":
            tls_model = 'For Input Structure with Anisotropic ADPs'
        else:
            tls_model = "Internal Error"
        x += '<td><b>%s</b></td>' % (tls_model)
        x += '</tr>'

        x += '<tr>'
        x += '<td>Least-Squares Weight Model</td>'
        if GLOBALS["WEIGHT_MODEL"]=="UNIT":
            weight = 'Unit Weights (All Weights 1.0)'
        elif GLOBALS["WEIGHT_MODEL"]=="IUISO":
            weight = 'Input Structure Atoms Weighted by '\
                     '<var>1.0/B<sub>iso</sub></var>'
        x += '<td><b>%s</b></td>' % (weight)
        x += '</tr>'


        x += '<tr>'
        x += '<td>Included Atoms</td>'
        x += '<td><b>%s</b></td>' % (GLOBALS["INCLUDE_ATOMS"])
        x += '</tr>'

        x += '<tr>'
        x += '<td>'
        x += 'Smallest Chain Subsegment Considered as a TLS Group'
        x += '</td>'
        x += '<td><b>%s Residues</b></td>' % (GLOBALS["MIN_SUBSEGMENT_SIZE"])
        x += '</tr>'

        x += '</table>'

        return x

    def write_tls_chain_optimization(self, chainopt):
        """Writes the HTML report analysis of a single TLS graphed chain.
        """
        begin_chain_timing(chainopt["chain_id"])
            
        path  = "%s_CHAIN%s_ANALYSIS.html" % (
            self.struct_id, chainopt["chain_id"])

        title = "Chain %s TLS Analysis" % (chainopt["chain_id"])

        self.pages_chain_motion_analysis.append(
            {"title": title,
             "href":  path })

        fil = open(path, "w")
        fil.write(self.html_tls_chain_optimization(chainopt))
        fil.close()

        end_chain_timing(chainopt["chain_id"])
        
    def html_tls_chain_optimization(self, chainopt):
        """Generates and returns the HTML string report analysis of a
        single TLS graphed chain.
        """
        title = "Chain %s TLS Analysis of %s" % (
            chainopt["chain_id"], self.struct_id)
        
        x  = self.html_head(title)
        x += self.html_title(title)
        x += '<center><a href="index.html">Back to Index</a></center>'
        x += '<br>\n'

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
        gp = LSQR_vs_TLS_Segments_Plot(chainopt)

        x  = ''
        x += '<center><h3>Chain %s Optimization Residual</h3></center>\n' % (
            chainopt["chain_id"])

        x += '<center>'
        x += '<table>'
        x += '<tr><td align="center">'
        x += '<img src="%s" alt="LSQR Plot">' % (gp.png_path)
        x += '</td></tr>'
        x += '<tr><td><p>%s</p></td></tr>' % (LSQR_CAPTION)
        x += '</table>'
        x += '</center>'
        
        return x

    def html_chain_alignment_plot(self, chainopt):
        """generate a plot comparing all segmentations
        """
        plot = TLSSegmentAlignmentPlot()
        
        for ntls, tlsopt in chainopt["ntls_list"]:
            plot.add_tls_segmentation(chainopt, ntls)

        ## create filename for plot PNG image file
        plot_path = "%s_CHAIN%s_ALIGN.png" % (
            self.struct_id, chainopt["chain_id"])
        
        plot.plot(plot_path)

        x  = ''
        x += '<center><h3>Chain %s TLS Segment Sequence '\
             'Alignment</h3></center>\n' % (chainopt["chain_id"])
        
        x += '<center>'
        x += '<table border="1">'
        x += '<tr><th>TLS Groups</th>'
        x += '<th>Chain %s Sequence Alignment</th></tr>'% (
            chainopt["chain_id"])
        x += '<tr>'
        
        x += '<td align="right">'
        x += '<table border="0" cellspacing="0" cellpadding="0">'

        for ntls, tlsopt in chainopt["ntls_list"]:
            x += '<tr><td align="right" valign="middle" height="20">'\
                 '<font size="-20">'\
                 '<a href="#NTLS%d">%d</a>'\
                 '</font></td></tr>' % (ntls, ntls)

        x += '</table>'
        x += '</td>'

        x += '<td><img src="%s" alt="Sequence Alignment Plot"></td>' % (
            plot_path)

        x += '</tr>'
        x += '</table>'
        x += '</center>'

        x += '<p>%s</p>' % (SEG_ALIGN_CAPTION)
        
        return x

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
        pml_path, png_path = self.raster3d_render_tls_graph_path(
            chainopt, tlsopt, ntls)

        ## JMol Viewer Page
        jmol_path = self.jmol_html(chainopt, tlsopt, ntls)
        jmol_animate_path = self.jmol_animate_html(chainopt, tlsopt)

        ## tlsout file
        tlsout_path = self.write_tlsout_file(chainopt, tlsopt, ntls)

        ## detailed analyisis of all TLS groups
        analysis_path = self.chain_ntls_analysis(chainopt, tlsopt)

        f1 = '<font size="-5">'
        f2 = '</font>'

        x  = ''
        x += '<hr>'
        x += '<center style="page-break-before: always">\n'
        x += '<h3><a name="NTLS%d">'\
             'Optimal TLS Group Partition using %d Groups</a></h3>\n' % (
            ntls, ntls)

        ## navigation links
        x += '<a href="%s">Motion and Error Analysis</a>' % (analysis_path)

        x += '&nbsp;&nbsp;&nbsp;&nbsp;'
         
        x += '<a href="." onClick="'\
             'window.open('\
             '&quot;%s&quot;,'\
             '&quot;&quot;,'\
             '&quot;width=%d,height=%d,screenX=10,'\
             'screenY=10,left=10,top=10&quot;);'\
             'return false;">View with JMol</a>' % (
            jmol_path, JMOL_SIZE, JMOL_SIZE)

        x += '&nbsp;&nbsp;&nbsp;&nbsp;'

        x += '<a href="." onClick="'\
             'window.open('\
             '&quot;%s&quot;,'\
             '&quot;&quot;,'\
             '&quot;width=%d,height=%d,screenX=10,'\
             'screenY=10,left=10,top=10&quot;);'\
             'return false;">Animate Screw Displacement with JMol</a>' % (
            jmol_animate_path, JMOL_SIZE, JMOL_SIZE)

        x += '<br>'
        x += '<a href="%s">Download TLSOUT File for TLSView</a>' % (
            tlsout_path)
        x += '&nbsp;&nbsp;&nbsp;&nbsp;'
        x += '<a href="%s">Generate PDBIN/TLSIN Files for REFMAC5</a>' % (
            "%s_REFINEMENT_PREP.html" % (self.struct_id))

        x += '</center>\n'

        ## raytraced image
        x += '<center><img src="%s" alt="iAlt"></center><br>\n' % (png_path)

        ## TLS group table
        x += '<table width="100%" border=1>\n'

        x += '<tr>\n'
        x += '<th align="center" colspan="12">'
        x += 'Analysis with %d TLS Groups' % (ntls)
        x += '</th>\n'
        x += '</tr>\n'

        x += '<tr>\n'
        x += '<th colspan="6">Input Structure</th>\n'
        x += '<th colspan="6">TLS Predictions</th>\n'
        x += '</tr>\n'

        x += ' <tr align="left">\n'
        x += '  <th>%sColor%s</th>\n' % (f1, f2)
        x += '  <th>%sSegment%s</th>\n' % (f1, f2)
        x += '  <th>%sResidues%s</th>\n' % (f1, f2)
        x += '  <th>%sAtoms%s</th>\n'  % (f1, f2)
        x += '  <th>%s&#60;B&#62;%s</th>\n'  % (f1, f2)
        x += '  <th>%s&#60;Aniso&#62;%s</th>\n' % (f1, f2)
        x += '  <th>%sLSQR%s</th>\n' % (f1, f2)
        x += '  <th>%sLSQR/Res%s</th>\n' % (f1, f2)
        x += '  <th>%seval(T<sup>r</sup>) <var>B</var>%s</th>\n' % (f1, f2)
        x += '  <th>%seval(L) <var>DEG<sup>2</sup></var>%s</th>\n' % (f1, f2)
        x += '  <th>%s&#60;B&#62;%s</th>\n'  % (f1, f2)
        x += '  <th>%s&#60;Aniso&#62;%s</th>\n' % (f1, f2)
        x += ' </tr>\n'

        for tls in tlsopt.tls_list:
            tls_group = tls["tls_group"]
            tls_info  = tls["tls_info"]

            L     = tls_group.L * RAD2DEG2
            L_ev  = [ev for ev in eigenvalues(L)]
            L_ev.sort()
            L_ev.reverse()
            
            T_red    = tls_info["rT'"]
            T_red_ev = [ev for ev in eigenvalues(T_red)]
            T_red_ev.sort()
            T_red_ev.reverse()

            x += '<tr>\n'

            x += '<td align="center" valign="middle">'\
                 '<img src="%s" alt="%s"></td>\n' % (
                tls["color"]["thumbnail_path"],
                tls["color"]["name"])

            x += '<td>%s%s-%s%s</td>\n' % (
                f1, tls["frag_id1"], tls["frag_id2"], f2)

            x += '<td>%s%d%s</td>\n'    % (
                f1, len(tls["segment"]), f2)

            x += '<td>%s%d%s</td>\n'    % (
                f1, len(tls_group), f2)

            x += '<td>%s%5.1f%s</td>\n' % (
                f1, tls_info["exp_mean_temp_factor"], f2)

            x += '<td>%s%4.2f%s</td>\n' % (
                f1, tls_info["exp_mean_anisotropy"], f2)

            x += '<td>%s%6.4f%s</td>\n' % (
                f1, tls["lsq_residual"], f2)

            x += '<td>%s%6.4f%s</td>\n' % (
                f1, tls["lsq_residual_per_res"], f2)

            x += '<td>%s%5.1f<br>%5.1f<br>%5.1f%s</td>\n' % (
                f1,
                T_red_ev[0]*U2B,
                T_red_ev[1]*U2B,
                T_red_ev[2]*U2B,
                f2)

            x += '<td>%s%5.2f<br>%5.2f<br>%5.2f%s</td>\n' % (
                f1, L_ev[0], L_ev[1], L_ev[2], f2)

            x += '<td>%s%5.1f%s</td>\n' % (
                f1, tls_info["tls_mean_temp_factor"], f2)
            
            x += '<td>%s%4.2f%s</td>\n' % (
                f1, tls_info["tls_mean_anisotropy"], f2)

            x += '</tr>\n'

        x += '</table>\n'
        x += '<br clear="all">\n'
        
        return x

    def raster3d_render_tls_graph_path(self, chainopt, tlsopt, ntls):
        """Render TLS visualizations using Raster3D.
        """
        basename = "%s_CHAIN%s_NTLS%d" % (
            self.struct_id, chainopt["chain_id"], ntls)

        png_path = "%s.png"   % (basename)

        start_timing()
        print "Raster3D: rendering %s..." % (basename)

        struct_id = self.struct_id
        chain     = chainopt["chain"]
        chain_id  = chainopt["chain_id"]

        driver = Raster3DDriver()

        ## XXX: Size hack: some structures have too many chains,
        ## or are just too large
        show_chain = {}
        for chx in self.struct.iter_chains():
            if chx.chain_id==chain_id:
                show_chain[chx.chain_id] = True
                continue
            
            if chx.count_fragments()>=500:
                show_chain[chx.chain_id] = False
                continue
            
            show_chain[chx.chain_id] = True
        ## end size hack

        viewer = GLViewer()
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
        gl_struct.glo_update_properties_path(
            "gl_axes/visible", False)
        gl_struct.glo_update_properties_path(
            "gl_unit_cell/visible", False)

        ## setup base structural visualization
        for gl_chain in gl_struct.glo_iter_children():
            if not isinstance(gl_chain, GLChain):
                continue

            ## chain is hidden
            if show_chain.get(gl_chain.chain.chain_id, False)==False:
                gl_chain.properties.update(visible=False)
                continue
            
            gl_chain.properties.update(
                oatm_visible       = False,
                side_chain_visible = False,
                hetatm_visible     = True,
                color              = "0.20,0.20,0.20",
                lines              = False,
                ball_stick         = True,
                ball_radius        = 0.20,
                stick_radius       = 0.20 )

            ## make chains other than the one we are analyizing visible,
            ## but pale
            if gl_chain.chain.chain_id!=chain_id:
                gl_chain.properties.update(
                    ball_radius  = 0.40,
                    stick_radius = 0.40,
                    color        = "0.9,0.9,0.9")

        ## add the TLS group visualizations
        for tls in tlsopt.tls_list:
            if tls["method"]!="TLS":
                continue
            
            tls_name = "TLS_%s_%s" % (
                tls["frag_id1"], tls["frag_id2"])
            
            gl_tls_group = GLTLSGroup(
                oatm_visible       = False,
                side_chain_visible = False,
                hetatm_visible     = True,
                adp_prob           = ADP_PROB,
                L1_visible         = True,
                L2_visible         = True,
                L3_visible         = True,
                L_axis_scale       = 2.0,
		both_phases        = True,
                tls_group          = tls["tls_group"],
                tls_info           = tls["tls_info"],
                tls_name           = tls_name,
                tls_color          = tls["color"]["name"])

            gl_struct.glo_add_child(gl_tls_group)

        ## set visualization: TLS traced surface
        for gl_tls_group in gl_struct.glo_iter_children():
            if not isinstance(gl_tls_group, GLTLSGroup):
                continue

            gl_tls_group.gl_atom_list.properties.update(
                trace_radius = 0.075)
            
            gl_tls_group.glo_update_properties(
                time = 0.25)

        driver.glr_set_render_png_path(png_path)
        viewer.glv_render_one(driver)
        print end_timing()

        return "", png_path

    def write_tls_pdb_file(self, chainopt, tlsopt, ntls):
        """Write out a PDB file with the TLS predicted anisotropic ADPs for
        this segmentation.
        """
        basename = "%s_CHAIN%s_NTLS%d_UTLS"  % (
            self.struct_id, chainopt["chain_id"], ntls)
        pdb_path = "%s.pdb" % (basename)

        ## temporarily set the atom temp_factor and U tensor to the Utls value
        old_temp_factor = {}
        old_U = {}
        for tls in tlsopt.tls_list:
            tls_group = tls["tls_group"]
            
            for atm, Utls in tls_group.iter_atm_Utls():
                old_temp_factor[atm] = atm.temp_factor
                old_U[atm] = atm.U

                atm.temp_factor = U2B * (trace(Utls)/3.0)
                atm.U = Utls

        SaveStructure(fil=pdb_path, struct=self.struct)

        ## restore atom temp_factor and U
        for atm, temp_factor in old_temp_factor.items():
            atm.temp_factor = temp_factor
            atm.U = old_U[atm]

    def write_tlsout_file(self, chainopt, tlsopt, ntls):
        """Writes the TLSOUT file for the segmentation.
        """
        basename = "%s_CHAIN%s_NTLS%d" % (
            self.struct_id, chainopt["chain_id"], ntls)
        tlsout_path = "%s.tlsout" % (basename)

        struct_id = self.struct_id
        chain_id  = chainopt["chain_id"]

        tls_file = TLSFile()
        tls_file.set_file_format(TLSFileFormatTLSOUT())

        for tls in tlsopt.tls_list:
            ## don't write out bypass edges
            if tls["method"]!="TLS":
                continue
            
            tls_desc = TLSGroupDesc()
            tls_file.tls_desc_list.append(tls_desc)
            
            tls_desc.set_tls_group(tls["tls_group"])
            tls_desc.add_range(
                chain_id, tls["frag_id1"],
                chain_id, tls["frag_id2"], "ALL")

        tls_file.save(open(tlsout_path, "w"))

        return tlsout_path

    def chain_ntls_analysis(self, chainopt, tlsopt):
        """Generate ntls optimization constraint report and free memory.
        """
        report = ChainNTLSAnalysisReport(chainopt, tlsopt, tlsopt.ntls)
        url = report.url
        report = None

        import gc
        gc.collect()

        return url

    def jmol_html(self, chainopt, tlsopt, ntls):
        """Writes out the HTML page which will display the
        structure using the JMol Applet.
        """
        jmol_path = "%s_CHAIN%s_NTLS%d_JMOL.html"  % (
            self.struct_id, chainopt["chain_id"], ntls)

        ## create the JMol script using cartoons and consisant
        ## coloring to represent the TLS groups
        js  = ''
        js += 'load %s;' % (self.struct_path)
        js += 'select *;'
        js += 'cpk off;'
        js += 'wireframe off;'
        js += 'select protein;'
        js += 'cartoon on;'

        ## loop over TLS groups and color
        for tls in tlsopt.tls_list:
            js += 'select %s-%s:%s;' % (
                tls["frag_id1"], tls["frag_id2"], tls["chain_id"])
            js += 'color [%d,%d,%d];' % (tls["color"]["rgbi"])

        ## select non-protein non-solvent and display
        js += 'select not protein and not solvent;'
        js += 'color CPK;'
        js += 'wireframe on; wireframe 0.5;'
        js += 'spacefill 80%;'
        js += 'spacefill on;'
        
        ## write the HTML page to render the script in
        x  = ''
        x += '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" '\
             '"http://www.w3.org/TR/html4/strict.dtd">'
        x += '<html>'
        x += '<head>'
        x += '<title>Chain %s using %d TLS Groups</title>' % (
            chainopt["chain_id"], ntls)
        x += '<script type="text/javascript" src="%s/Jmol.js">' % (JMOL_DIR)
        x += '</script>'
        x += '</head>'
        x += '<body>'
        x += '<script type="text/javascript">'
        x += 'jmolInitialize("%s");' % (JMOL_DIR)
        x += 'jmolSetAppletColor("white");'
        x += 'jmolApplet(%d, "%s");' % (JMOL_SIZE, js)
        x += '</script>'
        x += '</body>'
        x += '</html>'

        open(jmol_path, "w").write(x)
        return jmol_path

    def jmol_animate_html(self, chainopt, tlsopt):
        """Writes out the HTML page which will display the
        structure using the JMol Applet.
        """
        basename = "%s_CHAIN%s_NTLS%d_ANIMATE" % (
            self.struct_id, chainopt["chain_id"], tlsopt.ntls)

        html_path = "%s.html" % (basename)
        pdb_path  = "%s.pdb" % (basename)

        ## gerate animation PDB file

        try:
            print "TLSAnimate: creating animation PDB file..."
            start_timing()
            tlsa = TLSAnimate(self.struct, chainopt, tlsopt)
            tlsa.construct_animation(pdb_path)
            print end_timing()
        except TLSAnimateFailure:
            pass
        
        ## create the JMol script using cartoons and consisant
        ## coloring to represent the TLS groups
        js  = ''
        js += 'load %s;' % (pdb_path)
        js += 'select *;'
        js += 'cpk off;'
        js += 'wireframe off;'
        js += 'select protein;'
        js += 'trace on;'

        ## loop over TLS groups and color
        for tls in tlsopt.tls_list:
            chain_ids = [tlsa.L1_chain.chain_id,
                         tlsa.L2_chain.chain_id,
                         tlsa.L3_chain.chain_id]

            for chain_id in chain_ids:
                js += 'select %s-%s:%s;' % (
                    tls["frag_id1"], tls["frag_id2"], chain_id)
                js += 'color [%d,%d,%d];' % (tls["color"]["rgbi"])

        ## select non-protein non-solvent and display
        js += 'select not protein and not solvent;'
        js += 'color CPK;'
        js += 'wireframe on;'
        js += 'wireframe 0.5;'
        js += 'spacefill 80%;'
        js += 'spacefill on;'

        js += 'anim fps 2;'
        js += 'anim mode loop 0 0;'
        js += 'anim on;'

        ## write the HTML page to render the script in
        x  = ''
        x += '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" '\
             '"http://www.w3.org/TR/html4/strict.dtd">'
        x += '<html>'
        x += '<head>'
        x += '<title>Chain %s using %d TLS Groups</title>' % (
            chainopt["chain_id"], tlsopt.ntls)
        x += '<script type="text/javascript" src="%s/Jmol.js">' % (JMOL_DIR)
        x += '</script>'
        x += '</head>'
        x += '<body>'
        x += '<script type="text/javascript">'
        x += 'jmolInitialize("%s");' % (JMOL_DIR)
        x += 'jmolSetAppletColor("white");'
        x += 'jmolApplet(%d, "%s");' % (JMOL_SIZE, js)
        x += '</script>'
        x += '</body>'
        x += '</html>'

        ## manually free memory
        tlsa = None
        import gc
        gc.collect()

        open(html_path, "w").write(x)
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
        x += '<br>\n'

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
                        atm.temp_factor = U2B * (
                            Utls[0,0] + Utls[1,1] + Utls[2,2]) / 3.0
                        atm.U = Utls

            ## save the structure file
            SaveStructure(fil=pdb_path, struct=self.struct)

            ## restore atom temp_factor and U
            for atm, temp_factor in old_temp_factor.items():
                atm.temp_factor = temp_factor
                atm.U = old_U[atm]

            ## generate the TLS segmentation alignment plot for all chains
            plot = TLSSegmentAlignmentPlot()
            chain_id_list = []

            for chainopt in chainopt_list:
                if not chainopt["tlsopt"].has_key(ntls):
                    continue
                
                chain_id_list.append(chainopt["chain_id"])
                plot.add_tls_segmentation(chainopt, ntls)

            plot_path = "%s_ALIGN.png" % (basename)
            plot.plot(plot_path)

            ## write HTML
            x += '<h3>Chains Alignment using %d TLS Groups</h3>\n' % (
                ntls)

            ## plot table
            x += '<table border="1">'
            x += '<tr><th>Chain</th><th>Chain Alignment</th></tr>'
            x += '<tr>'
        
            x += '<td align="center">'
            x += '<table border="0" cellspacing="0" cellpadding="0">'

            for chain_id in chain_id_list:
                x += '<tr><td align="right" valign="middle" height="20">'\
                     '<font size="-20">%s</font></td></tr>' % (chain_id)

            x += '</table>'
            x += '</td>'

            x += '<td><img src="%s" alt="Segmentation Plot"></td>' % (
                plot_path)

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
        x += '<br>\n'

        x += '<form enctype="multipart/form-data" '\
             'action="%s" method="post">' % (
            GLOBALS["REFINEPREP_URL"])
        x += '<input type="hidden" name="job_id" value="%s">' % (
            GLOBALS["JOB_ID"])

        x += '<p>%s</p>' % (REFINEMENT_PREP_INFO)
        
        x += '<center><table><tr><td>'
        
        plot = LSQR_vs_TLS_Segments_All_Chains_Plot(chainopt_list)
        x += '<img src="%s" alt="LSQR Residual">' % (plot.png_path)

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
        self.dir   = "%s_CHAIN%s_NTLS%d"  % (
            self.struct_id, self.chain_id, self.ntls)
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
        title = "Chain %s Analysis using %d TLS Groups" % (
            self.chain_id, self.ntls)
        
        x  = self.html_head(title)
        x += self.html_title(title)

        x += '<center>'
        x += '<a href="../index.html">Back to Index</a>'
        x += '&nbsp;&nbsp;&nbsp;&nbsp;'
        path = "%s_CHAIN%s_ANALYSIS.html" % (
            self.struct_id, self.chain_id)
        x += '<a href="../%s">Back to Chain %s Analysis</a>' % (
            path, self.chain_id)
        x += '</center>'
        
        x += '<br>\n'

        x += self.html_translation_analysis()
        x += self.html_libration_analysis()
        x += self.html_fit_analysis()
        
        for tls in self.tlsopt.tls_list:
            ## don't write out bypass edges
            if tls["method"]!="TLS":
                continue
            x += self.html_tls_fit_histogram(tls)

        ## write out the HTML page
        x += self.html_foot()
        open(self.index, "w").write(x)

    def html_translation_analysis(self):
        """Perform a translation analysis of the protein chain as
        spanned by the tlsopt TLS groups.
        """
        x  = ''
        x += '<center>'
        x += '<h3>Translation Analysis of T<sup>r</sup></h3>'
        x += '</center>\n'

        tanalysis = TranslationAnalysis(self.chainopt, self.tlsopt)
        
        x += '<center>'
        x += '<img src="%s" alt="Translation Analysis">' % (tanalysis.png_path)
        x += '</center>\n'
        x += '<p>%s</p>' % (TRANSLATION_GRAPH_CAPTION)
        
        return x

    def html_libration_analysis(self):
        """Perform a libration analysis of the protein chain as
        spanned by the tlsopt TLS groups.
        """
        x  = ''
        x += '<center><h3>Screw Displacment Analysis</h3></center>\n'

        libration_analysis = LibrationAnalysis(self.chainopt, self.tlsopt)

        x += '<center>'
        x += '<img src="%s" alt="Libration Analysis">' % (
            libration_analysis.png_path)
        x += '</center>\n'
        x += '<p>%s</p>' % (LIBRATION_GRAPH_CAPTION)
        
        return x

    def html_fit_analysis(self):
        """Perform a fit analysis of the protein chain as
        spanned by the tlsopt TLS groups.
        """
        x  = ''
        x += '<center><h3>Main Chain TLS Fit Analysis</h3></center>\n'

        fit_analysis = FitAnalysis(self.chainopt, self.tlsopt)
        
        x += '<center>'
        x += '<img src="%s" alt="Fit Analysis">' % (
            fit_analysis.png_path)
        x += '</center>\n'
        x += '<p>%s</p>' % (FIT_GRAPH_CAPTION)

        return x

    def html_tls_fit_histogram(self, tls):
        """A complete analysis of a single TLS group output as HTML.
        """
        x  = ''
        x += '<center>'
        x += '<h3>Distribution Histogram of TLS Group %s%s-%s%s</h3>' % (
            self.chain_id, tls["frag_id1"], self.chain_id, tls["frag_id2"])
        x += '</center>'

        ## histogrm of atomic U_ISO - U_TLS_ISO
        his = UIso_vs_UtlsIso_Hisotgram(self.chainopt, self.tlsopt, tls)

        x += '<center><img src="%s" alt="iAlt"></center>\n' % (
            his.png_path)

        return x



