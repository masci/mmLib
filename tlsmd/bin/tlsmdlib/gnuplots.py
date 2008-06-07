## TLS Motion Determination (TLSMD)
## Copyright 2002-2008 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import subprocess
import itertools
import numpy

from mmLib import Constants, Colors, Gaussian, AtomMath, TLS

import console
import table
import misc
import conf
import tls_calcs


def FormatFigureHTML(title, caption, figure_html):

    if caption:
        l = ['<tr><td align="center">',
             '<p style="padding:2%%;background-color:#eeeeee;border-style:dashed;border-width:thin;border-color:black">%s</p>' % (caption),
             '</td></tr>']
        cap_html = "".join(l)
    else:
        cap_html = ""

    l = ['<table style="width:80%">',
         '<tr><th style="font-size:large">%s</th></tr>' % (title),
         '<tr><td align="center">',
         figure_html,
         '</td></tr>',
         cap_html,
         '</table>']

    return "".join(l)


class GNUPlot(object):
    """Provides useful methods for subclasses which need to run gnuplot.
    """
    def __init__(self, **args):
        self.gnuplot_path = args.get("gnuplot_path", conf.GNUPLOT)
        self.font_path    = args.get("font_path", conf.GNUPLOT_FONT)
        self.font_size    = args.get("fontsize", conf.GNUPLOT_FONT_SIZE)
        self.width        = args.get("width", conf.GNUPLOT_WIDTH)
        self.height       = args.get("height", conf.GNUPLOT_HEIGHT)

        self.basename     = None
        self.txt_path     = None
        self.plot_path    = None
        self.png_path     = None
        self.svg_path     = None

    def set_basename(self, basename):
        """Sets all output file names based on the plot's basename.
        """
        self.basename = basename
        self.txt_path = "%s.txt" % (basename)
        self.plot_path = "%s.plot" % (basename)
        self.png_path = "%s.png" % (basename)
        self.svg_path = "%s.svg" % (basename)

    def make_script(self):
        """Implement me.
        """
        pass

    def run_gnuplot(self, script):
        """Execute GNUPlot with the given script.
        """
        try:
            pobj = subprocess.Popen([self.gnuplot_path],
                                    stdin = subprocess.PIPE,
                                    stdout = subprocess.PIPE,
                                    stderr = subprocess.STDOUT,
                                    close_fds = True,
                                    bufsize = 8192)
        except OSError:
            console.stderrln("gnuplot failed to execute from path: %s" % (self.gnuplot_path))
            return
            
        pobj.stdin.write(script)
        pobj.stdin.close()
	pobj.stdout.read()
        pobj.stdout.close()
        pobj.wait()
        
    def output_png(self):
        """Runs gnuplot.  Expects self.plot_path and self.png_path to be set.
        """
        script0 = self.make_script()
        
        ## if a basename is given, then write the GnuPlot script as a file
        console.stdoutln("GNUPLot: Saving %s" % (self.png_path))

        ## set output size
        l = ['set term png font "%s" %d size %d,%d enhanced' % (self.font_path, self.font_size, self.width, self.height),
             'set output "%s"' % (self.png_path),
             '']
        
        script_png = "\n".join(l) + script0

        ## write a gnuplot script
        open(self.plot_path, "w").write(script_png)
        
        ## run gnuplot
        self.run_gnuplot(script_png)

        ## XXX: hack svg output
        if conf.globalconf.use_svg == True:
            l = ['set term svg size %d %d dynamic fsize 12 enhanced' % (self.width, self.height),
                 'set output "%s"' % (self.svg_path),
                 '']
            script_svg = "\n".join(l) + script0
            self.run_gnuplot(script_svg)

    def html_link(self, alt_text=None):
        if not alt_text: alt_text = self.basename

        if conf.globalconf.use_svg == True:
            l = ['<object type="image/svg+xml" data="./%s" width="%d" height="%d">' % (self.svg_path, self.width, self.height),
                 '<img src="%s" alt="%s">' % (self.png_path, alt_text),
                 '</object>']

            return "".join(l)
        else:
            return '<img src="%s" alt="%s">' % (self.png_path, alt_text)
        
    def html_markup(self, title, caption, alt_text = None):
        return FormatFigureHTML(title, caption, self.html_link(alt_text))
        

_LSQR_VS_TLS_SEGMENTS_TEMPLATE = """\
set xlabel "Number of TLS Segments"
set xrange [1:<nparts>]
set xtics 1
set ylabel "Residual"
set format y "%5.2f"
set style line 1 lw 3
set title "<title>"
plot "<txtfile>" using 1:2 title "Minimization (Weighted) Residual" ls 1 with linespoints
"""
        
class LSQR_vs_TLS_Segments_Plot(GNUPlot):
    ## TLSMD selects the optimal partition of a chain into 1 to 20 TLS groups
    ## by minimizing an overall residual function. This plot shows the value
    ## of the residual as a function of the number of TLS groups allowed in
    ## the partition. Adding additional TLS groups will always make this
    ## residual lower, but there is an issue of diminishing returns as you go
    ## to larger numbers of groups.
    ##
    ## NOTE (by Christoph):
    ## This class is called by def html_chain_lsq_residual_plot(self, chain) in html.py
    def __init__(self, chain,  **args):
        GNUPlot.__init__(self, **args)
        self.chain = chain
        self.output_png()

    def make_script(self):
        ## generate data and png paths
        basename = "%s_CHAIN%s_RESID" % (
            self.chain.partition_collection.struct.structure_id,
            self.chain.chain_id)
        self.set_basename(basename)

        tbl = table.StringTable(0, 3, "?",
                                column_titles = ["Number of TLS Groups", "RMSD B", "Residual"])
        for ntls, cpartition in self.chain.partition_collection.iter_ntls_chain_partitions():
            tbl.append_row(ntls, cpartition.rmsd_b(), cpartition.residual())

        open(self.txt_path, "w").write(str(tbl))

        ## modify script template
        script = _LSQR_VS_TLS_SEGMENTS_TEMPLATE
        script = script.replace("<nparts>", str(conf.globalconf.nparts))
        script = script.replace("<txtfile>", self.txt_path)
        script = script.replace("<title>", "Least Squares Residual vs. Number of TLS Segments for %s Chain %s " % (
            self.chain.partition_collection.struct.structure_id,
            self.chain.chain_id))

        return script


_LSQR_VS_TLS_SEGMENTS_ALL_CHAINS_TEMPLATE = """\
set xlabel "Number of TLS Segments"
set xrange [1:<nparts>]
set xtics 1
set ylabel "Minimization (Weighted) LSQR Residual"
set format y "%5.2f"
set title "<title>"
"""

class LSQR_vs_TLS_Segments_All_Chains_Plot(GNUPlot):
    def __init__(self, tlsmd_analysis, **args):
        GNUPlot.__init__(self, **args)
        self.tlsmd_analysis = tlsmd_analysis
        self.output_png()

    def make_script(self):
        struct_id = self.tlsmd_analysis.struct.structure_id
        
        ## generate data and png paths
        basename = "%s_RESID" % (struct_id)
        self.set_basename(basename)

        ## prepare gnuplot script
        script = _LSQR_VS_TLS_SEGMENTS_ALL_CHAINS_TEMPLATE
        script = script.replace("<nparts>", str(conf.globalconf.nparts))
        script = script.replace("<title>", "Least Squares Residual vs. Number of TLS Segments of %s" % (struct_id))

        ## re-use the data files of LSQRvsNTLS from the individual
        ## graphs; to do this the filenames have to be re-constructed
        plist = []
        for chain in self.tlsmd_analysis.iter_chains():
            filename = "%s_CHAIN%s_RESID.txt" % (struct_id, chain.chain_id)
            x = '"%s" using 1:2 title "Chain %s" lw 3 with linespoints' % (filename, chain.chain_id)
            plist.append(x)
            
        script += "plot " + ",\\\n    ".join(plist) + "\n"

	## EAM Feb 2008
	## Make a thumbnail (half-size) version for the summary page
	script += "set term png font '%s' 8 size 400,320 linewidth 0.5\n" % conf.GNUPLOT_FONT
	script += "set output 'summary.png'\n"
	script += "unset title; set ylabel 'Residual' offset 1; replot\n"

        console.stdoutln("GNUPLot: Saving summary.png") ## LOGLINE

        return script


_TRANSLATION_ANALYSIS_TEMPLATE = """\
set xlabel "Residue"
set xrange [<xrng1>:<xrng2>]
set ylabel "Angstroms Displacement"
set format y "%4.2f"
set title "<title>"
"""

class TranslationAnalysis(GNUPlot):
    def __init__(self, chain, cpartition, **args):
        GNUPlot.__init__(self, **args)
        self.chain = chain
        self.cpartition = cpartition
        self.output_png()

    def make_script(self):
        basename = "%s_CHAIN%s_NTLS%s_TRANSLATION" % (
            self.chain.struct.structure_id,
            self.chain.chain_id,
            self.cpartition.num_tls_segments())

        self.set_basename(basename)

        self.write_data_file()

        script = _TRANSLATION_ANALYSIS_TEMPLATE
        script = script.replace("<xrng1>", self.cpartition.first_frag_id())
        script = script.replace("<xrng2>", self.cpartition.last_frag_id())
        script = script.replace("<title>", "Translation Displacement Analysis of Atoms for %d TLS Groups" % (self.cpartition.num_tls_segments()))

        ## line style
        ls = 0
        for tls in self.cpartition.iter_tls_segments():
            ls += 1
            script += 'set style line %d lc rgb "%s" lw 3\n' % (ls, tls.color.rgbs)

        ## plot list
        plist = []
        ls = 0
        for itls in xrange(self.cpartition.num_tls_segments()):
            ls += 1
            
            for n in (0,1,2):
                col = 2 + 3*itls + n
                x = '"%s" using 1:%d smooth bezier notitle ls %d with lines' % (self.txt_path, col , ls)
                plist.append(x)

        script += "plot " + ",\\\n    ".join(plist) + "\n"
           
        return script

    def write_data_file(self):
        """Generate the data file and return the filename.
        """
        nrows = len(self.cpartition.chain)
        ncols = 1 + 3 * self.cpartition.num_tls_segments()
        tbl = table.StringTable(nrows, ncols, "?",
                                title = "TLS Model Translation Tensor RMSD Values of  Principal Components",
                                column_titles = ["Residue"])

        frag_id_iter = itertools.imap(lambda frag: frag.fragment_id, self.cpartition.chain.iter_fragments())
        tbl.set_column(0, 0, frag_id_iter)

        for itls, tls in enumerate(self.cpartition.iter_tls_segments()):
            tls_group = tls.tls_group
            tls_info = tls.model_tls_info
            O = tls_info["COR"]

            ## determine Tr translational eigenvalues
            t1 = Gaussian.GAUSS3C[conf.ADP_PROB] * tls_info["Tr1_rmsd"]
            t2 = Gaussian.GAUSS3C[conf.ADP_PROB] * tls_info["Tr2_rmsd"]
            t3 = Gaussian.GAUSS3C[conf.ADP_PROB] * tls_info["Tr3_rmsd"]

            for frag in tls.iter_fragments():
                i = frag.ifrag
                if t1 > 0.0: tbl[i, 1 + 3*itls] = "%6.4f" % (t1)
                if t2 > 0.0: tbl[i, 2 + 3*itls] = "%6.4f" % (t2)
                if t3 > 0.0: tbl[i, 3 + 3*itls] = "%6.4f" % (t3)
        
        open(self.txt_path, "w").write(str(tbl))


_LIBRATION_ANALYSIS_TEMPLATE = """\
set xlabel "Residue"
set xrange [<xrng1>:<xrng2>]
set ylabel "Angstroms Displacement"
set format y "%4.2f"
set title "<title>"
set datafile missing "?"
"""

class LibrationAnalysis(GNUPlot):
    def __init__(self, chain, cpartition, **args):
        GNUPlot.__init__(self, **args)
        self.chain = chain
        self.cpartition = cpartition
        self.output_png()

    def make_script(self):        
        basename = "%s_CHAIN%s_NTLS%s_LIBRATION" % (
            self.chain.struct.structure_id,
            self.chain.chain_id,
            self.cpartition.num_tls_segments())
        
        self.set_basename(basename)

        self.write_data_file()

        script = _LIBRATION_ANALYSIS_TEMPLATE
        script = script.replace("<xrng1>", self.cpartition.first_frag_id())
        script = script.replace("<xrng2>", self.cpartition.last_frag_id())
        script = script.replace("<title>", "Screw displacement analysis of backbone atoms using %d TLS Groups" % (self.cpartition.num_tls_segments()))

        ## line style
        ls = 0
        for tls in self.cpartition.iter_tls_segments():
            ls += 1
            script += 'set style line %d lc rgb "%s" lw 3\n' % (ls, tls.color.rgbs)

        ## plot list
        plist = []
        ls = 0
        for itls in xrange(self.cpartition.num_tls_segments()):
            ls += 1
            
            for n in (0,1,2):
                col = 2 + 3*itls + n
                x = '"%s" using 1:%d smooth bezier notitle ls %d with lines' % (self.txt_path, col , ls)
                plist.append(x)

        script += "plot " + ",\\\n    ".join(plist) + "\n"
           
        return script

    def write_data_file(self):
        """Generate the data file and return the filename.
        """
        nrows = len(self.cpartition.chain)
        ncols = 1 + 3 * self.cpartition.num_tls_segments()
        tbl = table.StringTable(nrows, ncols, "?")

        mpred = lambda f: f.fragment_id
        frag_id_iter = itertools.imap(mpred, self.cpartition.chain.iter_fragments())
        tbl.set_column(0, 0, frag_id_iter)

        for itls, tls in enumerate(self.cpartition.iter_tls_segments()):
            tls_group = tls.tls_group
            tls_info = tls.model_tls_info
            O = tls_info["COR"]

            for frag in tls.iter_fragments():
                atm = frag.get_atom("CA")
                if atm is None:
                    continue

                i = frag.ifrag

                for n, Lx_val, Lx_vec, Lx_rho, Lx_pitch in [
                    (0, "L1_eigen_val", "L1_eigen_vec", "L1_rho", "L1_pitch"),
                    (1, "L2_eigen_val", "L2_eigen_vec", "L2_rho", "L2_pitch"),
                    (2, "L3_eigen_val", "L3_eigen_vec", "L3_rho", "L3_pitch") ]:

                    Lval   = tls_info[Lx_val]
                    Lvec   = tls_info[Lx_vec]
                    Lrho   = tls_info[Lx_rho]
                    Lpitch = tls_info[Lx_pitch]

                    if numpy.allclose(Lval, 0.0):
                        continue

                    dvec = TLS.calc_LS_displacement(O, Lval, Lvec, Lrho, Lpitch, atm.position, conf.ADP_PROB)
                    tbl[i, 1 + 3*itls + n] = AtomMath.length(dvec)

        open(self.txt_path, "w").write(str(tbl))


_CA_TLS_DIFFERANCE_TEMPLATE = """\
set xlabel "Residue"
set xrange [<xrng1>:<xrng2>]
set ylabel "B_{obs} - B_{tls}"
set format y "%5.2f"
set title "<title>"
set datafile missing "?"
"""

class CA_TLS_Differance_Plot(GNUPlot):
    def __init__(self, chain, cpartition, **args):
        GNUPlot.__init__(self, **args)
        self.chain = chain
        self.cpartition = cpartition
        self.output_png()

    def make_script(self):
        basename = "%s_CHAIN%s_NTLS%s_CADIFF" % (
            self.chain.struct.structure_id,
            self.chain.chain_id,
            self.cpartition.num_tls_segments())

        self.set_basename(basename)

        self.write_data_file()

        script = _CA_TLS_DIFFERANCE_TEMPLATE
        script = script.replace("<xrng1>", self.cpartition.first_frag_id())
        script = script.replace("<xrng2>", self.cpartition.last_frag_id())
        script = script.replace("<title>", "Deviation of Observed CA B Factors from TLS Model for %d Group Partition" % (self.cpartition.num_tls_segments()))

        ## line style
        script += 'set style line 1 lc rgb "#000000" lw 1\n'
        
        ls = 1
        for itls, tls in enumerate(self.cpartition.iter_tls_segments()):
            ls += 1
            script += 'set style line %d lc rgb "%s" lw 3\n' % (ls, tls.color.rgbs)

        ## plot list
        plist = []
        plist.append("0.0 notitle ls 1 with lines")

        ls = 1
        for itls in xrange(self.cpartition.num_tls_segments()):
            ls += 1
            x = '"%s" using 1:%d notitle ls %d with lines' % (self.txt_path, itls+2, ls)
            plist.append(x)

        script += "plot " + ",\\\n    ".join(plist) + "\n"
           
        return script

    def write_data_file(self):
        nrows = len(self.cpartition.chain)
        ncols = self.cpartition.num_tls_segments() + 1
        tbl = table.StringTable(nrows, ncols, "?")

        frag_id_iter = itertools.imap(lambda frag: frag.fragment_id, self.cpartition.chain.iter_fragments())
        tbl.set_column(0, 0, frag_id_iter)
        
        for itls, tls in enumerate(self.cpartition.iter_tls_segments()):
            tls_group = tls.tls_group

            T = tls_group.itls_T
            L = tls_group.itls_L
            S = tls_group.itls_S
            O = tls_group.origin

            for frag in tls.iter_fragments():
                atm = frag.get_atom("CA")
                if atm is None:
                    continue
                i = frag.ifrag
                b_tls = Constants.U2B * TLS.calc_itls_uiso(T, L, S, atm.position - O)
                tbl[i, itls + 1] = atm.temp_factor - b_tls

        open(self.txt_path, "w").write(str(tbl))
        

_UISO_VS_UTLSISO_HISTOGRAM_TEMPLATE = """\
set xlabel "B_{obs} - B_{tls}"
set ylabel "Number of Atoms"
set style line 1 lc rgb "<rgb>" lw 3
set title "<title>"
plot "<txtfile>" using 1:2 ls 1 notitle with histeps
"""

class UIso_vs_UtlsIso_Histogram(GNUPlot):
    def __init__(self, chain, cpartition, tls, **args):
        GNUPlot.__init__(self, **args)
        self.chain = chain
        self.cpartition = cpartition
        self.tls = tls
        self.output_png()

    def make_script(self):
        tls = self.tls

        ## generate data and png paths
        basename  = "%s_CHAIN%s_TLS%s_BoBc" % (
            self.chain.struct.structure_id,
            self.chain.chain_id,
            tls.filename_label())

        self.set_basename(basename)

        ## write out the data file
        tls_group = tls.tls_group

        T = tls_group.itls_T
        L = tls_group.itls_L
        S = tls_group.itls_S
        O = tls_group.origin

        ## create a histogram of (Uiso - Utls_iso)
        bdiff_min = 0.0
        bdiff_max = 0.0

        for atm in tls_group:
            b_iso_tls = Constants.U2B * TLS.calc_itls_uiso(T, L, S, atm.position - O)
            bdiff = atm.temp_factor - b_iso_tls

            bdiff_min = min(bdiff_min, bdiff)
            bdiff_max = max(bdiff_max, bdiff)

        ## compute the bin width and range to bin over
        brange    = (bdiff_max - bdiff_min) + 2.0
        num_bins  = int(brange)
        bin_width = brange / float(num_bins)
        bins      = [0 for n in xrange(num_bins)]

        ## name the bins with their mean value
        bin_names = []
        for n in xrange(num_bins):
            bin_mean = bdiff_min + (float(n) * bin_width) + (bin_width / 2.0)
            bin_names.append(bin_mean)

        ## count the bins
        for atm in tls_group:
            b_iso_tls = Constants.U2B * TLS.calc_itls_uiso(T, L, S, atm.position - O)
            bdiff = atm.temp_factor - b_iso_tls
            bin = int((bdiff - bdiff_min)/ bin_width)
            bins[bin] += 1

        ## write out the gnuplot input file
        fil = open(self.txt_path, "w")
        fil.write("## Histogram of atoms in the TLS group binned by\n")
        fil.write("## the difference of their isotropic temperature factors\n")
        fil.write("## from the isotropic values predicted from the TLS model.\n")
        fil.write("##\n")
        fil.write("## Structure ----------------: %s\n" % (self.chain.struct.structure_id))
        fil.write("## Chain --------------------: %s\n" % (self.chain.chain_id))
        fil.write("## TLS Group Residue Range --: %s\n" % (tls.display_label()))

        for i in xrange(len(bins)):
            fil.write("%f %d\n" % (bin_names[i], bins[i]))

        fil.close()

        ## modify script template
        script = _UISO_VS_UTLSISO_HISTOGRAM_TEMPLATE
        script = script.replace("<txtfile>", self.txt_path)

        title = "Histogram of Observed B_{iso} - B_{tls} for TLS Group %s" % (tls.display_label())
        script = script.replace("<title>", title)
        script = script.replace("<rgb>", tls.color.rgbs)

        return script

     
_BMEAN_PLOT_TEMPLATE = """\
set xlabel "Residue"
set xrange [<xrng1>:<xrng2>]
set ylabel "<B>"
set yrange [0.0:]
set format y "%5.2f"
set title "<title>"
set datafile missing "?"
"""

class BMeanPlot(GNUPlot):
    def __init__(self, chain, cpartition, **args):
        GNUPlot.__init__(self, **args)
        self.chain = chain
        self.cpartition = cpartition
        self.tls_group_titles = True
        self.output_png()

    def make_script(self):
        basename = "%s_CHAIN%s_NTLS%s_BMEAN" % (
            self.chain.struct.structure_id,
            self.chain.chain_id,
            self.cpartition.num_tls_segments())

        self.set_basename(basename)

        ## write data file
        BISO_OBS = tls_calcs.calc_mean_biso_obs(self.chain)
        BISO = tls_calcs.calc_mean_biso_tls(self.chain, self.cpartition)

        nrows = len(self.cpartition.chain)
        ncols = 3 + self.cpartition.num_tls_segments()
        tbl = table.StringTable(nrows, ncols, "?")
        frag_id_iter = itertools.imap(lambda frag: frag.fragment_id, self.cpartition.chain.iter_fragments())
        tbl.set_column(0, 0, frag_id_iter)
        ifrag_iter = itertools.imap(lambda frag: frag.ifrag, self.cpartition.chain.iter_fragments())
        tbl.set_column(0, 1, ifrag_iter)
        tbl.set_column(0, 2, BISO_OBS)
        
        for itls, tls in enumerate(self.cpartition.iter_tls_segments()):
            for frag in tls.iter_fragments():
                i = frag.ifrag
                tbl[i, itls + 3] = BISO[i]

        open(self.txt_path, "w").write(str(tbl))

        ## Gnuplot Script
        script = _BMEAN_PLOT_TEMPLATE
        script = script.replace("<xrng1>", self.cpartition.first_frag_id())
        script = script.replace("<xrng2>", self.cpartition.last_frag_id())
        script = script.replace("<title>", "Observed and TLS Calculated Mean B Factor Per Residue")

        ## line style
        line_titles = []
        ls = 0

        ls += 1
        script += 'set style line %d lc rgb "#000000" lw 3\n' % (ls)
        
        for tls in self.cpartition.iter_tls_segments():
            ls += 1
            script += 'set style line %d lc rgb "%s" lw 3\n' % (ls, tls.color.rgbs)

            if self.tls_group_titles:
                title = 'title "%s TLS"' % (tls.display_label())
            else:
                title = 'notitle'
                
            line_titles.append(title)

        ## plot list
        plist = []
        ls = 0

        ## first the observed mean bfactors
        ls += 1
        x = '"%s" using 1:%d title "Observed" ls %d with lines' % (self.txt_path, 2+ls, ls)
        plist.append(x)

        ## second the TLS calculated bfactors
        for i in xrange(self.cpartition.num_tls_segments()):
            ls += 1
            x = '"%s" using 1:%d %s ls %d with lines' % (self.txt_path, 2+ls, line_titles[i], ls)
            plist.append(x)

        script += "plot " + ",\\\n    ".join(plist) + "\n"
        
        return script


_RMSD_PLOT_TEMPLATE = """\
set xlabel "Residue"
set xrange [<xrng1>:<xrng2>]
set ylabel "RMSD B"
set format y "%5.2f"
set title "<title>"
set datafile missing "?"
"""

class RMSDPlot(GNUPlot):
    def __init__(self, chain, cpartition, **args):
        GNUPlot.__init__(self, **args)
        self.chain = chain
        self.cpartition = cpartition
        self.output_png()

    def write_data_file(self):
        nrows = len(self.cpartition.chain)
        ncols = self.cpartition.num_tls_segments() + 2

        tbl = table.StringTable(nrows, ncols, "?")
        frag_id_iter = itertools.imap(lambda frag: frag.fragment_id, self.cpartition.chain.iter_fragments())
        tbl.set_column(0, 0, frag_id_iter)
                
        if self.cpartition.num_tls_segments() > 1:
            CMTX1 = tls_calcs.calc_residue_mean_rmsd(
                self.chain, self.chain.partition_collection.get_chain_partition(1))
            tbl.set_column(0, 1, CMTX1[0])

        CMTX = tls_calcs.calc_residue_mean_rmsd(self.chain, self.cpartition)
        for itls, tls in enumerate(self.cpartition.iter_tls_segments()):
            for frag in tls.iter_fragments():
                ifrag = frag.ifrag
                tbl[ifrag, itls + 2] = CMTX[itls, ifrag]

        open(self.txt_path, "w").write(str(tbl))

    def make_script(self):
        basename = "%s_CHAIN%s_NTLS%s_RMSD" % (
            self.chain.struct.structure_id,
            self.chain.chain_id,
            self.cpartition.num_tls_segments())

        self.set_basename(basename)

        self.write_data_file()

        ## Gnuplot Script
        script = _RMSD_PLOT_TEMPLATE
        script = script.replace("<xrng1>", self.cpartition.first_frag_id())
        script = script.replace("<xrng2>", self.cpartition.last_frag_id())
        script = script.replace("<title>", "TLS Model RMSD per Residue ")

        ## line style
        line_titles = ["notitle" for x in xrange(self.cpartition.num_tls_segments())]

        if self.cpartition.num_tls_segments() > 1:
            script += 'set style line 1 lc rgb "#000000" lw 3\n'

        ls = 1
        for itls, tls in enumerate(self.cpartition.iter_tls_segments()):
            ls += 1
            script += 'set style line %d lc rgb "%s" lw 3\n' % (ls, tls.color.rgbs)

            title = 'title "%s TLS"' % (tls.display_label())
            line_titles.append(title)

        if self.cpartition.num_tls_segments() > 1:
            ls += 1
            script += 'set style line %d lc rgb "#000000" lw 3\n' % (ls)
            
        ## plot list
        plist = []

        if self.cpartition.num_tls_segments() > 1:
            plist.append('"%s" using 1:2 title "Rigid Chain" ls 1 with lines' % (self.txt_path))

        ls = 1
        for itls in xrange(self.cpartition.num_tls_segments()):
            ls += 1
            col = itls + 3
            x = '"%s" using 1:%d %s ls %d with lines' % (self.txt_path, col, line_titles[itls], ls)
            plist.append(x)

        script += "plot " + ",\\\n    ".join(plist) + "\n"
        
        return script

