## TLS Minimized Domains (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import popen2
import string

## mmLib
from mmLib.Colors         import *
## tlsmdlib
from misc                 import *


class GNUPlot(object):
    """Provides useful methods for subclasses which need to run gnuplot.
    """
    def __init__(self, **args):
        self.gnuplot_path = args.get("gnuplot_path", "gnuplot")
        self.font_path    = args.get("font_path", GNUPLOT_FONT)
        self.font_size    = args.get("fontsize", 10)
        self.width        = args.get("width", 640)
        self.height       = args.get("height", 480)
        self.plot_path    = None
        self.png_path     = None

    def make_script(self):
        """Implement me.
        """
        pass

    def output_png(self):
        """Runs gnuplot.  Expects self.plot_path and self.png_path to be set.
        """
        script = self.make_script()
        
        ## if a basename is given, then write the GnuPlot script as a file
        print "GNUPLot: Saving %s" % (self.png_path)

        ## set output size
        x  = ''
        x += 'set term png font "%s" %d size %d,%d enhanced\n' % (self.font_path, self.font_size, self.width, self.height)
        x += 'set output "%s"\n' % (self.png_path)
        
        script = x + script

        ## optionally write a gnuplot script
        if self.plot_path!=None:
            open(self.plot_path, "w").write(script)
        
        ## run gnuplot
        stdout, stdin, stderr = popen2.popen3((self.gnuplot_path, ), 32768)
        stdin.write(script)
        stdin.close()
        stdout.close()
        stderr.close()



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
    def __init__(self, chainopt, **args):
        GNUPlot.__init__(self, **args)
        self.chainopt = chainopt
        self.output_png()

    def make_script(self):
        chainopt = self.chainopt

        ## generate data and png paths
        basename = "%s_CHAIN%s_RESID" % (chainopt["struct_id"] , chainopt["chain_id"])
        self.txt_path = "%s.txt" % (basename)
        self.png_path = "%s.png" % (basename)
        self.plot_path = "%s.plot" % (basename)

        fil = open(self.txt_path, "w")
        for h, tlsopt in chainopt["ntls_list"]:
            fil.write("%10d %f\n" % (h, tlsopt.residual))
        fil.close()

        ## modify script template
        script = _LSQR_VS_TLS_SEGMENTS_TEMPLATE
        script = script.replace("<nparts>", str(NPARTS))
        script = script.replace("<txtfile>", self.txt_path)
        script = script.replace("<title>", "Least Squares Residual vs. Number of TLS Segments for %s Chain %s " % (
            chainopt["struct_id"], chainopt["chain_id"]))

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
    def __init__(self, chainopt_list, **args):
        GNUPlot.__init__(self, **args)
        self.chainopt_list = chainopt_list
        self.output_png()

    def make_script(self):
        chainopt_list = self.chainopt_list
        
        struct_id = chainopt_list[0]["struct_id"]
        
        ## generate data and png paths
        basename = "%s_RESID" % (struct_id)
        self.png_path = "%s.png" % (basename)
        self.plot_path = "%s.plot" % (basename)

        ## prepare gnuplot script
        script = _LSQR_VS_TLS_SEGMENTS_ALL_CHAINS_TEMPLATE
        script = script.replace("<nparts>", str(NPARTS))
        script = script.replace("<title>", "Least Squares Residual vs. Number of TLS Segments of %s" % (struct_id))

        ## re-use the data files of LSQRvsNTLS from the individual
        ## graphs; to do this the filenames have to be re-constructed
        plist = []
        for chainopt in chainopt_list:
            chain_id = chainopt["chain_id"]
            filename = "%s_CHAIN%s_RESID.txt" % (struct_id, chain_id)
            x = '"%s" using 1:2 title "Chain %s" lw 3 with linespoints' % (filename, chain_id)
            plist.append(x)
            
        script += "plot " + string.join(plist, ",\\\n\t") + "\n"
        return script


_TRANSLATION_ANALYSIS_TEMPLATE = """\
set xlabel "Residue"
set xrange [<xrng1>:<xrng2>]
set ylabel "Angstroms Displacement"
set format y "%4.2f"
set title "<title>"
"""

class TranslationAnalysis(GNUPlot):
    def __init__(self, chainopt, tlsopt, **args):
        GNUPlot.__init__(self, **args)
        self.chainopt = chainopt
        self.tlsopt = tlsopt
        self.output_png()

    def make_script(self):
        chainopt = self.chainopt
        tlsopt = self.tlsopt
        
        basename = "%s_CHAIN%s_NTLS%s_TRANSLATION" % (chainopt["struct_id"], chainopt["chain_id"], tlsopt.ntls)
        self.png_path = "%s.png" % (basename)
        self.txt_path = "%s.txt" % (basename)

        self.write_data_file(tlsopt)

        script = _TRANSLATION_ANALYSIS_TEMPLATE
        script = script.replace("<xrng1>", tlsopt.tls_list[0]["frag_id1"])
        script = script.replace("<xrng2>", tlsopt.tls_list[-1]["frag_id2"])
        script = script.replace("<title>", "Translation Displacement Analysis of Atoms for %d TLS Groups" % (tlsopt.ntls))

        ## line style
        ls = 0
        for tls in tlsopt.tls_list:
            ls += 1
            script += 'set style line %d lc rgb "%s" lw 3\n' % (ls, tls["color"]["rgbs"])

        ## plot list
        plist = []
        ls = 0
        for itls in range(len(tlsopt.tls_list)):
            ls += 1
            
            for n in (0,1,2):
                col = 2 + 3*itls + n
                x = '"%s" using 1:%d smooth bezier notitle ls %d with lines' % (self.txt_path, col , ls)
                plist.append(x)

        script += "plot " + string.join(plist, ",\\\n\t") + "\n"
           
        return script

    def write_data_file(self, tlsopt):
        """Generate the data file and return the filename.
        """
        fil = open(self.txt_path, "w")

        ncols = 1 + 3*len(tlsopt.tls_list)

        for itls in range(len(tlsopt.tls_list)):
            tls = tlsopt.tls_list[itls]

            tls_group = tls["tls_group"]
            tls_info  = tls["model_tls_info"]
            O         = tls_info["COR"]

            for atm in tls_group:
                if atm.name!="CA": continue

                try:
                    ifrag = int(str(atm.fragment_id))
                except ValueError:
                    continue
                
                cols = ["?" for x in range(ncols)]
                cols[0] = str(ifrag)

                ## determine Tr translational eigenvalues
                t1 = GAUSS3C[ADP_PROB] * tls_info["Tr1_rmsd"]
                t2 = GAUSS3C[ADP_PROB] * tls_info["Tr2_rmsd"]
                t3 = GAUSS3C[ADP_PROB] * tls_info["Tr3_rmsd"]

                if t1>0.0: cols[1 + 3*itls] = "%6.4f" % (t1)
                if t2>0.0: cols[2 + 3*itls] = "%6.4f" % (t2)
                if t3>0.0: cols[3 + 3*itls] = "%6.4f" % (t3)
                
                fil.write(string.join(cols) + "\n")

        fil.close()


_LIBRATION_ANALYSIS_TEMPLATE = """\
set xlabel "Residue"
set xrange [<xrng1>:<xrng2>]
set ylabel "Angstroms Displacement"
set format y "%4.2f"
set title "<title>"
set datafile missing "?"
"""

class LibrationAnalysis(GNUPlot):
    def __init__(self, chainopt, tlsopt, **args):
        GNUPlot.__init__(self, **args)
        self.chainopt = chainopt
        self.tlsopt = tlsopt
        self.output_png()

    def make_script(self):
        chainopt = self.chainopt
        tlsopt = self.tlsopt
        
        basename = "%s_CHAIN%s_NTLS%s_LIBRATION" % (chainopt["struct_id"], chainopt["chain_id"], tlsopt.ntls)
        self.png_path = "%s.png" % (basename)
        self.txt_path = "%s.txt" % (basename)

        self.write_data_file(tlsopt)

        script = _LIBRATION_ANALYSIS_TEMPLATE
        script = script.replace("<xrng1>", tlsopt.tls_list[0]["frag_id1"])
        script = script.replace("<xrng2>", tlsopt.tls_list[-1]["frag_id2"])
        script = script.replace("<title>", "Screw displacement analysis of backbone atoms using %d TLS Groups" % (tlsopt.ntls))

        ## line style
        ls = 0
        for tls in tlsopt.tls_list:
            ls += 1
            script += 'set style line %d lc rgb "%s" lw 3\n' % (ls, tls["color"]["rgbs"])

        ## plot list
        plist = []
        ls = 0
        for itls in range(len(tlsopt.tls_list)):
            ls += 1
            
            for n in (0,1,2):
                col = 2 + 3*itls + n
                x = '"%s" using 1:%d smooth bezier notitle ls %d with lines' % (self.txt_path, col , ls)
                plist.append(x)

        script += "plot " + string.join(plist, ",\\\n\t") + "\n"
           
        return script

    def write_data_file(self, tlsopt):
        """Generate the data file and return the filename.
        """
        fil = open(self.txt_path, "w")

        ncols = 1 + 3*len(tlsopt.tls_list)

        for itls in range(len(tlsopt.tls_list)):
            tls = tlsopt.tls_list[itls]

            tls_group = tls["tls_group"]
            tls_info  = tls["model_tls_info"]
            O         = tls_info["COR"]

            for atm in tls_group:
                if atm.name!="CA": continue

                try:
                    ifrag = int(str(atm.fragment_id))
                except ValueError:
                    continue
                
                cols = ["?" for x in range(ncols)]
                cols[0] = str(ifrag)

                for n, Lx_val, Lx_vec, Lx_rho, Lx_pitch in [
                    (0, "L1_eigen_val", "L1_eigen_vec", "L1_rho", "L1_pitch"),
                    (1, "L2_eigen_val", "L2_eigen_vec", "L2_rho", "L2_pitch"),
                    (2, "L3_eigen_val", "L3_eigen_vec", "L3_rho", "L3_pitch") ]:

                    Lval   = tls_info[Lx_val]
                    Lvec   = tls_info[Lx_vec]
                    Lrho   = tls_info[Lx_rho]
                    Lpitch = tls_info[Lx_pitch]

                    if allclose(Lval, 0.0): continue

                    dvec = calc_LS_displacement(O, Lval, Lvec, Lrho, Lpitch, atm.position, ADP_PROB)
                    d = length(dvec)
                    cols[1 + 3*itls + n] = "%8.3f" % (d)

                fil.write(string.join(cols) + "\n")

        fil.close()


_CA_TLS_DIFFERANCE_TEMPLATE = """\
set xlabel "Residue"
set xrange [<xrng1>:<xrng2>]
set ylabel "B_{obs} - B_{tls}"
set format y "%5.2f"
set title "<title>"
set datafile missing "?"
"""

class CA_TLS_Differance_Plot(GNUPlot):
    def __init__(self, chainopt, tlsopt, **args):
        GNUPlot.__init__(self, **args)
        self.chainopt = chainopt
        self.tlsopt = tlsopt
        self.output_png()

    def make_script(self):
        chainopt = self.chainopt
        tlsopt = self.tlsopt

        basename = "%s_CHAIN%s_NTLS%s_CADIFF" % (chainopt["struct_id"], chainopt["chain_id"], tlsopt.ntls)
        self.png_path = "%s.png" % (basename)
        self.txt_path = "%s.txt" % (basename)

        self.write_data_file(tlsopt)

        script = _CA_TLS_DIFFERANCE_TEMPLATE
        script = script.replace("<xrng1>", tlsopt.tls_list[0]["frag_id1"])
        script = script.replace("<xrng2>", tlsopt.tls_list[-1]["frag_id2"])
        script = script.replace("<title>", "Deviation of Observed CA B Factors from TLS Model for %d Group Partition" % (tlsopt.ntls))

        ## line style
        script += 'set style line 1 lc rgb "#000000" lw 1\n'
        
        ls = 1
        for itls in range(len(tlsopt.tls_list)):
            tls = tlsopt.tls_list[itls]
            ls += 1
            script += 'set style line %d lc rgb "%s" lw 3\n' % (ls, tls["color"]["rgbs"])

        ## plot list
        plist = []
        plist.append("0.0 notitle ls 1 with lines")

        ls = 1
        for itls in range(len(tlsopt.tls_list)):
            tls = tlsopt.tls_list[itls]
            ls += 1
            x = '"%s" using 1:%d notitle ls %d with lines' % (self.txt_path, itls+2, ls)
            plist.append(x)

        script += "plot " + string.join(plist, ",\\\n\t") + "\n"
           
        return script

    def write_data_file(self, tlsopt):
        fil = open(self.txt_path, "w")
        ncols = len(tlsopt.tls_list) + 1
        
        for itls in range(len(tlsopt.tls_list)):
            tls = tlsopt.tls_list[itls]
            tls_group = tls["tls_group"]

            T = tls_group.itls_T
            L = tls_group.itls_L
            S = tls_group.itls_S
            O = tls_group.origin

            for atm in tls_group:
                if atm.name!="CA": continue

                try:
                    ifrag = int(str(atm.fragment_id))
                except ValueError:
                    continue
            
                b_tls = U2B * calc_itls_uiso(T, L, S, atm.position - O)
                bdiff = atm.temp_factor - b_tls

                ltmp = ["?" for x in range(ncols)]
                ltmp[0] = str(ifrag)
                ltmp[itls+1] = "%6.2f" % (bdiff)

                fil.write(string.join(ltmp) + "\n")
        
        fil.close()



_UISO_VS_UTLSISO_HISTOGRAM_TEMPLATE = """\
set xlabel "B_{obs} - B_{tls}"
set ylabel "Number of Atoms"
set style line 1 lc rgb "<rgb>" lw 3
set title "<title>"
plot "<txtfile>" using 1:2 ls 1 notitle with histeps
"""

class UIso_vs_UtlsIso_Histogram(GNUPlot):
    def __init__(self, chainopt, tlsopt, tls, **args):
        GNUPlot.__init__(self, **args)
        self.chainopt = chainopt
        self.tlsopt = tlsopt
        self.tls = tls
        self.output_png()

    def make_script(self):
        chainopt = self.chainopt
        tlsopt = self.tlsopt
        tls = self.tls

        ## generate data and png paths
        basename  = "%s_CHAIN%s_TLS%s_%s_BoBc" % (
            chainopt["struct_id"], chainopt["chain_id"], tls["frag_id1"], tls["frag_id2"])
        
        self.txt_path = "%s.txt" % (basename)
        self.png_path = "%s.png" % (basename)

        ## write out the data file
        tls_group = tls["tls_group"]

        T = tls_group.itls_T
        L = tls_group.itls_L
        S = tls_group.itls_S
        O = tls_group.origin

        ## create a histogram of (Uiso - Utls_iso)
        bdiff_min = 0.0
        bdiff_max = 0.0

        for atm in tls_group:
            b_iso_tls = U2B * calc_itls_uiso(T, L, S, atm.position - O)
            bdiff = atm.temp_factor - b_iso_tls

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
        for atm in tls_group:
            b_iso_tls = U2B * calc_itls_uiso(T, L, S, atm.position - O)
            bdiff = atm.temp_factor - b_iso_tls
            bin = int((bdiff - bdiff_min)/ bin_width)
            bins[bin] += 1

        ## write out the gnuplot input file
        fil = open(self.txt_path, "w")
        fil.write("## Histogram of atoms in the TLS group binned by\n")
        fil.write("## the difference of their isotropic temperature factors\n")
        fil.write("## from the isotropic values predicted from the TLS model.\n")
        fil.write("##\n")
        fil.write("## Structure ----------------: %s\n" % (chainopt["struct_id"]))
        fil.write("## Chain --------------------: %s\n" % (chainopt["chain_id"]))
        fil.write("## Number of TLS Groups -----: %d\n" % (tlsopt.ntls))
        fil.write("## TLS Group Residue Range --: %s-%s\n" % (tls["frag_id1"], tls["frag_id2"]))

        for i in range(len(bins)):
            fil.write("%f %d\n" % (bin_names[i], bins[i]))

        fil.close()

        ## modify script template
        script = _UISO_VS_UTLSISO_HISTOGRAM_TEMPLATE
        script = script.replace("<txtfile>", self.txt_path)

        title = "Histogram of Observed B_{iso} - B_{tls} for TLS Group %s%s-%s%s" % (tls["chain_id"], tls["frag_id1"], tls["chain_id"], tls["frag_id2"])
        script = script.replace("<title>", title)
        script = script.replace("<rgb>", tls["color"]["rgbs"])

        return script


def calc_mean_biso_obs(chainopt):
    """Calculates the mean B value per residue in the chain (as observed in the input structure).
    """
    chain = chainopt["chain"]
    num_res = chain.count_fragments()
    biso = zeros(num_res, Float)

    for frag in chain.iter_fragments():
        n = 0
        b_sum_obs = 0.0

        for atm in frag.iter_all_atoms():
            if atm.include==False: continue

            n += 1
            b_sum_obs += atm.temp_factor

        mean_b_obs = b_sum_obs / n
        biso[frag.ichain] = mean_b_obs

    return biso

def calc_mean_biso_tls(chainopt, tlsopt):
    """Calculated the mean B value per residue in the chain as calculated in the chain optimization.
    """
    chain = chainopt["chain"]
    num_res = chain.count_fragments()
    biso = zeros(num_res, Float)

    for i in range(len(tlsopt.tls_list)):
        tls       = tlsopt.tls_list[i]
        tls_group = tls["tls_group"]
        segment   = tls["segment"]

        T = tls_group.itls_T
        L = tls_group.itls_L
        S = tls_group.itls_S
        O = tls_group.origin

        for frag in segment.iter_fragments():
            n = 0
            b_sum_tls = 0.0

            for atm in frag.iter_all_atoms():
                if atm.include==False: continue

                n += 1
                b_sum_tls += U2B * calc_itls_uiso(T, L, S, atm.position - O)

            biso[frag.ichain] =  b_sum_tls / n

    return biso

def calc_accounted_biso(chainopt, tlsopt):
    chain = chainopt["chain"]

    num_res = chain.count_fragments()
    biso = zeros(num_res, Float)

    for i in range(len(tlsopt.tls_list)):
        tls       = tlsopt.tls_list[i]
        tls_group = tls["tls_group"]
        segment   = tls["segment"]

        T = tls_group.itls_T
        L = tls_group.itls_L
        S = tls_group.itls_S
        O = tls_group.origin

        for frag in segment:

            n = 0
            b_sum_obs = 0.0
            b_sum_tls = 0.0

            for atm in frag.iter_all_atoms():
                if atm.include==False: continue

                n += 1
                b_sum_obs += atm.temp_factor
                b_sum_tls += U2B * calc_itls_uiso(T, L, S, atm.position - O)

            mean_b_obs = b_sum_obs / n
            mean_b_tls = b_sum_tls / n

            ## set the cross prediction matrix
            biso[frag.ichain] = mean_b_obs - mean_b_tls

    return biso
     
     
_BMEAN_PLOT_TEMPLATE = """\
set xlabel "Residue"
set xrange [<xrng1>:<xrng2>]
set ylabel "<B>"
set yrange [0.0]
set format y "%5.2f"
set title "<title>"
set datafile missing "?"
"""

class BMeanPlot(GNUPlot):
    def __init__(self, chainopt, tlsopt, **args):
        GNUPlot.__init__(self, **args)
        self.chainopt = chainopt
        self.tlsopt = tlsopt
        self.tls_group_titles = True
        self.output_png()

    def make_script(self):
        chainopt = self.chainopt
        tlsopt = self.tlsopt
        
        basename = "%s_CHAIN%s_NTLS%s_BMEAN" % (chainopt["struct_id"], chainopt["chain_id"], tlsopt.ntls)
        self.png_path = "%s.png" % (basename)
        self.plot_path = "%s.plot" % (basename)

        ## write data file
        BISO1 = calc_mean_biso_obs(chainopt)
        BISO = calc_mean_biso_tls(chainopt, tlsopt)
        chain = chainopt["chain"]

        data_path = "%s.txt" % (basename)
        fil = open(data_path, "w")

        for itls in range(len(tlsopt.tls_list)):
            tls       = tlsopt.tls_list[itls]
            tls_group = tls["tls_group"]
            segment   = tls["segment"]

            for frag in segment:
                ifrag = frag.ichain

                x = ''
                x += '%5s' % (frag.fragment_id)
                x += '  %5d' % (ifrag)

                x += '  %6.2f' % (BISO1[ifrag])
                for i in range(tlsopt.ntls):
                    if i==itls:
                        x += '  %6.2f' % (BISO[ifrag])
                    else:
                        x += '      ?'
            
                x += '\n'
                fil.write(x)
                
        fil.close()

        ## Gnuplot Script
        script = _BMEAN_PLOT_TEMPLATE
        script = script.replace("<xrng1>", tlsopt.tls_list[0]["frag_id1"])
        script = script.replace("<xrng2>", tlsopt.tls_list[-1]["frag_id2"])
        script = script.replace("<title>", "Observed and TLS Calculated Mean B Factor Per Residue")

        ## line style
        line_titles = []
        ls = 0

        ls += 1
        script += 'set style line %d lc rgb "#000000" lw 3\n' % (ls)
        
        for tls in tlsopt.tls_list:
            ls += 1
            script += 'set style line %d lc rgb "%s" lw 3\n' % (ls, tls["color"]["rgbs"])

            if self.tls_group_titles:
                title = 'title "%s-%s TLS"' % (tls["frag_id1"], tls["frag_id2"])
            else:
                title = 'notitle'
                
            line_titles.append(title)

        ## plot list
        plist = []
        ls = 0

        ## first the observed mean bfactors
        ls += 1
        x = '"%s" using 1:%d title "Observed" ls %d with lines' % (data_path, 2+ls, ls)
        plist.append(x)

        ## second the TLS calculated bfactors
        for i in range(len(tlsopt.tls_list)):
            ls += 1
            x = '"%s" using 1:%d %s ls %d with lines' % (data_path, 2+ls, line_titles[i], ls)
            plist.append(x)

        script += "plot " + string.join(plist, ",\\\n\t") + "\n"
        
        return script


def calc_cross_prediction_matrix_rmsd(chainopt, tlsopt):
    chain = chainopt["chain"]

    num_tls = len(tlsopt.tls_list)
    num_res = chain.count_fragments()

    cmtx = zeros((num_tls, num_res), Float)

    for i in range(len(tlsopt.tls_list)):
        tls = tlsopt.tls_list[i]
        tls_group = tls["tls_group"]

        T = tls_group.itls_T
        L = tls_group.itls_L
        S = tls_group.itls_S
        O = tls_group.origin

        for j in range(len(chain)):
            frag = chain[j]

            ## calculate a atom-normalized rmsd deviation for each residue
            n = 0
            delta2 = 0.0
            for atm in frag.iter_all_atoms():
                if atm.include==False: continue

                n += 1
                b_iso_tls = U2B * calc_itls_uiso(T, L, S, atm.position - O)
                delta = atm.temp_factor - b_iso_tls
                delta2 += delta**2

            msd = delta2 / n
            rmsd = math.sqrt(msd)

            ## set the cross prediction matrix
            cmtx[i,j] = rmsd

    return cmtx
     
_RMSD_PLOT_TEMPLATE = """\
set xlabel "Residue"
set xrange [<xrng1>:<xrng2>]
set ylabel "RMSD B"
set format y "%5.2f"
set title "<title>"
set datafile missing "?"
"""

class RMSDPlot(GNUPlot):
    def __init__(self, chainopt, tlsopt, **args):
        GNUPlot.__init__(self, **args)
        self.chainopt = chainopt
        self.tlsopt = tlsopt
        self.output_png()

    def write_data_file(self, chainopt, tlsopt):
        ncols = len(tlsopt.tls_list) + 2
        
        if tlsopt.ntls>1:
            CMTX1 = calc_cross_prediction_matrix_rmsd(chainopt, chainopt["tlsopt"][1])
        else:
            CMTX1 = None

        CMTX = calc_cross_prediction_matrix_rmsd(chainopt, tlsopt)
        m, n = shape(CMTX)

        chain = chainopt["chain"]

        fil = open(self.txt_path, "w")
        
        for j in range(n):
            frag = chain[j]
            try:
                ifrag = int(str(frag.fragment_id))
            except ValueError:
                continue

            cols = ["?" for x in range(ncols)]
            cols[0] = str(ifrag)

            if CMTX1!=None:
                cols[1] = "%6.2f" % (CMTX1[0,j])

            for itls in range(len(tlsopt.tls_list)):
                tls = tlsopt.tls_list[itls]
                segment = tls["segment"]
                if frag in segment:
                    cols[itls+2] = "%6.2f" % (CMTX[itls,j])
                    break

            fil.write(string.join(cols) + "\n")
                
        fil.close()

    def make_script(self):
        chainopt = self.chainopt
        tlsopt = self.tlsopt
        
        basename = "%s_CHAIN%s_NTLS%s_RMSD" % (chainopt["struct_id"], chainopt["chain_id"], tlsopt.ntls)
        self.png_path = "%s.png" % (basename)
        self.plot_path = "%s.plot" % (basename)
        self.txt_path = "%s.txt" % (basename)

        self.write_data_file(chainopt, tlsopt)

        ## Gnuplot Script
        script = _RMSD_PLOT_TEMPLATE
        script = script.replace("<xrng1>", tlsopt.tls_list[0]["frag_id1"])
        script = script.replace("<xrng2>", tlsopt.tls_list[-1]["frag_id2"])
        script = script.replace("<title>", "TLS Model RMSD per Residue ")

        ## line style
        line_titles = ["notitle" for x in range(len(tlsopt.tls_list))]

        if tlsopt.ntls>1:
            script += 'set style line 1 lc rgb "#000000" lw 3\n'

        ls = 1
        for itls in range(len(tlsopt.tls_list)):
            tls = tlsopt.tls_list[itls]
            
            ls += 1
            script += 'set style line %d lc rgb "%s" lw 3\n' % (ls, tls["color"]["rgbs"])

            title = 'title "%s-%s TLS"' % (tls["frag_id1"], tls["frag_id2"])
            line_titles.append(title)

        if tlsopt.ntls>1:
            ls += 1
            script += 'set style line %d lc rgb "#000000" lw 3\n' % (ls)
            
        ## plot list
        plist = []

        if tlsopt.ntls>1:
            plist.append('"%s" using 1:2 title "Rigid Chain" ls 1 with lines' % (self.txt_path))

        ls = 1
        for itls in range(len(tlsopt.tls_list)):
            ls += 1
            col = itls + 3
            x = '"%s" using 1:%d %s ls %d with lines' % (self.txt_path, col, line_titles[itls], ls)
            plist.append(x)

        script += "plot " + string.join(plist, ",\\\n\t") + "\n"
        
        return script

    


