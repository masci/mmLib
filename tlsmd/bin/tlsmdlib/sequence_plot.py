## TLS Minimized Domains (TLSMD)
## Copyright 2002-2010 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

## Python Imaging Library imports
import Image
import ImageDraw
import ImageFont

## TLSMD
import misc
import conf

## target pixel width, height, and spacing of sequence
## alignment plots
ALIGN_TARGET_WIDTH = 500
ALIGN_HEIGHT       = 15
ALIGN_SPACING      = 5


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


class TLSSegmentAlignmentPlot(object):
    """Step 1: add all chains, generate unique list of ordered fragment ids,
               and hash all chain+frag_id->tls color
       Step 2: generate graphs
    """
    def __init__(self, **args):
        ## border pixels
        self.border_width = args.get("border_width", 10)

        ## bars are 15 pixels height
        self.pheight    = ALIGN_HEIGHT

        ## spacing pixels between stacked bars
        self.spacing    = ALIGN_SPACING

        ## background color
        self.bg_color   = misc.rgb_f2i((0.90, 0.90, 0.90))

        ## height of the x-axis scale
        self.xscale_height = 20

        self.frag_list      = []
        self.configurations = []

        self.chain_partition_list = []

    def add_tls_segmentation(self, cpartition):
        """Add a TLS optimization to the alignment plot.
        """
        self.chain_partition_list.append(cpartition)
        self.__update_frag_list(cpartition.chain)

    def __update_frag_list(self, chain):
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
        ## Force an ignore of deprecated functions
        import warnings
        warnings.filterwarnings("ignore")

        if len(self.frag_list) == 0 or len(self.chain_partition_list) == 0:
            return False

        nfrag = len(self.frag_list)
        target_width = 500
        fw = int(round(float(ALIGN_TARGET_WIDTH) / nfrag))
        one_frag_width = max(2, fw)

        ## calculate with pixel width/fragment
        ## adjust the width of the graph as necessary
        pheight = self.pheight
        img_width = (2 * self.border_width) + (one_frag_width * len(self.frag_list))

        ## calculate the total height of the image
        num_plots = len(self.chain_partition_list)
        img_height = (2 * self.border_width) + (pheight * num_plots) + (self.spacing * (num_plots-1)) + self.xscale_height

        ## create new image and 2D drawing object
        assert img_width > 0 and img_height > 0
        image = Image.new("RGBA", (img_width, img_height), self.bg_color)
        idraw = ImageDraw.Draw(image)

        idraw.setfill(True)
        idraw.setfont(ImageFont.truetype(conf.GNUPLOT_FONT_PATH, 12))

        x_start = self.border_width

        ## draw plots
        for i in xrange(num_plots):
            cpartition = self.chain_partition_list[i]
            yo = self.border_width + (i * pheight) + (i * self.spacing)
            self.__plot_segmentation(idraw, img_width - 2*self.border_width, one_frag_width, (x_start, yo), cpartition)

        ## draw labels draw some ruler lines
        x_start = self.border_width
        y_xscale = self.border_width + (num_plots * pheight) + ((num_plots+1) * self.spacing)

        y_ruler_top = self.border_width - 1
        y_ruler_bottom = self.border_width + (num_plots * (pheight + self.spacing))

        ## extra tick marks in the plot
        extra_tics = []
        extra_tics.append(0)
        extra_tics.append(len(self.frag_list)-1)

        idraw.setink((0,0,0))
        for i, fid in enumerate(self.frag_list):
            try:
                seq_res = int(str(fid))
            except ValueError:
                continue

            if seq_res%100==0:
                w1 = -1
                w2 = 1
            elif seq_res%50==0:
                w1 = 0
                w2 = 1
            elif seq_res%25==0:
                w1 = 0
                w2 = 0
            elif seq_res%5==0:
                xtic = x_start + i*one_frag_width + (one_frag_width/2)
                idraw.rectangle((xtic, y_ruler_bottom - self.spacing, xtic, y_ruler_bottom - 1))
                continue
            elif i in extra_tics:
                w1 = 0
                w2 = 0
            else:
                continue

            ## draw ruler tick
            xtic = x_start + i*one_frag_width + (one_frag_width/2)
            x1 = xtic + w1
            x2 = xtic + w2
            idraw.rectangle((x1, y_ruler_top, x2, y_ruler_bottom))

            ## draw text label
            sw, sh = idraw.textsize(str(seq_res))
            idraw.text((xtic - (sw/2), y_xscale), str(seq_res))

        ## save image and return
        image.save(path, "png")
        return True

    def __plot_segmentation(self, idraw, pwidth, fwidth, offset, cpartition):
        pheight = self.pheight
        nfrag   = len(self.frag_list)

        ## x/y offsets
        xo, yo = offset

        ## draw a gray background for the length of the chain;
        ## the colored TLS segments will be drawn on top of this
        fid1 = FragmentID(cpartition.chain[0].fragment_id)
        fid2 = FragmentID(cpartition.chain[-1].fragment_id)

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
        for tls in cpartition.iter_tls_segments():
            for frag_id1, frag_id2 in tls.iter_segment_ranges():
                fid1 = FragmentID(frag_id1)
                fid2 = FragmentID(frag_id2)

                i1 = self.frag_list.index(fid1)
                i2 = self.frag_list.index(fid2)

                x1 = i1       * fwidth
                x2 = (i2 + 1) * fwidth

                idraw.setink(tls.color.rgbi)
                idraw.rectangle((x1+xo, yo, x2+xo, pheight+yo))
