## Copyright 2002 by mmPython Development Group (see AUTHORS file)
## This code is part of the mmPython distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

PenultimateRotamerMap = {}
PenultimateRotamerMap["ARG"] = []
PenultimateRotamerMap["LYS"] = []
PenultimateRotamerMap["MET"] = []
PenultimateRotamerMap["GLU"] = []
PenultimateRotamerMap["GLN"] = []
PenultimateRotamerMap["ASP"] = []
PenultimateRotamerMap["ASN"] = []
PenultimateRotamerMap["ILE"] = []
PenultimateRotamerMap["LEU"] = []
PenultimateRotamerMap["HIS"] = []
PenultimateRotamerMap["TRP"] = []
PenultimateRotamerMap["TYR"] = []
PenultimateRotamerMap["PHE"] = []
PenultimateRotamerMap["PRO"] = []
PenultimateRotamerMap["THR"] = []
PenultimateRotamerMap["VAL"] = []
PenultimateRotamerMap["SER"] = []
PenultimateRotamerMap["CYS"] = []
PenultimateRotamerMap["DIS"] = []


def angdist(a1, a2):
    if (a1>=0 and a2>=0) or (a1<0 and a2<0):
        return abs(a1 - a2)

    else:
        a1_180 = 180 - abs(a1)
        a2_180 = 180 - abs(a2)

        diff0   = abs(a1 - a2)
        diff180 = a1_180 + a2_180

        return min(diff0, diff180)

def FindBestRotamer(aa_res):
    """Compare the chi1 angle with those of the standard rotamer
    library of the residue and return a list of rotamers which have the
    nearest values."""

    try:
        rotamer_list = PenultimateRotamerMap[aa_res.res_name]
    except KeyError:
        return None, None

    (chi1, chi2, chi3, chi4) = aa_res.calcTorsionChi()

    best_rotamer = None
    min_sum      = None

    for rotamer in rotamer_list:
        sum = 0.0

        if chi1 and rotamer.chi1_comm:
            d1 = angdist(chi1, rotamer.chi1_comm)
            sum += d1*d1

        if chi2 and rotamer.chi2_comm:
            d2 = angdist(chi2, rotamer.chi2_comm)
            sum += d2*d2

        if chi3 and rotamer.chi3_comm:
            d3 = angdist(chi3, rotamer.chi3_comm)
            sum += d3*d3

        if chi4 and rotamer.chi4_comm:
            d4 = angdist(chi4, rotamer.chi4_comm)
            sum += d4*d4

        if not min_sum or sum < min_sum:
            min_sum      = sum
            best_rotamer = rotamer
        elif sum == min_sum:
            pass

##     print "    chi1 = (%s, %s)" % (chi1, best_rotamer.chi1_comm)
##     print "    chi2 = (%s, %s)" % (chi2, best_rotamer.chi2_comm)
##     print "    chi3 = (%s, %s)" % (chi3, best_rotamer.chi3_comm)
##     print "    chi4 = (%s, %s)" % (chi4, best_rotamer.chi4_comm)
    
    return min_sum, best_rotamer


class PenultimateRotamer:
    def __init__(self,
                 name               = None,
                 observations       = None,
                 percent            = None,
                 percent_alpha      = None,
                 percent_beta       = None,
                 percent_other      = None,
                 attri              = None,

                 chi1_mode          = None,
                 chi1_comm          = None,
                 chi1_range         = None,
                 chi1_hwhhfa        = None,
                 chi1_hwhhfc        = None,

                 chi2_mode          = None,
                 chi2_comm          = None,
                 chi2_range         = None,
                 chi2_hwhhfa        = None,
                 chi2_hwhhfc        = None,
                 
                 chi3_mode          = None,
                 chi3_comm          = None,
                 chi3_range         = None,
                 chi3_hwhhfa        = None,
                 chi3_hwhhfc        = None,
                 
                 chi4_mode          = None,
                 chi4_comm          = None,
                 chi4_range         = None,
                 chi4_hwhhfa        = None,
                 chi4_hwhhfc        = None):
                 
        self.name             = name
        self.observations     = observations
        self.percent          = percent
        self.percent_alpha    = percent_alpha
        self.percent_beta     = percent_beta
        self.percent_other    = percent_other
        self.attri            = attri

        self.chi1_mode        = chi1_mode
        self.chi1_comm        = chi1_comm
        self.chi1_range       = chi1_range
        self.chi1_hwhhfa      = chi1_hwhhfa
        self.chi1_hwhhfc      = chi1_hwhhfc

        self.chi2_mode        = chi2_mode
        self.chi2_comm        = chi2_comm
        self.chi2_range       = chi2_range
        self.chi2_hwhhfa      = chi2_hwhhfa
        self.chi2_hwhhfc      = chi2_hwhhfc

        self.chi3_mode        = chi3_mode
        self.chi3_comm        = chi3_comm
        self.chi3_range       = chi3_range
        self.chi3_hwhhfa      = chi3_hwhhfa
        self.chi3_hwhhfc      = chi3_hwhhfc

        self.chi4_mode        = chi4_mode
        self.chi4_comm        = chi4_comm
        self.chi4_range       = chi4_range
        self.chi4_hwhhfa      = chi4_hwhhfa
        self.chi4_hwhhfc      = chi4_hwhhfc

    def __str__(self):
        return "%s: %s" % (self.name, str(self.percent))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "ptp85",
    observations      = 3,
    percent           = 0.010000,
    percent_alpha     = 0.0,
    percent_beta      = 0.01,
    percent_other     = 0.0,
    attri             = "minor",

    chi1_mode         = None,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = None,

    chi2_mode         = None,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = 65,
    chi3_range        = (35, 95),
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = 85,
    chi4_range        = (55, 115),
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "ptp180",
    observations      = 11,
    percent           = 0.010000,
    percent_alpha     = 0.0,
    percent_beta      = 0.02,
    percent_other     = 0.02,
    attri             = "major",

    chi1_mode         = 71,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = 14,

    chi2_mode         = 171,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (7, -16),
    chi2_hwhhfc       = 17,

    chi3_mode         = 65,
    chi3_comm         = 65,
    chi3_range        = (35, 95),
    chi3_hwhhfa       = (5, -5),
    chi3_hwhhfc       = 10,

    chi4_mode         = -161,
    chi4_comm         = -175,
    chi4_range        = (155, -145),
    chi4_hwhhfa       = (8, -8),
    chi4_hwhhfc       = 13))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "ptt85",
    observations      = 16,
    percent           = 0.020000,
    percent_alpha     = 0.01,
    percent_beta      = 0.02,
    percent_other     = 0.02,
    attri             = "major",

    chi1_mode         = 65,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = 13,

    chi2_mode         = -178,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (9, -9),
    chi2_hwhhfc       = 14,

    chi3_mode         = -179,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (7, -7),
    chi3_hwhhfc       = 13,

    chi4_mode         = 88,
    chi4_comm         = 85,
    chi4_range        = (55, 115),
    chi4_hwhhfa       = (11, -11),
    chi4_hwhhfc       = 17))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "ptt180",
    observations      = 16,
    percent           = 0.020000,
    percent_alpha     = 0.01,
    percent_beta      = 0.02,
    percent_other     = 0.02,
    attri             = "major",

    chi1_mode         = 59,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (10, -10),
    chi1_hwhhfc       = None,

    chi2_mode         = 176,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (7, -7),
    chi2_hwhhfc       = None,

    chi3_mode         = -178,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (10, -10),
    chi3_hwhhfc       = None,

    chi4_mode         = -177,
    chi4_comm         = 180,
    chi4_range        = (150, -150),
    chi4_hwhhfa       = (19, -8),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "ptt-85",
    observations      = 15,
    percent           = 0.020000,
    percent_alpha     = 0.01,
    percent_beta      = 0.02,
    percent_other     = 0.02,
    attri             = "major",

    chi1_mode         = 66,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = None,

    chi2_mode         = -176,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = None,

    chi3_mode         = -178,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (7, -7),
    chi3_hwhhfc       = None,

    chi4_mode         = -83,
    chi4_comm         = -85,
    chi4_range        = (-115, -55),
    chi4_hwhhfa       = (9, -9),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "ptm180",
    observations      = 6,
    percent           = 0.010000,
    percent_alpha     = 0.0,
    percent_beta      = 0.01,
    percent_other     = 0.01,
    attri             = "minor",

    chi1_mode         = None,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = None,

    chi2_mode         = None,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = -65,
    chi3_range        = (-95, -35),
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = 175,
    chi4_range        = (145, -155),
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "ptm-85",
    observations      = 5,
    percent           = 0.010000,
    percent_alpha     = 0.0,
    percent_beta      = 0.0,
    percent_other     = 0.01,
    attri             = "minor",

    chi1_mode         = None,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = None,

    chi2_mode         = None,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = -65,
    chi3_range        = (-95, -35),
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = -85,
    chi4_range        = (-115, -55),
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "tpp85",
    observations      = 11,
    percent           = 0.010000,
    percent_alpha     = 0.03,
    percent_beta      = 0.01,
    percent_other     = 0.0,
    attri             = "major",

    chi1_mode         = -178,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = None,

    chi2_mode         = 57,
    chi2_comm         = 65,
    chi2_range        = (35, 95),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = None,

    chi3_mode         = 57,
    chi3_comm         = 65,
    chi3_range        = (35, 95),
    chi3_hwhhfa       = (7, -7),
    chi3_hwhhfc       = None,

    chi4_mode         = 85,
    chi4_comm         = 85,
    chi4_range        = (55, 115),
    chi4_hwhhfa       = (9, -9),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "tpp180",
    observations      = 8,
    percent           = 0.010000,
    percent_alpha     = 0.01,
    percent_beta      = 0.0,
    percent_other     = 0.01,
    attri             = "minor",

    chi1_mode         = None,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = None,

    chi2_mode         = None,
    chi2_comm         = 65,
    chi2_range        = (35, 95),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = 65,
    chi3_range        = (35, 95),
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = -175,
    chi4_range        = (155, -145),
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "tpt85",
    observations      = 20,
    percent           = 0.020000,
    percent_alpha     = 0.03,
    percent_beta      = 0.02,
    percent_other     = 0.02,
    attri             = "major",

    chi1_mode         = 177,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = None,

    chi2_mode         = 64,
    chi2_comm         = 65,
    chi2_range        = (35, 95),
    chi2_hwhhfa       = (7, -7),
    chi2_hwhhfc       = None,

    chi3_mode         = 180,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (9, -9),
    chi3_hwhhfc       = None,

    chi4_mode         = 86,
    chi4_comm         = 85,
    chi4_range        = (55, 115),
    chi4_hwhhfa       = (9, -9),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "tpt180",
    observations      = 15,
    percent           = 0.020000,
    percent_alpha     = 0.03,
    percent_beta      = 0.01,
    percent_other     = 0.01,
    attri             = "major",

    chi1_mode         = 179,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = None,

    chi2_mode         = 60,
    chi2_comm         = 65,
    chi2_range        = (35, 95),
    chi2_hwhhfa       = (8, -15),
    chi2_hwhhfc       = None,

    chi3_mode         = 178,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (8, -8),
    chi3_hwhhfc       = None,

    chi4_mode         = 163,
    chi4_comm         = 180,
    chi4_range        = (150, -150),
    chi4_hwhhfa       = (11, -11),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "ttp85",
    observations      = 33,
    percent           = 0.040000,
    percent_alpha     = 0.05,
    percent_beta      = 0.03,
    percent_other     = 0.03,
    attri             = "major",

    chi1_mode         = -179,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = None,

    chi2_mode         = 177,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (12, -12),
    chi2_hwhhfc       = None,

    chi3_mode         = 65,
    chi3_comm         = 65,
    chi3_range        = (35, 95),
    chi3_hwhhfa       = (8, -8),
    chi3_hwhhfc       = None,

    chi4_mode         = 83,
    chi4_comm         = 85,
    chi4_range        = (55, 115),
    chi4_hwhhfa       = (9, -9),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "ttp180",
    observations      = 25,
    percent           = 0.030000,
    percent_alpha     = 0.05,
    percent_beta      = 0.03,
    percent_other     = 0.01,
    attri             = "major",

    chi1_mode         = -178,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = None,

    chi2_mode         = -178,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (8, -14),
    chi2_hwhhfc       = None,

    chi3_mode         = 65,
    chi3_comm         = 65,
    chi3_range        = (35, 95),
    chi3_hwhhfa       = (8, -8),
    chi3_hwhhfc       = None,

    chi4_mode         = -162,
    chi4_comm         = -175,
    chi4_range        = (155, -145),
    chi4_hwhhfa       = (30, -12),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "ttp-105",
    observations      = 9,
    percent           = 0.010000,
    percent_alpha     = 0.01,
    percent_beta      = 0.01,
    percent_other     = 0.01,
    attri             = "minor",

    chi1_mode         = None,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = None,

    chi2_mode         = None,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = 65,
    chi3_range        = (35, 95),
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = -105,
    chi4_range        = (-135, -75),
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "ttt85",
    observations      = 19,
    percent           = 0.020000,
    percent_alpha     = 0.02,
    percent_beta      = 0.02,
    percent_other     = 0.02,
    attri             = "major",

    chi1_mode         = -175,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = None,

    chi2_mode         = 176,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (9, -9),
    chi2_hwhhfc       = None,

    chi3_mode         = 179,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (7, -7),
    chi3_hwhhfc       = None,

    chi4_mode         = 83,
    chi4_comm         = 85,
    chi4_range        = (55, 115),
    chi4_hwhhfa       = (9, -9),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "ttt180",
    observations      = 33,
    percent           = 0.040000,
    percent_alpha     = 0.03,
    percent_beta      = 0.07,
    percent_other     = 0.03,
    attri             = "major",

    chi1_mode         = -179,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (10, -10),
    chi1_hwhhfc       = None,

    chi2_mode         = 177,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = None,

    chi3_mode         = -179,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (7, -7),
    chi3_hwhhfc       = None,

    chi4_mode         = 170,
    chi4_comm         = 180,
    chi4_range        = (150, -150),
    chi4_hwhhfa       = (10, -32),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "ttt-85",
    observations      = 26,
    percent           = 0.030000,
    percent_alpha     = 0.03,
    percent_beta      = 0.03,
    percent_other     = 0.02,
    attri             = "major",

    chi1_mode         = -179,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (10, -10),
    chi1_hwhhfc       = None,

    chi2_mode         = 179,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (9, -9),
    chi2_hwhhfc       = None,

    chi3_mode         = 180,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (9, -9),
    chi3_hwhhfc       = None,

    chi4_mode         = -86,
    chi4_comm         = -85,
    chi4_range        = (-115, -55),
    chi4_hwhhfa       = (10, -10),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "ttm105",
    observations      = 10,
    percent           = 0.010000,
    percent_alpha     = 0.02,
    percent_beta      = 0.01,
    percent_other     = 0.0,
    attri             = "major",

    chi1_mode         = -178,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (10, -10),
    chi1_hwhhfc       = None,

    chi2_mode         = 170,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (8, -13),
    chi2_hwhhfc       = None,

    chi3_mode         = -66,
    chi3_comm         = -65,
    chi3_range        = (-95, -35),
    chi3_hwhhfa       = (9, -9),
    chi3_hwhhfc       = None,

    chi4_mode         = 107,
    chi4_comm         = 105,
    chi4_range        = (75, 135),
    chi4_hwhhfa       = (9, -9),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "ttm180",
    observations      = 13,
    percent           = 0.010000,
    percent_alpha     = 0.0,
    percent_beta      = 0.04,
    percent_other     = 0.01,
    attri             = "major",

    chi1_mode         = 180,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = None,

    chi2_mode         = -178,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (7, -7),
    chi2_hwhhfc       = None,

    chi3_mode         = -67,
    chi3_comm         = -65,
    chi3_range        = (-95, -35),
    chi3_hwhhfa       = (6, -6),
    chi3_hwhhfc       = None,

    chi4_mode         = 176,
    chi4_comm         = 175,
    chi4_range        = (145, -155),
    chi4_hwhhfa       = (10, -10),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "ttm-85",
    observations      = 28,
    percent           = 0.030000,
    percent_alpha     = 0.03,
    percent_beta      = 0.03,
    percent_other     = 0.03,
    attri             = "major",

    chi1_mode         = -175,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = None,

    chi2_mode         = -178,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (11, -11),
    chi2_hwhhfc       = None,

    chi3_mode         = -65,
    chi3_comm         = -65,
    chi3_range        = (-95, -35),
    chi3_hwhhfa       = (10, -10),
    chi3_hwhhfc       = None,

    chi4_mode         = -84,
    chi4_comm         = -85,
    chi4_range        = (-115, -55),
    chi4_hwhhfa       = (8, -8),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "mtp85",
    observations      = 22,
    percent           = 0.020000,
    percent_alpha     = 0.02,
    percent_beta      = 0.03,
    percent_other     = 0.02,
    attri             = "major",

    chi1_mode         = -69,
    chi1_comm         = -67,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = None,

    chi2_mode         = 177,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (12, -12),
    chi2_hwhhfc       = None,

    chi3_mode         = 64,
    chi3_comm         = 65,
    chi3_range        = (35, 95),
    chi3_hwhhfa       = (8, -8),
    chi3_hwhhfc       = None,

    chi4_mode         = 84,
    chi4_comm         = 85,
    chi4_range        = (55, 115),
    chi4_hwhhfa       = (8, -8),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "mtp180",
    observations      = 45,
    percent           = 0.050000,
    percent_alpha     = 0.04,
    percent_beta      = 0.03,
    percent_other     = 0.06,
    attri             = "major",

    chi1_mode         = -65,
    chi1_comm         = -67,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = None,

    chi2_mode         = 176,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (14, -14),
    chi2_hwhhfc       = None,

    chi3_mode         = 64,
    chi3_comm         = 65,
    chi3_range        = (35, 95),
    chi3_hwhhfa       = (8, -8),
    chi3_hwhhfc       = None,

    chi4_mode         = -174,
    chi4_comm         = -175,
    chi4_range        = (155, -145),
    chi4_hwhhfa       = (14, -14),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "mtp-105",
    observations      = 7,
    percent           = 0.010000,
    percent_alpha     = 0.0,
    percent_beta      = 0.02,
    percent_other     = 0.01,
    attri             = "major",

    chi1_mode         = -62,
    chi1_comm         = -67,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (6, -6),
    chi1_hwhhfc       = None,

    chi2_mode         = 179,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (10, -10),
    chi2_hwhhfc       = None,

    chi3_mode         = 67,
    chi3_comm         = 65,
    chi3_range        = (35, 95),
    chi3_hwhhfa       = (8, -8),
    chi3_hwhhfc       = None,

    chi4_mode         = -113,
    chi4_comm         = -105,
    chi4_range        = (-135, -75),
    chi4_hwhhfa       = (6, -14),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "mtt85",
    observations      = 34,
    percent           = 0.040000,
    percent_alpha     = 0.04,
    percent_beta      = 0.04,
    percent_other     = 0.03,
    attri             = "major",

    chi1_mode         = -67,
    chi1_comm         = -67,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = None,

    chi2_mode         = 178,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (18, -10),
    chi2_hwhhfc       = None,

    chi3_mode         = 179,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (8, -8),
    chi3_hwhhfc       = None,

    chi4_mode         = 83,
    chi4_comm         = 85,
    chi4_range        = (55, 115),
    chi4_hwhhfa       = (10, -18),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "mtt180",
    observations      = 89,
    percent           = 0.090000,
    percent_alpha     = 0.09,
    percent_beta      = 0.05,
    percent_other     = 0.12,
    attri             = "major",

    chi1_mode         = -67,
    chi1_comm         = -67,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = None,

    chi2_mode         = -178,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = None,

    chi3_mode         = -177,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (8, -8),
    chi3_hwhhfc       = None,

    chi4_mode         = 174,
    chi4_comm         = 180,
    chi4_range        = (150, -150),
    chi4_hwhhfa       = (16, -16),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "mtt-85",
    observations      = 53,
    percent           = 0.060000,
    percent_alpha     = 0.04,
    percent_beta      = 0.07,
    percent_other     = 0.06,
    attri             = "major",

    chi1_mode         = -66,
    chi1_comm         = -67,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = None,

    chi2_mode         = -177,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = None,

    chi3_mode         = -179,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (8, -8),
    chi3_hwhhfc       = None,

    chi4_mode         = -83,
    chi4_comm         = -85,
    chi4_range        = (-115, -55),
    chi4_hwhhfa       = (8, -8),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "mtm105",
    observations      = 15,
    percent           = 0.020000,
    percent_alpha     = 0.01,
    percent_beta      = 0.01,
    percent_other     = 0.02,
    attri             = "major",

    chi1_mode         = -68,
    chi1_comm         = -67,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = None,

    chi2_mode         = -179,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = None,

    chi3_mode         = -65,
    chi3_comm         = -65,
    chi3_range        = (-95, -35),
    chi3_hwhhfa       = (8, -8),
    chi3_hwhhfc       = None,

    chi4_mode         = 103,
    chi4_comm         = 105,
    chi4_range        = (75, 135),
    chi4_hwhhfa       = (10, -10),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "mtm180",
    observations      = 48,
    percent           = 0.050000,
    percent_alpha     = 0.01,
    percent_beta      = 0.04,
    percent_other     = 0.08,
    attri             = "major",

    chi1_mode         = -68,
    chi1_comm         = -67,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = None,

    chi2_mode         = 173,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (12, -12),
    chi2_hwhhfc       = None,

    chi3_mode         = -64,
    chi3_comm         = -65,
    chi3_range        = (-95, -35),
    chi3_hwhhfa       = (8, -8),
    chi3_hwhhfc       = None,

    chi4_mode         = 180,
    chi4_comm         = 175,
    chi4_range        = (145, -155),
    chi4_hwhhfa       = (34, -12),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "mtm-85",
    observations      = 54,
    percent           = 0.060000,
    percent_alpha     = 0.13,
    percent_beta      = 0.02,
    percent_other     = 0.03,
    attri             = "major",

    chi1_mode         = -69,
    chi1_comm         = -67,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = None,

    chi2_mode         = -167,
    chi2_comm         = -167,
    chi2_range        = (165, -135),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = None,

    chi3_mode         = -63,
    chi3_comm         = -65,
    chi3_range        = (-95, -35),
    chi3_hwhhfa       = (8, -8),
    chi3_hwhhfc       = None,

    chi4_mode         = -86,
    chi4_comm         = -85,
    chi4_range        = (-115, -55),
    chi4_hwhhfa       = (8, -8),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "mmt85",
    observations      = 7,
    percent           = 0.010000,
    percent_alpha     = 0.01,
    percent_beta      = 0.01,
    percent_other     = 0.01,
    attri             = "minor",

    chi1_mode         = None,
    chi1_comm         = -62,
    chi1_range        = (-90, -30),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = None,

    chi2_mode         = None,
    chi2_comm         = -68,
    chi2_range        = (-100, -40),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = 85,
    chi4_range        = (55, 115),
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "mmt180",
    observations      = 18,
    percent           = 0.020000,
    percent_alpha     = 0.01,
    percent_beta      = 0.03,
    percent_other     = 0.02,
    attri             = "major",

    chi1_mode         = -63,
    chi1_comm         = -62,
    chi1_range        = (-90, -30),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = None,

    chi2_mode         = -66,
    chi2_comm         = -68,
    chi2_range        = (-100, -40),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = None,

    chi3_mode         = -179,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (6, -6),
    chi3_hwhhfc       = None,

    chi4_mode         = -168,
    chi4_comm         = 180,
    chi4_range        = (150, -150),
    chi4_hwhhfa       = (30, -18),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "mmt-85",
    observations      = 22,
    percent           = 0.020000,
    percent_alpha     = 0.0,
    percent_beta      = 0.04,
    percent_other     = 0.03,
    attri             = "major",

    chi1_mode         = -60,
    chi1_comm         = -62,
    chi1_range        = (-90, -30),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = None,

    chi2_mode         = -72,
    chi2_comm         = -68,
    chi2_range        = (-100, -40),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = None,

    chi3_mode         = -178,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (10, -10),
    chi3_hwhhfc       = None,

    chi4_mode         = -92,
    chi4_comm         = -85,
    chi4_range        = (-115, -55),
    chi4_hwhhfa       = (8, -8),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "mmm180",
    observations      = 11,
    percent           = 0.010000,
    percent_alpha     = 0.0,
    percent_beta      = 0.02,
    percent_other     = 0.02,
    attri             = "major",

    chi1_mode         = -64,
    chi1_comm         = -62,
    chi1_range        = (-90, -30),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = None,

    chi2_mode         = -74,
    chi2_comm         = -68,
    chi2_range        = (-100, -40),
    chi2_hwhhfa       = (10, -10),
    chi2_hwhhfc       = None,

    chi3_mode         = -67,
    chi3_comm         = -65,
    chi3_range        = (-95, -35),
    chi3_hwhhfa       = (6, -6),
    chi3_hwhhfc       = None,

    chi4_mode         = 172,
    chi4_comm         = 175,
    chi4_range        = (145, -155),
    chi4_hwhhfa       = (8, -8),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ARG"].append(PenultimateRotamer(
    name              = "mmm-85",
    observations      = 22,
    percent           = 0.020000,
    percent_alpha     = 0.02,
    percent_beta      = 0.03,
    percent_other     = 0.03,
    attri             = "major",

    chi1_mode         = -62,
    chi1_comm         = -62,
    chi1_range        = (-90, -30),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = None,

    chi2_mode         = -64,
    chi2_comm         = -68,
    chi2_range        = (-100, -40),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = None,

    chi3_mode         = -61,
    chi3_comm         = -65,
    chi3_range        = (-95, -35),
    chi3_hwhhfa       = (10, -10),
    chi3_hwhhfc       = None,

    chi4_mode         = -82,
    chi4_comm         = -85,
    chi4_range        = (-115, -55),
    chi4_hwhhfa       = (8, -8),
    chi4_hwhhfc       = None))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "ptpt",
    observations      = 7,
    percent           = 0.010000,
    percent_alpha     = 0.0,
    percent_beta      = 0.02,
    percent_other     = 0.0,
    attri             = "minor",

    chi1_mode         = None,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = None,

    chi2_mode         = None,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = 68,
    chi3_range        = (40, 100),
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = 180,
    chi4_range        = (150, -150),
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "pttp",
    observations      = 13,
    percent           = 0.010000,
    percent_alpha     = 0.0,
    percent_beta      = 0.01,
    percent_other     = 0.02,
    attri             = "major",

    chi1_mode         = 63,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = 13,

    chi2_mode         = 177,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (9, -9),
    chi2_hwhhfc       = 14,

    chi3_mode         = 172,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (8, -8),
    chi3_hwhhfc       = 14,

    chi4_mode         = 61,
    chi4_comm         = 65,
    chi4_range        = (35, 95),
    chi4_hwhhfa       = (9, -9),
    chi4_hwhhfc       = 11))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "pttt",
    observations      = 29,
    percent           = 0.020000,
    percent_alpha     = 0.0,
    percent_beta      = 0.04,
    percent_other     = 0.03,
    attri             = "major",

    chi1_mode         = 64,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = 13,

    chi2_mode         = 179,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (10, -10),
    chi2_hwhhfc       = 13,

    chi3_mode         = 180,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (9, -9),
    chi3_hwhhfc       = 13,

    chi4_mode         = 180,
    chi4_comm         = 180,
    chi4_range        = (150, -150),
    chi4_hwhhfa       = (8, -8),
    chi4_hwhhfc       = 10))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "pttm",
    observations      = 8,
    percent           = 0.010000,
    percent_alpha     = 0.0,
    percent_beta      = 0.01,
    percent_other     = 0.01,
    attri             = "minor",

    chi1_mode         = None,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = None,

    chi2_mode         = None,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = -65,
    chi4_range        = (-95, -35),
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "ptmt",
    observations      = 5,
    percent           = 0.000000,
    percent_alpha     = 0.0,
    percent_beta      = 0.01,
    percent_other     = 0.0,
    attri             = "minor",

    chi1_mode         = None,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = None,

    chi2_mode         = None,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = 68,
    chi3_range        = (-100, -40),
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = 180,
    chi4_range        = (150, -150),
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "tptp",
    observations      = 11,
    percent           = 0.010000,
    percent_alpha     = 0.01,
    percent_beta      = 0.01,
    percent_other     = 0.01,
    attri             = "major",

    chi1_mode         = 179,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = 13,

    chi2_mode         = 59,
    chi2_comm         = 68,
    chi2_range        = (40, 100),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = 12,

    chi3_mode         = 163,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = 10,

    chi4_mode         = 60,
    chi4_comm         = 65,
    chi4_range        = (35, 95),
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = 11))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "tptt",
    observations      = 32,
    percent           = 0.030000,
    percent_alpha     = 0.05,
    percent_beta      = 0.01,
    percent_other     = 0.02,
    attri             = "major",

    chi1_mode         = 179,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (4, -4),
    chi1_hwhhfc       = 10,

    chi2_mode         = 62,
    chi2_comm         = 68,
    chi2_range        = (40, 100),
    chi2_hwhhfa       = (3, -3),
    chi2_hwhhfc       = 10,

    chi3_mode         = 173,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (1, -1),
    chi3_hwhhfc       = 13,

    chi4_mode         = 171,
    chi4_comm         = 180,
    chi4_range        = (150, -150),
    chi4_hwhhfa       = (1, -12),
    chi4_hwhhfc       = 14))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "tptm",
    observations      = 7,
    percent           = 0.010000,
    percent_alpha     = 0.01,
    percent_beta      = 0.01,
    percent_other     = 0.0,
    attri             = "major",

    chi1_mode         = None,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = 14,

    chi2_mode         = None,
    chi2_comm         = 68,
    chi2_range        = (40, 100),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = 9,

    chi3_mode         = None,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = 12,

    chi4_mode         = None,
    chi4_comm         = -65,
    chi4_range        = (-95, -35),
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = 10))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "ttpp",
    observations      = 12,
    percent           = 0.010000,
    percent_alpha     = 0.01,
    percent_beta      = 0.0,
    percent_other     = 0.01,
    attri             = "minor",

    chi1_mode         = None,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = None,

    chi2_mode         = None,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = 68,
    chi3_range        = (40, 100),
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = 65,
    chi4_range        = (35, 95),
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "ttpt",
    observations      = 25,
    percent           = 0.020000,
    percent_alpha     = 0.02,
    percent_beta      = 0.05,
    percent_other     = 0.01,
    attri             = "major",

    chi1_mode         = -177,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = 14,

    chi2_mode         = 169,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (9, -9),
    chi2_hwhhfc       = 12,

    chi3_mode         = 71,
    chi3_comm         = 68,
    chi3_range        = (40, 100),
    chi3_hwhhfa       = (13, -13),
    chi3_hwhhfc       = 14,

    chi4_mode         = 172,
    chi4_comm         = 180,
    chi4_range        = (150, -150),
    chi4_hwhhfa       = (8, -8),
    chi4_hwhhfc       = 14))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "tttp",
    observations      = 49,
    percent           = 0.040000,
    percent_alpha     = 0.05,
    percent_beta      = 0.05,
    percent_other     = 0.03,
    attri             = "major",

    chi1_mode         = -177,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = 14,

    chi2_mode         = -170,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (9, -9),
    chi2_hwhhfc       = 13,

    chi3_mode         = -177,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (10, -10),
    chi3_hwhhfc       = 12,

    chi4_mode         = 69,
    chi4_comm         = 65,
    chi4_range        = (35, 95),
    chi4_hwhhfa       = (8, -8),
    chi4_hwhhfc       = 12))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "tttt",
    observations      = 162,
    percent           = 0.130000,
    percent_alpha     = 0.17,
    percent_beta      = 0.19,
    percent_other     = 0.1,
    attri             = "major",

    chi1_mode         = -177,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = 13,

    chi2_mode         = 178,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (10, -10),
    chi2_hwhhfc       = 13,

    chi3_mode         = 179,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (9, -9),
    chi3_hwhhfc       = 15,

    chi4_mode         = 180,
    chi4_comm         = 180,
    chi4_range        = (150, -150),
    chi4_hwhhfa       = (7, -7),
    chi4_hwhhfc       = 13))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "tttm",
    observations      = 37,
    percent           = 0.030000,
    percent_alpha     = 0.04,
    percent_beta      = 0.02,
    percent_other     = 0.03,
    attri             = "major",

    chi1_mode         = -176,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = 12,

    chi2_mode         = -175,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (9, -9),
    chi2_hwhhfc       = 13,

    chi3_mode         = -167,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (9, -9),
    chi3_hwhhfc       = 15,

    chi4_mode         = -59,
    chi4_comm         = -65,
    chi4_range        = (-95, -35),
    chi4_hwhhfa       = (8, -8),
    chi4_hwhhfc       = 13))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "ttmt",
    observations      = 20,
    percent           = 0.020000,
    percent_alpha     = 0.02,
    percent_beta      = 0.04,
    percent_other     = 0.01,
    attri             = "major",

    chi1_mode         = -174,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (10, -10),
    chi1_hwhhfc       = 14,

    chi2_mode         = -160,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = 14,

    chi3_mode         = -54,
    chi3_comm         = -68,
    chi3_range        = (-100, -40),
    chi3_hwhhfa       = (7, -7),
    chi3_hwhhfc       = 10,

    chi4_mode         = -174,
    chi4_comm         = 180,
    chi4_range        = (150, -150),
    chi4_hwhhfa       = (9, -9),
    chi4_hwhhfc       = 15))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "ttmm",
    observations      = 5,
    percent           = 0.000000,
    percent_alpha     = 0.01,
    percent_beta      = 0.0,
    percent_other     = 0.0,
    attri             = "minor",

    chi1_mode         = None,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = None,

    chi2_mode         = None,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = -68,
    chi3_range        = (-100, -40),
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = -65,
    chi4_range        = (-95, -35),
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "mptt",
    observations      = 4,
    percent           = 0.000000,
    percent_alpha     = 0.0,
    percent_beta      = 0.0,
    percent_other     = 0.01,
    attri             = "minor",

    chi1_mode         = None,
    chi1_comm         = -90,
    chi1_range        = (-120, -60),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = None,

    chi2_mode         = None,
    chi2_comm         = 68,
    chi2_range        = (40, 100),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = 180,
    chi4_range        = (150, -150),
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "mtpp",
    observations      = 12,
    percent           = 0.010000,
    percent_alpha     = 0.01,
    percent_beta      = 0.01,
    percent_other     = 0.01,
    attri             = "major",

    chi1_mode         = -62,
    chi1_comm         = -67,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = 10,

    chi2_mode         = -177,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (10, -10),
    chi2_hwhhfc       = 9,

    chi3_mode         = 71,
    chi3_comm         = 68,
    chi3_range        = (40, 100),
    chi3_hwhhfa       = (8, -8),
    chi3_hwhhfc       = 10,

    chi4_mode         = 68,
    chi4_comm         = 65,
    chi4_range        = (35, 95),
    chi4_hwhhfa       = (6, -6),
    chi4_hwhhfc       = 13))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "mtpt",
    observations      = 38,
    percent           = 0.030000,
    percent_alpha     = 0.04,
    percent_beta      = 0.02,
    percent_other     = 0.03,
    attri             = "major",

    chi1_mode         = -69,
    chi1_comm         = -67,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = 12,

    chi2_mode         = 174,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (10, -10),
    chi2_hwhhfc       = 13,

    chi3_mode         = 65,
    chi3_comm         = 68,
    chi3_range        = (40, 100),
    chi3_hwhhfa       = (10, -10),
    chi3_hwhhfc       = 11,

    chi4_mode         = 175,
    chi4_comm         = 180,
    chi4_range        = (150, -150),
    chi4_hwhhfa       = (11, -11),
    chi4_hwhhfc       = 9))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "mttp",
    observations      = 42,
    percent           = 0.030000,
    percent_alpha     = 0.02,
    percent_beta      = 0.04,
    percent_other     = 0.04,
    attri             = "major",

    chi1_mode         = -67,
    chi1_comm         = -67,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = 13,

    chi2_mode         = -176,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (7, -7),
    chi2_hwhhfc       = 13,

    chi3_mode         = 175,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (10, -10),
    chi3_hwhhfc       = 14,

    chi4_mode         = 71,
    chi4_comm         = 65,
    chi4_range        = (35, 95),
    chi4_hwhhfa       = (11, -11),
    chi4_hwhhfc       = 14))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "mttt",
    observations      = 244,
    percent           = 0.200000,
    percent_alpha     = 0.23,
    percent_beta      = 0.14,
    percent_other     = 0.21,
    attri             = "major",

    chi1_mode         = -67,
    chi1_comm         = -67,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = 14,

    chi2_mode         = 179,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (9, -9),
    chi2_hwhhfc       = 13,

    chi3_mode         = 179,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (9, -9),
    chi3_hwhhfc       = 12,

    chi4_mode         = 179,
    chi4_comm         = 180,
    chi4_range        = (150, -150),
    chi4_hwhhfa       = (10, -10),
    chi4_hwhhfc       = 14))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "mttm",
    observations      = 56,
    percent           = 0.050000,
    percent_alpha     = 0.03,
    percent_beta      = 0.05,
    percent_other     = 0.06,
    attri             = "major",

    chi1_mode         = -65,
    chi1_comm         = -67,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = 13,

    chi2_mode         = 180,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (9, -9),
    chi2_hwhhfc       = 12,

    chi3_mode         = -178,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (10, -10),
    chi3_hwhhfc       = 13,

    chi4_mode         = -63,
    chi4_comm         = -65,
    chi4_range        = (-95, -35),
    chi4_hwhhfa       = (10, -10),
    chi4_hwhhfc       = 14))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "mtmt",
    observations      = 40,
    percent           = 0.030000,
    percent_alpha     = 0.06,
    percent_beta      = 0.02,
    percent_other     = 0.03,
    attri             = "major",

    chi1_mode         = -68,
    chi1_comm         = -67,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = 12,

    chi2_mode         = -171,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (9, -9),
    chi2_hwhhfc       = 13,

    chi3_mode         = -68,
    chi3_comm         = -68,
    chi3_range        = (-100, -40),
    chi3_hwhhfa       = (9, -9),
    chi3_hwhhfc       = 14,

    chi4_mode         = -175,
    chi4_comm         = 180,
    chi4_range        = (150, -150),
    chi4_hwhhfa       = (9, -9),
    chi4_hwhhfc       = 13))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "mtmm",
    observations      = 12,
    percent           = 0.010000,
    percent_alpha     = 0.0,
    percent_beta      = 0.01,
    percent_other     = 0.01,
    attri             = "major",

    chi1_mode         = -68,
    chi1_comm         = -67,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = 12,

    chi2_mode         = -178,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = 12,

    chi3_mode         = -68,
    chi3_comm         = -68,
    chi3_range        = (-100, -40),
    chi3_hwhhfa       = (9, -9),
    chi3_hwhhfc       = 12,

    chi4_mode         = -63,
    chi4_comm         = -65,
    chi4_range        = (-95, -35),
    chi4_hwhhfa       = (8, -8),
    chi4_hwhhfc       = 11))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "mmtp",
    observations      = 9,
    percent           = 0.010000,
    percent_alpha     = 0.0,
    percent_beta      = 0.0,
    percent_other     = 0.01,
    attri             = "minor",

    chi1_mode         = None,
    chi1_comm         = -62,
    chi1_range        = (-90, -30),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = None,

    chi2_mode         = None,
    chi2_comm         = -68,
    chi2_range        = (-100, -40),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = 65,
    chi4_range        = (35, 95),
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "mmtt",
    observations      = 77,
    percent           = 0.060000,
    percent_alpha     = 0.03,
    percent_beta      = 0.05,
    percent_other     = 0.08,
    attri             = "major",

    chi1_mode         = -61,
    chi1_comm         = -62,
    chi1_range        = (-90, -30),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = 12,

    chi2_mode         = -64,
    chi2_comm         = -68,
    chi2_range        = (-100, -40),
    chi2_hwhhfa       = (9, -9),
    chi2_hwhhfc       = 13,

    chi3_mode         = -177,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (9, -9),
    chi3_hwhhfc       = 13,

    chi4_mode         = -179,
    chi4_comm         = 180,
    chi4_range        = (150, -150),
    chi4_hwhhfa       = (9, -9),
    chi4_hwhhfc       = 13))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "mmtm",
    observations      = 18,
    percent           = 0.010000,
    percent_alpha     = 0.01,
    percent_beta      = 0.01,
    percent_other     = 0.02,
    attri             = "major",

    chi1_mode         = -64,
    chi1_comm         = -62,
    chi1_range        = (-90, -30),
    chi1_hwhhfa       = (10, -10),
    chi1_hwhhfc       = 14,

    chi2_mode         = -68,
    chi2_comm         = -68,
    chi2_range        = (-100, -40),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = 12,

    chi3_mode         = -176,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (8, -8),
    chi3_hwhhfc       = 10,

    chi4_mode         = -71,
    chi4_comm         = -65,
    chi4_range        = (-95, -35),
    chi4_hwhhfa       = (10, -10),
    chi4_hwhhfc       = 15))

PenultimateRotamerMap["LYS"].append(PenultimateRotamer(
    name              = "mmmt",
    observations      = 10,
    percent           = 0.010000,
    percent_alpha     = 0.0,
    percent_beta      = 0.01,
    percent_other     = 0.01,
    attri             = "major",

    chi1_mode         = -59,
    chi1_comm         = -62,
    chi1_range        = (-90, -30),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = 12,

    chi2_mode         = -59,
    chi2_comm         = -68,
    chi2_range        = (-100, -40),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = 13,

    chi3_mode         = -70,
    chi3_comm         = -68,
    chi3_range        = (-100, -40),
    chi3_hwhhfa       = (8, -8),
    chi3_hwhhfc       = 10,

    chi4_mode         = -173,
    chi4_comm         = 180,
    chi4_range        = (150, -150),
    chi4_hwhhfa       = (10, -10),
    chi4_hwhhfc       = 15))

PenultimateRotamerMap["MET"].append(PenultimateRotamer(
    name              = "ptp",
    observations      = 12,
    percent           = 0.020000,
    percent_alpha     = 0.01,
    percent_beta      = 0.03,
    percent_other     = 0.03,
    attri             = "major",

    chi1_mode         = 68,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (11, -5),
    chi1_hwhhfc       = 11,

    chi2_mode         = -167,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (21, -6),
    chi2_hwhhfc       = 17,

    chi3_mode         = 88,
    chi3_comm         = 75,
    chi3_range        = (45, 105),
    chi3_hwhhfa       = (9, -9),
    chi3_hwhhfc       = 12,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["MET"].append(PenultimateRotamer(
    name              = "ptm",
    observations      = 17,
    percent           = 0.030000,
    percent_alpha     = 0.01,
    percent_beta      = 0.06,
    percent_other     = 0.04,
    attri             = "major",

    chi1_mode         = 67,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (5, -5),
    chi1_hwhhfc       = 9,

    chi2_mode         = 174,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (6, -6),
    chi2_hwhhfc       = 10,

    chi3_mode         = -78,
    chi3_comm         = -75,
    chi3_range        = (-105, -45),
    chi3_hwhhfa       = (6, -6),
    chi3_hwhhfc       = 9,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["MET"].append(PenultimateRotamer(
    name              = "tpp",
    observations      = 30,
    percent           = 0.050000,
    percent_alpha     = 0.08,
    percent_beta      = 0.02,
    percent_other     = 0.05,
    attri             = "major",

    chi1_mode         = -177,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = 10,

    chi2_mode         = 66,
    chi2_comm         = 65,
    chi2_range        = (35, 95),
    chi2_hwhhfa       = (14, -9),
    chi2_hwhhfc       = 15,

    chi3_mode         = 75,
    chi3_comm         = 75,
    chi3_range        = (45, 105),
    chi3_hwhhfa       = (14, -10),
    chi3_hwhhfc       = 15,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["MET"].append(PenultimateRotamer(
    name              = "tpt",
    observations      = 9,
    percent           = 0.020000,
    percent_alpha     = 0.01,
    percent_beta      = 0.04,
    percent_other     = 0.01,
    attri             = "major",

    chi1_mode         = 179,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (6, -6),
    chi1_hwhhfc       = 9,

    chi2_mode         = 67,
    chi2_comm         = 65,
    chi2_range        = (35, 95),
    chi2_hwhhfa       = (4, -4),
    chi2_hwhhfc       = 8,

    chi3_mode         = -179,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (6, -6),
    chi3_hwhhfc       = 9,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["MET"].append(PenultimateRotamer(
    name              = "ttp",
    observations      = 28,
    percent           = 0.050000,
    percent_alpha     = 0.07,
    percent_beta      = 0.07,
    percent_other     = 0.02,
    attri             = "major",

    chi1_mode         = 176,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = 10,

    chi2_mode         = 178,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (7, -7),
    chi2_hwhhfc       = 11,

    chi3_mode         = 73,
    chi3_comm         = 75,
    chi3_range        = (45, 105),
    chi3_hwhhfa       = (8, -8),
    chi3_hwhhfc       = 11,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["MET"].append(PenultimateRotamer(
    name              = "ttt",
    observations      = 17,
    percent           = 0.030000,
    percent_alpha     = 0.05,
    percent_beta      = 0.02,
    percent_other     = 0.02,
    attri             = "major",

    chi1_mode         = 180,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (5, -5),
    chi1_hwhhfc       = 9,

    chi2_mode         = 171,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (6, -6),
    chi2_hwhhfc       = 9,

    chi3_mode         = 174,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (12, -18),
    chi3_hwhhfc       = 19,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["MET"].append(PenultimateRotamer(
    name              = "ttm",
    observations      = 36,
    percent           = 0.070000,
    percent_alpha     = 0.03,
    percent_beta      = 0.1,
    percent_other     = 0.08,
    attri             = "major",

    chi1_mode         = -177,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = 10,

    chi2_mode         = 176,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (7, -7),
    chi2_hwhhfc       = 10,

    chi3_mode         = -78,
    chi3_comm         = -75,
    chi3_range        = (-105, -45),
    chi3_hwhhfa       = (10, -10),
    chi3_hwhhfc       = 13,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["MET"].append(PenultimateRotamer(
    name              = "mtp",
    observations      = 92,
    percent           = 0.170000,
    percent_alpha     = 0.22,
    percent_beta      = 0.1,
    percent_other     = 0.17,
    attri             = "major",

    chi1_mode         = -68,
    chi1_comm         = -67,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (6, -6),
    chi1_hwhhfc       = 10,

    chi2_mode         = 177,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = 12,

    chi3_mode         = 72,
    chi3_comm         = 75,
    chi3_range        = (45, 105),
    chi3_hwhhfa       = (10, -10),
    chi3_hwhhfc       = 14,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["MET"].append(PenultimateRotamer(
    name              = "mtt",
    observations      = 43,
    percent           = 0.080000,
    percent_alpha     = 0.09,
    percent_beta      = 0.08,
    percent_other     = 0.07,
    attri             = "major",

    chi1_mode         = -67,
    chi1_comm         = -67,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = 10,

    chi2_mode         = 177,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (10, -10),
    chi2_hwhhfc       = 13,

    chi3_mode         = -178,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (8, -14),
    chi3_hwhhfc       = 15,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["MET"].append(PenultimateRotamer(
    name              = "mtm",
    observations      = 58,
    percent           = 0.110000,
    percent_alpha     = 0.12,
    percent_beta      = 0.11,
    percent_other     = 0.09,
    attri             = "major",

    chi1_mode         = -67,
    chi1_comm         = -67,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = 12,

    chi2_mode         = -177,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = 11,

    chi3_mode         = -76,
    chi3_comm         = -75,
    chi3_range        = (-105, -45),
    chi3_hwhhfa       = (12, -12),
    chi3_hwhhfc       = 16,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["MET"].append(PenultimateRotamer(
    name              = "mmp",
    observations      = 15,
    percent           = 0.030000,
    percent_alpha     = 0.03,
    percent_beta      = 0.01,
    percent_other     = 0.04,
    attri             = "major",

    chi1_mode         = -64,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (6, -6),
    chi1_hwhhfc       = 9,

    chi2_mode         = -63,
    chi2_comm         = -65,
    chi2_range        = (-95, -35),
    chi2_hwhhfa       = (7, -7),
    chi2_hwhhfc       = 10,

    chi3_mode         = 103,
    chi3_comm         = 103,
    chi3_range        = (75, 135),
    chi3_hwhhfa       = (7, -7),
    chi3_hwhhfc       = 10,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["MET"].append(PenultimateRotamer(
    name              = "mmt",
    observations      = 10,
    percent           = 0.020000,
    percent_alpha     = 0.0,
    percent_beta      = 0.02,
    percent_other     = 0.03,
    attri             = "major",

    chi1_mode         = -63,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = 12,

    chi2_mode         = -64,
    chi2_comm         = -65,
    chi2_range        = (-95, -35),
    chi2_hwhhfa       = (11, -11),
    chi2_hwhhfc       = 14,

    chi3_mode         = 180,
    chi3_comm         = 180,
    chi3_range        = (150, -150),
    chi3_hwhhfa       = (15, -15),
    chi3_hwhhfc       = 19,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["MET"].append(PenultimateRotamer(
    name              = "mmm",
    observations      = 105,
    percent           = 0.190000,
    percent_alpha     = 0.21,
    percent_beta      = 0.16,
    percent_other     = 0.19,
    attri             = "major",

    chi1_mode         = -66,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = 11,

    chi2_mode         = -60,
    chi2_comm         = -65,
    chi2_range        = (-95, -35),
    chi2_hwhhfa       = (9, -9),
    chi2_hwhhfc       = 13,

    chi3_mode         = -67,
    chi3_comm         = -70,
    chi3_range        = (-100, -40),
    chi3_hwhhfa       = (15, -10),
    chi3_hwhhfc       = 16,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["GLU"].append(PenultimateRotamer(
    name              = "pt-20 ",
    observations      = 80,
    percent           = 0.050000,
    percent_alpha     = 0.01,
    percent_beta      = 0.09,
    percent_other     = 0.07,
    attri             = "major",

    chi1_mode         = 63,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (10, -10),
    chi1_hwhhfc       = 14,

    chi2_mode         = -175,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (10, -10),
    chi2_hwhhfc       = 13,

    chi3_mode         = -18,
    chi3_comm         = -20,
    chi3_range        = (-90, 90),
    chi3_hwhhfa       = (15, -23),
    chi3_hwhhfc       = 23,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["GLU"].append(PenultimateRotamer(
    name              = "pm0",
    observations      = 32,
    percent           = 0.020000,
    percent_alpha     = 0.0,
    percent_beta      = 0.0,
    percent_other     = 0.04,
    attri             = "major",

    chi1_mode         = 71,
    chi1_comm         = 70,
    chi1_range        = (40, 100),
    chi1_hwhhfa       = (11, -11),
    chi1_hwhhfc       = 14,

    chi2_mode         = -79,
    chi2_comm         = -80,
    chi2_range        = (-110, -50),
    chi2_hwhhfa       = (10, -10),
    chi2_hwhhfc       = 13,

    chi3_mode         = 5,
    chi3_comm         = 0,
    chi3_range        = (-50, 50),
    chi3_hwhhfa       = (13, -13),
    chi3_hwhhfc       = 17,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["GLU"].append(PenultimateRotamer(
    name              = "tp10",
    observations      = 91,
    percent           = 0.060000,
    percent_alpha     = 0.1,
    percent_beta      = 0.02,
    percent_other     = 0.06,
    attri             = "major",

    chi1_mode         = -177,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (10, -10),
    chi1_hwhhfc       = 14,

    chi2_mode         = 65,
    chi2_comm         = 65,
    chi2_range        = (35, 95),
    chi2_hwhhfa       = (10, -10),
    chi2_hwhhfc       = 13,

    chi3_mode         = 13,
    chi3_comm         = 10,
    chi3_range        = (-10, 90),
    chi3_hwhhfa       = (10, -16),
    chi3_hwhhfc       = 17,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["GLU"].append(PenultimateRotamer(
    name              = "tt 0 ",
    observations      = 350,
    percent           = 0.240000,
    percent_alpha     = 0.25,
    percent_beta      = 0.42,
    percent_other     = 0.18,
    attri             = "major",

    chi1_mode         = -177,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (10, -10),
    chi1_hwhhfc       = 14,

    chi2_mode         = 178,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (10, -10),
    chi2_hwhhfc       = 14,

    chi3_mode         = 2,
    chi3_comm         = 0,
    chi3_range        = (-90, 90),
    chi3_hwhhfa       = (25, -25),
    chi3_hwhhfc       = 30,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["GLU"].append(PenultimateRotamer(
    name              = "tm-20",
    observations      = 17,
    percent           = 0.010000,
    percent_alpha     = 0.01,
    percent_beta      = 0.01,
    percent_other     = 0.01,
    attri             = "major",

    chi1_mode         = None,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = 13,

    chi2_mode         = None,
    chi2_comm         = -80,
    chi2_range        = (-115, -55),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = 13,

    chi3_mode         = None,
    chi3_comm         = -25,
    chi3_range        = (-50, 10),
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = 15,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["GLU"].append(PenultimateRotamer(
    name              = "mp0",
    observations      = 88,
    percent           = 0.060000,
    percent_alpha     = 0.0,
    percent_beta      = 0.02,
    percent_other     = 0.1,
    attri             = "major",

    chi1_mode         = -65,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (10, -10),
    chi1_hwhhfc       = 14,

    chi2_mode         = 85,
    chi2_comm         = 85,
    chi2_range        = (55, 115),
    chi2_hwhhfa       = (9, -9),
    chi2_hwhhfc       = 13,

    chi3_mode         = -3,
    chi3_comm         = 0,
    chi3_range        = (-60, 60),
    chi3_hwhhfa       = (16, -25),
    chi3_hwhhfc       = 25,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["GLU"].append(PenultimateRotamer(
    name              = "mt-10 ",
    observations      = 484,
    percent           = 0.330000,
    percent_alpha     = 0.36,
    percent_beta      = 0.29,
    percent_other     = 0.32,
    attri             = "major",

    chi1_mode         = -67,
    chi1_comm         = -67,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = 13,

    chi2_mode         = 177,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (12, -12),
    chi2_hwhhfc       = 16,

    chi3_mode         = -10,
    chi3_comm         = -10,
    chi3_range        = (-90, 90),
    chi3_hwhhfa       = (23, -18),
    chi3_hwhhfc       = 25,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["GLU"].append(PenultimateRotamer(
    name              = "mm-40",
    observations      = 197,
    percent           = 0.130000,
    percent_alpha     = 0.19,
    percent_beta      = 0.07,
    percent_other     = 0.12,
    attri             = "major",

    chi1_mode         = -65,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (10, -10),
    chi1_hwhhfc       = 14,

    chi2_mode         = -58,
    chi2_comm         = -65,
    chi2_range        = (-95, -35),
    chi2_hwhhfa       = (10, -10),
    chi2_hwhhfc       = 14,

    chi3_mode         = -40,
    chi3_comm         = -40,
    chi3_range        = (-90, 30),
    chi3_hwhhfa       = (15, -26),
    chi3_hwhhfc       = 25,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["GLN"].append(PenultimateRotamer(
    name              = "pt 20",
    observations      = 37,
    percent           = 0.040000,
    percent_alpha     = 0.01,
    percent_beta      = 0.05,
    percent_other     = 0.06,
    attri             = "major",

    chi1_mode         = 64,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (6, -6),
    chi1_hwhhfc       = 13,

    chi2_mode         = 180,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (6, -6),
    chi2_hwhhfc       = 14,

    chi3_mode         = 20,
    chi3_comm         = 20,
    chi3_range        = (-90, 90),
    chi3_hwhhfa       = (8, -8),
    chi3_hwhhfc       = 16,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["GLN"].append(PenultimateRotamer(
    name              = "pm0",
    observations      = 15,
    percent           = 0.020000,
    percent_alpha     = 0.0,
    percent_beta      = 0.01,
    percent_other     = 0.03,
    attri             = "minor",

    chi1_mode         = None,
    chi1_comm         = 70,
    chi1_range        = (45, 105),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = None,

    chi2_mode         = None,
    chi2_comm         = -75,
    chi2_range        = (-105, -45),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = 0,
    chi3_range        = (-60, 60),
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["GLN"].append(PenultimateRotamer(
    name              = "tp-100",
    observations      = 14,
    percent           = 0.020000,
    percent_alpha     = 0.04,
    percent_beta      = 0.02,
    percent_other     = 0.0,
    attri             = "minor",

    chi1_mode         = None,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = None,

    chi2_mode         = None,
    chi2_comm         = 65,
    chi2_range        = (35, 95),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = -100,
    chi3_range        = (-150, 0),
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["GLN"].append(PenultimateRotamer(
    name              = "tp60",
    observations      = 78,
    percent           = 0.090000,
    percent_alpha     = 0.13,
    percent_beta      = 0.09,
    percent_other     = 0.07,
    attri             = "major",

    chi1_mode         = -175,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (6, -6),
    chi1_hwhhfc       = 14,

    chi2_mode         = 64,
    chi2_comm         = 65,
    chi2_range        = (35, 95),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = 15,

    chi3_mode         = 60,
    chi3_comm         = 60,
    chi3_range        = (0, 90),
    chi3_hwhhfa       = (11, -21),
    chi3_hwhhfc       = 24,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["GLN"].append(PenultimateRotamer(
    name              = "tt 0",
    observations      = 140,
    percent           = 0.160000,
    percent_alpha     = 0.16,
    percent_beta      = 0.29,
    percent_other     = 0.12,
    attri             = "major",

    chi1_mode         = -174,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (6, -6),
    chi1_hwhhfc       = 14,

    chi2_mode         = 173,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (6, -6),
    chi2_hwhhfc       = 13,

    chi3_mode         = -5,
    chi3_comm         = 0,
    chi3_range        = (-90, 90),
    chi3_hwhhfa       = (52, -11),
    chi3_hwhhfc       = 40,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["GLN"].append(PenultimateRotamer(
    name              = "mp0",
    observations      = 24,
    percent           = 0.030000,
    percent_alpha     = 0.0,
    percent_beta      = 0.01,
    percent_other     = 0.05,
    attri             = "minor",

    chi1_mode         = None,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = None,

    chi2_mode         = None,
    chi2_comm         = 85,
    chi2_range        = (55, 115),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = 0,
    chi3_range        = (-60, 60),
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["GLN"].append(PenultimateRotamer(
    name              = "mt -30",
    observations      = 304,
    percent           = 0.350000,
    percent_alpha     = 0.4,
    percent_beta      = 0.26,
    percent_other     = 0.36,
    attri             = "major",

    chi1_mode         = -67,
    chi1_comm         = -67,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = 16,

    chi2_mode         = 177,
    chi2_comm         = 180,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (7, -7),
    chi2_hwhhfc       = 15,

    chi3_mode         = -25,
    chi3_comm         = -25,
    chi3_range        = (-90, 90),
    chi3_hwhhfa       = (29, -29),
    chi3_hwhhfc       = 37,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["GLN"].append(PenultimateRotamer(
    name              = "mm-40",
    observations      = 127,
    percent           = 0.150000,
    percent_alpha     = 0.12,
    percent_beta      = 0.13,
    percent_other     = 0.17,
    attri             = "major",

    chi1_mode         = -66,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = 16,

    chi2_mode         = -60,
    chi2_comm         = -65,
    chi2_range        = (-95, -35),
    chi2_hwhhfa       = (10, -10),
    chi2_hwhhfc       = 18,

    chi3_mode         = -40,
    chi3_comm         = -40,
    chi3_range        = (-95, 0),
    chi3_hwhhfa       = (14, -19),
    chi3_hwhhfc       = 26,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["GLN"].append(PenultimateRotamer(
    name              = "mm100",
    observations      = 22,
    percent           = 0.030000,
    percent_alpha     = 0.04,
    percent_beta      = 0.01,
    percent_other     = 0.02,
    attri             = "minor",

    chi1_mode         = None,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = None,

    chi2_mode         = None,
    chi2_comm         = -65,
    chi2_range        = (-95, -35),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = 100,
    chi3_range        = (0, 150),
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ASP"].append(PenultimateRotamer(
    name              = "p-10",
    observations      = 203,
    percent           = 0.100000,
    percent_alpha     = 0.01,
    percent_beta      = 0.02,
    percent_other     = 0.13,
    attri             = "major",

    chi1_mode         = 61,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (6, -6),
    chi1_hwhhfc       = 9,

    chi2_mode         = -4,
    chi2_comm         = -10,
    chi2_range        = (-90, 0),
    chi2_hwhhfa       = (13, -21),
    chi2_hwhhfc       = 19,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ASP"].append(PenultimateRotamer(
    name              = "p30",
    observations      = 194,
    percent           = 0.090000,
    percent_alpha     = 0.01,
    percent_beta      = 0.05,
    percent_other     = 0.12,
    attri             = "major",

    chi1_mode         = 65,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = 8,

    chi2_mode         = 9,
    chi2_comm         = 30,
    chi2_range        = (0, 90),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = 14,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ASP"].append(PenultimateRotamer(
    name              = "t 0 ",
    observations      = 438,
    percent           = 0.210000,
    percent_alpha     = 0.08,
    percent_beta      = 0.44,
    percent_other     = 0.2,
    attri             = "major",

    chi1_mode         = -176,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (10, -10),
    chi1_hwhhfc       = 12,

    chi2_mode         = 1,
    chi2_comm         = 0,
    chi2_range        = (-50, 50),
    chi2_hwhhfa       = (16, -16),
    chi2_hwhhfc       = 30,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ASP"].append(PenultimateRotamer(
    name              = "t70",
    observations      = 118,
    percent           = 0.060000,
    percent_alpha     = 0.11,
    percent_beta      = 0.07,
    percent_other     = 0.04,
    attri             = "major",

    chi1_mode         = -179,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = 12,

    chi2_mode         = 65,
    chi2_comm         = 65,
    chi2_range        = (50, 90),
    chi2_hwhhfa       = (13, -13),
    chi2_hwhhfc       = 18,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ASP"].append(PenultimateRotamer(
    name              = "m -20 ",
    observations      = 1088,
    percent           = 0.510000,
    percent_alpha     = 0.77,
    percent_beta      = 0.38,
    percent_other     = 0.47,
    attri             = "major",

    chi1_mode         = -71,
    chi1_comm         = -70,
    chi1_range        = (-100, -40),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = 10,

    chi2_mode         = -15,
    chi2_comm         = -15,
    chi2_range        = (-90, 20),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = 16,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ASN"].append(PenultimateRotamer(
    name              = "p-10",
    observations      = 103,
    percent           = 0.070000,
    percent_alpha     = 0.0,
    percent_beta      = 0.01,
    percent_other     = 0.1,
    attri             = "major",

    chi1_mode         = 63,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (8, -3),
    chi1_hwhhfc       = 8,

    chi2_mode         = -13,
    chi2_comm         = -10,
    chi2_range        = (-90, 0),
    chi2_hwhhfa       = (10, -4),
    chi2_hwhhfc       = 9,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ASN"].append(PenultimateRotamer(
    name              = "p30",
    observations      = 132,
    percent           = 0.090000,
    percent_alpha     = 0.0,
    percent_beta      = 0.07,
    percent_other     = 0.12,
    attri             = "major",

    chi1_mode         = 64,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (3, -3),
    chi1_hwhhfc       = 6,

    chi2_mode         = 34,
    chi2_comm         = 30,
    chi2_range        = (0, 90),
    chi2_hwhhfa       = (5, -5),
    chi2_hwhhfc       = 7,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ASN"].append(PenultimateRotamer(
    name              = "t-20",
    observations      = 177,
    percent           = 0.120000,
    percent_alpha     = 0.05,
    percent_beta      = 0.21,
    percent_other     = 0.12,
    attri             = "major",

    chi1_mode         = -174,
    chi1_comm         = -174,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (3, -3),
    chi1_hwhhfc       = 5,

    chi2_mode         = -20,
    chi2_comm         = -20,
    chi2_range        = (-120, 0),
    chi2_hwhhfa       = (11, -26),
    chi2_hwhhfc       = 21,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ASN"].append(PenultimateRotamer(
    name              = "t30",
    observations      = 228,
    percent           = 0.150000,
    percent_alpha     = 0.13,
    percent_beta      = 0.18,
    percent_other     = 0.15,
    attri             = "major",

    chi1_mode         = -168,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (11, -11),
    chi1_hwhhfc       = 14,

    chi2_mode         = 31,
    chi2_comm         = 30,
    chi2_range        = (0, 80),
    chi2_hwhhfa       = (23, -17),
    chi2_hwhhfc       = 22,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ASN"].append(PenultimateRotamer(
    name              = "m-20",
    observations      = 580,
    percent           = 0.390000,
    percent_alpha     = 0.65,
    percent_beta      = 0.28,
    percent_other     = 0.33,
    attri             = "major",

    chi1_mode         = -71,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = 10,

    chi2_mode         = -23,
    chi2_comm         = -20,
    chi2_range        = (-60, 10),
    chi2_hwhhfa       = (17, -17),
    chi2_hwhhfc       = 20,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ASN"].append(PenultimateRotamer(
    name              = "m-80",
    observations      = 118,
    percent           = 0.080000,
    percent_alpha     = 0.08,
    percent_beta      = 0.09,
    percent_other     = 0.08,
    attri             = "major",

    chi1_mode         = -71,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = 9,

    chi2_mode         = -76,
    chi2_comm         = -75,
    chi2_range        = (-100, -60),
    chi2_hwhhfa       = (7, -7),
    chi2_hwhhfc       = 9,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ASN"].append(PenultimateRotamer(
    name              = "m120",
    observations      = 58,
    percent           = 0.040000,
    percent_alpha     = 0.03,
    percent_beta      = 0.03,
    percent_other     = 0.04,
    attri             = "major",

    chi1_mode         = -64,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = 9,

    chi2_mode         = 132,
    chi2_comm         = 120,
    chi2_range        = (60, 160),
    chi2_hwhhfa       = (7, -24),
    chi2_hwhhfc       = 18,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ILE"].append(PenultimateRotamer(
    name              = "pp",
    observations      = 10,
    percent           = 0.010000,
    percent_alpha     = 0.0,
    percent_beta      = 0.01,
    percent_other     = 0.0,
    attri             = "minor",

    chi1_mode         = None,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = None,

    chi2_mode         = None,
    chi2_comm         = 100,
    chi2_range        = (70, 130),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ILE"].append(PenultimateRotamer(
    name              = "pt",
    observations      = 216,
    percent           = 0.130000,
    percent_alpha     = 0.04,
    percent_beta      = 0.13,
    percent_other     = 0.22,
    attri             = "major",

    chi1_mode         = 61,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = 10,

    chi2_mode         = 171,
    chi2_comm         = 170,
    chi2_range        = (140, -160),
    chi2_hwhhfa       = (7, -7),
    chi2_hwhhfc       = 10,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ILE"].append(PenultimateRotamer(
    name              = "tp",
    observations      = 36,
    percent           = 0.020000,
    percent_alpha     = 0.02,
    percent_beta      = 0.01,
    percent_other     = 0.04,
    attri             = "major",

    chi1_mode         = -169,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (10, -10),
    chi1_hwhhfc       = 13,

    chi2_mode         = 66,
    chi2_comm         = 66,
    chi2_range        = (35, 95),
    chi2_hwhhfa       = (9, -9),
    chi2_hwhhfc       = 11,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ILE"].append(PenultimateRotamer(
    name              = "tt",
    observations      = 127,
    percent           = 0.080000,
    percent_alpha     = 0.01,
    percent_beta      = 0.08,
    percent_other     = 0.14,
    attri             = "major",

    chi1_mode         = -174,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (13, -8),
    chi1_hwhhfc       = 13,

    chi2_mode         = 167,
    chi2_comm         = 170,
    chi2_range        = (135, -165),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = 11,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ILE"].append(PenultimateRotamer(
    name              = "mp",
    observations      = 19,
    percent           = 0.010000,
    percent_alpha     = 0.0,
    percent_beta      = 0.02,
    percent_other     = 0.01,
    attri             = "minor",

    chi1_mode         = None,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (2, -2),
    chi1_hwhhfc       = None,

    chi2_mode         = None,
    chi2_comm         = 100,
    chi2_range        = (70, 130),
    chi2_hwhhfa       = (2, -2),
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ILE"].append(PenultimateRotamer(
    name              = "mt",
    observations      = 993,
    percent           = 0.600000,
    percent_alpha     = 0.81,
    percent_beta      = 0.58,
    percent_other     = 0.41,
    attri             = "major",

    chi1_mode         = -66,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = 10,

    chi2_mode         = 169,
    chi2_comm         = 170,
    chi2_range        = (140, -160),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = 10,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["ILE"].append(PenultimateRotamer(
    name              = "mm",
    observations      = 242,
    percent           = 0.150000,
    percent_alpha     = 0.1,
    percent_beta      = 0.16,
    percent_other     = 0.17,
    attri             = "major",

    chi1_mode         = -57,
    chi1_comm         = -57,
    chi1_range        = (-85, -25),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = 10,

    chi2_mode         = -59,
    chi2_comm         = -59,
    chi2_range        = (-90, -30),
    chi2_hwhhfa       = (7, -7),
    chi2_hwhhfc       = 10,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["LEU"].append(PenultimateRotamer(
    name              = "pp",
    observations      = 21,
    percent           = 0.010000,
    percent_alpha     = 0.0,
    percent_beta      = 0.02,
    percent_other     = 0.01,
    attri             = "minor",

    chi1_mode         = None,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = None,
    chi1_hwhhfc       = None,

    chi2_mode         = None,
    chi2_comm         = 80,
    chi2_range        = (50, 110),
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["LEU"].append(PenultimateRotamer(
    name              = "tp",
    observations      = 750,
    percent           = 0.290000,
    percent_alpha     = 0.3,
    percent_beta      = 0.36,
    percent_other     = 0.23,
    attri             = "major",

    chi1_mode         = 177,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = 10,

    chi2_mode         = 63,
    chi2_comm         = 65,
    chi2_range        = (35, 95),
    chi2_hwhhfa       = (7, -7),
    chi2_hwhhfc       = 10,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["LEU"].append(PenultimateRotamer(
    name              = "tt",
    observations      = 49,
    percent           = 0.020000,
    percent_alpha     = 0.01,
    percent_beta      = 0.03,
    percent_other     = 0.01,
    attri             = "major",

    chi1_mode         = -172,
    chi1_comm         = -172,
    chi1_range        = (160, -140),
    chi1_hwhhfa       = (6, -6),
    chi1_hwhhfc       = 9,

    chi2_mode         = 147,
    chi2_comm         = 147,
    chi2_range        = (120, 180),
    chi2_hwhhfa       = (6, -6),
    chi2_hwhhfc       = 9,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["LEU"].append(PenultimateRotamer(
    name              = "mp",
    observations      = 63,
    percent           = 0.020000,
    percent_alpha     = 0.01,
    percent_beta      = 0.05,
    percent_other     = 0.02,
    attri             = "major",

    chi1_mode         = -85,
    chi1_comm         = -85,
    chi1_range        = (-115, -55),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = 11,

    chi2_mode         = 66,
    chi2_comm         = 65,
    chi2_range        = (45, 105),
    chi2_hwhhfa       = (8, -12),
    chi2_hwhhfc       = 14,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["LEU"].append(PenultimateRotamer(
    name              = "mt",
    observations      = 1548,
    percent           = 0.590000,
    percent_alpha     = 0.62,
    percent_beta      = 0.46,
    percent_other     = 0.66,
    attri             = "major",

    chi1_mode         = -65,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = 11,

    chi2_mode         = 174,
    chi2_comm         = 174,
    chi2_range        = (150, -150),
    chi2_hwhhfa       = (9, -9),
    chi2_hwhhfc       = 11,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["HIS"].append(PenultimateRotamer(
    name              = "p-80",
    observations      = 51,
    percent           = 0.090000,
    percent_alpha     = 0.0,
    percent_beta      = 0.06,
    percent_other     = 0.13,
    attri             = "major",

    chi1_mode         = 60,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = 10,

    chi2_mode         = -75,
    chi2_comm         = -75,
    chi2_range        = (-120, -50),
    chi2_hwhhfa       = (9, -9),
    chi2_hwhhfc       = 12,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["HIS"].append(PenultimateRotamer(
    name              = "p80",
    observations      = 26,
    percent           = 0.040000,
    percent_alpha     = 0.0,
    percent_beta      = 0.04,
    percent_other     = 0.06,
    attri             = "major",

    chi1_mode         = 61,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (13, -13),
    chi1_hwhhfc       = 13,

    chi2_mode         = 78,
    chi2_comm         = 78,
    chi2_range        = (50, 120),
    chi2_hwhhfa       = (7, -7),
    chi2_hwhhfc       = 10,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["HIS"].append(PenultimateRotamer(
    name              = "t-160",
    observations      = 31,
    percent           = 0.050000,
    percent_alpha     = 0.05,
    percent_beta      = 0.14,
    percent_other     = 0.01,
    attri             = "major",

    chi1_mode         = -178,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (10, -10),
    chi1_hwhhfc       = 12,

    chi2_mode         = -163,
    chi2_comm         = -163,
    chi2_range        = (150, -120),
    chi2_hwhhfa       = (27, -9),
    chi2_hwhhfc       = 20,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["HIS"].append(PenultimateRotamer(
    name              = "t-80",
    observations      = 64,
    percent           = 0.110000,
    percent_alpha     = 0.17,
    percent_beta      = 0.09,
    percent_other     = 0.09,
    attri             = "major",

    chi1_mode         = -173,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = 10,

    chi2_mode         = -81,
    chi2_comm         = -81,
    chi2_range        = (-120, -50),
    chi2_hwhhfa       = (12, -28),
    chi2_hwhhfc       = 22,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["HIS"].append(PenultimateRotamer(
    name              = "t60",
    observations      = 94,
    percent           = 0.160000,
    percent_alpha     = 0.24,
    percent_beta      = 0.17,
    percent_other     = 0.12,
    attri             = "major",

    chi1_mode         = -178,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (10, -10),
    chi1_hwhhfc       = 13,

    chi2_mode         = 62,
    chi2_comm         = 62,
    chi2_range        = (50, 120),
    chi2_hwhhfa       = (25, -8),
    chi2_hwhhfc       = 19,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["HIS"].append(PenultimateRotamer(
    name              = "m-70",
    observations      = 174,
    percent           = 0.290000,
    percent_alpha     = 0.26,
    percent_beta      = 0.3,
    percent_other     = 0.3,
    attri             = "major",

    chi1_mode         = -60,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = 11,

    chi2_mode         = -69,
    chi2_comm         = -69,
    chi2_range        = (-120, -30),
    chi2_hwhhfa       = (20, -20),
    chi2_hwhhfc       = 23,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["HIS"].append(PenultimateRotamer(
    name              = "m170",
    observations      = 44,
    percent           = 0.070000,
    percent_alpha     = 0.09,
    percent_beta      = 0.03,
    percent_other     = 0.09,
    attri             = "major",

    chi1_mode         = -63,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = 10,

    chi2_mode         = 165,
    chi2_comm         = 165,
    chi2_range        = (120, -160),
    chi2_hwhhfa       = (8, -20),
    chi2_hwhhfc       = 16,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["HIS"].append(PenultimateRotamer(
    name              = "m80",
    observations      = 78,
    percent           = 0.130000,
    percent_alpha     = 0.14,
    percent_beta      = 0.1,
    percent_other     = 0.14,
    attri             = "major",

    chi1_mode         = -66,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = 11,

    chi2_mode         = 83,
    chi2_comm         = 83,
    chi2_range        = (50, 120),
    chi2_hwhhfa       = (21, -11),
    chi2_hwhhfc       = 18,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["TRP"].append(PenultimateRotamer(
    name              = "p-90",
    observations      = 67,
    percent           = 0.110000,
    percent_alpha     = 0.02,
    percent_beta      = 0.13,
    percent_other     = 0.14,
    attri             = "major",

    chi1_mode         = 58,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (12, -7),
    chi1_hwhhfc       = 12,

    chi2_mode         = -87,
    chi2_comm         = -90,
    chi2_range        = (-130, -60),
    chi2_hwhhfa       = (7, -7),
    chi2_hwhhfc       = 10,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["TRP"].append(PenultimateRotamer(
    name              = "p90",
    observations      = 34,
    percent           = 0.060000,
    percent_alpha     = 0.01,
    percent_beta      = 0.09,
    percent_other     = 0.06,
    attri             = "major",

    chi1_mode         = 60,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (13, -7),
    chi1_hwhhfc       = 12,

    chi2_mode         = 92,
    chi2_comm         = 90,
    chi2_range        = (60, 130),
    chi2_hwhhfa       = (6, -6),
    chi2_hwhhfc       = 8,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["TRP"].append(PenultimateRotamer(
    name              = "t-105",
    observations      = 100,
    percent           = 0.160000,
    percent_alpha     = 0.27,
    percent_beta      = 0.1,
    percent_other     = 0.14,
    attri             = "major",

    chi1_mode         = 178,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (16, -12),
    chi1_hwhhfc       = 16,

    chi2_mode         = -105,
    chi2_comm         = -105,
    chi2_range        = (-130, -60),
    chi2_hwhhfa       = (12, -12),
    chi2_hwhhfc       = 14,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["TRP"].append(PenultimateRotamer(
    name              = "t90",
    observations      = 109,
    percent           = 0.180000,
    percent_alpha     = 0.28,
    percent_beta      = 0.14,
    percent_other     = 0.15,
    attri             = "major",

    chi1_mode         = -178,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = 10,

    chi2_mode         = 88,
    chi2_comm         = 90,
    chi2_range        = (0, 100),
    chi2_hwhhfa       = (9, -9),
    chi2_hwhhfc       = 11,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["TRP"].append(PenultimateRotamer(
    name              = "m-90",
    observations      = 31,
    percent           = 0.050000,
    percent_alpha     = 0.0,
    percent_beta      = 0.07,
    percent_other     = 0.07,
    attri             = "major",

    chi1_mode         = -70,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (6, -6),
    chi1_hwhhfc       = 9,

    chi2_mode         = -87,
    chi2_comm         = -90,
    chi2_range        = (-130, -60),
    chi2_hwhhfa       = (10, -10),
    chi2_hwhhfc       = 12,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["TRP"].append(PenultimateRotamer(
    name              = "m0",
    observations      = 48,
    percent           = 0.080000,
    percent_alpha     = 0.15,
    percent_beta      = 0.02,
    percent_other     = 0.08,
    attri             = "major",

    chi1_mode         = -66,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = 9,

    chi2_mode         = -4,
    chi2_comm         = -4,
    chi2_range        = (-40, 20),
    chi2_hwhhfa       = (9, -26),
    chi2_hwhhfc       = 20,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["TRP"].append(PenultimateRotamer(
    name              = "m95",
    observations      = 195,
    percent           = 0.320000,
    percent_alpha     = 0.22,
    percent_beta      = 0.43,
    percent_other     = 0.29,
    attri             = "major",

    chi1_mode         = -69,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = 11,

    chi2_mode         = 95,
    chi2_comm         = 95,
    chi2_range        = (60, 130),
    chi2_hwhhfa       = (16, -16),
    chi2_hwhhfc       = 19,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["TYR"].append(PenultimateRotamer(
    name              = "p90",
    observations      = 182,
    percent           = 0.130000,
    percent_alpha     = 0.01,
    percent_beta      = 0.21,
    percent_other     = 0.12,
    attri             = "major",

    chi1_mode         = 63,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (10, -10),
    chi1_hwhhfc       = 13,

    chi2_mode         = 89,
    chi2_comm         = 90,
    chi2_range        = (60, 90),
    chi2_hwhhfa       = (10, -10),
    chi2_hwhhfc       = 13,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["TYR"].append(PenultimateRotamer(
    name              = "t80",
    observations      = 486,
    percent           = 0.340000,
    percent_alpha     = 0.55,
    percent_beta      = 0.25,
    percent_other     = 0.3,
    attri             = "major",

    chi1_mode         = 176,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = 11,

    chi2_mode         = 77,
    chi2_comm         = 80,
    chi2_range        = (20, 90),
    chi2_hwhhfa       = (12, -12),
    chi2_hwhhfc       = 14,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["TYR"].append(PenultimateRotamer(
    name              = "m-85",
    observations      = 618,
    percent           = 0.430000,
    percent_alpha     = 0.26,
    percent_beta      = 0.5,
    percent_other     = 0.45,
    attri             = "major",

    chi1_mode         = -65,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = 11,

    chi2_mode         = -87,
    chi2_comm         = -85,
    chi2_range        = (50, 90),
    chi2_hwhhfa       = (16, -16),
    chi2_hwhhfc       = 21,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["TYR"].append(PenultimateRotamer(
    name              = "m -30 ",
    observations      = 124,
    percent           = 0.090000,
    percent_alpha     = 0.15,
    percent_beta      = 0.04,
    percent_other     = 0.09,
    attri             = "major",

    chi1_mode         = -64,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = 11,

    chi2_mode         = -42,
    chi2_comm         = -30,
    chi2_range        = (-50, 0),
    chi2_hwhhfa       = (31, -6),
    chi2_hwhhfc       = 18,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["PHE"].append(PenultimateRotamer(
    name              = "p90",
    observations      = 202,
    percent           = 0.130000,
    percent_alpha     = 0.01,
    percent_beta      = 0.24,
    percent_other     = 0.11,
    attri             = "major",

    chi1_mode         = 59,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = 11,

    chi2_mode         = 88,
    chi2_comm         = 90,
    chi2_range        = (60, 90),
    chi2_hwhhfa       = (9, -9),
    chi2_hwhhfc       = 11,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["PHE"].append(PenultimateRotamer(
    name              = "t80",
    observations      = 522,
    percent           = 0.330000,
    percent_alpha     = 0.57,
    percent_beta      = 0.18,
    percent_other     = 0.29,
    attri             = "major",

    chi1_mode         = 177,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (11, -11),
    chi1_hwhhfc       = 13,

    chi2_mode         = 80,
    chi2_comm         = 80,
    chi2_range        = (20, 90),
    chi2_hwhhfa       = (11, -18),
    chi2_hwhhfc       = 17,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["PHE"].append(PenultimateRotamer(
    name              = "m-85",
    observations      = 697,
    percent           = 0.440000,
    percent_alpha     = 0.29,
    percent_beta      = 0.51,
    percent_other     = 0.47,
    attri             = "major",

    chi1_mode         = -64,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = 12,

    chi2_mode         = -83,
    chi2_comm         = -85,
    chi2_range        = (50, 90),
    chi2_hwhhfa       = (14, -14),
    chi2_hwhhfc       = 17,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["PHE"].append(PenultimateRotamer(
    name              = "m -30 ",
    observations      = 149,
    percent           = 0.090000,
    percent_alpha     = 0.12,
    percent_beta      = 0.05,
    percent_other     = 0.11,
    attri             = "major",

    chi1_mode         = -64,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = 9,

    chi2_mode         = -19,
    chi2_comm         = -30,
    chi2_range        = (-50, 0),
    chi2_hwhhfa       = (14, -21),
    chi2_hwhhfc       = 20,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["PRO"].append(PenultimateRotamer(
    name              = "Cg_endo",
    observations      = 379,
    percent           = 0.440000,
    percent_alpha     = 0.23,
    percent_beta      = 0.54,
    percent_other     = 0.43,
    attri             = "major",

    chi1_mode         = 30,
    chi1_comm         = 30,
    chi1_range        = (15, 60),
    chi1_hwhhfa       = (5, -5),
    chi1_hwhhfc       = 7,

    chi2_mode         = None,
    chi2_comm         = None,
    chi2_range        = None,
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["PRO"].append(PenultimateRotamer(
    name              = "Cg_exo",
    observations      = 372,
    percent           = 0.430000,
    percent_alpha     = 0.68,
    percent_beta      = 0.28,
    percent_other     = 0.44,
    attri             = "major",

    chi1_mode         = -29,
    chi1_comm         = -30,
    chi1_range        = (-60, -15),
    chi1_hwhhfa       = (5, -5),
    chi1_hwhhfc       = 6,

    chi2_mode         = None,
    chi2_comm         = None,
    chi2_range        = None,
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["PRO"].append(PenultimateRotamer(
    name              = "cis_Cg_endo",
    observations      = 56,
    percent           = 0.060000,
    percent_alpha     = 0.0,
    percent_beta      = 0.01,
    percent_other     = 0.07,
    attri             = "major",

    chi1_mode         = 31,
    chi1_comm         = 30,
    chi1_range        = (15, 60),
    chi1_hwhhfa       = (4, -4),
    chi1_hwhhfc       = 5,

    chi2_mode         = None,
    chi2_comm         = None,
    chi2_range        = None,
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["THR"].append(PenultimateRotamer(
    name              = "p",
    observations      = 1200,
    percent           = 0.490000,
    percent_alpha     = 0.25,
    percent_beta      = 0.31,
    percent_other     = 0.65,
    attri             = "major",

    chi1_mode         = 59,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = 10,

    chi2_mode         = None,
    chi2_comm         = None,
    chi2_range        = None,
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["THR"].append(PenultimateRotamer(
    name              = "t",
    observations      = 169,
    percent           = 0.070000,
    percent_alpha     = 0.0,
    percent_beta      = 0.13,
    percent_other     = 0.06,
    attri             = "major",

    chi1_mode         = -171,
    chi1_comm         = -175,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (10, -10),
    chi1_hwhhfc       = 6,

    chi2_mode         = None,
    chi2_comm         = None,
    chi2_range        = None,
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["THR"].append(PenultimateRotamer(
    name              = "m",
    observations      = 1062,
    percent           = 0.430000,
    percent_alpha     = 0.74,
    percent_beta      = 0.55,
    percent_other     = 0.29,
    attri             = "major",

    chi1_mode         = -61,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = 7,

    chi2_mode         = None,
    chi2_comm         = None,
    chi2_range        = None,
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["VAL"].append(PenultimateRotamer(
    name              = "p",
    observations      = 169,
    percent           = 0.060000,
    percent_alpha     = 0.02,
    percent_beta      = 0.08,
    percent_other     = 0.08,
    attri             = "major",

    chi1_mode         = 63,
    chi1_comm         = 63,
    chi1_range        = (35, 95),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = 8,

    chi2_mode         = None,
    chi2_comm         = -177,
    chi2_range        = None,
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["VAL"].append(PenultimateRotamer(
    name              = "t",
    observations      = 1931,
    percent           = 0.730000,
    percent_alpha     = 0.9,
    percent_beta      = 0.72,
    percent_other     = 0.63,
    attri             = "major",

    chi1_mode         = 175,
    chi1_comm         = 175,
    chi1_range        = (145, -155),
    chi1_hwhhfa       = (4, -4),
    chi1_hwhhfc       = 8,

    chi2_mode         = None,
    chi2_comm         = -65,
    chi2_range        = None,
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["VAL"].append(PenultimateRotamer(
    name              = "m",
    observations      = 526,
    percent           = 0.200000,
    percent_alpha     = 0.07,
    percent_beta      = 0.2,
    percent_other     = 0.28,
    attri             = "major",

    chi1_mode         = -64,
    chi1_comm         = -60,
    chi1_range        = (-90, -30),
    chi1_hwhhfa       = (5, -5),
    chi1_hwhhfc       = 7,

    chi2_mode         = None,
    chi2_comm         = 60,
    chi2_range        = None,
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["SER"].append(PenultimateRotamer(
    name              = "p",
    observations      = 1201,
    percent           = 0.480000,
    percent_alpha     = 0.33,
    percent_beta      = 0.36,
    percent_other     = 0.55,
    attri             = "major",

    chi1_mode         = 64,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = 10,

    chi2_mode         = None,
    chi2_comm         = None,
    chi2_range        = None,
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["SER"].append(PenultimateRotamer(
    name              = "t",
    observations      = 541,
    percent           = 0.220000,
    percent_alpha     = 0.22,
    percent_beta      = 0.34,
    percent_other     = 0.18,
    attri             = "major",

    chi1_mode         = 178,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (7, -7),
    chi1_hwhhfc       = 11,

    chi2_mode         = None,
    chi2_comm         = None,
    chi2_range        = None,
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["SER"].append(PenultimateRotamer(
    name              = "m",
    observations      = 714,
    percent           = 0.290000,
    percent_alpha     = 0.44,
    percent_beta      = 0.29,
    percent_other     = 0.25,
    attri             = "major",

    chi1_mode         = -65,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (6, -6),
    chi1_hwhhfc       = 9,

    chi2_mode         = None,
    chi2_comm         = None,
    chi2_range        = None,
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["CYS"].append(PenultimateRotamer(
    name              = "p",
    observations      = 64,
    percent           = 0.230000,
    percent_alpha     = 0.05,
    percent_beta      = 0.23,
    percent_other     = 0.34,
    attri             = "major",

    chi1_mode         = 55,
    chi1_comm         = 62,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (22, -3),
    chi1_hwhhfc       = 14,

    chi2_mode         = None,
    chi2_comm         = None,
    chi2_range        = None,
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["CYS"].append(PenultimateRotamer(
    name              = "t",
    observations      = 74,
    percent           = 0.260000,
    percent_alpha     = 0.2,
    percent_beta      = 0.45,
    percent_other     = 0.21,
    attri             = "major",

    chi1_mode         = -177,
    chi1_comm         = -177,
    chi1_range        = (155, -145),
    chi1_hwhhfa       = (12, -5),
    chi1_hwhhfc       = 10,

    chi2_mode         = None,
    chi2_comm         = None,
    chi2_range        = None,
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["CYS"].append(PenultimateRotamer(
    name              = "m",
    observations      = 142,
    percent           = 0.500000,
    percent_alpha     = 0.75,
    percent_beta      = 0.32,
    percent_other     = 0.43,
    attri             = "major",

    chi1_mode         = -65,
    chi1_comm         = -65,
    chi1_range        = (-95, -35),
    chi1_hwhhfa       = (15, -5),
    chi1_hwhhfc       = 11,

    chi2_mode         = None,
    chi2_comm         = None,
    chi2_range        = None,
    chi2_hwhhfa       = None,
    chi2_hwhhfc       = None,

    chi3_mode         = None,
    chi3_comm         = None,
    chi3_range        = None,
    chi3_hwhhfa       = None,
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["DIS"].append(PenultimateRotamer(
    name              = "mmm",
    observations      = 70,
    percent           = 0.360000,
    percent_alpha     = None,
    percent_beta      = None,
    percent_other     = None,
    attri             = "major",

    chi1_mode         = -61,
    chi1_comm         = None,
    chi1_range        = (-95, -30),
    chi1_hwhhfa       = (11, -11),
    chi1_hwhhfc       = None,

    chi2_mode         = -81,
    chi2_comm         = None,
    chi2_range        = (-110, -50),
    chi2_hwhhfa       = (9, -9),
    chi2_hwhhfc       = None,

    chi3_mode         = -75,
    chi3_comm         = None,
    chi3_range        = (-120, -40),
    chi3_hwhhfa       = (16, -16),
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["DIS"].append(PenultimateRotamer(
    name              = "ppp",
    observations      = 15,
    percent           = 0.080000,
    percent_alpha     = None,
    percent_beta      = None,
    percent_other     = None,
    attri             = "major",

    chi1_mode         = 63,
    chi1_comm         = None,
    chi1_range        = (30, 90),
    chi1_hwhhfa       = (10, -10),
    chi1_hwhhfc       = None,

    chi2_mode         = 85,
    chi2_comm         = None,
    chi2_range        = (55, 115),
    chi2_hwhhfa       = (14, -9),
    chi2_hwhhfc       = None,

    chi3_mode         = 85,
    chi3_comm         = None,
    chi3_range        = (60, 115),
    chi3_hwhhfa       = (11, -11),
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["DIS"].append(PenultimateRotamer(
    name              = "mpp",
    observations      = 33,
    percent           = 0.170000,
    percent_alpha     = None,
    percent_beta      = None,
    percent_other     = None,
    attri             = "major",

    chi1_mode         = -65,
    chi1_comm         = None,
    chi1_range        = (-100, -35),
    chi1_hwhhfa       = (12, -12),
    chi1_hwhhfc       = None,

    chi2_mode         = 100,
    chi2_comm         = None,
    chi2_range        = (70, 130),
    chi2_hwhhfa       = (9, -9),
    chi2_hwhhfc       = None,

    chi3_mode         = 85,
    chi3_comm         = None,
    chi3_range        = (50, 115),
    chi3_hwhhfa       = (14, -14),
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["DIS"].append(PenultimateRotamer(
    name              = "pmm",
    observations      = 6,
    percent           = 0.030000,
    percent_alpha     = None,
    percent_beta      = None,
    percent_other     = None,
    attri             = "minor",

    chi1_mode         = 90,
    chi1_comm         = None,
    chi1_range        = (60, 120),
    chi1_hwhhfa       = (10, -18),
    chi1_hwhhfc       = None,

    chi2_mode         = -91,
    chi2_comm         = None,
    chi2_range        = (-120, -60),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = None,

    chi3_mode         = -64,
    chi3_comm         = None,
    chi3_range        = (-95, -35),
    chi3_hwhhfa       = (15, -15),
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["DIS"].append(PenultimateRotamer(
    name              = "mpm",
    observations      = 11,
    percent           = 0.060000,
    percent_alpha     = None,
    percent_beta      = None,
    percent_other     = None,
    attri             = "major",

    chi1_mode         = -86,
    chi1_comm         = None,
    chi1_range        = (-120, -40),
    chi1_hwhhfa       = (11, -11),
    chi1_hwhhfc       = None,

    chi2_mode         = 102,
    chi2_comm         = None,
    chi2_range        = (70, 130),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = None,

    chi3_mode         = -102,
    chi3_comm         = None,
    chi3_range        = (-160, -90),
    chi3_hwhhfa       = (9, -9),
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["DIS"].append(PenultimateRotamer(
    name              = "mmt",
    observations      = 19,
    percent           = 0.100000,
    percent_alpha     = None,
    percent_beta      = None,
    percent_other     = None,
    attri             = "major",

    chi1_mode         = -92,
    chi1_comm         = None,
    chi1_range        = (-120, -30),
    chi1_hwhhfa       = (15, -10),
    chi1_hwhhfc       = None,

    chi2_mode         = -90,
    chi2_comm         = None,
    chi2_range        = (-120, -60),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = None,

    chi3_mode         = -149,
    chi3_comm         = None,
    chi3_range        = (150, -120),
    chi3_hwhhfa       = (12, -12),
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["DIS"].append(PenultimateRotamer(
    name              = "ppt",
    observations      = 16,
    percent           = 0.080000,
    percent_alpha     = None,
    percent_beta      = None,
    percent_other     = None,
    attri             = "major",

    chi1_mode         = 52,
    chi1_comm         = None,
    chi1_range        = (30, 95),
    chi1_hwhhfa       = (10, -10),
    chi1_hwhhfc       = None,

    chi2_mode         = 82,
    chi2_comm         = None,
    chi2_range        = (50, 110),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = None,

    chi3_mode         = 180,
    chi3_comm         = None,
    chi3_range        = (150, -120),
    chi3_hwhhfa       = (9, -9),
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["DIS"].append(PenultimateRotamer(
    name              = "mpt",
    observations      = 5,
    percent           = 0.030000,
    percent_alpha     = None,
    percent_beta      = None,
    percent_other     = None,
    attri             = "minor",

    chi1_mode         = -68,
    chi1_comm         = None,
    chi1_range        = (-100, -40),
    chi1_hwhhfa       = (13, -8),
    chi1_hwhhfc       = None,

    chi2_mode         = 96,
    chi2_comm         = None,
    chi2_range        = (65, 125),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = None,

    chi3_mode         = 147,
    chi3_comm         = None,
    chi3_range        = (115, 175),
    chi3_hwhhfa       = (14, -14),
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["DIS"].append(PenultimateRotamer(
    name              = "tmt",
    observations      = 6,
    percent           = 0.030000,
    percent_alpha     = None,
    percent_beta      = None,
    percent_other     = None,
    attri             = "minor",

    chi1_mode         = 172,
    chi1_comm         = None,
    chi1_range        = (140, -160),
    chi1_hwhhfa       = (9, -9),
    chi1_hwhhfc       = None,

    chi2_mode         = -83,
    chi2_comm         = None,
    chi2_range        = (-115, -55),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = None,

    chi3_mode         = -168,
    chi3_comm         = None,
    chi3_range        = (160, -140),
    chi3_hwhhfa       = (11, -11),
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

PenultimateRotamerMap["DIS"].append(PenultimateRotamer(
    name              = "tpt",
    observations      = 1,
    percent           = 0.010000,
    percent_alpha     = None,
    percent_beta      = None,
    percent_other     = None,
    attri             = "minor",

    chi1_mode         = 122,
    chi1_comm         = None,
    chi1_range        = (95, 150),
    chi1_hwhhfa       = (8, -8),
    chi1_hwhhfc       = None,

    chi2_mode         = 87,
    chi2_comm         = None,
    chi2_range        = (55, 155),
    chi2_hwhhfa       = (8, -8),
    chi2_hwhhfc       = None,

    chi3_mode         = 163,
    chi3_comm         = None,
    chi3_range        = (130, -170),
    chi3_hwhhfa       = (8, -8),
    chi3_hwhhfc       = None,

    chi4_mode         = None,
    chi4_comm         = None,
    chi4_range        = None,
    chi4_hwhhfa       = None,
    chi4_hwhhfc       = None))

### <TESTING>
if __name__ == "__main__":
    for (res_name, rot_list) in PenultimateRotamerMap.items():
        print res_name
        for rot in rot_list:
            print "    ",rot
### </TESTING>
