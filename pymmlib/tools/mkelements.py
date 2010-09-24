#!/usr/bin/env python
## Copyright 2002-2010 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
"""This program creates "elements.cif" using the data in this file, and the
CCP4 "atomsf.lib" library file.
"""

## Python
import os
import sys

## pymmlib
from mmLib.mmCIF import *

class Element(object):
    """Class for holding the properties of a atomic element.
    """
    def __init__(self,
                 name                 = "",
                 symbol               = "",
                 group                = "",
                 period               = "",
                 atomic_number        = 0,
                 atomic_weight        = 0.0,
                 atomic_radius        = 0.0,
                 covalent_radius      = 0.0,
                 van_der_waals_radius = 2.0,
                 electronegativity    = 0.0,
                 color                = (1.0, 1.0, 1.0)):

        self.name                 = name
        self.symbol               = symbol
        self.group                = group
        self.period               = period
        self.atomic_number        = atomic_number
        self.atomic_weight        = atomic_weight
        self.atomic_radius        = atomic_radius
        self.covalent_radius      = covalent_radius
        self.van_der_waals_radius = van_der_waals_radius
        self.electronegativity    = electronegativity
        self.color                = color

    def __str__(self):
        return "Element=%s" % (self.name)

H  = Element(
    name                 = "Hydrogen",
    symbol               = "H",
    atomic_number        = 1,
    atomic_weight        = 1.007940,
    van_der_waals_radius = 1.20)

He = Element(
    name                 = "Helium",
    symbol               = "He",
    atomic_number        = 2,
    atomic_weight        = 4.002602,
    van_der_waals_radius = 1.40)

Li = Element(
    name                 = "Lithium",
    symbol               = "Li",
    atomic_number        = 3,
    atomic_weight        = 6.941000,
    van_der_waals_radius = 1.82)

Be = Element(
    name                 = "Beryllium",
    symbol               = "Be",
    atomic_number        = 4,
    atomic_weight        = 9.012182)

B  = Element(
    name                 = "Boron",
    symbol               = "B",
    atomic_number        = 5,
    atomic_weight        = 10.811000)

C  = Element(
    name                 = "Carbon",
    symbol               = "C",
    atomic_number        = 6,
    atomic_weight        = 12.010700,
    van_der_waals_radius = 1.70,
    color                = (1.0, 1.0, 0.0))

N  = Element(
    name                 = "Nitrogen",
    symbol               = "N",
    atomic_number        = 7,
    atomic_weight        = 14.006700,
    van_der_waals_radius = 1.55,
    color                = (0.0, 0.0, 1.0))

O  = Element(
    name                 = "Oxygen",
    symbol               = "O",
    atomic_number        = 8,
    atomic_weight        = 15.999400,
    van_der_waals_radius = 1.52,
    color                = (1.0, 0.0, 0.0))

F  = Element(
    name                 = "Fluorine",
    symbol               = "F",
    atomic_number        = 9,
    atomic_weight        = 18.998403,
    van_der_waals_radius = 1.47)

Ne = Element(
    name                 = "Neon",
    symbol               = "Ne",
    atomic_number        = 10,
    atomic_weight        = 20.179700,
    van_der_waals_radius = 1.54)

Na = Element(
    name                 = "Sodium",
    symbol               = "Na",
    atomic_number        = 11,
    atomic_weight        = 22.989770,
    van_der_waals_radius = 2.27)

Mg = Element(
    name                 = "Magnesium",
    symbol               = "Mg",
    atomic_number        = 12,
    atomic_weight        = 24.305000,
    van_der_waals_radius = 1.73)

Al = Element(
    name                 = "Aluminium",
    symbol               = "Al",
    atomic_number        = 13,
    atomic_weight        = 26.981538)

Si = Element(
    name                 = "Silicon",
    symbol               = "Si",
    atomic_number        = 14,
    atomic_weight        = 28.085500,
    van_der_waals_radius = 2.10)

P  = Element(
    name                 = "Phosphorus",
    symbol               = "P",
    atomic_number        = 15,
    atomic_weight        = 30.973761,
    van_der_waals_radius = 1.80)

S  = Element(
    name                 = "Sulfur",
    symbol               = "S",
    atomic_number        = 16,
    atomic_weight        = 32.065000,
    van_der_waals_radius = 1.80,
    color                = (0.0, 1.0, 0.0))

Cl = Element(
    name                 = "Chlorine",
    symbol               = "Cl",
    atomic_number        = 17,
    atomic_weight        = 35.453000,
    van_der_waals_radius = 1.75)

Ar = Element(
    name                 = "Argon",
    symbol               = "Ar",
    atomic_number        = 18,
    atomic_weight        = 39.948000,
    van_der_waals_radius = 1.88)

K  = Element(
    name                 = "Potassium",
    symbol               = "K",
    atomic_number        = 19,
    atomic_weight        = 39.098300,
    van_der_waals_radius = 2.75)

Ca = Element(
    name                 = "Calcium",
    symbol               = "Ca",
    atomic_number        = 20,
    atomic_weight        = 40.078000)

Sc = Element(
    name                 = "Scandium",
    symbol               = "Sc",
    atomic_number        = 21,
    atomic_weight        = 44.955910)

Ti = Element(
    name                 = "Titanium",
    symbol               = "Ti",
    atomic_number        = 22,
    atomic_weight        = 47.867000)

V  = Element(
    name                 = "Vanadium",
    symbol               = "V",
    atomic_number        = 23,
    atomic_weight        = 50.941500)

Cr = Element(
    name                 = "Chromium",
    symbol               = "Cr",
    atomic_number        = 24,
    atomic_weight        = 51.996100)

Mn = Element(
    name                 = "Manganese",
    symbol               = "Mn",
    atomic_number        = 25,
    atomic_weight        = 54.938049)

Fe = Element(
    name                 = "Iron",
    symbol               = "Fe",
    atomic_number        = 26,
    atomic_weight        = 55.845000)

Co = Element(
    name                 = "Cobalt",
    symbol               = "Co",
    atomic_number        = 27,
    atomic_weight        = 58.933200)

Ni = Element(
    name                 = "Nickel",
    symbol               = "Ni",
    atomic_number        = 28,
    atomic_weight        = 58.693400,
    van_der_waals_radius = 1.63)

Cu = Element(
    name                 = "Copper",
    symbol               = "Cu",
    atomic_number        = 29,
    atomic_weight        = 63.546000,
    van_der_waals_radius = 1.40)

Zn = Element(
    name                 = "Zinc",
    symbol               = "Zn",
    atomic_number        = 30,
    atomic_weight        = 65.409000,
    van_der_waals_radius = 1.39)

Ga = Element(
    name                 = "Gallium",
    symbol               = "Ga",
    atomic_number        = 31,
    atomic_weight        = 69.723000,
    van_der_waals_radius = 1.87)

Ge = Element(
    name                 = "Germanium",
    symbol               = "Ge",
    atomic_number        = 32,
    atomic_weight        = 72.640000)

Ge = Element(
    name                 = "Germanium",
    symbol               = "Ge",
    atomic_number        = 32,
    atomic_weight        = 72.640000)

As = Element(
    name                 = "Arsenic",
    symbol               = "As",
    atomic_number        = 33,
    atomic_weight        = 74.921600,
    van_der_waals_radius = 1.85)

Se = Element(
    name                 = "Selenium",
    symbol               = "Se",
    atomic_number        = 34,
    atomic_weight        = 78.960000,
    van_der_waals_radius = 1.90)

Br = Element(
    name                 = "Bromine",
    symbol               = "Br",
    atomic_number        = 35,
    atomic_weight        = 79.904000,
    van_der_waals_radius = 1.85)

Kr = Element(
    name                 = "Krypton",
    symbol               = "Kr",
    atomic_number        = 36,
    atomic_weight        = 83.798000,
    van_der_waals_radius = 2.02)

Rb = Element(
    name                 = "Rubidium",
    symbol               = "Rb",
    atomic_number        = 37,
    atomic_weight        = 85.467800)

Rb = Element(
    name                 = "Rubidium",
    symbol               = "Rb",
    atomic_number        = 37,
    atomic_weight        = 85.467800)

Sr = Element(
    name                 = "Strontium",
    symbol               = "Sr",
    atomic_number        = 38,
    atomic_weight        = 87.620000)

Sr = Element(
    name                 = "Strontium",
    symbol               = "Sr",
    atomic_number        = 38,
    atomic_weight        = 87.620000)

Y  = Element(
    name                 = "Yttrium",
    symbol               = "Y",
    atomic_number        = 39,
    atomic_weight        = 88.905850)

Zr = Element(
    name                 = "Zirconium",
    symbol               = "Zr",
    atomic_number        = 40,
    atomic_weight        = 91.224000)

Zr = Element(
    name                 = "Zirconium",
    symbol               = "Zr",
    atomic_number        = 40,
    atomic_weight        = 91.224000)

Nb = Element(
    name                 = "Niobium",
    symbol               = "Nb",
    atomic_number        = 41,
    atomic_weight        = 92.906380)

Nb = Element(
    name                 = "Niobium",
    symbol               = "Nb",
    atomic_number        = 41,
    atomic_weight        = 92.906380)

Mo = Element(
    name                 = "Molybdenum",
    symbol               = "Mo",
    atomic_number        = 42,
    atomic_weight        = 95.940000)

Mo = Element(
    name                 = "Molybdenum",
    symbol               = "Mo",
    atomic_number        = 42,
    atomic_weight        = 95.940000)

Tc = Element(
    name                 = "Technetium",
    symbol               = "Tc",
    atomic_number        = 43,
    atomic_weight        = 98.000000)

Tc = Element(
    name                 = "Technetium",
    symbol               = "Tc",
    atomic_number        = 43,
    atomic_weight        = 98.000000)

Ru = Element(
    name                 = "Ruthenium",
    symbol               = "Ru",
    atomic_number        = 44,
    atomic_weight        = 101.070000)

Ru = Element(
    name                 = "Ruthenium",
    symbol               = "Ru",
    atomic_number        = 44,
    atomic_weight        = 101.070000)

Rh = Element(
    name                 = "Rhodium",
    symbol               = "Rh",
    atomic_number        = 45,
    atomic_weight        = 102.905500)

Pd = Element(
    name                 = "Palladium",
    symbol               = "Pd",
    atomic_number        = 46,
    atomic_weight        = 106.420000,
    van_der_waals_radius = 1.63)

Ag = Element(
    name                 = "Silver",
    symbol               = "Ag",
    atomic_number        = 47,
    atomic_weight        = 107.868200,
    van_der_waals_radius = 1.72)

Cd = Element(
    name                 = "Cadmium",
    symbol               = "Cd",
    atomic_number        = 48,
    atomic_weight        = 112.411000,
    van_der_waals_radius = 1.58)

In = Element(
    name                 = "Indium",
    symbol               = "In",
    atomic_number        = 49,
    atomic_weight        = 114.818000)

In = Element(
    name                 = "Indium",
    symbol               = "In",
    atomic_number        = 49,
    atomic_weight        = 114.818000,
    van_der_waals_radius = 1.93)

Sn = Element(
    name                 = "Tin",
    symbol               = "Sn",
    atomic_number        = 50,
    atomic_weight        = 118.710000,
    van_der_waals_radius = 2.17)

Sb = Element(
    name                 = "Antimony",
    symbol               = "Sb",
    atomic_number        = 51,
    atomic_weight        = 121.760000)

Te = Element(
    name                 = "Tellurium",
    symbol               = "Te",
    atomic_number        = 52,
    atomic_weight        = 127.600000,
    van_der_waals_radius = 2.06)

I  = Element(
    name                 = "Iodine",
    symbol               = "I",
    atomic_number        = 53,
    atomic_weight        = 126.904470,
    van_der_waals_radius = 1.98)

Xe = Element(
    name                 = "Xenon",
    symbol               = "Xe",
    atomic_number        = 54,
    atomic_weight        = 131.293000,
    van_der_waals_radius = 2.16)

Cs = Element(
    name                 = "Caesium",
    symbol               = "Cs",
    atomic_number        = 55,
    atomic_weight        = 132.905450)

Ba = Element(
    name                 = "Barium",
    symbol               = "Ba",
    atomic_number        = 56,
    atomic_weight        = 137.327000)

La = Element(
    name                 = "Lanthanum",
    symbol               = "La",
    atomic_number        = 57,
    atomic_weight        = 138.905500)

Ce = Element(
    name                 = "Cerium",
    symbol               = "Ce",
    atomic_number        = 58,
    atomic_weight        = 140.116000)

Pr = Element(
    name                 = "Praseodymium",
    symbol               = "Pr",
    atomic_number        = 59,
    atomic_weight        = 140.907650)

Nd = Element(
    name                 = "Neodymium",
    symbol               = "Nd",
    atomic_number        = 60,
    atomic_weight        = 144.240000)

Pm = Element(
    name                 = "Promethium",
    symbol               = "Pm",
    atomic_number        = 61,
    atomic_weight        = 145.000000)

Sm = Element(
    name                 = "Samarium",
    symbol               = "Sm",
    atomic_number        = 62,
    atomic_weight        = 150.360000)

Eu = Element(
    name                 = "Europium",
    symbol               = "Eu",
    atomic_number        = 63,
    atomic_weight        = 151.964000)

Gd = Element(
    name                 = "Gadolinium",
    symbol               = "Gd",
    atomic_number        = 64,
    atomic_weight        = 157.250000)

Tb = Element(
    name                 = "Terbium",
    symbol               = "Tb",
    atomic_number        = 65,
    atomic_weight        = 158.925340)

Dy = Element(
    name                 = "Dysprosium",
    symbol               = "Dy",
    atomic_number        = 66,
    atomic_weight        = 162.500000)

Ho = Element(
    name                 = "Holmium",
    symbol               = "Ho",
    atomic_number        = 67,
    atomic_weight        = 164.930320)

Er = Element(
    name                 = "Erbium",
    symbol               = "Er",
    atomic_number        = 68,
    atomic_weight        = 167.259000)

Tm = Element(
    name                 = "Thulium",
    symbol               = "Tm",
    atomic_number        = 69,
    atomic_weight        = 168.934210)

Yb = Element(
    name                 = "Ytterbium",
    symbol               = "Yb",
    atomic_number        = 70,
    atomic_weight        = 173.040000)

Lu = Element(
    name                 = "Lutetium",
    symbol               = "Lu",
    atomic_number        = 71,
    atomic_weight        = 174.967000)

Hf = Element(
    name                 = "Hafnium",
    symbol               = "Hf",
    atomic_number        = 72,
    atomic_weight        = 178.490000)

Ta = Element(
    name                 = "Tantalum",
    symbol               = "Ta",
    atomic_number        = 73,
    atomic_weight        = 180.947900)

W  = Element(
    name                 = "Tungsten",
    symbol               = "W",
    atomic_number        = 74,
    atomic_weight        = 183.840000)

Re = Element(
    name                 = "Rhenium",
    symbol               = "Re",
    atomic_number        = 75,
    atomic_weight        = 186.207000)

Os = Element(
    name                 = "Osmium",
    symbol               = "Os",
    atomic_number        = 76,
    atomic_weight        = 190.230000)

Ir = Element(
    name                 = "Iridium",
    symbol               = "Ir",
    atomic_number        = 77,
    atomic_weight        = 192.217000)

Pt = Element(
    name                 = "Platinum",
    symbol               = "Pt",
    atomic_number        = 78,
    atomic_weight        = 195.078000,
    van_der_waals_radius = 1.72)

Au = Element(
    name                 = "Gold",
    symbol               = "Au",
    atomic_number        = 79,
    atomic_weight        = 196.966550,
    van_der_waals_radius = 1.66)

Hg = Element(
    name                 = "Mercury",
    symbol               = "Hg",
    atomic_number        = 80,
    atomic_weight        = 200.590000,
    van_der_waals_radius = 1.55)

Tl = Element(
    name                 = "Thallium",
    symbol               = "Tl",
    atomic_number        = 81,
    atomic_weight        = 204.383300,
    van_der_waals_radius = 1.96)

Pb = Element(
    name                 = "Lead",
    symbol               = "Pb",
    atomic_number        = 82,
    atomic_weight        = 207.200000,
    van_der_waals_radius = 2.02)

Bi = Element(
    name                 = "Bismuth",
    symbol               = "Bi",
    atomic_number        = 83,
    atomic_weight        = 208.980380)

Po = Element(
    name                 = "Polonium",
    symbol               = "Po",
    atomic_number        = 84,
    atomic_weight        = 209.000000)

At = Element(
    name                 = "Astatine",
    symbol               = "At",
    atomic_number        = 85,
    atomic_weight        = 210.000000)

Rn = Element(
    name                 = "Radon",
    symbol               = "Rn",
    atomic_number        = 86,
    atomic_weight        = 222.000000)

Fr = Element(
    name                 = "Francium",
    symbol               = "Fr",
    atomic_number        = 87,
    atomic_weight        = 223.000000)

Ra = Element(
    name                 = "Radium",
    symbol               = "Ra",
    atomic_number        = 88,
    atomic_weight        = 226.000000)

Ac = Element(
    name                 = "Actinium",
    symbol               = "Ac",
    atomic_number        = 89,
    atomic_weight        = 227.000000)

Th = Element(
    name                 = "Thorium",
    symbol               = "Th",
    atomic_number        = 90,
    atomic_weight        = 232.038100)

Pa = Element(
    name                 = "Protactinium",
    symbol               = "Pa",
    atomic_number        = 91,
    atomic_weight        = 231.035880)

U  = Element(
    name                 = "Uranium",
    symbol               = "U",
    atomic_number        = 92,
    atomic_weight        = 238.028910,
    van_der_waals_radius = 1.86)


## this map includes upper-case versions of the element strings
ElementMap = {
    "H"   : H,
    "He"  : He,
    "HE"  : He,
    "Li"  : Li,
    "LI"  : Li,
    "Be"  : Be,
    "BE"  : Be,
    "B"   : B,
    "C"   : C,
    "N"   : N,
    "O"   : O,
    "F"   : F,
    "Ne"  : Ne,
    "NE"  : Ne,
    "Na"  : Na,
    "NA"  : Na,
    "Mg"  : Mg,
    "MG"  : Mg,
    "Al"  : Al,
    "AL"  : Al,
    "Si"  : Si,
    "SI"  : Si,
    "P"   : P,
    "S"   : S,
    "Cl"  : Cl,
    "CL"  : Cl,
    "Ar"  : Ar,
    "AR"  : Ar,
    "K"   : K,
    "Ca"  : Ca,
    "CA"  : Ca,
    "Sc"  : Sc,
    "SC"  : Sc,
    "Ti"  : Ti,
    "TI"  : Ti,
    "V"   : V,
    "Cr"  : Cr,
    "CR"  : Cr,
    "Mn"  : Mn,
    "MN"  : Mn,
    "Fe"  : Fe,
    "FE"  : Fe,
    "Co"  : Co,
    "CO"  : Co,
    "Ni"  : Ni,
    "NI"  : Ni,
    "Cu"  : Cu,
    "CU"  : Cu,
    "Zn"  : Zn,
    "ZN"  : Zn,
    "Ga"  : Ga,
    "GA"  : Ga,
    "Ge"  : Ge,
    "GE"  : Ge,
    "As"  : As,
    "AS"  : As,
    "Se"  : Se,
    "SE"  : Se,
    "Br"  : Br,
    "BR"  : Br,
    "Kr"  : Kr,
    "KR"  : Kr,
    "Rb"  : Rb,
    "RB"  : Rb,
    "Sr"  : Sr,
    "SR"  : Sr,
    "Y"   : Y,
    "Zr"  : Zr,
    "ZR"  : Zr,
    "Nb"  : Nb,
    "NB"  : Nb,
    "Mo"  : Mo,
    "MO"  : Mo,
    "Tc"  : Tc,
    "TC"  : Tc,
    "Ru"  : Ru,
    "RU"  : Ru,
    "Rh"  : Rh,
    "RH"  : Rh,
    "Pd"  : Pd,
    "PD"  : Pd,
    "Ag"  : Ag,
    "AG"  : Ag,
    "Cd"  : Cd,
    "CD"  : Cd,
    "In"  : In,
    "IN"  : In,
    "Sn"  : Sn,
    "SN"  : Sn,
    "Sb"  : Sb,
    "SB"  : Sb,
    "Te"  : Te,
    "TE"  : Te,
    "I"   : I,
    "Xe"  : Xe,
    "XE"  : Xe,
    "Cs"  : Cs,
    "CS"  : Cs,
    "Ba"  : Ba,
    "BA"  : Ba,
    "La"  : La,
    "LA"  : La,
    "Ce"  : Ce,
    "CE"  : Ce,
    "Pr"  : Pr,
    "PR"  : Pr,
    "Nd"  : Nd,
    "ND"  : Nd,
    "Pm"  : Pm,
    "PM"  : Pm,
    "Sm"  : Sm,
    "SM"  : Sm,
    "Eu"  : Eu,
    "EU"  : Eu,
    "Gd"  : Gd,
    "GD"  : Gd,
    "Tb"  : Tb,
    "TB"  : Tb,
    "Dy"  : Dy,
    "DY"  : Dy,
    "Ho"  : Ho,
    "HO"  : Ho,
    "Er"  : Er,
    "ER"  : Er,
    "Tm"  : Tm,
    "TM"  : Tm,
    "Yb"  : Yb,
    "YB"  : Yb,
    "Lu"  : Lu,
    "LU"  : Lu,
    "Hf"  : Hf,
    "HF"  : Hf,
    "Ta"  : Ta,
    "TA"  : Ta,
    "W"   : W,
    "Re"  : Re,
    "RE"  : Re,
    "Os"  : Os,
    "OS"  : Os,
    "Ir"  : Ir,
    "IR"  : Ir,
    "Pt"  : Pt,
    "PT"  : Pt,
    "Au"  : Au,
    "AU"  : Au,
    "Hg"  : Hg,
    "HG"  : Hg,
    "Tl"  : Tl,
    "TL"  : Tl,
    "Pb"  : Pb,
    "PB"  : Pb,
    "Bi"  : Bi,
    "BI"  : Bi,
    "Po"  : Po,
    "PO"  : Po,
    "At"  : At,
    "AT"  : At,
    "Rn"  : Rn,
    "RN"  : Rn,
    "Fr"  : Fr,
    "FR"  : Fr,
    "Ra"  : Ra,
    "RA"  : Ra,
    "Ac"  : Ac,
    "AC"  : Ac,
    "Th"  : Th,
    "TH"  : Th,
    "Pa"  : Pa,
    "PA"  : Pa,
    "U"   : U }


## converted to mmCIF
def new_doc():
    cif_file = mmCIFFile()

    dx    = {}
    elist = []

    for e in ElementMap.values():
        if not dx.has_key(e.name):
            dx[e.name] = True
            elist.append((e.atomic_number, e))

    elist.sort()

    for (n, e) in elist:
        cif_data = mmCIFData(e.symbol)
        cif_file.append(cif_data)

        cif_table = mmCIFTable(
            "element",
            ["symbol", "name", "number", "atomic_weight",
             "van_der_walls_radius", "color_rgb"])
        cif_data.append(cif_table)

        cif_table["symbol"] = e.symbol
        cif_table["name"] = e.name
        cif_table["number"] = str(e.atomic_number)
        cif_table["atomic_weight"] = "%.6f"%(e.atomic_weight)
        cif_table["van_der_walls_radius"] = "%.6f"%(e.van_der_waals_radius)

        if e.color != (1.0, 1.0, 1.0):
            cx = "#"
            for x in e.color:
                if x == 1.0: cx += "FF"
                else:        cx += "00"

            cif_table["color_rgb"] = cx
        else:
            cif_table["color_rgb"] = "#FFFFFF"

    return cif_file


def add_CCP4_atomsf(cif_file, path):
    """Adds the scattering factor data from the CCP4 file lib/data/atomsf.lib.
    """
    def get_element(tagName, attr, val):
        for x in doc.getElementsByTagName(tagName):
            if x.getAttribute(attr) == val:
                return x
        sys.stderr.write("Can't find %s:%s:%s\n" % (tagName, attr, val))
        return None

    state = -1
    ename = ""
    eion  = ""
    eelm  = None
    elsf  = None

    charge = "0"

    fil = open(path, "r")
    for ln in fil.readlines():
        ln = ln.rstrip()
        if ln == "":
            continue
        if ln[:2] == "AD":
            continue

        listx = ln.split()
        if len(listx) == 1:
            state = 0

        if state == 0:
            ename = listx[0]
            ## check if the element is a ion
            i = max(ename.find("+"), ename.find("-"))
            if i >= 0:
                charge = ename[i:]
                ename  = ename[:i]

                try:
                    elem_cif_data = cif_file[ename]
                except KeyError:
                    state = -1
                    continue

            else:
                charge = "0"

                try:
                    elem_cif_data = cif_file[ename]
                except KeyError:
                    state = -1
                    continue

            state += 1

        elif state == 1:
            (eiwt, eielec, ec) = listx

            ## add the <ScatteringFactors> section
            try:
                sf_cif_table = elem_cif_data["scattering_factors"]
            except KeyError:
                sf_cif_table = mmCIFTable(
                    "scattering_factors",
                    ["wavelength", "charge", "a1","a2","a3","a4",
                     "b1","b2","b3","b4", "c"])
                elem_cif_data.append(sf_cif_table)

            sf_cif_row = mmCIFRow()
            sf_cif_table.append(sf_cif_row)
            sf_cif_row["wavelength"] = "CuKa"
            sf_cif_row["charge"] = charge
            sf_cif_row["c"] = str(ec)

            state += 1

        elif state == 2:
            ## add A1, A2, A3, A4
            for i in range(len(listx)):
                column = "a%d" % (i+1)
                sf_cif_row[column] = str(listx[i])

            state += 1

        elif state == 3:
            ## add B1, B2, B3, B4
            for i in range(len(listx)):
                column = "b%d" % (i+1)
                sf_cif_row[column] = str(listx[i])

            state += 1

        elif state == 4:
            try:
                asf_cif_table = elem_cif_data["anomalous_scattering_factors"]
            except KeyError:
                asf_cif_table = mmCIFTable(
                    "anomalous_scattering_factors",
                    ["wavelength", "charge", "dF", "d2F"])
                elem_cif_data.append(asf_cif_table)

            asf_cif_row = mmCIFRow()
            asf_cif_table.append(asf_cif_row)
            asf_cif_row["wavelength"] = "Cu"
            asf_cif_row["charge"] = charge
            asf_cif_row["dF"] = str(listx[0])
            asf_cif_row["d2F"] = str(listx[1])

            asf_cif_row = mmCIFRow()
            asf_cif_table.append(asf_cif_row)
            asf_cif_row["wavelength"] = "Mo"
            asf_cif_row["charge"] = charge
            asf_cif_row["dF"] = str(listx[2])
            asf_cif_row["d2F"] = str(listx[3])

            state += 1


if __name__ == "__main__":
    try:
        path = sys.argv[1]
    except IndexError:
        doc = new_doc()
    else:
        doc = parse(open(sys.argv[1], "r"))

    ## add scattering factors
    ## NOTE: The 'atomsf.lib' file can be found in the lib/data/ directory of
    ## the latest CCP4 source tarball
    add_CCP4_atomsf(doc, "atomsf.lib")

    doc.save_file(fil=sys.stdout)
