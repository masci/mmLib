## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.


ElementNames = [
    'H', 'He', 'HE', 'Li', 'LI', 'Be', 'BE', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'NE', 'Na', 'NA', 'Mg', 'MG', 'Al', 'AL', 'Si', 'SI', 'P', 'S', 'Cl',
    'CL', 'Ar', 'AR', 'K', 'Ca', 'CA', 'Sc', 'SC', 'Ti', 'TI', 'V', 'Cr',
    'CR', 'Mn', 'MN', 'Fe', 'FE', 'Co', 'CO', 'Ni', 'NI', 'Cu', 'CU', 'Zn',
    'ZN', 'Ga', 'GA', 'Ge', 'GE', 'As', 'AS', 'Se', 'SE', 'Br', 'BR', 'Kr',
    'KR', 'Rb', 'RB', 'Sr', 'SR', 'Y', 'Zr', 'ZR', 'Nb', 'NB', 'Mo', 'MO',
    'Tc', 'TC', 'Ru', 'RU', 'Rh', 'RH', 'Pd', 'PD', 'Ag', 'AG', 'Cd', 'CD',
    'In', 'IN', 'Sn', 'SN', 'Sb', 'SB', 'Te', 'TE', 'I', 'Xe', 'XE', 'Cs',
    'CS', 'Ba', 'BA', 'La', 'LA', 'Ce', 'CE', 'Pr', 'PR', 'Nd', 'ND', 'Pm',
    'PM', 'Sm', 'SM', 'Eu', 'EU', 'Gd', 'GD', 'Tb', 'TB', 'Dy', 'DY', 'Ho',
    'HO', 'Er', 'ER', 'Tm', 'TM', 'Yb', 'YB', 'Lu', 'LU', 'Hf', 'HF', 'Ta',
    'TA', 'W', 'Re', 'RE', 'Os', 'OS', 'Ir', 'IR', 'Pt', 'PT', 'Au', 'AU',
    'Hg', 'HG', 'Tl', 'TL', 'Pb', 'PB', 'Bi', 'BI', 'Po', 'PO', 'At', 'AT',
    'Rn', 'RN', 'Fr', 'FR', 'Ra', 'RA', 'Ac', 'AC', 'Th', 'TH', 'Pa', 'PA',
    'U', 'Np', 'NP', 'Pu', 'PU', 'Am', 'AM', 'Cm', 'CM', 'Bk', 'BK', 'Cf',
    'CF', 'Es', 'ES', 'Fm', 'FM', 'Md', 'MD', 'No', 'NO', 'Lr', 'LR', 'Rf',
    'RF', 'Db', 'DB', 'Sg', 'SG', 'Bh', 'BH', 'Hs', 'HS', 'Mt', 'MT'
    ]


class Element:
    def __init__(self):
        self.name = None
        self.symbol = None
        self.group = None
        self.period = None
        self.atomic_number = None
        self.atomic_weight = None
        self.atomic_radius = None
        self.covalent_radius = None
        self.van_der_waals_radius = None
        self.electronegativity = None

    def __str__(self):
        return "Element=%s" % (self.name)


H = Element()
H.name = "Hydrogen"
H.symbol = "H"
H.atomic_number = 1
H.atomic_mass = 1.007940

He = Element()
He.name = "Helium"
He.symbol = "He"
He.atomic_number = 2
He.atomic_mass = 4.002602

Li = Element()
Li.name = "Lithium"
Li.symbol = "Li"
Li.atomic_number = 3
Li.atomic_mass = 6.941000

Be = Element()
Be.name = "Beryllium"
Be.symbol = "Be"
Be.atomic_number = 4
Be.atomic_mass = 9.012182

B = Element()
B.name = "Boron"
B.symbol = "B"
B.atomic_number = 5
B.atomic_mass = 10.811000

C = Element()
C.name = "Carbon"
C.symbol = "C"
C.atomic_number = 6
C.atomic_mass = 12.010700

N = Element()
N.name = "Nitrogen"
N.symbol = "N"
N.atomic_number = 7
N.atomic_mass = 14.006700

O = Element()
O.name = "Oxygen"
O.symbol = "O"
O.atomic_number = 8
O.atomic_mass = 15.999400

F = Element()
F.name = "Fluorine"
F.symbol = "F"
F.atomic_number = 9
F.atomic_mass = 18.998403

Ne = Element()
Ne.name = "Neon"
Ne.symbol = "Ne"
Ne.atomic_number = 10
Ne.atomic_mass = 20.179700

Na = Element()
Na.name = "Sodium"
Na.symbol = "Na"
Na.atomic_number = 11
Na.atomic_mass = 22.989770

Mg = Element()
Mg.name = "Magnesium"
Mg.symbol = "Mg"
Mg.atomic_number = 12
Mg.atomic_mass = 24.305000

Al = Element()
Al.name = "Aluminium"
Al.symbol = "Al"
Al.atomic_number = 13
Al.atomic_mass = 26.981538

Si = Element()
Si.name = "Silicon"
Si.symbol = "Si"
Si.atomic_number = 14
Si.atomic_mass = 28.085500

P = Element()
P.name = "Phosphorus"
P.symbol = "P"
P.atomic_number = 15
P.atomic_mass = 30.973761

S = Element()
S.name = "Sulfur"
S.symbol = "S"
S.atomic_number = 16
S.atomic_mass = 32.065000

Cl = Element()
Cl.name = "Chlorine"
Cl.symbol = "Cl"
Cl.atomic_number = 17
Cl.atomic_mass = 35.453000

Ar = Element()
Ar.name = "Argon"
Ar.symbol = "Ar"
Ar.atomic_number = 18
Ar.atomic_mass = 39.948000

K = Element()
K.name = "Potassium"
K.symbol = "K"
K.atomic_number = 19
K.atomic_mass = 39.098300

Ca = Element()
Ca.name = "Calcium"
Ca.symbol = "Ca"
Ca.atomic_number = 20
Ca.atomic_mass = 40.078000

Sc = Element()
Sc.name = "Scandium"
Sc.symbol = "Sc"
Sc.atomic_number = 21
Sc.atomic_mass = 44.955910

Ti = Element()
Ti.name = "Titanium"
Ti.symbol = "Ti"
Ti.atomic_number = 22
Ti.atomic_mass = 47.867000

V = Element()
V.name = "Vanadium"
V.symbol = "V"
V.atomic_number = 23
V.atomic_mass = 50.941500

Cr = Element()
Cr.name = "Chromium"
Cr.symbol = "Cr"
Cr.atomic_number = 24
Cr.atomic_mass = 51.996100

Mn = Element()
Mn.name = "Manganese"
Mn.symbol = "Mn"
Mn.atomic_number = 25
Mn.atomic_mass = 54.938049

Fe = Element()
Fe.name = "Iron"
Fe.symbol = "Fe"
Fe.atomic_number = 26
Fe.atomic_mass = 55.845000

Co = Element()
Co.name = "Cobalt"
Co.symbol = "Co"
Co.atomic_number = 27
Co.atomic_mass = 58.933200

Ni = Element()
Ni.name = "Nickel"
Ni.symbol = "Ni"
Ni.atomic_number = 28
Ni.atomic_mass = 58.693400

Cu = Element()
Cu.name = "Copper"
Cu.symbol = "Cu"
Cu.atomic_number = 29
Cu.atomic_mass = 63.546000

Zn = Element()
Zn.name = "Zinc"
Zn.symbol = "Zn"
Zn.atomic_number = 30
Zn.atomic_mass = 65.409000

Ga = Element()
Ga.name = "Gallium"
Ga.symbol = "Ga"
Ga.atomic_number = 31
Ga.atomic_mass = 69.723000

Ge = Element()
Ge.name = "Germanium"
Ge.symbol = "Ge"
Ge.atomic_number = 32
Ge.atomic_mass = 72.640000

As = Element()
As.name = "Arsenic"
As.symbol = "As"
As.atomic_number = 33
As.atomic_mass = 74.921600

Se = Element()
Se.name = "Selenium"
Se.symbol = "Se"
Se.atomic_number = 34
Se.atomic_mass = 78.960000

Br = Element()
Br.name = "Bromine"
Br.symbol = "Br"
Br.atomic_number = 35
Br.atomic_mass = 79.904000

Kr = Element()
Kr.name = "Krypton"
Kr.symbol = "Kr"
Kr.atomic_number = 36
Kr.atomic_mass = 83.798000

Rb = Element()
Rb.name = "Rubidium"
Rb.symbol = "Rb"
Rb.atomic_number = 37
Rb.atomic_mass = 85.467800

Sr = Element()
Sr.name = "Strontium"
Sr.symbol = "Sr"
Sr.atomic_number = 38
Sr.atomic_mass = 87.620000

Y = Element()
Y.name = "Yttrium"
Y.symbol = "Y"
Y.atomic_number = 39
Y.atomic_mass = 88.905850

Zr = Element()
Zr.name = "Zirconium"
Zr.symbol = "Zr"
Zr.atomic_number = 40
Zr.atomic_mass = 91.224000

Nb = Element()
Nb.name = "Niobium"
Nb.symbol = "Nb"
Nb.atomic_number = 41
Nb.atomic_mass = 92.906380

Mo = Element()
Mo.name = "Molybdenum"
Mo.symbol = "Mo"
Mo.atomic_number = 42
Mo.atomic_mass = 95.940000

Tc = Element()
Tc.name = "Technetium"
Tc.symbol = "Tc"
Tc.atomic_number = 43
Tc.atomic_mass = 98.000000

Ru = Element()
Ru.name = "Ruthenium"
Ru.symbol = "Ru"
Ru.atomic_number = 44
Ru.atomic_mass = 101.070000

Rh = Element()
Rh.name = "Rhodium"
Rh.symbol = "Rh"
Rh.atomic_number = 45
Rh.atomic_mass = 102.905500

Pd = Element()
Pd.name = "Palladium"
Pd.symbol = "Pd"
Pd.atomic_number = 46
Pd.atomic_mass = 106.420000

Ag = Element()
Ag.name = "Silver"
Ag.symbol = "Ag"
Ag.atomic_number = 47
Ag.atomic_mass = 107.868200

Cd = Element()
Cd.name = "Cadmium"
Cd.symbol = "Cd"
Cd.atomic_number = 48
Cd.atomic_mass = 112.411000

In = Element()
In.name = "Indium"
In.symbol = "In"
In.atomic_number = 49
In.atomic_mass = 114.818000

Sn = Element()
Sn.name = "Tin"
Sn.symbol = "Sn"
Sn.atomic_number = 50
Sn.atomic_mass = 118.710000

Sb = Element()
Sb.name = "Antimony"
Sb.symbol = "Sb"
Sb.atomic_number = 51
Sb.atomic_mass = 121.760000

Te = Element()
Te.name = "Tellurium"
Te.symbol = "Te"
Te.atomic_number = 52
Te.atomic_mass = 127.600000

I = Element()
I.name = "Iodine"
I.symbol = "I"
I.atomic_number = 53
I.atomic_mass = 126.904470

Xe = Element()
Xe.name = "Xenon"
Xe.symbol = "Xe"
Xe.atomic_number = 54
Xe.atomic_mass = 131.293000

Cs = Element()
Cs.name = "Caesium"
Cs.symbol = "Cs"
Cs.atomic_number = 55
Cs.atomic_mass = 132.905450

Ba = Element()
Ba.name = "Barium"
Ba.symbol = "Ba"
Ba.atomic_number = 56
Ba.atomic_mass = 137.327000

La = Element()
La.name = "Lanthanum"
La.symbol = "La"
La.atomic_number = 57
La.atomic_mass = 138.905500

Ce = Element()
Ce.name = "Cerium"
Ce.symbol = "Ce"
Ce.atomic_number = 58
Ce.atomic_mass = 140.116000

Pr = Element()
Pr.name = "Praseodymium"
Pr.symbol = "Pr"
Pr.atomic_number = 59
Pr.atomic_mass = 140.907650

Nd = Element()
Nd.name = "Neodymium"
Nd.symbol = "Nd"
Nd.atomic_number = 60
Nd.atomic_mass = 144.240000

Pm = Element()
Pm.name = "Promethium"
Pm.symbol = "Pm"
Pm.atomic_number = 61
Pm.atomic_mass = 145.000000

Sm = Element()
Sm.name = "Samarium"
Sm.symbol = "Sm"
Sm.atomic_number = 62
Sm.atomic_mass = 150.360000

Eu = Element()
Eu.name = "Europium"
Eu.symbol = "Eu"
Eu.atomic_number = 63
Eu.atomic_mass = 151.964000

Gd = Element()
Gd.name = "Gadolinium"
Gd.symbol = "Gd"
Gd.atomic_number = 64
Gd.atomic_mass = 157.250000

Tb = Element()
Tb.name = "Terbium"
Tb.symbol = "Tb"
Tb.atomic_number = 65
Tb.atomic_mass = 158.925340

Dy = Element()
Dy.name = "Dysprosium"
Dy.symbol = "Dy"
Dy.atomic_number = 66
Dy.atomic_mass = 162.500000

Ho = Element()
Ho.name = "Holmium"
Ho.symbol = "Ho"
Ho.atomic_number = 67
Ho.atomic_mass = 164.930320

Er = Element()
Er.name = "Erbium"
Er.symbol = "Er"
Er.atomic_number = 68
Er.atomic_mass = 167.259000

Tm = Element()
Tm.name = "Thulium"
Tm.symbol = "Tm"
Tm.atomic_number = 69
Tm.atomic_mass = 168.934210

Yb = Element()
Yb.name = "Ytterbium"
Yb.symbol = "Yb"
Yb.atomic_number = 70
Yb.atomic_mass = 173.040000

Lu = Element()
Lu.name = "Lutetium"
Lu.symbol = "Lu"
Lu.atomic_number = 71
Lu.atomic_mass = 174.967000

Hf = Element()
Hf.name = "Hafnium"
Hf.symbol = "Hf"
Hf.atomic_number = 72
Hf.atomic_mass = 178.490000

Ta = Element()
Ta.name = "Tantalum"
Ta.symbol = "Ta"
Ta.atomic_number = 73
Ta.atomic_mass = 180.947900

W = Element()
W.name = "Tungsten"
W.symbol = "W"
W.atomic_number = 74
W.atomic_mass = 183.840000

Re = Element()
Re.name = "Rhenium"
Re.symbol = "Re"
Re.atomic_number = 75
Re.atomic_mass = 186.207000

Os = Element()
Os.name = "Osmium"
Os.symbol = "Os"
Os.atomic_number = 76
Os.atomic_mass = 190.230000

Ir = Element()
Ir.name = "Iridium"
Ir.symbol = "Ir"
Ir.atomic_number = 77
Ir.atomic_mass = 192.217000

Pt = Element()
Pt.name = "Platinum"
Pt.symbol = "Pt"
Pt.atomic_number = 78
Pt.atomic_mass = 195.078000

Au = Element()
Au.name = "Gold"
Au.symbol = "Au"
Au.atomic_number = 79
Au.atomic_mass = 196.966550

Hg = Element()
Hg.name = "Mercury"
Hg.symbol = "Hg"
Hg.atomic_number = 80
Hg.atomic_mass = 200.590000

Tl = Element()
Tl.name = "Thallium"
Tl.symbol = "Tl"
Tl.atomic_number = 81
Tl.atomic_mass = 204.383300

Pb = Element()
Pb.name = "Lead"
Pb.symbol = "Pb"
Pb.atomic_number = 82
Pb.atomic_mass = 207.200000

Bi = Element()
Bi.name = "Bismuth"
Bi.symbol = "Bi"
Bi.atomic_number = 83
Bi.atomic_mass = 208.980380

Po = Element()
Po.name = "Polonium"
Po.symbol = "Po"
Po.atomic_number = 84
Po.atomic_mass = 209.000000

At = Element()
At.name = "Astatine"
At.symbol = "At"
At.atomic_number = 85
At.atomic_mass = 210.000000

Rn = Element()
Rn.name = "Radon"
Rn.symbol = "Rn"
Rn.atomic_number = 86
Rn.atomic_mass = 222.000000

Fr = Element()
Fr.name = "Francium"
Fr.symbol = "Fr"
Fr.atomic_number = 87
Fr.atomic_mass = 223.000000

Ra = Element()
Ra.name = "Radium"
Ra.symbol = "Ra"
Ra.atomic_number = 88
Ra.atomic_mass = 226.000000

Ac = Element()
Ac.name = "Actinium"
Ac.symbol = "Ac"
Ac.atomic_number = 89
Ac.atomic_mass = 227.000000

Th = Element()
Th.name = "Thorium"
Th.symbol = "Th"
Th.atomic_number = 90
Th.atomic_mass = 232.038100

Pa = Element()
Pa.name = "Protactinium"
Pa.symbol = "Pa"
Pa.atomic_number = 91
Pa.atomic_mass = 231.035880

U = Element()
U.name = "Uranium"
U.symbol = "U"
U.atomic_number = 92
U.atomic_mass = 238.028910

Np = Element()
Np.name = "Neptunium"
Np.symbol = "Np"
Np.atomic_number = 93
Np.atomic_mass = 237.000000

Pu = Element()
Pu.name = "Plutonium"
Pu.symbol = "Pu"
Pu.atomic_number = 94
Pu.atomic_mass = 244.000000

Am = Element()
Am.name = "Americium"
Am.symbol = "Am"
Am.atomic_number = 95
Am.atomic_mass = 243.000000

Cm = Element()
Cm.name = "Curium"
Cm.symbol = "Cm"
Cm.atomic_number = 96
Cm.atomic_mass = 247.000000

Bk = Element()
Bk.name = "Berkelium"
Bk.symbol = "Bk"
Bk.atomic_number = 97
Bk.atomic_mass = 247.000000

Cf = Element()
Cf.name = "Californium"
Cf.symbol = "Cf"
Cf.atomic_number = 98
Cf.atomic_mass = 251.000000

Es = Element()
Es.name = "Einsteinium"
Es.symbol = "Es"
Es.atomic_number = 99
Es.atomic_mass = 252.000000

Fm = Element()
Fm.name = "Fermium"
Fm.symbol = "Fm"
Fm.atomic_number = 100
Fm.atomic_mass = 257.000000

Md = Element()
Md.name = "Mendelevium"
Md.symbol = "Md"
Md.atomic_number = 101
Md.atomic_mass = 258.000000

No = Element()
No.name = "Nobelium"
No.symbol = "No"
No.atomic_number = 102
No.atomic_mass = 259.000000

Lr = Element()
Lr.name = "Lawrencium"
Lr.symbol = "Lr"
Lr.atomic_number = 103
Lr.atomic_mass = 262.000000

Rf = Element()
Rf.name = "Rutherfordium"
Rf.symbol = "Rf"
Rf.atomic_number = 104
Rf.atomic_mass = 261.000000

Db = Element()
Db.name = "Dubnium"
Db.symbol = "Db"
Db.atomic_number = 105
Db.atomic_mass = 262.000000

Sg = Element()
Sg.name = "Seaborgium"
Sg.symbol = "Sg"
Sg.atomic_number = 106
Sg.atomic_mass = 266.000000

Bh = Element()
Bh.name = "Bohrium"
Bh.symbol = "Bh"
Bh.atomic_number = 107
Bh.atomic_mass = 264.000000

Hs = Element()
Hs.name = "Hassium"
Hs.symbol = "Hs"
Hs.atomic_number = 108
Hs.atomic_mass = 277.000000

Mt = Element()
Mt.name = "Meitnerium"
Mt.symbol = "Mt"
Mt.atomic_number = 109
Mt.atomic_mass = 268.000000


## this map includes upper-case versions of the element strings
ElementMap = {
    "H" : H,
    "He" : He,
    "HE" : He,
    "Li" : Li,
    "LI" : Li,
    "Be" : Be,
    "BE" : Be,
    "B" : B,
    "C" : C,
    "N" : N,
    "O" : O,
    "F" : F,
    "Ne" : Ne,
    "NE" : Ne,
    "Na" : Na,
    "NA" : Na,
    "Mg" : Mg,
    "MG" : Mg,
    "Al" : Al,
    "AL" : Al,
    "Si" : Si,
    "SI" : Si,
    "P" : P,
    "S" : S,
    "Cl" : Cl,
    "CL" : Cl,
    "Ar" : Ar,
    "AR" : Ar,
    "K" : K,
    "Ca" : Ca,
    "CA" : Ca,
    "Sc" : Sc,
    "SC" : Sc,
    "Ti" : Ti,
    "TI" : Ti,
    "V" : V,
    "Cr" : Cr,
    "CR" : Cr,
    "Mn" : Mn,
    "MN" : Mn,
    "Fe" : Fe,
    "FE" : Fe,
    "Co" : Co,
    "CO" : Co,
    "Ni" : Ni,
    "NI" : Ni,
    "Cu" : Cu,
    "CU" : Cu,
    "Zn" : Zn,
    "ZN" : Zn,
    "Ga" : Ga,
    "GA" : Ga,
    "Ge" : Ge,
    "GE" : Ge,
    "As" : As,
    "AS" : As,
    "Se" : Se,
    "SE" : Se,
    "Br" : Br,
    "BR" : Br,
    "Kr" : Kr,
    "KR" : Kr,
    "Rb" : Rb,
    "RB" : Rb,
    "Sr" : Sr,
    "SR" : Sr,
    "Y" : Y,
    "Zr" : Zr,
    "ZR" : Zr,
    "Nb" : Nb,
    "NB" : Nb,
    "Mo" : Mo,
    "MO" : Mo,
    "Tc" : Tc,
    "TC" : Tc,
    "Ru" : Ru,
    "RU" : Ru,
    "Rh" : Rh,
    "RH" : Rh,
    "Pd" : Pd,
    "PD" : Pd,
    "Ag" : Ag,
    "AG" : Ag,
    "Cd" : Cd,
    "CD" : Cd,
    "In" : In,
    "IN" : In,
    "Sn" : Sn,
    "SN" : Sn,
    "Sb" : Sb,
    "SB" : Sb,
    "Te" : Te,
    "TE" : Te,
    "I" : I,
    "Xe" : Xe,
    "XE" : Xe,
    "Cs" : Cs,
    "CS" : Cs,
    "Ba" : Ba,
    "BA" : Ba,
    "La" : La,
    "LA" : La,
    "Ce" : Ce,
    "CE" : Ce,
    "Pr" : Pr,
    "PR" : Pr,
    "Nd" : Nd,
    "ND" : Nd,
    "Pm" : Pm,
    "PM" : Pm,
    "Sm" : Sm,
    "SM" : Sm,
    "Eu" : Eu,
    "EU" : Eu,
    "Gd" : Gd,
    "GD" : Gd,
    "Tb" : Tb,
    "TB" : Tb,
    "Dy" : Dy,
    "DY" : Dy,
    "Ho" : Ho,
    "HO" : Ho,
    "Er" : Er,
    "ER" : Er,
    "Tm" : Tm,
    "TM" : Tm,
    "Yb" : Yb,
    "YB" : Yb,
    "Lu" : Lu,
    "LU" : Lu,
    "Hf" : Hf,
    "HF" : Hf,
    "Ta" : Ta,
    "TA" : Ta,
    "W" : W,
    "Re" : Re,
    "RE" : Re,
    "Os" : Os,
    "OS" : Os,
    "Ir" : Ir,
    "IR" : Ir,
    "Pt" : Pt,
    "PT" : Pt,
    "Au" : Au,
    "AU" : Au,
    "Hg" : Hg,
    "HG" : Hg,
    "Tl" : Tl,
    "TL" : Tl,
    "Pb" : Pb,
    "PB" : Pb,
    "Bi" : Bi,
    "BI" : Bi,
    "Po" : Po,
    "PO" : Po,
    "At" : At,
    "AT" : At,
    "Rn" : Rn,
    "RN" : Rn,
    "Fr" : Fr,
    "FR" : Fr,
    "Ra" : Ra,
    "RA" : Ra,
    "Ac" : Ac,
    "AC" : Ac,
    "Th" : Th,
    "TH" : Th,
    "Pa" : Pa,
    "PA" : Pa,
    "U" : U,
    "Np" : Np,
    "NP" : Np,
    "Pu" : Pu,
    "PU" : Pu,
    "Am" : Am,
    "AM" : Am,
    "Cm" : Cm,
    "CM" : Cm,
    "Bk" : Bk,
    "BK" : Bk,
    "Cf" : Cf,
    "CF" : Cf,
    "Es" : Es,
    "ES" : Es,
    "Fm" : Fm,
    "FM" : Fm,
    "Md" : Md,
    "MD" : Md,
    "No" : No,
    "NO" : No,
    "Lr" : Lr,
    "LR" : Lr,
    "Rf" : Rf,
    "RF" : Rf,
    "Db" : Db,
    "DB" : Db,
    "Sg" : Sg,
    "SG" : Sg,
    "Bh" : Bh,
    "BH" : Bh,
    "Hs" : Hs,
    "HS" : Hs,
    "Mt" : Mt,
    "MT" : Mt }
