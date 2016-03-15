# -*- coding: utf-8 -*-
'''
This module contains atomic data and conversion factors.

Most of this data is from PyQuante.
https://github.com/rpmuller/pyquante2
'''
bohr2angstrom = 0.529177249
angstrom2meter = 1e-10
wavenumber2ghz = 29.9792458
# Planck's constant
# http://physics.nist.gov/cgi-bin/cuu/Value?h
planck = 6.626070040e-34  # Js

# Atomic mass unit
amu = 1.6605402e-27  # kg

mass = [
    0.00,
    1.0008, 4.0026,
    6.941, 9.0122,
    10.811, 12.011, 14.007, 15.999, 18.998, 20.179,
    22.990, 24.305,
    26.982, 28.086, 30.974, 32.066, 35.453, 39.948,
    39.098, 40.078,
    44.9559, 47.867, 50.9415, 51.9961, 54.938, 55.845,
    58.9332, 58.6934, 63.546, 65.39,
    69.723, 72.61, 74.9216, 78.96, 79.904, 83.80,
    85.4678, 87.62,
    88.90686, 91.224, 92.90638, 95.94, 98, 101.07,
    102.90550, 106.42, 107.8682, 112.411,
    114.818, 118.710, 121.760, 127.60, 126.90447, 131.29,
    132.90545, 137.327, 138.9055, 140.11, 140.90765, 144.24,
    145.0, 150.36, 151.964,
    157.25, 158.92534, 162.5, 164.93, 167.259, 168.934, 173.04, 174.967,
    178.49, 180.9479, 183.84, 186.207, 190.23, 192.217, 195.078, 196.96655,
    200.59]

symbol = [
    "X", "H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",  "Eu",
    "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn"]

sym2no = {}
for i in range(len(symbol)):
    sym2no[symbol[i]] = i
    sym2no[symbol[i].lower()] = i

name = [
    "dummy",
    "hydrogen", "helium",
    "lithium", "beryllium", "boron", "carbon", "nitrogen",
    "oxygen", "fluorine", "neon", "sodium", "magnesium",
    "aluminum", "silicon", "phosphorus", "sulfur", "chlorine",
    "argon", "potassium", "calcium", "scandium", "titanium",
    "vanadium", "chromium", "manganese", "iron",
    "cobalt", "nickel", "copper", "zinc",
    "gallium", "germanium", "arsenic", "selenium", "bromine",
    "krypton", "rubidium", "strontium", "yttrium", "zirconium",
    "niobium", "molybdenum", "technetium", "ruthenium", "rhodium",
    "palladium", "silver", "cadmium", "indium", "tin", "antimony",
    "tellerium", "iodine", "xenon", "cesium", "barium",
    "lanthanum", "cerium", "praseodymium", "neodymium", "promethium",
    "samarium", "europium", "gadolinium", "terbium", "dysprosium",
    "holmium", "erbium", "thulium", "ytterbium", "lutetium",
    "halfnium", "tantalum", "tungsten", "rhenium", "osmium", "iridium",
    "platinum", "gold", "mercury", "thallium", "lead", "bismuth",
    "polonium", "astatine", "radon"]

color = [
    (0, 0, 0),        # dummy
    (250, 235, 215),  # H,  1
    (255, 192, 203),  # He, 2
    (178, 34, 34),    # Li, 3
    (34, 139, 34),    # Be, 4
    (0, 255, 0),      # B,  5
    (112, 128, 144),  # C,  6
    (0, 191, 255),    # N,  7
    (255, 0, 0),      # O,  8
    (218, 165, 32),   # F,  9
    (255, 105, 180),  # Ne, 10
    (0, 0, 255),      # Na, 11
    (34, 139, 34),    # Mg, 12
    (190, 190, 190),  # Al, 13
    (218, 165, 32),   # Si, 14
    (255, 165, 0),    # P,  15
    (255, 255, 0),    # S,  16
    (0, 255, 0),      # Cl, 17
    (255, 192, 203),  # Ar, 18
    (255, 20, 147),   # K,  19
    (128, 128, 128),  # Ca, 20
    (190, 190, 190),  # Sc, 21
    (190, 190, 190),  # Ti, 22
    (190, 190, 190),  # V,  23
    (190, 190, 190),  # Cr, 24
    (190, 190, 190),  # Mn, 25
    (255, 165, 0),    # Fe, 26
    (165, 42, 42),    # Co, 27
    (165, 42, 42),    # Ni, 28
    (165, 42, 42),    # Cu, 29
    (165, 42, 42),    # Zn, 30
    (165, 42, 42),    # Ga, 31
    (85, 107, 47),    # Ge, 32
    (253, 245, 230),  # As, 33
    (152, 251, 152),  # Se, 34
    (165, 42, 42),    # Br, 35
    (50, 205, 50),    # Kr, 36
    (165, 42, 42),    # Rb, 37
    (190, 190, 190),  # Sr, 38
    (190, 190, 190),  # Y,  39
    (190, 190, 190),  # Zr, 40
    (190, 190, 190),  # Nb, 41
    (255, 127, 80),   # Mo, 42
    (190, 190, 190),  # Tc, 43
    (190, 190, 190),  # Ru, 44
    (190, 190, 190),  # Rh, 45
    (190, 190, 190),  # Pd, 46
    (190, 190, 190),  # Ag, 47
    (255, 140, 0),    # Cd, 48
    (190, 190, 190),  # In, 49
    (190, 190, 190),  # Sn, 50
    (190, 190, 190),  # Sb, 51
    (190, 190, 190),  # Te, 52
    (160, 32, 240),   # I,  53
    (255, 105, 180),  # Xe, 54
    (165, 42, 42),    # Cs, 55
    (190, 190, 190),  # Ba, 56
    (190, 190, 190),  # La, 57
    (190, 190, 190),  # Ce, 58
    (190, 190, 190),  # Pr, 59
    (190, 190, 190),  # Nd, 60
    (190, 190, 190),  # Pm, 61
    (190, 190, 190),  # Sm, 62
    (190, 190, 190),  # Eu, 63
    (190, 190, 190),  # Gd, 64
    (190, 190, 190),  # Tb, 65
    (190, 190, 190),  # Dy, 66
    (190, 190, 190),  # Ho, 67
    (190, 190, 190),  # Er, 68
    (190, 190, 190),  # Tm, 69
    (190, 190, 190),  # Yb, 70
    (190, 190, 190),  # Lu, 71
    (190, 190, 190),  # Hf, 72
    (190, 190, 190),  # Ta, 73
    (64, 224, 208),   # W,  74
    (190, 190, 190),  # Re, 75
    (190, 190, 190),  # Os, 76
    (190, 190, 190),  # Ir, 77
    (190, 190, 190),  # Pt, 78
    (255, 215, 0),    # Au, 79
    (190, 190, 190),  # Hg, 80
    (190, 190, 190),  # Tl, 81
    (190, 190, 190),  # Pb, 82
    (255, 181, 197),  # Bi, 83
    (190, 190, 190),  # Po, 84
    (190, 190, 190),  # At, 85
    (190, 190, 190),  # Rn, 86
    (190, 190, 190),  # Fr, 87
    (190, 190, 190),  # Ra, 88
    (190, 190, 190),  # Ac, 89
    (190, 190, 190),  # Th, 90
    (190, 190, 190),  # Pa, 91
    (90, 90, 90),     # U,  92
    (190, 190, 190),  # Np, 93
    (190, 190, 190),  # Pu, 94
    (190, 190, 190),  # Am, 95
    (190, 190, 190),  # Cm, 96
    (190, 190, 190),  # Bk, 97
    (190, 190, 190),  # Cf, 98
    (190, 190, 190),  # Es, 99
    (190, 190, 190),  # Fm, 100
    (190, 190, 190),  # Md, 101
    (190, 190, 190),  # No, 102
    (190, 190, 190),  # Lr, 103
    (190, 190, 190),  # Rf, 104
    (190, 190, 190),  # Db, 105
    (190, 190, 190),  # Sg, 106
    (190, 190, 190),  # Bh, 107
    (190, 190, 190),  # Hs, 108
    (190, 190, 190),  # Mt, 109
    (190, 190, 190),  # Ds, 110
    (190, 190, 190)]  # Rg, 111

floatcolor = [(col[0] / 255., col[1] / 255., col[2] / 255.) for col in color]

rvdw = [
    1.0000,  # dummy
    1.2000,  # H,  1
    1.4000,  # He, 2
    1.8200,  # Li, 3
    1.3725,  # Be, 4
    0.7950,  # B,  5
    1.7000,  # C,  6
    1.5500,  # N,  7
    1.5200,  # O,  8
    1.4700,  # F,  9
    1.5400,  # Ne, 10
    2.2700,  # Na, 11
    1.7300,  # Mg, 12
    1.7000,  # Al, 13
    2.1000,  # Si, 14
    1.8000,  # P,  15
    1.8000,  # S,  16
    1.7500,  # Cl, 17
    1.8800,  # Ar, 18
    2.7500,  # K,  19
    2.4500,  # Ca, 20
    1.3700,  # Sc, 21
    1.3700,  # Ti, 22
    1.3700,  # V,  23
    1.3700,  # Cr, 24
    1.3700,  # Mn, 25
    1.4560,  # Fe, 26
    0.8800,  # Co, 27
    0.6900,  # Ni, 28
    0.7200,  # Cu, 29
    0.7400,  # Zn, 30
    1.3700,  # Ga, 31
    1.9500,  # Ge, 32
    1.8500,  # As, 33
    1.9000,  # Se, 34
    1.8500,  # Br, 35
    2.0200,  # Kr, 36
    1.5800,  # Rb, 37
    2.1510,  # Sr, 38
    1.8010,  # Y,  39
    1.6020,  # Zr, 40
    1.4680,  # Nb, 41
    1.5260,  # Mo, 42
    1.3600,  # Tc, 43
    1.3390,  # Ru, 44
    1.3450,  # Rh, 45
    1.3760,  # Pd, 46
    1.2700,  # Ag, 47
    1.4240,  # Cd, 48
    1.6630,  # In, 49
    2.1000,  # Sn, 50
    2.0500,  # Sb, 51
    2.0600,  # Te, 52
    1.9800,  # I,  53
    2.0000,  # Xe, 54
    1.8400,  # Cs, 55
    2.2430,  # Ba, 56
    1.8770,  # La, 57
    None,   # Ce, 58
    None,   # Pr, 59
    None,   # Nd, 60
    None,   # Pm, 61
    None,   # Sm, 62
    None,   # Eu, 63
    None,   # Gd, 64
    None,   # Tb, 65
    None,   # Dy, 66
    None,   # Ho, 67
    None,   # Er, 68
    None,   # Tm, 69
    None,   # Yb, 70
    2.1700,  # Lu, 71
    1.5800,  # Hf, 72
    1.4670,  # Ta, 73
    1.5340,  # W,  74
    1.3750,  # Re, 75
    1.3530,  # Os, 76
    1.3570,  # Ir, 77
    1.7500,  # Pt, 78
    1.6600,  # Au, 79
    1.5500,  # Hg, 80
    1.9600,  # Tl, 81
    2.0200,  # Pb, 82
    2.1500]  # Bi, 83
# Po, 84
# At, 85
# Rn, 86
# Fr, 87
# Ra, 88
# Ac, 89
# Th, 90
# Pa, 91
# U,  92
# Np, 93
# Pu, 94
# Am, 95
# Cm, 96
# Bk, 97
# Cf, 98
# Es, 99
# Fm, 100
# Md, 101
# No, 102
# Lr, 103
# Rf, 104
# Db, 105
# Sg, 106
# Bh, 107
# Hs, 108
# Mt, 109
# Ds, 110
# Rg, 111

radius = rvdw

spectrumcolors = [
    '#800000',
    '#822828',
    '#8D533B',
    '#996675',
    '#9966A9',
    '#800080',
    '#65009B',
    '#4800E1',
    '#0400D0',
    '#0044DC',
    '#0172E2',
    '#019FE8',
    '#0BAFA2',
    '#17B34D',
    '#00D41C',
    '#00FF00',
    '#80FF00',
    '#C8FF00',
    '#FFFF00',
    '#FFDB00',
    '#FFB600',
    '#FF9200',
    '#FF6D00',
    '#FF4900',
    '#FF0000',
    '#FF0080',
    '#FF69B4',
    '#FF00FF',
    '#A800B9']
