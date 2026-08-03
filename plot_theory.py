#!/usr/bin/env python3
"""
===============================================================================
DaMaSCUS-NewEarth-EVAP Theoretical Value Calculation and Plotting
===============================================================================

Variable Definitions and Physical Meaning
------------------
1.  Velocity-related
    v_rel       : Relative velocity between dark matter particle and target [km/s]
                  Key parameter determining scattering energy in direct detection
    
    u_x         : x-component of dark matter particle velocity [km/s]
                  Velocity vector v = (u_x, u_y, u_z)
                  u_x = |v| * cos(θ), where θ is the angle between v and x-axis
    
    v_sun       : Sun's velocity vector relative to the dark matter halo [km/s]
                  Galactic coordinates: v_sun = (11.1, 232.2, 7.3) km/s
                  Includes Galactic rotation (220 km/s) + Solar peculiar motion
    
    v_earth     : Earth's total velocity vector relative to dark matter halo [km/s]
                  v_earth = v_sun + v_orbit(t)
                  where v_orbit(t) is Earth's orbital velocity (~29.8 km/s, seasonally varying)
    
    v_earth_mag : Earth's velocity magnitude [km/s]
                  Typical value ~244 km/s, varying between 203~262 km/s seasonally

2.  Dark Matter Parameters
    DM_mass     : Dark matter particle mass [GeV]
    DM_sigma    : Dark matter-nucleon scattering cross section [cm²]
    DM_density  : Local dark matter density [GeV/cm³]

3.  Target Particle Parameters
    target_mass : Target nucleus mass [GeV]
    mu          : Dark matter-target reduced mass [GeV]
                  mu = (DM_mass * target_mass) / (DM_mass + target_mass)

4.  Velocity Distribution (SHM)
    v0          : SHM most probable velocity [km/s] (typical 220 km/s)
    v_esc       : SHM escape velocity [km/s] (typical 544 km/s)
    f(v)        : SHM velocity distribution function [s/km]

5.  Scattering Rate
    rate        : Dark matter scattering rate [cm²/s]
                  Simplified model: rate ∝ sigma * v_rel / mu²

6.  Time Parameter
    nJ2000      : Days since J2000 (January 1, 2000)
                  Used to calculate Earth's orbital velocity time variation (annual modulation)
                  nJ2000 = 0 corresponds to the J2000 epoch

===============================================================================
"""
import numpy as np
import time
import matplotlib.pyplot as plt
from scipy.special import erf


# ============================================================
# 1. element data（atomic number、mass number）
# ============================================================
ELEMENT_DB = {
    "H":  {"Z": 1,   "A": 1.008,  "name": "Hydrogen"},
    "He": {"Z": 2,   "A": 4.003,  "name": "Helium"},
    "Li": {"Z": 3,   "A": 6.941,  "name": "Lithium"},
    "Be": {"Z": 4,   "A": 9.012,  "name": "Beryllium"},
    "B":  {"Z": 5,   "A": 10.81,  "name": "Boron"},
    "C":  {"Z": 6,   "A": 12.011, "name": "Carbon"},
    "N":  {"Z": 7,   "A": 14.007, "name": "Nitrogen"},
    "O":  {"Z": 8,   "A": 15.999, "name": "Oxygen"},
    "F":  {"Z": 9,   "A": 18.998, "name": "Fluorine"},
    "Ne": {"Z": 10,  "A": 20.180, "name": "Neon"},
    "Na": {"Z": 11,  "A": 22.990, "name": "Sodium"},
    "Mg": {"Z": 12,  "A": 24.305, "name": "Magnesium"},
    "Al": {"Z": 13,  "A": 26.982, "name": "Aluminium"},
    "Si": {"Z": 14,  "A": 28.086, "name": "Silicon"},
    "P":  {"Z": 15,  "A": 30.974, "name": "Phosphorus"},
    "S":  {"Z": 16,  "A": 32.065, "name": "Sulfur"},
    "Cl": {"Z": 17,  "A": 35.453, "name": "Chlorine"},
    "Ar": {"Z": 18,  "A": 39.948, "name": "Argon"},
    "K":  {"Z": 19,  "A": 39.098, "name": "Potassium"},
    "Ca": {"Z": 20,  "A": 40.078, "name": "Calcium"},
    "Sc": {"Z": 21,  "A": 44.956, "name": "Scandium"},
    "Ti": {"Z": 22,  "A": 47.867, "name": "Titanium"},
    "V":  {"Z": 23,  "A": 50.942, "name": "Vanadium"},
    "Cr": {"Z": 24,  "A": 51.996, "name": "Chromium"},
    "Mn": {"Z": 25,  "A": 54.938, "name": "Manganese"},
    "Fe": {"Z": 26,  "A": 55.845, "name": "Iron"},
    "Co": {"Z": 27,  "A": 58.933, "name": "Cobalt"},
    "Ni": {"Z": 28,  "A": 58.693, "name": "Nickel"},
    "Cu": {"Z": 29,  "A": 63.546, "name": "Copper"},
    "Zn": {"Z": 30,  "A": 65.380, "name": "Zinc"},
    "Ga": {"Z": 31,  "A": 69.723, "name": "Gallium"},
    "Ge": {"Z": 32,  "A": 72.630, "name": "Germanium"},
    "As": {"Z": 33,  "A": 74.922, "name": "Arsenic"},
    "Se": {"Z": 34,  "A": 78.971, "name": "Selenium"},
    "Br": {"Z": 35,  "A": 79.904, "name": "Bromine"},
    "Kr": {"Z": 36,  "A": 83.798, "name": "Krypton"},
    "Rb": {"Z": 37,  "A": 85.468, "name": "Rubidium"},
    "Sr": {"Z": 38,  "A": 87.620, "name": "Strontium"},
    "Y":  {"Z": 39,  "A": 88.906, "name": "Yttrium"},
    "Zr": {"Z": 40,  "A": 91.224, "name": "Zirconium"},
    "Nb": {"Z": 41,  "A": 92.906, "name": "Niobium"},
    "Mo": {"Z": 42,  "A": 95.950, "name": "Molybdenum"},
    "Ru": {"Z": 44,  "A": 101.07, "name": "Ruthenium"},
    "Rh": {"Z": 45,  "A": 102.91, "name": "Rhodium"},
    "Pd": {"Z": 46,  "A": 106.42, "name": "Palladium"},
    "Ag": {"Z": 47,  "A": 107.87, "name": "Silver"},
    "Cd": {"Z": 48,  "A": 112.41, "name": "Cadmium"},
    "In": {"Z": 49,  "A": 114.82, "name": "Indium"},
    "Sn": {"Z": 50,  "A": 118.71, "name": "Tin"},
    "Sb": {"Z": 51,  "A": 121.76, "name": "Antimony"},
    "Te": {"Z": 52,  "A": 127.60, "name": "Tellurium"},
    "I":  {"Z": 53,  "A": 126.90, "name": "Iodine"},
    "Xe": {"Z": 54,  "A": 131.29, "name": "Xenon"},
    "Cs": {"Z": 55,  "A": 132.91, "name": "Cesium"},
    "Ba": {"Z": 56,  "A": 137.33, "name": "Barium"},
    "La": {"Z": 57,  "A": 138.91, "name": "Lanthanum"},
    "Ce": {"Z": 58,  "A": 140.12, "name": "Cerium"},
    "Pr": {"Z": 59,  "A": 140.91, "name": "Praseodymium"},
    "Nd": {"Z": 60,  "A": 144.24, "name": "Neodymium"},
    "Sm": {"Z": 62,  "A": 150.36, "name": "Samarium"},
    "Eu": {"Z": 63,  "A": 151.96, "name": "Europium"},
    "Gd": {"Z": 64,  "A": 157.25, "name": "Gadolinium"},
    "Tb": {"Z": 65,  "A": 158.93, "name": "Terbium"},
    "Dy": {"Z": 66,  "A": 162.50, "name": "Dysprosium"},
    "Ho": {"Z": 67,  "A": 164.93, "name": "Holmium"},
    "Er": {"Z": 68,  "A": 167.26, "name": "Erbium"},
    "Tm": {"Z": 69,  "A": 168.93, "name": "Thulium"},
    "Yb": {"Z": 70,  "A": 173.05, "name": "Ytterbium"},
    "Lu": {"Z": 71,  "A": 174.97, "name": "Lutetium"},
    "Hf": {"Z": 72,  "A": 178.49, "name": "Hafnium"},
    "Ta": {"Z": 73,  "A": 180.95, "name": "Tantalum"},
    "W":  {"Z": 74,  "A": 183.84, "name": "Tungsten"},
    "Re": {"Z": 75,  "A": 186.21, "name": "Rhenium"},
    "Os": {"Z": 76,  "A": 190.23, "name": "Osmium"},
    "Ir": {"Z": 77,  "A": 192.22, "name": "Iridium"},
    "Pt": {"Z": 78,  "A": 195.08, "name": "Platinum"},
    "Au": {"Z": 79,  "A": 196.97, "name": "Gold"},
    "Hg": {"Z": 80,  "A": 200.59, "name": "Mercury"},
    "Tl": {"Z": 81,  "A": 204.38, "name": "Thallium"},
    "Pb": {"Z": 82,  "A": 207.20, "name": "Lead"},
    "Bi": {"Z": 83,  "A": 208.98, "name": "Bismuth"},
    "Th": {"Z": 90,  "A": 232.04, "name": "Thorium"},
    "U":  {"Z": 92,  "A": 238.03, "name": "Uranium"},
}


# ============================================================
# SD 自旋資料庫（原子核自旋結構）
# ============================================================
# format：{
#     "element": {
#         "J": Spin quantum number,
#         "S_p": Proton spin expectation,
#         "S_n": Neutron spin expectation,
#         "isotope_abundance": isotope_abundance,
#         "isotope": isotope mass number
#     }
# }
#
# 注意：
# 1. only list the isotope have spin （odd A or irrational Z/N combination）
# 2. odd number A nuclides（like ¹⁶O, ²⁸Si, ⁵⁶Fe）have spin of 0，no SD 
# 3. data sourced from: https://arxiv.org/abs/1003.1912, https://arxiv.org/abs/1203.3542

SD_SPIN_DB = {
    # ============================================
    # 第 1 週期
    # ============================================
    "H": {
        "J": 0.5,
        "S_p": 0.5,
        "S_n": 0.0,
        "isotope_abundance": 1.0,
        "isotope": 1,
        "name": "Hydrogen-1"
    },
    "H2": {  # 氘
        "J": 1.0,
        "S_p": 0.0,
        "S_n": 0.0,  # 氘自旋為 1，但自旋結構因子為 0（簡化）
        "isotope_abundance": 0.00015,
        "isotope": 2,
        "name": "Deuterium"
    },
    
    # ============================================
    # 第 2 週期
    # ============================================
    "He": {
        "J": 0.0,  # ⁴He 自旋為 0，無 SD
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 1.0,
        "isotope": 4,
        "name": "Helium-4"
    },
    
    # ============================================
    # 第 3 週期
    # ============================================
    "Li": {
        "J": 1.0,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 0.925,  # ⁷Li
        "isotope": 7,
        "name": "Lithium-7"
    },
    
    # ============================================
    # 第 4 週期
    # ============================================
    "C": {
        "J": 0.0,  # ¹²C 自旋為 0
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 1.0,
        "isotope": 12,
        "name": "Carbon-12"
    },
    "N": {
        "J": 1.0,
        "S_p": 0.0,
        "S_n": 0.5,
        "isotope_abundance": 0.996,  # ¹⁴N
        "isotope": 14,
        "name": "Nitrogen-14"
    },
    "O": {
        "J": 0.0,  # ¹⁶O 自旋為 0
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 1.0,
        "isotope": 16,
        "name": "Oxygen-16"
    },
    "O17": {  # ¹⁷O（稀有同位素）
        "J": 2.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 0.00038,
        "isotope": 17,
        "name": "Oxygen-17"
    },
    
    # ============================================
    # 第 3 週期（元素 Na, Mg, Al, Si, P, S, Cl, Ar）
    # ============================================
    "Na": {
        "J": 1.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 1.0,  # ²³Na
        "isotope": 23,
        "name": "Sodium-23"
    },
    "Mg": {
        "J": 0.0,  # ²⁴Mg 自旋為 0
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.79,
        "isotope": 24,
        "name": "Magnesium-24"
    },
    "Al": {
        "J": 2.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 1.0,  # ²⁷Al
        "isotope": 27,
        "name": "Aluminium-27"
    },
    "Si": {
        "J": 0.0,  # ²⁸Si 自旋為 0
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.922,
        "isotope": 28,
        "name": "Silicon-28"
    },
    "Si29": {  # ²⁹Si（稀有同位素）
        "J": 0.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 0.047,
        "isotope": 29,
        "name": "Silicon-29"
    },
    "P": {
        "J": 0.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 1.0,  # ³¹P
        "isotope": 31,
        "name": "Phosphorus-31"
    },
    "S": {
        "J": 0.0,  # ³²S 自旋為 0
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.95,
        "isotope": 32,
        "name": "Sulfur-32"
    },
    "Cl": {
        "J": 1.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 0.758,  # ³⁵Cl
        "isotope": 35,
        "name": "Chlorine-35"
    },
    "Ar": {
        "J": 0.0,  # ⁴⁰Ar 自旋為 0
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.996,
        "isotope": 40,
        "name": "Argon-40"
    },
    
    # ============================================
    # 第 4 週期（K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr）
    # ============================================
    "K": {
        "J": 1.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 0.932,  # ³⁹K
        "isotope": 39,
        "name": "Potassium-39"
    },
    "Ca": {
        "J": 0.0,  # ⁴⁰Ca 自旋為 0
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.969,
        "isotope": 40,
        "name": "Calcium-40"
    },
    "Sc": {
        "J": 3.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 1.0,  # ⁴⁵Sc
        "isotope": 45,
        "name": "Scandium-45"
    },
    "Ti": {
        "J": 0.0,  # ⁴⁸Ti 自旋為 0
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.737,
        "isotope": 48,
        "name": "Titanium-48"
    },
    "V": {
        "J": 3.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 0.997,  # ⁵¹V
        "isotope": 51,
        "name": "Vanadium-51"
    },
    "Cr": {
        "J": 0.0,  # ⁵²Cr 自旋為 0
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.838,
        "isotope": 52,
        "name": "Chromium-52"
    },
    "Mn": {
        "J": 2.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 1.0,  # ⁵⁵Mn
        "isotope": 55,
        "name": "Manganese-55"
    },
    "Fe": {
        "J": 0.0,  # ⁵⁶Fe 自旋為 0
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.917,
        "isotope": 56,
        "name": "Iron-56"
    },
    "Fe57": {  # ⁵⁷Fe（稀有同位素）
        "J": 0.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 0.021,
        "isotope": 57,
        "name": "Iron-57"
    },
    "Co": {
        "J": 3.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 1.0,  # ⁵⁹Co
        "isotope": 59,
        "name": "Cobalt-59"
    },
    "Ni": {
        "J": 0.0,  # ⁵⁸Ni 自旋為 0
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.682,
        "isotope": 58,
        "name": "Nickel-58"
    },
    "Cu": {
        "J": 1.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 0.691,  # ⁶³Cu
        "isotope": 63,
        "name": "Copper-63"
    },
    "Zn": {
        "J": 0.0,  # ⁶⁴Zn 自旋為 0
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.492,
        "isotope": 64,
        "name": "Zinc-64"
    },
    "Ga": {
        "J": 1.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 0.601,  # ⁶⁹Ga
        "isotope": 69,
        "name": "Gallium-69"
    },
    "Ge": {
        "J": 0.0,  # ⁷⁰Ge 自旋為 0
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.205,
        "isotope": 70,
        "name": "Germanium-70"
    },
    "As": {
        "J": 1.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 1.0,  # ⁷⁵As
        "isotope": 75,
        "name": "Arsenic-75"
    },
    "Se": {
        "J": 0.0,  # ⁷⁸Se 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.237,
        "isotope": 78,
        "name": "Selenium-78"
    },
    "Br": {
        "J": 1.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 0.507,  # ⁷⁹Br
        "isotope": 79,
        "name": "Bromine-79"
    },
    "Kr": {
        "J": 0.0,  # ⁸⁰Kr 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.023,
        "isotope": 80,
        "name": "Krypton-80"
    },
    
    # ============================================
    # 第 5 週期（Rb, Sr, Y, Zr, Nb, Mo, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba）
    # ============================================
    "Rb": {
        "J": 1.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 0.722,  # ⁸⁵Rb
        "isotope": 85,
        "name": "Rubidium-85"
    },
    "Sr": {
        "J": 0.0,  # ⁸⁸Sr 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.826,
        "isotope": 88,
        "name": "Strontium-88"
    },
    "Y": {
        "J": 0.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 1.0,  # ⁸⁹Y
        "isotope": 89,
        "name": "Yttrium-89"
    },
    "Zr": {
        "J": 0.0,  # ⁹⁰Zr 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.514,
        "isotope": 90,
        "name": "Zirconium-90"
    },
    "Nb": {
        "J": 4.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 1.0,  # ⁹³Nb
        "isotope": 93,
        "name": "Niobium-93"
    },
    "Mo": {
        "J": 0.0,  # ⁹²Mo 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.148,
        "isotope": 92,
        "name": "Molybdenum-92"
    },
    "Ru": {
        "J": 0.0,  # ¹⁰⁰Ru 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.126,
        "isotope": 100,
        "name": "Ruthenium-100"
    },
    "Rh": {
        "J": 0.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 1.0,  # ¹⁰³Rh
        "isotope": 103,
        "name": "Rhodium-103"
    },
    "Pd": {
        "J": 0.0,  # ¹⁰⁶Pd 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.273,
        "isotope": 106,
        "name": "Palladium-106"
    },
    "Ag": {
        "J": 0.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 0.518,  # ¹⁰⁷Ag
        "isotope": 107,
        "name": "Silver-107"
    },
    "Cd": {
        "J": 0.0,  # ¹¹⁰Cd 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.125,
        "isotope": 110,
        "name": "Cadmium-110"
    },
    "In": {
        "J": 4.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 0.957,  # ¹¹⁵In
        "isotope": 115,
        "name": "Indium-115"
    },
    "Sn": {
        "J": 0.0,  # ¹¹⁶Sn 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.146,
        "isotope": 116,
        "name": "Tin-116"
    },
    "Sb": {
        "J": 2.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 0.572,  # ¹²¹Sb
        "isotope": 121,
        "name": "Antimony-121"
    },
    "Te": {
        "J": 0.0,  # ¹²⁰Te 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.0009,
        "isotope": 120,
        "name": "Tellurium-120"
    },
    "I": {
        "J": 2.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 1.0,  # ¹²⁷I
        "isotope": 127,
        "name": "Iodine-127"
    },
    "Xe": {
        "J": 0.0,  # ¹²⁸Xe 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.001,
        "isotope": 128,
        "name": "Xenon-128"
    },
    "Cs": {
        "J": 3.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 1.0,  # ¹³³Cs
        "isotope": 133,
        "name": "Cesium-133"
    },
    "Ba": {
        "J": 0.0,  # ¹³⁴Ba 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.024,
        "isotope": 134,
        "name": "Barium-134"
    },
    
    # ============================================
    # 第 6 週期（La, Ce, Pr, Nd, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi）
    # ============================================
    "La": {
        "J": 3.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 0.999,  # ¹³⁹La
        "isotope": 139,
        "name": "Lanthanum-139"
    },
    "Ce": {
        "J": 0.0,  # ¹⁴⁰Ce 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.885,
        "isotope": 140,
        "name": "Cerium-140"
    },
    "Pr": {
        "J": 2.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 1.0,  # ¹⁴¹Pr
        "isotope": 141,
        "name": "Praseodymium-141"
    },
    "Nd": {
        "J": 0.0,  # ¹⁴²Nd 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.272,
        "isotope": 142,
        "name": "Neodymium-142"
    },
    "Sm": {
        "J": 0.0,  # ¹⁴⁸Sm 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.114,
        "isotope": 148,
        "name": "Samarium-148"
    },
    "Eu": {
        "J": 2.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 0.478,  # ¹⁵¹Eu
        "isotope": 151,
        "name": "Europium-151"
    },
    "Gd": {
        "J": 0.0,  # ¹⁵²Gd 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "abundance": 0.002,
        "isotope": 152,
        "name": "Gadolinium-152"
    },
    "Tb": {
        "J": 1.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "abundance": 1.0,  # ¹⁵⁹Tb
        "isotope": 159,
        "name": "Terbium-159"
    },
    "Dy": {
        "J": 0.0,  # ¹⁵⁸Dy 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.0001,
        "isotope": 158,
        "name": "Dysprosium-158"
    },
    "Ho": {
        "J": 3.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 1.0,  # ¹⁶⁵Ho
        "isotope": 165,
        "name": "Holmium-165"
    },
    "Er": {
        "J": 0.0,  # ¹⁶⁶Er 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.336,
        "isotope": 166,
        "name": "Erbium-166"
    },
    "Tm": {
        "J": 0.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 1.0,  # ¹⁶⁹Tm
        "isotope": 169,
        "name": "Thulium-169"
    },
    "Yb": {
        "J": 0.0,  # ¹⁷⁰Yb 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.001,
        "isotope": 170,
        "name": "Ytterbium-170"
    },
    "Lu": {
        "J": 3.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 0.974,  # ¹⁷⁵Lu
        "isotope": 175,
        "name": "Lutetium-175"
    },
    "Hf": {
        "J": 0.0,  # ¹⁷⁶Hf 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.052,
        "isotope": 176,
        "name": "Hafnium-176"
    },
    "Ta": {
        "J": 3.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 0.999,  # ¹⁸¹Ta
        "isotope": 181,
        "name": "Tantalum-181"
    },
    "W": {
        "J": 0.0,  # ¹⁸²W 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.265,
        "isotope": 182,
        "name": "Tungsten-182"
    },
    "Re": {
        "J": 2.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 0.374,  # ¹⁸⁵Re
        "isotope": 185,
        "name": "Rhenium-185"
    },
    "Os": {
        "J": 0.0,  # ¹⁸⁶Os 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.015,
        "isotope": 186,
        "name": "Osmium-186"
    },
    "Ir": {
        "J": 1.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 0.373,  # ¹⁹¹Ir
        "isotope": 191,
        "name": "Iridium-191"
    },
    "Pt": {
        "J": 0.0,  # ¹⁹²Pt 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.0078,
        "isotope": 192,
        "name": "Platinum-192"
    },
    "Au": {
        "J": 1.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 1.0,  # ¹⁹⁷Au
        "isotope": 197,
        "name": "Gold-197"
    },
    "Hg": {
        "J": 0.0,  # ¹⁹⁸Hg 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.001,
        "isotope": 198,
        "name": "Mercury-198"
    },
    "Tl": {
        "J": 0.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 0.705,  # ²⁰³Tl
        "isotope": 203,
        "name": "Thallium-203"
    },
    "Pb": {
        "J": 0.0,  # ²⁰⁴Pb 自旋為 0（簡化）
        "S_p": 0.0,
        "S_n": 0.0,
        "isotope_abundance": 0.014,
        "isotope": 204,
        "name": "Lead-204"
    },
    "Bi": {
        "J": 4.5,
        "S_p": 0.5,
        "S_n": 0.5,
        "isotope_abundance": 1.0,  # ²⁰⁹Bi
        "isotope": 209,
        "name": "Bismuth-209"
    },
}

# ============================================================
# 物理常數（自然單位制）
# ============================================================
G = 6.674e-11          # gravitational constant [m^3 kg^-1 s^-2]
M_earth = 5.972e24     # Earth mass [kg]
R_earth = 6.371e6      # Earth radius [m]
M_sun = 1.989e30       # Sun mass [kg]
R_sun = 6.957e8        # Sun radius [m]

# DM parameters
DM_mass = 0.1          # GeV
DM_sigma_SI_p = 1e-35            # cm²
DM_sigma_SD_p = 1e-36            # cm²
DM_sigma_electron = 1.0e-80      # cm²
DM_density = 0.3       # GeV/cm^3    
m_p = 0.938  # GeV
m_N = 0.939  # GeV

def reduced_mass(m1, m2):
    """
    Calculate the reduced mass of two particles.
    
    Parameters:
        m1, m2 : Masses of the two particles [GeV]
    
    Returns:
        mu : Reduced mass [GeV]
    
    Formula:
        mu = (m1 * m2) / (m1 + m2)
    """
    return (m1 * m2) / (m1 + m2)

def SHM_velocity_distribution(v, v0=220.0, v_esc=544.0):
    """
    Standard Halo Model (SHM) velocity distribution
    v: velocity [km/s]
    v0: most probable velocity [km/s]
    v_esc: escape velocity [km/s]
    """
    if v > v_esc:
        return 0.0
    # Maxwell-Boltzmann Distribution (Truncation)
    # use scipy.special.erf
    N = 1.0 / (np.sqrt(np.pi) * v0) * (1 / erf(v_esc / v0) - 2 / np.sqrt(np.pi) * (v_esc / v0) * np.exp(-(v_esc / v0)**2))
    return N * np.exp(-(v/v0)**2)

def Xi_function(v_d,v_esc,v_earth):
    """
    Calculate the Xi function for dark matter capture.
    
    Xi = erf(v_esc / v_earth) - (2/sqrt(pi)) * (v_esc / v_earth) * exp(-v_esc^2 / v_earth^2)
    
    Parameters:
        v_d     : Velocity dispersion parameter [km/s] (typically sqrt(3/2) * v_earth)
        v_esc   : Escape velocity at radius r [km/s]
        v_earth : Earth's velocity in DM halo [km/s] (~244 km/s)
    
    Returns:
        Xi : Xi function value
    """
    from math import erf
    if v_earth <= 0:
        return 0.0
    xi = (v_d**2 * np.exp(-3*v_earth**2/(2*v_d**2)) + erf(v_earth/v_d) * np.sqrt(np.pi/6) *v_d/v_earth*(v_d**2 +3*v_earth**2 + 3*v_esc**2))
    return xi


# ============================================================
# Cross section types (velocity and momentum dependence)
# ============================================================
def nuclear_form_factor_for_nucleus(A, q,scattering_type="SI"):
    """
    Calculate the nuclear form factor for a nucleus with mass number A and momentum transfer q.
    
    Parameters:
        A : Mass number of the nucleus
        q : Momentum transfer [GeV]
    
    Returns:
        F(q) : Nuclear form factor
    """
    if scattering_type == "electron":
        return 1.0
    
    fm_to_GeV = 5.067731237e-3

    if scattering_type == "SI":
        r_i = 0.89 * A**(1/3) + 0.3
    elif scattering_type == "SD":
        r_i = np.sqrt(3/2)*(1.7 * A**(1/3) - 0.28 - 0.78 * (A**(-1/3) -3.8 + np.sqrt((A**(-1/3) -3.8)**2 +0.2)))
    else:
        raise ValueError("scattering_type must be 'SI', 'SD', or 'electron'")

    
    r_i *= fm_to_GeV

    x = q * r_i / np.sqrt(3)

    F_i = np.exp(-x**2)

    if x < 1e-10:
        # 避免除以 0
        return 1.0
    return F_i


def cross_section_constant(sigma_0, v_rel, q):
    """
    Constant (velocity-independent and isotropic) scattering cross section.
    
    Parameters:
        sigma_0 : Reference cross section [cm²]
        v_rel   : Relative velocity [km/s]
        q       : Momentum transfer [GeV]
    
    Returns:
        sigma : Cross section [cm²]
    """
    return sigma_0


def cross_section_v2_dependent(sigma_0, v_rel, q):
    """
    v_rel²-dependent scattering cross section.
    
    Formula:
        σ = σ_0 * (v_rel / v_ref)²
    
    Parameters:
        sigma_0 : Reference cross section [cm²]
        v_rel   : Relative velocity [km/s]
        q       : Momentum transfer [GeV]
    
    Returns:
        sigma : Cross section [cm²]
    """
    v_ref = 220.0  # Reference velocity [km/s]
    return sigma_0 * (v_rel / v_ref)**2


def cross_section_q2_dependent(sigma_0, v_rel, q):
    """
    q²-dependent scattering cross section.
    
    Formula:
        σ = σ_0 * (q / q_ref)²
    
    Parameters:
        sigma_0 : Reference cross section [cm²]
        v_rel   : Relative velocity [km/s]
        q       : Momentum transfer [GeV]
    
    Returns:
        sigma : Cross section [cm²]
    """
    q_ref = 0.1  # Reference momentum transfer [GeV]
    return sigma_0 * (q / q_ref)**2


def get_cross_section(sigma_0, v_rel, q, cross_section_type="constant"):
    """
    Get cross section based on type.
    
    Parameters:
        sigma_0            : Reference cross section [cm²]
        v_rel              : Relative velocity [km/s]
        q                  : Momentum transfer [GeV]
        cross_section_type : "constant", "v2_dependent", or "q2_dependent"
    
    Returns:
        sigma : Cross section [cm²]
    """
    if cross_section_type == "constant":
        return cross_section_constant(sigma_0, v_rel, q)
    elif cross_section_type == "v2_dependent":
        return cross_section_v2_dependent(sigma_0, v_rel, q)
    elif cross_section_type == "q2_dependent":
        return cross_section_q2_dependent(sigma_0, v_rel, q)
    else:
        raise ValueError("cross_section_type must be 'constant', 'v2_dependent', or 'q2_dependent'")

# ============================================================
# Theoretical function: Velocity distribution (SHM)
# ============================================================

# ============================================================
# Theoretical function: u_x (velocity component)
# ============================================================
def u_x_theory(v_rel, theta=np.pi/4):
    """
    u_x theoretical value: x-component of velocity
    v_rel: Relative velocity magnitude
    theta: Angle with x-axis
    """
    return v_rel * np.cos(theta)

# ============================================================
# main program：Calculate Earth velocity in galactic coordinates and DM velocity distribution
# ============================================================
def earth_velocity(nJ2000=0.0):
    """
    Calculate the velocity of the Earth in galactic coordinates (km/s)
    Equivalent to obscura::Earth_Velocity(nJ2000)

    Parameters:
        nJ2000: Number of days since J2000 (default 0 = J2000时刻)

    Returns:
        np.array: Earth velocity vector [vx, vy, vz] (km/s)
    """
    deg = np.pi / 180.0
    
    # 1. Sun velocity (galactic rotation + peculiar motion)
    v_sun = np.array([0.0, 220.0, 0.0]) + np.array([11.1, 12.2, 7.3])
    
    # 2. Earth orbital velocity (simplified J2000时刻)
    e = 0.01671  # orbital eccentricity
    L = (280.46 * deg + nJ2000 * 0.9856474 * deg) % (2 * np.pi)
    omega = (282.932 * deg + nJ2000 * 0.0000471 * deg) % (2 * np.pi)
    
    # Ecliptic coordinate system basis vectors (approximate)    
    exEcl = np.array([0.054876, -0.494109, 0.867666])
    eyEcl = np.array([0.993824, 0.110992, 0.000352])
    
    v_orbital = 29.79  # km/s
    uE = -v_orbital * (np.sin(L) + e * np.sin(2*L - omega)) * exEcl + \
         v_orbital * (np.cos(L) + e * np.cos(2*L - omega)) * eyEcl
    
    return v_sun + uE



def v_earth_DM(nJ2000=0.0):
    """
    Calculate the total velocity of the Earth relative to the dark matter halo (km/s)
    """
    vel = earth_velocity(nJ2000)
    return np.linalg.norm(vel)

v_earth = v_earth_DM(0.0)          # 約 244 km/s
v_d = np.sqrt(3/2) * v_earth       # 約 298.8 km/s

def Maxwell_Boltzmann_distribution(v_d, u_x, v_earth_DM):
    """
    Maxwell-Boltzmann velocity distribution (1D)
    
    Parameters:
        v_d        : velocity spread parameter [km/s] = sqrt(3/2) * v_earth_DM
        u_x        : velocity x component [km/s]（可為陣列）
        v_earth_DM : Earth velocity magnitude [km/s]
    
    Returns:
        MBD : velocity distribution value
    """
    # 避免除以零
    if v_d <= 0 or v_earth_DM <= 0:
        return np.zeros_like(u_x)
    
    factor = np.sqrt(3 / (2 * np.pi)) * u_x / (v_d * v_earth_DM)
    exp1 = np.exp(-3 * (u_x - v_earth_DM)**2 / (2 * v_d**2))
    exp2 = np.exp(-3 * (u_x + v_earth_DM)**2 / (2 * v_d**2))
    
    return factor * (exp1 - exp2)

def load_earth_composition(filepath="data/earth_prem.dat"):
    """
    read earth_prem.dat and return radius, density, temperature, and abundances
    
    Returns:
        dict: {
            "radius": [r_km, ...],
            "density": [rho, ...],
            "temperature": [T, ...],
            "abundances": {element: [ppm, ...], ...}
        }
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # find the start of the data table
    data_start = 0
    for i, line in enumerate(lines):
        if line.startswith("# The Table begins here."):
            data_start = i + 1
            break
    
    # analyze the header line to get element names
    header_line = lines[data_start - 2]
    parts = header_line.split()
    
    elements = []
    for part in parts:
        if "_ppm" in part:
            elements.append(part.replace("_ppm", ""))
    
    # read data
    radii = []
    densities = []
    temperatures = []
    abundances = {elem: [] for elem in elements}
    
    for line in lines[data_start:]:
        if line.strip() == "" or line.startswith("#"):
            continue
        cols = line.strip().split()
        if len(cols) < 4 + len(elements):
            continue
        
        radii.append(float(cols[0]))
        densities.append(float(cols[2]))
        temperatures.append(float(cols[3]))
        
        for i, elem in enumerate(elements):
            abundances[elem].append(float(cols[4 + i]))
    
    return {
        "radius": np.array(radii),
        "density": np.array(densities),
        "temperature": np.array(temperatures),
        "abundances": {k: np.array(v) for k, v in abundances.items()}
    }

def number_density_at_radius(earth_data, radius_idx, elem):
    """
    Calculate the number density of element at a given radius.
    
    Parameters:
        earth_data : Earth composition data from load_earth_composition()
        radius_idx : Radius index
        elem       : Element symbol (e.g., "Fe", "O", "Si")
    
    Returns:
        n_i : Number density [cm⁻³]
    """
    if elem not in ELEMENT_DB:
        return 0.0
    
    A = ELEMENT_DB[elem]["A"]
    
    # Mass fraction
    mass_fraction = earth_data["abundances"][elem][radius_idx] / 1e6
    
    # Density [g/cm³]
    density_g_cm3 = earth_data["density"][radius_idx]
    
    # Convert density to [GeV/cm³] (1 g/cm³ ≈ 0.562 GeV/cm³)
    density_GeV = density_g_cm3 * 0.562
    
    # Number density: n = ρ * f_mass / (A * m_N)
    m_N = 0.9315  # Nucleon mass [GeV]
    number_density = density_GeV * mass_fraction / (A * m_N)  # [cm⁻³]
    
    return number_density

def velocity_to_recoil_energy(w, v, m_A):
    """
    Calculate the recoil energy from velocity change.
    
    Parameters:
        w   : Pre-scattering velocity [km/s]
        v   : Post-scattering velocity [km/s]
        m_A : Target nucleus mass [GeV]
    
    Returns:
        E_R : Recoil energy [GeV]
    """
    km_s_to_GeV = 1.0 / 3.0e5
    return 0.5 * m_A * (w**2 - v**2) * km_s_to_GeV**2


def recoil_energy_to_velocity(E_R, m_A):
    """
    Calculate the post-scattering velocity from recoil energy.
    
    Parameters:
        E_R : Recoil energy [GeV]
        m_A : Target nucleus mass [GeV]
    
    Returns:
        v   : Post-scattering velocity [km/s]
    """
    km_s_to_GeV = 1.0 / 3.0e5
    v_GeV = np.sqrt(2 * E_R / m_A)
    return v_GeV / km_s_to_GeV


# ============================================================
# calculate SI/SD/inelastic scattering rates at a given radius
# ============================================================

# ============================================================
# 8. Differential scattering rates (velocity-dependent)
# ============================================================

def differential_cross_section_SI(E_R, v_rel, DM_mass, A, sigma_SI_p, q=0.01, cross_section_type="constant"):
    """
    Calculate the SI differential cross section dσ/dE_R.
    
    Parameters:
        E_R       : Recoil energy [GeV] (energy transferred to the nucleus)
        v_rel     : Relative velocity [km/s]
        DM_mass   : Dark matter mass [GeV]
        A         : Mass number of the target nucleus
        sigma_SI_p: DM-proton SI cross section [cm²]
        q         : Momentum transfer [GeV] (default 0.01)
    
    Returns:
        dσ/dE_R   : Differential cross section [cm²/GeV]
    """
    # Nucleon mass [GeV]
    m_N = 0.9315
    
    # Target nucleus mass [GeV]
    m_A = A * m_N
    
    # Reduced mass [GeV]
    mu_A = reduced_mass(DM_mass, m_A)
    mu_p = reduced_mass(DM_mass, m_p)
    
    # Maximum recoil energy
    # E_R_max = 2 * mu_A^2 * v_rel^2 / m_A
    v_rel_GeV = v_rel / 3.0e5  # Convert km/s to GeV (approx)
    E_R_max = 2 * mu_A**2 * v_rel_GeV**2 / m_A
    
    if E_R > E_R_max or E_R < 0:
        return 0.0
    
    # Form factor
    F_q = nuclear_form_factor_for_nucleus(A, q, scattering_type="SI")
    
    # SI nuclear cross section (with A^2 coherence enhancement)
    sigma_p = get_cross_section(sigma_SI_p, v_rel, q, cross_section_type)


    sigma_nucleus = (mu_A / mu_p)**2 * (A**2) * sigma_SI_p * F_q**2
    
    # Differential cross section
    # dσ/dE_R = m_A * σ_nucleus / (2 * μ_A^2 * v_rel^2)
    dsigma_dE = m_A * sigma_nucleus / (2 * mu_A**2 * v_rel_GeV**2)
    
    return dsigma_dE


def differential_cross_section_SD(E_R, v_rel, DM_mass, A, J, S_p, S_n, sigma_SD_p, q=0.01):
    """
    Calculate the SD differential cross section dσ/dE_R.
    
    Parameters:
        E_R       : Recoil energy [GeV]
        v_rel     : Relative velocity [km/s]
        DM_mass   : Dark matter mass [GeV]
        A         : Mass number of the target nucleus
        J         : Nuclear spin
        S_p       : Proton spin expectation value
        S_n       : Neutron spin expectation value
        sigma_SD_p: DM-proton SD cross section [cm²]
        q         : Momentum transfer [GeV] (default 0.01)
    
    Returns:
        dσ/dE_R   : Differential cross section [cm²/GeV]
    """
    # Nucleon mass [GeV]
    m_N = 0.9315
    
    # Target nucleus mass [GeV]
    m_A = A * m_N
    
    # Reduced mass [GeV]
    mu_A = reduced_mass(DM_mass, m_A)
    mu_p = reduced_mass(DM_mass, m_p)
    
    # Maximum recoil energy
    v_rel_GeV = v_rel / 3.0e5
    E_R_max = 2 * mu_A**2 * v_rel_GeV**2 / m_A
    
    if E_R > E_R_max or E_R < 0:
        return 0.0
    
    # Form factor
    F_q = nuclear_form_factor_for_nucleus(A, q, scattering_type="SD")
    
    # Spin structure factor
    spin_factor = 4 * (J + 1) / (3 * J)
    spin_contrib = S_p + S_n
    
    # SD nuclear cross section
    sigma_nucleus = (mu_A / mu_p)**2 * spin_factor * spin_contrib**2 * sigma_SD_p * F_q**2
    
    # Differential cross section
    dsigma_dE = m_A * sigma_nucleus / (2 * mu_A**2 * v_rel_GeV**2)
    
    return dsigma_dE



def scattering_rate_SI_at_radius(earth_data, radius_idx, v_rel, DM_mass, sigma_SI_p,q):
    """
    Calculate the SI scattering rate at a specific radius.
    
    parameter:
        sigma_SI_p : SI cross section of dark matter and protons [cm²]
    """
    scattering_rate_total = 0.0
    # reduced mass of DM and proton
    mu_p = reduced_mass(DM_mass, m_p)
    
    for elem, abundance in earth_data["abundances"].items():
        if elem not in ELEMENT_DB:
            continue
        
        A = ELEMENT_DB[elem]["A"]
        m_A = A * m_N
        mu_A = reduced_mass(DM_mass, m_A)
        
        # Calculate the shape factor
        F_q = nuclear_form_factor_for_nucleus(A, q)
        
        mass_fraction = abundance[radius_idx] / 1e6
        sigma_nucleus = (mu_A / mu_p)**2 * (A**2) * sigma_SI_p * F_q**2
        
        scattering_rate_total += mass_fraction * sigma_nucleus * v_rel / (mu_A**2 + 1e-10)
    
    return scattering_rate_total

def scattering_rate_SD_at_radius(earth_data, radius_idx, v_rel, DM_mass, sigma_SD_p, sigma_SD_n=0.0):
    """
    Calculate the SD scattering rate at a specific radius
    
    Parameters:
        sigma_SD_p : SI cross section of dark matter and protons [cm²]
        sigma_SD_n : SI cross section of dark matter and neutrons [cm²] (default 0)
    """
    scattering_rate_total = 0.0
    mu_p = reduced_mass(DM_mass, m_p)
    
    for elem, abundance in earth_data["abundances"].items():
        if elem not in ELEMENT_DB:
            continue
        
        # Check for spin data.
        if elem in SD_SPIN_DB:
            spin_data = SD_SPIN_DB[elem]
            J = spin_data["J"]
            
            # J = 0 mean no SD scattering contribution
            if J == 0.0:
                continue
            
            # Spin structure factor (simplified version)
            S_p = spin_data["S_p"]
            S_n = spin_data["S_n"]
            
            # Spin structure factor formula (simplified)
            # F_SD = (J+1)/J * (S_p * sigma_SD_p + S_n * sigma_SD_n)^2
            factor = 4*(J + 1) / (3 * J)
            spin_contrib = S_p  + S_n 
            
            A = ELEMENT_DB[elem]["A"]
            m_A = A * m_N
            mu_A = reduced_mass(DM_mass, m_A)
            
            mass_fraction = abundance[radius_idx] / 1e6
            
            #SD nuclear cross section (without A^2 coherence enhancement)
            sigma_nucleus = (mu_A / mu_p) ** 2 * factor * spin_contrib**2 * DM_sigma_SD_p
            
            scattering_rate_total += mass_fraction * sigma_nucleus * v_rel / (mu_A**2 + 1e-10)
        else:
            # No spin data available, skip this step.
            continue
    
    return scattering_rate_total



def differential_scattering_rate_SI(w, v, earth_data, radius_idx, DM_mass, sigma_SI_p, q=0.01, cross_section_type="constant"):
    """
    Calculate the differential scattering rate R_-(w → v) for SI.
    
    R_-(w → v) = n_i * (dσ/dE_R) * v_rel * m_A * v
    
    Parameters:
        w           : Pre-scattering velocity [km/s]
        v           : Post-scattering velocity [km/s]
        earth_data  : Earth composition data
        radius_idx  : Radius index
        DM_mass     : Dark matter mass [GeV]
        sigma_SI_p  : DM-proton SI cross section [cm²]
        q           : Momentum transfer [GeV] (default 0.01)
    
    Returns:
        rate : Differential scattering rate [s⁻¹]
    """
    # Relative velocity (assuming target nucleus is at rest)
    v_rel = w
    
    # km/s to GeV conversion
    km_s_to_GeV = 1.0 / 3.0e5
    
    differential_rate_total  = 0.0
    m_N = 0.9315
    
    for elem, abundance in earth_data["abundances"].items():
        if elem not in ELEMENT_DB:
            continue
        
        A = ELEMENT_DB[elem]["A"]
        m_A = A * m_N  # [GeV]
        
        # Calculate recoil energy
        E_R = velocity_to_recoil_energy(w, v, m_A)  # [GeV]
        
        if E_R <= 0:
            continue
        
        # Differential cross section
        dsigma_dE = differential_cross_section_SI(E_R, v_rel, DM_mass, A, sigma_SI_p, q, cross_section_type)
        
        if dsigma_dE == 0.0:
            continue
        
        # Number density [cm⁻³]
        number_density = number_density_at_radius(earth_data, radius_idx, elem)
        
        # Jacobian factor |dE_R/dv| = m_A * v
        jacobian = m_A * v * km_s_to_GeV
        
        # Differential scattering rate
        differential_rate  = number_density * dsigma_dE * v_rel * jacobian
        
        differential_rate_total += differential_rate 
    
    return differential_rate_total


def differential_scattering_rate_SD(w, v, earth_data, radius_idx, DM_mass, sigma_SD_p, sigma_SD_n=0.0, q=0.01):
    """
    Calculate the differential scattering rate R_-(w → v) for SD.
    
    Parameters:
        w           : Pre-scattering velocity [km/s]
        v           : Post-scattering velocity [km/s]
        earth_data  : Earth composition data
        radius_idx  : Radius index
        DM_mass     : Dark matter mass [GeV]
        sigma_SD_p  : DM-proton SD cross section [cm²]
        sigma_SD_n  : DM-neutron SD cross section [cm²] (default 0)
        q           : Momentum transfer [GeV] (default 0.01)
    
    Returns:
        rate : Differential scattering rate [s⁻¹]
    """
    v_rel = w
    km_s_to_GeV = 1.0 / 3.0e5
    
    differential_rate_total = 0.0
    m_N = 0.9315
    
    for elem, abundance in earth_data["abundances"].items():
        if elem not in ELEMENT_DB:
            continue
        
        # Check if spin data exists
        if elem not in SD_SPIN_DB:
            continue
        
        spin_data = SD_SPIN_DB[elem]
        J = spin_data["J"]
        
        if J == 0.0:
            continue
        
        S_p = spin_data["S_p"]
        S_n = spin_data["S_n"]
        
        A = ELEMENT_DB[elem]["A"]
        m_A = A * m_N
        
        # Calculate recoil energy
        E_R = velocity_to_recoil_energy(w, v, m_A)
        
        if E_R <= 0:
            continue
        
        # Differential cross section
        dsigma_dE = differential_cross_section_SD(E_R, v_rel, DM_mass, A, J, S_p, S_n, sigma_SD_p, q)
        
        if dsigma_dE == 0.0:
            continue
        
        # Number density [cm⁻³]
        number_density = number_density_at_radius(earth_data, radius_idx, elem)
        
        # Jacobian factor
        jacobian = m_A * v * km_s_to_GeV
        
        # Differential scattering rate
        differential_rate  = number_density * dsigma_dE * v_rel * jacobian
        
        differential_rate_total += differential_rate 
    
    return differential_rate_total


def capture_probability(w, earth_data, radius_idx, DM_mass, sigma_SI_p, sigma_SD_p=0.0, scattering_type="SI", cross_section_type="constant"):
    """
    Calculate the capture probability for a DM particle with incoming velocity w.
    
    Parameters:
        w                  : Incoming velocity [km/s]
        earth_data         : Earth composition data
        radius_idx         : Radius index
        DM_mass            : Dark matter mass [GeV]
        sigma_SI_p         : DM-proton SI cross section [cm²]
        sigma_SD_p         : DM-proton SD cross section [cm²] (default 0)
        scattering_type    : "SI" or "SD" (default "SI")
        cross_section_type : "constant", "v2_dependent", or "q2_dependent" (default "constant")
    
    Returns:
        capture_probability : Probability (dimensionless, between 0 and 1)
    """
    # Escape velocity at this radius [km/s]
    r_km = earth_data["radius"][radius_idx]
    if r_km == 0:
        return 0.0
    r_m = r_km * 1e3  # km → m
    v_esc = np.sqrt(2 * G * M_earth / r_m) / 1000.0  # [km/s]
    
    if w <= v_esc:
        return 1.0  # Already bound
    
    # Use the same number of points for both integrals
    n_points = 200
    
    # Integrate over all possible v (0 to w)
    v_all = np.linspace(0, w, n_points)
    dv_all = v_all[1] - v_all[0]
    total_diff_rate = 0.0
    
    for v in v_all:
        if scattering_type == "SI":
            diff_rate = differential_scattering_rate_SI(
                w, v, earth_data, radius_idx, DM_mass, sigma_SI_p, 
                q=0.01, cross_section_type=cross_section_type
            )
        elif scattering_type == "SD":
            diff_rate = differential_scattering_rate_SD(
                w, v, earth_data, radius_idx, DM_mass, sigma_SD_p,
                sigma_SD_n=0.0, q=0.01, cross_section_type=cross_section_type
            )
        else:
            raise ValueError("scattering_type must be 'SI' or 'SD'")
        total_diff_rate += diff_rate * dv_all
    
    if total_diff_rate == 0.0:
        return 0.0
    
    # Integrate over v < v_esc (capture region)
    v_max_cap = min(w, v_esc)
    v_cap = np.linspace(0, v_max_cap, n_points)
    dv_cap = v_cap[1] - v_cap[0]
    capture_diff_rate = 0.0
    
    for v in v_cap:
        if scattering_type == "SI":
            diff_rate = differential_scattering_rate_SI(
                w, v, earth_data, radius_idx, DM_mass, sigma_SI_p,
                q=0.01, cross_section_type=cross_section_type
            )
        elif scattering_type == "SD":
            diff_rate = differential_scattering_rate_SD(
                w, v, earth_data, radius_idx, DM_mass, sigma_SD_p,
                sigma_SD_n=0.0, q=0.01, cross_section_type=cross_section_type
            )
        else:
            raise ValueError("scattering_type must be 'SI' or 'SD'")
        capture_diff_rate += diff_rate * dv_cap
    
    # Ensure probability is between 0 and 1 (cap at 1.0 due to numerical errors)
    prob = capture_diff_rate / total_diff_rate
    if prob > 1.0:
        prob = 1.0
    elif prob < 0.0:
        prob = 0.0
    
    return prob


def capture_rate_total(earth_data, DM_mass, sigma_SI_p, sigma_SD_p=0.0, scattering_type="SI", u_max=800.0, n_u=200, cross_section_type="constant"):
    """
    Complete version of capture rate calculation (with velocity distribution integration + gravitational focusing)
    
    Parameters:
        earth_data         : Earth composition data
        DM_mass            : Dark matter mass [GeV]
        sigma_SI_p         : DM-proton SI cross section [cm²]
        sigma_SD_p         : DM-proton SD cross section [cm²] (default 0)
        scattering_type    : "SI" or "SD" (default "SI")
        u_max              : Maximum velocity integration range [km/s] (default 800)
        n_u                : Number of velocity bins (default 200)
        cross_section_type : "constant", "v2_dependent", or "q2_dependent" (default "constant")
    """
    # Dark matter number density [cm⁻³]
    rho_chi = 0.3
    n_chi = rho_chi / DM_mass
    
    # Velocity integration range (SHM distribution, v_esc ~ 544 km/s)
    u_range = np.linspace(0, u_max, n_u)
    du = u_range[1] - u_range[0]
    
    total_capture_rate = 0.0
    
    # Progress tracking
    print(f"=== Computing total capture rate ===")
    print(f"DM mass: {DM_mass:.2f} GeV")
    print(f"SI cross section: {sigma_SI_p:.2e} cm²")
    print(f"Scattering type: {scattering_type}")
    print(f"Cross section type: {cross_section_type}")
    print(f"Local DM density: {rho_chi:.2f} GeV/cm³")
    print(f"Number of velocity bins: {n_u}")
    print(f"Number of radius layers: {len(earth_data['radius'])}")
    print(f"n_chi = {n_chi:.6e} cm⁻³")
    print("")
    
    # Debug counters
    debug_count = 0
    max_debug = 3
    
    # --- Loop over DM velocity distribution ---
    for idx_u, u in enumerate(u_range):
        # SHM velocity distribution
        f_v = SHM_velocity_distribution(u)
        
        if f_v == 0.0:
            continue
        
        # --- Loop over radius layers ---
        for r_idx in range(len(earth_data["radius"])):
            r_km = earth_data["radius"][r_idx]
            
            # Escape velocity at this radius [km/s]
            if r_km == 0:
                continue  # Skip center to avoid division by zero
            r_m = r_km * 1e3  # km → m
            v_esc = np.sqrt(2 * G * M_earth / r_m) / 1000.0   # [km/s]
            
            # Gravitational focusing: w(r) = sqrt(u² + v_esc²)
            w = np.sqrt(u**2 + v_esc**2)
            
            if w <= v_esc:
                continue
            
            # Compute capture probability at this radius
            # ✅ 加入 cross_section_type 傳遞
            capture_prob = capture_probability(
                w, earth_data, r_idx, DM_mass, 
                sigma_SI_p, sigma_SD_p, scattering_type, cross_section_type
            )
            
            if capture_prob == 0.0:
                continue
            
            # Volume element [cm³]
            r_cm = r_km * 1e5  # km → cm
            
            # Radial step [cm]
            if r_idx == 0:
                dr_km = earth_data["radius"][1] - earth_data["radius"][0]
            else:
                dr_km = earth_data["radius"][r_idx] - earth_data["radius"][r_idx - 1]
            dr_cm = dr_km * 1e5
            
            dV = 4.0 * np.pi * r_cm**2 * dr_cm  # [cm³]
            
            # Contribution to capture rate [s⁻¹]
            contribution = n_chi * f_v * u * w * dV * capture_prob * du
            
            # Debug output (only for first few iterations)
            if debug_count < max_debug:
                print(f"\n[DEBUG {debug_count+1}]")
                print(f"  u = {u:.2f} km/s, r = {r_km:.1f} km")
                print(f"  v_esc = {v_esc:.2f} km/s, w = {w:.2f} km/s")
                print(f"  f_v = {f_v:.6e} s/km")
                print(f"  capture_prob = {capture_prob:.6e}")
                print(f"  dV = {dV:.6e} cm³")
                print(f"  contribution = {contribution:.6e}")
                debug_count += 1
            
            total_capture_rate += contribution
    
    print(f"\n=== Total capture rate: {total_capture_rate:.6e} s⁻¹ ===")
    print(f"Equivalent to {total_capture_rate * 3.156e7:.2e} yr⁻¹")
    
    return total_capture_rate

def capture_rate_rest_limit(earth_data, DM_mass, sigma_SI_p, sigma_SD_p=0.0, 
                            scattering_type="SI", cross_section_type="constant"):
    """
    Capture rate with target nuclei at rest (T(r) = 0).
    
    This is the limit where the thermal motion of target nuclei is neglected.
    """
    # Simplified: use capture_rate_total with a flag
    # For now, just return the same as capture_rate_total
    # (In a full implementation, this would use the T(r) = 0 approximation)
    return capture_rate_total(earth_data, DM_mass, sigma_SI_p, sigma_SD_p, 
                             scattering_type, u_max=800.0, n_u=50, 
                             cross_section_type=cross_section_type)


def Capture_rate_geometric(xi, v_esc, v_d, DM_mass, R_earth, rho_chi=0.3):
    """
    Calculate the geometrical capture rate.
    
    Formula:
    C_geo = π R² (ρ/m) sqrt(8/(3π)) v_d * xi
    
    Parameters:
        xi      : Xi function value (dimensionless)
        v_esc   : Escape velocity at surface [km/s]
        v_d     : Velocity dispersion parameter [km/s]
        DM_mass : Dark matter mass [GeV]
        R_earth : Earth radius [cm]
        rho_chi : Local DM density [GeV/cm³] (default 0.3)
    """
    n_chi = rho_chi / DM_mass
    sigma_geo = np.pi * R_earth**2
    v_d_cm_s = v_d * 1e5
    coeff = np.sqrt(8 / (3 * np.pi))
    
    C_geo = sigma_geo * n_chi * coeff * v_d_cm_s * xi
    return C_geo
    
def final_capture_rate(C_geo, C_weak):
    """
    Calculate the final capture rate considering both geometric and weak scattering contributions.
    
    Parameters:
        C_geo  : Geometrical capture rate [s⁻¹]
        C_weak : Weak scattering capture rate [s⁻¹]
    """
    ratio = C_geo / C_weak
    
    if ratio > 50:  # exp(-50) ≈ 1.9e-22，已經非常接近 0
        # 完全飽和，捕捉率 = C_weak
        C_earth = C_weak
    elif ratio < 1e-10:  # 極小比值，exp(-small) ≈ 1 - small
        # 泰勒展開，避免數值問題
        C_earth = C_weak * (ratio - ratio**2 / 2 + ratio**3 / 6)
    else:
        C_earth = C_weak * (1.0 - np.exp(-ratio))
    
    return C_earth

def plot_capture_rates(earth_data, DM_masses, sigma_SI_p, sigma_SD_p, sigma_electron, 
                        rho_chi=0.3, v_earth=244.0, v_esc_earth=11.2, 
                        u_max=600.0, n_u=30, save_path="capture_rate_vs_dm_mass.png"):
    """
    Plot capture rate vs DM mass for SI, SD, and electron scattering.
    
    Parameters:
        earth_data      : Earth composition data from load_earth_composition()
        DM_masses       : Array of DM masses [GeV]
        sigma_SI_p      : SI cross section [cm²]
        sigma_SD_p      : SD cross section [cm²]
        sigma_electron  : Electron scattering cross section [cm²]
        rho_chi         : Local DM density [GeV/cm³] (default 0.3)
        v_earth         : Earth velocity in DM halo [km/s] (default 244)
        v_esc_earth     : Earth escape velocity [km/s] (default 11.2)
        u_max           : Max velocity for integration [km/s] (default 600)
        n_u             : Number of velocity bins (default 30)
        save_path       : Output file path (default "capture_rate_vs_dm_mass.png")
    
    Returns:
        C_SI, C_SD, C_electron : Lists of capture rates [s⁻¹]
    """
    
    print("\n" + "=" * 70)
    print("Computing Capture Rates vs DM Mass")
    print("=" * 70)
    
    # Initialize lists to store results
    C_SI = []
    C_SD = []
    C_electron = []
    
    start_time = time.time()
    
    print(f"\nDM mass range: {DM_masses[0]:.3f} to {DM_masses[-1]:.3f} GeV")
    print(f"Number of mass points: {len(DM_masses)}")
    print(f"Velocity bins: {n_u}")
    print("-" * 70)
    
    # Loop over DM masses
    for i, m in enumerate(DM_masses):
        print(f"Progress: {i+1:3d}/{len(DM_masses):3d}  m = {m:8.3f} GeV", end="", flush=True)
        
        # 1. SI capture rate (dot-dashed blue)
        try:
            C_SI_val = capture_rate_total(
                earth_data, m, sigma_SI_p, 
                scattering_type="SI",
                u_max=u_max, n_u=n_u,
                cross_section_type="constant"
            )
        except Exception as e:
            print(f"  [SI error: {e}]")
            C_SI_val = 0.0
        C_SI.append(C_SI_val)
        
        # 2. SD capture rate (dashed green)
        try:
            C_SD_val = capture_rate_total(
                earth_data, m, sigma_SD_p, 
                scattering_type="SD",
                u_max=u_max, n_u=n_u,
                cross_section_type="constant"
            )
        except Exception as e:
            print(f"  [SD error: {e}]")
            C_SD_val = 0.0
        C_SD.append(C_SD_val)
        
        # 3. Electron capture rate (solid red)
        try:
            C_e_val = capture_rate_total(
                earth_data, m, sigma_electron, 
                scattering_type="SI",
                u_max=u_max, n_u=n_u,
                cross_section_type="constant"
            )
            # Empirical scaling for electron scattering
            C_e_val = C_e_val * 0.01
        except Exception as e:
            print(f"  [Electron error: {e}]")
            C_e_val = 0.0
        C_electron.append(C_e_val)
        
        print("  [Done]")
    
    end_time = time.time()
    print(f"\nComputation completed in {end_time - start_time:.2f} seconds")
    
    # ============================================================
    # Generate Plot
    # ============================================================
    print("\n=== Generating Plot ===")
    
    # Calculate geometrical capture rate for reference
    v_d = np.sqrt(3/2) * v_earth
    R_earth_cm = earth_data["radius"][-1] * 1e5
    
    xi = Xi_function(v_d, v_esc_earth, v_earth)
    
    C_geo = Capture_rate_geometric(xi, v_esc_earth, v_d, DM_masses[0], R_earth_cm, rho_chi)
    
    # Create figure
    plt.figure(figsize=(10, 8))
    
    # SI: dot-dashed blue
    plt.plot(DM_masses, C_SI, 'b-.', linewidth=2.5, label='SI (Spin-Independent)', markersize=6)
    
    # SD: dashed green
    plt.plot(DM_masses, C_SD, 'g--', linewidth=2.5, label='SD (Spin-Dependent)', markersize=6)
    
    # Electron: solid red
    plt.plot(DM_masses, C_electron, 'r-', linewidth=2.5, label='DM-Electron', markersize=6)
    
    # Add C_geo as a horizontal reference line
    plt.axhline(y=C_geo, color='gray', linestyle=':', linewidth=1.5, 
                label=f'$C_{{\\mathrm{{geo}}}}$ = {C_geo:.2e} s$^{{-1}}$')
    
    # Plot formatting
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$m_\chi$ [GeV]', fontsize=14)
    plt.ylabel(r'Capture Rate $C$ [s$^{-1}$]', fontsize=14)
    plt.title('Dark Matter Capture Rate vs Mass for Earth', fontsize=16)
    plt.grid(True, alpha=0.3, linestyle='--')
    plt.legend(fontsize=12, loc='best')
    
    # Add parameter info box
    plt.text(0.02, 0.98, 
             f'$\\rho_\\chi$ = {rho_chi:.1f} GeV/cm³\n'
             f'$\\sigma_{SI}$ = {sigma_SI_p:.0e} cm²\n'
             f'$\\sigma_{SD}$ = {sigma_SD_p:.0e} cm²\n'
             f'$\\sigma_{e}$ = {sigma_electron:.0e} cm²',
             transform=plt.gca().transAxes,
             fontsize=10,
             verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved as '{save_path}'")
    plt.show()
    
    # ============================================================
    # Print numerical results
    # ============================================================
    print("\n=== Numerical Results ===")
    print("-" * 80)
    print(f"{'DM Mass [GeV]':<15} {'C_SI [s⁻¹]':<20} {'C_SD [s⁻¹]':<20} {'C_e [s⁻¹]':<20}")
    print("-" * 80)
    for i, m in enumerate(DM_masses):
        print(f"{m:<15.3f} {C_SI[i]:<20.6e} {C_SD[i]:<20.6e} {C_electron[i]:<20.6e}")
    print("-" * 80)
    
    return C_SI, C_SD, C_electron
# ============================================================
# Complete plot with 3x2 panels
# ============================================================

def plot_capture_rates_complete(earth_data, DM_masses, sigma_values, 
                               cross_section_types, save_path="complete_capture_rates.png"):
    """
    Complete capture rate plot with 3 rows x 2 columns.
    
    Left panels: Capture rates for SI, SD, and electron scattering.
    Right panels: Ratio of capture rates with respect to T(r) = 0 limit.
    
    Parameters:
        earth_data          : Earth composition data
        DM_masses           : Array of DM masses [GeV]
        sigma_values        : List of [sigma_SI, sigma_SD, sigma_electron] for each row
        cross_section_types : List of cross section types for each row
        save_path           : Output file path
    """
    fig, axes = plt.subplots(3, 2, figsize=(14, 18))
    
    # Load Earth data for rest limit
    R_earth_cm = earth_data["radius"][-1] * 1e5
    v_earth = 244.0
    v_esc_earth = 11.2
    v_d = np.sqrt(3/2) * v_earth
    xi = Xi_function(v_esc_earth, v_earth)
    
    for row, (sigma, cross_type) in enumerate(zip(sigma_values, cross_section_types)):
        
        sigma_SI_p, sigma_SD_p, sigma_electron = sigma
        
        C_SI = []
        C_SD = []
        C_e = []
        C_rest = []
        
        print(f"\n=== Row {row+1}: {cross_type}, σ = {sigma_SI_p:.0e} cm² ===")
        
        for i, m in enumerate(DM_masses):
            print(f"  Mass {i+1}/{len(DM_masses)}: {m:.3f} GeV", end="", flush=True)
            
            # SI
            C_SI_val = capture_rate_total(
                earth_data, m, sigma_SI_p, 
                scattering_type="SI",
                u_max=600.0, n_u=20,
                cross_section_type=cross_type
            )
            C_SI.append(C_SI_val)
            
            # SD
            C_SD_val = capture_rate_total(
                earth_data, m, sigma_SD_p, 
                scattering_type="SD",
                u_max=600.0, n_u=20,
                cross_section_type=cross_type
            )
            C_SD.append(C_SD_val)
            
            # Electron
            C_e_val = capture_rate_total(
                earth_data, m, sigma_electron, 
                scattering_type="SI",  # Use SI function as placeholder
                u_max=600.0, n_u=20,
                cross_section_type=cross_type
            )
            C_e.append(C_e_val)
            
            # Rest limit
            C_rest_val = capture_rate_rest_limit(
                earth_data, m, sigma_SI_p, 
                scattering_type="SI",
                cross_section_type=cross_type
            )
            C_rest.append(C_rest_val)
            
            print(" [Done]")
        
        # Calculate C_geo for reference
        C_geo = Capture_rate_geometric(xi, v_esc_earth, v_d, DM_masses[0], R_earth_cm, 0.3)
        
        # Left panel: Capture rates
        ax_left = axes[row, 0]
        ax_left.loglog(DM_masses, C_SI, 'b-.', linewidth=2.5, label='SI')
        ax_left.loglog(DM_masses, C_SD, 'g--', linewidth=2.5, label='SD')
        ax_left.loglog(DM_masses, C_e, 'r-', linewidth=2.5, label='Electron')
        ax_left.axhline(y=C_geo, color='k', linestyle=':', linewidth=1.5, label='C_geo')
        ax_left.set_xlabel(r'$m_\chi$ [GeV]', fontsize=12)
        ax_left.set_ylabel(r'$C$ [s$^{-1}$]', fontsize=12)
        ax_left.legend(fontsize=10, loc='best')
        ax_left.grid(True, alpha=0.3)
        ax_left.set_title(f'{cross_type}, σ = {sigma_SI_p:.0e} cm²', fontsize=12)
        
        # Right panel: Ratio to rest limit
        ax_right = axes[row, 1]
        ratio_SI = [c/c_rest if c_rest > 1e-50 else 0 for c, c_rest in zip(C_SI, C_rest)]
        ratio_SD = [c/c_rest if c_rest > 1e-50 else 0 for c, c_rest in zip(C_SD, C_rest)]
        ratio_e = [c/c_rest if c_rest > 1e-50 else 0 for c, c_rest in zip(C_e, C_rest)]
        
        ax_right.semilogx(DM_masses, ratio_SI, 'b-.', linewidth=2, label='SI')
        ax_right.semilogx(DM_masses, ratio_SD, 'g--', linewidth=2, label='SD')
        ax_right.semilogx(DM_masses, ratio_e, 'r-', linewidth=2, label='Electron')
        ax_right.axhline(y=1.0, color='k', linestyle=':', linewidth=1, label='T(r)=0 limit')
        ax_right.set_xlabel(r'$m_\chi$ [GeV]', fontsize=12)
        ax_right.set_ylabel(r'$C / C_{T=0}$', fontsize=12)
        ax_right.legend(fontsize=10, loc='best')
        ax_right.grid(True, alpha=0.3)
        ax_right.set_title(f'Ratio to T(r)=0 limit', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved as '{save_path}'")
    plt.show()

# ============================================================
# Main program
# ============================================================
if __name__ == "__main__":
    
    # 1. Load Earth composition data
    print("Loading Earth data...")
    earth_data = load_earth_composition("data/earth_prem.dat")
    print(f"Loaded {len(earth_data['radius'])} layers, {len(earth_data['abundances'])} elements")
    
    # 2. Set dark matter parameters
    DM_mass = 10.0              # GeV
    sigma_SI_p = 1e-35          # cm²
    sigma_SD_p = 1e-36          # cm²
    sigma_electron = 1e-40      # cm²
    
    # 3. Compute C_weak for a single mass point
    print("\n=== Computing C_weak for DM_mass = 10.0 GeV ===")
    C_weak = capture_rate_total(earth_data, DM_mass, sigma_SI_p)
    print(f"C_weak = {C_weak:.6e} s⁻¹")
    
    # 4. Define DM mass range for plotting
    DM_masses = np.logspace(-2, 3, 15)  # 0.01 to 1000 GeV
    
    # 5. Plot capture rates vs DM mass
    C_SI, C_SD, C_e = plot_capture_rates(
        earth_data=earth_data,
        DM_masses=DM_masses,
        sigma_SI_p=sigma_SI_p,
        sigma_SD_p=sigma_SD_p,
        sigma_electron=sigma_electron,
        rho_chi=0.3,
        v_earth=244.0,
        v_esc_earth=11.2,
        u_max=600.0,
        n_u=30,
        save_path="capture_rate_vs_dm_mass.png"
    )
# ============================================================
# main function
# ============================================================