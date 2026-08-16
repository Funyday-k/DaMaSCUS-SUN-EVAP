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
import os
from concurrent.futures import ProcessPoolExecutor, as_completed


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
# SD spin dataset
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
# notice：
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
# Physical constants
# ============================================================
G = 6.67430e-11                 # m^3 kg^-1 s^-2
M_earth = 5.972e24              # kg
R_earth = 6.371e6               # m

C_LIGHT_KM_S = 299792.458       # km/s
M_U_GEV = 0.93149410242         # atomic mass unit in GeV
M_U_G = 1.66053906660e-24       # atomic mass unit in g
M_P_GEV = 0.9382720813          # proton mass in GeV
M_E_GEV = 0.00051099895         # electron mass in GeV

FM_TO_GEV_INV = 5.067731237     # 1 fm = 5.067731237 GeV^-1

# Halo parameters
RHO_CHI_DEFAULT = 0.3           # GeV/cm^3
V0_DEFAULT = 220.0              # km/s
VESC_HALO_DEFAULT = 544.0       # km/s

K_B_J = 1.380649e-23
GEV_TO_KG = 1.78266192e-27   # 1 GeV/c^2 in kg

_GL_CACHE = {}

# ============================================================
# Basic helpers
# ============================================================

def leggauss_cached(n):
    if n not in _GL_CACHE:
        _GL_CACHE[n] = np.polynomial.legendre.leggauss(n)
    return _GL_CACHE[n]

def thermal_speed_scale_km_s(T_K, m_target_GeV):
    """
    Maxwell-Boltzmann thermal speed parameter:
        f(v) ~ exp(-v^2 / v_th^2)
    where
        v_th = sqrt(2 k_B T / m)
    """
    if T_K <= 0.0 or m_target_GeV <= 0.0:
        return 0.0
    m_kg = m_target_GeV * GEV_TO_KG
    v_th_m_s = np.sqrt(2.0 * K_B_J * T_K / m_kg)
    return v_th_m_s / 1000.0

def maxwell_speed_distribution(u, v_th):
    """
    3D Maxwell speed distribution:
        f(u) du = (4/sqrt(pi)) (u^2/v_th^3) exp(-u^2/v_th^2) du
    normalized on [0, inf)
    """
    u = np.asarray(u, dtype=float)
    if v_th <= 0.0:
        out = np.zeros_like(u)
        out[0] = 1.0
        return out
    return (4.0 / np.sqrt(np.pi)) * (u**2 / v_th**3) * np.exp(-(u**2) / v_th**2)

def reduced_mass(m1, m2):
    return (m1 * m2) / (m1 + m2)

def shm_speed_distribution(u, v0=V0_DEFAULT, vesc=VESC_HALO_DEFAULT):
    """
    Truncated Standard Halo Model speed distribution f(u), normalized so that:
        ∫ f(u) du = 1   over [0, vesc]
    Units:
        if u is in km/s, then f(u) is in (km/s)^-1
    """
    u = np.asarray(u, dtype=float)
    z = vesc / v0
    N = erf(z) - (2.0 / np.sqrt(np.pi)) * z * np.exp(-z**2)

    f = np.zeros_like(u)
    mask = (u >= 0.0) & (u < vesc)
    um = u[mask]

    f[mask] = (4.0 / np.sqrt(np.pi)) * (um**2 / v0**3) * np.exp(-(um**2) / v0**2) / N
    return f

# ============================================================
# Earth composition loader
# ============================================================
def load_earth_composition(filepath="data/earth_prem.dat"):
    """
    Read earth_prem.dat and return:
    {
        "radius": np.array([...])          # km
        "density": np.array([...])         # g/cm^3
        "temperature": np.array([...])     # K
        "abundances": {elem: np.array([...])}  # ppm by mass
    }
    """
    with open(filepath, "r", encoding="utf-8") as f:
        lines = f.readlines()

    data_start = None
    for i, line in enumerate(lines):
        if line.startswith("# The Table begins here."):
            data_start = i + 1
            break

    if data_start is None:
        raise ValueError("Cannot find '# The Table begins here.' in earth_prem.dat")

    header_line = lines[data_start - 2]
    parts = header_line.split()

    elements = []
    for part in parts:
        if "_ppm" in part:
            elements.append(part.replace("_ppm", ""))

    radii = []
    densities = []
    temperatures = []
    abundances = {elem: [] for elem in elements}

    for line in lines[data_start:]:
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.strip().split()
        if len(cols) < 4 + len(elements):
            continue

        radii.append(float(cols[0]))
        densities.append(float(cols[2]))
        temperatures.append(float(cols[3]))

        for i, elem in enumerate(elements):
            abundances[elem].append(float(cols[4 + i]))

    earth_data = {
        "radius": np.array(radii, dtype=float),
        "density": np.array(densities, dtype=float),
        "temperature": np.array(temperatures, dtype=float),
        "abundances": {k: np.array(v, dtype=float) for k, v in abundances.items()}
    }

    precompute_earth_shells(earth_data)
    build_active_shell_data(earth_data)
    return earth_data

def build_active_shell_data(earth_data, min_mass_fraction=1e-10):
    active_shells = []

    for i in range(len(earth_data["radius"])):
        shell_info = {
            "SI": [],
            "SD": [],
            "electron": {
                "n_e": electron_density_at_radius(earth_data, i),
                "T": earth_data["temperature"][i],
                "vesc": earth_data["v_esc_profile"][i],
                "dV": earth_data["shell_volume_cm3"][i],
            }
        }

        for elem in earth_data["abundances"]:
            if elem not in ELEMENT_DB:
                continue

            mass_fraction = earth_data["abundances"][elem][i] / 1e6
            if mass_fraction < min_mass_fraction:
                continue

            n_i = number_density_at_radius(earth_data, i, elem)
            if n_i <= 0.0:
                continue

            A_eff = int(round(ELEMENT_DB[elem]["A"]))
            m_t = A_eff * M_U_GEV

            shell_info["SI"].append({
                "elem": elem,
                "n_i": n_i,
                "m_t": m_t,
            })

            if elem in SD_SPIN_DB:
                spin_data = SD_SPIN_DB[elem]
                J = spin_data.get("J", 0.0)
                if J > 0.0:
                    iso_ab = spin_data.get("isotope_abundance", spin_data.get("abundance", 1.0))
                    A_sd = int(spin_data.get("isotope", round(ELEMENT_DB[elem]["A"])))
                    m_t_sd = A_sd * M_U_GEV

                    shell_info["SD"].append({
                        "elem": elem,
                        "n_i_eff": n_i * iso_ab,
                        "m_t": m_t_sd,
                    })

        active_shells.append(shell_info)

    earth_data["active_shells"] = active_shells
    return earth_data
# ============================================================
# Earth shell / escape velocity profile
# ============================================================
def precompute_earth_shells(earth_data):
    """
    Build shell volumes, enclosed mass, and a finite escape-velocity profile
    inside the Earth.
    """
    r_km = np.asarray(earth_data["radius"], dtype=float)
    rho_g_cm3 = np.asarray(earth_data["density"], dtype=float)

    r_m = r_km * 1e3
    rho_kg_m3 = rho_g_cm3 * 1000.0

    n = len(r_m)
    if n < 2:
        raise ValueError("Need at least 2 radius points in earth_prem.dat")

    # Build shell edges from radius grid
    edges = np.zeros(n + 1)
    edges[0] = 0.0
    edges[1:-1] = 0.5 * (r_m[:-1] + r_m[1:])
    edges[-1] = r_m[-1]

    shell_r_inner = edges[:-1]
    shell_r_outer = edges[1:]
    shell_r_center = 0.5 * (shell_r_inner + shell_r_outer)
    shell_dr = shell_r_outer - shell_r_inner

    shell_vol_m3 = (4.0 / 3.0) * np.pi * (shell_r_outer**3 - shell_r_inner**3)
    shell_mass_kg = rho_kg_m3 * shell_vol_m3
    M_enc = np.cumsum(shell_mass_kg)

    # Potential inside spherical body:
    # Phi(r) = -G [ M_enc(r)/r + ∫_r^R 4π rho(r') r' dr' ]
    outer_integrand = 4.0 * np.pi * rho_kg_m3 * shell_r_center
    outer_term = np.zeros(n)
    running = 0.0
    for i in range(n - 1, -1, -1):
        running += outer_integrand[i] * shell_dr[i]
        outer_term[i] = running

    r_safe = np.maximum(shell_r_center, 1.0)
    phi = -G * (M_enc / r_safe + outer_term)
    v_esc_m_s = np.sqrt(np.maximum(-2.0 * phi, 0.0))
    v_esc_km_s = v_esc_m_s / 1000.0

    earth_data["shell_r_center_m"] = shell_r_center
    earth_data["shell_r_center_cm"] = shell_r_center * 100.0
    earth_data["shell_dr_cm"] = shell_dr * 100.0
    earth_data["shell_volume_cm3"] = shell_vol_m3 * 1e6
    earth_data["shell_mass_kg"] = shell_mass_kg
    earth_data["M_enc_kg"] = M_enc
    earth_data["v_esc_profile"] = v_esc_km_s

# ============================================================
# Number densities
# ============================================================
def number_density_at_radius(earth_data, radius_idx, elem):
    """
    Number density of nuclei [cm^-3]
    """
    if elem not in ELEMENT_DB:
        return 0.0
    if elem not in earth_data["abundances"]:
        return 0.0

    A_mass = ELEMENT_DB[elem]["A"]
    mass_fraction = earth_data["abundances"][elem][radius_idx] / 1e6
    mass_fraction = max(mass_fraction, 0.0)

    rho_g_cm3 = earth_data["density"][radius_idx]
    n_i = rho_g_cm3 * mass_fraction / (A_mass * M_U_G)
    return n_i

def electron_density_at_radius(earth_data, radius_idx):
    """
    Electron number density [cm^-3], assuming neutrality.
    """
    n_e = 0.0
    for elem in earth_data["abundances"]:
        if elem not in ELEMENT_DB:
            continue
        Z = ELEMENT_DB[elem]["Z"]
        n_i = number_density_at_radius(earth_data, radius_idx, elem)
        n_e += Z * n_i
    return n_e

# ============================================================
# Form factor and cross sections
# ============================================================
def nuclear_form_factor_for_nucleus(A_eff, q_GeV, scattering_type="SI"):
    """
    Simple Gaussian form factor.
    q_GeV in GeV
    """
    if scattering_type == "electron":
        return 1.0

    A_eff = max(float(A_eff), 1.0)

    if scattering_type == "SI":
        r_fm = 1.2 * A_eff**(1.0 / 3.0)
    elif scattering_type == "SD":
        r_fm = 1.0 * A_eff**(1.0 / 3.0)
    else:
        raise ValueError("scattering_type must be 'SI', 'SD', or 'electron'")

    x = q_GeV * r_fm * FM_TO_GEV_INV / np.sqrt(3.0)
    return float(np.exp(-x**2))

def cross_section_constant(sigma_0, v_rel, q):
    return sigma_0

def cross_section_v2_dependent(sigma_0, v_rel, q):
    v_ref = 220.0  # km/s
    return sigma_0 * (v_rel / v_ref)**2

def cross_section_q2_dependent(sigma_0, v_rel, q):
    q_ref = 0.04  # GeV
    return sigma_0 * (q / q_ref)**2

def get_cross_section(sigma_0, v_rel, q, cross_section_type="constant"):
    if cross_section_type == "constant":
        return cross_section_constant(sigma_0, v_rel, q)
    elif cross_section_type == "v2_dependent":
        return cross_section_v2_dependent(sigma_0, v_rel, q)
    elif cross_section_type == "q2_dependent":
        return cross_section_q2_dependent(sigma_0, v_rel, q)
    else:
        raise ValueError("cross_section_type must be 'constant', 'v2_dependent', or 'q2_dependent'")

def sigma_nucleus_SI(DM_mass, elem, w_km_s, sigma_SI_p, cross_section_type="constant"):
    A_mass = ELEMENT_DB[elem]["A"]
    A_eff = int(round(A_mass))
    m_A = A_mass * M_U_GEV

    mu_A = reduced_mass(DM_mass, m_A)
    mu_p = reduced_mass(DM_mass, M_P_GEV)

    q_typ = 2.0 * mu_A * (w_km_s / C_LIGHT_KM_S)  # GeV
    F_q = nuclear_form_factor_for_nucleus(A_eff, q_typ, "SI")
    sigma_eff_p = get_cross_section(sigma_SI_p, w_km_s, q_typ, cross_section_type)

    sigma_A = (mu_A / mu_p)**2 * (A_eff**2) * sigma_eff_p * F_q**2
    return sigma_A, m_A

def sigma_nucleus_SD(DM_mass, elem, w_km_s, sigma_SD_p, sigma_SD_n=0.0, cross_section_type="constant"):
    if elem not in SD_SPIN_DB:
        return 0.0, None, 0.0

    spin_data = SD_SPIN_DB[elem]
    J = spin_data.get("J", 0.0)
    if J <= 0.0:
        return 0.0, None, 0.0

    iso_ab = spin_data.get("isotope_abundance", spin_data.get("abundance", 1.0))
    S_p = spin_data.get("S_p", 0.0)
    S_n = spin_data.get("S_n", 0.0)

    A_mass = ELEMENT_DB[elem]["A"]
    A_eff = int(round(A_mass))
    m_A = A_mass * M_U_GEV

    mu_A = reduced_mass(DM_mass, m_A)
    mu_p = reduced_mass(DM_mass, M_P_GEV)

    q_typ = 2.0 * mu_A * (w_km_s / C_LIGHT_KM_S)  # GeV
    F_q = nuclear_form_factor_for_nucleus(A_eff, q_typ, "SD")
    sigma_eff_p = get_cross_section(sigma_SD_p, w_km_s, q_typ, cross_section_type)

    spin_factor = 4.0 * (J + 1.0) / (3.0 * J)
    spin_contrib = (S_p + S_n)

    sigma_A = (mu_A / mu_p)**2 * spin_factor * (spin_contrib**2) * sigma_eff_p * F_q**2
    return sigma_A, m_A, iso_ab

def sigma_electron_effective(DM_mass, w_km_s, sigma_electron, cross_section_type="constant"):
    mu_e = reduced_mass(DM_mass, M_E_GEV)
    q_typ = 2.0 * mu_e * (w_km_s / C_LIGHT_KM_S)
    sigma_eff = get_cross_section(sigma_electron, w_km_s, q_typ, cross_section_type)
    return sigma_eff, M_E_GEV

def sigma_nucleus_SI_q(DM_mass, elem, g_km_s, q_GeV, sigma_SI_p, cross_section_type="constant"):
    """
    SI nuclear cross section at a given relative speed g and momentum transfer q.
    """
    A_mass = ELEMENT_DB[elem]["A"]
    A_eff = int(round(A_mass))
    m_A = A_eff * M_U_GEV

    mu_A = reduced_mass(DM_mass, m_A)
    mu_p = reduced_mass(DM_mass, M_P_GEV)

    F_q = nuclear_form_factor_for_nucleus(A_eff, q_GeV, scattering_type="SI")
    sigma_p_eff = get_cross_section(sigma_SI_p, g_km_s, q_GeV, cross_section_type)

    sigma_A = (mu_A / mu_p)**2 * (A_eff**2) * sigma_p_eff * F_q**2
    return sigma_A, m_A

def sigma_nucleus_SD_q(DM_mass, elem, g_km_s, q_GeV, sigma_SD_p, sigma_SD_n=0.0, cross_section_type="constant"):
    """
    SD nuclear cross section at a given relative speed g and momentum transfer q.
    """
    if elem not in SD_SPIN_DB:
        return 0.0, None, 0.0

    spin_data = SD_SPIN_DB[elem]
    J = spin_data.get("J", 0.0)
    if J <= 0.0:
        return 0.0, None, 0.0

    iso_ab = spin_data.get("isotope_abundance", spin_data.get("abundance", 1.0))
    S_p = spin_data.get("S_p", 0.0)
    S_n = spin_data.get("S_n", 0.0)

    A_eff = int(spin_data.get("isotope", round(ELEMENT_DB[elem]["A"])))
    m_A = A_eff * M_U_GEV

    mu_A = reduced_mass(DM_mass, m_A)
    mu_p = reduced_mass(DM_mass, M_P_GEV)

    F_q = nuclear_form_factor_for_nucleus(A_eff, q_GeV, scattering_type="SD")
    sigma_p_eff = get_cross_section(sigma_SD_p, g_km_s, q_GeV, cross_section_type)

    spin_factor = 4.0 * (J + 1.0) / (3.0 * J)
    spin_contrib = (S_p + S_n)

    sigma_A = (mu_A / mu_p)**2 * spin_factor * (spin_contrib**2) * sigma_p_eff * F_q**2
    return sigma_A, m_A, iso_ab

def sigma_electron_q(DM_mass, g_km_s, q_GeV, sigma_electron, cross_section_type="constant"):
    sigma_e = get_cross_section(sigma_electron, g_km_s, q_GeV, cross_section_type)
    return sigma_e, M_E_GEV

def sigma_capture_average_over_final_state(
    w_km_s,
    u_t_km_s,
    mu_in,
    m_chi,
    m_t,
    vesc_km_s,
    sigma_eval,
    n_scatter_mu=8,
    n_scatter_phi=12
):
    """
    Average over CM-frame final-state directions:
        < sigma(g,q) * Theta(v_esc - v'_chi) >_{Omega_*}

    Parameters
    ----------
    w_km_s : DM speed in lab
    u_t_km_s : target thermal speed in lab
    mu_in : cos(angle between target velocity and DM velocity)
    m_chi : DM mass [GeV]
    m_t : target mass [GeV]
    vesc_km_s : local escape speed [km/s]
    sigma_eval : function(g_km_s, q_GeV) -> sigma [cm^2]
    """
    # choose DM direction along +z
    sin_in = np.sqrt(max(0.0, 1.0 - mu_in**2))
    w_vec = np.array([0.0, 0.0, w_km_s])
    u_vec = np.array([u_t_km_s * sin_in, 0.0, u_t_km_s * mu_in])

    g_vec = w_vec - u_vec
    g_km_s = np.linalg.norm(g_vec)
    if g_km_s < 1e-12:
        return 0.0

    Mtot = m_chi + m_t
    alpha = m_t / Mtot
    mu_red = reduced_mass(m_chi, m_t)

    V_vec = (m_chi * w_vec + m_t * u_vec) / Mtot
    V_km_s = np.linalg.norm(V_vec)

    mu_s_nodes, mu_s_weights = leggauss_cached(n_scatter_mu)
    phi_nodes = 2.0 * np.pi * (np.arange(n_scatter_phi) + 0.5) / n_scatter_phi
    cosphi = np.cos(phi_nodes)[None, :]

    # momentum transfer for a given CM scattering angle theta*
    # q = mu_red * (g/c) * sqrt(2(1-cos theta*))
    q_vals = mu_red * (g_km_s / C_LIGHT_KM_S) * np.sqrt(2.0 * np.maximum(0.0, 1.0 - mu_s_nodes))
    sigma_vals = np.array([sigma_eval(g_km_s, q) for q in q_vals])

    # If V = 0, final speed is angle-independent
    if V_km_s < 1e-12:
        vprime2 = (alpha * g_km_s)**2
        captured = 1.0 if vprime2 <= vesc_km_s**2 else 0.0
        return captured * 0.5 * np.sum(mu_s_weights * sigma_vals)

    cos_beta = np.dot(g_vec, V_vec) / (g_km_s * V_km_s)
    cos_beta = np.clip(cos_beta, -1.0, 1.0)
    sin_beta = np.sqrt(max(0.0, 1.0 - cos_beta**2))

    mu_s = mu_s_nodes[:, None]
    sin_s = np.sqrt(np.maximum(0.0, 1.0 - mu_s**2))

    # cos(psi) = cos(beta) cos(theta*) + sin(beta) sin(theta*) cos(phi*)
    cos_psi = cos_beta * mu_s + sin_beta * sin_s * cosphi

    vprime2 = V_km_s**2 + (alpha * g_km_s)**2 + 2.0 * alpha * g_km_s * V_km_s * cos_psi
    cap_frac_vs_mu_s = np.mean(vprime2 <= vesc_km_s**2, axis=1)

    return 0.5 * np.sum(mu_s_weights * sigma_vals * cap_frac_vs_mu_s)

def thermal_average_gsigma_capture(
    w_km_s,
    vesc_km_s,
    T_K,
    m_chi,
    m_t,
    sigma_eval,
    include_thermal_targets=True,
    n_t_speed=6,
    n_t_mu=6,
    n_scatter_mu=8,
    n_scatter_phi=12,
    u_t_max_mult=5.0
):
    """
    Compute:
        < g * <sigma * Theta(capture)>_{Omega_*} >_{thermal target velocities}

    Returns
    -------
    avg_gsigma_capture : [cm^3 / s]
    """
    # T=0 limit
    if (not include_thermal_targets) or T_K <= 0.0:
        sigma_cap = sigma_capture_average_over_final_state(
            w_km_s=w_km_s,
            u_t_km_s=0.0,
            mu_in=1.0,
            m_chi=m_chi,
            m_t=m_t,
            vesc_km_s=vesc_km_s,
            sigma_eval=sigma_eval,
            n_scatter_mu=n_scatter_mu,
            n_scatter_phi=n_scatter_phi
        )
        return (w_km_s * 1e5) * sigma_cap

    v_th = thermal_speed_scale_km_s(T_K, m_t)
    if v_th < 1e-10:
        sigma_cap = sigma_capture_average_over_final_state(
            w_km_s=w_km_s,
            u_t_km_s=0.0,
            mu_in=1.0,
            m_chi=m_chi,
            m_t=m_t,
            vesc_km_s=vesc_km_s,
            sigma_eval=sigma_eval,
            n_scatter_mu=n_scatter_mu,
            n_scatter_phi=n_scatter_phi
        )
        return (w_km_s * 1e5) * sigma_cap

    u_nodes, u_weights = leggauss_cached(n_t_speed)
    mu_nodes, mu_weights = leggauss_cached(n_t_mu)

    u_max = u_t_max_mult * v_th
    u_grid = 0.5 * (u_nodes + 1.0) * u_max
    du_weights = 0.5 * u_max * u_weights
    f_u = maxwell_speed_distribution(u_grid, v_th)

    total = 0.0

    for uk, wuk, fuk in zip(u_grid, du_weights, f_u):
        inner_mu = 0.0
        for mu_in, wmu in zip(mu_nodes, mu_weights):
            sigma_cap = sigma_capture_average_over_final_state(
                w_km_s=w_km_s,
                u_t_km_s=uk,
                mu_in=mu_in,
                m_chi=m_chi,
                m_t=m_t,
                vesc_km_s=vesc_km_s,
                sigma_eval=sigma_eval,
                n_scatter_mu=n_scatter_mu,
                n_scatter_phi=n_scatter_phi
            )

            g_km_s = np.sqrt(w_km_s**2 + uk**2 - 2.0 * w_km_s * uk * mu_in)
            inner_mu += 0.5 * wmu * (g_km_s * 1e5) * sigma_cap

        total += wuk * fuk * inner_mu

    return total

# ============================================================
# Single-scatter capture probability
# ============================================================
def single_scatter_capture_probability(u_km_s, w_km_s, m_chi, m_target):
    """
    For an elastic scatter off a target initially at rest:
    max fractional energy loss:
        beta_plus = 4 m_chi m_T / (m_chi + m_T)^2

    capture condition:
        Delta >= u^2 / w^2
    """
    if w_km_s <= 0.0:
        return 0.0

    beta_plus = 4.0 * m_chi * m_target / (m_chi + m_target)**2
    if beta_plus <= 0.0:
        return 0.0

    delta_req = (u_km_s**2) / (w_km_s**2)

    if delta_req >= beta_plus:
        return 0.0

    p = 1.0 - delta_req / beta_plus
    return float(np.clip(p, 0.0, 1.0))

# ============================================================
# Local capture kernel Ω^-(w)
# ============================================================
def omega_minus_at_shell(
    earth_data,
    radius_idx,
    u_km_s,
    w_km_s,
    DM_mass,
    sigma_SI_p=0.0,
    sigma_SD_p=0.0,
    sigma_electron=0.0,
    scattering_type="SI",
    cross_section_type="constant",
    include_thermal_targets=True,
    n_t_speed=6,
    n_t_mu=6,
    n_scatter_mu=8,
    n_scatter_phi=12,
    min_mass_fraction=1e-10
):
    """
    Local capture kernel Ω^-(w,r,T) [s^-1]

    Parameters
    ----------
    earth_data : dict
        Earth model data.
    radius_idx : int
        Shell index.
    u_km_s : float
        Halo DM speed at infinity [km/s].
    w_km_s : float
        Local DM speed at radius r [km/s].
    DM_mass : float
        Dark matter mass [GeV].
    sigma_SI_p, sigma_SD_p, sigma_electron : float
        Reference cross sections [cm^2].
    scattering_type : str
        "SI", "SD", or "electron".
    cross_section_type : str
        "constant", "v2_dependent", or "q2_dependent".
    include_thermal_targets : bool
        True  -> finite-T target thermal motion included
        False -> T=0 limit
    n_t_speed, n_t_mu : int
        Quadrature order for target thermal averaging.
    n_scatter_mu, n_scatter_phi : int
        Quadrature order for final-state angular average.
    min_mass_fraction : float
        Skip tiny abundances below this threshold.

    Returns
    -------
    omega : float
        Local capture kernel [s^-1]
    """
    omega = 0.0

    # ------------------------------------------------------------
    # Fast path: use precomputed active_shells if available
    # ------------------------------------------------------------
    use_active_shells = "active_shells" in earth_data and earth_data["active_shells"] is not None

    if use_active_shells:
        shell_info = earth_data["active_shells"][radius_idx]
        T_loc = shell_info["electron"]["T"]
        vesc_loc = shell_info["electron"]["vesc"]

        if scattering_type == "SI":
            for info in shell_info["SI"]:
                elem = info["elem"]
                n_i = info["n_i"]
                m_t = info["m_t"]

                if n_i <= 0.0:
                    continue

                def sigma_eval(g_km_s, q_GeV, elem=elem):
                    sigma_A, _ = sigma_nucleus_SI_q(
                        DM_mass=DM_mass,
                        elem=elem,
                        g_km_s=g_km_s,
                        q_GeV=q_GeV,
                        sigma_SI_p=sigma_SI_p,
                        cross_section_type=cross_section_type
                    )
                    return sigma_A

                avg_gsigma_cap = thermal_average_gsigma_capture(
                    w_km_s=w_km_s,
                    vesc_km_s=vesc_loc,
                    T_K=T_loc,
                    m_chi=DM_mass,
                    m_t=m_t,
                    sigma_eval=sigma_eval,
                    include_thermal_targets=include_thermal_targets,
                    n_t_speed=n_t_speed,
                    n_t_mu=n_t_mu,
                    n_scatter_mu=n_scatter_mu,
                    n_scatter_phi=n_scatter_phi
                )

                omega += n_i * avg_gsigma_cap

        elif scattering_type == "SD":
            for info in shell_info["SD"]:
                elem = info["elem"]
                n_i_eff = info["n_i_eff"]
                m_t = info["m_t"]

                if n_i_eff <= 0.0:
                    continue

                def sigma_eval(g_km_s, q_GeV, elem=elem):
                    sigma_A, _, _ = sigma_nucleus_SD_q(
                        DM_mass=DM_mass,
                        elem=elem,
                        g_km_s=g_km_s,
                        q_GeV=q_GeV,
                        sigma_SD_p=sigma_SD_p,
                        sigma_SD_n=0.0,
                        cross_section_type=cross_section_type
                    )
                    return sigma_A

                avg_gsigma_cap = thermal_average_gsigma_capture(
                    w_km_s=w_km_s,
                    vesc_km_s=vesc_loc,
                    T_K=T_loc,
                    m_chi=DM_mass,
                    m_t=m_t,
                    sigma_eval=sigma_eval,
                    include_thermal_targets=include_thermal_targets,
                    n_t_speed=n_t_speed,
                    n_t_mu=n_t_mu,
                    n_scatter_mu=n_scatter_mu,
                    n_scatter_phi=n_scatter_phi
                )

                omega += n_i_eff * avg_gsigma_cap

        elif scattering_type == "electron":
            n_e = shell_info["electron"]["n_e"]
            if n_e > 0.0:
                m_t = M_E_GEV

                def sigma_eval(g_km_s, q_GeV):
                    sigma_e, _ = sigma_electron_q(
                        DM_mass=DM_mass,
                        g_km_s=g_km_s,
                        q_GeV=q_GeV,
                        sigma_electron=sigma_electron,
                        cross_section_type=cross_section_type
                    )
                    return sigma_e

                avg_gsigma_cap = thermal_average_gsigma_capture(
                    w_km_s=w_km_s,
                    vesc_km_s=vesc_loc,
                    T_K=T_loc,
                    m_chi=DM_mass,
                    m_t=m_t,
                    sigma_eval=sigma_eval,
                    include_thermal_targets=include_thermal_targets,
                    n_t_speed=n_t_speed,
                    n_t_mu=n_t_mu,
                    n_scatter_mu=n_scatter_mu,
                    n_scatter_phi=n_scatter_phi,
                    u_t_max_mult=6.0
                )

                omega += n_e * avg_gsigma_cap

        else:
            raise ValueError("scattering_type must be 'SI', 'SD', or 'electron'")

        return omega

    # ------------------------------------------------------------
    # Fallback path: no active_shells, scan original abundances
    # ------------------------------------------------------------
    T_loc = earth_data["temperature"][radius_idx]
    vesc_loc = earth_data["v_esc_profile"][radius_idx]

    if scattering_type == "SI":
        for elem in earth_data["abundances"]:
            if elem not in ELEMENT_DB:
                continue

            mass_fraction = earth_data["abundances"][elem][radius_idx] / 1e6
            if mass_fraction < min_mass_fraction:
                continue

            n_i = number_density_at_radius(earth_data, radius_idx, elem)
            if n_i <= 0.0:
                continue

            A_eff = int(round(ELEMENT_DB[elem]["A"]))
            m_t = A_eff * M_U_GEV

            def sigma_eval(g_km_s, q_GeV, elem=elem):
                sigma_A, _ = sigma_nucleus_SI_q(
                    DM_mass=DM_mass,
                    elem=elem,
                    g_km_s=g_km_s,
                    q_GeV=q_GeV,
                    sigma_SI_p=sigma_SI_p,
                    cross_section_type=cross_section_type
                )
                return sigma_A

            avg_gsigma_cap = thermal_average_gsigma_capture(
                w_km_s=w_km_s,
                vesc_km_s=vesc_loc,
                T_K=T_loc,
                m_chi=DM_mass,
                m_t=m_t,
                sigma_eval=sigma_eval,
                include_thermal_targets=include_thermal_targets,
                n_t_speed=n_t_speed,
                n_t_mu=n_t_mu,
                n_scatter_mu=n_scatter_mu,
                n_scatter_phi=n_scatter_phi
            )

            omega += n_i * avg_gsigma_cap

    elif scattering_type == "SD":
        for elem in earth_data["abundances"]:
            if elem not in ELEMENT_DB or elem not in SD_SPIN_DB:
                continue

            mass_fraction = earth_data["abundances"][elem][radius_idx] / 1e6
            if mass_fraction < min_mass_fraction:
                continue

            spin_data = SD_SPIN_DB[elem]
            J = spin_data.get("J", 0.0)
            if J <= 0.0:
                continue

            n_i = number_density_at_radius(earth_data, radius_idx, elem)
            if n_i <= 0.0:
                continue

            iso_ab = spin_data.get("isotope_abundance", spin_data.get("abundance", 1.0))
            n_i_eff = n_i * iso_ab

            A_eff = int(spin_data.get("isotope", round(ELEMENT_DB[elem]["A"])))
            m_t = A_eff * M_U_GEV

            def sigma_eval(g_km_s, q_GeV, elem=elem):
                sigma_A, _, _ = sigma_nucleus_SD_q(
                    DM_mass=DM_mass,
                    elem=elem,
                    g_km_s=g_km_s,
                    q_GeV=q_GeV,
                    sigma_SD_p=sigma_SD_p,
                    sigma_SD_n=0.0,
                    cross_section_type=cross_section_type
                )
                return sigma_A

            avg_gsigma_cap = thermal_average_gsigma_capture(
                w_km_s=w_km_s,
                vesc_km_s=vesc_loc,
                T_K=T_loc,
                m_chi=DM_mass,
                m_t=m_t,
                sigma_eval=sigma_eval,
                include_thermal_targets=include_thermal_targets,
                n_t_speed=n_t_speed,
                n_t_mu=n_t_mu,
                n_scatter_mu=n_scatter_mu,
                n_scatter_phi=n_scatter_phi
            )

            omega += n_i_eff * avg_gsigma_cap

    elif scattering_type == "electron":
        n_e = electron_density_at_radius(earth_data, radius_idx)
        if n_e > 0.0:
            m_t = M_E_GEV

            def sigma_eval(g_km_s, q_GeV):
                sigma_e, _ = sigma_electron_q(
                    DM_mass=DM_mass,
                    g_km_s=g_km_s,
                    q_GeV=q_GeV,
                    sigma_electron=sigma_electron,
                    cross_section_type=cross_section_type
                )
                return sigma_e

            avg_gsigma_cap = thermal_average_gsigma_capture(
                w_km_s=w_km_s,
                vesc_km_s=vesc_loc,
                T_K=T_loc,
                m_chi=DM_mass,
                m_t=m_t,
                sigma_eval=sigma_eval,
                include_thermal_targets=include_thermal_targets,
                n_t_speed=n_t_speed,
                n_t_mu=n_t_mu,
                n_scatter_mu=n_scatter_mu,
                n_scatter_phi=n_scatter_phi,
                u_t_max_mult=6.0
            )

            omega += n_e * avg_gsigma_cap

    else:
        raise ValueError("scattering_type must be 'SI', 'SD', or 'electron'")

    return omega
    
# ============================================================
# Geometric capture reference
# ============================================================
def capture_rate_geometric(DM_mass, R_earth_cm, rho_chi=RHO_CHI_DEFAULT, v0=V0_DEFAULT, vesc_surface=11.2):
    """
    Simple geometric reference:
        C_geo ~ pi R^2 n_chi <u> (1 + vesc^2 / <u>^2)
    """
    n_chi = rho_chi / DM_mass
    u_bar = 2.0 * v0 / np.sqrt(np.pi)      # km/s
    focus = 1.0 + (vesc_surface / u_bar)**2
    return np.pi * R_earth_cm**2 * n_chi * (u_bar * 1e5) * focus

# ============================================================
# Total capture rate
# ============================================================
def capture_rate_total(
    earth_data,
    DM_mass,
    sigma_SI_p=0.0,
    sigma_SD_p=0.0,
    sigma_electron=0.0,
    scattering_type="SI",
    cross_section_type="constant",
    rho_chi=RHO_CHI_DEFAULT,
    u_max=800.0,
    n_u=80,
    v0=V0_DEFAULT,
    vesc_halo=VESC_HALO_DEFAULT,
    include_thermal_targets=True,
    n_t_speed=6,
    n_t_mu=6,
    n_scatter_mu=8,
    n_scatter_phi=12,
    min_mass_fraction=1e-10
):
    """
    Total capture rate:
        C = n_chi ∫ dV ∫ du f(u) (w/u) Ω^-(w,r,T)

    If include_thermal_targets=False, this becomes the T(r)=0 limit.
    """
    n_chi = rho_chi / DM_mass
    u_grid = np.linspace(1e-3, u_max, n_u)
    f_u = shm_speed_distribution(u_grid, v0=v0, vesc=vesc_halo)

    total_C = 0.0

    for i in range(len(earth_data["radius"])):
        dV = earth_data["shell_volume_cm3"][i]
        vesc_loc = earth_data["v_esc_profile"][i]

        integrand = np.zeros_like(u_grid)

        for j, u in enumerate(u_grid):
            if f_u[j] <= 0.0:
                continue

            w = np.sqrt(u**2 + vesc_loc**2)

            omega = omega_minus_at_shell(
                earth_data=earth_data,
                radius_idx=i,
                u_km_s=u,
                w_km_s=w,
                DM_mass=DM_mass,
                sigma_SI_p=sigma_SI_p,
                sigma_SD_p=sigma_SD_p,
                sigma_electron=sigma_electron,
                scattering_type=scattering_type,
                cross_section_type=cross_section_type,
                include_thermal_targets=include_thermal_targets,
                n_t_speed=n_t_speed,
                n_t_mu=n_t_mu,
                n_scatter_mu=n_scatter_mu,
                n_scatter_phi=n_scatter_phi,
                min_mass_fraction=min_mass_fraction
            )

            integrand[j] = f_u[j] * (w / u) * omega

        total_C += n_chi * dV * np.trapezoid(integrand, u_grid)

    return total_C

def compute_one_mass_point(task):
    (
        earth_data,
        m,
        sigma_SI_p,
        sigma_SD_p,
        sigma_electron,
        cross_type,
        rho_chi,
        u_max,
        n_u,
        v0,
        n_t_speed,
        n_t_mu,
        n_scatter_mu,
        n_scatter_phi
    ) = task

    c_si_T = capture_rate_total(
        earth_data=earth_data,
        DM_mass=m,
        sigma_SI_p=sigma_SI_p,
        scattering_type="SI",
        cross_section_type=cross_type,
        rho_chi=rho_chi,
        u_max=u_max,
        n_u=n_u,
        v0=v0,
        include_thermal_targets=True,
        n_t_speed=n_t_speed,
        n_t_mu=n_t_mu,
        n_scatter_mu=n_scatter_mu,
        n_scatter_phi=n_scatter_phi
    )

    c_si_0 = capture_rate_total(
        earth_data=earth_data,
        DM_mass=m,
        sigma_SI_p=sigma_SI_p,
        scattering_type="SI",
        cross_section_type=cross_type,
        rho_chi=rho_chi,
        u_max=u_max,
        n_u=n_u,
        v0=v0,
        include_thermal_targets=False,
        n_t_speed=n_t_speed,
        n_t_mu=n_t_mu,
        n_scatter_mu=n_scatter_mu,
        n_scatter_phi=n_scatter_phi
    )

    c_sd_T = capture_rate_total(
        earth_data=earth_data,
        DM_mass=m,
        sigma_SD_p=sigma_SD_p,
        scattering_type="SD",
        cross_section_type=cross_type,
        rho_chi=rho_chi,
        u_max=u_max,
        n_u=n_u,
        v0=v0,
        include_thermal_targets=True,
        n_t_speed=n_t_speed,
        n_t_mu=n_t_mu,
        n_scatter_mu=n_scatter_mu,
        n_scatter_phi=n_scatter_phi
    )

    c_sd_0 = capture_rate_total(
        earth_data=earth_data,
        DM_mass=m,
        sigma_SD_p=sigma_SD_p,
        scattering_type="SD",
        cross_section_type=cross_type,
        rho_chi=rho_chi,
        u_max=u_max,
        n_u=n_u,
        v0=v0,
        include_thermal_targets=False,
        n_t_speed=n_t_speed,
        n_t_mu=n_t_mu,
        n_scatter_mu=n_scatter_mu,
        n_scatter_phi=n_scatter_phi
    )

    c_e_T = capture_rate_total(
        earth_data=earth_data,
        DM_mass=m,
        sigma_electron=sigma_electron,
        scattering_type="electron",
        cross_section_type=cross_type,
        rho_chi=rho_chi,
        u_max=u_max,
        n_u=n_u,
        v0=v0,
        include_thermal_targets=True,
        n_t_speed=n_t_speed,
        n_t_mu=n_t_mu,
        n_scatter_mu=n_scatter_mu,
        n_scatter_phi=n_scatter_phi
    )

    c_e_0 = capture_rate_total(
        earth_data=earth_data,
        DM_mass=m,
        sigma_electron=sigma_electron,
        scattering_type="electron",
        cross_section_type=cross_type,
        rho_chi=rho_chi,
        u_max=u_max,
        n_u=n_u,
        v0=v0,
        include_thermal_targets=False,
        n_t_speed=n_t_speed,
        n_t_mu=n_t_mu,
        n_scatter_mu=n_scatter_mu,
        n_scatter_phi=n_scatter_phi
    )

    return {
        "m": m,
        "C_SI_T": c_si_T,
        "C_SI_0": c_si_0,
        "C_SD_T": c_sd_T,
        "C_SD_0": c_sd_0,
        "C_e_T": c_e_T,
        "C_e_0": c_e_0,
    }

# ============================================================
# Plotting
# ============================================================
def plot_capture_rates_complete(
    earth_data,
    DM_masses,
    sigma_values,
    cross_section_types,
    save_path="complete_capture_rates_fixed.png",
    u_max=800.0,
    n_u=250,
    rho_chi=RHO_CHI_DEFAULT,
    v0=V0_DEFAULT
):
    """
    3x2 plot:
      Left  : C(m_chi)
      Right : C / C_geo(m_chi)

    注意：
    右欄不再畫假的 "T(r)=0 limit"，
    因為原程式沒有 thermal target model。
    """
    fig, axes = plt.subplots(3, 2, figsize=(14, 18))
    R_earth_cm = earth_data["radius"][-1] * 1e5

    for row, (sigma_triplet, cross_type) in enumerate(zip(sigma_values, cross_section_types)):
        sigma_SI_p, sigma_SD_p, sigma_electron = sigma_triplet

        C_SI = []
        C_SD = []
        C_e = []
        C_geo = []

        print(f"\n=== Row {row+1}: {cross_type} ===")
        print(f"σ_SI = {sigma_SI_p:.2e} cm², σ_SD = {sigma_SD_p:.2e} cm², σ_e = {sigma_electron:.2e} cm²")

        t0 = time.time()

        for k, m in enumerate(DM_masses):
            print(f"  [{k+1:2d}/{len(DM_masses):2d}] m_chi = {m:10.4g} GeV", flush=True)

            c_si = capture_rate_total(
                earth_data=earth_data,
                DM_mass=m,
                sigma_SI_p=sigma_SI_p,
                scattering_type="SI",
                cross_section_type=cross_type,
                rho_chi=rho_chi,
                u_max=u_max,
                n_u=n_u,
                v0=v0
            )

            c_sd = capture_rate_total(
                earth_data=earth_data,
                DM_mass=m,
                sigma_SD_p=sigma_SD_p,
                scattering_type="SD",
                cross_section_type=cross_type,
                rho_chi=rho_chi,
                u_max=u_max,
                n_u=n_u,
                v0=v0
            )

            c_e = capture_rate_total(
                earth_data=earth_data,
                DM_mass=m,
                sigma_electron=sigma_electron,
                scattering_type="electron",
                cross_section_type=cross_type,
                rho_chi=rho_chi,
                u_max=u_max,
                n_u=n_u,
                v0=v0
            )

            c_geo = capture_rate_geometric(
                DM_mass=m,
                R_earth_cm=R_earth_cm,
                rho_chi=rho_chi,
                v0=v0,
                vesc_surface=11.2
            )

            C_SI.append(c_si)
            C_SD.append(c_sd)
            C_e.append(c_e)
            C_geo.append(c_geo)

        print(f"  Row {row+1} done in {time.time() - t0:.2f} s")

        C_SI = np.array(C_SI)
        C_SD = np.array(C_SD)
        C_e = np.array(C_e)
        C_geo = np.array(C_geo)

        # ---------------- Left panel ----------------
        axL = axes[row, 0]
        axL.loglog(DM_masses, C_SI, "b-.", lw=2.5, label="SI")
        axL.loglog(DM_masses, C_SD, "g--", lw=2.5, label="SD")
        axL.loglog(DM_masses, C_e, "r-", lw=2.5, label="Electron")
        axL.loglog(DM_masses, C_geo, "k:", lw=1.8, label=r"$C_{\rm geo}$")

        axL.set_xlabel(r"$m_\chi$ [GeV]", fontsize=12)
        axL.set_ylabel(r"$C$ [s$^{-1}$]", fontsize=12)
        axL.set_title(f"{cross_type}, sigma = {sigma_SI_p:.0e} cm$^2$", fontsize=12)
        axL.grid(True, alpha=0.3)
        axL.legend(fontsize=10, loc="best")

        # ---------------- Right panel ----------------
        axR = axes[row, 1]
        ratio_SI = C_SI / np.maximum(C_geo, 1e-300)
        ratio_SD = C_SD / np.maximum(C_geo, 1e-300)
        ratio_e = C_e / np.maximum(C_geo, 1e-300)

        axR.semilogx(DM_masses, ratio_SI, "b-.", lw=2.0, label="SI / C_geo")
        axR.semilogx(DM_masses, ratio_SD, "g--", lw=2.0, label="SD / C_geo")
        axR.semilogx(DM_masses, ratio_e, "r-", lw=2.0, label="Electron / C_geo")
        axR.axhline(1.0, color="k", linestyle=":", lw=1.2, label=r"$C=C_{\rm geo}$")

        axR.set_xlabel(r"$m_\chi$ [GeV]", fontsize=12)
        axR.set_ylabel(r"$C / C_{\rm geo}$", fontsize=12)
        axR.set_title("Ratio to geometric limit", fontsize=12)
        axR.grid(True, alpha=0.3)
        axR.legend(fontsize=10, loc="best")

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"\nPlot saved as: {save_path}")
    plt.show()

def plot_capture_rates_complete_thermal(
    earth_data,
    DM_masses,
    sigma_values,
    cross_section_types,
    save_path="complete_capture_rates_thermal.png",
    u_max=800.0,
    n_u=80,
    rho_chi=RHO_CHI_DEFAULT,
    v0=V0_DEFAULT,
    n_t_speed=6,
    n_t_mu=6,
    n_scatter_mu=8,
    n_scatter_phi=12
):
    """
    3x2 plot:
    Left  : C(T != 0)
    Right : C(T != 0) / C(T = 0)
    """
    fig, axes = plt.subplots(3, 2, figsize=(14, 18))

    for row, (sigma_triplet, cross_type) in enumerate(zip(sigma_values, cross_section_types)):
        sigma_SI_p, sigma_SD_p, sigma_electron = sigma_triplet

        C_SI_T, C_SD_T, C_e_T = [], [], []
        C_SI_0, C_SD_0, C_e_0 = [], [], []

        print(f"\n=== Row {row+1}: {cross_type} ===")

        for k, m in enumerate(DM_masses):
            print(f"  [{k+1:2d}/{len(DM_masses):2d}] m_chi = {m:10.4g} GeV", flush=True)

            # SI
            c_si_T = capture_rate_total(
                earth_data=earth_data,
                DM_mass=m,
                sigma_SI_p=sigma_SI_p,
                scattering_type="SI",
                cross_section_type=cross_type,
                rho_chi=rho_chi,
                u_max=u_max,
                n_u=n_u,
                v0=v0,
                include_thermal_targets=True,
                n_t_speed=n_t_speed,
                n_t_mu=n_t_mu,
                n_scatter_mu=n_scatter_mu,
                n_scatter_phi=n_scatter_phi
            )

            c_si_0 = capture_rate_total(
                earth_data=earth_data,
                DM_mass=m,
                sigma_SI_p=sigma_SI_p,
                scattering_type="SI",
                cross_section_type=cross_type,
                rho_chi=rho_chi,
                u_max=u_max,
                n_u=n_u,
                v0=v0,
                include_thermal_targets=False,
                n_t_speed=n_t_speed,
                n_t_mu=n_t_mu,
                n_scatter_mu=n_scatter_mu,
                n_scatter_phi=n_scatter_phi
            )

            # SD
            c_sd_T = capture_rate_total(
                earth_data=earth_data,
                DM_mass=m,
                sigma_SD_p=sigma_SD_p,
                scattering_type="SD",
                cross_section_type=cross_type,
                rho_chi=rho_chi,
                u_max=u_max,
                n_u=n_u,
                v0=v0,
                include_thermal_targets=True,
                n_t_speed=n_t_speed,
                n_t_mu=n_t_mu,
                n_scatter_mu=n_scatter_mu,
                n_scatter_phi=n_scatter_phi
            )

            c_sd_0 = capture_rate_total(
                earth_data=earth_data,
                DM_mass=m,
                sigma_SD_p=sigma_SD_p,
                scattering_type="SD",
                cross_section_type=cross_type,
                rho_chi=rho_chi,
                u_max=u_max,
                n_u=n_u,
                v0=v0,
                include_thermal_targets=False,
                n_t_speed=n_t_speed,
                n_t_mu=n_t_mu,
                n_scatter_mu=n_scatter_mu,
                n_scatter_phi=n_scatter_phi
            )

            # Electron
            c_e_T = capture_rate_total(
                earth_data=earth_data,
                DM_mass=m,
                sigma_electron=sigma_electron,
                scattering_type="electron",
                cross_section_type=cross_type,
                rho_chi=rho_chi,
                u_max=u_max,
                n_u=n_u,
                v0=v0,
                include_thermal_targets=True,
                n_t_speed=n_t_speed,
                n_t_mu=n_t_mu,
                n_scatter_mu=n_scatter_mu,
                n_scatter_phi=n_scatter_phi
            )

            c_e_0 = capture_rate_total(
                earth_data=earth_data,
                DM_mass=m,
                sigma_electron=sigma_electron,
                scattering_type="electron",
                cross_section_type=cross_type,
                rho_chi=rho_chi,
                u_max=u_max,
                n_u=n_u,
                v0=v0,
                include_thermal_targets=False,
                n_t_speed=n_t_speed,
                n_t_mu=n_t_mu,
                n_scatter_mu=n_scatter_mu,
                n_scatter_phi=n_scatter_phi
            )

            C_SI_T.append(c_si_T)
            C_SI_0.append(c_si_0)
            C_SD_T.append(c_sd_T)
            C_SD_0.append(c_sd_0)
            C_e_T.append(c_e_T)
            C_e_0.append(c_e_0)

        C_SI_T = np.array(C_SI_T)
        C_SI_0 = np.array(C_SI_0)
        C_SD_T = np.array(C_SD_T)
        C_SD_0 = np.array(C_SD_0)
        C_e_T = np.array(C_e_T)
        C_e_0 = np.array(C_e_0)

        # Left panel: finite-T rates
        axL = axes[row, 0]
        axL.loglog(DM_masses, C_SI_T, "b-.", lw=2.5, label="SI, T(r)≠0")
        axL.loglog(DM_masses, C_SD_T, "g--", lw=2.5, label="SD, T(r)≠0")

        axL.loglog(DM_masses, C_SI_0, color="blue", alpha=0.25, ls=":", lw=2.0, label="SI, T=0")
        axL.loglog(DM_masses, C_SD_0, color="green", alpha=0.25, ls=":", lw=2.0, label="SD, T=0")
        axL.loglog(DM_masses, C_e_0, color="red", alpha=0.25, ls=":", lw=2.0, label="Electron, T=0")

        axL.set_xlabel(r"$m_\chi$ [GeV]", fontsize=12)
        axL.set_ylabel(r"$C$ [s$^{-1}$]", fontsize=12)
        axL.set_title(f"{cross_type}, sigma = {sigma_SI_p:.0e} cm$^2$", fontsize=12)
        axL.grid(True, alpha=0.3)
        axL.legend(fontsize=9, loc="best")

        # Right panel: thermal enhancement/suppression
        axR = axes[row, 1]
        ratio_SI = np.where(C_SI_0 > 0, C_SI_T / C_SI_0, np.nan)
        ratio_SD = np.where(C_SD_0 > 0, C_SD_T / C_SD_0, np.nan)
        ratio_e = np.where(C_e_0 > 0, C_e_T / C_e_0, np.nan)

        axR.semilogx(DM_masses, ratio_SI, "b-.", lw=2.0, label="SI")
        axR.semilogx(DM_masses, ratio_SD, "g--", lw=2.0, label="SD")
        axR.semilogx(DM_masses, ratio_e, "r-", lw=2.0, label="Electron")
        axR.axhline(1.0, color="k", linestyle=":", lw=1.2, label=r"$T(r)=0$ limit")

        axR.set_xlabel(r"$m_\chi$ [GeV]", fontsize=12)
        axR.set_ylabel(r"$C(T) / C(T=0)$", fontsize=12)
        axR.set_title("Thermal correction factor", fontsize=12)
        axR.grid(True, alpha=0.3)
        axR.legend(fontsize=10, loc="best")

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"\nPlot saved as: {save_path}")
    plt.show()

# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    print("Loading Earth data...")
    earth_data = load_earth_composition("data/earth_prem.dat")
    print(f"Loaded {len(earth_data['radius'])} layers, {len(earth_data['abundances'])} elements")

    DM_masses = np.logspace(-2, 3, 10)

    plot_capture_rates_complete_thermal(
        earth_data=earth_data,
        DM_masses=DM_masses,
        sigma_values=[
            (1e-40, 1e-40, 1e-40),
            (1e-42, 1e-42, 1e-42),
            (1e-42, 1e-42, 1e-42),
        ],
        cross_section_types=[
            "constant",
            "v2_dependent",
            "q2_dependent",
        ],
        save_path="complete_capture_rates_thermal.png",
        u_max=800.0,
        n_u=60,
        n_t_speed=4,
        n_t_mu=4,
        n_scatter_mu=6,
        n_scatter_phi=8
    )