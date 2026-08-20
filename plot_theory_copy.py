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
import csv



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
    # ============================================================
    # Verified isotopes (with literature values) - enabled=True
    # ============================================================
    "H": {
        "element_name": "Hydrogen",
        "isotopes": [
            {
                "A": 1,
                "name": "Hydrogen-1",
                "abundance": 1.0,
                "J": 0.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "verified",
                "enabled": True,
                "source": "proton benchmark",
                "notes": ""
            },
            {
                "A": 2,
                "name": "Deuterium",
                "abundance": 0.00015,
                "J": 1.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Deuterium spin=1, spin structure factor ≈ 0 (simplified)"
            },
        ]
    },
    "He": {
        "element_name": "Helium",
        "isotopes": [
            {
                "A": 4,
                "name": "Helium-4",
                "abundance": 1.0,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Li": {
        "element_name": "Lithium",
        "isotopes": [
            {
                "A": 7,
                "name": "Lithium-7",
                "abundance": 0.925,
                "J": 1.5,
                "S_p": 0.5,
                "S_n": 0.5,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: Sp=Sn=0.5, awaiting literature values"
            },
        ]
    },
    "Be": {
        "element_name": "Beryllium",
        "isotopes": [
            {
                "A": 9,
                "name": "Beryllium-9",
                "abundance": 1.0,
                "J": 1.5,
                "S_p": 0.5,
                "S_n": 0.5,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "B": {
        "element_name": "Boron",
        "isotopes": [
            {
                "A": 11,
                "name": "Boron-11",
                "abundance": 0.801,
                "J": 1.5,
                "S_p": 0.5,
                "S_n": 0.5,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "C": {
        "element_name": "Carbon",
        "isotopes": [
            {
                "A": 12,
                "name": "Carbon-12",
                "abundance": 1.0,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "N": {
        "element_name": "Nitrogen",
        "isotopes": [
            {
                "A": 14,
                "name": "Nitrogen-14",
                "abundance": 0.996,
                "J": 1.0,
                "S_p": 0.0,
                "S_n": 0.5,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "O": {
        "element_name": "Oxygen",
        "isotopes": [
            {
                "A": 16,
                "name": "Oxygen-16",
                "abundance": 1.0,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
            {
                "A": 17,
                "name": "Oxygen-17",
                "abundance": 0.00038,
                "J": 2.5,
                "S_p": 0.0,
                "S_n": 0.5,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "F": {
        "element_name": "Fluorine",
        "isotopes": [
            {
                "A": 19,
                "name": "Fluorine-19",
                "abundance": 1.0,
                "J": 0.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Ne": {
        "element_name": "Neon",
        "isotopes": [
            {
                "A": 20,
                "name": "Neon-20",
                "abundance": 0.904,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Na": {
        "element_name": "Sodium",
        "isotopes": [
            {
                "A": 23,
                "name": "Sodium-23",
                "abundance": 1.0,
                "J": 1.5,
                "S_p": 0.224,
                "S_n": 0.024,
                "status": "verified",
                "enabled": True,
                "source": "Klos et al. 2013 / Hu et al. 2022",
                "notes": "23Na, odd-Z, S_p dominates"
            },
        ]
    },
    "Mg": {
        "element_name": "Magnesium",
        "isotopes": [
            {
                "A": 24,
                "name": "Magnesium-24",
                "abundance": 0.79,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Al": {
        "element_name": "Aluminium",
        "isotopes": [
            {
                "A": 27,
                "name": "Aluminium-27",
                "abundance": 1.0,
                "J": 2.5,
                "S_p": 0.326,
                "S_n": 0.038,
                "status": "verified",
                "enabled": True,
                "source": "Klos et al. 2013 / Hu et al. 2022",
                "notes": "27Al, odd-Z, S_p dominates"
            },
        ]
    },
    "Si": {
        "element_name": "Silicon",
        "isotopes": [
            {
                "A": 28,
                "name": "Silicon-28",
                "abundance": 0.922,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
            {
                "A": 29,
                "name": "Silicon-29",
                "abundance": 0.047,
                "J": 0.5,
                "S_p": 0.016,
                "S_n": 0.156,
                "status": "verified",
                "enabled": True,
                "source": "Klos et al. 2013 / Hu et al. 2022",
                "notes": "29Si, odd-N, S_n dominates"
            },
        ]
    },
    "P": {
        "element_name": "Phosphorus",
        "isotopes": [
            {
                "A": 31,
                "name": "Phosphorus-31",
                "abundance": 1.0,
                "J": 0.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "S": {
        "element_name": "Sulfur",
        "isotopes": [
            {
                "A": 32,
                "name": "Sulfur-32",
                "abundance": 0.95,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Cl": {
        "element_name": "Chlorine",
        "isotopes": [
            {
                "A": 35,
                "name": "Chlorine-35",
                "abundance": 0.758,
                "J": 1.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Ar": {
        "element_name": "Argon",
        "isotopes": [
            {
                "A": 40,
                "name": "Argon-40",
                "abundance": 0.996,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "K": {
        "element_name": "Potassium",
        "isotopes": [
            {
                "A": 39,
                "name": "Potassium-39",
                "abundance": 0.932,
                "J": 1.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Ca": {
        "element_name": "Calcium",
        "isotopes": [
            {
                "A": 40,
                "name": "Calcium-40",
                "abundance": 0.969,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Sc": {
        "element_name": "Scandium",
        "isotopes": [
            {
                "A": 45,
                "name": "Scandium-45",
                "abundance": 1.0,
                "J": 3.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Ti": {
        "element_name": "Titanium",
        "isotopes": [
            {
                "A": 48,
                "name": "Titanium-48",
                "abundance": 0.737,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "V": {
        "element_name": "Vanadium",
        "isotopes": [
            {
                "A": 51,
                "name": "Vanadium-51",
                "abundance": 0.997,
                "J": 3.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Cr": {
        "element_name": "Chromium",
        "isotopes": [
            {
                "A": 52,
                "name": "Chromium-52",
                "abundance": 0.838,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Mn": {
        "element_name": "Manganese",
        "isotopes": [
            {
                "A": 55,
                "name": "Manganese-55",
                "abundance": 1.0,
                "J": 2.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Fe": {
        "element_name": "Iron",
        "isotopes": [
            {
                "A": 56,
                "name": "Iron-56",
                "abundance": 0.917,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
            {
                "A": 57,
                "name": "Iron-57",
                "abundance": 0.021,
                "J": 0.5,
                "S_p": 0.0,
                "S_n": 0.5,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Co": {
        "element_name": "Cobalt",
        "isotopes": [
            {
                "A": 59,
                "name": "Cobalt-59",
                "abundance": 1.0,
                "J": 3.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Ni": {
        "element_name": "Nickel",
        "isotopes": [
            {
                "A": 58,
                "name": "Nickel-58",
                "abundance": 0.682,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Cu": {
        "element_name": "Copper",
        "isotopes": [
            {
                "A": 63,
                "name": "Copper-63",
                "abundance": 0.691,
                "J": 1.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Zn": {
        "element_name": "Zinc",
        "isotopes": [
            {
                "A": 64,
                "name": "Zinc-64",
                "abundance": 0.492,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Ga": {
        "element_name": "Gallium",
        "isotopes": [
            {
                "A": 69,
                "name": "Gallium-69",
                "abundance": 0.601,
                "J": 1.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Ge": {
        "element_name": "Germanium",
        "isotopes": [
            {
                "A": 70,
                "name": "Germanium-70",
                "abundance": 0.205,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "As": {
        "element_name": "Arsenic",
        "isotopes": [
            {
                "A": 75,
                "name": "Arsenic-75",
                "abundance": 1.0,
                "J": 1.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Se": {
        "element_name": "Selenium",
        "isotopes": [
            {
                "A": 78,
                "name": "Selenium-78",
                "abundance": 0.237,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Br": {
        "element_name": "Bromine",
        "isotopes": [
            {
                "A": 79,
                "name": "Bromine-79",
                "abundance": 0.507,
                "J": 1.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Kr": {
        "element_name": "Krypton",
        "isotopes": [
            {
                "A": 80,
                "name": "Krypton-80",
                "abundance": 0.023,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Rb": {
        "element_name": "Rubidium",
        "isotopes": [
            {
                "A": 85,
                "name": "Rubidium-85",
                "abundance": 0.722,
                "J": 2.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Sr": {
        "element_name": "Strontium",
        "isotopes": [
            {
                "A": 88,
                "name": "Strontium-88",
                "abundance": 0.826,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Y": {
        "element_name": "Yttrium",
        "isotopes": [
            {
                "A": 89,
                "name": "Yttrium-89",
                "abundance": 1.0,
                "J": 0.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Zr": {
        "element_name": "Zirconium",
        "isotopes": [
            {
                "A": 90,
                "name": "Zirconium-90",
                "abundance": 0.514,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Nb": {
        "element_name": "Niobium",
        "isotopes": [
            {
                "A": 93,
                "name": "Niobium-93",
                "abundance": 1.0,
                "J": 4.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Mo": {
        "element_name": "Molybdenum",
        "isotopes": [
            {
                "A": 92,
                "name": "Molybdenum-92",
                "abundance": 0.148,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Ru": {
        "element_name": "Ruthenium",
        "isotopes": [
            {
                "A": 100,
                "name": "Ruthenium-100",
                "abundance": 0.126,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Rh": {
        "element_name": "Rhodium",
        "isotopes": [
            {
                "A": 103,
                "name": "Rhodium-103",
                "abundance": 1.0,
                "J": 0.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Pd": {
        "element_name": "Palladium",
        "isotopes": [
            {
                "A": 106,
                "name": "Palladium-106",
                "abundance": 0.273,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Ag": {
        "element_name": "Silver",
        "isotopes": [
            {
                "A": 107,
                "name": "Silver-107",
                "abundance": 0.518,
                "J": 0.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Cd": {
        "element_name": "Cadmium",
        "isotopes": [
            {
                "A": 110,
                "name": "Cadmium-110",
                "abundance": 0.125,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "In": {
        "element_name": "Indium",
        "isotopes": [
            {
                "A": 115,
                "name": "Indium-115",
                "abundance": 0.957,
                "J": 4.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Sn": {
        "element_name": "Tin",
        "isotopes": [
            {
                "A": 116,
                "name": "Tin-116",
                "abundance": 0.146,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Sb": {
        "element_name": "Antimony",
        "isotopes": [
            {
                "A": 121,
                "name": "Antimony-121",
                "abundance": 0.572,
                "J": 2.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Te": {
        "element_name": "Tellurium",
        "isotopes": [
            {
                "A": 120,
                "name": "Tellurium-120",
                "abundance": 0.0009,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "I": {
        "element_name": "Iodine",
        "isotopes": [
            {
                "A": 127,
                "name": "Iodine-127",
                "abundance": 1.0,
                "J": 2.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Xe": {
        "element_name": "Xenon",
        "isotopes": [
            {
                "A": 128,
                "name": "Xenon-128",
                "abundance": 0.001,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Cs": {
        "element_name": "Cesium",
        "isotopes": [
            {
                "A": 133,
                "name": "Cesium-133",
                "abundance": 1.0,
                "J": 3.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Ba": {
        "element_name": "Barium",
        "isotopes": [
            {
                "A": 134,
                "name": "Barium-134",
                "abundance": 0.024,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "La": {
        "element_name": "Lanthanum",
        "isotopes": [
            {
                "A": 139,
                "name": "Lanthanum-139",
                "abundance": 0.999,
                "J": 3.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Ce": {
        "element_name": "Cerium",
        "isotopes": [
            {
                "A": 140,
                "name": "Cerium-140",
                "abundance": 0.885,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Pr": {
        "element_name": "Praseodymium",
        "isotopes": [
            {
                "A": 141,
                "name": "Praseodymium-141",
                "abundance": 1.0,
                "J": 2.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Nd": {
        "element_name": "Neodymium",
        "isotopes": [
            {
                "A": 142,
                "name": "Neodymium-142",
                "abundance": 0.272,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Sm": {
        "element_name": "Samarium",
        "isotopes": [
            {
                "A": 148,
                "name": "Samarium-148",
                "abundance": 0.114,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Eu": {
        "element_name": "Europium",
        "isotopes": [
            {
                "A": 151,
                "name": "Europium-151",
                "abundance": 0.478,
                "J": 2.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Gd": {
        "element_name": "Gadolinium",
        "isotopes": [
            {
                "A": 152,
                "name": "Gadolinium-152",
                "abundance": 0.002,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Tb": {
        "element_name": "Terbium",
        "isotopes": [
            {
                "A": 159,
                "name": "Terbium-159",
                "abundance": 1.0,
                "J": 1.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Dy": {
        "element_name": "Dysprosium",
        "isotopes": [
            {
                "A": 158,
                "name": "Dysprosium-158",
                "abundance": 0.0001,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Ho": {
        "element_name": "Holmium",
        "isotopes": [
            {
                "A": 165,
                "name": "Holmium-165",
                "abundance": 1.0,
                "J": 3.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Er": {
        "element_name": "Erbium",
        "isotopes": [
            {
                "A": 166,
                "name": "Erbium-166",
                "abundance": 0.336,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Tm": {
        "element_name": "Thulium",
        "isotopes": [
            {
                "A": 169,
                "name": "Thulium-169",
                "abundance": 1.0,
                "J": 0.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Yb": {
        "element_name": "Ytterbium",
        "isotopes": [
            {
                "A": 170,
                "name": "Ytterbium-170",
                "abundance": 0.001,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Lu": {
        "element_name": "Lutetium",
        "isotopes": [
            {
                "A": 175,
                "name": "Lutetium-175",
                "abundance": 0.974,
                "J": 3.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Hf": {
        "element_name": "Hafnium",
        "isotopes": [
            {
                "A": 176,
                "name": "Hafnium-176",
                "abundance": 0.052,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Ta": {
        "element_name": "Tantalum",
        "isotopes": [
            {
                "A": 181,
                "name": "Tantalum-181",
                "abundance": 0.999,
                "J": 3.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "W": {
        "element_name": "Tungsten",
        "isotopes": [
            {
                "A": 182,
                "name": "Tungsten-182",
                "abundance": 0.265,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Re": {
        "element_name": "Rhenium",
        "isotopes": [
            {
                "A": 185,
                "name": "Rhenium-185",
                "abundance": 0.374,
                "J": 2.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Os": {
        "element_name": "Osmium",
        "isotopes": [
            {
                "A": 186,
                "name": "Osmium-186",
                "abundance": 0.015,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Ir": {
        "element_name": "Iridium",
        "isotopes": [
            {
                "A": 191,
                "name": "Iridium-191",
                "abundance": 0.373,
                "J": 1.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Pt": {
        "element_name": "Platinum",
        "isotopes": [
            {
                "A": 192,
                "name": "Platinum-192",
                "abundance": 0.0078,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Au": {
        "element_name": "Gold",
        "isotopes": [
            {
                "A": 197,
                "name": "Gold-197",
                "abundance": 1.0,
                "J": 1.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Hg": {
        "element_name": "Mercury",
        "isotopes": [
            {
                "A": 198,
                "name": "Mercury-198",
                "abundance": 0.001,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Tl": {
        "element_name": "Thallium",
        "isotopes": [
            {
                "A": 203,
                "name": "Thallium-203",
                "abundance": 0.705,
                "J": 0.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
    },
    "Pb": {
        "element_name": "Lead",
        "isotopes": [
            {
                "A": 204,
                "name": "Lead-204",
                "abundance": 0.014,
                "J": 0.0,
                "S_p": 0.0,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": False,
                "source": "",
                "notes": "J=0, no SD scattering"
            },
        ]
    },
    "Bi": {
        "element_name": "Bismuth",
        "isotopes": [
            {
                "A": 209,
                "name": "Bismuth-209",
                "abundance": 1.0,
                "J": 4.5,
                "S_p": 0.5,
                "S_n": 0.0,
                "status": "placeholder",
                "enabled": True,
                "source": "",
                "notes": "Placeholder: awaiting literature values"
            },
        ]
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


# ============================================================
# Basic helpers
# ============================================================
def print_sd_db_status_summary(sd_mode="include_placeholders", show_all=True):
    """
    Print a summary of the current SD_SPIN_DB status.

    Parameters
    ----------
    sd_mode : str
        "verified_only" or "include_placeholders"
        This affects the "active under current mode" summary.
    show_all : bool
        If True, print every isotope entry.
        If False, only print compact element-level summaries.
    """
    if sd_mode not in ("verified_only", "include_placeholders"):
        raise ValueError("sd_mode must be 'verified_only' or 'include_placeholders'")

    print("\n" + "=" * 90)
    print(f"SD_SPIN_DB status summary  (sd_mode = {sd_mode})")
    print("=" * 90)

    total_elements = 0
    total_isotopes = 0
    total_enabled = 0
    total_verified = 0
    total_placeholder = 0
    total_active_under_mode = 0

    active_labels_under_mode = []

    for elem in sorted(SD_SPIN_DB.keys()):
        elem_block = SD_SPIN_DB[elem]
        isotopes = elem_block.get("isotopes", [])

        total_elements += 1
        total_isotopes += len(isotopes)

        n_enabled = 0
        n_verified = 0
        n_placeholder = 0
        n_active_mode = 0

        element_active_labels = []

        for iso in isotopes:
            A = iso.get("A", None)
            label = f"{elem}{A}" if A is not None and A != int(round(ELEMENT_DB[elem]["A"])) else elem

            enabled = iso.get("enabled", True)
            status = iso.get("status", "verified")
            J = iso.get("J", None)
            S_p = iso.get("S_p", None)
            S_n = iso.get("S_n", None)

            if enabled:
                n_enabled += 1
                total_enabled += 1

            if status == "verified":
                n_verified += 1
                total_verified += 1
            elif status == "placeholder":
                n_placeholder += 1
                total_placeholder += 1

            # same logic as get_sd_isotopes_for_element(...)
            active_now = True
            if not enabled:
                active_now = False
            if sd_mode == "verified_only" and status != "verified":
                active_now = False
            if A is None or J is None or S_p is None or S_n is None:
                active_now = False
            elif J <= 0.0:
                active_now = False

            if active_now:
                n_active_mode += 1
                total_active_under_mode += 1
                element_active_labels.append(label)
                active_labels_under_mode.append(label)

        print(
            f"{elem:>3s} : "
            f"isotopes={len(isotopes):2d}, "
            f"enabled={n_enabled:2d}, "
            f"verified={n_verified:2d}, "
            f"placeholder={n_placeholder:2d}, "
            f"active_under_mode={n_active_mode:2d}"
        )

        if element_active_labels:
            print(f"      active labels -> {element_active_labels}")

        if show_all and isotopes:
            for iso in isotopes:
                A = iso.get("A", None)
                label = f"{elem}{A}" if A is not None and A != int(round(ELEMENT_DB[elem]["A"])) else elem

                print(
                    "      "
                    f"label={label:<8s} "
                    f"A={str(iso.get('A', None)):<4s} "
                    f"abundance={iso.get('abundance', None)!s:<10s} "
                    f"J={iso.get('J', None)!s:<8s} "
                    f"S_p={iso.get('S_p', None)!s:<10s} "
                    f"S_n={iso.get('S_n', None)!s:<10s} "
                    f"status={iso.get('status', 'verified'):<12s} "
                    f"enabled={iso.get('enabled', True)}"
                )

    print("-" * 90)
    print(f"Total elements in SD_SPIN_DB           : {total_elements}")
    print(f"Total isotope entries                  : {total_isotopes}")
    print(f"Total enabled entries                  : {total_enabled}")
    print(f"Total verified entries                 : {total_verified}")
    print(f"Total placeholder entries              : {total_placeholder}")
    print(f"Total active entries under '{sd_mode}' : {total_active_under_mode}")
    print("-" * 90)

    unique_active_labels = sorted(set(active_labels_under_mode))
    print("Active labels under current mode:")
    print(unique_active_labels)
    print(f"Count = {len(unique_active_labels)}")
    print("=" * 90)

def get_sd_isotopes_for_element(elem, sd_mode="include_placeholders"):
    """
    Return enabled SD isotopes for one base element from the new SD_SPIN_DB format.

    Parameters
    ----------
    elem : str
        Base element symbol, e.g. 'Si', 'Fe', 'Al'
    sd_mode : str
        'verified_only' or 'include_placeholders'

    Returns
    -------
    list of dict
        Each dict contains one usable SD isotope entry.
    """
    if elem not in SD_SPIN_DB:
        return []

    if sd_mode not in ("verified_only", "include_placeholders"):
        raise ValueError("sd_mode must be 'verified_only' or 'include_placeholders'")

    elem_block = SD_SPIN_DB[elem]
    isotopes = elem_block.get("isotopes", [])

    out = []

    for iso in isotopes:
        if not iso.get("enabled", True):
            continue

        status = iso.get("status", "verified")
        if sd_mode == "verified_only" and status != "verified":
            continue

        A = iso.get("A", None)
        J = iso.get("J", None)
        S_p = iso.get("S_p", None)
        S_n = iso.get("S_n", None)

        if A is None or J is None or S_p is None or S_n is None:
            continue

        if J <= 0.0:
            continue

        A = int(A)
        label = f"{elem}{A}" if A != int(round(ELEMENT_DB[elem]["A"])) else elem

        out.append({
            "label": label,                      # e.g. 'Si29', 'Fe57', or 'Al'
            "base_elem": elem,                  # e.g. 'Si'
            "A": A,
            "name": iso.get("name", f"{elem}-{A}"),
            "abundance": float(iso.get("abundance", 1.0)),
            "J": float(J),
            "S_p": float(S_p),
            "S_n": float(S_n),
            "m_t": A * M_U_GEV,
            "status": status,
            "source": iso.get("source", ""),
            "notes": iso.get("notes", ""),
        })

    return out

def validate_sd_spin_db():
    """
    Validate SD_SPIN_DB and print entries with missing fields.
    """
    problems = []

    for elem, elem_data in SD_SPIN_DB.items():
        isotopes = elem_data.get("isotopes", [])
        if not isotopes:
            problems.append(f"{elem}: no isotope entries")
            continue

        for iso in isotopes:
            A = iso.get("A", None)
            name = iso.get("name", f"{elem}-unknown")

            for field in ["A", "abundance", "J", "S_p", "S_n"]:
                if field not in iso:
                    problems.append(f"{elem} {name}: missing field '{field}'")

            if iso.get("enabled", True):
                if iso.get("J", None) is None:
                    problems.append(f"{elem} {name}: enabled but J is None")
                if iso.get("S_p", None) is None:
                    problems.append(f"{elem} {name}: enabled but S_p is None")
                if iso.get("S_n", None) is None:
                    problems.append(f"{elem} {name}: enabled but S_n is None")

    if problems:
        print("SD_SPIN_DB validation problems:")
        for p in problems:
            print(" -", p)
    else:
        print("SD_SPIN_DB validation passed.")

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
def load_earth_composition(
    filepath="data/earth_prem.dat",
    sd_mode="include_placeholders",
    min_mass_fraction=1e-10
):
    """
    Read earth_prem.dat and return:
    {
        "radius": np.array([...])          # km
        "density": np.array([...])         # g/cm^3
        "temperature": np.array([...])     # K
        "abundances": {elem: np.array([...])}  # ppm by mass
    }

    Parameters
    ----------
    sd_mode : str
        'verified_only' or 'include_placeholders'
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

        for j, elem in enumerate(elements):
            abundances[elem].append(float(cols[4 + j]))

    earth_data = {
        "radius": np.array(radii, dtype=float),
        "density": np.array(densities, dtype=float),
        "temperature": np.array(temperatures, dtype=float),
        "abundances": {k: np.array(v, dtype=float) for k, v in abundances.items()}
    }

    precompute_earth_shells(earth_data)
    build_active_shell_data(
        earth_data,
        min_mass_fraction=min_mass_fraction,
        sd_mode=sd_mode
    )
    return earth_data

def get_base_element_from_sd_key(sd_key):
    """
    Convert SD key like 'Si29', 'Fe57', 'O17', 'H2' -> base element
    'Si', 'Fe', 'O', 'H'. Exact element keys like 'Al' are returned unchanged.
    """
    import re

    m = re.fullmatch(r"([A-Z][a-z]?)(\d+)?", sd_key)
    if m is None:
        raise ValueError(f"Cannot parse SD key: {sd_key}")

    base_elem = m.group(1)
    if base_elem not in ELEMENT_DB:
        raise KeyError(f"Base element '{base_elem}' from SD key '{sd_key}' not found in ELEMENT_DB")

    return base_elem


def get_sd_spinful_isotopes_for_element(elem):
    """
    Return all spinful SD isotope entries associated with a base element.

    Examples
    --------
    elem = 'Si' -> ['Si29']
    elem = 'Fe' -> ['Fe57']
    elem = 'O'  -> ['O17']
    elem = 'Al' -> ['Al']
    """
    if elem not in ELEMENT_DB:
        return []

    out = []

    for sd_key, spin_data in SD_SPIN_DB.items():
        # exact key: 'Al'
        # isotope-style key: 'Si29', 'Fe57', ...
        if sd_key == elem:
            pass
        elif sd_key.startswith(elem) and sd_key[len(elem):].isdigit():
            pass
        else:
            continue

        J = spin_data.get("J", 0.0)
        if J <= 0.0:
            continue

        iso_ab = spin_data.get("isotope_abundance", spin_data.get("abundance", 1.0))
        A_sd = int(spin_data.get("isotope", round(ELEMENT_DB[elem]["A"])))

        out.append({
            "sd_key": sd_key,      # e.g. 'Si29', 'Fe57', 'Al'
            "J": J,
            "iso_ab": iso_ab,
            "A_sd": A_sd,
            "m_t_sd": A_sd * M_U_GEV,
        })

    return out

def build_active_shell_data(earth_data, min_mass_fraction=1e-10, sd_mode="include_placeholders"):
    """
    Precompute active elements/isotopes in each shell.

    SI:
        stored by element

    SD:
        for each base element in Earth composition, expand to all enabled
        spinful isotopes listed in SD_SPIN_DB.

    Electron:
        stored as before.

    Notes
    -----
    For SD isotopes, use isotope-specific effective number density:
        n_i_eff = rho * X_elem * f_iso / (A_iso * m_u)
    """
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

        rho_g_cm3 = earth_data["density"][i]

        for elem in earth_data["abundances"]:
            if elem not in ELEMENT_DB:
                continue

            mass_fraction = earth_data["abundances"][elem][i] / 1e6
            if mass_fraction < min_mass_fraction:
                continue

            # ---------------- SI ----------------
            n_i_elem = number_density_at_radius(earth_data, i, elem)
            if n_i_elem > 0.0:
                A_eff = int(round(ELEMENT_DB[elem]["A"]))
                m_t = A_eff * M_U_GEV

                shell_info["SI"].append({
                    "elem": elem,
                    "n_i": n_i_elem,
                    "m_t": m_t,
                })

            # ---------------- SD ----------------
            sd_isotopes = get_sd_isotopes_for_element(elem, sd_mode=sd_mode)

            for iso in sd_isotopes:
                A_sd = iso["A"]
                iso_ab = iso["abundance"]

                n_i_eff = rho_g_cm3 * mass_fraction * iso_ab / (A_sd * M_U_G)

                if n_i_eff <= 0.0:
                    continue

                shell_info["SD"].append({
                    "elem": iso["label"],          # e.g. 'Si29', 'Fe57', or 'Al'
                    "base_elem": iso["base_elem"], # e.g. 'Si'
                    "A": A_sd,
                    "iso_name": iso["name"],
                    "n_i_eff": n_i_eff,
                    "m_t": iso["m_t"],
                    "J": iso["J"],
                    "S_p": iso["S_p"],
                    "S_n": iso["S_n"],
                    "abundance": iso_ab,
                    "status": iso["status"],
                    "source": iso["source"],
                    "notes": iso["notes"],
                })

        active_shells.append(shell_info)

    earth_data["active_shells"] = active_shells
    earth_data["sd_mode"] = sd_mode
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


def sigma_electron_q(DM_mass, g_km_s, q_GeV, sigma_electron, cross_section_type="constant"):
    sigma_e = get_cross_section(sigma_electron, g_km_s, q_GeV, cross_section_type)
    return sigma_e, M_E_GEV

def sigma_nucleus_SD_from_info(
    DM_mass,
    sd_info,
    w_km_s,
    sigma_SD_p,
    sigma_SD_n=0.0,
    cross_section_type="constant"
):
    """
    SD nuclear cross section using one isotope record from shell_info['SD'].
    """
    J = sd_info["J"]
    if J <= 0.0:
        return 0.0, None, 0.0

    A_eff = int(sd_info["A"])
    m_A = float(sd_info["m_t"])
    S_p = float(sd_info["S_p"])
    S_n = float(sd_info["S_n"])

    mu_A = reduced_mass(DM_mass, m_A)
    mu_p = reduced_mass(DM_mass, M_P_GEV)

    q_typ = 2.0 * mu_A * (w_km_s / C_LIGHT_KM_S)  # GeV
    F_q = nuclear_form_factor_for_nucleus(A_eff, q_typ, "SD")
    sigma_eff_p = get_cross_section(sigma_SD_p, w_km_s, q_typ, cross_section_type)

    spin_factor = 4.0 * (J + 1.0) / (3.0 * J)
    spin_contrib = (S_p + S_n)   # keep same simplified convention as your current code

    sigma_A = (mu_A / mu_p)**2 * spin_factor * (spin_contrib**2) * sigma_eff_p * F_q**2
    return sigma_A, m_A, sd_info.get("abundance", 1.0)

def sigma_nucleus_SD_q_from_info(
    DM_mass,
    sd_info,
    g_km_s,
    q_GeV,
    sigma_SD_p,
    sigma_SD_n=0.0,
    cross_section_type="constant"
):
    """
    SD nuclear cross section at given g and q using one isotope record from shell_info['SD'].
    """
    J = sd_info["J"]
    if J <= 0.0:
        return 0.0, None, 0.0

    A_eff = int(sd_info["A"])
    m_A = float(sd_info["m_t"])
    S_p = float(sd_info["S_p"])
    S_n = float(sd_info["S_n"])

    mu_A = reduced_mass(DM_mass, m_A)
    mu_p = reduced_mass(DM_mass, M_P_GEV)

    F_q = nuclear_form_factor_for_nucleus(A_eff, q_GeV, scattering_type="SD")
    sigma_p_eff = get_cross_section(sigma_SD_p, g_km_s, q_GeV, cross_section_type)

    spin_factor = 4.0 * (J + 1.0) / (3.0 * J)
    spin_contrib = (S_p + S_n)   # keep same simplified convention as your current code

    sigma_A = (mu_A / mu_p)**2 * spin_factor * (spin_contrib**2) * sigma_p_eff * F_q**2
    return sigma_A, m_A, sd_info.get("abundance", 1.0)

_GL_CACHE = {}

def leggauss_cached(n):
    if n not in _GL_CACHE:
        _GL_CACHE[n] = np.polynomial.legendre.leggauss(n)
    return _GL_CACHE[n]

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
    w_km_s : float
        DM speed in lab frame [km/s]
    u_t_km_s : float
        target speed in lab frame [km/s]
    mu_in : float
        cos(angle between target velocity and DM velocity)
    m_chi : float
        DM mass [GeV]
    m_t : float
        target mass [GeV]
    vesc_km_s : float
        local escape velocity [km/s]
    sigma_eval : callable
        function sigma_eval(g_km_s, q_GeV) -> sigma [cm^2]
    n_scatter_mu : int
        quadrature order in cos(theta*)
    n_scatter_phi : int
        sampling number in phi*

    Returns
    -------
    float
        angle-averaged capture-weighted cross section [cm^2]
    """
    # Choose DM direction along +z
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

    # Center-of-mass velocity in lab
    V_vec = (m_chi * w_vec + m_t * u_vec) / Mtot
    V_km_s = np.linalg.norm(V_vec)

    # Quadrature in cos(theta*)
    mu_s_nodes, mu_s_weights = leggauss_cached(n_scatter_mu)

    # Uniform phi* sampling
    phi_nodes = 2.0 * np.pi * (np.arange(n_scatter_phi) + 0.5) / n_scatter_phi
    cosphi = np.cos(phi_nodes)[None, :]

    # Momentum transfer:
    # q = mu_red * (g/c) * sqrt(2(1-cos theta*))
    q_vals = mu_red * (g_km_s / C_LIGHT_KM_S) * np.sqrt(
        2.0 * np.maximum(0.0, 1.0 - mu_s_nodes)
    )
    sigma_vals = np.array([sigma_eval(g_km_s, q) for q in q_vals])

    # If COM velocity is ~0, final speed is angle-independent
    if V_km_s < 1e-12:
        vprime2 = (alpha * g_km_s)**2
        captured = 1.0 if vprime2 <= vesc_km_s**2 else 0.0
        return captured * 0.5 * np.sum(mu_s_weights * sigma_vals)

    # Angle between g_vec and V_vec
    cos_beta = np.dot(g_vec, V_vec) / (g_km_s * V_km_s)
    cos_beta = np.clip(cos_beta, -1.0, 1.0)
    sin_beta = np.sqrt(max(0.0, 1.0 - cos_beta**2))

    mu_s = mu_s_nodes[:, None]
    sin_s = np.sqrt(np.maximum(0.0, 1.0 - mu_s**2))

    # cos(psi) = cos(beta) cos(theta*) + sin(beta) sin(theta*) cos(phi*)
    cos_psi = cos_beta * mu_s + sin_beta * sin_s * cosphi

    # DM final speed squared in lab
    vprime2 = (
        V_km_s**2
        + (alpha * g_km_s)**2
        + 2.0 * alpha * g_km_s * V_km_s * cos_psi
    )

    # Fraction of phi* angles that lead to capture
    cap_frac_vs_mu_s = np.mean(vprime2 <= vesc_km_s**2, axis=1)

    # Average over cos(theta*) with factor 1/2
    return 0.5 * np.sum(mu_s_weights * sigma_vals * cap_frac_vs_mu_s) 

def thermal_average_gsigma_capture(
    w_km_s,
    vesc_km_s,
    T_K,
    m_chi,
    m_t,
    sigma_eval,
    include_thermal_targets=True,
    n_t_speed=16,
    n_t_mu=12,
    n_scatter_mu=8,
    n_scatter_phi=12,
    u_t_max_mult=8.0,
    thermal_grid_mode="linear"
):
    """
    Compute:
        < g * <sigma * Theta(capture)>_{Omega_*} >_{thermal target velocities}

    Returns
    -------
    avg_gsigma_capture : [cm^3 / s]

    Notes
    -----
    This version uses explicit grids + trapezoid integration for the target
    thermal speed and angle average. It is slower than the Gauss-Legendre
    version, but much more stable for light DM and threshold-sensitive capture.
    """
    # ------------------------------------------------------------
    # T = 0 limit
    # ------------------------------------------------------------
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

    if v_th < 1e-12:
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

    # Safety floors for stability
    n_t_speed = max(int(n_t_speed), 8)
    n_t_mu = max(int(n_t_mu), 6)

    # ------------------------------------------------------------
    # Build thermal-speed grid in x = u_t / v_th
    # ------------------------------------------------------------
    x_max = float(u_t_max_mult)

    if thermal_grid_mode == "linear":
        x_grid = np.linspace(0.0, x_max, n_t_speed)

    elif thermal_grid_mode == "hybrid":
        # More resolution near x ~ 0..1, then linear to the tail
        n_low = max(4, n_t_speed // 3)
        n_high = max(4, n_t_speed - n_low + 1)

        x_low = np.linspace(0.0, 1.0, n_low, endpoint=False)
        x_high = np.linspace(1.0, x_max, n_high)

        x_grid = np.unique(np.concatenate([x_low, x_high]))

    else:
        raise ValueError("thermal_grid_mode must be 'linear' or 'hybrid'")

    u_grid = x_grid * v_th
    f_u = maxwell_speed_distribution(u_grid, v_th)

    # ------------------------------------------------------------
    # Build mu grid for target-direction average
    # ------------------------------------------------------------
    mu_grid = np.linspace(-1.0, 1.0, n_t_mu)

    avg_over_mu = np.zeros_like(u_grid)

    for i, uk in enumerate(u_grid):
        mu_vals = np.zeros_like(mu_grid)

        for j, mu_in in enumerate(mu_grid):
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
            mu_vals[j] = (g_km_s * 1e5) * sigma_cap

        # average over mu_in uniformly on [-1, 1]
        avg_over_mu[i] = 0.5 * np.trapezoid(mu_vals, mu_grid)

    total = np.trapezoid(f_u * avg_over_mu, u_grid)
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
    """
    omega = 0.0
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
                    n_scatter_phi=n_scatter_phi,
                    thermal_grid_mode="hybrid"
                )
                omega += n_i * avg_gsigma_cap

        elif scattering_type == "SD":
            for info in shell_info["SD"]:
                n_i_eff = info["n_i_eff"]
                m_t = info["m_t"]

                if n_i_eff <= 0.0:
                    continue

                def sigma_eval(g_km_s, q_GeV, info=info):
                    sigma_A, _, _ = sigma_nucleus_SD_q_from_info(
                        DM_mass=DM_mass,
                        sd_info=info,
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
                    n_t_speed=max(n_t_speed, 10),
                    n_t_mu=max(n_t_mu, 10),
                    n_scatter_mu=n_scatter_mu,
                    n_scatter_phi=n_scatter_phi,
                    u_t_max_mult=10.0,
                    thermal_grid_mode="hybrid"
                )
                omega += n_e * avg_gsigma_cap
        else:
            raise ValueError("scattering_type must be 'SI', 'SD', or 'electron'")

        return omega

    # ------------------------------------------------------------
    # Fallback path
    # ------------------------------------------------------------
    T_loc = earth_data["temperature"][radius_idx]
    vesc_loc = earth_data["v_esc_profile"][radius_idx]
    rho_g_cm3 = earth_data["density"][radius_idx]
    sd_mode = earth_data.get("sd_mode", "include_placeholders")

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
                n_scatter_phi=n_scatter_phi,
                thermal_grid_mode="hybrid"
            )
            omega += n_i * avg_gsigma_cap

    elif scattering_type == "SD":
        for elem in earth_data["abundances"]:
            if elem not in ELEMENT_DB:
                continue

            mass_fraction = earth_data["abundances"][elem][radius_idx] / 1e6
            if mass_fraction < min_mass_fraction:
                continue

            sd_isotopes = get_sd_isotopes_for_element(elem, sd_mode=sd_mode)
            for iso in sd_isotopes:
                A_sd = iso["A"]
                n_i_eff = rho_g_cm3 * mass_fraction * iso["abundance"] / (A_sd * M_U_G)

                if n_i_eff <= 0.0:
                    continue

                info = {
                    "elem": iso["label"],
                    "base_elem": iso["base_elem"],
                    "A": iso["A"],
                    "iso_name": iso["name"],
                    "n_i_eff": n_i_eff,
                    "m_t": iso["m_t"],
                    "J": iso["J"],
                    "S_p": iso["S_p"],
                    "S_n": iso["S_n"],
                    "abundance": iso["abundance"],
                    "status": iso["status"],
                    "source": iso["source"],
                    "notes": iso["notes"],
                }

                def sigma_eval(g_km_s, q_GeV, info=info):
                    sigma_A, _, _ = sigma_nucleus_SD_q_from_info(
                        DM_mass=DM_mass,
                        sd_info=info,
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
                    m_t=info["m_t"],
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
                n_t_speed=max(n_t_speed, 10),
                n_t_mu=max(n_t_mu, 10),
                n_scatter_mu=n_scatter_mu,
                n_scatter_phi=n_scatter_phi,
                u_t_max_mult=10.0,
                thermal_grid_mode="hybrid"
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
    min_mass_fraction=1e-10,
    u_min=1e-3,
    u_grid_mode="log",
    shell_step=1

):
    """
    Total capture rate:
        C = n_chi ∫ dV ∫ du f(u) (w/u) Ω^-(w,r,T)

    Parameters
    ----------
    u_grid_mode : str
        "linear", "log", or "hybrid"
        - linear : uniform in u
        - log    : logarithmic in u
        - hybrid : log at low u + linear at high u
    """
    n_chi = rho_chi / DM_mass

    # ------------------------------------------------------------
    # Build halo-speed grid
    # ------------------------------------------------------------
    u_min = max(float(u_min), 1e-6)

    if u_grid_mode == "linear":
        u_grid = np.linspace(u_min, u_max, n_u)

    elif u_grid_mode == "log":
        u_grid = np.geomspace(u_min, u_max, n_u)

    elif u_grid_mode == "hybrid":
        n_low = max(10, n_u // 2)
        u_split = min(50.0, 0.5 * u_max)

        u_low = np.geomspace(u_min, u_split, n_low, endpoint=False)
        n_high = max(2, n_u - len(u_low) + 1)
        u_high = np.linspace(u_split, u_max, n_high)

        u_grid = np.unique(np.concatenate([u_low, u_high]))

    else:
        raise ValueError("u_grid_mode must be 'linear', 'log', or 'hybrid'")

    f_u = shm_speed_distribution(u_grid, v0=v0, vesc=vesc_halo)

    total_C = 0.0

    for i in range(0, len(earth_data["radius"]), shell_step):
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

    print(f"[worker] start m = {m:.4g} GeV", flush=True)

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

    c_e_ref = capture_rate_total(
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

    print(f"[worker] done  m = {m:.4g} GeV", flush=True)

    return {
        "m": m,
        "C_SI_T": c_si_T,
        "C_SI_0": c_si_0,
        "C_SD_T": c_sd_T,
        "C_SD_0": c_sd_0,
        "C_e_T": c_e_ref,
        "C_e_0": c_e_ref,
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
    n_scatter_phi=12,
    max_workers=None
):
    """
    3x2 plot:
      Left  : finite-T nuclear capture rates + electron curve (electron thermal motion disabled)
      Right : C(T)/C(T=0) for SI and SD only
    """
    import os
    from concurrent.futures import ProcessPoolExecutor, as_completed

    fig, axes = plt.subplots(3, 2, figsize=(14, 18))

    if max_workers is None:
        max_workers = max(1, os.cpu_count() - 1)

    for row, (sigma_triplet, cross_type) in enumerate(zip(sigma_values, cross_section_types)):
        sigma_SI_p, sigma_SD_p, sigma_electron = sigma_triplet

        print(f"\n=== Row {row+1}: {cross_type} ===")
        print(f"Using {max_workers} worker processes")

        tasks = []
        for m in DM_masses:
            tasks.append((
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
            ))

        results = []
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(compute_one_mass_point, task) for task in tasks]

            for i, fut in enumerate(as_completed(futures), 1):
                res = fut.result()
                print(f"  done {i}/{len(futures)} : m_chi = {res['m']:.4g} GeV")
                results.append(res)

        results.sort(key=lambda x: x["m"])

        C_SI_T = np.array([r["C_SI_T"] for r in results])
        C_SI_0 = np.array([r["C_SI_0"] for r in results])
        C_SD_T = np.array([r["C_SD_T"] for r in results])
        C_SD_0 = np.array([r["C_SD_0"] for r in results])
        C_e_T  = np.array([r["C_e_T"]  for r in results])
        C_e_0  = np.array([r["C_e_0"]  for r in results])

        # Left panel
        axL = axes[row, 0]

        axL.loglog(DM_masses, C_SI_T, "b-.", lw=2.5, label="SI, T(r)≠0")
        #axL.loglog(DM_masses, C_SD_T, "g--", lw=2.5, label="SD, T(r)≠0")

        axL.loglog(DM_masses, C_SI_0, color="blue",  alpha=0.25, ls=":", lw=2.0, label="SI, T=0")
       # axL.loglog(DM_masses, C_SD_0, color="green", alpha=0.25, ls=":", lw=2.0, label="SD, T=0")

        axL.loglog(DM_masses, C_e_T, "r-", lw=2.5, label="Electron (T=0 used)")

        axL.set_xlabel(r"$m_\chi$ [GeV]", fontsize=12)
        axL.set_ylabel(r"$C$ [s$^{-1}$]", fontsize=12)
        axL.set_title(f"{cross_type}, sigma = {sigma_SI_p:.0e} cm$^2$", fontsize=12)
        axL.grid(True, alpha=0.3)
        axL.legend(fontsize=9, loc="best")

        # Right panel
        # Right panel
        axR = axes[row, 1]

        si_floor = max(1e-6 * np.nanmax(C_SI_0), 1e-300)
        sd_floor = max(1e-6 * np.nanmax(C_SD_0), 1e-300)

        ratio_SI = np.where(C_SI_0 > si_floor, C_SI_T / C_SI_0, np.nan)
        ratio_SD = np.where(C_SD_0 > sd_floor, C_SD_T / C_SD_0, np.nan)

        mask_mass = DM_masses >= 1e-3

        axR.semilogx(DM_masses[mask_mass], ratio_SI[mask_mass], "b-.", lw=2.0, label="SI")
        #axR.semilogx(DM_masses[mask_mass], ratio_SD[mask_mass], "g--", lw=2.0, label="SD")
        axR.axhline(1.0, color="k", linestyle=":", lw=1.2, label=r"$T(r)=0$ limit")

        axR.set_xlabel(r"$m_\chi$ [GeV]", fontsize=12)
        axR.set_ylabel(r"$C(T)/C(T=0)$", fontsize=12)
        axR.set_title("Thermal correction factor", fontsize=12)
        axR.set_yscale("log")
        axR.grid(True, alpha=0.3)
        axR.legend(fontsize=10, loc="best")

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"\nPlot saved as: {save_path}")
    plt.show()

def test_convergence_at_selected_masses(
    earth_data,
    masses,
    sigma_SI_p=1e-40,
    sigma_SD_p=1e-40,
    cross_section_type="constant",
    rho_chi=RHO_CHI_DEFAULT,
    v0=V0_DEFAULT
):
    """
    Simple convergence test for selected masses.
    Compare capture rates under several integration settings.
    """
    settings = [
        {"n_u": 12, "n_t_speed": 2, "n_t_mu": 2, "n_scatter_mu": 2, "n_scatter_phi": 3},
        {"n_u": 20, "n_t_speed": 3, "n_t_mu": 3, "n_scatter_mu": 3, "n_scatter_phi": 4},
        {"n_u": 30, "n_t_speed": 4, "n_t_mu": 4, "n_scatter_mu": 4, "n_scatter_phi": 6},
    ]

    print("\n" + "=" * 80)
    print("Convergence test")
    print("=" * 80)

    for m in masses:
        print(f"\n--- m_chi = {m:.4g} GeV ---")
        for s in settings:
            c_si_T = capture_rate_total(
                earth_data=earth_data,
                DM_mass=m,
                sigma_SI_p=sigma_SI_p,
                scattering_type="SI",
                cross_section_type=cross_section_type,
                rho_chi=rho_chi,
                u_max=800.0,
                n_u=s["n_u"],
                v0=v0,
                include_thermal_targets=True,
                n_t_speed=s["n_t_speed"],
                n_t_mu=s["n_t_mu"],
                n_scatter_mu=s["n_scatter_mu"],
                n_scatter_phi=s["n_scatter_phi"]
            )

            c_si_0 = capture_rate_total(
                earth_data=earth_data,
                DM_mass=m,
                sigma_SI_p=sigma_SI_p,
                scattering_type="SI",
                cross_section_type=cross_section_type,
                rho_chi=rho_chi,
                u_max=800.0,
                n_u=s["n_u"],
                v0=v0,
                include_thermal_targets=False,
                n_t_speed=s["n_t_speed"],
                n_t_mu=s["n_t_mu"],
                n_scatter_mu=s["n_scatter_mu"],
                n_scatter_phi=s["n_scatter_phi"]
            )

            c_sd_T = capture_rate_total(
                earth_data=earth_data,
                DM_mass=m,
                sigma_SD_p=sigma_SD_p,
                scattering_type="SD",
                cross_section_type=cross_section_type,
                rho_chi=rho_chi,
                u_max=800.0,
                n_u=s["n_u"],
                v0=v0,
                include_thermal_targets=True,
                n_t_speed=s["n_t_speed"],
                n_t_mu=s["n_t_mu"],
                n_scatter_mu=s["n_scatter_mu"],
                n_scatter_phi=s["n_scatter_phi"]
            )

            c_sd_0 = capture_rate_total(
                earth_data=earth_data,
                DM_mass=m,
                sigma_SD_p=sigma_SD_p,
                scattering_type="SD",
                cross_section_type=cross_section_type,
                rho_chi=rho_chi,
                u_max=800.0,
                n_u=s["n_u"],
                v0=v0,
                include_thermal_targets=False,
                n_t_speed=s["n_t_speed"],
                n_t_mu=s["n_t_mu"],
                n_scatter_mu=s["n_scatter_mu"],
                n_scatter_phi=s["n_scatter_phi"]
            )

            ratio_si = c_si_T / c_si_0 if c_si_0 > 0 else np.nan
            ratio_sd = c_sd_T / c_sd_0 if c_sd_0 > 0 else np.nan

            print(
                f"n_u={s['n_u']:2d}, "
                f"n_t_speed={s['n_t_speed']}, "
                f"n_t_mu={s['n_t_mu']}, "
                f"n_scatter_mu={s['n_scatter_mu']}, "
                f"n_scatter_phi={s['n_scatter_phi']} | "
                f"SI_T={c_si_T:.4e}, SI_0={c_si_0:.4e}, SI_ratio={ratio_si:.4e} | "
                f"SD_T={c_sd_T:.4e}, SD_0={c_sd_0:.4e}, SD_ratio={ratio_sd:.4e}"
            )


def test_si_convergence_1d(
    earth_data,
    masses,
    sigma_SI_p=1e-40,
    cross_section_type="constant",
    rho_chi=RHO_CHI_DEFAULT,
    v0=V0_DEFAULT
):
    """
    Focused one-parameter-at-a-time convergence test for SI capture.

    After the previous debug run, we already know:
      - n_t_speed is not the dominant source of error
      - n_t_mu is not the dominant source of error
      - n_scatter_phi is relatively mild

    So now we focus on:
      - n_u
      - n_scatter_mu
      - n_scatter_phi
      - shell_step
    """
    import time

    print("\n" + "=" * 80)
    print("SI focused convergence test")
    print("=" * 80)

    # Fixed baseline from previous results
    base = {
        "n_u": 40,
        "n_t_speed": 4,
        "n_t_mu": 4,
        "n_scatter_mu": 4,
        "n_scatter_phi": 6,
        "shell_step": 4,
    }

    scans = {
        "n_u": [20, 30, 40, 60],
        "n_scatter_mu": [3, 4, 5, 6],
        "n_scatter_phi": [4, 6, 8],
        "shell_step": [4, 2, 1],
    }

    for m in masses:
        print(f"\n--- SI focused convergence at m_chi = {m:.4g} GeV ---", flush=True)

        for param_name, values in scans.items():
            print(f"\nScan: {param_name}", flush=True)

            for val in values:
                s = dict(base)
                s[param_name] = val

                print(
                    f"[start] m={m:.4g}, scan={param_name}, val={val}, "
                    f"settings=(n_u={s['n_u']}, n_t_speed={s['n_t_speed']}, "
                    f"n_t_mu={s['n_t_mu']}, n_scatter_mu={s['n_scatter_mu']}, "
                    f"n_scatter_phi={s['n_scatter_phi']}, shell_step={s['shell_step']})",
                    flush=True
                )

                t0 = time.time()

                c_si_T = capture_rate_total(
                    earth_data=earth_data,
                    DM_mass=m,
                    sigma_SI_p=sigma_SI_p,
                    scattering_type="SI",
                    cross_section_type=cross_section_type,
                    rho_chi=rho_chi,
                    u_max=800.0,
                    n_u=s["n_u"],
                    v0=v0,
                    include_thermal_targets=True,
                    n_t_speed=s["n_t_speed"],
                    n_t_mu=s["n_t_mu"],
                    n_scatter_mu=s["n_scatter_mu"],
                    n_scatter_phi=s["n_scatter_phi"],
                    u_grid_mode="log",
                    shell_step=s["shell_step"]
                )

                c_si_0 = capture_rate_total(
                    earth_data=earth_data,
                    DM_mass=m,
                    sigma_SI_p=sigma_SI_p,
                    scattering_type="SI",
                    cross_section_type=cross_section_type,
                    rho_chi=rho_chi,
                    u_max=800.0,
                    n_u=s["n_u"],
                    v0=v0,
                    include_thermal_targets=False,
                    n_t_speed=s["n_t_speed"],
                    n_t_mu=s["n_t_mu"],
                    n_scatter_mu=s["n_scatter_mu"],
                    n_scatter_phi=s["n_scatter_phi"],
                    u_grid_mode="log",
                    shell_step=s["shell_step"]
                )

                dt = time.time() - t0
                ratio_si = c_si_T / c_si_0 if c_si_0 > 0 else np.nan

                print(
                    f"[done ] {param_name}={val:2d} | "
                    f"n_u={s['n_u']:2d}, "
                    f"n_t_speed={s['n_t_speed']:2d}, "
                    f"n_t_mu={s['n_t_mu']:2d}, "
                    f"n_scatter_mu={s['n_scatter_mu']:2d}, "
                    f"n_scatter_phi={s['n_scatter_phi']:2d}, "
                    f"shell_step={s['shell_step']:1d} | "
                    f"SI_T={c_si_T:.4e}, SI_0={c_si_0:.4e}, SI_ratio={ratio_si:.4e} | "
                    f"time={dt:.1f}s",
                    flush=True
                )

def plot_capture_rates_si_only(
    earth_data,
    DM_masses,
    sigma_SI_p=1e-40,
    cross_section_type="constant",
    save_path="capture_rates_SI_only.png",
    save_csv=True,
    csv_path="si_capture_results.csv",
    u_max=800.0,
    n_u=40,
    rho_chi=RHO_CHI_DEFAULT,
    v0=V0_DEFAULT,
    n_t_speed=4,
    n_t_mu=4,
    n_scatter_mu=5,
    n_scatter_phi=6,
    max_workers=None
):
    """
    SI-only plotting function.

    Left panel:
        C_SI(T != 0) and C_SI(T = 0)

    Right panel:
        C_SI(T != 0) / C_SI(T = 0)

    Also optionally saves the SI results to CSV.
    """
    import os
    from concurrent.futures import ProcessPoolExecutor, as_completed

    if max_workers is None:
        max_workers = max(1, os.cpu_count() - 1)

    print("\n" + "=" * 80)
    print("Generating SI-only capture-rate plot")
    print("=" * 80)
    print(f"Using {max_workers} worker processes")

    DM_masses = np.asarray(DM_masses, dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    tasks = []
    for m in DM_masses:
        tasks.append((
            earth_data,
            m,
            sigma_SI_p,     # sigma_SI_p
            0.0,            # sigma_SD_p dummy
            0.0,            # sigma_electron dummy
            cross_section_type,
            rho_chi,
            u_max,
            n_u,
            v0,
            n_t_speed,
            n_t_mu,
            n_scatter_mu,
            n_scatter_phi
        ))

    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(compute_one_mass_point, task) for task in tasks]
        for i, fut in enumerate(as_completed(futures), 1):
            res = fut.result()
            print(f"  done {i}/{len(futures)} : m_chi = {res['m']:.4g} GeV")
            results.append(res)

    results.sort(key=lambda x: x["m"])

    masses_sorted = np.array([r["m"] for r in results], dtype=float)
    C_SI_T = np.array([r["C_SI_T"] for r in results], dtype=float)
    C_SI_0 = np.array([r["C_SI_0"] for r in results], dtype=float)

    # Raw ratio for saving
    ratio_SI_raw = np.divide(
        C_SI_T,
        C_SI_0,
        out=np.full_like(C_SI_T, np.nan, dtype=float),
        where=(C_SI_0 != 0.0)
    )

    # Ratio for plotting: keep all physically nonzero points
    ratio_SI_plot = np.divide(
        C_SI_T,
        C_SI_0,
        out=np.full_like(C_SI_T, np.nan, dtype=float),
        where=(C_SI_0 > 0.0)
    )
    # ------------------------------------------------------------
    # Save CSV
    # ------------------------------------------------------------
    if save_csv:
        csv_dir = os.path.dirname(csv_path)
        if csv_dir:
            os.makedirs(csv_dir, exist_ok=True)

        save_si_results_to_csv(
            DM_masses=masses_sorted,
            C_SI_T=C_SI_T,
            C_SI_0=C_SI_0,
            output_csv=csv_path
        )

    # ------------------------------------------------------------
    # Left panel
    # ------------------------------------------------------------
    axL = axes[0]
    axL.loglog(masses_sorted, C_SI_T, "b-.", lw=2.5, label="SI, T(r)≠0")
    axL.loglog(masses_sorted, C_SI_0, color="blue", alpha=0.35, ls=":", lw=2.2, label="SI, T=0")
    axL.set_xlabel(r"$m_\chi$ [GeV]", fontsize=12)
    axL.set_ylabel(r"$C$ [s$^{-1}$]", fontsize=12)
    axL.set_title(f"SI only, {cross_section_type}, sigma = {sigma_SI_p:.0e} cm$^2$", fontsize=12)
    axL.grid(True, alpha=0.3)
    axL.legend(fontsize=10, loc="best")

    # ------------------------------------------------------------
    # Right panel
    # ------------------------------------------------------------
    axR = axes[1]
    mask_mass = (masses_sorted >= 1e-2) & np.isfinite(ratio_SI_plot)    
    axR.semilogx(masses_sorted[mask_mass], ratio_SI_plot[mask_mass], "b-.", lw=2.2, label="SI")
    axR.axhline(1.0, color="k", linestyle=":", lw=1.2, label=r"$T(r)=0$ limit")
    axR.set_xlabel(r"$m_\chi$ [GeV]", fontsize=12)
    axR.set_ylabel(r"$C(T)/C(T=0)$", fontsize=12)
    axR.set_title("SI thermal correction factor", fontsize=12)
    axR.set_yscale("log")
    axR.grid(True, alpha=0.3)
    axR.legend(fontsize=10, loc="best")

    plt.tight_layout()

    fig_dir = os.path.dirname(save_path)
    if fig_dir:
        os.makedirs(fig_dir, exist_ok=True)

    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"[INFO] Plot saved as: {os.path.abspath(save_path)}")

    plt.close(fig)

    return {
        "DM_masses": masses_sorted,
        "C_SI_T": C_SI_T,
        "C_SI_0": C_SI_0,
        "ratio_SI": ratio_SI_raw
    }


def save_si_results_to_csv(
    DM_masses,
    C_SI_T,
    C_SI_0,
    output_csv="si_capture_results.csv"
):
    """
    Save SI-only capture results to CSV.

    Columns:
        m_chi_GeV
        C_SI_T
        C_SI_0
        ratio_SI

    Parameters
    ----------
    DM_masses : array-like
        Dark matter masses in GeV.
    C_SI_T : array-like
        SI capture rates with finite temperature, T(r) != 0.
    C_SI_0 : array-like
        SI capture rates in the T = 0 limit.
    output_csv : str
        Output CSV filename.

    Returns
    -------
    ratio_SI : np.ndarray
        The ratio C_SI_T / C_SI_0.
    """
    DM_masses = np.asarray(DM_masses, dtype=float)
    C_SI_T = np.asarray(C_SI_T, dtype=float)
    C_SI_0 = np.asarray(C_SI_0, dtype=float)

    if not (len(DM_masses) == len(C_SI_T) == len(C_SI_0)):
        raise ValueError(
            "DM_masses, C_SI_T, and C_SI_0 must have the same length."
        )

    ratio_SI = np.divide(
        C_SI_T,
        C_SI_0,
        out=np.full_like(C_SI_T, np.nan, dtype=float),
        where=(C_SI_0 != 0.0)
    )

    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["m_chi_GeV", "C_SI_T", "C_SI_0", "ratio_SI"])

        for m, c_t, c_0, r in zip(DM_masses, C_SI_T, C_SI_0, ratio_SI):
            writer.writerow([
                f"{m:.8e}",
                f"{c_t:.8e}",
                f"{c_0:.8e}",
                f"{r:.8e}" if np.isfinite(r) else "nan"
            ])

    print(f"[INFO] SI results saved to: {output_csv}")
    return ratio_SI

def save_si_run_info(
    output_txt,
    cross_section_type,
    sigma_SI_p,
    DM_masses,
    results,
    u_max,
    n_u,
    n_t_speed,
    n_t_mu,
    n_scatter_mu,
    n_scatter_phi,
    rho_chi,
    v0,
    max_workers
):
    """
    Save a simple text summary for one SI-only run.
    """
    ratio = np.asarray(results["ratio_SI"], dtype=float)
    masses = np.asarray(results["DM_masses"], dtype=float)
    C_SI_T = np.asarray(results["C_SI_T"], dtype=float)
    C_SI_0 = np.asarray(results["C_SI_0"], dtype=float)

    finite_mask = np.isfinite(ratio)
    if np.any(finite_mask):
        idx_max = np.nanargmax(ratio)
        max_ratio = ratio[idx_max]
        max_ratio_mass = masses[idx_max]
    else:
        max_ratio = np.nan
        max_ratio_mass = np.nan

    out_dir = os.path.dirname(output_txt)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(output_txt, "w", encoding="utf-8") as f:
        f.write("SI-only Earth capture run info\n")
        f.write("=" * 60 + "\n")
        f.write(f"cross_section_type = {cross_section_type}\n")
        f.write(f"sigma_SI_p         = {sigma_SI_p:.6e} cm^2\n")
        f.write(f"rho_chi            = {rho_chi:.6e} GeV/cm^3\n")
        f.write(f"v0                 = {v0:.6f} km/s\n")
        f.write(f"u_max              = {u_max:.6f} km/s\n")
        f.write(f"n_u                = {n_u}\n")
        f.write(f"n_t_speed          = {n_t_speed}\n")
        f.write(f"n_t_mu             = {n_t_mu}\n")
        f.write(f"n_scatter_mu       = {n_scatter_mu}\n")
        f.write(f"n_scatter_phi      = {n_scatter_phi}\n")
        f.write(f"max_workers        = {max_workers}\n")
        f.write(f"DM mass range      = {masses.min():.6e} -> {masses.max():.6e} GeV\n")
        f.write(f"Number of masses   = {len(masses)}\n")
        f.write("\n")
        f.write("Summary:\n")
        f.write(f"max ratio C(T)/C(T=0) = {max_ratio:.6e}\n")
        f.write(f"at m_chi              = {max_ratio_mass:.6e} GeV\n")
        f.write(f"C_SI_T(min,max)       = ({np.nanmin(C_SI_T):.6e}, {np.nanmax(C_SI_T):.6e})\n")
        f.write(f"C_SI_0(min,max)       = ({np.nanmin(C_SI_0):.6e}, {np.nanmax(C_SI_0):.6e})\n")

    print(f"[INFO] Run info saved to: {output_txt}")

def run_si_only_all_cross_sections(
    earth_data,
    DM_masses,
    sigma_SI_p=1e-40,
    output_root=".",
    u_max=800.0,
    n_u=40,
    rho_chi=RHO_CHI_DEFAULT,
    v0=V0_DEFAULT,
    n_t_speed=4,
    n_t_mu=4,
    n_scatter_mu=5,
    n_scatter_phi=6,
    max_workers=None
):
    """
    Run SI-only capture for:
        1. constant
        2. v2_dependent
        3. q2_dependent

    For each case, save:
        - figure
        - CSV
        - run-info text
    """
    if max_workers is None:
        max_workers = max(1, os.cpu_count() - 1)

    DM_masses = np.asarray(DM_masses, dtype=float)

    cross_section_types = [
        "constant",
        "v2_dependent",
        "q2_dependent",
    ]

    short_tag = {
        "constant": "constant",
        "v2_dependent": "v2",
        "q2_dependent": "q2",
    }

    sigma_tag = f"{sigma_SI_p:.0e}".replace("+", "")

    all_results = {}

    figures_dir = os.path.join(output_root, "figures")
    results_dir = os.path.join(output_root, "results")
    logs_dir = os.path.join(output_root, "logs")

    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    print("\n" + "=" * 80)
    print("Running SI-only suite: constant / v2_dependent / q2_dependent")
    print("=" * 80)

    for i, cross_type in enumerate(cross_section_types, 1):
        tag = short_tag[cross_type]

        fig_path = os.path.join(
            figures_dir,
            f"capture_rates_SI_only_{tag}_sigma{sigma_tag}.png"
        )
        csv_path = os.path.join(
            results_dir,
            f"si_capture_results_{tag}_sigma{sigma_tag}.csv"
        )
        log_path = os.path.join(
            logs_dir,
            f"si_run_info_{tag}_sigma{sigma_tag}.txt"
        )

        print("\n" + "-" * 80)
        print(f"[{i}/{len(cross_section_types)}] Running SI-only case: {cross_type}")
        print("-" * 80)

        results = plot_capture_rates_si_only(
            earth_data=earth_data,
            DM_masses=DM_masses,
            sigma_SI_p=sigma_SI_p,
            cross_section_type=cross_type,
            save_path=fig_path,
            save_csv=True,
            csv_path=csv_path,
            u_max=u_max,
            n_u=n_u,
            rho_chi=rho_chi,
            v0=v0,
            n_t_speed=n_t_speed,
            n_t_mu=n_t_mu,
            n_scatter_mu=n_scatter_mu,
            n_scatter_phi=n_scatter_phi,
            max_workers=max_workers
        )

        save_si_run_info(
            output_txt=log_path,
            cross_section_type=cross_type,
            sigma_SI_p=sigma_SI_p,
            DM_masses=DM_masses,
            results=results,
            u_max=u_max,
            n_u=n_u,
            n_t_speed=n_t_speed,
            n_t_mu=n_t_mu,
            n_scatter_mu=n_scatter_mu,
            n_scatter_phi=n_scatter_phi,
            rho_chi=rho_chi,
            v0=v0,
            max_workers=max_workers
        )

        all_results[cross_type] = results

    print("\n" + "=" * 80)
    print("All SI-only runs finished.")
    print("=" * 80)

    return all_results

def redraw_q2_si_from_csv(
    csv_path="si_capture_results_q2_sigma1e-40.csv",
    save_path="capture_rates_SI_only_q2_sigma1e-40_redrawn.png",
    sigma_SI_p=1e-40
):
    """
    Redraw the SI-only q2-dependent figure from an existing CSV file.
    No capture-rate recalculation is performed.
    """
    import os
    import csv
    import numpy as np
    import matplotlib.pyplot as plt

    masses = []
    C_SI_T = []
    C_SI_0 = []

    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            masses.append(float(row["m_chi_GeV"]))
            C_SI_T.append(float(row["C_SI_T"]))
            C_SI_0.append(float(row["C_SI_0"]))

    masses = np.asarray(masses, dtype=float)
    C_SI_T = np.asarray(C_SI_T, dtype=float)
    C_SI_0 = np.asarray(C_SI_0, dtype=float)

    ratio_SI_plot = np.divide(
        C_SI_T,
        C_SI_0,
        out=np.full_like(C_SI_T, np.nan, dtype=float),
        where=(C_SI_0 > 0.0)
    )

    mask_mass = (masses >= 1e-2) & np.isfinite(ratio_SI_plot)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left panel
    axL = axes[0]
    axL.loglog(masses, C_SI_T, "b-.", lw=2.5, label="SI, T(r)≠0")
    axL.loglog(masses, C_SI_0, color="blue", alpha=0.35, ls=":", lw=2.2, label="SI, T=0")
    axL.set_xlabel(r"$m_\chi$ [GeV]", fontsize=12)
    axL.set_ylabel(r"$C$ [s$^{-1}$]", fontsize=12)
    axL.set_title(f"SI only, q2_dependent, sigma = {sigma_SI_p:.0e} cm$^2$", fontsize=12)
    axL.grid(True, alpha=0.3)
    axL.legend(fontsize=10, loc="best")

    # Right panel
    axR = axes[1]
    axR.semilogx(masses[mask_mass], ratio_SI_plot[mask_mass], "b-.", lw=2.2, label="SI")
    axR.axhline(1.0, color="k", linestyle=":", lw=1.2, label=r"$T(r)=0$ limit")
    axR.set_xlabel(r"$m_\chi$ [GeV]", fontsize=12)
    axR.set_ylabel(r"$C(T)/C(T=0)$", fontsize=12)
    axR.set_title("SI thermal correction factor", fontsize=12)
    axR.set_yscale("log")
    axR.grid(True, alpha=0.3)
    axR.legend(fontsize=10, loc="best")

    plt.tight_layout()

    out_dir = os.path.dirname(save_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"[INFO] Redrawn q2 figure saved as: {os.path.abspath(save_path)}")
    plt.close(fig)

def plot_si_thermal_ratio_comparison_from_csv(
    csv_constant="results/si_capture_results_constant_sigma1e-40.csv",
    csv_v2="results/si_capture_results_v2_sigma1e-40.csv",
    csv_q2="results/si_capture_results_q2_sigma1e-40.csv",
    save_path="figures/si_thermal_ratio_comparison_sigma1e-40.png"
):
    """
    Read three SI-only CSV files and plot a comparison of
        C(T) / C(T=0)
    for:
        - constant
        - v2_dependent
        - q2_dependent
    """
    import os
    import csv
    import numpy as np
    import matplotlib.pyplot as plt

    def load_ratio_from_csv(csv_path):
        masses = []
        C_SI_T = []
        C_SI_0 = []

        with open(csv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                masses.append(float(row["m_chi_GeV"]))
                C_SI_T.append(float(row["C_SI_T"]))
                C_SI_0.append(float(row["C_SI_0"]))

        masses = np.asarray(masses, dtype=float)
        C_SI_T = np.asarray(C_SI_T, dtype=float)
        C_SI_0 = np.asarray(C_SI_0, dtype=float)

        ratio = np.divide(
            C_SI_T,
            C_SI_0,
            out=np.full_like(C_SI_T, np.nan, dtype=float),
            where=(C_SI_0 > 0.0)
        )

        mask = (masses >= 1e-2) & np.isfinite(ratio)
        return masses, ratio, mask

    m_const, r_const, mask_const = load_ratio_from_csv(csv_constant)
    m_v2, r_v2, mask_v2 = load_ratio_from_csv(csv_v2)
    m_q2, r_q2, mask_q2 = load_ratio_from_csv(csv_q2)

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.semilogx(
        m_const[mask_const], r_const[mask_const],
        color="blue", linestyle="-", linewidth=2.4,
        label="constant"
    )
    ax.semilogx(
        m_v2[mask_v2], r_v2[mask_v2],
        color="green", linestyle="--", linewidth=2.4,
        label=r"$v^2$ dependent"
    )
    ax.semilogx(
        m_q2[mask_q2], r_q2[mask_q2],
        color="red", linestyle="-.", linewidth=2.4,
        label=r"$q^2$ dependent"
    )

    ax.axhline(
        1.0,
        color="black",
        linestyle=":",
        linewidth=1.3,
        label=r"$T(r)=0$ limit"
    )

    ax.set_xlabel(r"$m_\chi\ [{\rm GeV}]$", fontsize=13)
    ax.set_ylabel(r"$C(T)/C(T=0)$", fontsize=13)
    ax.set_title(r"SI thermal correction comparison", fontsize=14)
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=11, loc="best")

    plt.tight_layout()

    out_dir = os.path.dirname(save_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"[INFO] SI ratio comparison figure saved as: {os.path.abspath(save_path)}")
    plt.close(fig)

def compute_one_mass_point_electron(task):
    (
        earth_data,
        m,
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

    print(f"[worker:e] start m = {m:.4g} GeV", flush=True)

    c_e = capture_rate_total(
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

    print(f"[worker:e] done  m = {m:.4g} GeV", flush=True)

    return {
        "m": m,
        "C_electron": c_e,
    }
def save_electron_results_to_csv(
    DM_masses,
    C_electron,
    output_csv="electron_capture_results.csv"
):
    """
    Save electron-only baseline capture results to CSV.

    Columns:
        m_chi_GeV
        C_electron_T0
    """
    import csv
    import numpy as np

    DM_masses = np.asarray(DM_masses, dtype=float)
    C_electron = np.asarray(C_electron, dtype=float)

    if len(DM_masses) != len(C_electron):
        raise ValueError("DM_masses and C_electron must have the same length.")

    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["m_chi_GeV", "C_electron_T0"])

        for m, c in zip(DM_masses, C_electron):
            writer.writerow([
                f"{m:.8e}",
                f"{c:.8e}"
            ])

    print(f"[INFO] Electron results saved to: {output_csv}")

def save_electron_run_info(
    output_txt,
    cross_section_type,
    sigma_electron,
    results,
    u_max,
    n_u,
    n_t_speed,
    n_t_mu,
    n_scatter_mu,
    n_scatter_phi,
    rho_chi,
    v0,
    max_workers
):
    """
    Save a simple text summary for one electron-only run.
    """
    import os
    import numpy as np

    masses = np.asarray(results["DM_masses"], dtype=float)
    C_electron = np.asarray(results["C_electron"], dtype=float)

    out_dir = os.path.dirname(output_txt)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(output_txt, "w", encoding="utf-8") as f:
        f.write("Electron-only Earth capture baseline run info\n")
        f.write("=" * 60 + "\n")
        f.write("NOTE: electron thermal motion is disabled in current baseline.\n\n")
        f.write(f"cross_section_type = {cross_section_type}\n")
        f.write(f"sigma_electron     = {sigma_electron:.6e} cm^2\n")
        f.write(f"rho_chi            = {rho_chi:.6e} GeV/cm^3\n")
        f.write(f"v0                 = {v0:.6f} km/s\n")
        f.write(f"u_max              = {u_max:.6f} km/s\n")
        f.write(f"n_u                = {n_u}\n")
        f.write(f"n_t_speed          = {n_t_speed}\n")
        f.write(f"n_t_mu             = {n_t_mu}\n")
        f.write(f"n_scatter_mu       = {n_scatter_mu}\n")
        f.write(f"n_scatter_phi      = {n_scatter_phi}\n")
        f.write(f"max_workers        = {max_workers}\n")
        f.write(f"DM mass range      = {masses.min():.6e} -> {masses.max():.6e} GeV\n")
        f.write(f"Number of masses   = {len(masses)}\n")
        f.write("\n")
        f.write("Summary:\n")
        f.write(f"C_electron(min,max) = ({np.nanmin(C_electron):.6e}, {np.nanmax(C_electron):.6e})\n")

    print(f"[INFO] Electron run info saved to: {output_txt}")

def plot_capture_rates_electron_only(
    earth_data,
    DM_masses,
    sigma_electron=1e-40,
    cross_section_type="constant",
    save_path="capture_rates_electron_only.png",
    save_csv=True,
    csv_path="electron_capture_results.csv",
    u_max=800.0,
    n_u=40,
    rho_chi=RHO_CHI_DEFAULT,
    v0=V0_DEFAULT,
    n_t_speed=4,
    n_t_mu=4,
    n_scatter_mu=5,
    n_scatter_phi=6,
    max_workers=None
):
    """
    Electron-only baseline plotting function.

    Current status:
        electron thermal motion is disabled,
        so this plots the baseline T=0-equivalent electron capture rate only.
    """
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from concurrent.futures import ProcessPoolExecutor, as_completed

    if max_workers is None:
        max_workers = max(1, os.cpu_count() - 1)

    print("\n" + "=" * 80)
    print("Generating electron-only baseline capture-rate plot")
    print("=" * 80)
    print(f"Using {max_workers} worker processes")

    DM_masses = np.asarray(DM_masses, dtype=float)

    tasks = []
    for m in DM_masses:
        tasks.append((
            earth_data,
            m,
            sigma_electron,
            cross_section_type,
            rho_chi,
            u_max,
            n_u,
            v0,
            n_t_speed,
            n_t_mu,
            n_scatter_mu,
            n_scatter_phi
        ))

    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(compute_one_mass_point_electron, task) for task in tasks]
        for i, fut in enumerate(as_completed(futures), 1):
            res = fut.result()
            print(f"  done {i}/{len(futures)} : m_chi = {res['m']:.4g} GeV")
            results.append(res)

    results.sort(key=lambda x: x["m"])

    masses_sorted = np.array([r["m"] for r in results], dtype=float)
    C_electron = np.array([r["C_electron"] for r in results], dtype=float)

    if save_csv:
        csv_dir = os.path.dirname(csv_path)
        if csv_dir:
            os.makedirs(csv_dir, exist_ok=True)

        save_electron_results_to_csv(
            DM_masses=masses_sorted,
            C_electron=C_electron,
            output_csv=csv_path
        )

    fig, ax = plt.subplots(1, 1, figsize=(8.5, 6))

    ax.loglog(masses_sorted, C_electron, color="red", linestyle="-", linewidth=2.5, label="Electron baseline")
    ax.set_xlabel(r"$m_\chi\ [{\rm GeV}]$", fontsize=13)
    ax.set_ylabel(r"$C\ [{\rm s}^{-1}]$", fontsize=13)
    ax.set_title(f"Electron only, {cross_section_type}, sigma = {sigma_electron:.0e} cm$^2$", fontsize=13)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=11, loc="best")

    ax.text(
        0.03, 0.06,
        "Current baseline: electron thermal motion disabled",
        transform=ax.transAxes,
        fontsize=10,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
    )

    plt.tight_layout()

    fig_dir = os.path.dirname(save_path)
    if fig_dir:
        os.makedirs(fig_dir, exist_ok=True)

    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"[INFO] Electron figure saved as: {os.path.abspath(save_path)}")
    plt.close(fig)

    return {
        "DM_masses": masses_sorted,
        "C_electron": C_electron
    }


def run_electron_only_all_cross_sections(
    earth_data,
    DM_masses,
    sigma_electron=1e-40,
    output_root=".",
    u_max=800.0,
    n_u=40,
    rho_chi=RHO_CHI_DEFAULT,
    v0=V0_DEFAULT,
    n_t_speed=4,
    n_t_mu=4,
    n_scatter_mu=5,
    n_scatter_phi=6,
    max_workers=None
):
    """
    Run electron-only baseline capture for:
        1. constant
        2. v2_dependent
        3. q2_dependent

    For each case, save:
        - figure
        - CSV
        - run-info text
    """
    import os
    import numpy as np

    if max_workers is None:
        max_workers = max(1, os.cpu_count() - 1)

    DM_masses = np.asarray(DM_masses, dtype=float)

    cross_section_types = [
        "constant",
        "v2_dependent",
        "q2_dependent",
    ]

    short_tag = {
        "constant": "constant",
        "v2_dependent": "v2",
        "q2_dependent": "q2",
    }

    sigma_tag = f"{sigma_electron:.0e}".replace("+", "")

    figures_dir = os.path.join(output_root, "figures")
    results_dir = os.path.join(output_root, "results")
    logs_dir = os.path.join(output_root, "logs")

    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    all_results = {}

    print("\n" + "=" * 80)
    print("Running electron-only suite: constant / v2_dependent / q2_dependent")
    print("=" * 80)

    for i, cross_type in enumerate(cross_section_types, 1):
        tag = short_tag[cross_type]

        fig_path = os.path.join(
            figures_dir,
            f"capture_rates_electron_only_{tag}_sigma{sigma_tag}.png"
        )
        csv_path = os.path.join(
            results_dir,
            f"electron_capture_results_{tag}_sigma{sigma_tag}.csv"
        )
        log_path = os.path.join(
            logs_dir,
            f"electron_run_info_{tag}_sigma{sigma_tag}.txt"
        )

        print("\n" + "-" * 80)
        print(f"[{i}/{len(cross_section_types)}] Running electron-only case: {cross_type}")
        print("-" * 80)

        results = plot_capture_rates_electron_only(
            earth_data=earth_data,
            DM_masses=DM_masses,
            sigma_electron=sigma_electron,
            cross_section_type=cross_type,
            save_path=fig_path,
            save_csv=True,
            csv_path=csv_path,
            u_max=u_max,
            n_u=n_u,
            rho_chi=rho_chi,
            v0=v0,
            n_t_speed=n_t_speed,
            n_t_mu=n_t_mu,
            n_scatter_mu=n_scatter_mu,
            n_scatter_phi=n_scatter_phi,
            max_workers=max_workers
        )

        save_electron_run_info(
            output_txt=log_path,
            cross_section_type=cross_type,
            sigma_electron=sigma_electron,
            results=results,
            u_max=u_max,
            n_u=n_u,
            n_t_speed=n_t_speed,
            n_t_mu=n_t_mu,
            n_scatter_mu=n_scatter_mu,
            n_scatter_phi=n_scatter_phi,
            rho_chi=rho_chi,
            v0=v0,
            max_workers=max_workers
        )

        all_results[cross_type] = results

    print("\n" + "=" * 80)
    print("All electron-only runs finished.")
    print("=" * 80)

    return all_results

def plot_electron_comparison_from_csv(
    csv_constant="results/electron_capture_results_constant_sigma1e-40.csv",
    csv_v2="results/electron_capture_results_v2_sigma1e-40.csv",
    csv_q2="results/electron_capture_results_q2_sigma1e-40.csv",
    save_path="figures/electron_capture_comparison_sigma1e-40.png"
):
    """
    Read three electron-only baseline CSV files and plot a comparison of
    absolute capture rates C(m_chi) for:
        - constant
        - v2_dependent
        - q2_dependent
    """
    import os
    import csv
    import numpy as np
    import matplotlib.pyplot as plt

    def load_electron_csv(csv_path):
        masses = []
        C_electron = []

        with open(csv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                masses.append(float(row["m_chi_GeV"]))
                C_electron.append(float(row["C_electron_T0"]))

        masses = np.asarray(masses, dtype=float)
        C_electron = np.asarray(C_electron, dtype=float)

        mask = np.isfinite(masses) & np.isfinite(C_electron) & (masses > 0.0) & (C_electron > 0.0)
        return masses[mask], C_electron[mask]

    m_const, C_const = load_electron_csv(csv_constant)
    m_v2, C_v2 = load_electron_csv(csv_v2)
    m_q2, C_q2 = load_electron_csv(csv_q2)

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.loglog(
        m_const, C_const,
        color="red", linestyle="-", linewidth=2.5,
        label="constant"
    )
    ax.loglog(
        m_v2, C_v2,
        color="blue", linestyle="--", linewidth=2.5,
        label=r"$v^2$ dependent"
    )
    ax.loglog(
        m_q2, C_q2,
        color="green", linestyle="-.", linewidth=2.5,
        label=r"$q^2$ dependent"
    )

    ax.set_xlabel(r"$m_\chi\ [{\rm GeV}]$", fontsize=13)
    ax.set_ylabel(r"$C\ [{\rm s}^{-1}]$", fontsize=13)
    ax.set_title(r"Electron-only capture comparison", fontsize=14)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=11, loc="best")

    ax.text(
        0.03, 0.06,
        "Baseline only: electron thermal motion disabled",
        transform=ax.transAxes,
        fontsize=10,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85)
    )

    plt.tight_layout()

    out_dir = os.path.dirname(save_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"[INFO] Electron comparison figure saved as: {os.path.abspath(save_path)}")
    plt.close(fig)

def compute_one_mass_point_sd_baseline(task):
    (
        earth_data,
        m,
        sigma_SD_p,
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

    print(f"[worker:sd] start m = {m:.4g} GeV", flush=True)

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

    print(f"[worker:sd] done  m = {m:.4g} GeV", flush=True)

    return {
        "m": m,
        "C_SD_0": c_sd_0,
    }

def save_sd_baseline_results_to_csv(
    DM_masses,
    C_SD_0,
    output_csv="sd_capture_results_constant_T0.csv"
):
    """
    Save SD-only T=0 baseline capture results to CSV.

    Columns:
        m_chi_GeV
        C_SD_0
    """
    import csv
    import numpy as np

    DM_masses = np.asarray(DM_masses, dtype=float)
    C_SD_0 = np.asarray(C_SD_0, dtype=float)

    if len(DM_masses) != len(C_SD_0):
        raise ValueError("DM_masses and C_SD_0 must have the same length.")

    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["m_chi_GeV", "C_SD_0"])

        for m, c in zip(DM_masses, C_SD_0):
            writer.writerow([
                f"{m:.8e}",
                f"{c:.8e}"
            ])

    print(f"[INFO] SD baseline results saved to: {output_csv}")

def save_sd_baseline_run_info(
    output_txt,
    sigma_SD_p,
    results,
    u_max,
    n_u,
    n_t_speed,
    n_t_mu,
    n_scatter_mu,
    n_scatter_phi,
    rho_chi,
    v0,
    max_workers
):
    """
    Save a simple text summary for one SD-only T=0 baseline run.
    """
    import os
    import numpy as np

    masses = np.asarray(results["DM_masses"], dtype=float)
    C_SD_0 = np.asarray(results["C_SD_0"], dtype=float)

    out_dir = os.path.dirname(output_txt)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(output_txt, "w", encoding="utf-8") as f:
        f.write("SD-only Earth capture baseline run info\n")
        f.write("=" * 60 + "\n")
        f.write("NOTE: this is SD-only, constant cross section, T=0 baseline.\n")
        f.write("NOTE: current SD implementation should be treated as diagnostic/preliminary.\n\n")
        f.write(f"cross_section_type = constant\n")
        f.write(f"sigma_SD_p         = {sigma_SD_p:.6e} cm^2\n")
        f.write(f"rho_chi            = {rho_chi:.6e} GeV/cm^3\n")
        f.write(f"v0                 = {v0:.6f} km/s\n")
        f.write(f"u_max              = {u_max:.6f} km/s\n")
        f.write(f"n_u                = {n_u}\n")
        f.write(f"n_t_speed          = {n_t_speed}\n")
        f.write(f"n_t_mu             = {n_t_mu}\n")
        f.write(f"n_scatter_mu       = {n_scatter_mu}\n")
        f.write(f"n_scatter_phi      = {n_scatter_phi}\n")
        f.write(f"max_workers        = {max_workers}\n")
        f.write(f"DM mass range      = {masses.min():.6e} -> {masses.max():.6e} GeV\n")
        f.write(f"Number of masses   = {len(masses)}\n")
        f.write("\n")
        f.write("Summary:\n")
        f.write(f"C_SD_0(min,max) = ({np.nanmin(C_SD_0):.6e}, {np.nanmax(C_SD_0):.6e})\n")

    print(f"[INFO] SD baseline run info saved to: {output_txt}")

def plot_capture_rates_sd_baseline(
    earth_data,
    DM_masses,
    sigma_SD_p=1e-40,
    save_path="capture_rates_SD_only_baseline_constant.png",
    save_csv=True,
    csv_path="sd_capture_results_constant_T0.csv",
    u_max=800.0,
    n_u=40,
    rho_chi=RHO_CHI_DEFAULT,
    v0=V0_DEFAULT,
    n_t_speed=4,
    n_t_mu=4,
    n_scatter_mu=5,
    n_scatter_phi=6,
    max_workers=None
):
    """
    SD-only baseline plotting function.

    Current version:
        - SD only
        - constant cross section
        - T = 0 baseline
    """
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from concurrent.futures import ProcessPoolExecutor, as_completed

    if max_workers is None:
        max_workers = max(1, os.cpu_count() - 1)

    print("\n" + "=" * 80)
    print("Generating SD-only baseline capture-rate plot")
    print("=" * 80)
    print(f"Using {max_workers} worker processes")

    DM_masses = np.asarray(DM_masses, dtype=float)

    tasks = []
    for m in DM_masses:
        tasks.append((
            earth_data,
            m,
            sigma_SD_p,
            "constant",
            rho_chi,
            u_max,
            n_u,
            v0,
            n_t_speed,
            n_t_mu,
            n_scatter_mu,
            n_scatter_phi
        ))

    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(compute_one_mass_point_sd_baseline, task) for task in tasks]
        for i, fut in enumerate(as_completed(futures), 1):
            res = fut.result()
            print(f"  done {i}/{len(futures)} : m_chi = {res['m']:.4g} GeV")
            results.append(res)

    results.sort(key=lambda x: x["m"])

    masses_sorted = np.array([r["m"] for r in results], dtype=float)
    C_SD_0 = np.array([r["C_SD_0"] for r in results], dtype=float)

    if save_csv:
        csv_dir = os.path.dirname(csv_path)
        if csv_dir:
            os.makedirs(csv_dir, exist_ok=True)

        save_sd_baseline_results_to_csv(
            DM_masses=masses_sorted,
            C_SD_0=C_SD_0,
            output_csv=csv_path
        )

    fig, ax = plt.subplots(1, 1, figsize=(8.5, 6))

    ax.loglog(
        masses_sorted, C_SD_0,
        color="purple", linestyle="-", linewidth=2.5,
        label="SD baseline, T=0"
    )
    ax.set_xlabel(r"$m_\chi\ [{\rm GeV}]$", fontsize=13)
    ax.set_ylabel(r"$C\ [{\rm s}^{-1}]$", fontsize=13)
    ax.set_title(r"SD only, constant, T = 0 baseline", fontsize=13)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=11, loc="best")

    ax.text(
        0.03, 0.06,
        "Current SD baseline: constant cross section, thermal motion disabled",
        transform=ax.transAxes,
        fontsize=10,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
    )

    plt.tight_layout()

    fig_dir = os.path.dirname(save_path)
    if fig_dir:
        os.makedirs(fig_dir, exist_ok=True)

    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"[INFO] SD baseline figure saved as: {os.path.abspath(save_path)}")
    plt.close(fig)

    return {
        "DM_masses": masses_sorted,
        "C_SD_0": C_SD_0
    }

def run_sd_baseline_constant(
    earth_data,
    DM_masses,
    sigma_SD_p=1e-40,
    output_root=".",
    u_max=800.0,
    n_u=40,
    rho_chi=RHO_CHI_DEFAULT,
    v0=V0_DEFAULT,
    n_t_speed=4,
    n_t_mu=4,
    n_scatter_mu=5,
    n_scatter_phi=6,
    max_workers=None
):
    """
    Minimal SD baseline runner:
        - SD only
        - constant
        - T = 0
    """
    import os
    import numpy as np

    if max_workers is None:
        max_workers = max(1, os.cpu_count() - 1)

    DM_masses = np.asarray(DM_masses, dtype=float)

    sigma_tag = f"{sigma_SD_p:.0e}".replace("+", "")

    figures_dir = os.path.join(output_root, "figures")
    results_dir = os.path.join(output_root, "results")
    logs_dir = os.path.join(output_root, "logs")

    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    fig_path = os.path.join(
        figures_dir,
        f"capture_rates_SD_only_baseline_constant_sigma{sigma_tag}.png"
    )
    csv_path = os.path.join(
        results_dir,
        f"sd_capture_results_constant_T0_sigma{sigma_tag}.csv"
    )
    log_path = os.path.join(
        logs_dir,
        f"sd_run_info_baseline_constant_sigma{sigma_tag}.txt"
    )

    results = plot_capture_rates_sd_baseline(
        earth_data=earth_data,
        DM_masses=DM_masses,
        sigma_SD_p=sigma_SD_p,
        save_path=fig_path,
        save_csv=True,
        csv_path=csv_path,
        u_max=u_max,
        n_u=n_u,
        rho_chi=rho_chi,
        v0=v0,
        n_t_speed=n_t_speed,
        n_t_mu=n_t_mu,
        n_scatter_mu=n_scatter_mu,
        n_scatter_phi=n_scatter_phi,
        max_workers=max_workers
    )

def get_verified_only_active_sd_labels(earth_data):
    """
    Collect active SD labels actually present in earth_data["active_shells"].
    """
    active_labels = set()
    label_status = {}
    for shell in earth_data["active_shells"]:
        for info in shell["SD"]:
            label = info["elem"]
            status = info.get("status", "unknown")
            active_labels.add(label)
            label_status.setdefault(label, set()).add(status)
    return sorted(active_labels), label_status


def validate_verified_only_sd_baseline(
    earth_data,
    expected_labels=("Al", "H", "Na", "Si29")
):
    """
    Hard validation for the verified-only refined SD baseline.
    """
    if earth_data.get("sd_mode", None) != "verified_only":
        raise RuntimeError(
            "earth_data['sd_mode'] is not 'verified_only'. "
            "Reload Earth composition with sd_mode='verified_only'."
        )

    active_labels, label_status = get_verified_only_active_sd_labels(earth_data)

    print("\n" + "=" * 80)
    print("Validated active SD labels for verified-only baseline")
    print("=" * 80)
    print("Active labels under current mode:")
    print(active_labels)
    print(f"Count = {len(active_labels)}")

    if set(active_labels) != set(expected_labels):
        raise RuntimeError(
            "Verified-only SD baseline validation failed.\n"
            f"Expected labels: {sorted(expected_labels)}\n"
            f"Actual labels:   {active_labels}"
        )

    bad_status = {}
    for label, statuses in label_status.items():
        if statuses != {"verified"}:
            bad_status[label] = statuses

    if bad_status:
        raise RuntimeError(
            "Some active labels in verified_only mode are not purely verified:\n"
            f"{bad_status}"
        )

    print("Status check passed: all active SD labels are verified.")
    print("=" * 80)


def compute_one_mass_point_verified_only_sd_refined(task):
    """
    Worker for one DM mass point:
      - total SD capture
      - per-label contributions for verified-only labels
      - constant cross section
      - T = 0
    """
    (
        earth_data,
        m,
        sigma_SD_p,
        rho_chi,
        u_max,
        n_u,
        v0,
        n_t_speed,
        n_t_mu,
        n_scatter_mu,
        n_scatter_phi,
        label_order
    ) = task

    print(f"[worker:sd-refined] start m = {m:.4g} GeV", flush=True)

    out = {"m": m}

    # total SD, verified-only, constant, T=0
    c_total = capture_rate_total(
        earth_data=earth_data,
        DM_mass=m,
        sigma_SD_p=sigma_SD_p,
        scattering_type="SD",
        cross_section_type="constant",
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
    out["C_total"] = c_total

    # per-label SD contributions
    for label in label_order:
        earth_data_one = build_sd_single_element_earth_data(earth_data, label)
        c_label = capture_rate_total(
            earth_data=earth_data_one,
            DM_mass=m,
            sigma_SD_p=sigma_SD_p,
            scattering_type="SD",
            cross_section_type="constant",
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
        out[f"C_{label}"] = c_label

    print(f"[worker:sd-refined] done  m = {m:.4g} GeV", flush=True)
    return out


def save_verified_only_sd_refined_results_to_csv(
    DM_masses,
    C_total,
    C_H,
    C_Al,
    C_Si29,
    C_Na,
    output_csv="verified_only_refined_sd_baseline.csv"
):
    """
    Save verified-only refined SD baseline results to CSV.
    """
    DM_masses = np.asarray(DM_masses, dtype=float)
    C_total = np.asarray(C_total, dtype=float)
    C_H = np.asarray(C_H, dtype=float)
    C_Al = np.asarray(C_Al, dtype=float)
    C_Si29 = np.asarray(C_Si29, dtype=float)
    C_Na = np.asarray(C_Na, dtype=float)

    f_H = np.divide(C_H, C_total, out=np.zeros_like(C_H), where=(C_total > 0.0))
    f_Al = np.divide(C_Al, C_total, out=np.zeros_like(C_Al), where=(C_total > 0.0))
    f_Si29 = np.divide(C_Si29, C_total, out=np.zeros_like(C_Si29), where=(C_total > 0.0))
    f_Na = np.divide(C_Na, C_total, out=np.zeros_like(C_Na), where=(C_total > 0.0))

    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "m_chi_GeV",
            "C_total",
            "C_H",
            "C_Al",
            "C_Si29",
            "C_Na",
            "f_H",
            "f_Al",
            "f_Si29",
            "f_Na"
        ])
        for i in range(len(DM_masses)):
            writer.writerow([
                f"{DM_masses[i]:.8e}",
                f"{C_total[i]:.8e}",
                f"{C_H[i]:.8e}",
                f"{C_Al[i]:.8e}",
                f"{C_Si29[i]:.8e}",
                f"{C_Na[i]:.8e}",
                f"{f_H[i]:.8e}",
                f"{f_Al[i]:.8e}",
                f"{f_Si29[i]:.8e}",
                f"{f_Na[i]:.8e}",
            ])

    print(f"[INFO] Verified-only refined SD baseline CSV saved to: {output_csv}")


def print_verified_only_sd_probe_summary(
    DM_masses,
    C_total,
    C_H,
    C_Al,
    C_Si29,
    C_Na,
    probe_masses=(0.7, 2.0, 25.0, 100.0)
):
    """
    Print the same probe-mass summary you already used for interpretation.
    """
    DM_masses = np.asarray(DM_masses, dtype=float)
    C_total = np.asarray(C_total, dtype=float)
    C_H = np.asarray(C_H, dtype=float)
    C_Al = np.asarray(C_Al, dtype=float)
    C_Si29 = np.asarray(C_Si29, dtype=float)
    C_Na = np.asarray(C_Na, dtype=float)

    print("\n" + "=" * 80)
    print("Verified-only refined SD probe-mass summary")
    print("=" * 80)

    for mp in probe_masses:
        idx = np.argmin(np.abs(DM_masses - mp))

        total = C_total[idx]
        pieces = [
            ("H",    C_H[idx]),
            ("Al",   C_Al[idx]),
            ("Si29", C_Si29[idx]),
            ("Na",   C_Na[idx]),
        ]
        pieces = [(lab, val, (val / total if total > 0.0 else np.nan)) for lab, val in pieces]
        pieces.sort(key=lambda x: x[1], reverse=True)

        print(f"\nm_chi = {DM_masses[idx]:.4g} GeV")
        print(f"C_total = {total:.6e} s^-1")
        for lab, val, frac in pieces:
            print(f"  {lab:>4s} : {100.0*frac:8.4f}%   (C = {val:.6e})")


def plot_verified_only_refined_sd_baseline(
    earth_data,
    DM_masses,
    sigma_SD_p=1e-40,
    output_root=".",
    u_max=800.0,
    n_u=40,
    rho_chi=RHO_CHI_DEFAULT,
    v0=V0_DEFAULT,
    n_t_speed=4,
    n_t_mu=4,
    n_scatter_mu=5,
    n_scatter_phi=6,
    max_workers=None
):
    """
    Formal verified-only refined SD baseline:
        - SD only
        - constant cross section
        - T = 0
        - verified_only
        - explicitly decomposed into H, Al, Si29, Na
    """
    import os
    from concurrent.futures import ProcessPoolExecutor, as_completed

    expected_labels = ("Al", "H", "Na", "Si29")
    label_order_for_calc = ("H", "Al", "Si29", "Na")

    validate_verified_only_sd_baseline(
        earth_data=earth_data,
        expected_labels=expected_labels
    )

    if max_workers is None:
        max_workers = max(1, os.cpu_count() - 1)

    DM_masses = np.asarray(DM_masses, dtype=float)

    figures_dir = os.path.join(output_root, "figures")
    results_dir = os.path.join(output_root, "results")
    logs_dir = os.path.join(output_root, "logs")
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    sigma_tag = f"{sigma_SD_p:.0e}".replace("+", "")
    fig_png = os.path.join(
        figures_dir,
        f"verified_only_refined_sd_baseline_sigma{sigma_tag}.png"
    )
    fig_pdf = os.path.join(
        figures_dir,
        f"verified_only_refined_sd_baseline_sigma{sigma_tag}.pdf"
    )
    csv_path = os.path.join(
        results_dir,
        f"verified_only_refined_sd_baseline_sigma{sigma_tag}.csv"
    )
    log_path = os.path.join(
        logs_dir,
        f"verified_only_refined_sd_baseline_sigma{sigma_tag}.txt"
    )

    print("\n" + "=" * 80)
    print("Generating verified-only refined SD baseline")
    print("=" * 80)
    print(f"Using {max_workers} worker processes")

    tasks = []
    for m in DM_masses:
        tasks.append((
            earth_data,
            m,
            sigma_SD_p,
            rho_chi,
            u_max,
            n_u,
            v0,
            n_t_speed,
            n_t_mu,
            n_scatter_mu,
            n_scatter_phi,
            label_order_for_calc
        ))

    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(compute_one_mass_point_verified_only_sd_refined, task) for task in tasks]
        for i, fut in enumerate(as_completed(futures), 1):
            res = fut.result()
            print(f"  done {i}/{len(futures)} : m_chi = {res['m']:.4g} GeV")
            results.append(res)

    results.sort(key=lambda x: x["m"])

    masses_sorted = np.array([r["m"] for r in results], dtype=float)
    C_total = np.array([r["C_total"] for r in results], dtype=float)
    C_H = np.array([r["C_H"] for r in results], dtype=float)
    C_Al = np.array([r["C_Al"] for r in results], dtype=float)
    C_Si29 = np.array([r["C_Si29"] for r in results], dtype=float)
    C_Na = np.array([r["C_Na"] for r in results], dtype=float)

    f_H = np.divide(C_H, C_total, out=np.zeros_like(C_H), where=(C_total > 0.0))
    f_Al = np.divide(C_Al, C_total, out=np.zeros_like(C_Al), where=(C_total > 0.0))
    f_Si29 = np.divide(C_Si29, C_total, out=np.zeros_like(C_Si29), where=(C_total > 0.0))
    f_Na = np.divide(C_Na, C_total, out=np.zeros_like(C_Na), where=(C_total > 0.0))

    save_verified_only_sd_refined_results_to_csv(
        DM_masses=masses_sorted,
        C_total=C_total,
        C_H=C_H,
        C_Al=C_Al,
        C_Si29=C_Si29,
        C_Na=C_Na,
        output_csv=csv_path
    )

    print_verified_only_sd_probe_summary(
        DM_masses=masses_sorted,
        C_total=C_total,
        C_H=C_H,
        C_Al=C_Al,
        C_Si29=C_Si29,
        C_Na=C_Na,
        probe_masses=(0.7, 2.0, 25.0, 100.0)
    )

    # -----------------------------
    # Save a short run-info text
    # -----------------------------
    with open(log_path, "w", encoding="utf-8") as f:
        f.write("Verified-only refined SD baseline run info\n")
        f.write("=" * 60 + "\n")
        f.write("Definition:\n")
        f.write("  - SD only\n")
        f.write("  - constant cross section\n")
        f.write("  - T = 0\n")
        f.write("  - verified_only\n")
        f.write("  - active labels forced to {H, Al, Si29, Na}\n\n")
        f.write(f"sigma_SD_p    = {sigma_SD_p:.6e} cm^2\n")
        f.write(f"rho_chi       = {rho_chi:.6e} GeV/cm^3\n")
        f.write(f"v0            = {v0:.6f} km/s\n")
        f.write(f"u_max         = {u_max:.6f} km/s\n")
        f.write(f"n_u           = {n_u}\n")
        f.write(f"n_t_speed     = {n_t_speed}\n")
        f.write(f"n_t_mu        = {n_t_mu}\n")
        f.write(f"n_scatter_mu  = {n_scatter_mu}\n")
        f.write(f"n_scatter_phi = {n_scatter_phi}\n")
        f.write(f"max_workers   = {max_workers}\n")
        f.write(f"mass range    = {masses_sorted.min():.6e} -> {masses_sorted.max():.6e} GeV\n")
        f.write(f"n_masses      = {len(masses_sorted)}\n\n")
        f.write(f"C_total(min,max) = ({np.nanmin(C_total):.6e}, {np.nanmax(C_total):.6e})\n")
    print(f"[INFO] Run info saved to: {log_path}")

    # -----------------------------
    # Plot
    # -----------------------------
    fig, (ax1, ax2) = plt.subplots(
        2, 1,
        figsize=(8.0, 8.5),
        sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0]}
    )

    # top panel: absolute rates
    ax1.loglog(masses_sorted, C_total, color="black", lw=2.6, label="total")
    ax1.loglog(masses_sorted, C_H,    color="#2ca02c", lw=2.0, label="H")
    ax1.loglog(masses_sorted, C_Al,   color="#d62728", lw=2.0, label="Al")
    ax1.loglog(masses_sorted, C_Si29, color="#1f77b4", lw=2.0, label="Si29")
    ax1.loglog(masses_sorted, C_Na,   color="#ff7f0e", lw=2.0, label="Na")

    ax1.set_ylabel(r"$C\ [{\rm s}^{-1}]$", fontsize=12)
    ax1.set_title(r"Verified-only refined SD baseline: constant, $T=0$", fontsize=13)
    ax1.grid(True, which="both", alpha=0.3)
    ax1.legend(fontsize=10, frameon=False, ncol=2)

    # bottom panel: fractions
    ax2.semilogx(masses_sorted, f_H,    color="#2ca02c", lw=2.0, label="H")
    ax2.semilogx(masses_sorted, f_Al,   color="#d62728", lw=2.0, label="Al")
    ax2.semilogx(masses_sorted, f_Si29, color="#1f77b4", lw=2.0, label="Si29")
    ax2.semilogx(masses_sorted, f_Na,   color="#ff7f0e", lw=2.0, label="Na")

    ax2.set_xlabel(r"$m_\chi\ [{\rm GeV}]$", fontsize=12)
    ax2.set_ylabel("fraction", fontsize=12)
    ax2.set_ylim(0.0, 1.02)
    ax2.grid(True, which="both", alpha=0.3)

    plt.tight_layout()
    fig.savefig(fig_png, dpi=300, bbox_inches="tight")
    fig.savefig(fig_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"[INFO] Figure saved to: {fig_png}")
    print(f"[INFO] Figure saved to: {fig_pdf}")

    return {
        "DM_masses": masses_sorted,
        "C_total": C_total,
        "C_H": C_H,
        "C_Al": C_Al,
        "C_Si29": C_Si29,
        "C_Na": C_Na,
        "f_H": f_H,
        "f_Al": f_Al,
        "f_Si29": f_Si29,
        "f_Na": f_Na,
        "csv_path": csv_path,
        "fig_png": fig_png,
        "fig_pdf": fig_pdf,
        "log_path": log_path,
    }

    save_sd_baseline_run_info(
        output_txt=log_path,
        sigma_SD_p=sigma_SD_p,
        results=results,
        u_max=u_max,
        n_u=n_u,
        n_t_speed=n_t_speed,
        n_t_mu=n_t_mu,
        n_scatter_mu=n_scatter_mu,
        n_scatter_phi=n_scatter_phi,
        rho_chi=rho_chi,
        v0=v0,
        max_workers=max_workers
    )

    return results

def build_sd_single_element_earth_data(earth_data, sd_label):
    """
    Build a lightweight copy of earth_data whose SD active shells keep only one SD label,
    e.g. 'H', 'Al', 'Si29', 'Fe57', etc.
    """
    earth_data_elem = dict(earth_data)
    active_shells_new = []

    for shell in earth_data["active_shells"]:
        shell_new = {
            "SI": [],
            "SD": [],
            "electron": dict(shell["electron"]),
        }

        for info in shell["SD"]:
            if info["elem"] == sd_label:
                shell_new["SD"].append(dict(info))

        active_shells_new.append(shell_new)

    earth_data_elem["active_shells"] = active_shells_new
    earth_data_elem["sd_mode"] = earth_data.get("sd_mode", "include_placeholders")
    return earth_data_elem

def diagnose_sd_element_contributions(
    earth_data,
    DM_mass,
    sigma_SD_p=1e-40,
    cross_section_type="constant",
    rho_chi=RHO_CHI_DEFAULT,
    u_max=800.0,
    n_u=40,
    v0=V0_DEFAULT,
    n_t_speed=4,
    n_t_mu=4,
    n_scatter_mu=5,
    n_scatter_phi=6,
    top_n=10,
    save_csv=False,
    output_csv=None,
    verbose=True
):
    """
    Diagnose which SD isotope labels dominate the capture rate at one DM mass.

    Returns
    -------
    diagnostics : dict
        {
            "DM_mass": float,
            "total_sd_capture": float,
            "contributions": [
                {
                    "elem": str,
                    "C_SD_0": float,
                    "fraction": float or np.nan
                },
                ...
            ]
        }
    """
    import csv
    import numpy as np

    top_n = max(int(top_n), 1)

    # ------------------------------------------------------------
    # Total SD capture rate
    # ------------------------------------------------------------
    C_total = capture_rate_total(
        earth_data=earth_data,
        DM_mass=DM_mass,
        sigma_SD_p=sigma_SD_p,
        scattering_type="SD",
        cross_section_type=cross_section_type,
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
    C_total = float(C_total)

    # ------------------------------------------------------------
    # Collect active SD labels from current earth_data
    # ------------------------------------------------------------
    sd_labels = set()
    for shell in earth_data["active_shells"]:
        for info in shell["SD"]:
            sd_labels.add(info["elem"])

    sd_labels = sorted(sd_labels)

    # ------------------------------------------------------------
    # Compute one-label-at-a-time contributions
    # ------------------------------------------------------------
    contributions = []

    for sd_label in sd_labels:
        earth_data_elem = build_sd_single_element_earth_data(earth_data, sd_label)

        C_elem = capture_rate_total(
            earth_data=earth_data_elem,
            DM_mass=DM_mass,
            sigma_SD_p=sigma_SD_p,
            scattering_type="SD",
            cross_section_type=cross_section_type,
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
        C_elem = float(C_elem)

        if C_total > 0.0:
            frac = float(C_elem / C_total)
        else:
            frac = float("nan")

        contributions.append({
            "elem": str(sd_label),
            "C_SD_0": C_elem,
            "fraction": frac
        })

    contributions.sort(key=lambda x: x["C_SD_0"], reverse=True)

    diagnostics = {
        "DM_mass": float(DM_mass),
        "total_sd_capture": C_total,
        "contributions": contributions
    }

    # ------------------------------------------------------------
    # Verbose printout
    # ------------------------------------------------------------
    if verbose:
        print("\n" + "=" * 80)
        print(f"SD element contribution diagnosis at m_chi = {float(DM_mass):.6g} GeV")
        print("=" * 80)
        print(f"Total SD capture = {C_total:.6e} s^-1")
        print("-" * 80)
        print(f"{'rank':>4s}  {'elem':>8s}  {'C_SD_0 [s^-1]':>18s}  {'fraction':>12s}")
        print("-" * 80)

        n_show = min(top_n, len(contributions))
        for i, item in enumerate(contributions[:n_show], 1):
            c_val = float(item["C_SD_0"])
            frac_val = item["fraction"]

            if np.isfinite(frac_val):
                frac_str = f"{float(frac_val):12.6e}"
            else:
                frac_str = f"{'nan':>12s}"

            print(
                f"{i:4d}  "
                f"{item['elem']:>8s}  "
                f"{c_val:18.6e}  "
                f"{frac_str}"
            )

        frac_sum_top = float(np.nansum([x["fraction"] for x in contributions[:n_show]]))

        print("-" * 80)
        print(f"Top-{n_show} fraction sum = {frac_sum_top:.6e}")

    # ------------------------------------------------------------
    # Optional CSV output
    # ------------------------------------------------------------
    if save_csv:
        if output_csv is None:
            output_csv = f"sd_element_diagnostics_m{float(DM_mass):.3g}GeV.csv"

        with open(output_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["m_chi_GeV", "elem", "C_SD_0", "fraction"])

            for item in contributions:
                frac_val = item["fraction"]
                frac_out = f"{float(frac_val):.8e}" if np.isfinite(frac_val) else "nan"

                writer.writerow([
                    f"{float(DM_mass):.8e}",
                    item["elem"],
                    f"{float(item['C_SD_0']):.8e}",
                    frac_out
                ])

        print(f"[INFO] SD element diagnostic CSV saved to: {output_csv}")

    return diagnostics

def run_sd_diagnostics_for_selected_masses(
    earth_data,
    masses,
    sigma_SD_p=1e-40,
    cross_section_type="constant",
    rho_chi=RHO_CHI_DEFAULT,
    u_max=800.0,
    n_u=40,
    v0=V0_DEFAULT,
    n_t_speed=4,
    n_t_mu=4,
    n_scatter_mu=5,
    n_scatter_phi=6,
    top_n=10,
    save_csv=True,
    output_dir="results/sd_diagnostics"
):
    """
    Run SD element contribution diagnostics for a list of selected masses.
    """
    import os

    if save_csv:
        os.makedirs(output_dir, exist_ok=True)

    all_diag = {}

    for m in masses:
        output_csv = None
        if save_csv:
            output_csv = os.path.join(
                output_dir,
                f"sd_element_diagnostics_m{m:.3g}GeV.csv"
            )

        diag = diagnose_sd_element_contributions(
            earth_data=earth_data,
            DM_mass=m,
            sigma_SD_p=sigma_SD_p,
            cross_section_type=cross_section_type,
            rho_chi=rho_chi,
            u_max=u_max,
            n_u=n_u,
            v0=v0,
            n_t_speed=n_t_speed,
            n_t_mu=n_t_mu,
            n_scatter_mu=n_scatter_mu,
            n_scatter_phi=n_scatter_phi,
            top_n=top_n,
            save_csv=save_csv,
            output_csv=output_csv,
            verbose=True
        )

        all_diag[m] = diag

    return all_diag

def collect_active_sd_labels_from_earth_data(earth_data):
    """
    Collect all active SD labels from earth_data["active_shells"].
    """
    labels = set()
    for shell in earth_data["active_shells"]:
        for info in shell["SD"]:
            labels.add(info["elem"])
    return sorted(labels)


def compute_one_mass_point_sd_mode_total(task):
    """
    Worker for one SD total-capture mass point under a given earth_data mode.
    """
    (
        earth_data,
        m,
        sigma_SD_p,
        cross_section_type,
        rho_chi,
        u_max,
        n_u,
        v0,
        n_t_speed,
        n_t_mu,
        n_scatter_mu,
        n_scatter_phi
    ) = task

    mode_tag = earth_data.get("sd_mode", "unknown")
    print(f"[worker:sd:{mode_tag}] start m = {m:.4g} GeV", flush=True)

    c_sd_0 = capture_rate_total(
        earth_data=earth_data,
        DM_mass=m,
        sigma_SD_p=sigma_SD_p,
        scattering_type="SD",
        cross_section_type=cross_section_type,
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

    print(f"[worker:sd:{mode_tag}] done  m = {m:.4g} GeV", flush=True)

    return {
        "m": float(m),
        "C_SD_0": float(c_sd_0),
    }


def run_sd_total_curve_parallel(
    earth_data,
    DM_masses,
    sigma_SD_p=1e-40,
    cross_section_type="constant",
    rho_chi=RHO_CHI_DEFAULT,
    u_max=800.0,
    n_u=40,
    v0=V0_DEFAULT,
    n_t_speed=4,
    n_t_mu=4,
    n_scatter_mu=5,
    n_scatter_phi=6,
    max_workers=None
):
    """
    Compute one SD total capture curve in parallel for one earth_data mode.
    """
    import os
    from concurrent.futures import ProcessPoolExecutor, as_completed

    if max_workers is None:
        max_workers = max(1, os.cpu_count() - 1)

    DM_masses = np.asarray(DM_masses, dtype=float)

    tasks = []
    for m in DM_masses:
        tasks.append((
            earth_data,
            float(m),
            sigma_SD_p,
            cross_section_type,
            rho_chi,
            u_max,
            n_u,
            v0,
            n_t_speed,
            n_t_mu,
            n_scatter_mu,
            n_scatter_phi
        ))

    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(compute_one_mass_point_sd_mode_total, task) for task in tasks]
        for i, fut in enumerate(as_completed(futures), 1):
            res = fut.result()
            print(
                f"  [{earth_data.get('sd_mode', 'unknown')}] done {i}/{len(futures)} : "
                f"m_chi = {res['m']:.4g} GeV"
            )
            results.append(res)

    results.sort(key=lambda x: x["m"])

    masses_sorted = np.array([r["m"] for r in results], dtype=float)
    C_SD_0 = np.array([r["C_SD_0"] for r in results], dtype=float)

    return {
        "DM_masses": masses_sorted,
        "C_SD_0": C_SD_0
    }


def save_sd_verified_vs_placeholder_comparison_to_csv(
    DM_masses,
    C_verified,
    C_with_placeholders,
    output_csv="sd_verified_vs_include_placeholders_comparison.csv"
):
    """
    Save comparison of:
        verified_only
        include_placeholders
    """
    DM_masses = np.asarray(DM_masses, dtype=float)
    C_verified = np.asarray(C_verified, dtype=float)
    C_with_placeholders = np.asarray(C_with_placeholders, dtype=float)

    if not (len(DM_masses) == len(C_verified) == len(C_with_placeholders)):
        raise ValueError("Input arrays must have the same length.")

    ratio = np.divide(
        C_with_placeholders,
        C_verified,
        out=np.full_like(C_with_placeholders, np.nan, dtype=float),
        where=(C_verified > 0.0)
    )

    frac_excess = np.divide(
        C_with_placeholders - C_verified,
        C_verified,
        out=np.full_like(C_with_placeholders, np.nan, dtype=float),
        where=(C_verified > 0.0)
    )

    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "m_chi_GeV",
            "C_SD_verified_only",
            "C_SD_include_placeholders",
            "ratio_include_over_verified",
            "fractional_excess"
        ])

        for m, cv, cp, r, fe in zip(DM_masses, C_verified, C_with_placeholders, ratio, frac_excess):
            writer.writerow([
                f"{float(m):.8e}",
                f"{float(cv):.8e}",
                f"{float(cp):.8e}",
                f"{float(r):.8e}" if np.isfinite(r) else "nan",
                f"{float(fe):.8e}" if np.isfinite(fe) else "nan",
            ])

    print(f"[INFO] SD comparison CSV saved to: {output_csv}")


def plot_sd_verified_vs_include_placeholders_comparison(
    earth_prem_path="data/earth_prem.dat",
    DM_masses=None,
    sigma_SD_p=1e-40,
    output_root=".",
    u_max=800.0,
    n_u=40,
    rho_chi=RHO_CHI_DEFAULT,
    v0=V0_DEFAULT,
    n_t_speed=4,
    n_t_mu=4,
    n_scatter_mu=5,
    n_scatter_phi=6,
    max_workers=None,
    min_mass_fraction=1e-10
):
    """
    Compare SD total capture curves between:
        1. verified_only
        2. include_placeholders

    Fixed definition:
        - SD only
        - constant cross section
        - T = 0
    """
    import os

    if DM_masses is None:
        DM_masses = np.logspace(np.log10(0.3), np.log10(300.0), 161)

    DM_masses = np.asarray(DM_masses, dtype=float)

    if max_workers is None:
        max_workers = max(1, os.cpu_count() - 1)

    sigma_tag = f"{sigma_SD_p:.0e}".replace("+", "")

    figures_dir = os.path.join(output_root, "figures")
    results_dir = os.path.join(output_root, "results")
    logs_dir = os.path.join(output_root, "logs")
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    fig_png = os.path.join(
        figures_dir,
        f"sd_verified_vs_include_placeholders_comparison_sigma{sigma_tag}.png"
    )
    fig_pdf = os.path.join(
        figures_dir,
        f"sd_verified_vs_include_placeholders_comparison_sigma{sigma_tag}.pdf"
    )
    csv_path = os.path.join(
        results_dir,
        f"sd_verified_vs_include_placeholders_comparison_sigma{sigma_tag}.csv"
    )
    log_path = os.path.join(
        logs_dir,
        f"sd_verified_vs_include_placeholders_comparison_sigma{sigma_tag}.txt"
    )

    print("\n" + "=" * 90)
    print("Loading Earth data for SD comparison")
    print("=" * 90)

    earth_data_verified = load_earth_composition(
        filepath=earth_prem_path,
        sd_mode="verified_only",
        min_mass_fraction=min_mass_fraction
    )

    earth_data_placeholders = load_earth_composition(
        filepath=earth_prem_path,
        sd_mode="include_placeholders",
        min_mass_fraction=min_mass_fraction
    )

    labels_verified = collect_active_sd_labels_from_earth_data(earth_data_verified)
    labels_placeholders = collect_active_sd_labels_from_earth_data(earth_data_placeholders)

    expected_verified = {"Al", "H", "Na", "Si29"}
    if set(labels_verified) != expected_verified:
        raise RuntimeError(
            "verified_only active SD labels do not match expectation.\n"
            f"Expected: {sorted(expected_verified)}\n"
            f"Actual:   {labels_verified}"
        )

    extra_labels = sorted(set(labels_placeholders) - set(labels_verified))

    print("\n[verified_only] active labels:")
    print(labels_verified)
    print(f"Count = {len(labels_verified)}")

    print("\n[include_placeholders] active labels:")
    print(labels_placeholders)
    print(f"Count = {len(labels_placeholders)}")

    print("\nExtra labels contributed only by include_placeholders:")
    print(extra_labels)
    print(f"Count = {len(extra_labels)}")

    print("\n" + "=" * 90)
    print("Computing SD curve: verified_only")
    print("=" * 90)
    res_verified = run_sd_total_curve_parallel(
        earth_data=earth_data_verified,
        DM_masses=DM_masses,
        sigma_SD_p=sigma_SD_p,
        cross_section_type="constant",
        rho_chi=rho_chi,
        u_max=u_max,
        n_u=n_u,
        v0=v0,
        n_t_speed=n_t_speed,
        n_t_mu=n_t_mu,
        n_scatter_mu=n_scatter_mu,
        n_scatter_phi=n_scatter_phi,
        max_workers=max_workers
    )

    print("\n" + "=" * 90)
    print("Computing SD curve: include_placeholders")
    print("=" * 90)
    res_placeholders = run_sd_total_curve_parallel(
        earth_data=earth_data_placeholders,
        DM_masses=DM_masses,
        sigma_SD_p=sigma_SD_p,
        cross_section_type="constant",
        rho_chi=rho_chi,
        u_max=u_max,
        n_u=n_u,
        v0=v0,
        n_t_speed=n_t_speed,
        n_t_mu=n_t_mu,
        n_scatter_mu=n_scatter_mu,
        n_scatter_phi=n_scatter_phi,
        max_workers=max_workers
    )

    masses_verified = np.asarray(res_verified["DM_masses"], dtype=float)
    masses_placeholders = np.asarray(res_placeholders["DM_masses"], dtype=float)

    if len(masses_verified) != len(masses_placeholders) or not np.allclose(masses_verified, masses_placeholders):
        raise RuntimeError("Mass grids do not match between verified_only and include_placeholders runs.")

    masses = masses_verified
    C_verified = np.asarray(res_verified["C_SD_0"], dtype=float)
    C_placeholders = np.asarray(res_placeholders["C_SD_0"], dtype=float)

    ratio = np.divide(
        C_placeholders,
        C_verified,
        out=np.full_like(C_placeholders, np.nan, dtype=float),
        where=(C_verified > 0.0)
    )

    frac_excess = np.divide(
        C_placeholders - C_verified,
        C_verified,
        out=np.full_like(C_placeholders, np.nan, dtype=float),
        where=(C_verified > 0.0)
    )

    save_sd_verified_vs_placeholder_comparison_to_csv(
        DM_masses=masses,
        C_verified=C_verified,
        C_with_placeholders=C_placeholders,
        output_csv=csv_path
    )

    finite_ratio = np.isfinite(ratio)
    if np.any(finite_ratio):
        idx_max = int(np.nanargmax(ratio))
        idx_min = int(np.nanargmin(ratio))
        ratio_max = float(ratio[idx_max])
        ratio_min = float(ratio[idx_min])
        m_ratio_max = float(masses[idx_max])
        m_ratio_min = float(masses[idx_min])
    else:
        ratio_max = float("nan")
        ratio_min = float("nan")
        m_ratio_max = float("nan")
        m_ratio_min = float("nan")

    with open(log_path, "w", encoding="utf-8") as f:
        f.write("SD verified_only vs include_placeholders comparison\n")
        f.write("=" * 70 + "\n")
        f.write("Definition:\n")
        f.write("  - SD only\n")
        f.write("  - constant cross section\n")
        f.write("  - T = 0\n\n")
        f.write(f"sigma_SD_p         = {sigma_SD_p:.6e} cm^2\n")
        f.write(f"rho_chi            = {rho_chi:.6e} GeV/cm^3\n")
        f.write(f"v0                 = {v0:.6f} km/s\n")
        f.write(f"u_max              = {u_max:.6f} km/s\n")
        f.write(f"n_u                = {n_u}\n")
        f.write(f"n_t_speed          = {n_t_speed}\n")
        f.write(f"n_t_mu             = {n_t_mu}\n")
        f.write(f"n_scatter_mu       = {n_scatter_mu}\n")
        f.write(f"n_scatter_phi      = {n_scatter_phi}\n")
        f.write(f"max_workers        = {max_workers}\n")
        f.write(f"earth_prem_path    = {earth_prem_path}\n")
        f.write(f"mass range         = {masses.min():.6e} -> {masses.max():.6e} GeV\n")
        f.write(f"n_masses           = {len(masses)}\n\n")

        f.write("verified_only active labels:\n")
        f.write(f"  {labels_verified}\n\n")

        f.write("include_placeholders active labels:\n")
        f.write(f"  {labels_placeholders}\n\n")

        f.write("extra labels in include_placeholders:\n")
        f.write(f"  {extra_labels}\n\n")

        f.write("Summary:\n")
        f.write(f"max ratio(include/verified) = {ratio_max:.6e} at m_chi = {m_ratio_max:.6e} GeV\n")
        f.write(f"min ratio(include/verified) = {ratio_min:.6e} at m_chi = {m_ratio_min:.6e} GeV\n")
        f.write(f"C_verified(min,max)         = ({np.nanmin(C_verified):.6e}, {np.nanmax(C_verified):.6e})\n")
        f.write(f"C_placeholders(min,max)     = ({np.nanmin(C_placeholders):.6e}, {np.nanmax(C_placeholders):.6e})\n")

    print(f"[INFO] Comparison run info saved to: {log_path}")

    # ------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(
        2, 1,
        figsize=(8.2, 8.6),
        sharex=True,
        gridspec_kw={"height_ratios": [2.1, 1.0]}
    )

    ax1.loglog(
        masses, C_verified,
        color="black", lw=2.7,
        label="verified_only"
    )
    ax1.loglog(
        masses, C_placeholders,
        color="magenta", lw=2.3, ls="--",
        label="include_placeholders"
    )

    ax1.set_ylabel(r"$C\ [{\rm s}^{-1}]$", fontsize=12)
    ax1.set_title(r"SD comparison: verified-only vs include-placeholders, constant, $T=0$", fontsize=13)
    ax1.grid(True, which="both", alpha=0.3)
    ax1.legend(fontsize=10, frameon=False, loc="best")

    ratio_mask = np.isfinite(ratio) & (ratio > 0.0)
    ax2.semilogx(
        masses[ratio_mask],
        ratio[ratio_mask],
        color="magenta", lw=2.3,
        label=r"$C_{\rm include}/C_{\rm verified}$"
    )
    ax2.axhline(1.0, color="black", ls=":", lw=1.2)
    ax2.set_xlabel(r"$m_\chi\ [{\rm GeV}]$", fontsize=12)
    ax2.set_ylabel("ratio", fontsize=12)
    ax2.grid(True, which="both", alpha=0.3)
    ax2.legend(fontsize=10, frameon=False, loc="best")

    plt.tight_layout()
    fig.savefig(fig_png, dpi=300, bbox_inches="tight")
    fig.savefig(fig_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"[INFO] Comparison figure saved to: {fig_png}")
    print(f"[INFO] Comparison figure saved to: {fig_pdf}")

    return {
        "DM_masses": masses,
        "C_verified_only": C_verified,
        "C_include_placeholders": C_placeholders,
        "ratio_include_over_verified": ratio,
        "fractional_excess": frac_excess,
        "verified_labels": labels_verified,
        "placeholder_labels": labels_placeholders,
        "extra_placeholder_labels": extra_labels,
        "csv_path": csv_path,
        "fig_png": fig_png,
        "fig_pdf": fig_pdf,
        "log_path": log_path,
    }

def collect_sd_label_metadata_from_earth_data(earth_data):
    """
    Collect one metadata record per active SD label from earth_data["active_shells"].
    """
    meta = {}
    for shell in earth_data["active_shells"]:
        for info in shell["SD"]:
            label = info["elem"]
            if label not in meta:
                meta[label] = {
                    "elem": label,
                    "base_elem": info.get("base_elem", ""),
                    "A": info.get("A", None),
                    "iso_name": info.get("iso_name", ""),
                    "status": info.get("status", "unknown"),
                    "source": info.get("source", ""),
                    "notes": info.get("notes", ""),
                }
    return meta


def diagnose_sd_element_contributions_detailed(
    earth_data,
    DM_mass,
    sigma_SD_p=1e-40,
    cross_section_type="constant",
    rho_chi=RHO_CHI_DEFAULT,
    u_max=800.0,
    n_u=40,
    v0=V0_DEFAULT,
    n_t_speed=4,
    n_t_mu=4,
    n_scatter_mu=5,
    n_scatter_phi=6,
    top_n=15,
    save_csv=False,
    output_csv=None,
    verbose=True
):
    """
    Detailed SD contributor diagnosis at one DM mass.
    Includes label status (verified / placeholder) and summary fractions.
    """
    import csv
    import numpy as np

    top_n = max(int(top_n), 1)

    label_meta = collect_sd_label_metadata_from_earth_data(earth_data)
    sd_labels = sorted(label_meta.keys())

    C_total = capture_rate_total(
        earth_data=earth_data,
        DM_mass=DM_mass,
        sigma_SD_p=sigma_SD_p,
        scattering_type="SD",
        cross_section_type=cross_section_type,
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
    C_total = float(C_total)

    contributions = []
    for sd_label in sd_labels:
        earth_data_elem = build_sd_single_element_earth_data(earth_data, sd_label)

        C_elem = capture_rate_total(
            earth_data=earth_data_elem,
            DM_mass=DM_mass,
            sigma_SD_p=sigma_SD_p,
            scattering_type="SD",
            cross_section_type=cross_section_type,
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
        C_elem = float(C_elem)

        frac = float(C_elem / C_total) if C_total > 0.0 else float("nan")

        meta = label_meta[sd_label]
        contributions.append({
            "elem": str(sd_label),
            "base_elem": str(meta.get("base_elem", "")),
            "A": meta.get("A", None),
            "iso_name": str(meta.get("iso_name", "")),
            "status": str(meta.get("status", "unknown")),
            "source": str(meta.get("source", "")),
            "notes": str(meta.get("notes", "")),
            "C_SD_0": C_elem,
            "fraction": frac
        })

    contributions.sort(key=lambda x: x["C_SD_0"], reverse=True)

    verified_fraction_total = float(np.nansum([
        x["fraction"] for x in contributions if x["status"] == "verified"
    ]))
    placeholder_fraction_total = float(np.nansum([
        x["fraction"] for x in contributions if x["status"] == "placeholder"
    ]))

    n_show = min(top_n, len(contributions))
    verified_fraction_top = float(np.nansum([
        x["fraction"] for x in contributions[:n_show] if x["status"] == "verified"
    ]))
    placeholder_fraction_top = float(np.nansum([
        x["fraction"] for x in contributions[:n_show] if x["status"] == "placeholder"
    ]))
    top_fraction_sum = float(np.nansum([
        x["fraction"] for x in contributions[:n_show]
    ]))

    diagnostics = {
        "DM_mass": float(DM_mass),
        "total_sd_capture": C_total,
        "verified_fraction_total": verified_fraction_total,
        "placeholder_fraction_total": placeholder_fraction_total,
        "verified_fraction_top": verified_fraction_top,
        "placeholder_fraction_top": placeholder_fraction_top,
        "top_fraction_sum": top_fraction_sum,
        "contributions": contributions
    }

    if verbose:
        print("\n" + "=" * 100)
        print(f"Detailed SD contributor diagnosis at m_chi = {float(DM_mass):.6g} GeV")
        print(f"sd_mode = {earth_data.get('sd_mode', 'unknown')}")
        print("=" * 100)
        print(f"Total SD capture = {C_total:.6e} s^-1")
        print("-" * 100)
        print(
            f"{'rank':>4s}  "
            f"{'elem':>8s}  "
            f"{'status':>12s}  "
            f"{'C_SD_0 [s^-1]':>18s}  "
            f"{'fraction':>12s}"
        )
        print("-" * 100)

        for i, item in enumerate(contributions[:n_show], 1):
            frac_val = item["fraction"]
            frac_str = f"{float(frac_val):12.6e}" if np.isfinite(frac_val) else f"{'nan':>12s}"
            print(
                f"{i:4d}  "
                f"{item['elem']:>8s}  "
                f"{item['status']:>12s}  "
                f"{float(item['C_SD_0']):18.6e}  "
                f"{frac_str}"
            )

        print("-" * 100)
        print(f"Top-{n_show} fraction sum            = {top_fraction_sum:.6e}")
        print(f"Top-{n_show} verified fraction       = {verified_fraction_top:.6e}")
        print(f"Top-{n_show} placeholder fraction    = {placeholder_fraction_top:.6e}")
        print(f"All-label verified fraction          = {verified_fraction_total:.6e}")
        print(f"All-label placeholder fraction       = {placeholder_fraction_total:.6e}")

        top_placeholder = [x for x in contributions[:n_show] if x["status"] == "placeholder"]
        if top_placeholder:
            print("-" * 100)
            print("Top placeholder contributors within displayed ranks:")
            for item in top_placeholder:
                print(
                    f"  {item['elem']:>8s} : "
                    f"fraction = {float(item['fraction']):.6e}, "
                    f"C = {float(item['C_SD_0']):.6e}"
                )

    if save_csv:
        if output_csv is None:
            output_csv = f"sd_element_diagnostics_detailed_m{float(DM_mass):.3g}GeV.csv"

        with open(output_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow([
                "m_chi_GeV",
                "elem",
                "base_elem",
                "A",
                "iso_name",
                "status",
                "C_SD_0",
                "fraction",
                "source",
                "notes"
            ])

            for item in contributions:
                writer.writerow([
                    f"{float(DM_mass):.8e}",
                    item["elem"],
                    item["base_elem"],
                    item["A"],
                    item["iso_name"],
                    item["status"],
                    f"{float(item['C_SD_0']):.8e}",
                    f"{float(item['fraction']):.8e}" if np.isfinite(item["fraction"]) else "nan",
                    item["source"],
                    item["notes"]
                ])

        print(f"[INFO] Detailed SD diagnostic CSV saved to: {output_csv}")

    return diagnostics


def run_include_placeholders_top_contributor_diagnostics(
    earth_prem_path="data/earth_prem.dat",
    masses=(5.0, 25.0, 45.0, 70.0, 100.0, 300.0),
    sigma_SD_p=1e-40,
    cross_section_type="constant",
    rho_chi=RHO_CHI_DEFAULT,
    u_max=800.0,
    n_u=40,
    v0=V0_DEFAULT,
    n_t_speed=4,
    n_t_mu=4,
    n_scatter_mu=5,
    n_scatter_phi=6,
    top_n=15,
    save_csv=True,
    output_dir="results/sd_placeholder_diagnostics",
    summary_txt="logs/sd_placeholder_diagnostics_summary.txt",
    min_mass_fraction=1e-10
):
    """
    Run detailed contributor diagnostics in sd_mode='include_placeholders'
    for selected DM masses.
    """
    import os
    import numpy as np

    os.makedirs(output_dir, exist_ok=True)

    summary_dir = os.path.dirname(summary_txt)
    if summary_dir:
        os.makedirs(summary_dir, exist_ok=True)

    earth_data = load_earth_composition(
        filepath=earth_prem_path,
        sd_mode="include_placeholders",
        min_mass_fraction=min_mass_fraction
    )

    active_labels = collect_active_sd_labels_from_earth_data(earth_data)
    label_meta = collect_sd_label_metadata_from_earth_data(earth_data)

    verified_labels = sorted([
        lab for lab, meta in label_meta.items()
        if meta.get("status", "unknown") == "verified"
    ])
    placeholder_labels = sorted([
        lab for lab, meta in label_meta.items()
        if meta.get("status", "unknown") == "placeholder"
    ])

    print("\n" + "=" * 100)
    print("Running include_placeholders top-contributor SD diagnostics")
    print("=" * 100)
    print("Active labels:")
    print(active_labels)
    print(f"Count = {len(active_labels)}")
    print("\nVerified labels:")
    print(verified_labels)
    print(f"Count = {len(verified_labels)}")
    print("\nPlaceholder labels:")
    print(placeholder_labels)
    print(f"Count = {len(placeholder_labels)}")

    all_diag = {}

    for m in masses:
        output_csv = None
        if save_csv:
            output_csv = os.path.join(
                output_dir,
                f"sd_include_placeholders_detailed_m{float(m):.3g}GeV.csv"
            )

        diag = diagnose_sd_element_contributions_detailed(
            earth_data=earth_data,
            DM_mass=float(m),
            sigma_SD_p=sigma_SD_p,
            cross_section_type=cross_section_type,
            rho_chi=rho_chi,
            u_max=u_max,
            n_u=n_u,
            v0=v0,
            n_t_speed=n_t_speed,
            n_t_mu=n_t_mu,
            n_scatter_mu=n_scatter_mu,
            n_scatter_phi=n_scatter_phi,
            top_n=top_n,
            save_csv=save_csv,
            output_csv=output_csv,
            verbose=True
        )

        all_diag[float(m)] = diag

    # ------------------------------------------------------------
    # Save compact text summary
    # ------------------------------------------------------------
    with open(summary_txt, "w", encoding="utf-8") as f:
        f.write("include_placeholders SD top-contributor diagnostics summary\n")
        f.write("=" * 80 + "\n")
        f.write(f"earth_prem_path    = {earth_prem_path}\n")
        f.write(f"sigma_SD_p         = {sigma_SD_p:.6e} cm^2\n")
        f.write(f"cross_section_type = {cross_section_type}\n")
        f.write(f"rho_chi            = {rho_chi:.6e} GeV/cm^3\n")
        f.write(f"v0                 = {v0:.6f} km/s\n")
        f.write(f"u_max              = {u_max:.6f} km/s\n")
        f.write(f"n_u                = {n_u}\n")
        f.write(f"n_t_speed          = {n_t_speed}\n")
        f.write(f"n_t_mu             = {n_t_mu}\n")
        f.write(f"n_scatter_mu       = {n_scatter_mu}\n")
        f.write(f"n_scatter_phi      = {n_scatter_phi}\n")
        f.write(f"top_n              = {top_n}\n\n")

        f.write("Active labels:\n")
        f.write(f"{active_labels}\n\n")

        f.write("Verified labels:\n")
        f.write(f"{verified_labels}\n\n")

        f.write("Placeholder labels:\n")
        f.write(f"{placeholder_labels}\n\n")

        for m in masses:
            diag = all_diag[float(m)]
            f.write("-" * 80 + "\n")
            f.write(f"m_chi = {float(m):.6g} GeV\n")
            f.write(f"C_total = {float(diag['total_sd_capture']):.6e} s^-1\n")
            f.write(f"verified_fraction_total    = {float(diag['verified_fraction_total']):.6e}\n")
            f.write(f"placeholder_fraction_total = {float(diag['placeholder_fraction_total']):.6e}\n")
            f.write(f"top_fraction_sum           = {float(diag['top_fraction_sum']):.6e}\n")
            f.write(f"verified_fraction_top      = {float(diag['verified_fraction_top']):.6e}\n")
            f.write(f"placeholder_fraction_top   = {float(diag['placeholder_fraction_top']):.6e}\n")

            f.write("Top contributors:\n")
            for item in diag["contributions"][:top_n]:
                f.write(
                    f"  {item['elem']:>8s} | "
                    f"{item['status']:>11s} | "
                    f"C = {float(item['C_SD_0']):.6e} | "
                    f"f = {float(item['fraction']):.6e}\n"
                )
            f.write("\n")

    print(f"[INFO] Placeholder diagnostic summary saved to: {summary_txt}")

    return {
        "earth_data": earth_data,
        "active_labels": active_labels,
        "verified_labels": verified_labels,
        "placeholder_labels": placeholder_labels,
        "diagnostics": all_diag,
        "summary_txt": summary_txt,
        "output_dir": output_dir,
    }

def save_sd_placeholder_culprit_summary_csv(
    sd_placeholder_diag,
    output_csv="results/sd_placeholder_diagnostics/sd_placeholder_culprit_summary.csv",
    top_k_overall=3,
    top_k_placeholder=3,
    top_k_verified=3
):
    """
    Save a compact cross-mass summary table for include_placeholders SD diagnostics.

    Parameters
    ----------
    sd_placeholder_diag : dict
        Output from run_include_placeholders_top_contributor_diagnostics(...)
    output_csv : str
        Output CSV path
    top_k_overall : int
        Number of top overall contributors to store
    top_k_placeholder : int
        Number of top placeholder contributors to store
    top_k_verified : int
        Number of top verified contributors to store
    """
    import os
    import csv
    import numpy as np

    top_k_overall = max(int(top_k_overall), 1)
    top_k_placeholder = max(int(top_k_placeholder), 1)
    top_k_verified = max(int(top_k_verified), 1)

    diagnostics = sd_placeholder_diag["diagnostics"]

    masses = sorted(float(m) for m in diagnostics.keys())

    out_dir = os.path.dirname(output_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    header = [
        "m_chi_GeV",
        "C_total_include_placeholders",
        "verified_fraction_total",
        "placeholder_fraction_total",
        "estimated_include_over_verified",
    ]

    for i in range(1, top_k_overall + 1):
        header.extend([
            f"top{i}_elem",
            f"top{i}_status",
            f"top{i}_fraction",
            f"top{i}_C_SD_0",
        ])

    for i in range(1, top_k_placeholder + 1):
        header.extend([
            f"top_placeholder_{i}_elem",
            f"top_placeholder_{i}_fraction",
            f"top_placeholder_{i}_C_SD_0",
        ])

    for i in range(1, top_k_verified + 1):
        header.extend([
            f"top_verified_{i}_elem",
            f"top_verified_{i}_fraction",
            f"top_verified_{i}_C_SD_0",
        ])

    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)

        for m in masses:
            diag = diagnostics[m]
            contributions = diag["contributions"]

            C_total = float(diag["total_sd_capture"])
            verified_fraction_total = float(diag["verified_fraction_total"])
            placeholder_fraction_total = float(diag["placeholder_fraction_total"])

            est_ratio = (
                float(1.0 / verified_fraction_total)
                if verified_fraction_total > 0.0 else float("nan")
            )

            overall = contributions[:top_k_overall]
            placeholder_only = [x for x in contributions if x["status"] == "placeholder"][:top_k_placeholder]
            verified_only = [x for x in contributions if x["status"] == "verified"][:top_k_verified]

            row = [
                f"{m:.8e}",
                f"{C_total:.8e}",
                f"{verified_fraction_total:.8e}",
                f"{placeholder_fraction_total:.8e}",
                f"{est_ratio:.8e}" if np.isfinite(est_ratio) else "nan",
            ]

            for i in range(top_k_overall):
                if i < len(overall):
                    item = overall[i]
                    row.extend([
                        item["elem"],
                        item["status"],
                        f"{float(item['fraction']):.8e}",
                        f"{float(item['C_SD_0']):.8e}",
                    ])
                else:
                    row.extend(["", "", "", ""])

            for i in range(top_k_placeholder):
                if i < len(placeholder_only):
                    item = placeholder_only[i]
                    row.extend([
                        item["elem"],
                        f"{float(item['fraction']):.8e}",
                        f"{float(item['C_SD_0']):.8e}",
                    ])
                else:
                    row.extend(["", "", ""])

            for i in range(top_k_verified):
                if i < len(verified_only):
                    item = verified_only[i]
                    row.extend([
                        item["elem"],
                        f"{float(item['fraction']):.8e}",
                        f"{float(item['C_SD_0']):.8e}",
                    ])
                else:
                    row.extend(["", "", ""])

            writer.writerow(row)

    print(f"[INFO] Placeholder culprit summary CSV saved to: {output_csv}")


def print_sd_placeholder_culprit_summary(
    sd_placeholder_diag,
    top_k_overall=3,
    top_k_placeholder=3
):
    """
    Print a compact console summary from run_include_placeholders_top_contributor_diagnostics(...).
    """
    import numpy as np

    diagnostics = sd_placeholder_diag["diagnostics"]
    masses = sorted(float(m) for m in diagnostics.keys())

    print("\n" + "=" * 110)
    print("Compact SD placeholder culprit summary")
    print("=" * 110)

    for m in masses:
        diag = diagnostics[m]
        contributions = diag["contributions"]

        verified_fraction_total = float(diag["verified_fraction_total"])
        placeholder_fraction_total = float(diag["placeholder_fraction_total"])
        est_ratio = (
            float(1.0 / verified_fraction_total)
            if verified_fraction_total > 0.0 else float("nan")
        )

        overall = contributions[:top_k_overall]
        placeholder_only = [x for x in contributions if x["status"] == "placeholder"][:top_k_placeholder]

        print(f"\nm_chi = {m:.6g} GeV")
        print(f"  verified_fraction_total    = {verified_fraction_total:.6e}")
        print(f"  placeholder_fraction_total = {placeholder_fraction_total:.6e}")
        print(f"  estimated include/verified = {est_ratio:.6e}" if np.isfinite(est_ratio) else
              "  estimated include/verified = nan")

        print("  top overall contributors:")
        for item in overall:
            print(
                f"    {item['elem']:>8s} | "
                f"{item['status']:>11s} | "
                f"f = {float(item['fraction']):.6e} | "
                f"C = {float(item['C_SD_0']):.6e}"
            )

        print("  top placeholder contributors:")
        for item in placeholder_only:
            print(
                f"    {item['elem']:>8s} | "
                f"f = {float(item['fraction']):.6e} | "
                f"C = {float(item['C_SD_0']):.6e}"
            )

def compute_one_mass_point_verified_only_sd_thermal(task):
    """
    Worker for one mass point:
        verified-only
        SD only
        constant cross section
        finite-T and T=0
    """
    (
        earth_data,
        m,
        sigma_SD_p,
        rho_chi,
        u_max,
        n_u,
        v0,
        n_t_speed,
        n_t_mu,
        n_scatter_mu,
        n_scatter_phi
    ) = task

    print(f"[worker:sd-thermal] start m = {m:.4g} GeV", flush=True)

    c_sd_T = capture_rate_total(
        earth_data=earth_data,
        DM_mass=m,
        sigma_SD_p=sigma_SD_p,
        scattering_type="SD",
        cross_section_type="constant",
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
        cross_section_type="constant",
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

    print(f"[worker:sd-thermal] done  m = {m:.4g} GeV", flush=True)

    return {
        "m": float(m),
        "C_SD_T": float(c_sd_T),
        "C_SD_0": float(c_sd_0),
    }


def save_verified_only_sd_thermal_results_to_csv(
    DM_masses,
    C_SD_T,
    C_SD_0,
    output_csv="results/verified_only_sd_thermal_constant_sigma1e-40.csv"
):
    """
    Save verified-only SD thermal results to CSV.
    Columns:
        m_chi_GeV
        C_SD_T
        C_SD_0
        ratio_SD
    """
    DM_masses = np.asarray(DM_masses, dtype=float)
    C_SD_T = np.asarray(C_SD_T, dtype=float)
    C_SD_0 = np.asarray(C_SD_0, dtype=float)

    if not (len(DM_masses) == len(C_SD_T) == len(C_SD_0)):
        raise ValueError("DM_masses, C_SD_T, and C_SD_0 must have the same length.")

    ratio_SD = np.divide(
        C_SD_T,
        C_SD_0,
        out=np.full_like(C_SD_T, np.nan, dtype=float),
        where=(C_SD_0 > 0.0)
    )

    out_dir = os.path.dirname(output_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["m_chi_GeV", "C_SD_T", "C_SD_0", "ratio_SD"])
        for m, cT, c0, r in zip(DM_masses, C_SD_T, C_SD_0, ratio_SD):
            writer.writerow([
                f"{float(m):.8e}",
                f"{float(cT):.8e}",
                f"{float(c0):.8e}",
                f"{float(r):.8e}" if np.isfinite(r) else "nan"
            ])

    print(f"[INFO] Verified-only SD thermal CSV saved to: {output_csv}")


def plot_verified_only_sd_thermal_constant(
    earth_data,
    DM_masses,
    sigma_SD_p=1e-40,
    output_root=".",
    u_max=800.0,
    n_u=40,
    rho_chi=RHO_CHI_DEFAULT,
    v0=V0_DEFAULT,
    n_t_speed=4,
    n_t_mu=4,
    n_scatter_mu=5,
    n_scatter_phi=6,
    max_workers=None
):
    """
    Verified-only SD thermal comparison:
        - SD only
        - constant cross section
        - verified_only
        - compare finite-T vs T=0
    """
    import os
    from concurrent.futures import ProcessPoolExecutor, as_completed

    expected_labels = ("Al", "H", "Na", "Si29")
    validate_verified_only_sd_baseline(
        earth_data=earth_data,
        expected_labels=expected_labels
    )

    if max_workers is None:
        max_workers = max(1, os.cpu_count() - 1)

    DM_masses = np.asarray(DM_masses, dtype=float)

    sigma_tag = f"{sigma_SD_p:.0e}".replace("+", "")

    figures_dir = os.path.join(output_root, "figures")
    results_dir = os.path.join(output_root, "results")
    logs_dir = os.path.join(output_root, "logs")
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    fig_png = os.path.join(
        figures_dir,
        f"verified_only_sd_thermal_constant_sigma{sigma_tag}.png"
    )
    fig_pdf = os.path.join(
        figures_dir,
        f"verified_only_sd_thermal_constant_sigma{sigma_tag}.pdf"
    )
    csv_path = os.path.join(
        results_dir,
        f"verified_only_sd_thermal_constant_sigma{sigma_tag}.csv"
    )
    log_path = os.path.join(
        logs_dir,
        f"verified_only_sd_thermal_constant_sigma{sigma_tag}.txt"
    )

    print("\n" + "=" * 90)
    print("Generating verified-only SD thermal comparison")
    print("=" * 90)
    print(f"Using {max_workers} worker processes")

    tasks = []
    for m in DM_masses:
        tasks.append((
            earth_data,
            float(m),
            sigma_SD_p,
            rho_chi,
            u_max,
            n_u,
            v0,
            n_t_speed,
            n_t_mu,
            n_scatter_mu,
            n_scatter_phi
        ))

    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(compute_one_mass_point_verified_only_sd_thermal, task) for task in tasks]
        for i, fut in enumerate(as_completed(futures), 1):
            res = fut.result()
            print(f"  done {i}/{len(futures)} : m_chi = {res['m']:.4g} GeV")
            results.append(res)

    results.sort(key=lambda x: x["m"])

    masses_sorted = np.array([r["m"] for r in results], dtype=float)
    C_SD_T = np.array([r["C_SD_T"] for r in results], dtype=float)
    C_SD_0 = np.array([r["C_SD_0"] for r in results], dtype=float)

    ratio_SD_raw = np.divide(
        C_SD_T,
        C_SD_0,
        out=np.full_like(C_SD_T, np.nan, dtype=float),
        where=(C_SD_0 > 0.0)
    )

    # Plot mask: suppress meaningless huge ratios if denominator is numerically tiny
    sd_floor = max(1e-6 * float(np.nanmax(C_SD_0)), 1e-300)
    ratio_SD_plot = np.where(C_SD_0 > sd_floor, ratio_SD_raw, np.nan)

    save_verified_only_sd_thermal_results_to_csv(
        DM_masses=masses_sorted,
        C_SD_T=C_SD_T,
        C_SD_0=C_SD_0,
        output_csv=csv_path
    )

    finite_mask = np.isfinite(ratio_SD_plot)
    if np.any(finite_mask):
        idx_max = int(np.nanargmax(ratio_SD_plot))
        idx_min = int(np.nanargmin(ratio_SD_plot))
        max_ratio = float(ratio_SD_plot[idx_max])
        min_ratio = float(ratio_SD_plot[idx_min])
        max_ratio_mass = float(masses_sorted[idx_max])
        min_ratio_mass = float(masses_sorted[idx_min])
    else:
        max_ratio = float("nan")
        min_ratio = float("nan")
        max_ratio_mass = float("nan")
        min_ratio_mass = float("nan")

    with open(log_path, "w", encoding="utf-8") as f:
        f.write("Verified-only SD thermal comparison run info\n")
        f.write("=" * 70 + "\n")
        f.write("Definition:\n")
        f.write("  - SD only\n")
        f.write("  - constant cross section\n")
        f.write("  - verified_only\n")
        f.write("  - compare T(r)!=0 vs T=0\n")
        f.write("  - validated active labels = {H, Al, Si29, Na}\n\n")
        f.write(f"sigma_SD_p         = {sigma_SD_p:.6e} cm^2\n")
        f.write(f"rho_chi            = {rho_chi:.6e} GeV/cm^3\n")
        f.write(f"v0                 = {v0:.6f} km/s\n")
        f.write(f"u_max              = {u_max:.6f} km/s\n")
        f.write(f"n_u                = {n_u}\n")
        f.write(f"n_t_speed          = {n_t_speed}\n")
        f.write(f"n_t_mu             = {n_t_mu}\n")
        f.write(f"n_scatter_mu       = {n_scatter_mu}\n")
        f.write(f"n_scatter_phi      = {n_scatter_phi}\n")
        f.write(f"max_workers        = {max_workers}\n")
        f.write(f"mass range         = {masses_sorted.min():.6e} -> {masses_sorted.max():.6e} GeV\n")
        f.write(f"n_masses           = {len(masses_sorted)}\n\n")
        f.write("Summary:\n")
        f.write(f"C_SD_T(min,max)    = ({np.nanmin(C_SD_T):.6e}, {np.nanmax(C_SD_T):.6e})\n")
        f.write(f"C_SD_0(min,max)    = ({np.nanmin(C_SD_0):.6e}, {np.nanmax(C_SD_0):.6e})\n")
        f.write(f"max ratio C_T/C_0  = {max_ratio:.6e} at m_chi = {max_ratio_mass:.6e} GeV\n")
        f.write(f"min ratio C_T/C_0  = {min_ratio:.6e} at m_chi = {min_ratio_mass:.6e} GeV\n")

    print(f"[INFO] Verified-only SD thermal run info saved to: {log_path}")

    # ------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    axL = axes[0]
    axL.loglog(
        masses_sorted,
        C_SD_T,
        color="purple",
        linestyle="-",
        linewidth=2.5,
        label="SD, verified_only, $T(r)\\neq 0$"
    )
    axL.loglog(
        masses_sorted,
        C_SD_0,
        color="purple",
        alpha=0.35,
        linestyle=":",
        linewidth=2.2,
        label="SD, verified_only, $T=0$"
    )
    axL.set_xlabel(r"$m_\chi\ [{\rm GeV}]$", fontsize=12)
    axL.set_ylabel(r"$C\ [{\rm s}^{-1}]$", fontsize=12)
    axL.set_title(r"Verified-only SD, constant", fontsize=13)
    axL.grid(True, which="both", alpha=0.3)
    axL.legend(fontsize=10, loc="best")

    axR = axes[1]
    mask_plot = np.isfinite(ratio_SD_plot)
    axR.semilogx(
        masses_sorted[mask_plot],
        ratio_SD_plot[mask_plot],
        color="purple",
        linestyle="-",
        linewidth=2.3,
        label="SD"
    )
    axR.axhline(
        1.0,
        color="black",
        linestyle=":",
        linewidth=1.2,
        label=r"$T=0$ limit"
    )
    axR.set_xlabel(r"$m_\chi\ [{\rm GeV}]$", fontsize=12)
    axR.set_ylabel(r"$C(T)/C(T=0)$", fontsize=12)
    axR.set_title(r"Verified-only SD thermal correction", fontsize=13)
    axR.set_yscale("log")
    axR.grid(True, which="both", alpha=0.3)
    axR.legend(fontsize=10, loc="best")

    plt.tight_layout()
    fig.savefig(fig_png, dpi=300, bbox_inches="tight")
    fig.savefig(fig_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"[INFO] Verified-only SD thermal figure saved to: {fig_png}")
    print(f"[INFO] Verified-only SD thermal figure saved to: {fig_pdf}")

    return {
        "DM_masses": masses_sorted,
        "C_SD_T": C_SD_T,
        "C_SD_0": C_SD_0,
        "ratio_SD": ratio_SD_raw,
        "csv_path": csv_path,
        "fig_png": fig_png,
        "fig_pdf": fig_pdf,
        "log_path": log_path,
    }

def test_verified_only_sd_thermal_convergence(
    earth_data,
    masses,
    sigma_SD_p=1e-40,
    rho_chi=RHO_CHI_DEFAULT,
    v0=V0_DEFAULT,
    u_max=800.0,
    top_note=True
):
    """
    One-parameter-at-a-time convergence test for verified-only SD thermal capture.

    Purpose
    -------
    Diagnose whether oscillations in
        C_SD(T) / C_SD(T=0)
    are driven by numerical resolution choices.

    Parameters scanned
    ------------------
    - n_u
    - n_t_speed
    - n_t_mu
    - n_scatter_mu
    - n_scatter_phi
    - shell_step

    Notes
    -----
    This assumes:
        - earth_data already loaded with sd_mode="verified_only"
        - SD only
        - constant cross section
    """
    import time
    import numpy as np

    expected_labels = ("Al", "H", "Na", "Si29")
    validate_verified_only_sd_baseline(
        earth_data=earth_data,
        expected_labels=expected_labels
    )

    masses = np.asarray(masses, dtype=float)

    if top_note:
        print("\n" + "=" * 100)
        print("Verified-only SD thermal convergence test")
        print("=" * 100)
        print("Mode assumptions:")
        print("  - SD only")
        print("  - constant cross section")
        print("  - verified_only")
        print("  - compare T(r)!=0 vs T=0")
        print("=" * 100)

    # ------------------------------------------------------------
    # Baseline settings
    # ------------------------------------------------------------
    base = {
        "n_u": 40,
        "n_t_speed": 4,
        "n_t_mu": 4,
        "n_scatter_mu": 5,
        "n_scatter_phi": 6,
        "shell_step": 1,
        "u_grid_mode": "log",
    }

    # ------------------------------------------------------------
    # One-parameter scans
    # ------------------------------------------------------------
    scans = {
        "n_u": [30, 40, 60, 80],
        "n_t_speed": [3, 4, 5, 6],
        "n_t_mu": [3, 4, 5, 6],
        "n_scatter_mu": [4, 5, 6, 8],
        "n_scatter_phi": [4, 6, 8, 12],
        "shell_step": [2, 1],
    }

    for m in masses:
        print("\n" + "-" * 100)
        print(f"Verified-only SD thermal convergence at m_chi = {float(m):.6g} GeV")
        print("-" * 100)

        for param_name, values in scans.items():
            print(f"\nScan: {param_name}")

            for val in values:
                s = dict(base)
                s[param_name] = val

                print(
                    f"[start] m={float(m):.6g}, "
                    f"scan={param_name}, val={val}, "
                    f"settings=("
                    f"n_u={s['n_u']}, "
                    f"n_t_speed={s['n_t_speed']}, "
                    f"n_t_mu={s['n_t_mu']}, "
                    f"n_scatter_mu={s['n_scatter_mu']}, "
                    f"n_scatter_phi={s['n_scatter_phi']}, "
                    f"shell_step={s['shell_step']}, "
                    f"u_grid_mode={s['u_grid_mode']}"
                    f")",
                    flush=True
                )

                t0 = time.time()

                c_sd_T = capture_rate_total(
                    earth_data=earth_data,
                    DM_mass=float(m),
                    sigma_SD_p=sigma_SD_p,
                    scattering_type="SD",
                    cross_section_type="constant",
                    rho_chi=rho_chi,
                    u_max=u_max,
                    n_u=s["n_u"],
                    v0=v0,
                    include_thermal_targets=True,
                    n_t_speed=s["n_t_speed"],
                    n_t_mu=s["n_t_mu"],
                    n_scatter_mu=s["n_scatter_mu"],
                    n_scatter_phi=s["n_scatter_phi"],
                    u_grid_mode=s["u_grid_mode"],
                    shell_step=s["shell_step"]
                )

                c_sd_0 = capture_rate_total(
                    earth_data=earth_data,
                    DM_mass=float(m),
                    sigma_SD_p=sigma_SD_p,
                    scattering_type="SD",
                    cross_section_type="constant",
                    rho_chi=rho_chi,
                    u_max=u_max,
                    n_u=s["n_u"],
                    v0=v0,
                    include_thermal_targets=False,
                    n_t_speed=s["n_t_speed"],
                    n_t_mu=s["n_t_mu"],
                    n_scatter_mu=s["n_scatter_mu"],
                    n_scatter_phi=s["n_scatter_phi"],
                    u_grid_mode=s["u_grid_mode"],
                    shell_step=s["shell_step"]
                )

                dt = time.time() - t0

                c_sd_T = float(c_sd_T)
                c_sd_0 = float(c_sd_0)
                ratio_sd = float(c_sd_T / c_sd_0) if c_sd_0 > 0.0 else float("nan")

                print(
                    f"[done ] {param_name}={val!s:>3s} | "
                    f"C_SD_T={c_sd_T:.6e}, "
                    f"C_SD_0={c_sd_0:.6e}, "
                    f"ratio={ratio_sd:.6e} | "
                    f"time={dt:.1f}s",
                    flush=True
                )

def plot_verified_only_sd_thermal_constant_refined(
    earth_data,
    DM_masses,
    sigma_SD_p=1e-40,
    output_root=".",
    u_max=800.0,
    n_u=80,
    rho_chi=RHO_CHI_DEFAULT,
    v0=V0_DEFAULT,
    n_t_speed=4,
    n_t_mu=4,
    n_scatter_mu=8,
    n_scatter_phi=12,
    max_workers=None
):
    """
    Refined verified-only SD thermal comparison:
        - SD only
        - constant cross section
        - verified_only
        - compare T(r)!=0 vs T=0
        - refined numerical settings motivated by convergence study

    Recommended defaults:
        n_u = 80
        n_t_speed = 4
        n_t_mu = 4
        n_scatter_mu = 8
        n_scatter_phi = 12
        shell_step is kept fixed at 1 inside capture_rate_total
    """
    import os
    from concurrent.futures import ProcessPoolExecutor, as_completed

    expected_labels = ("Al", "H", "Na", "Si29")
    validate_verified_only_sd_baseline(
        earth_data=earth_data,
        expected_labels=expected_labels
    )

    if max_workers is None:
        max_workers = max(1, os.cpu_count() - 1)

    DM_masses = np.asarray(DM_masses, dtype=float)
    sigma_tag = f"{sigma_SD_p:.0e}".replace("+", "")

    figures_dir = os.path.join(output_root, "figures")
    results_dir = os.path.join(output_root, "results")
    logs_dir = os.path.join(output_root, "logs")
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    fig_png = os.path.join(
        figures_dir,
        f"verified_only_sd_thermal_constant_refined_sigma{sigma_tag}.png"
    )
    fig_pdf = os.path.join(
        figures_dir,
        f"verified_only_sd_thermal_constant_refined_sigma{sigma_tag}.pdf"
    )
    csv_path = os.path.join(
        results_dir,
        f"verified_only_sd_thermal_constant_refined_sigma{sigma_tag}.csv"
    )
    log_path = os.path.join(
        logs_dir,
        f"verified_only_sd_thermal_constant_refined_sigma{sigma_tag}.txt"
    )

    print("\n" + "=" * 90)
    print("Generating refined verified-only SD thermal comparison")
    print("=" * 90)
    print(f"Using {max_workers} worker processes")
    print(
        f"Refined settings: "
        f"n_u={n_u}, n_t_speed={n_t_speed}, n_t_mu={n_t_mu}, "
        f"n_scatter_mu={n_scatter_mu}, n_scatter_phi={n_scatter_phi}, shell_step=1"
    )

    tasks = []
    for m in DM_masses:
        tasks.append((
            earth_data,
            float(m),
            sigma_SD_p,
            rho_chi,
            u_max,
            n_u,
            v0,
            n_t_speed,
            n_t_mu,
            n_scatter_mu,
            n_scatter_phi
        ))

    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(compute_one_mass_point_verified_only_sd_thermal, task) for task in tasks]
        for i, fut in enumerate(as_completed(futures), 1):
            res = fut.result()
            print(f"  done {i}/{len(futures)} : m_chi = {res['m']:.4g} GeV")
            results.append(res)

    results.sort(key=lambda x: x["m"])

    masses_sorted = np.array([r["m"] for r in results], dtype=float)
    C_SD_T = np.array([r["C_SD_T"] for r in results], dtype=float)
    C_SD_0 = np.array([r["C_SD_0"] for r in results], dtype=float)

    ratio_SD_raw = np.divide(
        C_SD_T,
        C_SD_0,
        out=np.full_like(C_SD_T, np.nan, dtype=float),
        where=(C_SD_0 > 0.0)
    )

    sd_floor = max(1e-6 * float(np.nanmax(C_SD_0)), 1e-300)
    ratio_SD_plot = np.where(C_SD_0 > sd_floor, ratio_SD_raw, np.nan)

    save_verified_only_sd_thermal_results_to_csv(
        DM_masses=masses_sorted,
        C_SD_T=C_SD_T,
        C_SD_0=C_SD_0,
        output_csv=csv_path
    )

    finite_mask = np.isfinite(ratio_SD_plot)
    if np.any(finite_mask):
        idx_max = int(np.nanargmax(ratio_SD_plot))
        idx_min = int(np.nanargmin(ratio_SD_plot))
        max_ratio = float(ratio_SD_plot[idx_max])
        min_ratio = float(ratio_SD_plot[idx_min])
        max_ratio_mass = float(masses_sorted[idx_max])
        min_ratio_mass = float(masses_sorted[idx_min])
    else:
        max_ratio = float("nan")
        min_ratio = float("nan")
        max_ratio_mass = float("nan")
        min_ratio_mass = float("nan")

    with open(log_path, "w", encoding="utf-8") as f:
        f.write("Refined verified-only SD thermal comparison run info\n")
        f.write("=" * 70 + "\n")
        f.write("Definition:\n")
        f.write("  - SD only\n")
        f.write("  - constant cross section\n")
        f.write("  - verified_only\n")
        f.write("  - compare T(r)!=0 vs T=0\n")
        f.write("  - refined settings from convergence study\n")
        f.write("  - validated active labels = {H, Al, Si29, Na}\n\n")
        f.write(f"sigma_SD_p         = {sigma_SD_p:.6e} cm^2\n")
        f.write(f"rho_chi            = {rho_chi:.6e} GeV/cm^3\n")
        f.write(f"v0                 = {v0:.6f} km/s\n")
        f.write(f"u_max              = {u_max:.6f} km/s\n")
        f.write(f"n_u                = {n_u}\n")
        f.write(f"n_t_speed          = {n_t_speed}\n")
        f.write(f"n_t_mu             = {n_t_mu}\n")
        f.write(f"n_scatter_mu       = {n_scatter_mu}\n")
        f.write(f"n_scatter_phi      = {n_scatter_phi}\n")
        f.write(f"shell_step         = 1\n")
        f.write(f"max_workers        = {max_workers}\n")
        f.write(f"mass range         = {masses_sorted.min():.6e} -> {masses_sorted.max():.6e} GeV\n")
        f.write(f"n_masses           = {len(masses_sorted)}\n\n")
        f.write("Summary:\n")
        f.write(f"C_SD_T(min,max)    = ({np.nanmin(C_SD_T):.6e}, {np.nanmax(C_SD_T):.6e})\n")
        f.write(f"C_SD_0(min,max)    = ({np.nanmin(C_SD_0):.6e}, {np.nanmax(C_SD_0):.6e})\n")
        f.write(f"max ratio C_T/C_0  = {max_ratio:.6e} at m_chi = {max_ratio_mass:.6e} GeV\n")
        f.write(f"min ratio C_T/C_0  = {min_ratio:.6e} at m_chi = {min_ratio_mass:.6e} GeV\n")

    print(f"[INFO] Refined verified-only SD thermal run info saved to: {log_path}")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    axL = axes[0]
    axL.loglog(
        masses_sorted,
        C_SD_T,
        color="purple",
        linestyle="-",
        linewidth=2.5,
        label="SD, verified_only, $T(r)\\neq 0$"
    )
    axL.loglog(
        masses_sorted,
        C_SD_0,
        color="purple",
        alpha=0.35,
        linestyle=":",
        linewidth=2.2,
        label="SD, verified_only, $T=0$"
    )
    axL.set_xlabel(r"$m_\chi\ [{\rm GeV}]$", fontsize=12)
    axL.set_ylabel(r"$C\ [{\rm s}^{-1}]$", fontsize=12)
    axL.set_title(r"Verified-only SD, constant, refined", fontsize=13)
    axL.grid(True, which="both", alpha=0.3)
    axL.legend(fontsize=10, loc="best")

    axR = axes[1]
    mask_plot = np.isfinite(ratio_SD_plot)
    axR.semilogx(
        masses_sorted[mask_plot],
        ratio_SD_plot[mask_plot],
        color="purple",
        linestyle="-",
        linewidth=2.3,
        label="SD"
    )
    axR.axhline(
        1.0,
        color="black",
        linestyle=":",
        linewidth=1.2,
        label=r"$T=0$ limit"
    )
    axR.set_xlabel(r"$m_\chi\ [{\rm GeV}]$", fontsize=12)
    axR.set_ylabel(r"$C(T)/C(T=0)$", fontsize=12)
    axR.set_title(r"Verified-only SD thermal correction, refined", fontsize=13)
    axR.set_yscale("log")
    axR.grid(True, which="both", alpha=0.3)
    axR.legend(fontsize=10, loc="best")

    plt.tight_layout()
    fig.savefig(fig_png, dpi=300, bbox_inches="tight")
    fig.savefig(fig_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"[INFO] Refined verified-only SD thermal figure saved to: {fig_png}")
    print(f"[INFO] Refined verified-only SD thermal figure saved to: {fig_pdf}")

    return {
        "DM_masses": masses_sorted,
        "C_SD_T": C_SD_T,
        "C_SD_0": C_SD_0,
        "ratio_SD": ratio_SD_raw,
        "csv_path": csv_path,
        "fig_png": fig_png,
        "fig_pdf": fig_pdf,
        "log_path": log_path,
    }

def compute_one_mass_point_verified_only_sd_operator(task):
    """
    Worker for one mass point:
        verified-only
        SD only
        T = 0
        compare:
            - constant
            - v2_dependent
            - q2_dependent
    """
    (
        earth_data,
        m,
        sigma_SD_p,
        rho_chi,
        u_max,
        n_u,
        v0,
        n_t_speed,
        n_t_mu,
        n_scatter_mu,
        n_scatter_phi,
        shell_step,
        u_grid_mode
    ) = task

    print(f"[worker:sd-operator] start m = {m:.4g} GeV", flush=True)

    c_constant = capture_rate_total(
        earth_data=earth_data,
        DM_mass=m,
        sigma_SD_p=sigma_SD_p,
        scattering_type="SD",
        cross_section_type="constant",
        rho_chi=rho_chi,
        u_max=u_max,
        n_u=n_u,
        v0=v0,
        include_thermal_targets=False,
        n_t_speed=n_t_speed,
        n_t_mu=n_t_mu,
        n_scatter_mu=n_scatter_mu,
        n_scatter_phi=n_scatter_phi,
        shell_step=shell_step,
        u_grid_mode=u_grid_mode
    )

    c_v2 = capture_rate_total(
        earth_data=earth_data,
        DM_mass=m,
        sigma_SD_p=sigma_SD_p,
        scattering_type="SD",
        cross_section_type="v2_dependent",
        rho_chi=rho_chi,
        u_max=u_max,
        n_u=n_u,
        v0=v0,
        include_thermal_targets=False,
        n_t_speed=n_t_speed,
        n_t_mu=n_t_mu,
        n_scatter_mu=n_scatter_mu,
        n_scatter_phi=n_scatter_phi,
        shell_step=shell_step,
        u_grid_mode=u_grid_mode
    )

    c_q2 = capture_rate_total(
        earth_data=earth_data,
        DM_mass=m,
        sigma_SD_p=sigma_SD_p,
        scattering_type="SD",
        cross_section_type="q2_dependent",
        rho_chi=rho_chi,
        u_max=u_max,
        n_u=n_u,
        v0=v0,
        include_thermal_targets=False,
        n_t_speed=n_t_speed,
        n_t_mu=n_t_mu,
        n_scatter_mu=n_scatter_mu,
        n_scatter_phi=n_scatter_phi,
        shell_step=shell_step,
        u_grid_mode=u_grid_mode
    )

    print(f"[worker:sd-operator] done  m = {m:.4g} GeV", flush=True)

    return {
        "m": float(m),
        "C_constant": float(c_constant),
        "C_v2": float(c_v2),
        "C_q2": float(c_q2),
    }


def save_verified_only_sd_operator_results_to_csv(
    DM_masses,
    C_constant,
    C_v2,
    C_q2,
    output_csv="results/verified_only_sd_operator_comparison_T0_sigma1e-40.csv"
):
    """
    Save verified-only SD operator-comparison results to CSV.

    Columns
    -------
    m_chi_GeV
    C_SD_constant_T0
    C_SD_v2_T0
    C_SD_q2_T0
    ratio_v2_over_constant
    ratio_q2_over_constant
    """
    DM_masses = np.asarray(DM_masses, dtype=float)
    C_constant = np.asarray(C_constant, dtype=float)
    C_v2 = np.asarray(C_v2, dtype=float)
    C_q2 = np.asarray(C_q2, dtype=float)

    if not (len(DM_masses) == len(C_constant) == len(C_v2) == len(C_q2)):
        raise ValueError("DM_masses, C_constant, C_v2, and C_q2 must have the same length.")

    ratio_v2 = np.divide(
        C_v2,
        C_constant,
        out=np.full_like(C_v2, np.nan, dtype=float),
        where=(C_constant > 0.0)
    )

    ratio_q2 = np.divide(
        C_q2,
        C_constant,
        out=np.full_like(C_q2, np.nan, dtype=float),
        where=(C_constant > 0.0)
    )

    out_dir = os.path.dirname(output_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "m_chi_GeV",
            "C_SD_constant_T0",
            "C_SD_v2_T0",
            "C_SD_q2_T0",
            "ratio_v2_over_constant",
            "ratio_q2_over_constant"
        ])
        for m, c0, cv2, cq2, rv2, rq2 in zip(DM_masses, C_constant, C_v2, C_q2, ratio_v2, ratio_q2):
            writer.writerow([
                f"{float(m):.8e}",
                f"{float(c0):.8e}",
                f"{float(cv2):.8e}",
                f"{float(cq2):.8e}",
                f"{float(rv2):.8e}" if np.isfinite(rv2) else "nan",
                f"{float(rq2):.8e}" if np.isfinite(rq2) else "nan",
            ])

    print(f"[INFO] Verified-only SD operator CSV saved to: {output_csv}")


def plot_verified_only_sd_operator_comparison(
    earth_data,
    DM_masses,
    sigma_SD_p=1e-40,
    output_root=".",
    u_max=800.0,
    n_u=80,
    rho_chi=RHO_CHI_DEFAULT,
    v0=V0_DEFAULT,
    n_t_speed=4,
    n_t_mu=4,
    n_scatter_mu=8,
    n_scatter_phi=12,
    shell_step=1,
    u_grid_mode="log",
    max_workers=None
):
    """
    Verified-only SD operator comparison:
        - SD only
        - verified_only
        - T = 0
        - compare:
            1. constant
            2. v2_dependent
            3. q2_dependent

    Notes
    -----
    For T = 0 runs, n_t_speed and n_t_mu are effectively irrelevant,
    but kept for interface consistency.
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed

    expected_labels = ("Al", "H", "Na", "Si29")
    validate_verified_only_sd_baseline(
        earth_data=earth_data,
        expected_labels=expected_labels
    )

    if max_workers is None:
        max_workers = max(1, os.cpu_count() - 1)

    DM_masses = np.asarray(DM_masses, dtype=float)
    sigma_tag = f"{sigma_SD_p:.0e}".replace("+", "")

    figures_dir = os.path.join(output_root, "figures")
    results_dir = os.path.join(output_root, "results")
    logs_dir = os.path.join(output_root, "logs")
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    fig_png = os.path.join(
        figures_dir,
        f"verified_only_sd_operator_comparison_T0_sigma{sigma_tag}.png"
    )
    fig_pdf = os.path.join(
        figures_dir,
        f"verified_only_sd_operator_comparison_T0_sigma{sigma_tag}.pdf"
    )
    csv_path = os.path.join(
        results_dir,
        f"verified_only_sd_operator_comparison_T0_sigma{sigma_tag}.csv"
    )
    log_path = os.path.join(
        logs_dir,
        f"verified_only_sd_operator_comparison_T0_sigma{sigma_tag}.txt"
    )

    print("\n" + "=" * 90)
    print("Generating verified-only SD operator comparison")
    print("=" * 90)
    print(f"Using {max_workers} worker processes")
    print(
        f"Settings: "
        f"n_u={n_u}, n_t_speed={n_t_speed}, n_t_mu={n_t_mu}, "
        f"n_scatter_mu={n_scatter_mu}, n_scatter_phi={n_scatter_phi}, "
        f"shell_step={shell_step}, u_grid_mode={u_grid_mode}"
    )

    tasks = []
    for m in DM_masses:
        tasks.append((
            earth_data,
            float(m),
            sigma_SD_p,
            rho_chi,
            u_max,
            n_u,
            v0,
            n_t_speed,
            n_t_mu,
            n_scatter_mu,
            n_scatter_phi,
            shell_step,
            u_grid_mode
        ))

    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(compute_one_mass_point_verified_only_sd_operator, task) for task in tasks]
        for i, fut in enumerate(as_completed(futures), 1):
            res = fut.result()
            print(f"  done {i}/{len(futures)} : m_chi = {res['m']:.4g} GeV")
            results.append(res)

    results.sort(key=lambda x: x["m"])

    masses_sorted = np.array([r["m"] for r in results], dtype=float)
    C_constant = np.array([r["C_constant"] for r in results], dtype=float)
    C_v2 = np.array([r["C_v2"] for r in results], dtype=float)
    C_q2 = np.array([r["C_q2"] for r in results], dtype=float)

    ratio_v2_raw = np.divide(
        C_v2,
        C_constant,
        out=np.full_like(C_v2, np.nan, dtype=float),
        where=(C_constant > 0.0)
    )
    ratio_q2_raw = np.divide(
        C_q2,
        C_constant,
        out=np.full_like(C_q2, np.nan, dtype=float),
        where=(C_constant > 0.0)
    )

    c_floor = max(1e-8 * float(np.nanmax(C_constant)), 1e-300)
    ratio_v2_plot = np.where(C_constant > c_floor, ratio_v2_raw, np.nan)
    ratio_q2_plot = np.where(C_constant > c_floor, ratio_q2_raw, np.nan)

    save_verified_only_sd_operator_results_to_csv(
        DM_masses=masses_sorted,
        C_constant=C_constant,
        C_v2=C_v2,
        C_q2=C_q2,
        output_csv=csv_path
    )

    finite_v2 = np.isfinite(ratio_v2_plot)
    finite_q2 = np.isfinite(ratio_q2_plot)

    if np.any(finite_v2):
        idx_v2_max = int(np.nanargmax(ratio_v2_plot))
        idx_v2_min = int(np.nanargmin(ratio_v2_plot))
        v2_ratio_max = float(ratio_v2_plot[idx_v2_max])
        v2_ratio_min = float(ratio_v2_plot[idx_v2_min])
        v2_mass_max = float(masses_sorted[idx_v2_max])
        v2_mass_min = float(masses_sorted[idx_v2_min])
    else:
        v2_ratio_max = float("nan")
        v2_ratio_min = float("nan")
        v2_mass_max = float("nan")
        v2_mass_min = float("nan")

    if np.any(finite_q2):
        idx_q2_max = int(np.nanargmax(ratio_q2_plot))
        idx_q2_min = int(np.nanargmin(ratio_q2_plot))
        q2_ratio_max = float(ratio_q2_plot[idx_q2_max])
        q2_ratio_min = float(ratio_q2_plot[idx_q2_min])
        q2_mass_max = float(masses_sorted[idx_q2_max])
        q2_mass_min = float(masses_sorted[idx_q2_min])
    else:
        q2_ratio_max = float("nan")
        q2_ratio_min = float("nan")
        q2_mass_max = float("nan")
        q2_mass_min = float("nan")

    with open(log_path, "w", encoding="utf-8") as f:
        f.write("Verified-only SD operator comparison run info\n")
        f.write("=" * 70 + "\n")
        f.write("Definition:\n")
        f.write("  - SD only\n")
        f.write("  - verified_only\n")
        f.write("  - T = 0\n")
        f.write("  - compare {constant, v2_dependent, q2_dependent}\n")
        f.write("  - validated active labels = {H, Al, Si29, Na}\n\n")
        f.write(f"sigma_SD_p         = {sigma_SD_p:.6e} cm^2\n")
        f.write(f"rho_chi            = {rho_chi:.6e} GeV/cm^3\n")
        f.write(f"v0                 = {v0:.6f} km/s\n")
        f.write(f"u_max              = {u_max:.6f} km/s\n")
        f.write(f"n_u                = {n_u}\n")
        f.write(f"n_t_speed          = {n_t_speed}\n")
        f.write(f"n_t_mu             = {n_t_mu}\n")
        f.write(f"n_scatter_mu       = {n_scatter_mu}\n")
        f.write(f"n_scatter_phi      = {n_scatter_phi}\n")
        f.write(f"shell_step         = {shell_step}\n")
        f.write(f"u_grid_mode        = {u_grid_mode}\n")
        f.write(f"max_workers        = {max_workers}\n")
        f.write(f"mass range         = {masses_sorted.min():.6e} -> {masses_sorted.max():.6e} GeV\n")
        f.write(f"n_masses           = {len(masses_sorted)}\n\n")
        f.write("Summary:\n")
        f.write(f"C_constant(min,max) = ({np.nanmin(C_constant):.6e}, {np.nanmax(C_constant):.6e})\n")
        f.write(f"C_v2(min,max)       = ({np.nanmin(C_v2):.6e}, {np.nanmax(C_v2):.6e})\n")
        f.write(f"C_q2(min,max)       = ({np.nanmin(C_q2):.6e}, {np.nanmax(C_q2):.6e})\n")
        f.write("\n")
        f.write(f"max(v2/constant)    = {v2_ratio_max:.6e} at m_chi = {v2_mass_max:.6e} GeV\n")
        f.write(f"min(v2/constant)    = {v2_ratio_min:.6e} at m_chi = {v2_mass_min:.6e} GeV\n")
        f.write(f"max(q2/constant)    = {q2_ratio_max:.6e} at m_chi = {q2_mass_max:.6e} GeV\n")
        f.write(f"min(q2/constant)    = {q2_ratio_min:.6e} at m_chi = {q2_mass_min:.6e} GeV\n")

    print(f"[INFO] Verified-only SD operator log saved to: {log_path}")

    # ------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    axL = axes[0]
    axL.loglog(
        masses_sorted,
        C_constant,
        color="black",
        linestyle="-",
        linewidth=2.6,
        label="verified-only SD constant"
    )
    axL.loglog(
        masses_sorted,
        C_v2,
        color="#1f77b4",
        linestyle="--",
        linewidth=2.4,
        label="verified-only SD v2"
    )
    axL.loglog(
        masses_sorted,
        C_q2,
        color="#d62728",
        linestyle="-.",
        linewidth=2.4,
        label="verified-only SD q2"
    )
    axL.set_xlabel(r"$m_\chi\ [{\rm GeV}]$", fontsize=12)
    axL.set_ylabel(r"$C\ [{\rm s}^{-1}]$", fontsize=12)
    axL.set_title(r"Verified-only SD operator comparison, $T=0$", fontsize=13)
    axL.grid(True, which="both", alpha=0.3)
    axL.legend(fontsize=10, loc="best")

    axR = axes[1]
    mask_v2 = np.isfinite(ratio_v2_plot)
    mask_q2 = np.isfinite(ratio_q2_plot)

    axR.semilogx(
        masses_sorted[mask_v2],
        ratio_v2_plot[mask_v2],
        color="#1f77b4",
        linestyle="--",
        linewidth=2.3,
        label=r"verified-only SD $v^2$ / constant"
    )
    axR.semilogx(
        masses_sorted[mask_q2],
        ratio_q2_plot[mask_q2],
        color="#d62728",
        linestyle="-.",
        linewidth=2.3,
        label=r"verified-only SD $q^2$ / constant"
    )
    axR.axhline(
        1.0,
        color="black",
        linestyle=":",
        linewidth=1.2,
        label="constant reference"
    )
    axR.set_xlabel(r"$m_\chi\ [{\rm GeV}]$", fontsize=12)
    axR.set_ylabel("ratio to constant", fontsize=12)
    axR.set_title(r"Operator dependence relative to constant", fontsize=13)
    axR.set_yscale("log")
    axR.grid(True, which="both", alpha=0.3)
    axR.legend(fontsize=10, loc="best")

    plt.tight_layout()
    fig.savefig(fig_png, dpi=300, bbox_inches="tight")
    fig.savefig(fig_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"[INFO] Verified-only SD operator figure saved to: {fig_png}")
    print(f"[INFO] Verified-only SD operator figure saved to: {fig_pdf}")

    return {
        "DM_masses": masses_sorted,
        "C_constant": C_constant,
        "C_v2": C_v2,
        "C_q2": C_q2,
        "ratio_v2_over_constant": ratio_v2_raw,
        "ratio_q2_over_constant": ratio_q2_raw,
        "csv_path": csv_path,
        "fig_png": fig_png,
        "fig_pdf": fig_pdf,
        "log_path": log_path,
    }

def compute_one_mass_point_verified_only_sd_operator_target_decomposition(task):
    """
    Worker for one mass point:
        verified-only
        SD only
        T = 0
        operators:
            - constant
            - v2_dependent
            - q2_dependent
        targets:
            - H
            - Al
            - Si29
            - Na
    """
    (
        earth_data,
        m,
        sigma_SD_p,
        rho_chi,
        u_max,
        n_u,
        v0,
        n_t_speed,
        n_t_mu,
        n_scatter_mu,
        n_scatter_phi,
        shell_step,
        u_grid_mode,
        label_order
    ) = task

    print(f"[worker:sd-op-target] start m = {m:.4g} GeV", flush=True)

    out = {"m": float(m)}

    earth_data_by_label = {
        label: build_sd_single_element_earth_data(earth_data, label)
        for label in label_order
    }

    operator_map = {
        "constant": "constant",
        "v2": "v2_dependent",
        "q2": "q2_dependent",
    }

    for op_tag, cross_type in operator_map.items():
        c_total = capture_rate_total(
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
            n_scatter_phi=n_scatter_phi,
            shell_step=shell_step,
            u_grid_mode=u_grid_mode
        )
        out[f"C_{op_tag}_total"] = float(c_total)

        for label in label_order:
            c_label = capture_rate_total(
                earth_data=earth_data_by_label[label],
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
                n_scatter_phi=n_scatter_phi,
                shell_step=shell_step,
                u_grid_mode=u_grid_mode
            )
            out[f"C_{op_tag}_{label}"] = float(c_label)

    print(f"[worker:sd-op-target] done  m = {m:.4g} GeV", flush=True)
    return out


def save_verified_only_sd_operator_target_decomposition_to_csv(
    DM_masses,
    C_constant_total,
    C_constant_H,
    C_constant_Al,
    C_constant_Si29,
    C_constant_Na,
    C_v2_total,
    C_v2_H,
    C_v2_Al,
    C_v2_Si29,
    C_v2_Na,
    C_q2_total,
    C_q2_H,
    C_q2_Al,
    C_q2_Si29,
    C_q2_Na,
    output_csv="results/verified_only_sd_operator_target_decomposition_T0_sigma1e-40.csv"
):
    """
    Save verified-only SD operator target decomposition to CSV.
    """
    DM_masses = np.asarray(DM_masses, dtype=float)

    C_constant_total = np.asarray(C_constant_total, dtype=float)
    C_constant_H = np.asarray(C_constant_H, dtype=float)
    C_constant_Al = np.asarray(C_constant_Al, dtype=float)
    C_constant_Si29 = np.asarray(C_constant_Si29, dtype=float)
    C_constant_Na = np.asarray(C_constant_Na, dtype=float)

    C_v2_total = np.asarray(C_v2_total, dtype=float)
    C_v2_H = np.asarray(C_v2_H, dtype=float)
    C_v2_Al = np.asarray(C_v2_Al, dtype=float)
    C_v2_Si29 = np.asarray(C_v2_Si29, dtype=float)
    C_v2_Na = np.asarray(C_v2_Na, dtype=float)

    C_q2_total = np.asarray(C_q2_total, dtype=float)
    C_q2_H = np.asarray(C_q2_H, dtype=float)
    C_q2_Al = np.asarray(C_q2_Al, dtype=float)
    C_q2_Si29 = np.asarray(C_q2_Si29, dtype=float)
    C_q2_Na = np.asarray(C_q2_Na, dtype=float)

    def frac(x, y):
        return np.divide(
            x,
            y,
            out=np.zeros_like(x, dtype=float),
            where=(y > 0.0)
        )

    f_constant_H = frac(C_constant_H, C_constant_total)
    f_constant_Al = frac(C_constant_Al, C_constant_total)
    f_constant_Si29 = frac(C_constant_Si29, C_constant_total)
    f_constant_Na = frac(C_constant_Na, C_constant_total)

    f_v2_H = frac(C_v2_H, C_v2_total)
    f_v2_Al = frac(C_v2_Al, C_v2_total)
    f_v2_Si29 = frac(C_v2_Si29, C_v2_total)
    f_v2_Na = frac(C_v2_Na, C_v2_total)

    f_q2_H = frac(C_q2_H, C_q2_total)
    f_q2_Al = frac(C_q2_Al, C_q2_total)
    f_q2_Si29 = frac(C_q2_Si29, C_q2_total)
    f_q2_Na = frac(C_q2_Na, C_q2_total)

    out_dir = os.path.dirname(output_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "m_chi_GeV",

            "C_constant_total",
            "C_constant_H",
            "C_constant_Al",
            "C_constant_Si29",
            "C_constant_Na",
            "f_constant_H",
            "f_constant_Al",
            "f_constant_Si29",
            "f_constant_Na",

            "C_v2_total",
            "C_v2_H",
            "C_v2_Al",
            "C_v2_Si29",
            "C_v2_Na",
            "f_v2_H",
            "f_v2_Al",
            "f_v2_Si29",
            "f_v2_Na",

            "C_q2_total",
            "C_q2_H",
            "C_q2_Al",
            "C_q2_Si29",
            "C_q2_Na",
            "f_q2_H",
            "f_q2_Al",
            "f_q2_Si29",
            "f_q2_Na",
        ])

        for i in range(len(DM_masses)):
            writer.writerow([
                f"{DM_masses[i]:.8e}",

                f"{C_constant_total[i]:.8e}",
                f"{C_constant_H[i]:.8e}",
                f"{C_constant_Al[i]:.8e}",
                f"{C_constant_Si29[i]:.8e}",
                f"{C_constant_Na[i]:.8e}",
                f"{f_constant_H[i]:.8e}",
                f"{f_constant_Al[i]:.8e}",
                f"{f_constant_Si29[i]:.8e}",
                f"{f_constant_Na[i]:.8e}",

                f"{C_v2_total[i]:.8e}",
                f"{C_v2_H[i]:.8e}",
                f"{C_v2_Al[i]:.8e}",
                f"{C_v2_Si29[i]:.8e}",
                f"{C_v2_Na[i]:.8e}",
                f"{f_v2_H[i]:.8e}",
                f"{f_v2_Al[i]:.8e}",
                f"{f_v2_Si29[i]:.8e}",
                f"{f_v2_Na[i]:.8e}",

                f"{C_q2_total[i]:.8e}",
                f"{C_q2_H[i]:.8e}",
                f"{C_q2_Al[i]:.8e}",
                f"{C_q2_Si29[i]:.8e}",
                f"{C_q2_Na[i]:.8e}",
                f"{f_q2_H[i]:.8e}",
                f"{f_q2_Al[i]:.8e}",
                f"{f_q2_Si29[i]:.8e}",
                f"{f_q2_Na[i]:.8e}",
            ])

    print(f"[INFO] Verified-only SD operator target decomposition CSV saved to: {output_csv}")


def build_verified_only_sd_operator_target_probe_summary_lines(
    DM_masses,
    C_constant_total,
    C_constant_H,
    C_constant_Al,
    C_constant_Si29,
    C_constant_Na,
    C_v2_total,
    C_v2_H,
    C_v2_Al,
    C_v2_Si29,
    C_v2_Na,
    C_q2_total,
    C_q2_H,
    C_q2_Al,
    C_q2_Si29,
    C_q2_Na,
    probe_masses=(0.7, 2.0, 25.0, 100.0)
):
    """
    Build human-readable probe-mass summary lines.
    """
    DM_masses = np.asarray(DM_masses, dtype=float)

    op_data = {
        "constant": {
            "total": np.asarray(C_constant_total, dtype=float),
            "H": np.asarray(C_constant_H, dtype=float),
            "Al": np.asarray(C_constant_Al, dtype=float),
            "Si29": np.asarray(C_constant_Si29, dtype=float),
            "Na": np.asarray(C_constant_Na, dtype=float),
        },
        "v2": {
            "total": np.asarray(C_v2_total, dtype=float),
            "H": np.asarray(C_v2_H, dtype=float),
            "Al": np.asarray(C_v2_Al, dtype=float),
            "Si29": np.asarray(C_v2_Si29, dtype=float),
            "Na": np.asarray(C_v2_Na, dtype=float),
        },
        "q2": {
            "total": np.asarray(C_q2_total, dtype=float),
            "H": np.asarray(C_q2_H, dtype=float),
            "Al": np.asarray(C_q2_Al, dtype=float),
            "Si29": np.asarray(C_q2_Si29, dtype=float),
            "Na": np.asarray(C_q2_Na, dtype=float),
        }
    }

    lines = []
    lines.append("=" * 90)
    lines.append("Verified-only SD operator target-decomposition probe summary")
    lines.append("=" * 90)

    for op_name, block in op_data.items():
        lines.append("")
        lines.append(f"Operator = {op_name}")
        lines.append("-" * 90)
        for mp in probe_masses:
            idx = int(np.argmin(np.abs(DM_masses - mp)))
            total = float(block["total"][idx])
            pieces = [
                ("H", float(block["H"][idx])),
                ("Al", float(block["Al"][idx])),
                ("Si29", float(block["Si29"][idx])),
                ("Na", float(block["Na"][idx])),
            ]
            pieces = [
                (lab, val, val / total if total > 0.0 else np.nan)
                for lab, val in pieces
            ]
            pieces.sort(key=lambda x: x[1], reverse=True)

            lines.append(f"m_chi = {DM_masses[idx]:.4g} GeV")
            lines.append(f"  C_total = {total:.6e} s^-1")
            for lab, val, frac in pieces:
                lines.append(f"  {lab:>4s} : {100.0 * frac:8.4f}%   (C = {val:.6e})")
        lines.append("-" * 90)

    return lines


def plot_verified_only_sd_operator_target_decomposition(
    earth_data,
    DM_masses,
    sigma_SD_p=1e-40,
    output_root=".",
    u_max=800.0,
    n_u=80,
    rho_chi=RHO_CHI_DEFAULT,
    v0=V0_DEFAULT,
    n_t_speed=4,
    n_t_mu=4,
    n_scatter_mu=8,
    n_scatter_phi=12,
    shell_step=1,
    u_grid_mode="log",
    max_workers=None
):
    """
    Verified-only SD per-target operator decomposition at T = 0.

    Produces:
        - CSV
        - PNG
        - PDF
        - log

    Figure layout:
        3 rows x 2 columns
            left  : absolute capture rates
            right : target fractions
        rows correspond to:
            constant, v2_dependent, q2_dependent
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed

    expected_labels = ("Al", "H", "Na", "Si29")
    label_order = ("H", "Al", "Si29", "Na")

    validate_verified_only_sd_baseline(
        earth_data=earth_data,
        expected_labels=expected_labels
    )

    if max_workers is None:
        max_workers = max(1, os.cpu_count() - 1)

    DM_masses = np.asarray(DM_masses, dtype=float)
    sigma_tag = f"{sigma_SD_p:.0e}".replace("+", "")

    figures_dir = os.path.join(output_root, "figures")
    results_dir = os.path.join(output_root, "results")
    logs_dir = os.path.join(output_root, "logs")
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    fig_png = os.path.join(
        figures_dir,
        f"verified_only_sd_operator_target_decomposition_T0_sigma{sigma_tag}.png"
    )
    fig_pdf = os.path.join(
        figures_dir,
        f"verified_only_sd_operator_target_decomposition_T0_sigma{sigma_tag}.pdf"
    )
    csv_path = os.path.join(
        results_dir,
        f"verified_only_sd_operator_target_decomposition_T0_sigma{sigma_tag}.csv"
    )
    log_path = os.path.join(
        logs_dir,
        f"verified_only_sd_operator_target_decomposition_T0_sigma{sigma_tag}.txt"
    )

    print("\n" + "=" * 100)
    print("Generating verified-only SD operator target decomposition at T = 0")
    print("=" * 100)
    print(f"Using {max_workers} worker processes")
    print(
        f"Settings: "
        f"n_u={n_u}, n_t_speed={n_t_speed}, n_t_mu={n_t_mu}, "
        f"n_scatter_mu={n_scatter_mu}, n_scatter_phi={n_scatter_phi}, "
        f"shell_step={shell_step}, u_grid_mode={u_grid_mode}"
    )

    tasks = []
    for m in DM_masses:
        tasks.append((
            earth_data,
            float(m),
            sigma_SD_p,
            rho_chi,
            u_max,
            n_u,
            v0,
            n_t_speed,
            n_t_mu,
            n_scatter_mu,
            n_scatter_phi,
            shell_step,
            u_grid_mode,
            label_order
        ))

    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(
                compute_one_mass_point_verified_only_sd_operator_target_decomposition,
                task
            )
            for task in tasks
        ]
        for i, fut in enumerate(as_completed(futures), 1):
            res = fut.result()
            print(f"  done {i}/{len(futures)} : m_chi = {res['m']:.4g} GeV")
            results.append(res)

    results.sort(key=lambda x: x["m"])

    masses_sorted = np.array([r["m"] for r in results], dtype=float)

    C_constant_total = np.array([r["C_constant_total"] for r in results], dtype=float)
    C_constant_H = np.array([r["C_constant_H"] for r in results], dtype=float)
    C_constant_Al = np.array([r["C_constant_Al"] for r in results], dtype=float)
    C_constant_Si29 = np.array([r["C_constant_Si29"] for r in results], dtype=float)
    C_constant_Na = np.array([r["C_constant_Na"] for r in results], dtype=float)

    C_v2_total = np.array([r["C_v2_total"] for r in results], dtype=float)
    C_v2_H = np.array([r["C_v2_H"] for r in results], dtype=float)
    C_v2_Al = np.array([r["C_v2_Al"] for r in results], dtype=float)
    C_v2_Si29 = np.array([r["C_v2_Si29"] for r in results], dtype=float)
    C_v2_Na = np.array([r["C_v2_Na"] for r in results], dtype=float)

    C_q2_total = np.array([r["C_q2_total"] for r in results], dtype=float)
    C_q2_H = np.array([r["C_q2_H"] for r in results], dtype=float)
    C_q2_Al = np.array([r["C_q2_Al"] for r in results], dtype=float)
    C_q2_Si29 = np.array([r["C_q2_Si29"] for r in results], dtype=float)
    C_q2_Na = np.array([r["C_q2_Na"] for r in results], dtype=float)

    def frac(x, y):
        return np.divide(
            x,
            y,
            out=np.zeros_like(x, dtype=float),
            where=(y > 0.0)
        )

    f_constant_H = frac(C_constant_H, C_constant_total)
    f_constant_Al = frac(C_constant_Al, C_constant_total)
    f_constant_Si29 = frac(C_constant_Si29, C_constant_total)
    f_constant_Na = frac(C_constant_Na, C_constant_total)

    f_v2_H = frac(C_v2_H, C_v2_total)
    f_v2_Al = frac(C_v2_Al, C_v2_total)
    f_v2_Si29 = frac(C_v2_Si29, C_v2_total)
    f_v2_Na = frac(C_v2_Na, C_v2_total)

    f_q2_H = frac(C_q2_H, C_q2_total)
    f_q2_Al = frac(C_q2_Al, C_q2_total)
    f_q2_Si29 = frac(C_q2_Si29, C_q2_total)
    f_q2_Na = frac(C_q2_Na, C_q2_total)

    save_verified_only_sd_operator_target_decomposition_to_csv(
        DM_masses=masses_sorted,

        C_constant_total=C_constant_total,
        C_constant_H=C_constant_H,
        C_constant_Al=C_constant_Al,
        C_constant_Si29=C_constant_Si29,
        C_constant_Na=C_constant_Na,

        C_v2_total=C_v2_total,
        C_v2_H=C_v2_H,
        C_v2_Al=C_v2_Al,
        C_v2_Si29=C_v2_Si29,
        C_v2_Na=C_v2_Na,

        C_q2_total=C_q2_total,
        C_q2_H=C_q2_H,
        C_q2_Al=C_q2_Al,
        C_q2_Si29=C_q2_Si29,
        C_q2_Na=C_q2_Na,

        output_csv=csv_path
    )

    summary_lines = build_verified_only_sd_operator_target_probe_summary_lines(
        DM_masses=masses_sorted,

        C_constant_total=C_constant_total,
        C_constant_H=C_constant_H,
        C_constant_Al=C_constant_Al,
        C_constant_Si29=C_constant_Si29,
        C_constant_Na=C_constant_Na,

        C_v2_total=C_v2_total,
        C_v2_H=C_v2_H,
        C_v2_Al=C_v2_Al,
        C_v2_Si29=C_v2_Si29,
        C_v2_Na=C_v2_Na,

        C_q2_total=C_q2_total,
        C_q2_H=C_q2_H,
        C_q2_Al=C_q2_Al,
        C_q2_Si29=C_q2_Si29,
        C_q2_Na=C_q2_Na,

        probe_masses=(0.7, 2.0, 25.0, 100.0)
    )

    for line in summary_lines:
        print(line)

    with open(log_path, "w", encoding="utf-8") as f:
        f.write("Verified-only SD operator target decomposition at T = 0\n")
        f.write("=" * 80 + "\n")
        f.write("Definition:\n")
        f.write("  - SD only\n")
        f.write("  - verified_only\n")
        f.write("  - T = 0\n")
        f.write("  - operators = {constant, v2_dependent, q2_dependent}\n")
        f.write("  - target decomposition = {H, Al, Si29, Na}\n")
        f.write("  - validated active labels = {H, Al, Si29, Na}\n\n")
        f.write(f"sigma_SD_p         = {sigma_SD_p:.6e} cm^2\n")
        f.write(f"rho_chi            = {rho_chi:.6e} GeV/cm^3\n")
        f.write(f"v0                 = {v0:.6f} km/s\n")
        f.write(f"u_max              = {u_max:.6f} km/s\n")
        f.write(f"n_u                = {n_u}\n")
        f.write(f"n_t_speed          = {n_t_speed}\n")
        f.write(f"n_t_mu             = {n_t_mu}\n")
        f.write(f"n_scatter_mu       = {n_scatter_mu}\n")
        f.write(f"n_scatter_phi      = {n_scatter_phi}\n")
        f.write(f"shell_step         = {shell_step}\n")
        f.write(f"u_grid_mode        = {u_grid_mode}\n")
        f.write(f"max_workers        = {max_workers}\n")
        f.write(f"mass range         = {masses_sorted.min():.6e} -> {masses_sorted.max():.6e} GeV\n")
        f.write(f"n_masses           = {len(masses_sorted)}\n\n")

        f.write("Absolute ranges:\n")
        f.write(f"  C_constant_total(min,max) = ({np.nanmin(C_constant_total):.6e}, {np.nanmax(C_constant_total):.6e})\n")
        f.write(f"  C_v2_total(min,max)       = ({np.nanmin(C_v2_total):.6e}, {np.nanmax(C_v2_total):.6e})\n")
        f.write(f"  C_q2_total(min,max)       = ({np.nanmin(C_q2_total):.6e}, {np.nanmax(C_q2_total):.6e})\n\n")

        for line in summary_lines:
            f.write(line + "\n")

    print(f"[INFO] Verified-only SD operator target decomposition log saved to: {log_path}")

    # ------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------
    fig, axes = plt.subplots(3, 2, figsize=(14, 16), sharex=True)

    color_total = "black"
    color_H = "#2ca02c"
    color_Al = "#d62728"
    color_Si29 = "#1f77b4"
    color_Na = "#ff7f0e"

    # ---------------- Row 1: constant ----------------
    axL = axes[0, 0]
    axL.loglog(masses_sorted, C_constant_total, color=color_total, lw=2.6, label="total")
    axL.loglog(masses_sorted, C_constant_H, color=color_H, lw=2.0, label="H")
    axL.loglog(masses_sorted, C_constant_Al, color=color_Al, lw=2.0, label="Al")
    axL.loglog(masses_sorted, C_constant_Si29, color=color_Si29, lw=2.0, label="Si29")
    axL.loglog(masses_sorted, C_constant_Na, color=color_Na, lw=2.0, label="Na")
    axL.set_ylabel(r"$C\ [{\rm s}^{-1}]$", fontsize=12)
    axL.set_title(r"Verified-only SD constant, $T=0$", fontsize=13)
    axL.grid(True, which="both", alpha=0.3)
    axL.legend(fontsize=9, frameon=False, ncol=2)

    axR = axes[0, 1]
    axR.semilogx(masses_sorted, f_constant_H, color=color_H, lw=2.0, label="H")
    axR.semilogx(masses_sorted, f_constant_Al, color=color_Al, lw=2.0, label="Al")
    axR.semilogx(masses_sorted, f_constant_Si29, color=color_Si29, lw=2.0, label="Si29")
    axR.semilogx(masses_sorted, f_constant_Na, color=color_Na, lw=2.0, label="Na")
    axR.set_ylabel("fraction", fontsize=12)
    axR.set_ylim(0.0, 1.02)
    axR.set_title(r"Target fractions: constant", fontsize=13)
    axR.grid(True, which="both", alpha=0.3)
    axR.legend(fontsize=9, frameon=False, loc="best")

    # ---------------- Row 2: v2 ----------------
    axL = axes[1, 0]
    axL.loglog(masses_sorted, C_v2_total, color=color_total, lw=2.6, label="total")
    axL.loglog(masses_sorted, C_v2_H, color=color_H, lw=2.0, label="H")
    axL.loglog(masses_sorted, C_v2_Al, color=color_Al, lw=2.0, label="Al")
    axL.loglog(masses_sorted, C_v2_Si29, color=color_Si29, lw=2.0, label="Si29")
    axL.loglog(masses_sorted, C_v2_Na, color=color_Na, lw=2.0, label="Na")
    axL.set_ylabel(r"$C\ [{\rm s}^{-1}]$", fontsize=12)
    axL.set_title(r"Verified-only SD $v^2$ dependent, $T=0$", fontsize=13)
    axL.grid(True, which="both", alpha=0.3)
    axL.legend(fontsize=9, frameon=False, ncol=2)

    axR = axes[1, 1]
    axR.semilogx(masses_sorted, f_v2_H, color=color_H, lw=2.0, label="H")
    axR.semilogx(masses_sorted, f_v2_Al, color=color_Al, lw=2.0, label="Al")
    axR.semilogx(masses_sorted, f_v2_Si29, color=color_Si29, lw=2.0, label="Si29")
    axR.semilogx(masses_sorted, f_v2_Na, color=color_Na, lw=2.0, label="Na")
    axR.set_ylabel("fraction", fontsize=12)
    axR.set_ylim(0.0, 1.02)
    axR.set_title(r"Target fractions: $v^2$ dependent", fontsize=13)
    axR.grid(True, which="both", alpha=0.3)
    axR.legend(fontsize=9, frameon=False, loc="best")

    # ---------------- Row 3: q2 ----------------
    axL = axes[2, 0]
    axL.loglog(masses_sorted, C_q2_total, color=color_total, lw=2.6, label="total")
    axL.loglog(masses_sorted, C_q2_H, color=color_H, lw=2.0, label="H")
    axL.loglog(masses_sorted, C_q2_Al, color=color_Al, lw=2.0, label="Al")
    axL.loglog(masses_sorted, C_q2_Si29, color=color_Si29, lw=2.0, label="Si29")
    axL.loglog(masses_sorted, C_q2_Na, color=color_Na, lw=2.0, label="Na")
    axL.set_xlabel(r"$m_\chi\ [{\rm GeV}]$", fontsize=12)
    axL.set_ylabel(r"$C\ [{\rm s}^{-1}]$", fontsize=12)
    axL.set_title(r"Verified-only SD $q^2$ dependent, $T=0$", fontsize=13)
    axL.grid(True, which="both", alpha=0.3)
    axL.legend(fontsize=9, frameon=False, ncol=2)

    axR = axes[2, 1]
    axR.semilogx(masses_sorted, f_q2_H, color=color_H, lw=2.0, label="H")
    axR.semilogx(masses_sorted, f_q2_Al, color=color_Al, lw=2.0, label="Al")
    axR.semilogx(masses_sorted, f_q2_Si29, color=color_Si29, lw=2.0, label="Si29")
    axR.semilogx(masses_sorted, f_q2_Na, color=color_Na, lw=2.0, label="Na")
    axR.set_xlabel(r"$m_\chi\ [{\rm GeV}]$", fontsize=12)
    axR.set_ylabel("fraction", fontsize=12)
    axR.set_ylim(0.0, 1.02)
    axR.set_title(r"Target fractions: $q^2$ dependent", fontsize=13)
    axR.grid(True, which="both", alpha=0.3)
    axR.legend(fontsize=9, frameon=False, loc="best")

    plt.tight_layout()
    fig.savefig(fig_png, dpi=300, bbox_inches="tight")
    fig.savefig(fig_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"[INFO] Verified-only SD operator target decomposition figure saved to: {fig_png}")
    print(f"[INFO] Verified-only SD operator target decomposition figure saved to: {fig_pdf}")

    return {
        "DM_masses": masses_sorted,

        "C_constant_total": C_constant_total,
        "C_constant_H": C_constant_H,
        "C_constant_Al": C_constant_Al,
        "C_constant_Si29": C_constant_Si29,
        "C_constant_Na": C_constant_Na,
        "f_constant_H": f_constant_H,
        "f_constant_Al": f_constant_Al,
        "f_constant_Si29": f_constant_Si29,
        "f_constant_Na": f_constant_Na,

        "C_v2_total": C_v2_total,
        "C_v2_H": C_v2_H,
        "C_v2_Al": C_v2_Al,
        "C_v2_Si29": C_v2_Si29,
        "C_v2_Na": C_v2_Na,
        "f_v2_H": f_v2_H,
        "f_v2_Al": f_v2_Al,
        "f_v2_Si29": f_v2_Si29,
        "f_v2_Na": f_v2_Na,

        "C_q2_total": C_q2_total,
        "C_q2_H": C_q2_H,
        "C_q2_Al": C_q2_Al,
        "C_q2_Si29": C_q2_Si29,
        "C_q2_Na": C_q2_Na,
        "f_q2_H": f_q2_H,
        "f_q2_Al": f_q2_Al,
        "f_q2_Si29": f_q2_Si29,
        "f_q2_Na": f_q2_Na,

        "csv_path": csv_path,
        "fig_png": fig_png,
        "fig_pdf": fig_pdf,
        "log_path": log_path,
    }

def compute_one_mass_point_combined_constant_T0(task):
    """
    Worker for one mass point:
        - SI, constant, T = 0
        - verified-only SD, constant, T = 0
        - electron, constant, T = 0
        - geometric reference
    """
    (
        earth_data,
        m,
        sigma_SI_p,
        sigma_SD_p,
        sigma_electron,
        rho_chi,
        u_max,
        n_u,
        v0,
        n_t_speed,
        n_t_mu,
        n_scatter_mu,
        n_scatter_phi,
        shell_step,
        u_grid_mode,
        R_earth_cm
    ) = task

    print(f"[worker:combined-constant-T0] start m = {m:.4g} GeV", flush=True)

    c_si_0 = capture_rate_total(
        earth_data=earth_data,
        DM_mass=m,
        sigma_SI_p=sigma_SI_p,
        scattering_type="SI",
        cross_section_type="constant",
        rho_chi=rho_chi,
        u_max=u_max,
        n_u=n_u,
        v0=v0,
        include_thermal_targets=False,
        n_t_speed=n_t_speed,
        n_t_mu=n_t_mu,
        n_scatter_mu=n_scatter_mu,
        n_scatter_phi=n_scatter_phi,
        shell_step=shell_step,
        u_grid_mode=u_grid_mode
    )

    c_sd_0 = capture_rate_total(
        earth_data=earth_data,
        DM_mass=m,
        sigma_SD_p=sigma_SD_p,
        scattering_type="SD",
        cross_section_type="constant",
        rho_chi=rho_chi,
        u_max=u_max,
        n_u=n_u,
        v0=v0,
        include_thermal_targets=False,
        n_t_speed=n_t_speed,
        n_t_mu=n_t_mu,
        n_scatter_mu=n_scatter_mu,
        n_scatter_phi=n_scatter_phi,
        shell_step=shell_step,
        u_grid_mode=u_grid_mode
    )

    c_e_0 = capture_rate_total(
        earth_data=earth_data,
        DM_mass=m,
        sigma_electron=sigma_electron,
        scattering_type="electron",
        cross_section_type="constant",
        rho_chi=rho_chi,
        u_max=u_max,
        n_u=n_u,
        v0=v0,
        include_thermal_targets=False,
        n_t_speed=n_t_speed,
        n_t_mu=n_t_mu,
        n_scatter_mu=n_scatter_mu,
        n_scatter_phi=n_scatter_phi,
        shell_step=shell_step,
        u_grid_mode=u_grid_mode
    )

    c_geo = capture_rate_geometric(
        DM_mass=m,
        R_earth_cm=R_earth_cm,
        rho_chi=rho_chi,
        v0=v0,
        vesc_surface=11.2
    )

    print(f"[worker:combined-constant-T0] done  m = {m:.4g} GeV", flush=True)

    return {
        "m": float(m),
        "C_SI_0": float(c_si_0),
        "C_SD_0": float(c_sd_0),
        "C_e_0": float(c_e_0),
        "C_geo": float(c_geo),
    }


def save_combined_constant_T0_results_to_csv(
    DM_masses,
    C_SI_0,
    C_SD_0,
    C_e_0,
    C_geo,
    output_csv="results/combined_SI_verifiedSD_electron_constant_T0_sigma1e-40.csv"
):
    """
    Save combined SI / verified-only SD / electron constant T=0 results to CSV.
    """
    DM_masses = np.asarray(DM_masses, dtype=float)
    C_SI_0 = np.asarray(C_SI_0, dtype=float)
    C_SD_0 = np.asarray(C_SD_0, dtype=float)
    C_e_0 = np.asarray(C_e_0, dtype=float)
    C_geo = np.asarray(C_geo, dtype=float)

    if not (len(DM_masses) == len(C_SI_0) == len(C_SD_0) == len(C_e_0) == len(C_geo)):
        raise ValueError("All input arrays must have the same length.")

    ratio_SI_geo = np.divide(
        C_SI_0,
        C_geo,
        out=np.full_like(C_SI_0, np.nan, dtype=float),
        where=(C_geo > 0.0)
    )
    ratio_SD_geo = np.divide(
        C_SD_0,
        C_geo,
        out=np.full_like(C_SD_0, np.nan, dtype=float),
        where=(C_geo > 0.0)
    )
    ratio_e_geo = np.divide(
        C_e_0,
        C_geo,
        out=np.full_like(C_e_0, np.nan, dtype=float),
        where=(C_geo > 0.0)
    )

    out_dir = os.path.dirname(output_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "m_chi_GeV",
            "C_SI_0",
            "C_SD_0_verified_only",
            "C_e_0",
            "C_geo",
            "ratio_SI_over_geo",
            "ratio_SD_over_geo",
            "ratio_e_over_geo",
        ])
        for m, csi, csd, ce, cg, rsi, rsd, re in zip(
            DM_masses, C_SI_0, C_SD_0, C_e_0, C_geo,
            ratio_SI_geo, ratio_SD_geo, ratio_e_geo
        ):
            writer.writerow([
                f"{m:.8e}",
                f"{csi:.8e}",
                f"{csd:.8e}",
                f"{ce:.8e}",
                f"{cg:.8e}",
                f"{rsi:.8e}" if np.isfinite(rsi) else "nan",
                f"{rsd:.8e}" if np.isfinite(rsd) else "nan",
                f"{re:.8e}" if np.isfinite(re) else "nan",
            ])

    print(f"[INFO] Combined constant T=0 CSV saved to: {output_csv}")


def plot_combined_SI_verified_only_SD_electron_operator_grid_polished(
    earth_data,
    DM_masses,
    sigma_SI_p=1e-40,
    sigma_SD_p=1e-40,
    sigma_electron=1e-40,
    output_root=".",
    u_max=800.0,
    n_u=80,
    rho_chi=RHO_CHI_DEFAULT,
    v0=V0_DEFAULT,
    n_t_speed=4,
    n_t_mu=4,
    n_scatter_mu=8,
    n_scatter_phi=12,
    shell_step=1,
    u_grid_mode="log",
    max_workers=None,
    show_geometric_top_left_only=True,
    add_figure_note=True
):
    """
    Final-style combined operator grid:
        rows   = constant / v2_dependent / q2_dependent
        left   = absolute capture rates
        right  = informative comparison panels

    New plotting policy
    -------------------
    1. top-right panel:
         constant baseline hierarchy relative to SI constant
    2. geometric line:
         by default shown only in the top-left panel
    3. legend text:
         simplified for paper-style readability
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed

    expected_labels = ("Al", "H", "Na", "Si29")
    validate_verified_only_sd_baseline(
        earth_data=earth_data,
        expected_labels=expected_labels
    )

    if max_workers is None:
        max_workers = max(1, os.cpu_count() - 1)

    DM_masses = np.asarray(DM_masses, dtype=float)
    R_earth_cm = earth_data["radius"][-1] * 1e5

    if np.isclose(sigma_SI_p, sigma_SD_p) and np.isclose(sigma_SI_p, sigma_electron):
        sigma_tag = f"{sigma_SI_p:.0e}".replace("+", "")
        file_tag = f"sigma{sigma_tag}"
    else:
        si_tag = f"{sigma_SI_p:.0e}".replace("+", "")
        sd_tag = f"{sigma_SD_p:.0e}".replace("+", "")
        e_tag = f"{sigma_electron:.0e}".replace("+", "")
        file_tag = f"si{si_tag}_sd{sd_tag}_e{e_tag}"

    figures_dir = os.path.join(output_root, "figures")
    results_dir = os.path.join(output_root, "results")
    logs_dir = os.path.join(output_root, "logs")
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    fig_png = os.path.join(
        figures_dir,
        f"combined_SI_verifiedSD_electron_operator_grid_{file_tag}.png"
    )
    fig_pdf = os.path.join(
        figures_dir,
        f"combined_SI_verifiedSD_electron_operator_grid_{file_tag}.pdf"
    )
    csv_path = os.path.join(
        results_dir,
        f"combined_SI_verifiedSD_electron_operator_grid_{file_tag}.csv"
    )
    log_path = os.path.join(
        logs_dir,
        f"combined_SI_verifiedSD_electron_operator_grid_{file_tag}.txt"
    )

    print("\n" + "=" * 100)
    print("Generating combined SI / verified-only SD / electron operator grid")
    print("=" * 100)
    print(f"Using {max_workers} worker processes")
    print(
        f"Settings: "
        f"n_u={n_u}, n_t_speed={n_t_speed}, n_t_mu={n_t_mu}, "
        f"n_scatter_mu={n_scatter_mu}, n_scatter_phi={n_scatter_phi}, "
        f"shell_step={shell_step}, u_grid_mode={u_grid_mode}"
    )

    tasks = []
    for m in DM_masses:
        tasks.append((
            earth_data,
            float(m),
            sigma_SI_p,
            sigma_SD_p,
            sigma_electron,
            rho_chi,
            u_max,
            n_u,
            v0,
            n_t_speed,
            n_t_mu,
            n_scatter_mu,
            n_scatter_phi,
            shell_step,
            u_grid_mode,
            R_earth_cm
        ))

    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(compute_one_mass_point_combined_operator_grid, task)
            for task in tasks
        ]
        for i, fut in enumerate(as_completed(futures), 1):
            res = fut.result()
            print(f"  done {i}/{len(futures)} : m_chi = {res['m']:.4g} GeV")
            results.append(res)

    results.sort(key=lambda x: x["m"])

    masses_sorted = np.array([r["m"] for r in results], dtype=float)
    C_geo = np.array([r["C_geo"] for r in results], dtype=float)

    C_SI_constant = np.array([r["C_SI_constant"] for r in results], dtype=float)
    C_SI_v2 = np.array([r["C_SI_v2"] for r in results], dtype=float)
    C_SI_q2 = np.array([r["C_SI_q2"] for r in results], dtype=float)

    C_SD_constant = np.array([r["C_SD_constant"] for r in results], dtype=float)
    C_SD_v2 = np.array([r["C_SD_v2"] for r in results], dtype=float)
    C_SD_q2 = np.array([r["C_SD_q2"] for r in results], dtype=float)

    C_e_constant = np.array([r["C_e_constant"] for r in results], dtype=float)
    C_e_v2 = np.array([r["C_e_v2"] for r in results], dtype=float)
    C_e_q2 = np.array([r["C_e_q2"] for r in results], dtype=float)

    # ------------------------------------------------------------
    # Ratios to same-channel constant baseline
    # ------------------------------------------------------------
    ratio_SI_v2 = np.divide(
        C_SI_v2,
        C_SI_constant,
        out=np.full_like(C_SI_v2, np.nan, dtype=float),
        where=(C_SI_constant > 0.0)
    )
    ratio_SI_q2 = np.divide(
        C_SI_q2,
        C_SI_constant,
        out=np.full_like(C_SI_q2, np.nan, dtype=float),
        where=(C_SI_constant > 0.0)
    )

    ratio_SD_v2 = np.divide(
        C_SD_v2,
        C_SD_constant,
        out=np.full_like(C_SD_v2, np.nan, dtype=float),
        where=(C_SD_constant > 0.0)
    )
    ratio_SD_q2 = np.divide(
        C_SD_q2,
        C_SD_constant,
        out=np.full_like(C_SD_q2, np.nan, dtype=float),
        where=(C_SD_constant > 0.0)
    )

    ratio_e_v2 = np.divide(
        C_e_v2,
        C_e_constant,
        out=np.full_like(C_e_v2, np.nan, dtype=float),
        where=(C_e_constant > 0.0)
    )
    ratio_e_q2 = np.divide(
        C_e_q2,
        C_e_constant,
        out=np.full_like(C_e_q2, np.nan, dtype=float),
        where=(C_e_constant > 0.0)
    )

    # ------------------------------------------------------------
    # New top-right informative baseline hierarchy ratios
    # relative to SI constant
    # ------------------------------------------------------------
    ratio_SDconst_over_SIconst = np.divide(
        C_SD_constant,
        C_SI_constant,
        out=np.full_like(C_SD_constant, np.nan, dtype=float),
        where=(C_SI_constant > 0.0)
    )
    ratio_econst_over_SIconst = np.divide(
        C_e_constant,
        C_SI_constant,
        out=np.full_like(C_e_constant, np.nan, dtype=float),
        where=(C_SI_constant > 0.0)
    )
    ratio_SIconst_over_SIconst = np.ones_like(C_SI_constant, dtype=float)

    data = {
        "DM_masses": masses_sorted,
        "C_geo": C_geo,
        "C_SI_constant": C_SI_constant,
        "C_SI_v2": C_SI_v2,
        "C_SI_q2": C_SI_q2,
        "C_SD_constant": C_SD_constant,
        "C_SD_v2": C_SD_v2,
        "C_SD_q2": C_SD_q2,
        "C_e_constant": C_e_constant,
        "C_e_v2": C_e_v2,
        "C_e_q2": C_e_q2,
    }

    save_combined_operator_grid_results_to_csv(
        data=data,
        output_csv=csv_path
    )

    with open(log_path, "w", encoding="utf-8") as f:
        f.write("Combined SI / verified-only SD / electron operator grid\n")
        f.write("=" * 80 + "\n")
        f.write("Definition:\n")
        f.write("  - rows = {constant, v2_dependent, q2_dependent}\n")
        f.write("  - left panels = absolute capture rates\n")
        f.write("  - top-right = constant baseline hierarchy relative to SI constant\n")
        f.write("  - middle/bottom-right = ratio to same-channel constant\n")
        f.write("  - electron thermal motion disabled in baseline\n")
        f.write("  - SD mode = verified_only\n")
        f.write(f"  - show_geometric_top_left_only = {show_geometric_top_left_only}\n\n")

        f.write(f"sigma_SI_p         = {sigma_SI_p:.6e} cm^2\n")
        f.write(f"sigma_SD_p         = {sigma_SD_p:.6e} cm^2\n")
        f.write(f"sigma_electron     = {sigma_electron:.6e} cm^2\n")
        f.write(f"rho_chi            = {rho_chi:.6e} GeV/cm^3\n")
        f.write(f"v0                 = {v0:.6f} km/s\n")
        f.write(f"u_max              = {u_max:.6f} km/s\n")
        f.write(f"n_u                = {n_u}\n")
        f.write(f"n_t_speed          = {n_t_speed}\n")
        f.write(f"n_t_mu             = {n_t_mu}\n")
        f.write(f"n_scatter_mu       = {n_scatter_mu}\n")
        f.write(f"n_scatter_phi      = {n_scatter_phi}\n")
        f.write(f"shell_step         = {shell_step}\n")
        f.write(f"u_grid_mode        = {u_grid_mode}\n")
        f.write(f"max_workers        = {max_workers}\n")
        f.write(f"mass range         = {masses_sorted.min():.6e} -> {masses_sorted.max():.6e} GeV\n")
        f.write(f"n_masses           = {len(masses_sorted)}\n\n")

        f.write("Absolute ranges:\n")
        f.write(f"  C_SI_constant(min,max) = ({np.nanmin(C_SI_constant):.6e}, {np.nanmax(C_SI_constant):.6e})\n")
        f.write(f"  C_SI_v2(min,max)       = ({np.nanmin(C_SI_v2):.6e}, {np.nanmax(C_SI_v2):.6e})\n")
        f.write(f"  C_SI_q2(min,max)       = ({np.nanmin(C_SI_q2):.6e}, {np.nanmax(C_SI_q2):.6e})\n")
        f.write(f"  C_SD_constant(min,max) = ({np.nanmin(C_SD_constant):.6e}, {np.nanmax(C_SD_constant):.6e})\n")
        f.write(f"  C_SD_v2(min,max)       = ({np.nanmin(C_SD_v2):.6e}, {np.nanmax(C_SD_v2):.6e})\n")
        f.write(f"  C_SD_q2(min,max)       = ({np.nanmin(C_SD_q2):.6e}, {np.nanmax(C_SD_q2):.6e})\n")
        f.write(f"  C_e_constant(min,max)  = ({np.nanmin(C_e_constant):.6e}, {np.nanmax(C_e_constant):.6e})\n")
        f.write(f"  C_e_v2(min,max)        = ({np.nanmin(C_e_v2):.6e}, {np.nanmax(C_e_v2):.6e})\n")
        f.write(f"  C_e_q2(min,max)        = ({np.nanmin(C_e_q2):.6e}, {np.nanmax(C_e_q2):.6e})\n\n")

        f.write("Ratio ranges:\n")
        f.write(f"  SI  v2/constant(min,max) = ({np.nanmin(ratio_SI_v2):.6e}, {np.nanmax(ratio_SI_v2):.6e})\n")
        f.write(f"  SI  q2/constant(min,max) = ({np.nanmin(ratio_SI_q2):.6e}, {np.nanmax(ratio_SI_q2):.6e})\n")
        f.write(f"  SD  v2/constant(min,max) = ({np.nanmin(ratio_SD_v2):.6e}, {np.nanmax(ratio_SD_v2):.6e})\n")
        f.write(f"  SD  q2/constant(min,max) = ({np.nanmin(ratio_SD_q2):.6e}, {np.nanmax(ratio_SD_q2):.6e})\n")
        f.write(f"  e   v2/constant(min,max) = ({np.nanmin(ratio_e_v2):.6e}, {np.nanmax(ratio_e_v2):.6e})\n")
        f.write(f"  e   q2/constant(min,max) = ({np.nanmin(ratio_e_q2):.6e}, {np.nanmax(ratio_e_q2):.6e})\n")
        f.write(f"  SD const / SI const(min,max) = ({np.nanmin(ratio_SDconst_over_SIconst):.6e}, {np.nanmax(ratio_SDconst_over_SIconst):.6e})\n")
        f.write(f"  e  const / SI const(min,max) = ({np.nanmin(ratio_econst_over_SIconst):.6e}, {np.nanmax(ratio_econst_over_SIconst):.6e})\n")

    print(f"[INFO] Combined operator-grid log saved to: {log_path}")

    fig, axes = plt.subplots(3, 2, figsize=(14, 16), sharex="col")
    fig.suptitle(
        r"Combined SI / verified-only SD / electron operator comparison",
        fontsize=15,
        y=0.995
    )

    color_e = "red"
    color_sd = "green"
    color_si = "blue"
    color_geo = "black"

    ls_e = "-"
    ls_sd = "--"
    ls_si = ":"
    ls_geo = (0, (3, 3))

    sigma_text = (
        rf"$\sigma_{{\mathrm{{SI}},p}}={sigma_SI_p:.0e}\,\mathrm{{cm}}^2$" + "\n" +
        rf"$\sigma_{{\mathrm{{SD}},p}}={sigma_SD_p:.0e}\,\mathrm{{cm}}^2$" + "\n" +
        rf"$\sigma_e={sigma_electron:.0e}\,\mathrm{{cm}}^2$"
    )

    def _mask_pos(y):
        y = np.asarray(y, dtype=float)
        return np.isfinite(y) & (y > 0.0)

    # ------------------------------------------------------------
    # Row 1: constant
    # ------------------------------------------------------------
    axL = axes[0, 0]
    axL.loglog(
        masses_sorted, C_e_constant,
        color=color_e, linestyle=ls_e, linewidth=2.2, label="Electrons"
    )
    axL.loglog(
        masses_sorted, C_SD_constant,
        color=color_sd, linestyle=ls_sd, linewidth=2.2, label="Nucleons-SD"
    )
    axL.loglog(
        masses_sorted, C_SI_constant,
        color=color_si, linestyle=ls_si, linewidth=2.6, label="Nucleons-SI"
    )
    axL.loglog(
        masses_sorted, C_geo,
        color=color_geo, linestyle=ls_geo, linewidth=1.5, alpha=0.75, label="Geometric"
    )
    axL.set_ylabel(r"$C\ [{\rm s}^{-1}]$", fontsize=12)
    axL.set_title(r"constant: absolute rates", fontsize=13)
    axL.grid(True, which="both", alpha=0.3)
    axL.legend(fontsize=9, loc="best")
    axL.text(
        0.98, 0.97,
        sigma_text,
        transform=axL.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85)
    )

    axR = axes[0, 1]
    mask_top_si = _mask_pos(ratio_SIconst_over_SIconst)
    mask_top_sd = _mask_pos(ratio_SDconst_over_SIconst)
    mask_top_e = _mask_pos(ratio_econst_over_SIconst)

    axR.semilogx(
        masses_sorted[mask_top_e], ratio_econst_over_SIconst[mask_top_e],
        color=color_e, linestyle=ls_e, linewidth=2.2, label="Electrons"
    )
    axR.semilogx(
        masses_sorted[mask_top_sd], ratio_SDconst_over_SIconst[mask_top_sd],
        color=color_sd, linestyle=ls_sd, linewidth=2.2, label="Nucleons-SD"
    )
    axR.semilogx(
        masses_sorted[mask_top_si], ratio_SIconst_over_SIconst[mask_top_si],
        color=color_si, linestyle=ls_si, linewidth=2.6, label="Nucleons-SI"
    )
    axR.axhline(1.0, color="black", linestyle=":", linewidth=1.1, alpha=0.8)
    axR.set_ylabel(r"ratio to $C_{\rm SI}^{\rm constant}$", fontsize=12)
    axR.set_title(r"constant baseline hierarchy", fontsize=13)
    axR.set_yscale("log")
    _set_ratio_panel_ylim(
        axR,
        [ratio_econst_over_SIconst, ratio_SDconst_over_SIconst, ratio_SIconst_over_SIconst],
        top_row=False
    )
    axR.grid(True, which="both", alpha=0.3)
    axR.legend(fontsize=9, loc="best")

    # ------------------------------------------------------------
    # Row 2: v2
    # ------------------------------------------------------------
    axL = axes[1, 0]
    axL.loglog(
        masses_sorted, C_e_v2,
        color=color_e, linestyle=ls_e, linewidth=2.2, label="Electrons"
    )
    axL.loglog(
        masses_sorted, C_SD_v2,
        color=color_sd, linestyle=ls_sd, linewidth=2.2, label="Nucleons-SD"
    )
    axL.loglog(
        masses_sorted, C_SI_v2,
        color=color_si, linestyle=ls_si, linewidth=2.6, label="Nucleons-SI"
    )
    if not show_geometric_top_left_only:
        axL.loglog(
            masses_sorted, C_geo,
            color=color_geo, linestyle=ls_geo, linewidth=1.2, alpha=0.35, label="Geometric"
        )
    axL.set_ylabel(r"$C\ [{\rm s}^{-1}]$", fontsize=12)
    axL.set_title(r"$v^2$ dependent: absolute rates", fontsize=13)
    axL.grid(True, which="both", alpha=0.3)
    axL.legend(fontsize=9, loc="best")

    axR = axes[1, 1]
    mask_v2_e = _mask_pos(ratio_e_v2)
    mask_v2_sd = _mask_pos(ratio_SD_v2)
    mask_v2_si = _mask_pos(ratio_SI_v2)

    axR.semilogx(
        masses_sorted[mask_v2_e], ratio_e_v2[mask_v2_e],
        color=color_e, linestyle=ls_e, linewidth=2.2, label="Electrons"
    )
    axR.semilogx(
        masses_sorted[mask_v2_sd], ratio_SD_v2[mask_v2_sd],
        color=color_sd, linestyle=ls_sd, linewidth=2.2, label="Nucleons-SD"
    )
    axR.semilogx(
        masses_sorted[mask_v2_si], ratio_SI_v2[mask_v2_si],
        color=color_si, linestyle=ls_si, linewidth=2.6, label="Nucleons-SI"
    )
    axR.axhline(
        1.0,
        color="black",
        linestyle=":",
        linewidth=1.1,
        alpha=0.8,
        label="constant reference"
    )
    axR.set_ylabel(r"$C_{v^2}/C_{\rm constant}$", fontsize=12)
    axR.set_title(r"$v^2$ relative to constant", fontsize=13)
    axR.set_yscale("log")
    _set_ratio_panel_ylim(axR, [ratio_e_v2, ratio_SD_v2, ratio_SI_v2], top_row=False)
    axR.grid(True, which="both", alpha=0.3)
    axR.legend(fontsize=9, loc="best")

    # ------------------------------------------------------------
    # Row 3: q2
    # ------------------------------------------------------------
    axL = axes[2, 0]
    axL.loglog(
        masses_sorted, C_e_q2,
        color=color_e, linestyle=ls_e, linewidth=2.2, label="Electrons"
    )
    axL.loglog(
        masses_sorted, C_SD_q2,
        color=color_sd, linestyle=ls_sd, linewidth=2.2, label="Nucleons-SD"
    )
    axL.loglog(
        masses_sorted, C_SI_q2,
        color=color_si, linestyle=ls_si, linewidth=2.6, label="Nucleons-SI"
    )
    if not show_geometric_top_left_only:
        axL.loglog(
            masses_sorted, C_geo,
            color=color_geo, linestyle=ls_geo, linewidth=1.2, alpha=0.35, label="Geometric"
        )
    axL.set_xlabel(r"$m_\chi\ [{\rm GeV}]$", fontsize=12)
    axL.set_ylabel(r"$C\ [{\rm s}^{-1}]$", fontsize=12)
    axL.set_title(r"$q^2$ dependent: absolute rates", fontsize=13)
    axL.grid(True, which="both", alpha=0.3)
    axL.legend(fontsize=9, loc="best")

    axR = axes[2, 1]
    mask_q2_e = _mask_pos(ratio_e_q2)
    mask_q2_sd = _mask_pos(ratio_SD_q2)
    mask_q2_si = _mask_pos(ratio_SI_q2)

    axR.semilogx(
        masses_sorted[mask_q2_e], ratio_e_q2[mask_q2_e],
        color=color_e, linestyle=ls_e, linewidth=2.2, label="Electrons"
    )
    axR.semilogx(
        masses_sorted[mask_q2_sd], ratio_SD_q2[mask_q2_sd],
        color=color_sd, linestyle=ls_sd, linewidth=2.2, label="Nucleons-SD"
    )
    axR.semilogx(
        masses_sorted[mask_q2_si], ratio_SI_q2[mask_q2_si],
        color=color_si, linestyle=ls_si, linewidth=2.6, label="Nucleons-SI"
    )
    axR.axhline(
        1.0,
        color="black",
        linestyle=":",
        linewidth=1.1,
        alpha=0.8,
        label="constant reference"
    )
    axR.set_xlabel(r"$m_\chi\ [{\rm GeV}]$", fontsize=12)
    axR.set_ylabel(r"$C_{q^2}/C_{\rm constant}$", fontsize=12)
    axR.set_title(r"$q^2$ relative to constant", fontsize=13)
    axR.set_yscale("log")
    _set_ratio_panel_ylim(axR, [ratio_e_q2, ratio_SD_q2, ratio_SI_q2], top_row=False)
    axR.grid(True, which="both", alpha=0.3)
    axR.legend(fontsize=9, loc="best")

    if add_figure_note:
        fig.text(
            0.5, 0.008,
            "SD curves use the verified-only SD dataset. Electron baseline uses T=0 (thermal motion disabled).",
            ha="center",
            va="bottom",
            fontsize=9
        )

    plt.tight_layout(rect=(0.0, 0.02, 1.0, 0.98))
    fig.savefig(fig_png, dpi=300, bbox_inches="tight")
    fig.savefig(fig_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"[INFO] Combined operator-grid figure saved to: {fig_png}")
    print(f"[INFO] Combined operator-grid figure saved to: {fig_pdf}")

    return {
        "DM_masses": masses_sorted,
        "C_geo": C_geo,
        "C_SI_constant": C_SI_constant,
        "C_SI_v2": C_SI_v2,
        "C_SI_q2": C_SI_q2,
        "C_SD_constant": C_SD_constant,
        "C_SD_v2": C_SD_v2,
        "C_SD_q2": C_SD_q2,
        "C_e_constant": C_e_constant,
        "C_e_v2": C_e_v2,
        "C_e_q2": C_e_q2,
        "ratio_SI_v2": ratio_SI_v2,
        "ratio_SI_q2": ratio_SI_q2,
        "ratio_SD_v2": ratio_SD_v2,
        "ratio_SD_q2": ratio_SD_q2,
        "ratio_e_v2": ratio_e_v2,
        "ratio_e_q2": ratio_e_q2,
        "ratio_SDconst_over_SIconst": ratio_SDconst_over_SIconst,
        "ratio_econst_over_SIconst": ratio_econst_over_SIconst,
        "csv_path": csv_path,
        "fig_png": fig_png,
        "fig_pdf": fig_pdf,
        "log_path": log_path,
    }
def compute_one_mass_point_combined_operator_grid(task):
    """
    Worker for one mass point:
        - SI: constant / v2 / q2, T = 0
        - verified-only SD: constant / v2 / q2, T = 0
        - electron: constant / v2 / q2, T = 0
        - geometric reference
    """
    (
        earth_data,
        m,
        sigma_SI_p,
        sigma_SD_p,
        sigma_electron,
        rho_chi,
        u_max,
        n_u,
        v0,
        n_t_speed,
        n_t_mu,
        n_scatter_mu,
        n_scatter_phi,
        shell_step,
        u_grid_mode,
        R_earth_cm
    ) = task

    print(f"[worker:combined-operator-grid] start m = {m:.4g} GeV", flush=True)

    operator_map = {
        "constant": "constant",
        "v2": "v2_dependent",
        "q2": "q2_dependent",
    }

    out = {"m": float(m)}

    for op_tag, cross_type in operator_map.items():
        c_si = capture_rate_total(
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
            n_scatter_phi=n_scatter_phi,
            shell_step=shell_step,
            u_grid_mode=u_grid_mode
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
            v0=v0,
            include_thermal_targets=False,
            n_t_speed=n_t_speed,
            n_t_mu=n_t_mu,
            n_scatter_mu=n_scatter_mu,
            n_scatter_phi=n_scatter_phi,
            shell_step=shell_step,
            u_grid_mode=u_grid_mode
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
            v0=v0,
            include_thermal_targets=False,
            n_t_speed=n_t_speed,
            n_t_mu=n_t_mu,
            n_scatter_mu=n_scatter_mu,
            n_scatter_phi=n_scatter_phi,
            shell_step=shell_step,
            u_grid_mode=u_grid_mode
        )

        out[f"C_SI_{op_tag}"] = float(c_si)
        out[f"C_SD_{op_tag}"] = float(c_sd)
        out[f"C_e_{op_tag}"] = float(c_e)

    c_geo = capture_rate_geometric(
        DM_mass=m,
        R_earth_cm=R_earth_cm,
        rho_chi=rho_chi,
        v0=v0,
        vesc_surface=11.2
    )
    out["C_geo"] = float(c_geo)

    print(f"[worker:combined-operator-grid] done  m = {m:.4g} GeV", flush=True)
    return out


def save_combined_operator_grid_results_to_csv(
    data,
    output_csv="results/combined_SI_verifiedSD_electron_operator_grid_sigma1e-40.csv"
):
    """
    Save combined operator-grid results to CSV.
    """
    DM_masses = np.asarray(data["DM_masses"], dtype=float)

    C_geo = np.asarray(data["C_geo"], dtype=float)

    C_SI_constant = np.asarray(data["C_SI_constant"], dtype=float)
    C_SI_v2 = np.asarray(data["C_SI_v2"], dtype=float)
    C_SI_q2 = np.asarray(data["C_SI_q2"], dtype=float)

    C_SD_constant = np.asarray(data["C_SD_constant"], dtype=float)
    C_SD_v2 = np.asarray(data["C_SD_v2"], dtype=float)
    C_SD_q2 = np.asarray(data["C_SD_q2"], dtype=float)

    C_e_constant = np.asarray(data["C_e_constant"], dtype=float)
    C_e_v2 = np.asarray(data["C_e_v2"], dtype=float)
    C_e_q2 = np.asarray(data["C_e_q2"], dtype=float)

    ratio_SI_v2 = np.divide(
        C_SI_v2,
        C_SI_constant,
        out=np.full_like(C_SI_v2, np.nan, dtype=float),
        where=(C_SI_constant > 0.0)
    )
    ratio_SI_q2 = np.divide(
        C_SI_q2,
        C_SI_constant,
        out=np.full_like(C_SI_q2, np.nan, dtype=float),
        where=(C_SI_constant > 0.0)
    )

    ratio_SD_v2 = np.divide(
        C_SD_v2,
        C_SD_constant,
        out=np.full_like(C_SD_v2, np.nan, dtype=float),
        where=(C_SD_constant > 0.0)
    )
    ratio_SD_q2 = np.divide(
        C_SD_q2,
        C_SD_constant,
        out=np.full_like(C_SD_q2, np.nan, dtype=float),
        where=(C_SD_constant > 0.0)
    )

    ratio_e_v2 = np.divide(
        C_e_v2,
        C_e_constant,
        out=np.full_like(C_e_v2, np.nan, dtype=float),
        where=(C_e_constant > 0.0)
    )
    ratio_e_q2 = np.divide(
        C_e_q2,
        C_e_constant,
        out=np.full_like(C_e_q2, np.nan, dtype=float),
        where=(C_e_constant > 0.0)
    )

    out_dir = os.path.dirname(output_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "m_chi_GeV",
            "C_geo",

            "C_SI_constant",
            "C_SI_v2",
            "C_SI_q2",

            "C_SD_constant",
            "C_SD_v2",
            "C_SD_q2",

            "C_e_constant",
            "C_e_v2",
            "C_e_q2",

            "ratio_SI_v2_over_constant",
            "ratio_SI_q2_over_constant",
            "ratio_SD_v2_over_constant",
            "ratio_SD_q2_over_constant",
            "ratio_e_v2_over_constant",
            "ratio_e_q2_over_constant",
        ])

        for i in range(len(DM_masses)):
            writer.writerow([
                f"{DM_masses[i]:.8e}",
                f"{C_geo[i]:.8e}",

                f"{C_SI_constant[i]:.8e}",
                f"{C_SI_v2[i]:.8e}",
                f"{C_SI_q2[i]:.8e}",

                f"{C_SD_constant[i]:.8e}",
                f"{C_SD_v2[i]:.8e}",
                f"{C_SD_q2[i]:.8e}",

                f"{C_e_constant[i]:.8e}",
                f"{C_e_v2[i]:.8e}",
                f"{C_e_q2[i]:.8e}",

                f"{ratio_SI_v2[i]:.8e}" if np.isfinite(ratio_SI_v2[i]) else "nan",
                f"{ratio_SI_q2[i]:.8e}" if np.isfinite(ratio_SI_q2[i]) else "nan",
                f"{ratio_SD_v2[i]:.8e}" if np.isfinite(ratio_SD_v2[i]) else "nan",
                f"{ratio_SD_q2[i]:.8e}" if np.isfinite(ratio_SD_q2[i]) else "nan",
                f"{ratio_e_v2[i]:.8e}" if np.isfinite(ratio_e_v2[i]) else "nan",
                f"{ratio_e_q2[i]:.8e}" if np.isfinite(ratio_e_q2[i]) else "nan",
            ])

    print(f"[INFO] Combined operator-grid CSV saved to: {output_csv}")


def _set_ratio_panel_ylim(ax, arrays, top_row=False):
    """
    Helper for ratio-panel y-limits.
    """
    if top_row:
        ax.set_ylim(1e-1, 1e1)
        return

    vals = []
    for arr in arrays:
        arr = np.asarray(arr, dtype=float)
        mask = np.isfinite(arr) & (arr > 0.0)
        if np.any(mask):
            vals.append(arr[mask])

    if not vals:
        ax.set_ylim(1e-6, 1e1)
        return

    vals = np.concatenate(vals)
    vmin = np.nanmin(vals)
    vmax = np.nanmax(vals)

    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin <= 0.0:
        ax.set_ylim(1e-6, 1e1)
        return

    y_min = 10.0 ** np.floor(np.log10(vmin) - 0.3)
    y_max = 10.0 ** np.ceil(np.log10(vmax) + 0.3)

    if y_min == y_max:
        y_min /= 10.0
        y_max *= 10.0

    ax.set_ylim(y_min, y_max)


def compute_one_mass_point_combined_thermal_grid(task):
    """
    Worker for one mass point:
        - SI / verified-only SD / electron
        - operators: constant / v2 / q2
        - compute both T!=0 and T=0
    """
    (
        earth_data,
        m,
        sigma_SI_p,
        sigma_SD_p,
        sigma_electron,
        rho_chi,
        u_max,
        n_u,
        v0,
        n_t_speed,
        n_t_mu,
        n_scatter_mu,
        n_scatter_phi,
        shell_step,
        u_grid_mode
    ) = task

    print(f"[worker:combined-thermal-grid] start m = {m:.4g} GeV", flush=True)

    operator_map = {
        "constant": "constant",
        "v2": "v2_dependent",
        "q2": "q2_dependent",
    }

    out = {"m": float(m)}

    for op_tag, cross_type in operator_map.items():
        # SI
        out[f"C_SI_{op_tag}_T"] = float(capture_rate_total(
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
            n_scatter_phi=n_scatter_phi,
            shell_step=shell_step,
            u_grid_mode=u_grid_mode
        ))
        out[f"C_SI_{op_tag}_0"] = float(capture_rate_total(
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
            n_scatter_phi=n_scatter_phi,
            shell_step=shell_step,
            u_grid_mode=u_grid_mode
        ))

        # SD
        out[f"C_SD_{op_tag}_T"] = float(capture_rate_total(
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
            n_scatter_phi=n_scatter_phi,
            shell_step=shell_step,
            u_grid_mode=u_grid_mode
        ))
        out[f"C_SD_{op_tag}_0"] = float(capture_rate_total(
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
            n_scatter_phi=n_scatter_phi,
            shell_step=shell_step,
            u_grid_mode=u_grid_mode
        ))

        # electron
        out[f"C_e_{op_tag}_T"] = float(capture_rate_total(
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
            n_scatter_phi=n_scatter_phi,
            shell_step=shell_step,
            u_grid_mode=u_grid_mode
        ))
        out[f"C_e_{op_tag}_0"] = float(capture_rate_total(
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
            n_scatter_phi=n_scatter_phi,
            shell_step=shell_step,
            u_grid_mode=u_grid_mode
        ))

    print(f"[worker:combined-thermal-grid] done  m = {m:.4g} GeV", flush=True)
    return out
def plot_combined_SI_verified_only_SD_electron_thermal_grid(
    earth_data,
    DM_masses,
    sigma_SI_p=1e-40,
    sigma_SD_p=1e-40,
    sigma_electron=1e-40,
    output_root=".",
    u_max=800.0,
    n_u=80,
    rho_chi=RHO_CHI_DEFAULT,
    v0=V0_DEFAULT,
    n_t_speed=10,
    n_t_mu=10,
    n_scatter_mu=8,
    n_scatter_phi=12,
    shell_step=1,
    u_grid_mode="log",
    max_workers=None
):
    from concurrent.futures import ProcessPoolExecutor, as_completed

    expected_labels = ("Al", "H", "Na", "Si29")
    validate_verified_only_sd_baseline(
        earth_data=earth_data,
        expected_labels=expected_labels
    )

    if max_workers is None:
        max_workers = max(1, os.cpu_count() - 1)

    DM_masses = np.asarray(DM_masses, dtype=float)

    sigma_tag = f"{sigma_SI_p:.0e}".replace("+", "")
    figures_dir = os.path.join(output_root, "figures")
    os.makedirs(figures_dir, exist_ok=True)

    fig_png = os.path.join(
        figures_dir,
        f"combined_SI_verifiedSD_electron_thermal_grid_sigma{sigma_tag}.png"
    )
    fig_pdf = os.path.join(
        figures_dir,
        f"combined_SI_verifiedSD_electron_thermal_grid_sigma{sigma_tag}.pdf"
    )

    tasks = []
    for m in DM_masses:
        tasks.append((
            earth_data,
            float(m),
            sigma_SI_p,
            sigma_SD_p,
            sigma_electron,
            rho_chi,
            u_max,
            n_u,
            v0,
            n_t_speed,
            n_t_mu,
            n_scatter_mu,
            n_scatter_phi,
            shell_step,
            u_grid_mode
        ))

    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(compute_one_mass_point_combined_thermal_grid, task)
            for task in tasks
        ]
        for i, fut in enumerate(as_completed(futures), 1):
            res = fut.result()
            print(f"  done {i}/{len(futures)} : m_chi = {res['m']:.4g} GeV")
            results.append(res)

    results.sort(key=lambda x: x["m"])
    masses = np.array([r["m"] for r in results], dtype=float)

    def arr(key):
        return np.array([r[key] for r in results], dtype=float)

    def ratio(num, den):
        return np.divide(
            num, den,
            out=np.full_like(num, np.nan, dtype=float),
            where=(den > 0.0)
        )

    C = {}
    for op in ("constant", "v2", "q2"):
        C[f"SI_{op}_T"] = arr(f"C_SI_{op}_T")
        C[f"SI_{op}_0"] = arr(f"C_SI_{op}_0")
        C[f"SD_{op}_T"] = arr(f"C_SD_{op}_T")
        C[f"SD_{op}_0"] = arr(f"C_SD_{op}_0")
        C[f"e_{op}_T"] = arr(f"C_e_{op}_T")
        C[f"e_{op}_0"] = arr(f"C_e_{op}_0")

    fig, axes = plt.subplots(3, 2, figsize=(14, 16), sharex="col")
    fig.suptitle(
        r"Combined SI / verified-only SD / electron thermal comparison",
        fontsize=15,
        y=0.995
    )

    color_e = "red"
    color_sd = "green"
    color_si = "blue"

    ls_e = "-"
    ls_sd = "--"
    ls_si = ":"

    row_meta = [
        ("constant", r"constant"),
        ("v2", r"$v_{\rm rel}^2$ dependent"),
        ("q2", r"$q^2$ dependent"),
    ]

    for irow, (op, title_text) in enumerate(row_meta):
        axL = axes[irow, 0]
        axR = axes[irow, 1]

        # inset-inline-start: C(T != 0)
        axL.loglog(masses, C[f"e_{op}_T"], color=color_e, linestyle=ls_e, linewidth=2.2, label="Electrons")
        axL.loglog(masses, C[f"SD_{op}_T"], color=color_sd, linestyle=ls_sd, linewidth=2.2, label="Nucleons-SD")
        axL.loglog(masses, C[f"SI_{op}_T"], color=color_si, linestyle=ls_si, linewidth=2.6, label="Nucleons-SI")
        axL.set_ylabel(r"$C\ [{\rm s}^{-1}]$", fontsize=12)
        axL.set_title(f"{title_text}: capture rates", fontsize=13)
        axL.grid(True, which="both", alpha=0.3)
        axL.legend(fontsize=9, loc="best")

        # inset-inline-end: C(T)/C(T=0)
        R_e = ratio(C[f"e_{op}_T"], C[f"e_{op}_0"])
        R_sd = ratio(C[f"SD_{op}_T"], C[f"SD_{op}_0"])
        R_si = ratio(C[f"SI_{op}_T"], C[f"SI_{op}_0"])

        axR.semilogx(masses, R_e, color=color_e, linestyle=ls_e, linewidth=2.2, label="Electrons")
        axR.semilogx(masses, R_sd, color=color_sd, linestyle=ls_sd, linewidth=2.2, label="Nucleons-SD")
        axR.semilogx(masses, R_si, color=color_si, linestyle=ls_si, linewidth=2.6, label="Nucleons-SI")
        axR.axhline(1.0, color="black", linestyle=":", linewidth=1.1, alpha=0.8, label=r"$T(r)=0$ reference")
        axR.set_ylabel(r"$C(T)/C(T=0)$", fontsize=12)
        axR.set_title(f"{title_text}: thermal correction", fontsize=13)
        axR.set_yscale("log")
        axR.grid(True, which="both", alpha=0.3)
        axR.legend(fontsize=9, loc="best")

    axes[2, 0].set_xlabel(r"$m_\chi\ [{\rm GeV}]$", fontsize=12)
    axes[2, 1].set_xlabel(r"$m_\chi\ [{\rm GeV}]$", fontsize=12)

    fig.text(
        0.5, 0.008,
        "SD curves use the verified-only SD dataset. Electron finite-T is included in the free thermal-electron approximation.",
        ha="center",
        va="bottom",
        fontsize=9
    )

    plt.tight_layout(rect=(0.0, 0.02, 1.0, 0.98))
    fig.savefig(fig_png, dpi=300, bbox_inches="tight")
    fig.savefig(fig_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"[INFO] Combined thermal-grid figure saved to: {fig_png}")
    print(f"[INFO] Combined thermal-grid figure saved to: {fig_pdf}")

    return {"fig_png": fig_png, "fig_pdf": fig_pdf}

def test_electron_thermal_convergence(
    earth_data,
    masses,
    sigma_electron=1e-40,
    cross_section_types=("constant", "v2_dependent", "q2_dependent"),
    rho_chi=RHO_CHI_DEFAULT,
    v0=V0_DEFAULT,
    u_max=800.0,
    top_note=True
):
    """
    One-parameter-at-a-time convergence test for electron thermal capture.

    Purpose
    -------
    Diagnose whether the electron thermal enhancement
        C_e(T) / C_e(T=0)
    is numerically stable.

    Parameters scanned
    ------------------
    - n_u
    - n_t_speed
    - n_t_mu
    - n_scatter_mu
    - n_scatter_phi
    - shell_step

    Notes
    -----
    This assumes:
        - electron finite-T has been enabled in omega_minus_at_shell(...)
        - scattering_type = "electron"
        - T(r) != 0 and T = 0 are both computed
    """
    import time
    import numpy as np

    masses = np.asarray(masses, dtype=float)

    if top_note:
        print("\n" + "=" * 100)
        print("Electron thermal convergence test")
        print("=" * 100)
        print("Mode assumptions:")
        print("  - electron only")
        print("  - compare T(r)!=0 vs T=0")
        print("  - free thermal-electron approximation")
        print(f"  - cross_section_types = {cross_section_types}")
        print("=" * 100)

    # ------------------------------------------------------------
    # Baseline settings
    # ------------------------------------------------------------
    base = {
        "n_u": 30,
        "n_t_speed": 4,
        "n_t_mu": 4,
        "n_scatter_mu": 4,
        "n_scatter_phi": 6,
        "shell_step": 2,
        "u_grid_mode": "log",
    }

    # ------------------------------------------------------------
    # One-parameter scans
    # ------------------------------------------------------------
    scans = {
        "n_u": [20, 30, 40, 60],
        "n_t_speed": [3, 4, 6, 8],
        "n_t_mu": [3, 4, 6, 8],
        "n_scatter_mu": [3, 4, 6, 8],
        "n_scatter_phi": [4, 6, 8, 12],
        "shell_step": [4, 2, 1],
    }

    all_results = {}

    for cross_type in cross_section_types:
        all_results[cross_type] = {}

        print("\n" + "#" * 100)
        print(f"Operator / cross_section_type = {cross_type}")
        print("#" * 100)

        for m in masses:
            m = float(m)
            all_results[cross_type][m] = {}

            print("\n" + "-" * 100)
            print(f"Electron thermal convergence at m_chi = {m:.6g} GeV")
            print("-" * 100)

            for param_name, values in scans.items():
                print(f"\nScan: {param_name}")

                scan_rows = []
                prev_ratio = None
                prev_cT = None
                prev_c0 = None

                for val in values:
                    s = dict(base)
                    s[param_name] = val

                    print(
                        f"[start] op={cross_type}, "
                        f"m={m:.6g}, "
                        f"scan={param_name}, val={val}, "
                        f"settings=("
                        f"n_u={s['n_u']}, "
                        f"n_t_speed={s['n_t_speed']}, "
                        f"n_t_mu={s['n_t_mu']}, "
                        f"n_scatter_mu={s['n_scatter_mu']}, "
                        f"n_scatter_phi={s['n_scatter_phi']}, "
                        f"shell_step={s['shell_step']}, "
                        f"u_grid_mode={s['u_grid_mode']}"
                        f")",
                        flush=True
                    )

                    t0 = time.time()

                    c_e_T = capture_rate_total(
                        earth_data=earth_data,
                        DM_mass=m,
                        sigma_electron=sigma_electron,
                        scattering_type="electron",
                        cross_section_type=cross_type,
                        rho_chi=rho_chi,
                        u_max=u_max,
                        n_u=s["n_u"],
                        v0=v0,
                        include_thermal_targets=True,
                        n_t_speed=s["n_t_speed"],
                        n_t_mu=s["n_t_mu"],
                        n_scatter_mu=s["n_scatter_mu"],
                        n_scatter_phi=s["n_scatter_phi"],
                        u_grid_mode=s["u_grid_mode"],
                        shell_step=s["shell_step"]
                    )

                    c_e_0 = capture_rate_total(
                        earth_data=earth_data,
                        DM_mass=m,
                        sigma_electron=sigma_electron,
                        scattering_type="electron",
                        cross_section_type=cross_type,
                        rho_chi=rho_chi,
                        u_max=u_max,
                        n_u=s["n_u"],
                        v0=v0,
                        include_thermal_targets=False,
                        n_t_speed=s["n_t_speed"],
                        n_t_mu=s["n_t_mu"],
                        n_scatter_mu=s["n_scatter_mu"],
                        n_scatter_phi=s["n_scatter_phi"],
                        u_grid_mode=s["u_grid_mode"],
                        shell_step=s["shell_step"]
                    )

                    dt = time.time() - t0

                    c_e_T = float(c_e_T)
                    c_e_0 = float(c_e_0)
                    ratio_e = float(c_e_T / c_e_0) if c_e_0 > 0.0 else float("nan")

                    if prev_ratio is None or not np.isfinite(prev_ratio) or prev_ratio == 0.0:
                        delta_ratio = np.nan
                    else:
                        delta_ratio = (ratio_e - prev_ratio) / prev_ratio

                    if prev_cT is None or prev_cT == 0.0:
                        delta_cT = np.nan
                    else:
                        delta_cT = (c_e_T - prev_cT) / prev_cT

                    if prev_c0 is None or prev_c0 == 0.0:
                        delta_c0 = np.nan
                    else:
                        delta_c0 = (c_e_0 - prev_c0) / prev_c0

                    print(
                        f"[done ] {param_name}={val!s:>3s} | "
                        f"C_e(T)={c_e_T:.6e}, "
                        f"C_e(0)={c_e_0:.6e}, "
                        f"ratio={ratio_e:.6e} | "
                        f"dC_T={delta_cT:+.3e}, "
                        f"dC_0={delta_c0:+.3e}, "
                        f"dratio={delta_ratio:+.3e} | "
                        f"time={dt:.1f}s",
                        flush=True
                    )

                    row = {
                        "param_name": param_name,
                        "param_value": val,
                        "settings": dict(s),
                        "C_e_T": c_e_T,
                        "C_e_0": c_e_0,
                        "ratio": ratio_e,
                        "delta_C_e_T": delta_cT,
                        "delta_C_e_0": delta_c0,
                        "delta_ratio": delta_ratio,
                        "time_sec": dt,
                    }
                    scan_rows.append(row)

                    prev_ratio = ratio_e
                    prev_cT = c_e_T
                    prev_c0 = c_e_0

                all_results[cross_type][m][param_name] = scan_rows

    return all_results
def summarize_electron_thermal_convergence(
    conv_results,
    ratio_tol=0.05,
    rate_tol=0.05,
    output_txt=None,
    verbose=True
):
    """
    Summarize the output of test_electron_thermal_convergence(...).

    Parameters
    ----------
    conv_results : dict
        Output from test_electron_thermal_convergence(...)

    ratio_tol : float
        Allowed relative deviation threshold for the thermal ratio
            C_e(T) / C_e(0)
        when comparing lower-resolution settings to the highest-resolution
        setting in each one-parameter scan.

    rate_tol : float
        Allowed relative deviation threshold for both
            C_e(T)
            C_e(0)

    output_txt : str or None
        If not None, save the summary to a text file.

    verbose : bool
        If True, print the summary.

    Returns
    -------
    summary : dict
        Structured summary with:
            - per_scan diagnostics
            - per_scan recommended setting
            - global recommended settings
    """
    import os
    import numpy as np

    def rel_dev(x, x_ref):
        x = float(x)
        x_ref = float(x_ref)
        if not np.isfinite(x) or not np.isfinite(x_ref):
            return np.nan
        if x_ref == 0.0:
            return np.nan
        return abs((x - x_ref) / x_ref)

    def classify_status(last_ratio_dev, last_cT_dev, last_c0_dev, ratio_tol, rate_tol):
        """
        Simple quality label based on the last-step deviation to the highest-resolution point.
        """
        vals = [last_ratio_dev, last_cT_dev, last_c0_dev]
        if any((not np.isfinite(v)) for v in vals):
            return "UNRESOLVED"
        if last_ratio_dev <= ratio_tol and last_cT_dev <= rate_tol and last_c0_dev <= rate_tol:
            return "OK"
        if last_ratio_dev <= 2.0 * ratio_tol and last_cT_dev <= 2.0 * rate_tol and last_c0_dev <= 2.0 * rate_tol:
            return "MARGINAL"
        return "UNSTABLE"

    def choose_recommended_index(rows, ratio_tol, rate_tol):
        """
        Choose the earliest scan index i such that all rows from i onward
        stay within tolerance of the highest-resolution row (the last row).
        """
        final_row = rows[-1]
        final_ratio = final_row["ratio"]
        final_cT = final_row["C_e_T"]
        final_c0 = final_row["C_e_0"]

        n = len(rows)
        ratio_dev_to_final = [rel_dev(r["ratio"], final_ratio) for r in rows]
        cT_dev_to_final = [rel_dev(r["C_e_T"], final_cT) for r in rows]
        c0_dev_to_final = [rel_dev(r["C_e_0"], final_c0) for r in rows]

        for i in range(n):
            ok_tail = True
            for j in range(i, n):
                r_ok = np.isfinite(ratio_dev_to_final[j]) and (ratio_dev_to_final[j] <= ratio_tol)
                t_ok = np.isfinite(cT_dev_to_final[j]) and (cT_dev_to_final[j] <= rate_tol)
                z_ok = np.isfinite(c0_dev_to_final[j]) and (c0_dev_to_final[j] <= rate_tol)
                if not (r_ok and t_ok and z_ok):
                    ok_tail = False
                    break
            if ok_tail:
                return i

        return n - 1

    lines = []
    lines.append("=" * 100)
    lines.append("Electron thermal convergence summary")
    lines.append("=" * 100)
    lines.append(f"ratio_tol = {ratio_tol:.3e}")
    lines.append(f"rate_tol  = {rate_tol:.3e}")
    lines.append("Criterion:")
    lines.append("  A setting is accepted if, within a one-parameter scan,")
    lines.append("  its C_e(T), C_e(0), and C_e(T)/C_e(0) remain within tolerance")
    lines.append("  of the highest-resolution setting for that scan.")
    lines.append("")

    summary = {
        "ratio_tol": float(ratio_tol),
        "rate_tol": float(rate_tol),
        "per_scan": {},
        "global_recommended_settings": {}
    }

    # Store scan value order for each parameter, so we can build a global recommendation
    param_orders = {}
    # Track the strictest required index across all masses/operators
    global_required_index = {}

    cross_types = list(conv_results.keys())

    for cross_type in cross_types:
        summary["per_scan"][cross_type] = {}
        lines.append("#" * 100)
        lines.append(f"Operator / cross_section_type = {cross_type}")
        lines.append("#" * 100)

        masses = sorted(conv_results[cross_type].keys(), key=float)

        for m in masses:
            summary["per_scan"][cross_type][m] = {}
            lines.append("")
            lines.append("-" * 100)
            lines.append(f"m_chi = {float(m):.6g} GeV")
            lines.append("-" * 100)

            scans_for_mass = conv_results[cross_type][m]

            for param_name, rows in scans_for_mass.items():
                if len(rows) == 0:
                    continue

                values_in_order = [r["param_value"] for r in rows]
                param_orders.setdefault(param_name, values_in_order)

                final_row = rows[-1]
                final_ratio = final_row["ratio"]
                final_cT = final_row["C_e_T"]
                final_c0 = final_row["C_e_0"]

                # Deviations of every scan point to the highest-resolution point
                ratio_dev_to_final = np.array([rel_dev(r["ratio"], final_ratio) for r in rows], dtype=float)
                cT_dev_to_final = np.array([rel_dev(r["C_e_T"], final_cT) for r in rows], dtype=float)
                c0_dev_to_final = np.array([rel_dev(r["C_e_0"], final_c0) for r in rows], dtype=float)

                rec_idx = choose_recommended_index(rows, ratio_tol, rate_tol)
                rec_row = rows[rec_idx]
                rec_val = rec_row["param_value"]

                global_required_index[param_name] = max(
                    global_required_index.get(param_name, 0),
                    rec_idx
                )

                # Last non-final point deviation to final, useful quick diagnostic
                if len(rows) >= 2:
                    last_ratio_dev = ratio_dev_to_final[-2]
                    last_cT_dev = cT_dev_to_final[-2]
                    last_c0_dev = c0_dev_to_final[-2]
                else:
                    last_ratio_dev = np.nan
                    last_cT_dev = np.nan
                    last_c0_dev = np.nan

                status = classify_status(
                    last_ratio_dev=last_ratio_dev,
                    last_cT_dev=last_cT_dev,
                    last_c0_dev=last_c0_dev,
                    ratio_tol=ratio_tol,
                    rate_tol=rate_tol
                )

                summary["per_scan"][cross_type][m][param_name] = {
                    "recommended_index": int(rec_idx),
                    "recommended_value": rec_val,
                    "final_ratio": float(final_ratio),
                    "final_C_e_T": float(final_cT),
                    "final_C_e_0": float(final_c0),
                    "max_ratio_dev_to_final": float(np.nanmax(ratio_dev_to_final)),
                    "max_C_e_T_dev_to_final": float(np.nanmax(cT_dev_to_final)),
                    "max_C_e_0_dev_to_final": float(np.nanmax(c0_dev_to_final)),
                    "last_ratio_dev_to_final": float(last_ratio_dev) if np.isfinite(last_ratio_dev) else np.nan,
                    "last_C_e_T_dev_to_final": float(last_cT_dev) if np.isfinite(last_cT_dev) else np.nan,
                    "last_C_e_0_dev_to_final": float(last_c0_dev) if np.isfinite(last_c0_dev) else np.nan,
                    "status": status,
                    "scan_values_in_order": list(values_in_order),
                }

                lines.append(
                    f"{param_name:>14s} | "
                    f"recommended = {str(rec_val):>6s} | "
                    f"status = {status:<10s} | "
                    f"final ratio = {final_ratio:.6e} | "
                    f"last devs: "
                    f"dratio = {last_ratio_dev:+.3e}, "
                    f"dC_T = {last_cT_dev:+.3e}, "
                    f"dC_0 = {last_c0_dev:+.3e}"
                )

    # ------------------------------------------------------------
    # Global recommended settings
    # ------------------------------------------------------------
    lines.append("")
    lines.append("=" * 100)
    lines.append("Global recommended settings")
    lines.append("=" * 100)

    global_settings = {}

    for param_name, ordered_values in param_orders.items():
        idx_req = global_required_index.get(param_name, len(ordered_values) - 1)
        idx_req = min(max(idx_req, 0), len(ordered_values) - 1)
        val_req = ordered_values[idx_req]

        global_settings[param_name] = {
            "required_index": int(idx_req),
            "recommended_value": val_req,
            "scan_order": list(ordered_values),
        }

        lines.append(
            f"{param_name:>14s} : recommended = {str(val_req):>6s}   "
            f"(scan order = {ordered_values}, required_index = {idx_req})"
        )

    summary["global_recommended_settings"] = global_settings

    lines.append("")
    lines.append("Interpretation note:")
    lines.append("  - For n_u, n_t_speed, n_t_mu, n_scatter_mu, n_scatter_phi,")
    lines.append("    larger values usually mean higher resolution.")
    lines.append("  - For shell_step, smaller values mean finer radial resolution;")
    lines.append("    the scan order is used exactly as provided by the original test.")
    lines.append("")

    # ------------------------------------------------------------
    # Save or print
    # ------------------------------------------------------------
    text = "\n".join(lines)

    if verbose:
        print(text)

    if output_txt is not None:
        out_dir = os.path.dirname(output_txt)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(output_txt, "w", encoding="utf-8") as f:
            f.write(text)
        print(f"[INFO] Electron thermal convergence summary saved to: {output_txt}")

    return summary
def plot_electron_thermal_convergence_summary(
    conv_results,
    ratio_tol=0.05,
    rate_tol=0.05,
    metric="ratio",
    output_root=".",
    save_prefix="electron_thermal_convergence_summary",
    show=False
):
    """
    Plot convergence-summary figures from test_electron_thermal_convergence(...).

    Parameters
    ----------
    conv_results : dict
        Output from test_electron_thermal_convergence(...)

    ratio_tol : float
        Tolerance for the thermal ratio
            R = C_e(T) / C_e(0)

    rate_tol : float
        Tolerance for the absolute rates
            C_e(T), C_e(0)

    metric : str
        Which convergence metric to visualize.
        Allowed:
            - "ratio" : relative deviation of C_e(T)/C_e(0)
            - "C_e_T" : relative deviation of C_e(T)
            - "C_e_0" : relative deviation of C_e(0)

    output_root : str
        Root output directory.

    save_prefix : str
        Prefix for saved figure names.

    show : bool
        If True, call plt.show() for each figure.

    Returns
    -------
    out : dict
        {
            cross_section_type: {
                "fig_png": ...,
                "fig_pdf": ...,
                "recommended": {param_name: recommended_value, ...}
            },
            ...
        }
    """
    import os
    import numpy as np
    import matplotlib.pyplot as plt

    allowed_metrics = {"ratio", "C_e_T", "C_e_0"}
    if metric not in allowed_metrics:
        raise ValueError(f"metric must be one of {allowed_metrics}, got {metric!r}")

    def rel_dev(x, x_ref):
        x = float(x)
        x_ref = float(x_ref)
        if not np.isfinite(x) or not np.isfinite(x_ref):
            return np.nan
        if x_ref == 0.0:
            return np.nan
        return abs((x - x_ref) / x_ref)

    def choose_recommended_index(rows, ratio_tol, rate_tol):
        """
        Choose earliest scan index i such that all rows from i onward
        remain within tolerance of the highest-resolution row.
        """
        final_row = rows[-1]
        final_ratio = final_row["ratio"]
        final_cT = final_row["C_e_T"]
        final_c0 = final_row["C_e_0"]

        n = len(rows)

        ratio_dev_to_final = [rel_dev(r["ratio"], final_ratio) for r in rows]
        cT_dev_to_final = [rel_dev(r["C_e_T"], final_cT) for r in rows]
        c0_dev_to_final = [rel_dev(r["C_e_0"], final_c0) for r in rows]

        for i in range(n):
            ok_tail = True
            for j in range(i, n):
                r_ok = np.isfinite(ratio_dev_to_final[j]) and (ratio_dev_to_final[j] <= ratio_tol)
                t_ok = np.isfinite(cT_dev_to_final[j]) and (cT_dev_to_final[j] <= rate_tol)
                z_ok = np.isfinite(c0_dev_to_final[j]) and (c0_dev_to_final[j] <= rate_tol)
                if not (r_ok and t_ok and z_ok):
                    ok_tail = False
                    break
            if ok_tail:
                return i

        return n - 1

    def build_dev_array(rows, metric):
        """
        Build relative-deviation array to the highest-resolution row.
        """
        final_row = rows[-1]

        if metric == "ratio":
            ref = final_row["ratio"]
            vals = [r["ratio"] for r in rows]
        elif metric == "C_e_T":
            ref = final_row["C_e_T"]
            vals = [r["C_e_T"] for r in rows]
        elif metric == "C_e_0":
            ref = final_row["C_e_0"]
            vals = [r["C_e_0"] for r in rows]
        else:
            raise ValueError("Invalid metric")

        return np.array([rel_dev(v, ref) for v in vals], dtype=float)

    scan_order = [
        "n_u",
        "n_t_speed",
        "n_t_mu",
        "n_scatter_mu",
        "n_scatter_phi",
        "shell_step",
    ]

    pretty_label = {
        "n_u": r"$n_u$",
        "n_t_speed": r"$n_{t,\mathrm{speed}}$",
        "n_t_mu": r"$n_{t,\mu}$",
        "n_scatter_mu": r"$n_{\mathrm{scatter},\mu}$",
        "n_scatter_phi": r"$n_{\mathrm{scatter},\phi}$",
        "shell_step": r"${\rm shell\_step}$",
    }

    if metric == "ratio":
        y_label = r"$\left|R-R_{\rm final}\right|/R_{\rm final}$"
        tol = float(ratio_tol)
        metric_title = r"ratio convergence: $R=C_e(T)/C_e(0)$"
    elif metric == "C_e_T":
        y_label = r"$\left|C_e(T)-C_{e,{\rm final}}(T)\right|/C_{e,{\rm final}}(T)$"
        tol = float(rate_tol)
        metric_title = r"finite-$T$ rate convergence: $C_e(T)$"
    else:  # C_e_0
        y_label = r"$\left|C_e(0)-C_{e,{\rm final}}(0)\right|/C_{e,{\rm final}}(0)$"
        tol = float(rate_tol)
        metric_title = r"$T=0$ rate convergence: $C_e(0)$"

    figures_dir = os.path.join(output_root, "figures")
    os.makedirs(figures_dir, exist_ok=True)

    out = {}

    for cross_type in conv_results.keys():
        masses = sorted(conv_results[cross_type].keys(), key=float)

        # ------------------------------------------------------------
        # Determine operator-level recommended settings
        # strictest required index across all masses
        # ------------------------------------------------------------
        recommended = {}
        for param_name in scan_order:
            max_idx = None
            values_ref = None

            for m in masses:
                if param_name not in conv_results[cross_type][m]:
                    continue
                rows = conv_results[cross_type][m][param_name]
                if len(rows) == 0:
                    continue
                idx = choose_recommended_index(rows, ratio_tol, rate_tol)
                if max_idx is None or idx > max_idx:
                    max_idx = idx
                    values_ref = [r["param_value"] for r in rows]

            if max_idx is not None and values_ref is not None:
                recommended[param_name] = values_ref[max_idx]

        # ------------------------------------------------------------
        # Plot
        # ------------------------------------------------------------
        fig, axes = plt.subplots(2, 3, figsize=(16, 9))
        axes = axes.flatten()

        fig.suptitle(
            f"Electron thermal convergence summary: {cross_type}, {metric_title}",
            fontsize=15,
            y=0.995
        )

        cmap = plt.get_cmap("viridis")
        colors = [cmap(i) for i in np.linspace(0.1, 0.9, len(masses))]

        for iax, param_name in enumerate(scan_order):
            ax = axes[iax]

            local_values = None
            recommended_x_idx = None
            recommended_value = recommended.get(param_name, None)

            for color, m in zip(colors, masses):
                if param_name not in conv_results[cross_type][m]:
                    continue

                rows = conv_results[cross_type][m][param_name]
                if len(rows) == 0:
                    continue

                scan_values = [r["param_value"] for r in rows]
                local_values = scan_values

                dev = build_dev_array(rows, metric)
                x = np.arange(len(scan_values))

                ax.plot(
                    x,
                    dev,
                    marker="o",
                    linewidth=2.0,
                    markersize=4,
                    color=color,
                    label=fr"$m_\chi={float(m):g}\,{{\rm GeV}}$"
                )

                if recommended_value is not None:
                    try:
                        recommended_x_idx = scan_values.index(recommended_value)
                    except ValueError:
                        pass

            ax.axhline(
                tol,
                color="black",
                linestyle=":",
                linewidth=1.2,
                label="tolerance" if iax == 0 else None
            )

            if recommended_x_idx is not None:
                ax.axvline(
                    recommended_x_idx,
                    color="tab:red",
                    linestyle="--",
                    linewidth=1.1,
                    alpha=0.8
                )

            ax.set_yscale("log")
            ax.set_title(pretty_label.get(param_name, param_name), fontsize=12)
            ax.set_ylabel(y_label, fontsize=11)

            if local_values is not None:
                ax.set_xticks(np.arange(len(local_values)))
                ax.set_xticklabels([str(v) for v in local_values])

            if iax >= 3:
                ax.set_xlabel("scan value", fontsize=11)

            ax.grid(True, which="both", alpha=0.3)

            if recommended_value is not None:
                ax.text(
                    0.03, 0.06,
                    f"recommended = {recommended_value}",
                    transform=ax.transAxes,
                    fontsize=9,
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
                )

        handles, labels = axes[0].get_legend_handles_labels()
        if handles:
            fig.legend(
                handles,
                labels,
                loc="upper center",
                ncol=min(5, len(labels)),
                fontsize=9,
                frameon=False,
                bbox_to_anchor=(0.5, 0.97)
            )

        plt.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))

        safe_metric = metric.replace("/", "_")
        fig_png = os.path.join(
            figures_dir,
            f"{save_prefix}_{cross_type}_{safe_metric}.png"
        )
        fig_pdf = os.path.join(
            figures_dir,
            f"{save_prefix}_{cross_type}_{safe_metric}.pdf"
        )

        fig.savefig(fig_png, dpi=300, bbox_inches="tight")
        fig.savefig(fig_pdf, bbox_inches="tight")

        if show:
            plt.show()
        else:
            plt.close(fig)

        print(f"[INFO] Electron convergence summary figure saved to: {fig_png}")
        print(f"[INFO] Electron convergence summary figure saved to: {fig_pdf}")

        out[cross_type] = {
            "fig_png": fig_png,
            "fig_pdf": fig_pdf,
            "recommended": recommended,
        }

    return out
def run_and_plot_electron_thermal_convergence_suite(
    earth_data,
    masses=(0.3, 1.0, 10.0, 100.0),
    sigma_electron=1e-40,
    cross_section_types=("constant", "v2_dependent", "q2_dependent"),
    rho_chi=RHO_CHI_DEFAULT,
    v0=V0_DEFAULT,
    u_max=800.0,
    ratio_tol=0.05,
    rate_tol=0.05,
    output_root=".",
    summary_txt="logs/electron_thermal_convergence_summary.txt",
    show_plots=False,
    top_note=True
):
    """
    Run the full electron thermal convergence suite:
        1. test_electron_thermal_convergence(...)
        2. summarize_electron_thermal_convergence(...)
        3. plot_electron_thermal_convergence_summary(...) for:
             - ratio
             - C_e_T
             - C_e_0

    Parameters
    ----------
    earth_data : dict
        Earth composition / shell data.

    masses : sequence of float
        DM masses [GeV] to test.

    sigma_electron : float
        Electron reference cross section [cm^2].

    cross_section_types : tuple/list of str
        e.g. ("constant", "v2_dependent", "q2_dependent")

    rho_chi, v0, u_max : float
        Standard capture inputs.

    ratio_tol : float
        Tolerance for
            C_e(T)/C_e(0)

    rate_tol : float
        Tolerance for
            C_e(T), C_e(0)

    output_root : str
        Root directory for figures/logs.

    summary_txt : str
        Path to text summary, relative or absolute.

    show_plots : bool
        If True, show figures interactively.

    top_note : bool
        Passed to test_electron_thermal_convergence(...)

    Returns
    -------
    out : dict
        {
            "conv_results": ...,
            "summary": ...,
            "fig_ratio": ...,
            "fig_C_e_T": ...,
            "fig_C_e_0": ...,
        }
    """
    import os
    import numpy as np

    print("\n" + "=" * 110)
    print("Running electron thermal convergence suite")
    print("=" * 110)

    masses = np.asarray(masses, dtype=float)

    # ------------------------------------------------------------
    # Resolve summary path
    # ------------------------------------------------------------
    if os.path.isabs(summary_txt):
        summary_txt_path = summary_txt
    else:
        summary_txt_path = os.path.join(output_root, summary_txt)

    # ------------------------------------------------------------
    # 1. Raw convergence test
    # ------------------------------------------------------------
    conv_results = test_electron_thermal_convergence(
        earth_data=earth_data,
        masses=masses,
        sigma_electron=sigma_electron,
        cross_section_types=cross_section_types,
        rho_chi=rho_chi,
        v0=v0,
        u_max=u_max,
        top_note=top_note
    )

    # ------------------------------------------------------------
    # 2. Text summary
    # ------------------------------------------------------------
    summary = summarize_electron_thermal_convergence(
        conv_results=conv_results,
        ratio_tol=ratio_tol,
        rate_tol=rate_tol,
        output_txt=summary_txt_path,
        verbose=True
    )

    # ------------------------------------------------------------
    # 3a. Plot ratio convergence summary
    # ------------------------------------------------------------
    fig_ratio = plot_electron_thermal_convergence_summary(
        conv_results=conv_results,
        ratio_tol=ratio_tol,
        rate_tol=rate_tol,
        metric="ratio",
        output_root=output_root,
        save_prefix="electron_thermal_conv_ratio",
        show=show_plots
    )

    # ------------------------------------------------------------
    # 3b. Plot C_e(T) convergence summary
    # ------------------------------------------------------------
    fig_C_e_T = plot_electron_thermal_convergence_summary(
        conv_results=conv_results,
        ratio_tol=ratio_tol,
        rate_tol=rate_tol,
        metric="C_e_T",
        output_root=output_root,
        save_prefix="electron_thermal_conv_CeT",
        show=show_plots
    )

    # ------------------------------------------------------------
    # 3c. Plot C_e(0) convergence summary
    # ------------------------------------------------------------
    fig_C_e_0 = plot_electron_thermal_convergence_summary(
        conv_results=conv_results,
        ratio_tol=ratio_tol,
        rate_tol=rate_tol,
        metric="C_e_0",
        output_root=output_root,
        save_prefix="electron_thermal_conv_Ce0",
        show=show_plots
    )

    # ------------------------------------------------------------
    # Print compact global recommendation
    # ------------------------------------------------------------
    print("\n" + "=" * 110)
    print("Compact global recommended settings")
    print("=" * 110)

    global_rec = summary.get("global_recommended_settings", {})
    for param_name, info in global_rec.items():
        print(
            f"{param_name:>14s} : "
            f"{str(info.get('recommended_value', None)):>6s}   "
            f"(required_index = {info.get('required_index', None)})"
        )

    print("=" * 110)
    print("Electron thermal convergence suite finished.")
    print("=" * 110)

    return {
        "conv_results": conv_results,
        "summary": summary,
        "fig_ratio": fig_ratio,
        "fig_C_e_T": fig_C_e_T,
        "fig_C_e_0": fig_C_e_0,
    }
def print_electron_thermal_recommended_production_block(
    summary_or_suite,
    function_name="plot_combined_SI_verified_only_SD_electron_thermal_grid",
    variable_name="combined_thermal_grid_production",
    earth_data_name="earth_data",
    masses_name="DM_masses_baseline",
    sigma_SI_p=1e-40,
    sigma_SD_p=1e-40,
    sigma_electron=1e-40,
    rho_chi_name="RHO_CHI_DEFAULT",
    v0_name="V0_DEFAULT",
    u_max=800.0,
    output_root=".",
    max_workers=10,
    print_result=True
):
    """
    Format a recommended production-call block from the output of
    summarize_electron_thermal_convergence(...) or
    run_and_plot_electron_thermal_convergence_suite(...).

    Parameters
    ----------
    summary_or_suite : dict
        Either:
          - the direct output of summarize_electron_thermal_convergence(...)
          - or the full suite output from
            run_and_plot_electron_thermal_convergence_suite(...),
            in which case this function will use summary_or_suite["summary"]

    function_name : str
        Target plotting function name to print.

    variable_name : str
        Variable name on the left-hand side of the printed assignment.

    earth_data_name : str
        Name of the earth_data variable to print.

    masses_name : str
        Name of the DM mass grid variable to print.

    sigma_SI_p, sigma_SD_p, sigma_electron : float
        Cross sections to print.

    rho_chi_name, v0_name : str
        Names of constants to print literally in the block.

    u_max : float
        u_max value to print.

    output_root : str
        output_root value to print.

    max_workers : int
        max_workers value to print.

    print_result : bool
        If True, print the generated block.

    Returns
    -------
    block : str
        The formatted production-call block.
    """
    # ------------------------------------------------------------
    # Resolve summary object
    # ------------------------------------------------------------
    if "global_recommended_settings" in summary_or_suite:
        summary = summary_or_suite
    elif "summary" in summary_or_suite and "global_recommended_settings" in summary_or_suite["summary"]:
        summary = summary_or_suite["summary"]
    else:
        raise ValueError(
            "Input must be either a summary dict from "
            "summarize_electron_thermal_convergence(...) "
            "or a suite dict containing ['summary']."
        )

    global_rec = summary.get("global_recommended_settings", {})

    def get_val(name, default):
        if name in global_rec:
            return global_rec[name].get("recommended_value", default)
        return default

    n_u = get_val("n_u", 40)
    n_t_speed = get_val("n_t_speed", 6)
    n_t_mu = get_val("n_t_mu", 6)
    n_scatter_mu = get_val("n_scatter_mu", 6)
    n_scatter_phi = get_val("n_scatter_phi", 8)
    shell_step = get_val("shell_step", 1)

    block = f'''{variable_name} = {function_name}(
    earth_data={earth_data_name},
    DM_masses={masses_name},
    sigma_SI_p={sigma_SI_p:.0e},
    sigma_SD_p={sigma_SD_p:.0e},
    sigma_electron={sigma_electron:.0e},
    output_root="{output_root}",
    u_max={u_max},
    n_u={n_u},
    rho_chi={rho_chi_name},
    v0={v0_name},
    n_t_speed={n_t_speed},
    n_t_mu={n_t_mu},
    n_scatter_mu={n_scatter_mu},
    n_scatter_phi={n_scatter_phi},
    shell_step={shell_step},
    u_grid_mode="log",
    max_workers={max_workers}
)

print("\\nElectron thermal production run finished.")
'''

    if print_result:
        print("\n" + "=" * 100)
        print("Recommended production block")
        print("=" * 100)
        print(block)
        print("=" * 100)

    return block

# ============================================================
# Main
if __name__ == "__main__":
    import os
    import numpy as np

    print("Loading Earth data...")
    earth_data = load_earth_composition(
        filepath="data/earth_prem.dat",
        sd_mode="verified_only",
        min_mass_fraction=1e-10
    )

    print(f"Loaded {len(earth_data['radius'])} layers, {len(earth_data['abundances'])} elements")
    print("[INFO] current working directory:", os.getcwd())

    # ------------------------------------------------------------
    # Sanity check: verified-only SD active labels
    # ------------------------------------------------------------
    print_sd_db_status_summary(sd_mode="verified_only", show_all=False)

    # ------------------------------------------------------------
    # Mass grid for SD baseline / thermal runs
    # ------------------------------------------------------------
    DM_masses_baseline = np.logspace(np.log10(0.3), np.log10(300.0), 161)

    # ------------------------------------------------------------
    # 1. Verified-only refined SD baseline
    #    SD only, constant, T=0
    # ------------------------------------------------------------
    sd_baseline_verified = plot_verified_only_refined_sd_baseline(
        earth_data=earth_data,
        DM_masses=DM_masses_baseline,
        sigma_SD_p=1e-40,
        output_root=".",
        u_max=800.0,
        n_u=40,
        rho_chi=RHO_CHI_DEFAULT,
        v0=V0_DEFAULT,
        n_t_speed=4,
        n_t_mu=4,
        n_scatter_mu=5,
        n_scatter_phi=6,
        max_workers=None
    )

    # ------------------------------------------------------------
    # 2. Refined verified-only SD thermal comparison
    #    SD only, constant, T(r)!=0 vs T=0
    #    Refined settings from convergence study:
    #      n_u = 80
    #      n_t_speed = 4
    #      n_t_mu = 4
    #      n_scatter_mu = 8
    #      n_scatter_phi = 12
    #      shell_step = 1
    # ------------------------------------------------------------
    sd_thermal_verified_refined = plot_verified_only_sd_thermal_constant_refined(
        earth_data=earth_data,
        DM_masses=DM_masses_baseline,
        sigma_SD_p=1e-40,
        output_root=".",
        u_max=800.0,
        n_u=80,
        rho_chi=RHO_CHI_DEFAULT,
        v0=V0_DEFAULT,
        n_t_speed=4,
        n_t_mu=4,
        n_scatter_mu=8,
        n_scatter_phi=12,
        max_workers=10
    )

    print("\nAll verified-only SD baseline and thermal-refined runs finished.")
        # ------------------------------------------------------------
    # 3. Verified-only SD operator comparison at T = 0
    #    Compare:
    #      - constant
    #      - v2_dependent
    #      - q2_dependent
    # ------------------------------------------------------------
    sd_operator_verified = plot_verified_only_sd_operator_comparison(
        earth_data=earth_data,
        DM_masses=DM_masses_baseline,
        sigma_SD_p=1e-40,
        output_root=".",
        u_max=800.0,
        n_u=80,
        rho_chi=RHO_CHI_DEFAULT,
        v0=V0_DEFAULT,
        n_t_speed=4,
        n_t_mu=4,
        n_scatter_mu=8,
        n_scatter_phi=12,
        shell_step=1,
        u_grid_mode="log",
        max_workers=10
    )

    print("\nAll verified-only SD baseline, thermal-refined, and operator-comparison runs finished.")

    # ------------------------------------------------------------
    # 4. Verified-only SD per-target operator decomposition at T = 0
    #    Operators:
    #      - constant
    #      - v2_dependent
    #      - q2_dependent
    #    Target decomposition:
    #      - H
    #      - Al
    #      - Si29
    #      - Na
    # ------------------------------------------------------------
    sd_operator_target_decomp = plot_verified_only_sd_operator_target_decomposition(
        earth_data=earth_data,
        DM_masses=DM_masses_baseline,
        sigma_SD_p=1e-40,
        output_root=".",
        u_max=800.0,
        n_u=80,
        rho_chi=RHO_CHI_DEFAULT,
        v0=V0_DEFAULT,
        n_t_speed=4,
        n_t_mu=4,
        n_scatter_mu=8,
        n_scatter_phi=12,
        shell_step=1,
        u_grid_mode="log",
        max_workers=10
    )

    print("\nAll verified-only SD baseline, thermal, operator comparison, and target decomposition runs finished.")
    # ------------------------------------------------------------
    # 5. Combined SI / verified-only SD / electron plot
    #    constant, T = 0
    # ------------------------------------------------------------
    print("\n[WARN] plot_combined_SI_verified_only_SD_electron_constant_T0 not found; skipping section 5.")

    # ------------------------------------------------------------
    # 6. Combined SI / verified-only SD / electron operator grid
    #    rows = constant / v2 / q2
    #    left = absolute capture rates
    #    right = informative comparison panels
    # ------------------------------------------------------------
    combined_operator_grid = plot_combined_SI_verified_only_SD_electron_operator_grid_polished(
        earth_data=earth_data,
        DM_masses=DM_masses_baseline,
        sigma_SI_p=1e-40,
        sigma_SD_p=1e-40,
        sigma_electron=1e-40,
        output_root=".",
        u_max=800.0,
        n_u=80,
        rho_chi=RHO_CHI_DEFAULT,
        v0=V0_DEFAULT,
        n_t_speed=4,
        n_t_mu=4,
        n_scatter_mu=8,
        n_scatter_phi=12,
        shell_step=1,
        u_grid_mode="log",
        max_workers=10
    )

    print("\nCombined SI / verified-only SD / electron operator grid run finished.")
    DM_masses_thermal_test = np.array([0.3, 1.0, 3.0, 10.0, 30.0, 100.0])

    combined_thermal_grid = plot_combined_SI_verified_only_SD_electron_thermal_grid(
        earth_data=earth_data,
        DM_masses=DM_masses_thermal_test,
        sigma_SI_p=1e-40,
        sigma_SD_p=1e-40,
        sigma_electron=1e-40,
        output_root=".",
        u_max=800.0,
        n_u=30,
        rho_chi=RHO_CHI_DEFAULT,
        v0=V0_DEFAULT,
        n_t_speed=4,
        n_t_mu=4,
        n_scatter_mu=4,
        n_scatter_phi=6,
        shell_step=2,
        u_grid_mode="log",
        max_workers=4
    )

    print("\nCombined SI / verified-only SD / electron thermal grid test run finished.")
    # ------------------------------------------------------------
    # 7. Combined SI / verified-only SD / electron thermal grid
    #    quick sanity run with reduced settings
    # ------------------------------------------------------------
    DM_masses_thermal_test = np.array([0.3, 1.0, 3.0, 10.0, 30.0, 100.0])

    combined_thermal_grid_test = plot_combined_SI_verified_only_SD_electron_thermal_grid(
        earth_data=earth_data,
        DM_masses=DM_masses_thermal_test,
        sigma_SI_p=1e-40,
        sigma_SD_p=1e-40,
        sigma_electron=1e-40,
        output_root=".",
        u_max=800.0,
        n_u=30,
        rho_chi=RHO_CHI_DEFAULT,
        v0=V0_DEFAULT,
        n_t_speed=4,
        n_t_mu=4,
        n_scatter_mu=4,
        n_scatter_phi=6,
        shell_step=2,
        u_grid_mode="log",
        max_workers=4
    )

    print("\nCombined SI / verified-only SD / electron thermal grid test run finished.")

    # ------------------------------------------------------------
    # 8. Electron thermal convergence suite
    #    This checks whether the large electron thermal enhancement
    #    is numerically stable.
    # ------------------------------------------------------------
    electron_conv_suite = run_and_plot_electron_thermal_convergence_suite(
        earth_data=earth_data,
        masses=(0.3, 1.0, 10.0, 100.0),
        sigma_electron=1e-40,
        cross_section_types=("constant", "v2_dependent", "q2_dependent"),
        rho_chi=RHO_CHI_DEFAULT,
        v0=V0_DEFAULT,
        u_max=800.0,
        ratio_tol=0.05,
        rate_tol=0.05,
        output_root=".",
        summary_txt="logs/electron_thermal_convergence_summary.txt",
        show_plots=False,
        top_note=True
    )

    print("\nElectron thermal convergence suite finished.")

    # ------------------------------------------------------------
    # 9. Print recommended production block from convergence summary
    # ------------------------------------------------------------
    production_block = print_electron_thermal_recommended_production_block(
        electron_conv_suite,
        function_name="plot_combined_SI_verified_only_SD_electron_thermal_grid",
        variable_name="combined_thermal_grid_production",
        earth_data_name="earth_data",
        masses_name="DM_masses_baseline",
        sigma_SI_p=1e-40,
        sigma_SD_p=1e-40,
        sigma_electron=1e-40,
        rho_chi_name="RHO_CHI_DEFAULT",
        v0_name="V0_DEFAULT",
        u_max=800.0,
        output_root=".",
        max_workers=10,
        print_result=True
    )
    print("\nRecommended production block has been printed.")

    # ------------------------------------------------------------
    # 10. Optional production run
    #     Uncomment only after you are satisfied with the
    #     electron convergence summary.
    # ------------------------------------------------------------
    # combined_thermal_grid_production = plot_combined_SI_verified_only_SD_electron_thermal_grid(
    #     earth_data=earth_data,
    #     DM_masses=DM_masses_baseline,
    #     sigma_SI_p=1e-40,
    #     sigma_SD_p=1e-40,
    #     sigma_electron=1e-40,
    #     output_root=".",
    #     u_max=800.0,
    #     n_u=40,
    #     rho_chi=RHO_CHI_DEFAULT,
    #     v0=V0_DEFAULT,
    #     n_t_speed=6,
    #     n_t_mu=6,
    #     n_scatter_mu=6,
    #     n_scatter_phi=8,
    #     shell_step=1,
    #     u_grid_mode="log",
    #     max_workers=10
    # )
    #
    # print("\nCombined SI / verified-only SD / electron thermal grid production run finished.")
        # ------------------------------------------------------------
    # 11. Combined SI / verified-only SD / electron thermal grid
    #     production run
    # ------------------------------------------------------------
    combined_thermal_grid_production = plot_combined_SI_verified_only_SD_electron_thermal_grid(
        earth_data=earth_data,
        DM_masses=DM_masses_baseline,
        sigma_SI_p=1e-40,
        sigma_SD_p=1e-40,
        sigma_electron=1e-40,
        output_root=".",
        u_max=800.0,
        n_u=60,
        rho_chi=RHO_CHI_DEFAULT,
        v0=V0_DEFAULT,
        n_t_speed=4,
        n_t_mu=4,
        n_scatter_mu=6,
        n_scatter_phi=4,
        shell_step=1,
        u_grid_mode="log",
        max_workers=10
    )

    print("\nCombined SI / verified-only SD / electron thermal grid production run finished.")
