# -*- coding:utf-8 -*-
"""
Created on Nov 5, 2015

@author: Filip Lindau
"""

from numpy import *
from scipy.special import kv
from scipy.integrate import quad
from scipy.interpolate import interp1d

c = 299792458.0
qe = 1.602e-19
me = 9.11e-31
h = 6.626e-34
eps0 = 8.85418782e-12

hb = h / (2 * pi)

B = 0.5
E = 3e9 * qe
I = 500e-3

gamma = E / (me * c ** 2)

# From http://photon-science.desy.de/sites/site_photonscience/content/e62/e189219/e187240/e187241/e187242/infoboxContent187244/f2_eng.pdf
Ec = 1.5 * qe * hb * c ** 2 * B * E ** 2 / (me * c ** 2) ** 3
# Ec = 665.0255 * B * (E * 1e-9 / qe) ** 2 * qe
lc = h * c / Ec
omegac = 2 * pi * c / lc

l = linspace(200e-9, 1400e-9, 100)
# l = logspace(-10, -6, 100)
omega = flipud(2 * pi * c / l)
l = flipud(l)

omega = linspace(2 * pi * c / 1.4e-6, 2 * pi * c / 0.2e-6, 100)
l = 2 * pi * c / omega

# From Juan I. Larruquert, J. Opt. Soc. Am. A / Vol. 28, No. 11 / November 2011 / 2340
# Self-consistent optical constants of SiC thin films
SiC_l = array([199, 300, 400, 500, 600, 700, 800, 900, 1000 , 1100, 1200, 1300, 1400]) * 1e-9
SiC_n = array([2.2, 3.35, 3.45, 3.4, 3.4, 3.35, 3.3, 3.25, 3.25, 3.25, 3.25, 3.25, 3.25])
SiC_k = array([2.25, 1.25, 0.75, 0.5, 0.3, 0.25, 0.2, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15])

# From A. M. Hofmeister: The Astrophysical Journal, 696:1502–1516, 2009 May 10
# OPTICAL CONSTANTS OF SILICON CARBIDE FOR ASTROPHYSICAL APPLICATIONS. II. 
# EXTENDING OPTICAL FUNCTIONS FROM INFRARED TO ULTRAVIOLET USING SINGLE-CRYSTAL ABSORPTION SPECTRA
SiC_R = ((SiC_n - 1) ** 2 + SiC_k ** 2) / ((SiC_n + 1) ** 2 + SiC_k ** 2)

SiC_int = interp1d(SiC_l, SiC_R, fill_value=0.28, bounds_error=False)


Al_UVenh = array([219.8, 76.546,
227.3, 79.598,
236, 83.349,
247.1, 87.076,
262, 89.907,
280.6, 91.7,
301.7, 92.594,
330.2, 92.797,
361.2, 92.378,
398.4, 91.542,
438.1, 90.016,
477.8, 88.586,
518.8, 86.899,
554.7, 85.515,
598.1, 84.062,
629.1, 83.096,
666.4, 82.265,
698.6, 81.434,
725.9, 80.673,
749.5, 79.911,
768.1, 79.218,
786.7, 78.399,
796.6, 77.906,
807.8, 77.56,
821.4, 77.355,
835, 77.635,
849.9, 78.333,
869.8, 80.074,
894.6, 82.707,
916.9, 84.989,
940.5, 87.065,
965.3, 88.794,
996.3, 90.293])

Al_UVenh_l = Al_UVenh[0::2] * 1e-9
Al_UVenh_R = Al_UVenh[1::2] * 1e-2
Al_UVenh_int = interp1d(Al_UVenh_l, Al_UVenh_R, fill_value=0.95, bounds_error=False)


# Retinal response from ICNIRP GUIDELINES
# ON LIMITS OF EXPOSURE TO INCOHERENT VISIBLE
# AND INFRARED RADIATION
# PUBLISHED IN: HEALTH PHYSICS 105(1):74‐96; 2013
retina = array([200, 6.0, 0.01, 0,
300, 6.00, 0.01, 0,
305, 6.00, 0.01, 0,
310, 6.00, 0.01, 0,
315, 6.00, 0.01, 0,
320, 6.00, 0.01, 0,
330, 6.00, 0.01, 0,
335, 6.00, 0.01, 0,
340, 5.88, 0.01, 0,
345, 5.71, 0.01, 0,
350, 5.46, 0.01, 0,
355, 5.22, 0.01, 0,
360, 4.62, 0.01, 0,
365, 4.29, 0.01, 0,
370, 3.75, 0.01, 0,
375, 3.56, 0.01, 0,
380, 3.19, 0.01, 0.01,
385, 2.31, 0.0125, 0.0125,
390, 1.88, 0.025, 0.025,
395, 1.58, 0.050, 0.05,
400, 1.43, 0.100, 0.1,
405, 1.30, 0.200, 0.2,
410, 1.25, 0.400, 0.4,
415, 1.20, 0.800, 0.8,
420, 1.15, 0.900, 0.9,
425, 1.11, 0.950, 0.95,
430, 1.07, 0.980, 0.98,
435, 1.03, 1.000, 1.0,
440, 1.000, 1.000, 1.0,
445, 0.970, 0.970, 1.0,
450, 0.940, 0.940, 1.0,
455, 0.900, 0.900, 1.0,
460, 0.800, 0.800, 1.0,
465, 0.700, 0.700, 1.0,
470, 0.620, 0.620, 1.0,
475, 0.550, 0.550, 1.0,
480, 0.450, 0.450, 1.0,
485, 0.400, 0.400, 1.0,
490, 0.220, 0.220, 1.0,
495, 0.160, 0.160, 1.0,
500, 0.100, 0.100, 1.0,
505, 0.079, 0.079, 1.0,
510, 0.063, 0.063, 1.0,
515, 0.050, 0.050, 1.0,
520, 0.040, 0.040, 1.0,
525, 0.032, 0.032, 1.0,
530, 0.025, 0.025, 1.0,
535, 0.020, 0.020, 1.0,
540, 0.016, 0.016, 1.0,
545, 0.013, 0.013, 1.0,
550, 0.010, 0.010, 1.0,
555, 0.008, 0.008, 1.0,
560, 0.006, 0.006, 1.0,
565, 0.005, 0.005, 1.0,
570, 0.004, 0.004, 1.0,
575, 0.003, 0.003, 1.0,
580, 0.002, 0.002, 1.0,
585, 0.002, 0.002, 1.0,
590, 0.001, 0.001, 1.0,
595, 0.001, 0.001, 1.0,
600, 0.001, 0.001, 1.0,
700, 0.001, 0.001, 1.0,
1050, 0, 0, 1.0,
1200, 0, 0, 0.2,
1400, 0, 0, 0.02])
l_ret = arange(300, 1400, 5) * 1e-9  # Wavelengths for retinal response data
l_ret = retina[0::4] * 1e-9  # Wavelengths for retinal response data
B_ret = retina[2::4]  # Blue light hazard function 
th_ret = retina[3::4]  # Thermal hazard function
B_int = interp1d(l_ret, B_ret, fill_value=0.0, bounds_error=False)
th_int = interp1d(l_ret, th_ret, fill_value=0.0, bounds_error=False)

###########################################
# Bending magnet calculations
#-------------------------------

# From wikipedia, integrated over all angles:
K = lambda x: kv(5 / 3.0, x)
dWdw = []
G1 = []
for omegax in omega:
    ostart = omegax / omegac
    S = quad(K, ostart, inf)
    dWdw.append(sqrt(3) * qe ** 2 / (4 * pi * eps0 * c) * gamma * omegax / omegac * S[0] * I / qe)
    G1.append(ostart * S[0])
dWdw = array(dWdw)
G1 = array(G1)

# Total visible power
P = trapz(dWdw, omega)

# Calculation from http://www.eecs.berkeley.edu/~attwood/srms/2007/Lec08.pdf
d2Fdthdo = 2.46e13 * E / qe / 1e9 * I * G1
dth = 1 / gamma
Eph = h * omega / (2 * pi)
dPdw = d2Fdthdo * dth * 1e3 * Eph
dww = gradient(omega) / omega / 0.1e-2
# Visible power:
P2 = sum(dPdw * dww)

# Acceptance angle of the beamline:
d = 5e-2  # Aperture diameter
R = 4  # Distance to aperture


# From wikipedia, radiation integral
# Looking at theta=0 (in the plane of the electrons)
rho = gamma * me * c / (B * qe)
K23 = lambda x: kv(2 / 3.0, x)

theta = 0
xi = rho * omega / (3 * c * gamma ** 3) * (1 + gamma ** 2 * theta ** 2) ** 1.5

# Radiation power (scaled with ring current) per steradian and frequency interval:
dWdodw = qe ** 2 / (16 * pi ** 3 * eps0 * c) * (2 * omega * rho / (3 * c * gamma ** 2)) ** 2 * (1 + gamma ** 2 * theta ** 2) ** 2 * (kv(2.0 / 3, xi) ** 2 + gamma ** 2 * theta ** 2 / (1 + gamma ** 2 * theta ** 2) * kv(1.0 / 3, xi) ** 2) * I / qe
# The solid angle at exposure over 1 cm**2 at a distance R is 1e-4/R**2
dWdw_ret = dWdodw * 1e-4 / R ** 2

# Assume 1/gamma angular beam in y dir (not correct)
dOmega = 1 / gamma * d / R

# Assume beam filling the mirror in both directions:
dOmega = (d / 2) ** 2 * pi / R ** 2  # Solid angle of the beam
A = dOmega * R ** 2

P3 = trapz(dWdodw * SiC_int(l) * Al_UVenh_int(l) * dOmega, omega)
# P3 = trapz(dWdodw * dOmega, omega)

I = P3 / A * 1e-4  # Total W/cm**2

uvInd = l < 400e-9
Puv = trapz(dWdodw[uvInd] * SiC_int(l[uvInd]) * Al_UVenh_int(l[uvInd]) * dOmega, omega[uvInd])
Iuv = Puv / A * 1e-4  # UV W/cm**2
MPEuv = 3e-3

visInd = all([l < 400e-9, l > 700e-9], axis=0)
Pvis = trapz(dWdodw[visInd] * SiC_int(l[visInd]) * Al_UVenh_int(l[visInd]) * dOmega, omega[visInd])
Ivis = Pvis / A * 1e-4  # Visible W/cm**2
MPEvis = 0.75e-3

# Then there is the lens...
