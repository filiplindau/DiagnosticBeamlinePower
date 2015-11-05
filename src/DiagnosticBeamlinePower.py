# -*- coding:utf-8 -*-
"""
Created on Nov 5, 2015

@author: Filip Lindau
"""

from numpy import *
from scipy.special import kv
from scipy.integrate import quad

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

# From Juan I. Larruquert, J. Opt. Soc. Am. A / Vol. 28, No. 11 / November 2011 / 2340
# Self-consistent optical constants of SiC thin films
SiC_l = array([200, 300, 400, 500, 600, 700, 800, 900, 1000 , 1100, 1200, 1300, 1400])
SiC_n = array([2.2, 3.35, 3.45, 3.4, 3.4, 3.35, 3.3, 3.25, 3.25, 3.25, 3.25, 3.25, 3.25])
SiC_k = array([2.25, 1.25, 0.75, 0.5, 0.3, 0.25, 0.2, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15])

# From A. M. Hofmeister: The Astrophysical Journal, 696:1502–1516, 2009 May 10
# OPTICAL CONSTANTS OF SILICON CARBIDE FOR ASTROPHYSICAL APPLICATIONS. II. 
# EXTENDING OPTICAL FUNCTIONS FROM INFRARED TO ULTRAVIOLET USING SINGLE-CRYSTAL ABSORPTION SPECTRA
SiC_R = ((SiC_n - 1) ** 2 + SiC_k ** 2) / ((SiC_n + 1) ** 2 + SiC_k ** 2)

K = lambda x: kv(5 / 3.0, x)
dWdw = []
for omegax in omega:
    ostart = omegax / omegac
    S = quad(K, ostart, inf)
    dWdw.append(sqrt(3) * qe ** 2 / (4 * pi * eps0 * c) * gamma * omegax / omegac * S[0] * I / qe)
dWdw = array(dWdw)
