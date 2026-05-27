#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import time
from math import sin, cos, sqrt, pi, exp, atan2, hypot
from scipy.integrate import quad, nquad
from scipy.special import erf, i0
from constants import *

tex0 = 0.0001779/10
x0 = 0.0001686/10
tey0 = 2.23998e-06
y0 = 1.339297e-05/10
E0 = 2000.e6
w0 = h * c / 527e-9 / e
L = 50.
hole = 0.0005
N = int(1e6)
wmax = w_max(E0, w0)

def omega_max(E_0, w_0, alpha=pi):
    lam = (4 * E_0 * w_0 * sin(0.5*alpha)**2) / me**2
    return ((E_0 * lam) / (1 + lam))
def theta_c(E_0, w_0):
    return (me / E_0 * sqrt(1 + 4*E_0*w_0 / me**2))
def omega_theta(E_0, w_0, alpha=pi, theta):
    return w_max(E_0, w_0, alpha) / (1 + (theta / theta_c(E_0, w_0))**2)
def theta_omega(E_0, w_0, alpha=pi, w):
    return theta_c(E_0, w_0) * sqrt(omega_max(E_0, w_0, alpha) / w - 1)
def delta_ring(x, y, r):
    if x**2 + y**2 == r**2:
        return 1
    else:
        return 0
def ellipsoid_hole(x, y, a, b):
    if x**2 / a**2 + y**2 / b**2 == 1:
        return 1
    else:
        return 0
def gamma_func_KEDR(beta_func_KEDR):
    return 1 / beta_func_KEDR
#pdf -- particle distribution function
def pdf(r, phi, phi0, r0, )
