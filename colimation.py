#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import time
from math import sin, cos, sqrt, pi, exp, atan2, hypot
from scipy.integrate import quad, nquad
from scipy.special import erf, i0
import scipy.stats as stats
from constants import *

tex0 = 0.0001779/10
x0 = 0.0001686/10
tey0 = 2.23998e-06
y0 = 1.339297e-05/10
E0 = 2000.e6
w0 = h * c / 527e-9 / e

def omega_max(E_0, w_0, alpha=pi):
    lam = (4 * E_0 * w_0 * sin(0.5*alpha)**2) / me**2
    return ((E_0 * lam) / (1 + lam))
def theta_c(E_0, w_0):
    return (me / E_0 * sqrt(1 + 4*E_0*w_0 / me**2))
def omega_theta(E_0, w_0, theta, alpha=pi):
    return omega_max(E_0, w_0, alpha) / (1 + (theta / theta_c(E_0, w_0))**2)
def theta_omega(E_0, w_0, w, alpha=pi):
    return theta_c(E_0, w_0) * sqrt(omega_max(E_0, w_0, alpha) / w - 1)
def delta_ring(x, y, r):
    if x**2 + y**2 == r**2:
        return 1
    else:
        return 0

L = 50.
hole = 0.0005
N = int(1e6)
wmax = omega_max(E0, w0)

def ellipsoid_hole(x, y, a, b):
    if x**2 / a**2 + y**2 / b**2 == 1:
        return 1
    else:
        return 0
def gamma_func(alpha_func, beta_func):
    return (1 + alpha_func**2)/ beta_func

#pdf -- particle distribution function
def pdf(r, phi, phi0, r0, sx=0.5, sy=1):
    x = (r*cos(phi) - r0*cos(phi0))/sx
    y = (r*sin(phi) - r0*sin(phi0))/sy
    return r*exp(-0.5*x**2 - 0.5*y**2) / 2/pi/sx/sy

sx = np.hypot(tex0, x0/L)
sy = np.hypot(tey0, y0/L)

def beam_generator(alpha_x, alpha_y, beta_x, beta_y, eps_x, eps_y):
    """Function that generates beam distribution using np.random and other numpy tools. Should be equals to the pdf"""
    mu_x, sigma_x = 0, 1         # Mean and standard deviation
    lower_bound_x = -5 * sigma_x # -5 sigma
    upper_bound_x = 5 * sigma_x  # +5 sigma

    # Truncnorm requires standard bounds scaled by sigma
    a_x = (lower_bound_x - mu_x) / sigma_x
    b_x = (upper_bound_x - mu_x) / sigma_x

    mu_y, sigma_y = 0, 1.2         # Mean and standard deviation
    lower_bound_y = -5 * sigma_y # -5 sigma
    upper_bound_y = 5 * sigma_y  # +5 sigma

    # Truncnorm requires standard bounds scaled by sigma
    a_y = (lower_bound_y - mu_y) / sigma_y
    b_y = (upper_bound_y - mu_y) / sigma_y


    gf_x = gamma_func(alpha_x, beta_x)
    gf_y = gamma_func(alpha_y, beta_y)
    coord_disp_x = sqrt(eps_x * beta_x)
    coord_disp_y = sqrt(eps_y * beta_y)

    rng = np.random.default_rng()

    r_x = np.arange(-coord_disp_x, coord_disp_x, 0.01)
    gauss_x = stats.truncnorm.rvs(a_x, b_x, loc=mu_x, scale=sigma_x, size=10000)
    angle_distr_x = rng.uniform(0, 2*pi, 10000)
    angle_distr_x *= gauss_x
    #angle_distr_x *= rng.uniform(0, 7.2)
    gauss_y = stats.truncnorm.rvs(a_y, b_y, loc=mu_y, scale=sigma_y, size=3775)
    r_y = np.arange(-coord_disp_y, coord_disp_y, 0.01)
    angle_distr_y = rng.uniform(0, 2*pi, 3775)
    angle_distr_y *= gauss_y
    #angle_distr_y = rng.normal(0, 3.3)

    x_ref = sqrt(eps_x * beta_x) * r_x * np.cos(angle_distr_x)
    x_hatch_ref = sqrt(eps_x / beta_x) * alpha_x * r_x * np.cos(angle_distr_x) + sqrt(eps_x / beta_x) * r_x * np.sin(angle_distr_x)

    y_ref = sqrt(eps_y * beta_y) * r_y * np.cos(angle_distr_y)
    y_hatch_ref = sqrt(eps_y / beta_y) * alpha_y * r_y * np.cos(angle_distr_y) + sqrt(eps_y / beta_y) * r_y * np.sin(angle_distr_y)

    np.savetxt("2d_distr_coord_x.txt", x_ref)
    np.savetxt("2d_distr_angle_x.txt", x_hatch_ref)

    np.savetxt("2d_distr_coord_y.txt", y_ref)
    np.savetxt("2d_distr_angle_y.txt", y_hatch_ref)

if __name__ == "__main__":
    beta_x = 100 #cm
    beta_y = 15 #cm
    eps_x = 25 #nm*rad
    eps_y = 0.95 * eps_x #nm*rad
    alpha_x = 5 #cm
    alpha_y = 0.3 #cm
    beam_generator(alpha_x, alpha_y, beta_x, beta_y, eps_x, eps_y)

