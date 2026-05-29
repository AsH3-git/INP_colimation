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
def gamma_func(beta_func):
    return 1 / beta_func

#pdf -- particle distribution function
def pdf(r, phi, phi0, r0, sx=0.5, sy=1):
    x = (r*cos(phi) - r0*cos(phi0))/sx
    y = (r*sin(phi) - r0*sin(phi0))/sy
    return r*exp(-0.5*x**2 - 0.5*y**2) / 2/pi/sx/sy

sx = np.hypot(tex0, x0/L)
sy = np.hypot(tey0, y0/L)

def beam_generator(alpha_func, beta_func, epsilon):
    """Function that generates beam distribution using np.random and other numpy tools. Should be equals to the pdf"""
    gf = gamma_func(beta_func)
    coord_disp = sqrt(epsilon * beta_func)
    x = np.arange(-coord_disp, coord_disp, 0.5)
    x0_hatch = -x * alpha_func/beta_func
    x0 = -x * alpha_func/gf
    rng = np.random.default_rng()
    weights_arr = np.array([])
    total_distr_angle = np.array([])
    total_distr_coord = np.array([])
    angle_radius = sqrt(epsilon / beta_func)
    coord_radius = sqrt(epsilon / gf)
    for i in range(len(x)):
        j = abs(i - len(x)/2)
        weights_arr = np.append(weights_arr, j)
        low_angle = min(x0_hatch[i] - angle_radius, x0_hatch[i] + angle_radius)
        high_angle = max(x0_hatch[i] - angle_radius, x0_hatch[i] + angle_radius)
        angle_dot = rng.uniform(low=low_angle, high=high_angle, size=100)
        low_coord = min(x0[i] - coord_radius, x0[i] + coord_radius)
        high_coord = max(x0[i] - coord_radius, x0[i] + coord_radius)
        coord_dot = rng.uniform(low=low_coord, high=high_coord, size=100)
        total_distr_angle = np.concatenate((total_distr_angle, angle_dot))
        total_distr_coord = np.concatenate((total_distr_coord, coord_dot))

    j_repl_idx = np.where(weights_arr == 0)
    weights_arr[j_repl_idx] = 1
    for i, element in enumerate(weights_arr):
        if element == 1:
            continue
        weights_arr[i] = abs(element - (max(weights_arr) - 1)) / max(weights_arr)
    #print(weights_arr)
    #weights_arr1 = weights_arr[:len(weights_arr)//2]
    #weights_arr2 = weights_arr[len(weights_arr)//2:]
    #weights_arr = np.concatenate((weights_arr1, weights_arr2))
    
    """
    total_distr_angle = total_distr_angle * weights_arr
    print(total_distr_angle)
    total_distr_coord = total_distr_coord * weights_arr
    
    angle_diff = x0_hatch - angle_radius
    coord_diff = x0 - coord_radius
    angle_diff = np.repeat(angle_diff, 100)
    coord_diff = np.repeat(coord_diff, 100)
    weights_arr = np.repeat(weights_arr, 100)
    mask_angle = (element < (angle_diff * weights_arr))
    mask_coord = (element < (coord_diff * weights_arr))

    total_distr_angle = total_distr_angle[mask_angle]
    total_distr_coord = total_distr_coord[mask_coord]
    """
    
    np.savetxt("2d_distr_coord.txt", total_distr_coord)
    np.savetxt("2d_distr_angle.txt", total_distr_angle)

if __name__ == "__main__":
    beta_func = 100 #cm
    epsilon = 25 #nm*rad
    alpha_func = 5 #cm
    beam_generator(alpha_func, beta_func, epsilon)

