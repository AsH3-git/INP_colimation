#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import time
from math import sin, cos, sqrt, pi, exp, atan2, hypot
from scipy.integrate import quad, nquad
from scipy.special import erf, i0
import scipy.stats as stats
from constants import *

def lam(E_0, w_0, alpha=pi):
    return (4 * E_0 * w_0 * sin(0.5*alpha)**2) / me**2
def omega_max(E_0, w_0, alpha=pi):
    return ((E_0 * lam(E_0, w_0)) / (1 + lam(E_0, w_0)))
def theta_c(E_0, w_0):
    return (me / E_0 * np.sqrt(1 + 4*E_0*w_0 / me**2))
def omega_theta(E_0, w_0, theta, alpha=pi):
    return omega_max(E_0, w_0, alpha) / (1 + (theta / theta_c(E_0, w_0))**2)
def theta_omega(E_0, w_0, w, alpha=pi):
    return theta_c(E_0, w_0) * np.sqrt((omega_max(E_0, w_0) / w) - 1)
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
def gamma_func(alpha_func, beta_func):
    return (1 + alpha_func**2)/ beta_func

#pdf -- particle distribution function
def pdf(r, phi, phi0, r0, sx=0.5, sy=1):
    x = (r*cos(phi) - r0*cos(phi0))/sx
    y = (r*sin(phi) - r0*sin(phi0))/sy
    return r*exp(-0.5*x**2 - 0.5*y**2) / 2/pi/sx/sy

#sx = np.hypot(tex0, x0/L)
#sy = np.hypot(tey0, y0/L)

def beam_generator(alpha_x, alpha_y, beta_x, beta_y, eps_x, eps_y, N_gen):
    """Function that generates beam distribution using np.random and other numpy tools. Should be equals to the pdf"""
    gf_x = gamma_func(alpha_x, beta_x)
    #axis = np.linspace(-sqrt(eps_x * beta_x), sqrt(eps_x * beta_x), num=10000)

    mu_x, sigma_x = 0, 1    
    mu_y, sigma_y = 0, 1.2  

    gauss_x = np.random.normal(loc=mu_x, scale=sigma_x, size=N_gen)
    gauss_x_hatch = np.random.normal(loc=mu_x, scale=sigma_x, size=N_gen)

    mask = (gauss_x**2 + gauss_x_hatch**2) <= 25.0
    coord_distr_x = gauss_x[mask]
    angle_distr_x = gauss_x_hatch[mask]

    x = np.sqrt(eps_x * beta_x) * coord_distr_x
    x_prime = np.sqrt(eps_x / beta_x) * (-alpha_x * coord_distr_x + angle_distr_x)
    x_plane = np.vstack((x, x_prime))

    gauss_y = np.random.normal(loc=mu_y, scale=sigma_y, size=N_gen)
    gauss_y_hatch = np.random.normal(loc=mu_y, scale=sigma_y, size=N_gen)

    mask = (gauss_y**2 + gauss_y_hatch**2) <= 25.0
    coord_distr_y = gauss_y[mask]
    angle_distr_y = gauss_y_hatch[mask]

    y = np.sqrt(eps_y * beta_y) * coord_distr_y
    y_prime = np.sqrt(eps_y / beta_y) * (-alpha_y * coord_distr_y + angle_distr_y)
    if len(x) > len(y):
        for i in range(len(y), len(x)):
            y = np.append(y, 0)
            y_prime = np.append(y_prime, 0)
    elif len(x) < len(y):
        for i in range(len(x), len(y)):
            x = np.append(x, 0)
            x_prime = np.append(x_prime, 0)
        x_plane = np.vstack((x, x_prime))
    else:
        print("Arrays along planes (x, x`) and (y, y`) have same dimensions!")
    y_plane = np.vstack((y, y_prime))
    #print(np.shape(x_plane))
    #print(np.shape(y_plane))

    delta_E = np.random.normal(loc=0.0, scale=5e-3, size=len(y))

    four_dim_space = np.vstack((x_plane, y_plane))
    four_dim_space_energy = np.vstack((four_dim_space, delta_E))
    #print(np.shape(four_dim_space))
    four_dim_space_arr = four_dim_space.T

    np.savetxt("four_dim_space_distr.txt", four_dim_space_arr)

    #x0_prime = -axis * alpha_x/beta_x

    plt.figure(figsize=(12,8))
    plt.scatter(four_dim_space[0], four_dim_space[1], s=3, alpha=0.4)
    plt.scatter(four_dim_space[2], four_dim_space[3], s=3, alpha=0.4)
    #plt.plot(axis, x0_prime, color="black")
    #plt.show()

    return four_dim_space_energy, len(x)

def compton_backscattering(E_0, w_0, L, r_0):
    four_dim_space_with_energy, N_compton = beam_generator(alpha_x, alpha_y, beta_x, beta_y, eps_x, eps_y, N_gen=int(1e6))
    CBS_max_energy = omega_max(E_0, w_0)
    energy_distribution = (four_dim_space_with_energy[4] + 1) * E_0
    CBS_energy_distribution = omega_max(energy_distribution, w_0)
    w = np.random.uniform(0, CBS_energy_distribution, N_compton)
    CBS_angle = theta_omega(energy_distribution, w_0, w)
    w0_gauss = np.random.normal(loc=2.35, scale=1e-4, size=N_compton)
    CBS_energy_distribution_with_gauss_omega0 = omega_max(energy_distribution, w0_gauss)
    y = CBS_energy_distribution_with_gauss_omega0 / energy_distribution
    #epsilon = (CBS_energy_distribution_with_gauss_omega0) / (me * c**2)
    
    a_0 = 1 / (1 - y)
    a_1 = (4 * y) / (lam(energy_distribution, w0_gauss)*(1 - y))
    a_2 = (4 * y**2) / (lam(energy_distribution, w0_gauss)**2 * (1 - y)**2)
    Klein_Nishina_disp = (2*pi*r_0**2) / (energy_distribution * lam(energy_distribution, w0_gauss)) * (a_0 + 1 - y - a_1 + a_2)
    np.savetxt("spectra.txt", Klein_Nishina_disp)

if __name__ == "__main__":
    tex0 = 0.0001779/10
    x0 = 0.0001686/10
    tey0 = 2.23998e-06
    y0 = 1.339297e-05/10
    E_0 = 2.0e9 #eV
    w_0 = h * c / 527e-9 / e
    r_0 = 2.8179e-15

    L = 50. #m

    hole = 5.e-4 #m
    N_gen = int(1e6)

    beta_x = 1. #m
    beta_y = 0.15 #m
    eps_x = 0.25e-5 #m*rad
    eps_y = 0.95 * eps_x #m*rad
    alpha_x = 5.e-2 #m
    alpha_y = 3.e-3 #m

    compton_backscattering(E_0, w_0, L, r_0)
