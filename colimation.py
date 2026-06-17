#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import time
from math import sin, cos, sqrt, pi, exp, atan2, hypot
from scipy.integrate import quad, nquad
from scipy.special import erf, i0
import scipy.stats as stats
from constants import *
from numba import njit

@njit
def lam(E_0, w_0, alpha=pi):
    return (4 * E_0 * w_0 * sin(0.5*alpha)**2) / me**2
@njit
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
       
    return four_dim_space_energy, len(x)

@njit
def Klein_Nishina(E_0, w_0, w):
    y = w / E_0
    a_0 = 1 / (1 - y)
    a_1 = (4 * y) / (lam(E_0, w_0)*(1 - y))
    a_2 = (4 * y**2) / (lam(E_0, w_0)**2 * (1 - y)**2)
    Klein_Nishina_distr = (2*pi*r_0**2) / (E_0 * lam(E_0, w_0)) * (a_0 + 1 - y - a_1 + a_2)
    return Klein_Nishina_distr

@njit
def spectra_generator(energy_distribution, w0_gauss, num_samples):
    accepted_samples_E = np.empty(num_samples)
    accepted_samples_W = np.empty(num_samples)

    n_elements = len(energy_distribution)

    f_max = 0.0
    for i in range(n_elements):
        E_0 = energy_distribution[i]
        w_0 = w0_gauss[i]
        w_m = omega_max(E_0, w_0)

        val_start = Klein_Nishina(E_0, w_0, 0.0)
        val_end = Klein_Nishina(E_0, w_0, w_m)

        if val_start > f_max: f_max = val_start
        if val_end > f_max: f_max = val_end

    #f_max *= 1.12

    count = 0
    while count < num_samples:
        rand_idx = np.random.randint(0, n_elements)
        E_rand = energy_distribution[rand_idx]
        w0_rand = w0_gauss[rand_idx]

        w_max_local = omega_max(E_rand, w0_rand)
        w_rand = np.random.uniform(0.0, w_max_local)
        y_rand = np.random.uniform(0.0, f_max)
        f_val = Klein_Nishina(E_rand, w0_rand, w_rand)

        if y_rand <= f_val:
            accepted_samples_E[count] = E_rand
            accepted_samples_W[count] = w_rand
            count += 1
    return accepted_samples_E, accepted_samples_W

def compton_backscattering(E_0, w_0, L, r_0):
    four_dim_space_with_energy, N_compton = beam_generator(alpha_x, alpha_y, beta_x, beta_y, eps_x, eps_y, N_gen=int(1e5))
    CBS_max_energy = omega_max(E_0, w_0)
    energy_distribution = (four_dim_space_with_energy[4] + 1) * E_0
    CBS_energy_distribution = omega_max(energy_distribution, w_0)
    w0_gauss = np.random.normal(loc=2.35, scale=1e-4, size=N_compton)
    CBS_energy_distribution_with_gauss_omega0 = omega_max(energy_distribution, w0_gauss)

    Klein_Nishina_energy_distr, Klein_Nishina_photon_energy_distr = spectra_generator(energy_distribution, w0_gauss, N_compton)
    np.savetxt("Klein_Nishina.txt", Klein_Nishina_photon_energy_distr)

    CBS_angle = theta_omega(Klein_Nishina_energy_distr, w0_gauss, Klein_Nishina_photon_energy_distr)
    CBS_angle_with_energy = np.vstack((CBS_angle, Klein_Nishina_photon_energy_distr))
    CBS_angle_with_energy_transpose = CBS_angle_with_energy.T
    np.savetxt("angle.txt", CBS_angle_with_energy_transpose)

    return CBS_angle_with_energy, w0_gauss, energy_distribution, four_dim_space_with_energy, N_compton
"""
@njit
def bin_size_calculator(energy_array):
    start_val = np.cumsum(energy_array)
    mods = start_val % 1e6
    condition = (mods >= 9.9e5) | (mods <= 1.9e5)

    return np.sum(condition)
"""

def colimation_cycle(CBS_angle_with_energy, N_compton, four_dim_space_with_energy):
    theta = CBS_angle_with_energy[0][:N_compton]
    phi = np.random.uniform(0, 2 * np.pi, N_compton)
    
    theta_x = np.arcsin(np.cos(phi) * np.sin(theta))
    theta_y = np.arcsin(np.sin(phi) * np.sin(theta))
    """
    theta_rec = np.sqrt(theta_x**2 + theta_y**2)
    theta_recovered = np.vstack((theta_rec, CBS_angle_with_energy[1]))
    theta_recovered = theta_recovered.T
    np.savetxt("theta_rec.txt", theta_recovered)
    """
    photon_x = four_dim_space_with_energy[0][:N_compton]
    photon_x_prime = four_dim_space_with_energy[1][:N_compton] + theta_x[:N_compton]
    photon_y = four_dim_space_with_energy[2][:N_compton]
    photon_y_prime = four_dim_space_with_energy[3][:N_compton] + theta_y[:N_compton]
    
    return photon_x, photon_x_prime, photon_y, photon_y_prime

def hole_colimation(L, hole, largest_diameter, smallest_diameter, rect_x_min, rect_x_max, rect_y_min, rect_y_max):
    CBS_angle_with_energy, w0_gauss, energy_distribution, four_dim_space_with_energy, N_compton = compton_backscattering(E_0, w_0, L, r_0)
        
    photon_x, photon_x_prime, photon_y, photon_y_prime = colimation_cycle(CBS_angle_with_energy, N_compton, four_dim_space_with_energy)

    photon_x_plane = np.vstack((photon_x, photon_x_prime))
    photon_y_plane = np.vstack((photon_y, photon_y_prime))
    photon_4d_distr = np.vstack((photon_x_plane, photon_y_plane))
    photon_4d_distr_with_energy = np.vstack((photon_4d_distr, CBS_angle_with_energy[1]))
    photon_4d_distr_with_energy_transpose = photon_4d_distr_with_energy.T
    np.savetxt("photons_4d_distr.txt", photon_4d_distr_with_energy_transpose)
    
    empty_gap_matrix = np.array([[1, L, 0, 0], [0, 1, 0, 0], [0, 0, 1, L], [0, 0, 0, 1]])
    photons_4d_distr_after_gap = empty_gap_matrix @ photon_4d_distr

    mask_for_circle_colimation = (photons_4d_distr_after_gap[0]**2 + photons_4d_distr_after_gap[2]**2 < hole**2)
    x_after_circle_colimation = photons_4d_distr_after_gap[0][mask_for_circle_colimation]
    y_after_circle_colimation = photons_4d_distr_after_gap[2][mask_for_circle_colimation]
    x_prime_after_circle_colimation = photons_4d_distr_after_gap[1][mask_for_circle_colimation]
    y_prime_after_circle_colimation = photons_4d_distr_after_gap[3][mask_for_circle_colimation]

    x_plane_after_circle_colimation = np.vstack((x_after_circle_colimation, x_prime_after_circle_colimation))
    y_plane_after_circle_colimation = np.vstack((y_after_circle_colimation, y_prime_after_circle_colimation))
    photons_after_circle_colimation = np.vstack((x_plane_after_circle_colimation, y_plane_after_circle_colimation))
    CBS_energy_after_circle_colimation = CBS_angle_with_energy[1][mask_for_circle_colimation]
    
    mask_for_ellipsoidal_colimation = (photons_4d_distr_after_gap[0]**2 / largest_diameter**2 + photons_4d_distr_after_gap[2]**2 / smallest_diameter**2 < 1)
    x_after_ellipsoidal_colimation = photons_4d_distr_after_gap[0][mask_for_ellipsoidal_colimation]
    y_after_ellipsoidal_colimation = photons_4d_distr_after_gap[2][mask_for_ellipsoidal_colimation]
    x_prime_after_ellipsoidal_colimation = photons_4d_distr_after_gap[1][mask_for_ellipsoidal_colimation]
    y_prime_after_ellipsoidal_colimation = photons_4d_distr_after_gap[3][mask_for_ellipsoidal_colimation]

    x_plane_after_ellipsoidal_colimation = np.vstack((x_after_ellipsoidal_colimation, x_prime_after_ellipsoidal_colimation))
    y_plane_after_ellipsoidal_colimation = np.vstack((y_after_ellipsoidal_colimation, y_prime_after_ellipsoidal_colimation))
    photons_after_ellipsoidal_colimation = np.vstack((x_plane_after_ellipsoidal_colimation, y_plane_after_ellipsoidal_colimation))
    CBS_energy_after_ellipsoidal_colimation = CBS_angle_with_energy[1][mask_for_ellipsoidal_colimation]

    mask_for_rectangle_colimation = ((photons_4d_distr_after_gap[0] > rect_x_min) & (photons_4d_distr_after_gap[0] < rect_x_max) & (photons_4d_distr_after_gap[2] > rect_y_min) & (photons_4d_distr_after_gap[2] < rect_y_max))
    x_after_rectangle_colimation = photons_4d_distr_after_gap[0][mask_for_rectangle_colimation]
    y_after_rectangle_colimation = photons_4d_distr_after_gap[2][mask_for_rectangle_colimation]
    x_prime_after_rectangle_colimation = photons_4d_distr_after_gap[1][mask_for_rectangle_colimation]
    y_prime_after_rectangle_colimation = photons_4d_distr_after_gap[3][mask_for_rectangle_colimation]

    x_plane_after_rectangle_colimation = np.vstack((x_after_rectangle_colimation, x_prime_after_rectangle_colimation))
    y_plane_after_rectangle_colimation = np.vstack((y_after_rectangle_colimation, y_prime_after_rectangle_colimation))
    photons_after_rectangle_colimation = np.vstack((x_plane_after_rectangle_colimation, y_plane_after_rectangle_colimation))
    CBS_energy_after_rectangle_colimation = CBS_angle_with_energy[1][mask_for_rectangle_colimation]

    fig, axes = plt.subplots(3, 3, figsize=(12,8))
    
    axes[0][0].scatter(photon_4d_distr[0], photon_4d_distr[1], s=2, alpha=0.4)
    axes[0][0].scatter(photon_4d_distr[2], photon_4d_distr[3], s=2, alpha=0.4)
    axes[0][0].set_title("photons distribution")
    #axes[0][0].set_aspect('equal')

    axes[0][1].scatter(four_dim_space_with_energy[0], four_dim_space_with_energy[1], s=2, alpha=0.4)
    axes[0][1].scatter(four_dim_space_with_energy[2], four_dim_space_with_energy[3], s=2, alpha=0.4)
    axes[0][1].set_title("electrons distribution")
    #axes[0][1].set_aspect('equal')

    axes[0][2].scatter(photons_4d_distr_after_gap[0], photons_4d_distr_after_gap[2], s=2, alpha=0.4)
    axes[0][2].set_title("photons distribution before collimation, (x,y) plane")
    axes[0][2].set_aspect('equal', adjustable='box')

    axes[1][0].scatter(photons_after_circle_colimation[0], photons_after_circle_colimation[2], s=2, alpha=0.4)
    axes[1][0].set_title("photons distribution after circle collimation in (x,y) plane")
    axes[1][0].set_aspect('equal', adjustable='box')

    axes[1][1].scatter(photons_after_ellipsoidal_colimation[0], photons_after_ellipsoidal_colimation[2], s=2, alpha=0.4)
    axes[1][1].set_title("photons distribution after ellipsoidal collimation in (x,y) plane")
    axes[1][1].set_aspect('equal', adjustable='box')

    axes[1][2].scatter(photons_after_rectangle_colimation[0], photons_after_rectangle_colimation[2], s=2, alpha=0.4)
    axes[1][2].set_title("photons distribution after rectangle collimation in (x,y) plane")
    #axes[1][2].set_aspect('equal', adjustable='box')
    """
    bins_circle = bin_size_calculator(CBS_energy_after_circle_colimation)
    print(bins_circle)
    bins_ellipsoidal = bin_size_calculator(CBS_energy_after_ellipsoidal_colimation)
    bins_rectangle = bin_size_calculator(CBS_energy_after_rectangle_colimation)
    """
    bins_size = np.histogram_bin_edges(CBS_angle_with_energy[1], bins=100)
    axes[2][0].hist(CBS_angle_with_energy[1], bins=bins_size, histtype='step', fill=False)
    axes[2][0].hist(CBS_energy_after_circle_colimation, bins=bins_size)
    axes[2][0].set_title("CBS energy after circle collimation")

    axes[2][1].hist(CBS_angle_with_energy[1], bins=bins_size, histtype='step', fill=False)
    axes[2][1].hist(CBS_energy_after_ellipsoidal_colimation, bins=bins_size)
    axes[2][1].set_title("CBS energy after ellipsoidal collimation")

    axes[2][2].hist(CBS_angle_with_energy[1], bins=bins_size, histtype='step', fill=False)
    axes[2][2].hist(CBS_energy_after_rectangle_colimation, bins=bins_size)
    axes[2][2].set_title("CBS energy after rectangle collimation")

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    tex0 = 0.0001779/10
    x0 = 0.0001686/10
    tey0 = 2.23998e-06
    y0 = 1.339297e-05/10
    E_0 = 2.0e9 #eV
    w_0 = h * c / 527e-9 / e
    r_0 = 2.8179e-15

    L = 50. #m

    hole = 3.e-3 #m
    largest_diameter = 3.5e-3
    smallest_diameter = 2.e-3
    rect_x_min = -3.e-3
    rect_x_max = 3e-3
    rect_y_min = -2.e-3
    rect_y_max = 2.e-3
    N_gen = int(1e5)

    beta_x = 1. #m
    beta_y = 0.15 #m
    eps_x = 0.25e-13 #m*rad
    eps_y = 0.05 * eps_x #m*rad
    alpha_x = 5.e-2 #m
    alpha_y = 3.e-3 #m

    compton_backscattering(E_0, w_0, L, r_0)
    hole_colimation(L, hole, largest_diameter, smallest_diameter, rect_x_min, rect_x_max, rect_y_min, rect_y_max)
