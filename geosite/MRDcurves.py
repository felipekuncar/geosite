"""
geosite / MRDcurves.py

Functions used to create modulus reduction and damping curves.
"""

# =====================================================================================================================
# IMPORT LIBRARIES AND FUNCTIONS
# =====================================================================================================================

import numpy as np
from scipy.optimize import fsolve
import os.path
import pandas as pd
from math import sin, pi
import math
import os
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from GeoSite.class_PDMY02 import PDMY02

# Compute Dmin based on Darendeli (2001)
def getDminDarendeli(sigma_v0_eff, K0, PI, OCR, freq=1):
    """

    Purpose: Obtain Dmin array for a given site-response model

    Input
    ----------
    GWL: Freatic level (m)
    thickness: Thickness array (m)
    gamma: Specific weight array (kN/m3)
    K0: At-rest earth pressure coefficient array
    PI: Soil plasticity array (%)
    OCR: OCR array

    Output
    ----------
    Dmin: Small-starin material damping ratio (dimensionless)

    References
    ----------
    Darendeli M. (2001). Development of a new family of normalized modulus reduction and material damping curves.
    PhD Thesis, The University of Texas Austin.
    """

    PI = np.array(PI)
    K0 = np.array(K0)
    OCR = np.array(OCR)

    # Atmospheric pressure (kPa)
    pa = 101.325

    # In-situ mean effective stress (atm)
    sigma_0_eff = (1 + 2 * K0) * sigma_v0_eff / (3 * pa)

    ### Darendeli (2001)
    # Parameters
    phi_6 = 0.8005
    phi_7 = 0.0129
    phi_8 = -0.1069
    phi_9 = -0.2889
    phi_10 = 0.2919
    Dmin = (phi_6 + phi_7 * PI * (OCR ** phi_8)) * (sigma_0_eff ** phi_9) * (1 + phi_10 * np.log(freq))
    Dmin = Dmin / 100
    Dmin = np.array(Dmin)

    return Dmin

def Darendeli_MRDcurves(sigma_v0_eff, K0, PI, OCR, freq=1, N=10):
    '''
    Input parameters:
    sigma_v0_eff: In-situ vertical effective stress (kPa)
    K0: At-rest earth pressure coefficient
    PI: Plasticity index
    OCR: Overconsolidation ratio
    freq: Frequency (defaut = 1 Hz)
    N: Number of cycles (default = 10)

    Output:
    shearstrain: Array with 30 log-spaced shear strain values (%)
    G_ratio: Array with G/Gmax curve
    D: Array with damping curve (%)
    a: Curvature coefficient
    gamma_r: Pseudo-reference shear strain (%)
    Dmin: Small-starin material damping ratio (dimensionless)
    '''

    # Atmospheric pressure (kPa)
    pa = 101.325

    # shearStrain: Array with 30 shear strain values (%) for which the modulus reduction will be computed
    shearStrain = np.logspace(-5, 1, 30)

    # Model coefficients
    phi_1 = 0.0352
    phi_2 = 0.0010
    phi_3 = 0.3246
    phi_4 = 0.3483
    phi_5 = 0.919
    phi_6 = 0.8005
    phi_7 = 0.0129
    phi_8 = -0.1069
    phi_9 = -0.2889
    phi_10 = 0.2919
    phi_11 = 0.6329
    phi_12 = -0.0057
    a = phi_5
    c_1 = 0.2523 + 1.8618 * a - 1.1143 * (a ** 2)
    c_2 = -0.0095 - 0.0710 * a + 0.0805 * (a ** 2)
    c_3 = 0.0003 + 0.0002 * a - 0.0005 * (a ** 2)

    # Mean effective confining pressure (atm)
    sigma_0_eff = (1 + 2 * K0) * sigma_v0_eff / (3 * pa)

    # Pseudo-reference strain, gamma_r (%)
    gamma_r = (phi_1 + phi_2 * PI * (OCR ** phi_3)) * (sigma_0_eff ** phi_4)

    # G/Gmax
    G_ratio = 1 / (1 + ((shearStrain / gamma_r) ** a))

    # Minimum damping (%)
    Dmin = (phi_6 + phi_7 * PI * (OCR ** phi_8)) * (sigma_0_eff ** phi_9) * (1 + phi_10 * np.log(freq))

    # Masing damping for a = 1 (%)
    D_M_a1 = (100 / np.pi) * ((4 * (shearStrain - gamma_r * np.log((shearStrain + gamma_r) / gamma_r)) / ((shearStrain ** 2) / (shearStrain + gamma_r))) - 2)

    # Masing damping for a = phi_5 (%)
    D_M = c_1 * D_M_a1 + c_2 * (D_M_a1 ** 2) + c_3 * (D_M_a1 ** 3)

    # Damping curve (%)
    b = phi_11 + phi_12 * np.log(N)
    D = Dmin + b * D_M * (G_ratio ** 0.1)

    # Damping curve (dimensionless)
    D = D / 100

    return shearStrain, G_ratio, D, a, gamma_r, Dmin

def Yee_MRcurve(a, gamma_r, shearStrain, G_ratio_ref, gamma_1, G_max, tau_max):
    '''
    Purpose: Adjust modulus reduction curve to be consistent with shear strength at large strains, according to
    the procedure proposed by Yee et al. (2013).

    Reference:
    Yee, E., Stewart, J., and Tokimatsu, K. (2013). “Elastic and large-strain
    nonlinear seismic site response from analysis of vertical array recordings.”
    J. Geotech. Geoenviron. Eng., 139(10): 1789-1801.

    Input parameters:
    a: Curvature coefficient
    gamma_r: Pseudo-reference shear strain (%)
    shearStrain: Array of shear strain values (%) for which the modulus reduction will be computed
    G_ratio_ref: Array with reference (e.g., Darendeli) modulus reduction curve, computed for every strain value in shearStrain
    gamma_1: User-defined transitional shear strain (%)
    G_max: Small-strain shear modulus (MPa)
    tau_max: Shear strength (kPa)

    Output:
    G_ratio: Array with G/Gmax curve

    '''

    # Mpa to kPa
    G_max = G_max * 1000

    # % to unitless
    gamma_r = gamma_r / 100
    gamma_1 = gamma_1 / 100
    shearStrain = shearStrain / 100

    # Tangent Shear Modulus (G_gamma1) divided by G_max. Equation 6, Yee et al. (2013).
    G_gamma1_ratio = (1 + (1 - a) * ((gamma_1 / gamma_r) ** a)) / (1 + (gamma_1 / gamma_r) ** a) ** 2

    # Shear strees at gamma1 (tau_1). Page 1796, Yee et al. (2013).
    tau_1 =  (G_max * gamma_1) / (1 + (gamma_1 / gamma_r) ** a)

    # Warning
    if tau_1 <= 0.3 * tau_max:
        print('tau_1 <= 0.3 tau_max -> Proper value of gamma1')
    else:
        print('tau_1 > 0.3 tau_max -> It could be appropriate to modify (reduce) gamma1!')

    # Reference shear strain (unitless). Page 1796, Yee et al. (2013).
    gamma_ref_prime = (tau_max - tau_1) / (G_gamma1_ratio * G_max)

    # gamma_prime. Page 1796, Yee et al. (2013).
    gamma_prime = shearStrain - gamma_1 * np.ones(len(shearStrain))

    # gamma_prime divided by reference shear strain
    gamma_prime_ratio = gamma_prime / gamma_ref_prime

    # Initialize G_ratio
    G_ratio = np.zeros(len(G_ratio_ref))

    for i, strain in enumerate(shearStrain):

        if strain < gamma_1:

            G_ratio[i] = G_ratio_ref[i]

        else:

            # Equation 7, Yee et al. (2013).
            G_ratio[i] = ((gamma_1 / (1 + (gamma_1 / gamma_r) ** a)) + ((G_gamma1_ratio * gamma_prime[i]) / (1 + gamma_prime_ratio[i]))) / shearStrain[i]

    return G_ratio

def Darendeli_MRDcurves_forPySeismoSoil(sigma_v0_eff, K0, PI, OCR, freq=1, N=10):
    '''
    Input parameters:
    sigma_v0_eff: In-situ vertical effective stress (kPa)
    K0: At-rest earth pressure coefficient
    PI: Plasticity index
    OCR: Overconsolidation ratio
    freq: Frequency (defaut = 1 Hz)
    N: Number of cycles (default = 10)

    Output:
    shearstrain: Array with 30 log-spaced shear strain values (%)
    G_ratio: Array with G/Gmax curve
    D: Array with damping curve (%)
    a: Curvature coefficient
    gamma_r: Pseudo-reference shear strain (%)
    Dmin: Small-starin material damping ratio (dimensionless)
    '''

    # Atmospheric pressure (kPa)
    pa = 101.325

    # shearStrain: Array with shear strain values (%) for which the modulus reduction will be computed
    shearStrain = np.array([0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3])

    # Model coefficients
    phi_1 = 0.0352
    phi_2 = 0.0010
    phi_3 = 0.3246
    phi_4 = 0.3483
    phi_5 = 0.919
    phi_6 = 0.8005
    phi_7 = 0.0129
    phi_8 = -0.1069
    phi_9 = -0.2889
    phi_10 = 0.2919
    phi_11 = 0.6329
    phi_12 = -0.0057
    a = phi_5
    c_1 = 0.2523 + 1.8618 * a - 1.1143 * (a ** 2)
    c_2 = -0.0095 - 0.0710 * a + 0.0805 * (a ** 2)
    c_3 = 0.0003 + 0.0002 * a - 0.0005 * (a ** 2)

    # Mean effective confining pressure (atm)
    sigma_0_eff = (1 + 2 * K0) * sigma_v0_eff / (3 * pa)

    # Pseudo-reference strain, gamma_r (%)
    gamma_r = (phi_1 + phi_2 * PI * (OCR ** phi_3)) * (sigma_0_eff ** phi_4)

    # G/Gmax
    G_ratio = 1 / (1 + ((shearStrain / gamma_r) ** a))

    # Minimum damping (%)
    Dmin = (phi_6 + phi_7 * PI * (OCR ** phi_8)) * (sigma_0_eff ** phi_9) * (1 + phi_10 * np.log(freq))

    # Masing damping for a = 1 (%)
    D_M_a1 = (100 / np.pi) * ((4 * (shearStrain - gamma_r * np.log((shearStrain + gamma_r) / gamma_r)) / ((shearStrain ** 2) / (shearStrain + gamma_r))) - 2)

    # Masing damping for a = phi_5 (%)
    D_M = c_1 * D_M_a1 + c_2 * (D_M_a1 ** 2) + c_3 * (D_M_a1 ** 3)

    # Damping curve (%)
    b = phi_11 + phi_12 * np.log(N)
    D = Dmin + b * D_M * (G_ratio ** 0.1)

    # Damping curve (dimensionless)
    D = D / 100

    return shearStrain, G_ratio, D, a, gamma_r, Dmin