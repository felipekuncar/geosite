"""
geosite / class_Vs.py

Creates a Vs profile object.
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

# =====================================================================================================================
# DEFINE CONSTANTS
# =====================================================================================================================

# Atmospheric pressure (kPa)
pa = 101.325

# Unit weight of water (kN/m3)
gamma_w = 9.80665

# Gravity (m/s2)
g = 9.80665

# Poisson's ratio of soil
nu = 0.25

# =====================================================================================================================
# CLASS Vs
# =====================================================================================================================

class Vs(object):
    """
    A shear wave velocity (Vs) profile object

    Parameters
    ----------
    depth (array_like):
        Depth - in (m)
    Vs_values (array_like):
        Cone resistance - in (kPa)
    GWL (float (default = 0.00)):
        Groundwater level - in (m)
    """

    def __init__(self, depth, values, GWL=0.00):
        self.depth = np.array(depth)
        self.values = np.array(values)
        self.GWL = GWL if GWL > 0 else -GWL
        self._sigma_v0 = None
        self._u0 = None

    @property
    def gamma(self):
        """
        It estimates the soil total unit weight - in (kN/m3)

        According to Rix et al. (2019)

        References
        ----------
        Rix et al. (2019). Manual on subsurface investigations. Web-Only Document 258, NCHRP Project 21-10.

        """

        gamma = 7.83 * np.log10(self.values) - 0.125 * np.ones(len(self.depth))

        return gamma

    def gamma_corr(self, correlation=1):
        """
        Estimate the soil total unit weight - in (kN/m3)

        Parameters
        ----------
        correlation: str, optional (default = 1)
            If correlation = 1, it computes the soil total unit weight according to Rix et al. (2019).
            If correlation = 2, it computes the soil total unit weight according to Boore (2015).

        References
        ----------
        Rix et al. (2019). Manual on subsurface investigations. Web-Only Document 258, NCHRP Project 21-10.

        Boore D. (2015). Notes on relating density to velocity for use in site amplification calculations.

        """
        if correlation == 1:
            gamma_corr = 7.83 * np.log10(self.values) - 0.125 * np.ones(len(self.depth))

        elif correlation == 2:
            # conversion of Vs from (m/s) to (km/s)
            Vs = self.values / 1000
            # density in (g/cm3) - equation valid for Vs < 3550 (m/s)
            rho = np.where(Vs < 0.30, 1 + ((1.53 * (Vs ** 0.85))/(0.35 + 1.889 * (Vs ** 1.7))), 1.74 * (0.9409 + 2.0947 * Vs - 0.8206 * (Vs ** 2) + 0.2683 * (Vs ** 3) - 0.0251 * (Vs ** 4)) ** 0.25)
            # soil total unit weight (kN/m3)
            gamma_corr = rho * 9.81

        return gamma_corr

    @property
    def u0(self):
        """ In-situ equilibrium pore pressure - in (%) """
        if self._u0 is None:
            u0 = np.zeros(len(self.depth))
            mask = self.depth < self.GWL
            u0[mask] = 0.0
            u0[~mask] = gamma_w * (self.depth[~mask] - self.GWL)
            self._u0 = u0
        return self._u0

    @property
    def sigma_v0(self):
        """ In-situ vertical total stress - in (kPa) """
        if self._sigma_v0 is None:
            sigma_v0 = np.zeros(len(self.depth))
            sigma_v0[0] = sigma_v0[0] + self.gamma[0] * self.depth[0]
            for i in range(1, len(self.depth)):
                sigma_v0[i] = sigma_v0[i - 1] + self.gamma[i] * (self.depth[i] - self.depth[i - 1])
            self._sigma_v0 = sigma_v0
        return self._sigma_v0

    @property
    def sigma_v0_eff(self):
        """ In-situ vertical effective stress - in (kPa) """
        sigma_v0_eff = self.sigma_v0 - self.u0
        return sigma_v0_eff

    @property
    def N1_60(self):
        """
        Corrected SPT N-Values (N1)60 - dimensionless

        Parameters
        ----------
        correlation: str, optional (default = 1)
            If correlation = 1, it estimates the N60 according to Rix et al. (2019).

        References
        ----------
        Rix et al. (2019). Manual on subsurface investigations. Web-Only Document 258, NCHRP Project 21-10.
        """

        N_60 = (self.values/97) ** (1 / 0.314)
        C_N = np.where((pa/self.sigma_v0_eff) ** 0.5 <= 2, (pa/self.sigma_v0_eff) ** 0.5, 2)
        N1_60 = C_N * N_60

        N1_60_avg = np.zeros(len(N1_60))
        # Average values of each layer
        for i in range(0, len(N1_60)):
            if (i % 2) == 0:
                N1_60_avg[i] = (0.5 * N1_60[i] + 0.5 * N1_60[i+1])
            else:
                N1_60_avg[i] = (0.5 * N1_60[i] + 0.5 * N1_60[i-1])
        N1_60_avg[-1] = N1_60[-1]

        N1_60 = N1_60_avg

        return N1_60

    def relative_density(self, N1_60, correlation=1):
        """
        It estimates the relative density of the soil, Dr - in (%)

        Parameters
        ----------
        correlation: str, optional (default = 1)
            If correlation = 1, it computes Dr according to Idriss & Boulanger (2008)
            If correlation = 2, it computes Dr according to Kulhawy & Mayne (1990)

        References
        ----------
        Idriss, I. M., and Boulanger, R. W. (2008) "Soil Liquefaction during Earthquakes." MNO-12, Earthquake
        Engineering Research Institute, Oakland, CA.

        Kulhawy, F.H., Mayne, P.H., 1990. Manual on estimating soil properties for foundation design, Report EL-6800
        Power Research Institute, EPRI, August 1990.

        """
        if correlation == 1:
            relative_density = np.sqrt(N1_60 / 46) * 100

        return relative_density

    def friction_angle(self, N1_60, correlation=1):
        """
        It estimates the friction angle of the soil - in (°)

        Parameters
        ----------
        correlation: str, optional (default = 1)
            If correlation = 1, it estimates the friction angle according to Rix et al. (2019).

        References
        ----------
        Rix et al. (2019). Manual on subsurface investigations. Web-Only Document 258, NCHRP Project 21-10.
        """

        if correlation == 1:
            friction_angle = 20 + np.sqrt(15.4 * N1_60)

        return friction_angle
