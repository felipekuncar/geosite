"""
geosite / class_CPTu

Creates a CPT profile object.
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
# CLASS CPTu
# =====================================================================================================================

class CPTu(object):
    """
    A cone penetration test with pore pressure measurement (CPTu) object

    Parameters
    ----------
    depth (array_like):
        Depth - in (m)
    qc (array_like):
        Cone resistance - in (kPa)
    fs (array_like):
        Sleeve friction resistance - in (kPA)
    u2 (array_like):
        Pore pressure, measured just behind the cone - in (kPa)
    GWL (float (default = 0.00)):
        Groundwater level - in (m)
    a (float (default = 0.80)):
        Net area ratio for cone - dimensionless
    """

    def __init__(self, depth, qc, fs, u2, GWL=0.00, a=0.80):
        self.depth = np.array(depth)
        self.qc_orig = np.array(qc)
        self.qc = np.where(np.array(qc) > 0, np.array(qc), 0.01)
        self.fs_orig = np.array(fs)
        self.fs = np.where(np.array(fs) > 0, np.array(fs), 0.01)
        self.u2 = np.array(u2)
        self.GWL = GWL if GWL > 0 else -GWL
        self.a = a
        self._sigma_v0 = None
        self._u0 = None

    @property
    def qt(self):
        """ Corrected cone resistance (corrected for pore water effects) - in (kPa) """
        qt = self.qc + self.u2 * (1 - self.a)
        return qt

    @property
    def Rf(self):
        """ Friction ratio - in (%) """
        Rf = (self.fs / self.qt) * 100
        return Rf

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
    def gamma(self):
        """
        It estimates the soil total unit weight - in (kN/m3)

        According to Robertson & Cabal (2010)

        References
        ----------
        Robertson P.K., Cabal K.L. (2010). Estimating soil unit weight from CPT. 2nd International Symposium on Cone Penetration Test, CPT'10, Huntington Beach, CA, USA.
        """
        gamma = (0.27 * np.log10(self.Rf) + 0.36 * np.log10(self.qt / pa) + 1.236) * gamma_w
        # If the values of qc or fs are zero, negative or non-existent, then enforce gamma = 18 (kN/m3)
        gamma = np.where((self.qc == 0.01) | (self.fs == 0.01), 18, gamma)
        return gamma

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
    def Qt(self):
        """ Normalized cone penetration resistance - dimensionless """
        Qt = (self.qt - self.sigma_v0) / self.sigma_v0_eff
        return Qt

    @property
    def Fr(self):
        """ Normalized friction ratio - in (%) """
        Fr = (self.fs / (self.qt - self.sigma_v0)) * 100
        return Fr

    def Ic(self, version=1):
        """
        Soil behaviour type index - dimensionless

        Parameters
        ----------
        version: str, optional (default = 1)
            If version = 1, it computes the normalized soil behaviour type index,
            according to Robertson (2009)
            If version = 2, it computes the non-normalized soil behaviour type index,
            according to Robertson (2010)

        """
        if version == 1:
            Ic = np.zeros(len(self.depth))
            Qtn = np.zeros(len(self.depth))
            for i in range(0, len(self.depth)):
                Fr = self.Fr[i]
                qt = self.qt[i]
                sigma_v0 = self.sigma_v0[i]
                sigma_v0_eff = self.sigma_v0_eff[i]
                def equations(p):
                    Ic_, Qtn_, n_ = p
                    eq1 = Ic_ - (((3.47 - np.log10(Qtn_)) ** 2 + (np.log10(Fr) + 1.22) ** 2) ** 0.5)
                    eq2 = Qtn_ - (((qt - sigma_v0) / pa) * ((pa / sigma_v0_eff) ** n_))
                    eq3 = n_ - min(0.381 * Ic_ + 0.05 * (sigma_v0_eff / pa) - 0.15, 1)
                    return (eq1, eq2, eq3)
                Ic_, Qtn_, n_ = fsolve(equations, [1.80, 100.00, 0.50])
                Ic[i] = Ic_
                Qtn[i] = Qtn_

        if version == 2:
            Ic = ((3.47 - np.log10(self.qt / pa)) ** 2 + (np.log10(self.Rf) + 1.22) ** 2) ** 0.5
            Qtn = None

        # Only consider Ic if qc and fs are realistic values (existent values greater than 0)
        Ic = np.where((self.qc_orig > 0) & (self.fs_orig > 0), Ic, float('NaN'))

        return Ic, Qtn

    def Vs(self, Ic, correlation=1):
        """
        It estimates the shear-wave velocity of the soil, Vs - in (m/s)

        Parameters
        ----------
        correlation: str, optional (default = 1)
            If version = 1, it computes Vs according to McGann et al. (2015)
            If version = 2, it computes Vs according to Robertson (2009) / CPT Guide 6th Edition (2015)

        References
        ----------
        McGann, C. R., Bradley, B. A., Taylor, M. L., Wotherspoon, L. M., & Cubrinovski, M. (2015). Development of an
        empirical correlation for predicting shear wave velocity of Christchurch soils from cone penetration test data.
        Soil Dynamics and Earthquake Engineering, 75, 66–75.

        Robertson, P.K. (2009). Interpretation of cone penetration tests - a unified approach. Can. Geotech. J. 46: 1337-1355.

        Robertson, P.K., Cabal K,L. (2015). Guide to cone penetration testing for geotechnical engineering, 6th Edition,
        GREGG.
        """
        if correlation == 1:
            Vs = 18.4 * (self.qt ** 0.144) * (self.fs ** 0.0823) * (self.depth ** 0.278)

        elif correlation == 2:
            alpha = 10 ** (0.55 * Ic + 1.68)
            Vs = (alpha * (self.qt - self.sigma_v0) / pa) ** 0.5

        Vs = np.where((self.qc_orig < 0) | (self.fs_orig < 0), float('NaN'), Vs)

        return Vs

    def relative_density(self, Ic, Qtn, correlation=1):
        """
        It estimates the relative density of the soil, Dr - in (%)

        Parameters
        ----------
        correlation: str, optional (default = 1)
            If correlation = 1, it computes Dr according to Idriss & Boulanger (2008)
            If correlation = 2, it computes Dr according to Kulhawy & Mayne (1990) - simplified expression from CPT
            Guide 6th Edition (2015)

        References
        ----------
        Idriss, I. M., and Boulanger, R. W. (2008) "Soil Liquefaction during Earthquakes." MNO-12, Earthquake
        Engineering Research Institute, Oakland, CA.

        Kulhawy, F.H., Mayne, P.H., 1990. Manual on estimating soil properties for foundation design, Report EL-6800
        Power Research Institute, EPRI, August 1990.

        Robertson, P.K., Cabal K,L. (2015). Guide to cone penetration testing for geotechnical engineering, 6th Edition,
        GREGG.
        """
        if correlation == 1:
            relative_density = np.zeros(len(self.depth))
            for i in range(0, len(self.depth)):
                qt = self.qt[i]
                sigma_v0_eff = self.sigma_v0_eff[i]
                def equations(p):
                    Dr, Cn, m = p
                    eq1 = Dr - (0.478 * ((Cn * qt / pa) ** 0.264) - 1.063)
                    eq2 = Cn - (min((pa / sigma_v0_eff) ** m, 1.70))
                    eq3 = m - (0.784 - (0.521 * Dr))
                    return (eq1, eq2, eq3)
                Dr, Cn, m = fsolve(equations, [0.50, 1.50, 0.52])
                if Ic[i] <= 2.6:
                    relative_density[i] = min(Dr * 100, 100)
                else:
                    relative_density[i] = None

        if correlation == 2:
            relative_density = np.where(Ic <= 2.6, np.where(np.sqrt(Qtn / 350) <= 1.00, (np.sqrt(Qtn / 350)), 1.00) * 100, float('NaN'))

        # Only consider Dr if qc and fs are realistic values (existent values greater than 0)
        relative_density = np.where((self.qc_orig > 0) & (self.fs_orig > 0), relative_density, float('NaN'))

        return relative_density

    def friction_angle(self, Ic, Qtn, correlation=1):
        """
        It estimates the friction angle of the soil, Dr - in (°)

        Parameters
        ----------
        correlation: str, optional (default = 1)
            If correlation = 1, it computes the peak friction angle according to Robertson (2010)
            If correlation = 2, it computes the peak friction angle according to Kulhawy & Mayne (1990)
            If correlation = 3, it computes the peak friction angle according to Robertson & Campanella (1983)

        References
        ----------
        Kulhawy, F.H., Mayne, P.H. (1990). Manual on estimating soil properties for foundation design, Report EL-6800
        Power Research Institute, EPRI, August 1990.

        Robertson, P.K., Campanella, R.G. (1983). Interpretation of cone penetration tests – Part I (sand). Canadian
        Geotechnical Journal, 20(4):718-733.

        """
        if correlation == 1:
            # correction factor
            Kc = np.where(Ic <= 1.64, 1.00, 5.581 * (Ic ** 3) - 0.403 * (Ic ** 4) - 21.63 * (Ic ** 2) + 33.75 * Ic - 17.88)
            # equivalent clean sand normalized cone resistance
            Qtncs = Kc * Qtn
            # constant volume (or critical state) friction angle
            # typically about 33 (°) for quartz sands but can be as high as 40 (°) for felspathic sand
            phi_cv = 33 * np.ones(len(self.depth))
            # peak friction angle
            friction_angle = np.where(Ic < 2.60, phi_cv + 15.84 * np.log10(Qtncs) - 26.88 * np.ones(len(self.depth)), float('NaN'))

        elif correlation == 2:
            friction_angle = np.where(Ic < 2.60, 17.6 + 11 * np.log10(Qtn), float('NaN'))

        elif correlation == 3:
            friction_angle = np.where(Ic < 2.6, np.arctan((1/2.68) * (np.log10(self.qt / self.sigma_v0_eff) + 0.29)) * 360 / (2 * np.pi), float('NaN'))

        # Only consider friction_angle if qc and fs are realistic values (existent values greater than 0)
        friction_angle = np.where((self.qc_orig > 0) & (self.fs_orig > 0), friction_angle, float('NaN'))

        return friction_angle

    def plasticityIndex(self, Ic, correlation=1):
        """
        It estimates the plasticity index of the soil, PI - in (%)

        Parameters
        ----------
        correlation: str, optional (default = 1)
            If correlation = 1, it computes PI according to Cetin & Ozan (2009)
            If correlation = 2, it computes PI according to a simple rule based on engineering judgement.

        References
        ----------
        Cetin, K.O., Ozan, C. (2009). CPT-Based Probabilistic Soil Characterization and Classification. Journal of
        Geotechnical and Geoenvironmental Engineering, 135(1): 84-107.
        """

        if correlation == 1:

            '''
            # Power law stress normalization exponent, c, as recommended by Cetin & Isik (2007)
            # Observation: Iteration for c is not working correctly. A simpler implementation is adopted for now.
            c = np.zeros(len(self.depth))
            for i in range(0, len(self.depth)):
                def equations(p):
                    qt1net_, c_, R_ = p
                    eq1 = qt1net_ - (self.qt[i] - self.sigma_v0[i]) / ((self.sigma_v0_eff[i] / pa) ** c_)
                    eq2 = c_ - (R_ - 272.38) / (275.19 - 272.38)
                    eq3 = R_ - np.sqrt((np.log10(self.Fr[i]) + 243.91) ** 2 + (np.log10(qt1net_ / pa) - 126.24) ** 2)
                    return (eq1, eq2, eq3)
                qt1net_, c_, R_ = fsolve(equations, [2000, 0.75, 273.5])
                print('c_:', c_)
                c_ = max(0.25, c_)
                c_ = min(1.00, c_)
                c[i] = c_
                print('c_corrected:', c_)

            # Normalized net cone tip resistance, qt1net
            denominator = ((self.sigma_v0_eff / pa) ** c)
            denominator = np.where(denominator > 2, 2, denominator)
            qt1net = (self.qt - self.sigma_v0) / denominator

            # qt1net has to be entered in MPa, and equation was not correct in the original paper
            # (denominator 2.25 goes outside - see Rockscience manual)
            PI = 10 ** ((2.37 + 1.33 + np.log10(self.Fr) - np.log10(qt1net / 1000)) / 2.25)
            print(PI)

            plt.plot(self.depth, PI)
            plt.plot(self.depth, c)
            plt.show()
            '''

        # Alternatively, simpler implementation:
    
        if correlation == 1:

            # Power law stress normalization exponent, c, as recommended by Robertson (1999) - simpler implementation
            Ic = np.array(Ic)
            c = np.where(Ic >= 2.60, 1, np.where(Ic >= 2.05, 0.75, 0.50))

            # Normalized net cone tip resistance, qt1net
            denominator = ((self.sigma_v0_eff / pa) ** c)
            denominator = np.where(denominator > 2, 2, denominator)
            qt1net = (self.qt - self.sigma_v0) / denominator

            # qt1net has to be entered in MPa, and equation was not correct in the original paper
            # (denominator 2.25 goes outside - see Rockscience manual)
            PI = 10 ** ((2.37 + 1.33 + np.log10(self.Fr) - np.log10(qt1net / 1000)) / 2.25)

        elif correlation == 2:

            PI = np.zeros(len(Ic))
            for i in range(len(Ic)):
                if Ic[i] < 2.05:
                    PI[i] = 0
                elif 2.05 <= Ic[i] < 2.60:
                    PI[i] = 5
                elif 2.60 <= Ic[i] < 2.95:
                    PI[i] = 12
                else:
                    PI[i] = 20

        return PI

    def undrainedShearStrength(self, Ic):
        """
        Estimate the undrained shear strength (Su) - in (kPa)

        References
        ----------
        Robertson, P.K., Cabal K,L. (2022). Guide to cone penetration testing, 7th Edition,
        GREGG.

        """

        N_kt = 14 * np.ones(len(Ic))
        Su = np.where(Ic < 2.60, 0, (self.qt - self.sigma_v0) / N_kt)

        # Rate-effect correction factor (Vardanega & Bolton, 2013)
        Su = 1.2 * Su

        return Su


