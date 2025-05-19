"""
geosite / class_PDMY02.py

Defines the PDMY02 class, associated with the PDMY02 OpenSees constitutive model.
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
# CLASS PDMY02
# =====================================================================================================================

class PDMY02(object):
    """
    A PDMY02 object that allows the estimation of model-specific parameters from the relative density values, according
    to the calibration of Karimi & Dashti (2015) for Nevada sand

    Parameters
    ----------
    depth (array_like):
        Depth - in (m)
    Dr (array_like):
        Relative density - in (%)

    References
    ----------
    Karimi, Z., & Dashti, S. (2016). Numerical and Centrifuge Modeling of Seismic Soil–Foundation–Structure Interaction
    on Liquefiable Ground. In Journal of Geotechnical and Geoenvironmental Engineering (Vol. 142, Issue 1, p. 04015061).
    American Society of Civil Engineers (ASCE).
    """

    def __init__(self, depth, Dr):
        self.depth = np.array(depth)
        self.Dr = np.array(Dr)

    def phiAng(self):
        """
        Estimate the triaxial friction angle used by the model - in (°)

        """
        phiAng = np.where(self.Dr < 30.0,  31.0,
                np.where(self.Dr < 40.0,  31.0 + 0.1000 * (self.Dr - 30.0),
                np.where(self.Dr < 50.0,  32.0 + 0.1500 * (self.Dr - 40.0),
                np.where(self.Dr < 63.0,  33.5 + 0.0768 * (self.Dr - 50.0),
                np.where(self.Dr < 68.0,  34.5 + 0.3000 * (self.Dr - 63.0),
                np.where(self.Dr < 90.0,  36.0 + 0.1818 * (self.Dr - 68.0), 40.0))))))

        return phiAng

    def PTAng(self):
        """
        Estimate the phase transformation angle - in (°)

        """
        PTAng = np.where(self.Dr < 30.0,  31.0,
                np.where(self.Dr < 40.0,  31.0 + -0.1000 * (self.Dr - 30.0),
                np.where(self.Dr < 50.0,  30.0 + -0.4500 * (self.Dr - 40.0),
                np.where(self.Dr < 63.0,  25.5 +  0.0846 * (self.Dr - 50.0),
                np.where(self.Dr < 68.0,  26.5 + -0.1200 * (self.Dr - 63.0),
                np.where(self.Dr < 90.0,  26.0 +  0.0227 * (self.Dr - 68.0), 26.5))))))

        return PTAng

    def e(self):
        """
        It estimates the void ratio - dimensionless

        """
        e = np.where(self.Dr < 30.0,  0.76,
            np.where(self.Dr < 40.0,  0.76 + -0.0030 * (self.Dr - 30.0),
            np.where(self.Dr < 50.0,  0.73 + -0.0030 * (self.Dr - 40.0),
            np.where(self.Dr < 63.0,  0.70 + -0.0031 * (self.Dr - 50.0),
            np.where(self.Dr < 68.0,  0.66 + -0.0020 * (self.Dr - 63.0),
            np.where(self.Dr < 90.0,  0.65 + -0.0032 * (self.Dr - 68.0), 0.58))))))

        return e

    def contrac1(self):
        """
        It estimates the parameter contrac1

        """
        contrac1 = np.where(self.Dr < 30.0,  0.087,
                   np.where(self.Dr < 40.0,  0.087 + -0.0020 * (self.Dr - 30.0),
                   np.where(self.Dr < 50.0,  0.067 + -0.0017 * (self.Dr - 40.0),
                   np.where(self.Dr < 63.0,  0.050 + -0.0008 * (self.Dr - 50.0),
                   np.where(self.Dr < 68.0,  0.040 + -0.0040 * (self.Dr - 63.0),
                   np.where(self.Dr < 90.0,  0.020 + -0.0002 * (self.Dr - 68.0), 0.016))))))

        return contrac1

    def contrac2(self):
        """
        It estimates the parameter contrac2

        """
        contrac2 = np.where(self.Dr < 30.0,  5.00,
                   np.where(self.Dr < 40.0,  5.00 + -0.0500 * (self.Dr - 30.0),
                   np.where(self.Dr < 50.0,  4.50 + -0.0500 * (self.Dr - 40.0),
                   np.where(self.Dr < 63.0,  4.50 + -0.1154 * (self.Dr - 50.0),
                   np.where(self.Dr < 68.0,  2.50 + -0.2000 * (self.Dr - 63.0),
                   np.where(self.Dr < 90.0,  1.50 + -0.0023 * (self.Dr - 68.0), 1.40))))))

        return contrac2

    def contrac3(self):
        """
        It estimates the parameter contrac3

        """
        contrac3 = np.where(self.Dr < 30.0,  0.30,
                   np.where(self.Dr < 40.0,  0.30 + -0.0030 * (self.Dr - 30.0),
                   np.where(self.Dr < 50.0,  0.27 + -0.0020 * (self.Dr - 40.0),
                   np.where(self.Dr < 63.0,  0.25 + -0.0038 * (self.Dr - 50.0),
                   np.where(self.Dr < 68.0,  0.20 + -0.0100 * (self.Dr - 63.0),
                   np.where(self.Dr < 90.0,  0.15 + -0.0005 * (self.Dr - 68.0), 0.14))))))

        return contrac3

    def dilat1(self):
        """
        It estimates the parameter dilat1

        """
        dilat1 = np.where(self.Dr < 30.0,  0.01,
                 np.where(self.Dr < 40.0,  0.01 + 0.0010 * (self.Dr - 30.0),
                 np.where(self.Dr < 50.0,  0.02 + 0.0040 * (self.Dr - 40.0),
                 np.where(self.Dr < 63.0,  0.06 + 0.0008 * (self.Dr - 50.0),
                 np.where(self.Dr < 68.0,  0.07 + 0.0160 * (self.Dr - 63.0),
                 np.where(self.Dr < 90.0,  0.15 + 0.0045 * (self.Dr - 68.0), 0.25))))))

        return dilat1

    def dilat2(self):
        """
        It estimates the parameter dilat2

        """
        dilat2 = 3 * np.ones(len(self.Dr))

        return dilat2

    def dilat3(self):
        """
        It estimates the parameter dilat3

        """
        dilat3 = np.zeros(len(self.Dr))

        return dilat3