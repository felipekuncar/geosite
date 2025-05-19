"""
geosite / functions.py

Several functions used along with the geosite classes.
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

from geosite.class_PDMY02 import PDMY02

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
# FUNCTIONS FOR GENERATING TRIAL (OR INITIAL) MODELS
# =====================================================================================================================

def trialModelsFromCPTu(layerThickness, CPTu):

    # Obtain Ic
    Ic, Qtn = CPTu.Ic(version=1)

    # Obtain layersDepth from layerThickness
    layersDepth = np.zeros(1 + len(layerThickness))
    for i in range(1, len(layersDepth)):
        layersDepth[i] = layersDepth[i - 1] + layerThickness[i - 1]

    # Obtain counter that determine data within each layer
    j = 1
    counterLim = [0]
    for i in range(0, len(CPTu.depth)):
        if CPTu.depth[i] > layersDepth[j]:
            counterLim.append(i)
            j = j + 1

    # layer thicknesses considered
    layerThickness_trial = layerThickness[0:len(counterLim) - 1]

    # Layer numbers
    layerNumber_trial = np.arange(1, len(layerThickness_trial) + 1, 1)

    # Initialize arrays
    Ic_trial = []
    Qtn_trial = []
    Vs_trial1 = []
    Vs1 = CPTu.Vs(Ic, correlation=1)
    Vs_trial2 = []
    Vs2 = CPTu.Vs(Ic, correlation=2)
    gamma_trial = []
    gamma = CPTu.gamma
    Dr_trial1 = []
    Dr1 = CPTu.relative_density(Ic, Qtn, correlation=1)
    Dr_trial2 = []
    Dr2 = CPTu.relative_density(Ic, Qtn, correlation=2)
    PI_trial = []
    PI = CPTu.plasticityIndex(Ic, correlation=1)
    phi_trial1 = []
    phi1 = CPTu.friction_angle(Ic, Qtn, correlation=1)
    phi_trial2 = []
    phi2 = CPTu.friction_angle(Ic, Qtn, correlation=2)
    phi_trial3 = []
    phi3 = CPTu.friction_angle(Ic, Qtn, correlation=3)
    Su_trial = []
    Su = CPTu.undrainedShearStrength(Ic)

    for i in range(1, len(layerThickness_trial) + 1):
        Ic_trial.append(np.nanmean(Ic[counterLim[i-1]:counterLim[i]]))
        Qtn_trial.append(np.nanmean(Qtn[counterLim[i - 1]:counterLim[i]]))
        Vs_trial1.append(np.nanmean(Vs1[counterLim[i-1]:counterLim[i]]))
        Vs_trial2.append(np.nanmean(Vs2[counterLim[i - 1]:counterLim[i]]))
        gamma_trial.append(np.nanmean(gamma[counterLim[i - 1]:counterLim[i]]))
        Dr_trial1.append(np.nanmean(Dr1[counterLim[i - 1]:counterLim[i]]))
        Dr_trial2.append(np.nanmean(Dr2[counterLim[i - 1]:counterLim[i]]))
        PI_trial.append(np.nanmean(PI[counterLim[i - 1]:counterLim[i]]))
        phi_trial1.append(np.nanmean(phi1[counterLim[i - 1]:counterLim[i]]))
        phi_trial2.append(np.nanmean(phi2[counterLim[i - 1]:counterLim[i]]))
        phi_trial3.append(np.nanmean(phi3[counterLim[i - 1]:counterLim[i]]))
        Su_trial.append(np.nanmean(Su[counterLim[i - 1]:counterLim[i]]))

    # Logic tree
    Vs_trial1 = np.array(Vs_trial1)
    Vs_trial2 = np.array(Vs_trial2)
    Vs_trial = 0.75 * Vs_trial1 + 0.25 * Vs_trial2
    Dr_trial1 = np.array(Dr_trial1)
    Dr_trial2 = np.array(Dr_trial2)
    Dr_trial = 0.5 * Dr_trial1 + 0.5 * Dr_trial2
    phi_trial1 = np.array(phi_trial1)
    phi_trial2 = np.array(phi_trial2)
    phi_trial3 = np.array(phi_trial3)
    phi_trial = (1/3) * phi_trial1 + (1/3) * phi_trial2 + (1/3) * phi_trial3

    bottomDepth = np.zeros(len(layerThickness_trial))
    for i in range(len(layerThickness_trial)):
        if i == 0:
            bottomDepth[i] = layerThickness_trial[i]
        else:
            bottomDepth[i] = bottomDepth[i - 1] + layerThickness_trial[i]
    topDepth = bottomDepth - layerThickness_trial

    return layerNumber_trial, layerThickness_trial, topDepth, bottomDepth, Ic_trial, Qtn_trial, Vs_trial, gamma_trial, Dr_trial, phi_trial, Su_trial, PI_trial

# =====================================================================================================================
# FUNCTION readVsData
# =====================================================================================================================

def readVsData(format, PathVsSW, nonUniqueVsModel=None):

    # Format
    ##########################################
    # Format 1
    # .xlsx file
    # columns: Thickness (m) | Vs (m/s)
    ##########################################

    LayerThickness = []
    VsValue = []

    if format == 1:

        Vs_read = pd.read_excel(PathVsSW)
        # Assign class Vs
        Vs_thickness = Vs_read.iloc[:, 0]
        Vs_values = Vs_read.iloc[:, 1]

        Vs_depth = np.zeros(1 + len(Vs_thickness))
        for i in range(1, len(Vs_depth)):
            Vs_depth[i] = Vs_depth[i - 1] + Vs_thickness[i - 1]
        Vs_depth = np.repeat(Vs_depth, 2)
        Vs_depth = np.delete(Vs_depth, [0, len(Vs_depth) - 1])

        Vs_values = np.repeat(Vs_values, 2)

    return Vs_depth, Vs_values


# =====================================================================================================================
# FUNCTION PressDepVs
# =====================================================================================================================

# Created by Chris de la Torre
def PressDepVs(PathOutputs, thick, Vs, gamma, phi, siteID, pressDependCoe, refDepth, VsInversionTop, waterDepth):

    numLayers = len(thick)

    waterDepth = waterDepth if waterDepth > 0 else -waterDepth

    # For the effects of this computation lets consider a value of phi=40 for the bottom layer (half-space)
    # This is only considered for being able to plot the Vs profile
    # It doesn't matter the specific value of phi for this layer
    phi[-1] = 40

    # ------------------------------------------------------------------------------------------------------------------
    # Define the size of the arrays to be generated containing the reference parameters
    # ------------------------------------------------------------------------------------------------------------------

    # Define the size of the arrays to be generated
    # Size = number of layers
    vertStress = np.zeros(numLayers)  # 1) Effective vertical stress (kN/m2 = kPa)
    kNot = np.zeros(numLayers)        # 2) Coefficient of lateral earth pressure at rest
    refStress = np.zeros(numLayers)   # 3) Reference mean effective confining pressure (kN/m2 = kPa)
    Gref = np.zeros(numLayers)        # 4) Reference low-strain shear modulus (kN/m2 = kPa)

    # ------------------------------------------------------------------------------------------------------------------
    # Compute the reference parameters for the top layer
    # ------------------------------------------------------------------------------------------------------------------

    # 1) Effective vertical stress, vertStress (kN/m2)
    if thick[0] > waterDepth:  # just in the case where waterDepth = 0
        vertStress[0] = gamma[0] * thick[0] * refDepth[0] - gamma_w * (thick[0] * refDepth[0] - waterDepth)
    else:
        vertStress[0] = gamma[0] * thick[0] * refDepth[0]

    # 2) Coefficient of lateral earth pressure at rest, kNot
    kNot[0] = 1.0 - sin(phi[0] * 2.0 * pi / 360.0)

    # 3) Reference mean effective confining pressure, refStress (kN/m2)
    # p' = (sigma'_v + 2*sigma'_h)/3 = (sigma'_v + 2*KNot*sigma'_v)/3 = sigma'_v*(1 + 2*kNot)/3
    refStress[0] = vertStress[0] * (1.0 + 2.0 * kNot[0]) / 3.0

    # 4) Reference low-strain shear modulus, Gref (kN/m2 = kPa)
    # Vs_not (m/s) is the constant required for an exponential function of Vs to have equal travel time to a layer of
    # constant Vs
    Vs_not = Vs[0] / ((1.0 - pressDependCoe[0] / 2.0) * thick[0] ** (pressDependCoe[0] / 2))
    if VsInversionTop == 'No':
        Gref[0] = (gamma[0] / g) * (Vs_not ** 2.0) * (gamma[0] * (1.0 + 2.0 * kNot[0]) / (3.0 * refStress[0])) ** (-pressDependCoe[0])
        # Gref[0] = rho[0] * (Vs_not ** 2.0) * ((thick[0] * refDepth[0]) ** pressDependCoe[0])
    else:
        Gref[0] = (gamma[0] / g) * Vs[0] ** 2.0

    # ------------------------------------------------------------------------------------------------------------------
    # Compute the reference parameters for all other layers
    # ------------------------------------------------------------------------------------------------------------------

    # Compute depth to bottom and midpoint of each layer (m)
    bottomDepth = np.zeros(numLayers)
    midDepth = np.zeros(numLayers)
    bottomDepth[0] = thick[0]
    midDepth[0] = thick[0] / 2.0
    for i in range(1, numLayers):
        bottomDepth[i] = bottomDepth[i-1] + thick[i]
        midDepth[i] = (bottomDepth[i] + bottomDepth[i-1]) / 2.0

    # For each layer compute:
    # 1) Effective vertical stress, vertStress (kN/m2)
    # 2) Coefficient of lateral earth pressure at rest, kNot
    # 3) Reference mean effective confining pressure, refStress (kN/m2)
    for i in range(1, numLayers):
        if bottomDepth[i] > waterDepth and bottomDepth[i-1] > waterDepth:
            vertStress[i] = vertStress[i - 1] + ((gamma[i - 1] - gamma_w) * thick[i - 1] * (1 - refDepth[i - 1])) + ((gamma[i] - gamma_w) * thick[i] * refDepth[i])
        elif bottomDepth[i] > waterDepth:
            vertStress[i] = vertStress[i - 1] + (gamma[i - 1] * thick[i - 1] * (1 - refDepth[i - 1])) + ((gamma[i] - gamma_w) * thick[i] * refDepth[i])
        else:
            vertStress[i] = vertStress[i - 1] + (gamma[i - 1] * thick[i - 1] * (1 - refDepth[i - 1])) + (gamma[i] * thick[i] * refDepth[i])
        kNot[i] = 1.0 - sin(phi[i] * 2.0 * pi / 360.0)
        refStress[i] = vertStress[i] * (1.0 + 2.0 * kNot[i]) / 3.0

    # 4) Reference low-strain shear modulus, Gref (kN/m2)
    # For top layer function of Vs_not and refPress to match pressure independent (PI) travel time
    # For all other layers, directly computed from Vs (of constant profile)
    for i in range(1, numLayers):
        if i == 0:
            Gref[i] = (gamma[i] / g) * (Vs_not ** 2.0) * (gamma[0] * (1.0 + 2.0 * kNot[0]) / (3.0 * refStress[0])) ** (-pressDependCoe[0])
        else:
            Gref[i] = (gamma[i] / g) * Vs[i] ** 2.0

    # ------------------------------------------------------------------------------------------------------------------
    # Compute the pressure dependent shear modulus
    # ------------------------------------------------------------------------------------------------------------------

    # Define the depths to consider for evaluating and plotting
    totalDepth = np.sum(thick)  # Total depth of the deposit (m)
    depthIncr = 0.001  # Increment of depth (m)
    depths = np.arange(0, totalDepth + depthIncr, depthIncr) # Return evenly spaced values within the interval 'depthIcr' (m)
                                                             # It doesn't include the value 'totalDepth + depthIncr'

    # Effective vertical stress at the bottom of each layer, bottomVertStress (kN/m2)
    bottomVertStress = np.zeros(len(bottomDepth))  # Define the size of the array
    # Top layer
    if thick[0] > waterDepth:  # Just in the case where waterDepth = 0
        bottomVertStress[0] = (gamma[0] - gamma_w) * thick[0]
    else:
        bottomVertStress[0] = gamma[0] * thick[0]
    # All other layers
    for i in range(1, numLayers):
        if bottomDepth[i] > waterDepth:
            bottomVertStress[i] = bottomVertStress[i - 1] + (gamma[i] - gamma_w) * thick[i]
        else:
            bottomVertStress[i] = bottomVertStress[i - 1] + gamma[i] * thick[i]

    # Generate two arrays to be used in subsequent computations
    # np.hstack: stack arrays in sequence
    interfaceDepth = np.hstack(([0], bottomDepth)) # Array with depths of interfaces between layers, including 0 (m)
    interfaceVertStress = np.hstack(([0], bottomVertStress)) # Array with vertical effective stresses at interfaces between layers, including 0 (kPa) at surface (kN/m2 = kPa)

    # Generate an array containing the vertical effective stresses evaluated in 'depths' (kN/m2)
    # np.interp: one-dimensional linear interpolation for monotonically increasing sample points
    # np.iterp(x,xp,fp)
    # x: The x-coordinates at which to evaluate the interpolated values
    # xp: The x-coordinates of the data points. Must be increasing if argument period is not specified
    # fp: The y-coordinates of the data points. Same length as xp
    vertStress = np.interp(depths, interfaceDepth, interfaceVertStress) #lin interp between stresses at layer bottoms

    #### compute pressure dependent shear modulus
    # Define the size of the arrays to be generated
    G_independ = np.zeros(len(depths))           # Reference (pressure-independent) low-strain shear modulus (kN/m2)
    G_depend = np.zeros(len(depths))             # Pressure-dependent low-strain shear modulus (kN/m2)
    Vs_depend  = np.zeros(len(depths))           # Pressure-dependent shear wave velocity profile (m/s)
    VsArray = np.zeros(len(depths))              # Pressure-independent shear wave velocity profile (m/s)
    pressDependCoeArray = np.zeros(len(depths))  # Pressure dependent coefficient
    rhoArray = np.zeros(len(depths))             # Density (Mg/m3)
    refStressArray = np.zeros(len(depths))       # Reference mean effective confining pressure, p' (kN/m2)
    meanStress = np.zeros(len(depths))           # Mean effective confining pressure, p'r (kN/m2)

    limit = np.zeros(len(depths))

    # Fill some of the arrays with the information already known
    # Loop over layers
    for i in range(0, numLayers):
        # Loop over 'depths'
        for j in range(0, len(depths)):
            # Check if the depth depths[i] is in layer i
            if depths[j] <= bottomDepth[i] and G_independ[j] == 0.:
                G_independ[j] = Gref[i]
                pressDependCoeArray[j] = pressDependCoe[i]
                refStressArray[j] = refStress[i]
                rhoArray[j] = gamma[i] / g
                VsArray[j] = Vs[i]
                meanStress[j] = vertStress[j] * (1.0 + 2.0 *kNot[i]) / 3.0
                limit[j] = 20

    # Compute the pressure-dependent low-strain shear modulus and shear wave velocity for each depth in 'depths'
    for i in range(0, len(depths)):
        G_depend[i] = G_independ[i] * (meanStress[i] / refStressArray[i]) ** pressDependCoeArray[i]
        Vs_depend[i] = np.sqrt(G_depend[i] / rhoArray[i])

    # Compute the travel time of shear wave through both pressure-independent and -dependent profiles

    # Travel time (s) of the pressure-independent Vs profile
    travelTimeLayerIndep = np.zeros(numLayers)
    for i in range(0, numLayers):
        travelTimeLayerIndep[i] = thick[i]/Vs[i]
    travelTimeIndep = np.sum(travelTimeLayerIndep)

    # Travel time (s) of the pressure-dependent Vs profile
    travelTimeLayerDepend = np.zeros(numLayers)
    # integrate.cumtrapz(y,x): Cumulatively integrate y(x) using the composite trapezoidal rule
    travelTimeDependIncr = integrate.cumtrapz(1 / Vs_depend[1:], depths[1:])
    # integrate.trapz(y,x): Integrate y(x) using the composite trapezoidal rule.
    travelTimeDependTotal = integrate.trapz(1 / Vs_depend[1:], depths[1:])
    # Travel time (s) of the pressure-dependent Vs profile, by layer
    for i in range(0, numLayers):
        if i == 0:
            travelTimeLayerDepend[i] = travelTimeDependIncr[int(bottomDepth[i] / depthIncr - 2.0)] # Need to subtract 1 b/c leaving out first increment from 0 to depthIncr...
                                                                                              # and subtract 1 b/c 1st index is 0 index
        else:
            travelTimeLayerDepend[i] = travelTimeDependIncr[int(bottomDepth[i] / depthIncr - 2.0)] - sum(travelTimeLayerDepend[0:i])

    # ------------------------------------------------------------------------------------------------------------------
    # PLOTTING
    # ------------------------------------------------------------------------------------------------------------------

    plt.figure(figsize=(12, 14))

    plt.plot(VsArray[1:], depths[1:], 'k-')
    plt.plot(Vs_depend[1:], depths[1:], 'b-')
    plt.legend(['Pressure Independent (PI)', 'Pressure Dependent (PD)'], prop={'size':18}, fontsize=30, loc='lower left')
    plt.xlabel('Vs (m/s)', size=18)
    plt.ylabel('Depth (m)', size=18)
    #plt.title('%s' %siteID, y = 1.07, size=26)
    plt.xlim([0, None])
    plt.ylim([totalDepth, 0])

    # Plot the pressure dependent coefficients
    plt.text(Vs[0] - 50.0, midDepth[0], 'd=%.2f' %pressDependCoe[0], ha='right', va='top', fontsize=14)
    for i in range(1, numLayers):
        plt.text(Vs[i] - 50.0, midDepth[i], 'd=%.2f' %pressDependCoe[i], ha='right', va='center', fontsize=14)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    plt.tight_layout()

    plt.savefig(os.path.join(PathOutputs, '%s_PD_Vs.pdf' %siteID), format='pdf')

    print('Travel time PI Vs:', travelTimeIndep)
    print('Travel time PD Vs:', travelTimeDependTotal)

    return depths, VsArray, Vs_depend

# =====================================================================================================================
# FUNCTION getInSituStresses
# =====================================================================================================================

# Compute in situ stresses
def getInSituStresses(GWL, thickness, gamma):
    """
    Mean effective confining pressure (atm)

    Input parameters
    ----------
    GWL: Freatic level (m)
    thickness: Thickness vector (m)
    gamma: Specific weight vector (kN/m3)

    Output
    ----------
    sigma_v0: In-situ vertical total stress (kPa)
    sigma_v0_eff: In-situ vertical effective stress (kPa)
    sigma_mean_eff: In-situ mean effective stress (kPa)
    u0: In-situ equilibrium pore pressure (kPa)

    References
    ----------
    """

    # Unit weight of water (kN/m3)
    gamma_w = 9.80665

    # midDepth vector - Depth to the middle of each layer
    midDepth = np.zeros(len(thickness))
    midDepth[0] = thickness[0] / 2
    cumThickness = thickness[0]
    for i in range(1, len(midDepth)):
        midDepth[i] = cumThickness + thickness[i] / 2
        cumThickness = cumThickness + thickness[i]

    # In-situ equilibrium pore pressure (kPa)
    u0 = np.zeros(len(thickness))
    mask = midDepth < GWL
    u0[mask] = 0.0
    u0[~mask] = gamma_w * (midDepth[~mask] - GWL)

    # In-situ vertical total stress (kPa)
    sigma_v0 = np.zeros(len(thickness))
    sigma_v0[0] = gamma[0] * midDepth[0]
    for i in range(1, len(midDepth)):
        sigma_v0[i] = sigma_v0[i - 1] + gamma[i - 1] * (thickness[i - 1] / 2) + gamma[i] * (thickness[i] / 2)

    # In-situ vertical effective stress (kPa)
    sigma_v0_eff = sigma_v0 - u0

    return sigma_v0, sigma_v0_eff, u0

# =====================================================================================================================
# FUNCTION plotFormat
# =====================================================================================================================

### Process profiles to plot them
def plotFormat(thickness, value):

    topDepth = np.zeros(len(thickness) + 1)
    for i in range(1, len(topDepth)):
        topDepth[i] = topDepth[i - 1] + thickness[i - 1]

    # For plotting
    plotDepth = np.repeat(topDepth, 2)
    plotDepth = np.delete(plotDepth, [0, len(plotDepth)-1])
    plotValue = np.array(value)
    plotValue = np.repeat(plotValue, 2)

    return plotDepth, plotValue

### Process profiles to plot them and use the function SRImethod
def processProfile(thickness, value, unit):

    topDepth = np.zeros(len(thickness) + 1)
    for i in range(1, len(topDepth)):
        topDepth[i] = topDepth[i - 1] + thickness[i - 1]

    # For plotting
    plotDepth = np.repeat(topDepth, 2)
    plotDepth = np.delete(plotDepth, [0, len(plotDepth)-1])
    plotValue = np.array(value)
    plotValue = np.repeat(plotValue, 2)

    # For computing Site Amp
    if unit == 'km':
        # Finer resolution in the first 10 m
        interpDepth1 = np.arange(0, 0.01, 0.000005)
        # Coarser resolution for deeper layers
        interpDepth2 = np.arange(0.01, np.max(topDepth), 0.0005)
    elif unit == 'm':
        # Finer resolution in the first 10 m
        interpDepth1 = np.arange(0, 0.01 * 1000, 0.000005 * 1000)
        # Coarser resolution for deeper layers
        interpDepth2 = np.arange(0.01 * 1000, np.max(topDepth), 0.0005 * 1000)
    interpDepth = np.concatenate((interpDepth1, interpDepth2))
    interpValue = np.interp(interpDepth, plotDepth, plotValue)

    return plotDepth, plotValue, interpDepth, interpValue

# =====================================================================================================================
# FUNCTION getTmodel
# =====================================================================================================================

### Compute T_model
def getTmodel(layerThickness, Vs):

    # Compute average Vs
    numerator = 0
    denominator = 0
    for i in range(len(Vs)):
        numerator = numerator + layerThickness[i]
        denominator = denominator + layerThickness[i] / Vs[i]
    aveVs = numerator / denominator

    Tmodel = 4 * np.sum(layerThickness) / aveVs

    return Tmodel