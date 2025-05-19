"""
Example 1

Purpose: Generate site-response models based on Vs and CPT data.
Written by Felipe Kuncar (kuncarg@gmail.com)
"""

# =====================================================================================================================
# IMPORT LIBRARIES, FUNCTIONS, CLASSES
# =====================================================================================================================

import sys
from pathlib import Path

# Add the root folder of your repo to sys.path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd
import os.path
from pathlib import Path
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import time

from geosite.class_CPTu import CPTu
from geosite.class_Vs import Vs
from geosite.class_PDMY02 import PDMY02
from geosite.genOpenSeesModel import create_modelOpenSees
from geosite.functions import readVsData, trialModelsFromCPTu, PressDepVs, getTmodel, getInSituStresses, plotFormat
from geosite.MRDcurves import getDminDarendeli, Darendeli_MRDcurves, Yee_MRcurve, Darendeli_MRDcurves_forPySeismoSoil

start = time.time()

# =====================================================================================================================
# MODEL INFORMATION
# =====================================================================================================================

# Site ID
site = 'PRPC'

# Model Version
version = 'v2'

# Water table (positive value in m)
# (always at the interface of two model layers)
GWL = 2.20

# Damping model
# D1: Lab-based Dmin (Darendeli, 2001)
# D2: 3 x Lab-based Dmin (Darendeli, 2001)
# D3: Vs-based Dmin (Cambell, 2009)
dampingModel = 'D2'

# Model option for OpenSees
# modelOption = 1: Default MR curves
# modelOption = 2: User-defined MR curves
modelOption = 2

# =====================================================================================================================
# PATHS
# =====================================================================================================================

# Root directory: the outer directory called 'geosite'
rootDir = Path(__file__).parents[2]
print(rootDir)

# Example directory
exampleDir = Path(rootDir / 'examples' / 'example_1')

# Define the path to the directory containing the OpenSees templates
templateDir = Path(exampleDir / 'OpenSeesTemplates')

# Define the path to the directory that contains the site characterization data to read
inputPath = Path(exampleDir / 'Input')
# Define the path to the directory that will save the results
outputPath = Path(exampleDir / 'Output')
if not os.path.exists(outputPath):
    os.makedirs(outputPath)

# =====================================================================================================================
# INITIALIZE LISTS THAT WILL SAVE DATA
# =====================================================================================================================

CPTu_Path = []
CPTu_Data = []
SCPT_Path = []
SCPT_Data = []
SW_Path = []
SW_Data = []

# =====================================================================================================================
# FILES TO BE READ
# =====================================================================================================================

# CPTu
# Format (Excel): depth (m) | qc (kPa) | fs (kPa) | u2 (kPa)
CPTu_Path += [os.path.join(inputPath, '%s_CPTu.xlsx' % site)]

'''
# SCPT
# Format (Excel): depth (m) | Vs (m/s)
SCPT_Path += [os.path.join(inputPath, '%s_SCPT_Vs.xlsx' % site)]
'''

# SW testing
# Format (csv file): data_point_number | depth (m) | Vs (m/s)
SW_Path += [os.path.join(inputPath, '%s_SW_Vs.xlsx' % site)]

# =====================================================================================================================
# READ DATA AND DEFINE OBJECTS
# =====================================================================================================================

### Define CPTu objects

# CPTu 1
# Read data
CPTu_read = pd.read_excel(CPTu_Path[0])
# Assign class CPTu
CPTu_depth = CPTu_read.iloc[:, 0]
qc = CPTu_read.iloc[:, 1]
fs = CPTu_read.iloc[:, 2]
u2 = CPTu_read.iloc[:, 3]
CPTu_Data.append(CPTu(CPTu_depth, qc, fs, u2, GWL))

### Define Vs objects

'''
# From SCPT

# SCPT 1
# Read data
Vs_read = pd.read_excel(SCPT_Path[0])
# Assign class Vs
Vs_depth = Vs_read.iloc[:, 0]
Vs_values = Vs_read.iloc[:, 1]
SCPT_Data.append(Vs(Vs_depth, Vs_values, GWL))
'''

# From SW testing

# SW 1
# Read data
format = 1
LR = 0.0
Vs_depth, Vs_values = readVsData(format, SW_Path[0])
# Assign class Vs
SW_Data.append(Vs(Vs_depth, Vs_values, GWL))

# =====================================================================================================================
# ADJUST SITE-RESPONSE MODEL PARAMETERS
# =====================================================================================================================

# Read model parameters from Excel file and adjust them iteratively
siteResponseModel = pd.read_excel(Path(exampleDir / f'layers{site}.xlsx'))

# Arrays defining the site-response model parameters
# Thickness of model layers (including elastic half-space layer)
layerThickness_siteModel = np.array(siteResponseModel['layerThickness'])
# Constitutive model flag: 0 for Linear | 1 for Nonlinear - Sand-like | 2 for Nonlinear - Clay-like
constModelFlag_siteModel = np.array(siteResponseModel['constModelFlag'])
# Shear-wave velocity (m/s)
Vs_siteModel = np.array(siteResponseModel['Vs'])
# Specific weight (kN/m3)
gamma_siteModel = np.array(siteResponseModel['gamma'])
# Plasticity index
PI_siteModel = np.array(siteResponseModel['PI'])
# Relative density (%)
Dr_siteModel = np.array(siteResponseModel['Dr'])
# Peak friction angle (°) - Simple shear (for DEEPSOIL and PySeismoSoil)
phi_DSS_siteModel = np.array(siteResponseModel['phi_DSS'])
# Peak friction angle (°) - In octahedral space (for OpenSees) - Khosravifar et al. (2012, 2018)
phi_oct_siteModel = (180/np.pi) * np.arcsin((3 * np.tan(phi_DSS_siteModel * np.pi/180)) / (2 * np.sqrt(3) + np.tan(phi_DSS_siteModel * np.pi/180)))
# Undrained shear strength (kPa) - Simple shear (for DEEPSOIL and PySeismoSoil)
Su_DSS_siteModel = np.array(siteResponseModel['Su_DSS'])
# Cohesion (kPa) - In octahedral space (for OpenSees) - Khosravifar et al. (2012)
# Multiply by (np.sqrt(6) / 3) to transform Su_DSS to Su_oct
# Multiply by (3 / np.sqrt(8)) to transform Su_oct to cohesion_oct
cohesion_oct_siteModel = (np.sqrt(6) / 3) * Su_DSS_siteModel * (3 / np.sqrt(8))
# Coefficient of earth pressure at rest
K0_siteModel = 1 - np.sin(phi_DSS_siteModel * np.pi/180)
# Overconsolidation ratio, considered equal to 1 for all soils
OCR_siteModel = np.ones(len(layerThickness_siteModel))
# Layer numbers
layerNumber = np.arange(1, len(layerThickness_siteModel) + 1, 1)
# Elastic half-space
layerNumber[-1] = 0

# For OpenSees:
# Depth-dependent parameters
pressDependCoe_siteModel = np.array(siteResponseModel['pressDependCoe'])
refDepth_siteModel = np.array(siteResponseModel['refDepth'])
# Flags for water table and Vs inversion
# set VsInvTopLayer to "Yes" if there is a velocity inversion immediately below upper layer (else "No")
# set waterTopLayer to "Yes" if water table within upper most layer and layer was split in two (else "No")
# if waterTopLayer == "Yes", should set refDepth(numLayers) = 1.0 and refDepth(numLayers-1) = 0.0
VsInvTopLayer = 'No'
waterTopLayer = 'No'

# For the Yee et al. (2013) procedure
# User-defined transitional shear strain (%)
gamma1_siteModel = np.array(siteResponseModel['gamma1_Yee'])

# =====================================================================================================================
# CREATE OUTPUT EXCEL FILES
# =====================================================================================================================

### CPTu outputs
Ic, Qtn = CPTu_Data[0].Ic(version=1)
I_SBT, none = CPTu_Data[0].Ic(version=2)

# Basic CPT results
CPTu_basicResults = np.array([CPTu_Data[0].depth, (CPTu_Data[0].qc / 1000), CPTu_Data[0].fs, CPTu_Data[0].u2, CPTu_Data[0].u0, (CPTu_Data[0].qt / 1000), CPTu_Data[0].Rf, CPTu_Data[0].gamma, CPTu_Data[0].sigma_v0, I_SBT, Qtn, CPTu_Data[0].Fr, Ic])
CPTu_basicResults = pd.DataFrame(np.transpose(CPTu_basicResults), columns=['Depth (m)', 'qc (Mpa)', 'fs (kPa)', 'u2 (kPa)', 'u0 (kPa)', 'qt (MPa)', 'Rf (%)', 'Unit weight (kN/m3)', 'sigma_v0 [kPa]', 'I SBT', 'Qtn', 'Fr (%)', 'Ic'])
CPTu_basicResults.to_excel(os.path.join(outputPath, '%s_CPTu_basicResults.xlsx' % site))

'''
# Estimated parameters from CPTu_based correlations
CPTu_estimatedParameters = np.array([CPTu_Data[0].depth, CPTu_Data[0].Vs(Ic, correlation=1), CPTu_Data[0].gamma, CPTu_Data[0].relative_density(Ic, Qtn, correlation=1), CPTu_Data[0].friction_angle(Ic, Qtn, correlation=1), CPTu_Data[0].undrainedShearStrength(Ic)])
CPTu_estimatedParameters = pd.DataFrame(np.transpose(CPTu_estimatedParameters), columns=['Depth (m)', 'Vs (m/s)', 'Unit weight (kN/m3)', 'Dr (%)', 'Friction angle (°)', 'Undrained shear strength (kPa)'])
CPTu_estimatedParameters.to_excel(os.path.join(outputPath, '%s_CPTu_estimatedParameters.xlsx' % site))
'''

# Suggested or trial models obtained from CPTu-based correlations
layerNumber_trial, layerThickness_trial, topDepth, bottomDepth, Ic_trial, Qtn_trial, Vs_trial, gamma_trial, Dr_trial, phi_trial, Su_trial, PI_trial = trialModelsFromCPTu(layerThickness_siteModel, CPTu_Data[0])
CPTu_trialModels = np.array([layerNumber_trial, layerThickness_trial, topDepth, bottomDepth, Ic_trial, Qtn_trial, Vs_trial , gamma_trial, Dr_trial, phi_trial, Su_trial, PI_trial])
CPTu_trialModels = pd.DataFrame(np.transpose(CPTu_trialModels), columns=['Layer', 'Thickness (m)', 'Top Depth (m)', 'Bottom Depth (m)', 'Ic', 'Qtn', 'Vs (m/s)', 'Unit weight (kN/m3)', 'Dr (%)', 'Friction angle (°)', 'Su (kPa)', 'PI'])
CPTu_trialModels.to_excel(os.path.join(outputPath, '%s_CPTu_trialModels.xlsx' % site))

# Suggested or trial models obtained from Vs-based correlations
N1_60 = SW_Data[0].N1_60
Vs_trialModels = np.array([SW_Data[0].depth, SW_Data[0].values, 0.5 * SW_Data[0].gamma_corr(correlation=1) + 0.5 * SW_Data[0].gamma_corr(correlation=2), SW_Data[0].relative_density(N1_60, correlation=1), SW_Data[0].friction_angle(N1_60, correlation=1)])
Vs_trialModels = pd.DataFrame(np.transpose(Vs_trialModels), columns=['Depth (m)', 'Vs (m/s)', 'Unit weight (kN/m3)', 'Dr (%)', 'Friction angle (°)'])
Vs_trialModels.to_excel(os.path.join(outputPath, '%s_Vs_trialModels.xlsx' % site))

'''
# Suggested or trial models obtained from Dr-based correlations
PDMY02_Model = PDMY02(layerThickness_siteModel, Dr_siteModel)
Dr_trialModels = np.array([layerThickness_siteModel, PDMY02_Model.phiAng()])
Dr_trialModels = pd.DataFrame(np.transpose(Dr_trialModels), columns=['Thickness (m)', 'Friction angle (°)'])
Dr_trialModels.to_excel(os.path.join(outputPath, '%s_Dr_trialModels.xlsx' % site))
'''

# =====================================================================================================================
# COMPUTE MINIMUM DAMPING
# =====================================================================================================================

### Dmin (dimensionless) - D1 Formulation: Lab-based Dmin (Darendeli, 2001)
sigma_v0_siteModel, sigma_v0_eff_siteModel, u0_siteModel = getInSituStresses(GWL, layerThickness_siteModel, gamma_siteModel)
D1_siteModel = getDminDarendeli(sigma_v0_eff_siteModel, K0_siteModel, PI_siteModel, OCR_siteModel, freq=1)
print(PI_siteModel)
print(K0_siteModel)
print(D1_siteModel)

### Dmin (dimensionless) - D2 Formulation: 3 x Lab-based Dmin (Darendeli, 2001)
D2_siteModel = 3 * D1_siteModel

### Dim (dimensionless) - D3 Formulation: Vs-based Dmin
# Quality factor, Q (based on Cambell, 2009)
Q_siteModel = 7.17 + 0.0276 * Vs_siteModel
D3_siteModel = 1 / (2 * Q_siteModel)

# Dmin
if dampingModel == 'D1':
      Dmin_siteModel = D1_siteModel
elif dampingModel == 'D2':
      Dmin_siteModel = D2_siteModel
elif dampingModel == 'D3':
      Dmin_siteModel = D3_siteModel

soilType = []
shearStrength_DS = np.zeros(len(Dmin_siteModel))
for i in range(len(Dmin_siteModel)):
      if constModelFlag_siteModel[i] == 1:
            soilType.append('Sand-like soil')
            shearStrength_DS[i] = sigma_v0_eff_siteModel[i] * np.tan(phi_DSS_siteModel[i] * np.pi/180)
      elif constModelFlag_siteModel[i] == 2:
            soilType.append('Clay-like soil')
            # shearStrength[i] = 0.22 * sigma_v0_eff[i]
            shearStrength_DS[i] = Su_DSS_siteModel[i]

# Save results in Excel file
profile = np.array([layerNumber, layerThickness_siteModel, sigma_v0_siteModel, u0_siteModel, sigma_v0_eff_siteModel, Vs_siteModel, D1_siteModel * 100, D2_siteModel * 100, D3_siteModel * 100])
profile = pd.DataFrame(np.transpose(profile), columns=['Layer Number', 'Thickness (m)', 'sigma_v0 (kN/m2)', 'u_0 (kN/m2)', 'sigma_v0_eff (kN/m2)', 'Vs (m/s)', 'Dmin1 (%)', 'Dmin2 (%)', 'Dmin2 (%)'])
profile.to_excel(os.path.join(outputPath, '%s_Dmin_profile.xlsx' % site))

# =====================================================================================================================
# COMPUTE MODULUS REDUCTION CURVES
# =====================================================================================================================

# Density in kg/m3
rho_siteModel = gamma_siteModel * 1000 / 9.81
# G_max (MPa)
Gmax_siteModel = rho_siteModel * (Vs_siteModel) ** 2 / 1000000
# gamma1 parameter for Yee et al. (2013) procedure
gamma1 = 0.1

curve = []
for i, layer in enumerate(layerThickness_siteModel):
    # Compute Darendeli (2001) MR curve
    #shearStrain, G_ratio_ref, a, gamma_r = Darendeli_MRcurve(sigma_v0_eff_siteModel[i], K0_siteModel[i], PI_siteModel[i], OCR_siteModel[i], freq=1, N=10)
    shearStrain, G_ratio_ref, _, a, gamma_r, _ = Darendeli_MRDcurves(sigma_v0_eff_siteModel[i], K0_siteModel[i], PI_siteModel[i], OCR_siteModel[i], freq=1, N=10)
    # Apply Yee et al. (2013) procedure to adjust curve at large strains
    G_ratio = Yee_MRcurve(a, gamma_r, shearStrain, G_ratio_ref, gamma1_siteModel[i], Gmax_siteModel[i], shearStrength_DS[i])
    # Save curve in OpenSees format
    valuesOpenSees = np.zeros(len(shearStrain))
    k = 0
    values = []
    for j in range(len(shearStrain)):
        k = k + 1
        values.append(int(k))
        # gamma (unitless)
        values.append(shearStrain[j] / 100)
        k = k + 1
        values.append(int(k))
        # G/Gmax
        values.append(G_ratio[j])
    curve.append(values)
MRcurves = np.vstack(curve)

# =====================================================================================================================
# CREATE INPUT FILE FOR OPENSEES
# =====================================================================================================================

if not os.path.exists(Path(outputPath / 'inputForOpenSees')):
    os.makedirs(Path(outputPath / 'inputForOpenSees'))

numerator = 0
denominator = 0
for i in range(len(Dmin_siteModel)):
      numerator = numerator + layerThickness_siteModel[i]
      denominator = denominator + layerThickness_siteModel[i] / Dmin_siteModel[i]
damp_OS = numerator / denominator
print('Average value of Dmin is', damp_OS)

# Site period (s)
T_siteModel = getTmodel(layerThickness_siteModel[:-1], Vs_siteModel[:-1])
print('Tmodel is', T_siteModel, 's')

# Rayleigh frequencies (Hz)
freq1Rayleigh = 1/T_siteModel
freq2Rayleigh = 5 * freq1Rayleigh
print('Rayleigh frequencies:', freq1Rayleigh, '(Hz) and', freq2Rayleigh, '(Hz)')
omega1_OS = 2 * np.pi * freq1Rayleigh
print(omega1_OS)
omega2_OS = 2 * np.pi * freq2Rayleigh
print(omega2_OS)

### Write tcl files for OpenSees
# 1) Only input, for running analyses in PC, and
# 2) Complete model, for running analyses in workflow implemented in Mahuika/Maui
outputPathOpenSees = os.path.join(outputPath, 'inputForOpenSees')
layerThickness_OS = layerThickness_siteModel.copy()
# Very small value (0.10 m) assigned to layer 0 (so this layer is almost insigificant and a modification of the OpenSees script is not needed)
layerThickness_OS[-1] = 0.1
create_modelOpenSees(modelOption, exampleDir, templateDir, site, outputPathOpenSees, layerThickness_OS, GWL, gamma_siteModel, Vs_siteModel, constModelFlag_siteModel, phi_oct_siteModel, Dr_siteModel, Su_DSS_siteModel, MRcurves, refDepth_siteModel, pressDependCoe_siteModel, K0_siteModel, VsInvTopLayer, waterTopLayer, LR, damp_OS, omega1_OS, omega2_OS)

# =====================================================================================================================
# CREATE INPUT FILE FOR DEEPSOIL
# =====================================================================================================================

if not os.path.exists(Path(exampleDir / 'Output' / 'inputForDEEPSOIL')):
    os.makedirs(Path(exampleDir / 'Output' / 'inputForDEEPSOIL'))

profile = np.array([layerNumber, layerNumber, layerThickness_siteModel, gamma_siteModel, Vs_siteModel, Dmin_siteModel * 100, constModelFlag_siteModel, shearStrength_DS, PI_siteModel, K0_siteModel])
profile = pd.DataFrame(np.transpose(profile), columns=['Layer Number', 'Layer Name', 'Thickness (m)', 'Unit Weight (kN/m3)', 'Shear Wave Velocity (m/s)', 'Dmin (%)', 'Soil Type', 'Shear Strength (kPa)', 'PI', 'K0'])
profile.to_excel(os.path.join(outputPath, 'inputForDEEPSOIL', '%s_DEEPSOIL_profile.xlsx' % site))

# =====================================================================================================================
# CREATE INPUT FILES FOR PYSEISMOSOIL
# =====================================================================================================================

if not os.path.exists(Path(exampleDir / 'Output' / 'inputForPySeismoSoil')):
    os.makedirs(Path(exampleDir / 'Output' / 'inputForPySeismoSoil'))

layerNumber_PSS = np.arange(1, len(Dmin_siteModel) + 1, 1)
layerNumber_PSS[-1] = 0

# Elastic half-space
layerThickness_PSS = layerThickness_siteModel
layerThickness_PSS[-1] = 0

# Compute density (kg/m3)
rho_PSS = gamma_siteModel * 1000 / 9.806

profile_PSS = np.column_stack((layerThickness_PSS, Vs_siteModel, Dmin_siteModel, rho_PSS, layerNumber_PSS))

# Profile
sep = '\t'
precision=('%.2f', '%.2f', '%.4g', '%.5g', '%d')
np.savetxt(os.path.join(outputPath, 'inputForPySeismoSoil', '%s_PSSprofile.txt' % site), profile_PSS, fmt=precision, delimiter=sep)

# Shear strength (pa)
sep = '\t'
precision=('%.2f')
shearStrength_PSS = shearStrength_DS[:-1] * 1000
np.savetxt(os.path.join(outputPath, 'inputForPySeismoSoil', '%s_PSSshearStrength.txt' % site), shearStrength_PSS, fmt=precision)

# MRD curves for small strains (PySeismoSoil adjust them at large strains)
sep = '\t'
for i, layer in enumerate(layerThickness_siteModel[:-1]):
    shearStrain_PSS, G_ratio_PSS, D_PSS, _, _, _ = Darendeli_MRDcurves_forPySeismoSoil(sigma_v0_eff_siteModel[i], K0_siteModel[i], PI_siteModel[i], OCR_siteModel[i], freq=1, N=10)
    G_ratio_PSS = np.around(G_ratio_PSS, decimals=5)
    D_PSS = np.around(D_PSS, decimals=5) * 100
    if i == 0:
        MRDcurves_PSS = np.column_stack((shearStrain_PSS, G_ratio_PSS, shearStrain_PSS, D_PSS))
    else:
        MRDcurves_PSS = np.column_stack((MRDcurves_PSS, shearStrain_PSS, G_ratio_PSS, shearStrain_PSS, D_PSS))
np.savetxt(os.path.join(outputPath, 'inputForPySeismoSoil', '%s_PSScurves.txt' % site), MRDcurves_PSS, delimiter=sep, fmt='%.5f')

# =====================================================================================================================
# MODIFICATIONS OF VARIABLES TO PLOT
# =====================================================================================================================

# Modify the arrays to plot
depth_modelPlot = np.zeros(1 + len(layerThickness_siteModel))
for i in range(1, len(depth_modelPlot)):
    depth_modelPlot[i] = depth_modelPlot[i-1] + layerThickness_siteModel[i-1]
depth_modelPlot = np.repeat(depth_modelPlot, 2)
depth_modelPlot = np.delete(depth_modelPlot, [0, len(depth_modelPlot) - 1])

Ic_trial = np.repeat(Ic_trial, 2)
Ic_trialModel = np.empty(len(depth_modelPlot))
Ic_trialModel.fill(np.NaN)
Ic_trialModel[0:len(Ic_trial)] = Ic_trial
Qtn_trial = np.repeat(Qtn_trial, 2)
Qtn_trialModel = np.empty(len(depth_modelPlot))
Qtn_trialModel.fill(np.NaN)
Qtn_trialModel[0:len(Qtn_trial)] = Qtn_trial

gamma_modelPlot = np.repeat(gamma_siteModel, 2)
Vs_modelPlot = np.repeat(Vs_siteModel, 2)
Dr_modelPlot = np.repeat(Dr_siteModel, 2)
phi_DSS_modelPlot = np.repeat(phi_DSS_siteModel, 2)
phi_oct_modelPlot = np.repeat(phi_oct_siteModel, 2)
Su_DSS_modelPlot = np.repeat(Su_DSS_siteModel, 2)
cohesion_oct_modelPlot = np.repeat(cohesion_oct_siteModel, 2)
PI_modelPlot = np.repeat(PI_siteModel, 2)

### Define the model Vs object
# Assign class Vs
Vs_m = Vs(depth_modelPlot, Vs_modelPlot, GWL)

### Define PDMY02 object
# Assign class PDMY02
PDMY02_Model = PDMY02(depth_modelPlot, Dr_modelPlot)

### Vs
N1_60 = Vs_m.N1_60

# =====================================================================================================================
# CREATE FIGURES
# =====================================================================================================================

# x limit
x_lim_Vs = max(Vs_modelPlot) + 30
# y limit
y_lim_cpt = max(CPTu_Data[0].depth)
y_lim_model = max(depth_modelPlot) + 5

#----------------------------------------------------------------------------------------------------------------------
# PRESSURE-DEPENDENT Vs PROFILE
# ---------------------------------------------------------------------------------------------------------------------

PressDepVs(outputPath, layerThickness_siteModel, Vs_siteModel, gamma_siteModel, phi_oct_siteModel, site,
           pressDependCoe_siteModel, refDepth_siteModel, VsInvTopLayer, GWL)

# ---------------------------------------------------------------------------------------------------------------------
# SOIL PROFILE
# ---------------------------------------------------------------------------------------------------------------------

# Create a figure
fig = plt.figure(figsize=(16, 8))

# Title
#fig.suptitle("Soil Profile - %s" % (SiteID), fontsize=16)

#x_lim_Vs = max(Vs_Site * 1000) + 30
x_lim_Vs = 450
y_lim_soilProfile = 32

# Create the first gridspec for ploting the model biases
gs = fig.add_gridspec(nrows=1, ncols=4, width_ratios=[1, 1, 0.65, 0.65], left=0.038, bottom=0.062, right=0.995, top=0.97, wspace=0.25)

### Shear-wave velocity (m/s)
ax = fig.add_subplot(gs[0, 0])
# Direct measurements:
#ax3.plot(SCPT_Data[0].values, SCPT_Data[0].depth, 'k-', linewidth=1.2, label='SCPT test')
#ax.plot(SW_Data[0].values, SW_Data[0].depth, color='k', linewidth=1.8, label='Surface Wave Testing')
# CPTu-based correlations:
Ic, Qtn = CPTu_Data[0].Ic(version=1)
ax.plot(CPTu_Data[0].Vs(Ic, correlation=1), CPTu_Data[0].depth, 'b', linewidth=1.0, label='McGann et al. (2015) correlation')
ax.plot(CPTu_Data[0].Vs(Ic, correlation=2), CPTu_Data[0].depth, 'c', linewidth=1.0, label='Robertson (2009) correlation')
ax.plot(Vs_modelPlot, depth_modelPlot, color='k', linewidth=1.8, linestyle=(0, (1, 0)), label='Site-specific profile')
ax.legend(loc=3, fontsize=10)
ax.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax.set_title('Shear-wave velocity', size=14)
ax.set_xlabel('$V_S$ (m/s)', size=14)
ax.set_ylabel('Depth, z (m)', size=14)
ax.set_xlim([0, x_lim_Vs])
ax.set_ylim([0, y_lim_soilProfile])
plt.xticks(np.arange(0, x_lim_Vs, 100))
plt.yticks(np.arange(0, y_lim_soilProfile, 2))
plt.gca().invert_yaxis()

### Soil behaviour type index
ax = fig.add_subplot(gs[0, 1])
### Measured
ax.plot(Ic, CPTu_Data[0].depth, color='yellow', linewidth=1.2)
ax.set_title('Soil behavior type index', size=14)
ax.set_xlabel('$I_{c}$', size=14)
ax.set_ylabel('Depth, z (m)', size=14)
ax.set_xlim([1, 4])
ax.set_ylim([0, y_lim_soilProfile])
plt.yticks(np.arange(0, y_lim_soilProfile, 2))
ax.add_patch(patches.Rectangle((0.00, 0.00), 1.31, y_lim_soilProfile, facecolor='whitesmoke'))
ax.add_patch(patches.Rectangle((1.31, 0.00), 0.74, y_lim_soilProfile, facecolor='gainsboro'))
ax.add_patch(patches.Rectangle((2.05, 0.00), 0.55, y_lim_soilProfile, facecolor='silver'))
ax.add_patch(patches.Rectangle((2.60, 0.00), 0.35, y_lim_soilProfile, facecolor='darkgray'))
ax.add_patch(patches.Rectangle((2.95, 0.00), 0.65, y_lim_soilProfile, facecolor='gray'))
ax.add_patch(patches.Rectangle((3.60, 0.00), 0.40, y_lim_soilProfile, facecolor='dimgray'))
ax.text(1.00+(1.31-1.00)/2, y_lim_soilProfile/2, "Gravelly sand - Dense sand", ha='center', va='center', rotation=90, size=14, color="black")
ax.text(1.31+0.74/2, y_lim_soilProfile/2, " Sands: clean sand - silty sand", ha='center', va='center', rotation=90, size=14, color="black")
ax.text(2.05+0.55/2, y_lim_soilProfile/2, "Sand mixtures: silty sand - sandy silt", ha='center', va='center', rotation=90, size=14, color="black")
ax.text(2.60+0.35/2, y_lim_soilProfile/2, "Silt mixtures: clayey silt - silty clay", ha='center', va='center', rotation=90, size=14, color="black")
ax.text(2.95+0.65/2, y_lim_soilProfile/2, "Clays: silty clay - clay", ha='center', va='center', rotation=90, size=14, color="black")
ax.text(3.60+(4.00-3.60)/2, y_lim_soilProfile/2, "Organic soils: peats", ha='center', va='center', rotation=90, size=14, color="black")
ax.plot([1, 4], [GWL, GWL], 'b--', linewidth=1.0)
ax.text(3.0, GWL+1.0, 'Water table', color='blue', size=12)
plt.gca().invert_yaxis()

# Cone resistance (MPa)
ax = fig.add_subplot(gs[0, 2])
ax.plot(CPTu_Data[0].qt / 1000, CPTu_Data[0].depth, 'k-', linewidth=1.0)
ax.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax.set_title('Cone resistance', size=14)
ax.set_xlabel('Tip resistance $q_t$ (MPa)', size=14)
ax.set_ylabel('Depth, z (m)', size=14)
ax.set_xlim([min(CPTu_Data[0].qt) / 1000, 1.1 * max(CPTu_Data[0].qt) / 1000])
ax.set_xlim([0, None])
ax.set_ylim([0, 35])
plt.yticks(np.arange(0, y_lim_soilProfile, 2))
plt.gca().invert_yaxis()

# Friction ratio (%)
ax = fig.add_subplot(gs[0, 3])
ax.plot(CPTu_Data[0].Rf, CPTu_Data[0].depth, 'k-', linewidth=1.0)
ax.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax.set_title('Friction ratio', size=14)
ax.set_xlabel('$R_f$ (%)', size=14)
ax.set_ylabel('Depth, z (m)', size=14)
ax.set_xlim([min(CPTu_Data[0].Rf), max(CPTu_Data[0].Rf)])
ax.set_xlim([0, 6])
ax.set_ylim([0, 35])
plt.yticks(np.arange(0, y_lim_soilProfile, 2))
plt.gca().invert_yaxis()

# Save figure
plt.savefig(os.path.join(outputPath, '%s_soilProfile.pdf' % site))
plt.savefig(os.path.join(outputPath, '%s_soilProfile.png' % site), dpi=600)

# ---------------------------------------------------------------------------------------------------------------------
# SOIL PROFILE AND SITE-RESPONSE MODEL
# ---------------------------------------------------------------------------------------------------------------------

# Create a figure
fig = plt.figure(figsize=(14, 8))

y_lim_soilProfile = 30

# Create the first gridspec for ploting the model biases
gs = fig.add_gridspec(nrows=1, ncols=5, left=0.040, bottom=0.06, right=0.98, top=0.97, wspace=0.20, width_ratios=[1, 1, 0.65, 0.65, 0.65])

### Layering
deleteValues = []
j = 1
for i in range(1, len(bottomDepth) + 1):
    deleteValues.append(j)
    j = 1 + (2 * i)
bottomDepth = np.delete(depth_modelPlot, deleteValues)

### Shear-wave velocity (m/s)
ax = fig.add_subplot(gs[0, 0])
# Direct measurements:
ax.plot(SW_Data[0].values, SW_Data[0].depth, color='k', linewidth=1.5)
for i in range(1, len(bottomDepth)):
    if i == 1:
        ax.text(350, 0.5, 'Layer %s' % i, fontsize=8, color='b')
        ax.plot([0, 1000], [bottomDepth[i], bottomDepth[i]], color='k', linestyle='--', linewidth=0.5)
        ax.text(20, bottomDepth[i] - 0.1, '%.1f m' % bottomDepth[i], fontsize=8, color='gray')
    elif i == len(bottomDepth) - 1:
        ax.text(20,  bottomDepth[i-1] + 0.5, 'Elastic Half-Space', fontsize=8, color='b')
    else:
        ax.text(350, bottomDepth[i-1] + 0.5, 'Layer %s' % i, fontsize=8, color='b')
        ax.plot([0, 1000], [bottomDepth[i], bottomDepth[i]], color='k', linestyle='--', linewidth=0.5)
        ax.text(20, bottomDepth[i] - 0.1, '%.1f m' % bottomDepth[i], fontsize=8, color='gray')
ax.set_title('Shear-wave velocity', size=12)
ax.set_xlabel('$V_S$ (m/s)', size=12)
ax.set_ylabel('Depth, z (m)', size=12)
ax.set_xlim([0, 450])
ax.set_ylim([0, y_lim_soilProfile])
plt.xticks(np.arange(0, 500, 100))
plt.yticks(np.arange(0, y_lim_soilProfile + 2, 2))
plt.gca().invert_yaxis()

### Soil behaviour type index
ax = fig.add_subplot(gs[0, 1])
### Measured
ax.plot(Ic, CPTu_Data[0].depth, color='yellow', linewidth=1.2, label='Actual value')
# Layering:
for i in range(1, len(bottomDepth)):
    if i == len(bottomDepth) - 1:
        continue
    else:
        ax.plot([0, 1000], [bottomDepth[i], bottomDepth[i]], color='k', linestyle='--', linewidth=0.5)
### Site-response model
ax.plot(Ic_trialModel, depth_modelPlot, 'r', linewidth=1.2, linestyle=(0, (5, 3)), label='Average value per layer')
ax.legend(loc=3, fontsize=8)
ax.set_title('Soil behaviour type index', size=12)
ax.set_xlabel('$I_{c}$', size=12)
ax.set_xlim([1, 4])
ax.set_ylim([0, y_lim_soilProfile])
plt.yticks(np.arange(0, y_lim_soilProfile + 2, 2))
ax.add_patch(patches.Rectangle((0.00, 0.00), 1.31, y_lim_soilProfile, facecolor='whitesmoke'))
ax.add_patch(patches.Rectangle((1.31, 0.00), 0.74, y_lim_soilProfile, facecolor='gainsboro'))
ax.add_patch(patches.Rectangle((2.05, 0.00), 0.55, y_lim_soilProfile, facecolor='silver'))
ax.add_patch(patches.Rectangle((2.60, 0.00), 0.35, y_lim_soilProfile, facecolor='darkgray'))
ax.add_patch(patches.Rectangle((2.95, 0.00), 0.65, y_lim_soilProfile, facecolor='gray'))
ax.add_patch(patches.Rectangle((3.60, 0.00), 0.40, y_lim_soilProfile, facecolor='dimgray'))
ax.text(1.00+(1.31-1.00)/2, y_lim_soilProfile/2, "gravelly sand - dense sand", ha='center', va='center', rotation=90, size=10, color="black")
ax.text(1.31+0.74/2, y_lim_soilProfile/2, " sands: clean sand - silty sand", ha='center', va='center', rotation=90, size=10, color="black")
ax.text(2.05+0.55/2, y_lim_soilProfile/2, "sand mixtures: silty sand - sandy silt", ha='center', va='center', rotation=90, size=10, color="black")
ax.text(2.60+0.35/2, y_lim_soilProfile/2, "silt mixtures: clayey silt - silty clay", ha='center', va='center', rotation=90, size=10, color="black")
ax.text(2.95+0.65/2, y_lim_soilProfile/2, "clays: silty clay - clay", ha='center', va='center', rotation=90, size=10, color="black")
ax.text(3.60+(4.00-3.60)/2, y_lim_soilProfile/2, "organic soils: peats", ha='center', va='center', rotation=90, size=10, color="black")
ax.plot([1, 4], [GWL, GWL], 'b--', linewidth=1.0)
ax.text(3, GWL-0.2, 'Water table', color='b', fontsize=9)
plt.gca().invert_yaxis()

# Cone resistance (MPa)
ax = fig.add_subplot(gs[0, 2])
# Measured
ax.plot(CPTu_Data[0].qt / 1000, CPTu_Data[0].depth, 'k-', linewidth=1.0)
# Layering:
for i in range(1, len(bottomDepth)):
    if i == len(bottomDepth) - 1:
        continue
    else:
        ax.plot([0, 1000], [bottomDepth[i], bottomDepth[i]], color='k', linestyle='--', linewidth=0.5)
ax.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax.set_title('Cone resistance', size=12)
ax.set_xlabel('Tip resistance $q_t$ (MPa)', size=12)
ax.set_xlim([min(CPTu_Data[0].qt) / 1000, 1.1 * max(CPTu_Data[0].qt) / 1000])
ax.set_xlim([0, None])
ax.set_ylim([0, y_lim_soilProfile])
plt.yticks(np.arange(0, y_lim_soilProfile + 2, 2))
plt.gca().invert_yaxis()

# Friction ratio (%)
ax = fig.add_subplot(gs[0, 3])
# Measured
ax.plot(CPTu_Data[0].Rf, CPTu_Data[0].depth, 'k-', linewidth=1.0)
# Layering:
for i in range(1, len(bottomDepth)):
    if i == len(bottomDepth) - 1:
        continue
    else:
        ax.plot([0, 1000], [bottomDepth[i], bottomDepth[i]], color='k', linestyle='--', linewidth=0.5)
ax.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax.set_title('Friction ratio', size=12)
ax.set_xlabel('$R_f$ (%)', size=12)
ax.set_xlim([min(CPTu_Data[0].Rf), max(CPTu_Data[0].Rf)])
ax.set_xlim([0, None])
ax.set_ylim([0, y_lim_soilProfile])
plt.yticks(np.arange(0, y_lim_soilProfile + 2, 2))
plt.gca().invert_yaxis()

# Normalized cone penetration resistance, Qtn (dimensionless)
ax = fig.add_subplot(gs[0, 4])
### Measured
ax.plot(Qtn, CPTu_Data[0].depth, color='k', linewidth=1.2, label='Actual value')
### Site-response model
ax.plot(Qtn_trialModel, depth_modelPlot, 'r', linewidth=1.2, linestyle=(0, (5, 3)), label='Average value per layer')
# Layering:
for i in range(1, len(bottomDepth)):
    if i == len(bottomDepth) - 1:
        continue
    else:
        ax.plot([0, 1000], [bottomDepth[i], bottomDepth[i]], color='k', linestyle='--', linewidth=0.5)
ax.legend(loc=3, fontsize=8)
ax.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax.set_title('Normalized cone resistance', size=12)
ax.set_xlabel('$Q_{tn}$', size=12)
ax.set_xlim([0, 400])
ax.set_ylim([0, y_lim_soilProfile])
plt.yticks(np.arange(0, y_lim_soilProfile + 2, 2))
plt.gca().invert_yaxis()

# Save figure
plt.savefig(os.path.join(outputPath, '%s_soilProfile_model.pdf' % site))
plt.savefig(os.path.join(outputPath, '%s_soilProfile_model.png' % site), dpi=600)

# ---------------------------------------------------------------------------------------------------------------------
# BASIC CPTu FIGURE
# ---------------------------------------------------------------------------------------------------------------------

# Create a figure
fig1 = plt.figure(figsize=(15, 8))

# Title
fig1.suptitle("CPTu Information - CPTu1 - Site %s" % site, fontsize=16)

# Create the first gridspec for ploting the model biases
gs1 = fig1.add_gridspec(nrows=1, ncols=5, left=0.05, bottom=0.08, right=0.97, top=0.90, wspace=0.35)

# Cone resistance (MPa)
ax1 = fig1.add_subplot(gs1[0, 0])
ax1.plot(CPTu_Data[0].qc / 1000, CPTu_Data[0].depth, color='0.6', linewidth=1.0, label='$q_c$')
ax1.plot(CPTu_Data[0].qt / 1000, CPTu_Data[0].depth, 'k-', linewidth=1.0, label='$q_t$')
ax1.legend(loc=4)
ax1.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax1.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax1.set_title('Cone resistance', size=14)
ax1.set_xlabel('Tip resistance (MPa)', size=12)
ax1.set_ylabel('Depth (m)', size=12)
ax1.set_xlim([min(CPTu_Data[0].qt) / 1000, max(CPTu_Data[0].qt) / 1000])
ax1.set_xlim([0, None])
ax1.set_ylim([0, y_lim_cpt])
plt.gca().invert_yaxis()

# Sleeve friction (kPa)
ax2 = fig1.add_subplot(gs1[0, 1])
ax2.plot(CPTu_Data[0].fs, CPTu_Data[0].depth, 'k-', linewidth=1.0)
ax2.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax2.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax2.set_title('Sleeve friction', size=14)
ax2.set_xlabel('$f_s$ (kPa)', size=12)
ax2.set_ylabel('Depth (m)', size=12)
ax2.set_xlim([min(CPTu_Data[0].fs), max(CPTu_Data[0].fs)])
ax2.set_xlim([0, None])
ax2.set_ylim([0, y_lim_cpt])
plt.gca().invert_yaxis()

# Friction ratio (%)
ax3 = fig1.add_subplot(gs1[0, 2])
ax3.plot(CPTu_Data[0].Rf, CPTu_Data[0].depth, 'k-', linewidth=1.0)
ax3.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax3.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax3.set_title('Friction ratio', size=14)
ax3.set_xlabel('$R_f$ (%)', size=12)
ax3.set_ylabel('Depth (m)', size=12)
ax3.set_xlim([min(CPTu_Data[0].Rf), max(CPTu_Data[0].Rf)])
ax3.set_xlim([0, 10])
ax3.set_ylim([0, y_lim_cpt])
plt.gca().invert_yaxis()

# Pore pressure (kPa)
ax4 = fig1.add_subplot(gs1[0, 3])
ax4.plot([np.minimum(np.min(CPTu_Data[0].u2), np.min(CPTu_Data[0].u0)), np.maximum(np.max(CPTu_Data[0].u2), np.max(CPTu_Data[0].u0))], [CPTu_Data[0].GWL, CPTu_Data[0].GWL], color='cyan', linestyle=(0, (5, 5)), linewidth=1.5)
ax4.plot(CPTu_Data[0].u0, CPTu_Data[0].depth, 'b-', linewidth=1.0, label='$u_0$')
ax4.plot(CPTu_Data[0].u2, CPTu_Data[0].depth, 'k-', linewidth=1.0, label='$u_2$')
ax4.legend(loc=4)
ax4.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax4.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax4.set_title('Pore pressure (kPa)', size=14)
ax4.set_xlabel('Pressure (kPa)', size=12)
ax4.set_ylabel('Depth (m)', size=12)
ax4.set_xlim([-100, max(CPTu_Data[0].u2)])
ax4.set_ylim([0, y_lim_cpt])
plt.gca().invert_yaxis()

# Soil behaviour type index
ax5 = fig1.add_subplot(gs1[0, 4])
ax5.plot(Ic, CPTu_Data[0].depth, color='yellow', linewidth=1.5)
ax5.set_title('Soil Behaviour Type Index', size=14)
ax5.set_xlabel('$I_{c}$', size=12)
ax5.set_ylabel('Depth (m)', size=12)
ax5.set_xlim([1, 4])
ax5.set_ylim([0, y_lim_cpt])
ax5.add_patch(patches.Rectangle((0.00, 0.00), 1.31, y_lim_cpt, facecolor='whitesmoke'))
ax5.add_patch(patches.Rectangle((1.31, 0.00), 0.74, y_lim_cpt, facecolor='gainsboro'))
ax5.add_patch(patches.Rectangle((2.05, 0.00), 0.55, y_lim_cpt, facecolor='silver'))
ax5.add_patch(patches.Rectangle((2.60, 0.00), 0.35, y_lim_cpt, facecolor='darkgray'))
ax5.add_patch(patches.Rectangle((2.95, 0.00), 0.65, y_lim_cpt, facecolor='gray'))
ax5.add_patch(patches.Rectangle((3.60, 0.00), 0.40, y_lim_cpt, facecolor='dimgray'))
ax5.text(1.00+(1.31-1.00)/2, y_lim_cpt/2, "gravelly sand - dense sand", ha='center', va='center', rotation=90, size=10, color="black")
ax5.text(1.31+0.74/2, y_lim_cpt/2, " sands: clean sand - silty sand", ha='center', va='center', rotation=90, size=10, color="black")
ax5.text(2.05+0.55/2, y_lim_cpt/2, "sand mixtures: silty sand - sandy silt", ha='center', va='center', rotation=90, size=10, color="black")
ax5.text(2.60+0.35/2, y_lim_cpt/2, "silt mixtures: clayey silt - silty clay", ha='center', va='center', rotation=90, size=10, color="black")
ax5.text(2.95+0.65/2, y_lim_cpt/2, "clays: silty clay - clay", ha='center', va='center', rotation=90, size=10, color="black")
ax5.text(3.60+(4.00-3.60)/2, y_lim_cpt/2, "organic soils: peats", ha='center', va='center', rotation=90, size=10, color="black")
plt.gca().invert_yaxis()

# Save figure
plt.savefig(os.path.join(outputPath, '%s_CPTu_basic.pdf' % site))

# ---------------------------------------------------------------------------------------------------------------------
# ESTIMATED PROPERTIES
# ---------------------------------------------------------------------------------------------------------------------

# Create a figure
fig2 = plt.figure(figsize=(22, 8))

# Create the first gridspec for ploting the model biases
gs2 = fig2.add_gridspec(nrows=1, ncols=8, left=0.025, bottom=0.06, right=0.995, top=0.97, wspace=0.20)

# Soil behaviour type index
ax1 = fig2.add_subplot(gs2[0, 0])
### Measured
ax1.plot(Ic, CPTu_Data[0].depth, color='yellow', linewidth=1.2, label='Actual value')
### Site-response model
ax1.plot(Ic_trialModel, depth_modelPlot, 'r', linewidth=1.2, linestyle=(0, (5, 3)), label='Average value per layer')
ax1.legend(loc=3, fontsize=7)
# Layering:
for i in range(1, len(bottomDepth)):
    if i == len(bottomDepth) - 1:
        continue
    else:
        ax1.plot([0, 1000], [bottomDepth[i], bottomDepth[i]], color='k', linestyle='--', linewidth=0.5)
ax1.set_title('Soil behavior type index', size=12)
ax1.set_xlabel('$I_{c}$', size=12)
ax1.set_ylabel('Depth, z (m)', size=12)
ax1.set_xlim([1, 4])
ax1.set_ylim([0, y_lim_model])
plt.yticks(np.arange(0, y_lim_model, 2))
ax1.add_patch(patches.Rectangle((0.00, 0.00), 1.31, y_lim_model, facecolor='whitesmoke'))
ax1.add_patch(patches.Rectangle((1.31, 0.00), 0.74, y_lim_model, facecolor='gainsboro'))
ax1.add_patch(patches.Rectangle((2.05, 0.00), 0.55, y_lim_model, facecolor='silver'))
ax1.add_patch(patches.Rectangle((2.60, 0.00), 0.35, y_lim_model, facecolor='darkgray'))
ax1.add_patch(patches.Rectangle((2.95, 0.00), 0.65, y_lim_model, facecolor='gray'))
ax1.add_patch(patches.Rectangle((3.60, 0.00), 0.40, y_lim_model, facecolor='dimgray'))
ax1.text(1.00+(1.31-1.00)/2, y_lim_model/2, "gravelly sand - dense sand", ha='center', va='center', rotation=90, size=10, color="black")
ax1.text(1.31+0.74/2, y_lim_model/2, " sands: clean sand - silty sand", ha='center', va='center', rotation=90, size=10, color="black")
ax1.text(2.05+0.55/2, y_lim_model/2, "sand mixtures: silty sand - sandy silt", ha='center', va='center', rotation=90, size=10, color="black")
ax1.text(2.60+0.35/2, y_lim_model/2, "silt mixtures: clayey silt - silty clay", ha='center', va='center', rotation=90, size=10, color="black")
ax1.text(2.95+0.65/2, y_lim_model/2, "clays: silty clay - clay", ha='center', va='center', rotation=90, size=10, color="black")
ax1.text(3.60+(4.00-3.60)/2, y_lim_model/2, "organic soils: peats", ha='center', va='center', rotation=90, size=10, color="black")
ax1.plot([1, 4], [GWL, GWL], 'b--', linewidth=1.0)
plt.gca().invert_yaxis()

# Normalized cone penetration resistance, Qtn (dimensionless)
ax2 = fig2.add_subplot(gs2[0, 1])
### Measured
ax2.plot(Qtn, CPTu_Data[0].depth, color='yellow', linewidth=1.2, label='From measurments')
### Site-response model
ax2.plot(Qtn_trialModel, depth_modelPlot, 'r', linewidth=1.2, linestyle=(0, (5, 3)), label='Average value per layer')
ax2.legend(loc=3, fontsize=7)
# Layering
for i in range(1, len(bottomDepth)):
    if i == 1:
        ax2.text(290, 0.5, 'Layer %s' % i, fontsize=8, color='b')
        ax2.plot([0, 1000], [bottomDepth[i], bottomDepth[i]], color='k', linestyle='--', linewidth=0.5)
    elif i == len(bottomDepth) - 1:
        ax2.text(190,  bottomDepth[i-1] + 0.5, 'Elastic half-space', fontsize=8, color='b')
    else:
        ax2.text(290, bottomDepth[i-1] + 0.5, 'Layer %s' % i, fontsize=8, color='b')
        ax2.plot([0, 1000], [bottomDepth[i], bottomDepth[i]], color='k', linestyle='--', linewidth=0.5)
ax2.legend(loc=3, fontsize=7)
#ax2.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
#ax2.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax2.set_title('Norm. cone resistance', size=12)
ax2.set_xlabel('$Q_{tn}$', size=12)
ax2.set_xlim([0, 400])
ax2.set_ylim([0, y_lim_model])
plt.yticks(np.arange(0, y_lim_model, 2))
plt.gca().invert_yaxis()

# Shear-wave velocity (m/s)
ax3 = fig2.add_subplot(gs2[0, 2])
# Direct measurements:
ax3.plot(SW_Data[0].values, SW_Data[0].depth, color='grey', linewidth=1.4, label='SW testing')
# CPTu-based correlations:
ax3.plot(CPTu_Data[0].Vs(Ic, correlation=1), CPTu_Data[0].depth, 'b', linewidth=1.0, label='CPTu - McGann et al. (2015)')
ax3.plot(CPTu_Data[0].Vs(Ic, correlation=2), CPTu_Data[0].depth, 'c', linewidth=1.0, label='CPTu - CPT Guide (2015)')
# Site-response model:
ax3.plot(Vs_modelPlot, depth_modelPlot, 'r', linewidth=1.4, linestyle=(0, (5, 3)), label='Site-response model')
ax3.legend(loc=3, fontsize=7)
# Layering:
for i in range(1, len(bottomDepth)):
    if i == len(bottomDepth) - 1:
        continue
    else:
        ax3.plot([0, 1000], [bottomDepth[i], bottomDepth[i]], color='k', linestyle='--', linewidth=0.5)
#ax3.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
#ax3.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax3.set_title('Shear-wave velocity', size=12)
ax3.set_xlabel('$V_s$ (m/s)', size=12)
ax3.set_xlim([0, x_lim_Vs])
ax3.set_ylim([0, y_lim_model])
plt.xticks(np.arange(0, x_lim_Vs, 100))
plt.yticks(np.arange(0, y_lim_model, 2))
plt.gca().invert_yaxis()

# Soil total unit weight (kPa)
ax4 = fig2.add_subplot(gs2[0, 3])
# CPTu-based correlations:
ax4.plot(CPTu_Data[0].gamma, CPTu_Data[0].depth, 'b-', linewidth=1.0, label='CPTu - Robertson & Cabal (2010)')
# Vs-based correlations:
ax4.plot(Vs_m.gamma_corr(correlation=1), Vs_m.depth, 'green', linewidth=1.1, label='Vs - NCHRP (2019)')
ax4.plot(Vs_m.gamma_corr(correlation=2), Vs_m.depth, color='limegreen', linewidth=1.1, label='Vs - Boore (2015)')
# Site-response model:
ax4.plot(gamma_modelPlot, depth_modelPlot, 'r', linewidth=1.4, linestyle=(0, (5, 3)), label='Site-response model')
# Layering:
for i in range(1, len(bottomDepth)):
    if i == len(bottomDepth) - 1:
        continue
    else:
        ax4.plot([0, 1000], [bottomDepth[i], bottomDepth[i]], color='k', linestyle='--', linewidth=0.5)
ax4.legend(loc=3, fontsize=7)
#ax4.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.1)
#ax4.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.1)
ax4.set_title('Soil total unit weight', size=12)
ax4.set_xlabel('$\gamma$ (kN/m3)', size=12)
ax4.set_xlim([0, 25])
ax4.set_ylim([0, y_lim_model])
plt.yticks(np.arange(0, y_lim_model, 2))
plt.gca().invert_yaxis()

# Relative density (%)
ax5 = fig2.add_subplot(gs2[0, 4])
# CPTu-based correlations:
ax5.plot(CPTu_Data[0].relative_density(Ic, Qtn, correlation=1), CPTu_Data[0].depth, 'b-', linewidth=1.0, label='CPTu - Idriss & Boulanger (2008)')
ax5.plot(CPTu_Data[0].relative_density(Ic, Qtn, correlation=2), CPTu_Data[0].depth, 'c--', linewidth=1.0, label='CPTu - Kulhawy & Mayne (1990)')
# Vs-based correlations:
ax5.plot(Vs_m.relative_density(N1_60, correlation=1), Vs_m.depth, 'g-', linewidth=1.1, label='Vs - Idriss & Boulanger (2008)')
# Site-response model:
ax5.plot(Dr_modelPlot, depth_modelPlot, 'r--', linewidth=1.4, label='Site-response model', linestyle=(0, (5, 3)))
# Layering:
for i in range(1, len(bottomDepth)):
    if i == len(bottomDepth) - 1:
        continue
    else:
        ax5.plot([0, 1000], [bottomDepth[i], bottomDepth[i]], color='k', linestyle='--', linewidth=0.5)
ax5.legend(loc=3, fontsize=7)
#ax5.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.1)
#ax5.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.1)
ax5.set_title('Relative density', size=12)
ax5.set_xlabel('$D_r$ (%)', size=12)
ax5.set_xlim([0, 101])
ax5.set_ylim([0, y_lim_model])
plt.yticks(np.arange(0, y_lim_model, 2))
plt.gca().invert_yaxis()

# Friction angle (°)
ax6 = fig2.add_subplot(gs2[0, 5])
# CPTu-based correlations:
ax6.plot(CPTu_Data[0].friction_angle(Ic, Qtn, correlation=1), CPTu_Data[0].depth, 'b', linewidth=1.0, label='CPTu - Robertson (2010)')
ax6.plot(CPTu_Data[0].friction_angle(Ic, Qtn, correlation=2), CPTu_Data[0].depth, 'darkviolet', linewidth=1.0, label='CPTu - Kulhawy & Mayne (1990)')
ax6.plot(CPTu_Data[0].friction_angle(Ic, Qtn, correlation=3), CPTu_Data[0].depth, 'c', linewidth=1.0, label='CPTu - Robertson & Campanella (1983)')
ax6.plot(PDMY02_Model.phiAng(), PDMY02_Model.depth, 'gray', linewidth=1.4, linestyle=(2, (2, 2)), label='Dr - Karimi & Dashti (2015)')
# Vs-based correlations:
ax6.plot(Vs_m.friction_angle(N1_60, correlation=1), Vs_m.depth, 'g-', linewidth=1.1, label='Vs - NCHRP (2019)')
# Site-response model:
ax6.plot(phi_DSS_modelPlot, depth_modelPlot, 'r', linewidth=1.4, label='Site-response model (simple shear)', linestyle=(0, (5, 3)))
ax6.plot(phi_oct_modelPlot, depth_modelPlot, 'orange', linewidth=1.4, label='Site-response model (octahedral)', linestyle=(0, (2, 1)))
# Layering:
for i in range(1, len(bottomDepth)):
    if i == len(bottomDepth) - 1:
        continue
    else:
        ax6.plot([0, 1000], [bottomDepth[i], bottomDepth[i]], color='k', linestyle='--', linewidth=0.5)
ax6.legend(loc=3, fontsize=7)
#ax6.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
#ax6.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax6.set_title('Friction angle', size=12)
ax6.set_xlabel('$\phi$ (°)', size=12)
ax6.set_xlim([20, 50])
ax6.set_ylim([0, y_lim_model])
plt.yticks(np.arange(0, y_lim_model, 2))
plt.gca().invert_yaxis()

# Undrained shear strength (kPa)
ax7 = fig2.add_subplot(gs2[0, 6])
# CPTu-based correlations:
ax7.plot(CPTu_Data[0].undrainedShearStrength(Ic), CPTu_Data[0].depth, 'b', linewidth=1.0, label='CPTu - Robertson & Cabal (2022)')
# Site-response model:
ax7.plot(Su_DSS_modelPlot, depth_modelPlot, 'r', linewidth=1.4, label='Site-response model (simple shear)', linestyle=(0, (5, 3)))
ax7.plot(cohesion_oct_modelPlot, depth_modelPlot, 'orange', linewidth=1.4, label='Site-response model (octahedral cohesion)', linestyle=(0, (2, 1)))
# Layering:
for i in range(1, len(bottomDepth)):
    if i == len(bottomDepth) - 1:
        continue
    else:
        ax7.plot([0, 1000], [bottomDepth[i], bottomDepth[i]], color='k', linestyle='--', linewidth=0.5)
ax7.legend(loc=3, fontsize=7)
#ax7.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
#ax7.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax7.set_title('Undrained shear strength', size=12)
ax7.set_xlabel('$S_u$ (kPa)', size=12)
ax7.set_xlim([0, 120])
ax7.set_ylim([0, y_lim_model])
plt.yticks(np.arange(0, y_lim_model, 2))
plt.gca().invert_yaxis()

# Plasticity index
ax8 = fig2.add_subplot(gs2[0, 7])
# CPTu-based correlations:
ax8.plot(CPTu_Data[0].plasticityIndex(Ic, correlation=1), CPTu_Data[0].depth, 'b', linewidth=1.0, label='CPTu - Cetin & Ozan (2009)')
ax8.plot(CPTu_Data[0].plasticityIndex(Ic, correlation=2), CPTu_Data[0].depth, 'c', linewidth=1.0, label='Simple rule')
# Site-response model:
ax8.plot(PI_modelPlot, depth_modelPlot, 'r', linewidth=1.4, label='Site-response model', linestyle=(0, (5, 3)))
# Layering:
for i in range(1, len(bottomDepth)):
    if i == len(bottomDepth) - 1:
        continue
    else:
        ax8.plot([0, 1000], [bottomDepth[i], bottomDepth[i]], color='k', linestyle='--', linewidth=0.5)
ax8.legend(loc=3, fontsize=7)
#ax8.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
#ax8.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax8.set_title('Plasticity index', size=12)
ax8.set_xlabel('PI', size=12)
ax8.set_xlim([-1, 65])
ax8.set_ylim([0, y_lim_model])
plt.yticks(np.arange(0, y_lim_model, 2))
plt.gca().invert_yaxis()

# Save figure
plt.savefig(os.path.join(outputPath, '%s_parameterEstimate.pdf' % site))
plt.savefig(os.path.join(outputPath, '%s_parameterEstimate.png' % site), dpi=600)

# ---------------------------------------------------------------------------------------------------------------------
# OTHERS PDMY02 MODEL PARAMETERS
# ---------------------------------------------------------------------------------------------------------------------

# create a figure
fig3 = plt.figure(figsize=(15, 8))

# title
fig3.suptitle("PDMY02 Parameters Estimate - Site %s" % (site), fontsize=16)

# create the first gridspec for ploting the model biases
gs3 = fig3.add_gridspec(nrows=1, ncols=8, left=0.05, bottom=0.08, right=0.97, top=0.90, wspace=0.35)

# PTAng
ax1 = fig3.add_subplot(gs3[0, 0])
# site-response model:
ax1.plot(PDMY02_Model.PTAng(), PDMY02_Model.depth, 'r', linewidth=1.4, linestyle=(0, (5, 3)))
ax1.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax1.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax1.set_title('$PTAng (°)', size=14)
ax1.set_ylabel('Depth (m)', size=12)
ax1.set_xlim([20, 35])
ax1.set_ylim([0, y_lim_model])
plt.gca().invert_yaxis()

# e
ax2 = fig3.add_subplot(gs3[0, 1])
# site-response model:
ax2.plot(PDMY02_Model.e(), PDMY02_Model.depth, 'r', linewidth=1.4, linestyle=(0, (5, 3)))
ax2.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax2.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax2.set_title('$e', size=14)
ax2.set_xlim([0.50, 0.90])
ax2.set_ylim([0, y_lim_model])
plt.gca().invert_yaxis()

# contrac1
ax3 = fig3.add_subplot(gs3[0, 2])
# site-response model:
ax3.plot(PDMY02_Model.contrac1(), PDMY02_Model.depth, 'r', linewidth=1.4, linestyle=(0, (5, 3)))
ax3.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax3.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax3.set_title('$contrac1', size=14)
ax3.set_xlim([0.00, 0.10])
ax3.set_ylim([0, y_lim_model])
plt.gca().invert_yaxis()

# contrac2
ax4 = fig3.add_subplot(gs3[0, 3])
# site-response model:
ax4.plot(PDMY02_Model.contrac2(), PDMY02_Model.depth, 'r', linewidth=1.4, linestyle=(0, (5, 3)))
ax4.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax4.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax4.set_title('$contrac2', size=14)
ax4.set_xlim([0.00, 6.00])
ax4.set_ylim([0, y_lim_model])
plt.gca().invert_yaxis()

# contrac3
ax5 = fig3.add_subplot(gs3[0, 4])
# site-response model:
ax5.plot(PDMY02_Model.contrac3(), PDMY02_Model.depth, 'r', linewidth=1.4, linestyle=(0, (5, 3)))
ax5.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax5.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax5.set_title('$contrac3', size=14)
ax5.set_xlim([0.00, 0.40])
ax5.set_ylim([0, y_lim_model])
plt.gca().invert_yaxis()

# dilat1
ax6 = fig3.add_subplot(gs3[0, 5])
# site-response model:
ax6.plot(PDMY02_Model.dilat1(), PDMY02_Model.depth, 'r', linewidth=1.4, linestyle=(0, (5, 3)))
ax6.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax6.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax6.set_title('$dilat1', size=14)
ax6.set_xlim([0.00, 0.40])
ax6.set_ylim([0, y_lim_model])
plt.gca().invert_yaxis()

# dilat2
ax7 = fig3.add_subplot(gs3[0, 6])
# site-response model:
ax7.plot(PDMY02_Model.dilat2(), PDMY02_Model.depth, 'r', linewidth=1.4, linestyle=(0, (5, 3)))
ax7.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax7.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax7.set_title('$dilat2', size=14)
ax7.set_xlim([0.00, 4.00])
ax7.set_ylim([0, y_lim_model])
plt.gca().invert_yaxis()

# dilat3
ax8 = fig3.add_subplot(gs3[0, 7])
# site-response model:
ax8.plot(PDMY02_Model.dilat3(), PDMY02_Model.depth, 'r', linewidth=1.4, linestyle=(0, (5, 3)))
ax8.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax8.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax8.set_title('$dilat3', size=14)
ax8.set_xlim([-0.01, 0.50])
ax8.set_ylim([0, y_lim_model])
plt.gca().invert_yaxis()

# save figure
plt.savefig(os.path.join(outputPath, '%s_PDMY02_parameterEstimate.pdf' % site))
plt.savefig(os.path.join(outputPath, '%s_PDMY02_parameterEstimate.png' % site), dpi=600)

# ---------------------------------------------------------------------------------------------------------------------
# STRESS AND MIN. DAMPING PROFILES
# ---------------------------------------------------------------------------------------------------------------------

plotDepth, plot_Vs = plotFormat(layerThickness_siteModel, Vs_siteModel)
plotDepth, plot_sigma_v0 = plotFormat(layerThickness_siteModel, sigma_v0_siteModel)
plotDepth, plot_sigma_v0_eff = plotFormat(layerThickness_siteModel, sigma_v0_eff_siteModel)
plotDepth, plot_u0 = plotFormat(layerThickness_siteModel, u0_siteModel)
plotDepth, plot_D1 = plotFormat(layerThickness_siteModel, D1_siteModel)
plotDepth, plot_D2 = plotFormat(layerThickness_siteModel, D2_siteModel)
plotDepth, plot_D3 = plotFormat(layerThickness_siteModel, D3_siteModel)

# Create a figure
fig = plt.figure(figsize=(9, 5))
# Title
#fig.suptitle('Soil Profiles', fontsize=16)
# Create the first gridspec for ploting the model biases
gs = fig.add_gridspec(nrows=1, ncols=3, wspace=0.085)

# Vs Profile
ax = fig.add_subplot(gs[0, 0])
ax.plot(plot_Vs, plotDepth, color='k', linewidth=2.0, linestyle=(0, (1, 0)))
ax.plot([0, 3.5*1000], [100, 100], color='k', linewidth=2.0, linestyle=(4, (2, 2)))
ax.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax.set_xlabel('$V_s$ (m/s)', size=14)
ax.set_ylabel('Depth, z (m)', size=14)
ax.set_xlim([0, 450])
ax.set_ylim([0, 30])
plt.gca().invert_yaxis()

# Stressess
ax = fig.add_subplot(gs[0, 1])
ax.plot(plot_sigma_v0, plotDepth, color='k', linewidth=2.0, linestyle=(0, (1, 0)), label=r'$\sigma_{v0}$')
ax.plot(plot_u0, plotDepth, color='b', linewidth=2.0, linestyle=(2, (2, 2)), label='$u_0$')
ax.plot(plot_sigma_v0_eff, plotDepth, color='k', linewidth=2.0, linestyle=(2, (2, 2)), label=r'$\sigma_{v0,eff}$')
ax.legend(loc=1, fontsize=10)
ax.plot([0, 3.5*1000], [100, 100], color='k', linewidth=2.0, linestyle=(4, (2, 2)))
ax.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax.set_xlabel('Stress (kPa)', size=14)
ax.set_xlim([0, 650])
ax.set_ylim([0, 30])
plt.yticks(fontsize=4, color='white')
plt.gca().invert_yaxis()

# Damping
ax = fig.add_subplot(gs[0, 2])
#ax.text(580, 20, 'Site %s' % site, color='k', fontsize=14)
ax.plot(plot_D1 * 100, plotDepth, color='lightgray', linewidth=2.0, linestyle=(0, (1, 0)), label='Lab-based')
ax.plot(plot_D2 * 100, plotDepth, color='gray', linewidth=2.0, linestyle=(0, (1, 0)), label='3 x Lab-based')
ax.plot(plot_D3 * 100, plotDepth, color='k', linewidth=2.0, linestyle=(0, (1, 0)), label='Vs-Q correlation')
ax.legend(loc=4, fontsize=10)
ax.plot([0, 3.5*1000], [100, 100], color='k', linewidth=2.0, linestyle=(4, (2, 2)))
ax.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax.set_xlabel('$D_{min}$ (%)', size=14)
ax.set_xlim([0, 8])
ax.set_ylim([0, 30])
plt.yticks(fontsize=4, color='white')
plt.gca().invert_yaxis()

fig.subplots_adjust(top=0.985, bottom=0.10, left=0.085, right=0.98)

# save figure
plt.savefig(os.path.join(outputPath, '%s_dampingProfiles.pdf' % site))
plt.savefig(os.path.join(outputPath, '%s_dampingProfiles.png' % site), dpi=600)

end = time.time()
print("Done in {:10.1f} secs".format(end-start))

# ---------------------------------------------------------------------------------------------------------------------
# MR curves
# ---------------------------------------------------------------------------------------------------------------------

# Create a figure
fig = plt.figure(figsize=(11, 11))

# Create the first gridspec for ploting the model biases
gs = fig.add_gridspec(nrows=2, ncols=1, left=0.060, bottom=0.05, right=0.98, top=0.99, hspace=0.1)

from matplotlib.pyplot import cm
color = cm.rainbow(np.linspace(0, 1, len(layerThickness_siteModel)))

# G/Gmax
ax = fig.add_subplot(gs[0, 0])
ax.text(2.2, 0.98, 'Site %s' % site, fontsize=14, color='k')
for i in range(len(MRcurves) - 1):
    MRcurve = MRcurves[i]
    strain = MRcurve[1::4]
    G_ratio = MRcurve[3::4]
    plt.plot(strain * 100, G_ratio, color=color[i], label='Layer %s' % (i+1))
    plt.xscale('log')
    ax.legend(loc=3)
    plt.ylabel(r'$G$/$G_{max}$', fontsize=14)
ax.set_xlim([0.00001, 10])
ax.set_ylim([0, 1.05])

'''
# Alternatively:
# G/Gmax
ax = fig.add_subplot(gs[0, 0])
ax.text(2.2, 0.98, 'Site %s' % site, fontsize=14, color='k')
for i, layer in enumerate(layerThickness_siteModel[:-1]):
    shearStrain, G_ratio_ref, _, a, gamma_r, _ = Darendeli_MRDcurves(sigma_v0_eff_siteModel[i], K0_siteModel[i], PI_siteModel[i], OCR_siteModel[i], freq=1, N=10)
    G_ratio = Yee_MRcurve(a, gamma_r, shearStrain, G_ratio_ref, gamma1, Gmax_siteModel[i], shearStrength_DS[i])
    plt.plot(shearStrain, G_ratio, color=color[i], label='Layer %s' % (i + 1))
    plt.xscale('log')
    ax.legend(loc=3)
    plt.ylabel(r'$G$/$G_{max}$', fontsize=14)
ax.set_xlim([0.00001, 10])
ax.set_ylim([0, 1.05])
'''

# Shear strain, gamma (%)
ax = fig.add_subplot(gs[1, 0])
for i in range(len(MRcurves) - 1):
    MRcurve = MRcurves[i]
    strain = MRcurve[1::4]
    G_ratio = MRcurve[3::4]
    plt.plot(strain * 100, G_ratio * (Gmax_siteModel[i] * 1000) * strain, color=color[i])
    plt.plot([0.00001, 10], [shearStrength_DS[i], shearStrength_DS[i]], linewidth=1.0, linestyle=(0, (5, 10)), color=color[i])
    plt.xscale('log')
    plt.xlabel(r'Shear strain, $\gamma$ (%)', fontsize=14)
    plt.ylabel(r'Shear stress, $\tau$ (kPa)', fontsize=14)
ax.set_xlim([0.00001, 10])
ax.set_xlim([0, None])

'''
# Alternatively:
# Shear strain, gamma (%)
ax = fig.add_subplot(gs[1, 0])
for i, layer in enumerate(layerThickness_siteModel[:-1]):
    shearStrain, G_ratio_ref, _, a, gamma_r, _ = Darendeli_MRDcurves(sigma_v0_eff_siteModel[i], K0_siteModel[i], PI_siteModel[i], OCR_siteModel[i], freq=1, N=10)
    G_ratio = Yee_MRcurve(a, gamma_r, shearStrain, G_ratio_ref, gamma1, Gmax_siteModel[i], shearStrength_DS[i])
    plt.plot(shearStrain, G_ratio * (Gmax_siteModel[i] * 1000) * (shearStrain / 100), color=color[i])
    plt.plot([0.00001, 10], [shearStrength_DS[i], shearStrength_DS[i]], linewidth=1.0, linestyle=(0, (5, 10)), color=color[i])
    plt.xscale('log')
    plt.xlabel(r'Shear strain, $\gamma$ (%)', fontsize=14)
    plt.ylabel(r'Shear stress, $\tau$ (kPa)', fontsize=14)
ax.set_xlim([0.00001, 10])
ax.set_xlim([0, None])
'''

# Save figure
plt.savefig(os.path.join(outputPath, '%s_MRcurves.pdf' % site))
plt.savefig(os.path.join(outputPath, '%s_MRcurves.png' % site), dpi=600)

# ---------------------------------------------------------------------------------------------------------------------
# MRD curves
# ---------------------------------------------------------------------------------------------------------------------

# Create a figure
fig = plt.figure(figsize=(11, 11))

# Create the first gridspec for ploting the model biases
gs = fig.add_gridspec(nrows=2, ncols=1, left=0.060, bottom=0.05, right=0.98, top=0.99, hspace=0.1)

from matplotlib.pyplot import cm
color = cm.rainbow(np.linspace(0, 1, len(layerThickness_siteModel)))

# G/Gmax
ax = fig.add_subplot(gs[0, 0])
ax.text(2.2, 0.98, 'Site %s' % site, fontsize=14, color='k')
for i, layer in enumerate(layerThickness_siteModel[:-1]):
    shearStrain, G_ratio_ref, _, a, gamma_r, _ = Darendeli_MRDcurves(sigma_v0_eff_siteModel[i], K0_siteModel[i], PI_siteModel[i], OCR_siteModel[i], freq=1, N=10)
    print('Yee et al. (2013) procedure for Layer %s:' % (i + 1))
    G_ratio = Yee_MRcurve(a, gamma_r, shearStrain, G_ratio_ref, gamma1_siteModel[i], Gmax_siteModel[i], shearStrength_DS[i])
    plt.plot(shearStrain, G_ratio, color=color[i], label='Layer %s' % (i + 1))
    plt.xscale('log')
    ax.legend(loc=3)
    plt.ylabel(r'$G$/$G_{max}$', fontsize=14)
ax.set_xlim([0.00001, 10])
ax.set_ylim([0, 1.05])

# Damping (%)
ax = fig.add_subplot(gs[1, 0])
for i, layer in enumerate(layerThickness_siteModel[:-1]):
    shearStrain, _, D, _, _, _ = Darendeli_MRDcurves(sigma_v0_eff_siteModel[i], K0_siteModel[i], PI_siteModel[i], OCR_siteModel[i], freq=1, N=10)
    plt.plot(shearStrain, D * 100, color=color[i])
    plt.xscale('log')
    plt.xlabel(r'Shear strain, $\gamma$ (%)', fontsize=14)
    plt.ylabel(r'Damping ratio, D (%)', fontsize=14)
ax.set_xlim([0.00001, 10])
ax.set_xlim([0, None])

# Save figure
plt.savefig(os.path.join(outputPath, '%s_MRDcurves.pdf' % site))
plt.savefig(os.path.join(outputPath, '%s_MRDcurves.png' % site), dpi=600)