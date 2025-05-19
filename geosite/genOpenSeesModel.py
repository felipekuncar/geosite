"""
geosite / genOpenSeesModel.py

Provides functions for automatically creating a OpenSees site-response model.
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
# FUNCTION create_inputOpenSees
# =====================================================================================================================

### To run site-response analysis in personal computer/laptop

def create_inputOpenSees(modelOption, siteID, pathOutput, layerThick, GWL, gamma, Vs, constModelFlag, phi, Dr, cohesion, MRcurves, refDepth, pressCoeff, K0, VsInvTopLayer, waterTopLayer, damp, omega1, omega2):

    # Model option
    # modelOption = 1: Default MR curves
    # modelOption = 2: User-defined MR curves

    GWL = GWL if GWL > 0 else -GWL

    # Reverse order so that idx = 0 is deepest layer (i.e., model base)
    layerThick = layerThick[::-1]
    rho = gamma[::-1] / 9.81
    Vs = Vs[::-1]
    constModelFlag = constModelFlag[::-1]
    phi = phi[::-1]
    Dr = Dr[::-1]
    cohesion = cohesion[::-1]
    refDepth = refDepth[::-1]
    pressCoeff = pressCoeff[::-1]
    K0 = K0[::-1]
    MRcurves = MRcurves[::-1]

    # define PDMY02 object
    PDMY02model = PDMY02(None, Dr)

    # define some variables to be used
    numLayers = len(layerThick) - 1
    soilThick = np.cumsum(layerThick)[-1]

    # open the .tcl file and write
    inputFile = open(os.path.join(pathOutput, '%s_input.tcl' % siteID), 'w')

    inputFile.write('# Define gravity and pi \n')
    inputFile.write('set g 9.80665 \n')
    inputFile.write('set pi [expr atan(1)*4] \n')
    inputFile.write('\n')

    inputFile.write('# ----------------------------------------------------------------------------------------- \n')
    inputFile.write('#  1. DEFINE SOIL GEOMETERY AND MATERIAL PARAMETERS \n')
    inputFile.write('# ----------------------------------------------------------------------------------------- \n')

    inputFile.write('\n')
    inputFile.write('#---SOIL GEOMETRY \n')
    inputFile.write('\n')

    inputFile.write('# Thicknesses of profile (m) (counting layer 0, which is the elastic half-space) \n')
    inputFile.write('set soilThick %s \n' % soilThick)
    inputFile.write('\n')

    inputFile.write('# Number of soil layers (not counting layer 0, which is the elastic half-space) \n')
    inputFile.write('set numLayers %s \n' % numLayers)
    inputFile.write('\n')

    inputFile.write('# Layer thicknesses \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('set layerThick(%s) %s \n' % (ele, layerThick[ele]))
    inputFile.write('\n')

    inputFile.write('# Depth of water table (create a new layer at WT) \n')
    inputFile.write('# If water not present set waterTable anywhere below depth of model \n')
    inputFile.write('set waterTable %s \n' % GWL)
    inputFile.write('\n')

    inputFile.write('#---BASIC MATERIAL PROPERTIES \n')
    inputFile.write('\n')

    inputFile.write('# Shear-wave velocity, Vs (m/s) \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('set Vs(%s) %s \n' % (ele, Vs[ele]))
    inputFile.write('\n')

    inputFile.write('# Reference depth, refDepth (0 is ToL, 1 is BoL, 0.5 is in the middle) \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('set refDepth(%s) %s \n' % (ele, refDepth[ele]))
    inputFile.write('\n')

    inputFile.write('# Pressure dependency coefficient \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('set pressCoeff(%s) %s \n' % (ele, pressCoeff[ele]))
    inputFile.write('\n')

    inputFile.write('# Mass density, rho (Mg/m^3) \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('set rho(%s) %s \n' % (ele, rho[ele]))
    inputFile.write('set rhoWater 1.0 \n')
    inputFile.write('\n')

    inputFile.write('# Poisson ratio of the soil \n')
    inputFile.write('set nu 0.25 \n')
    inputFile.write('\n')

    inputFile.write('# Rock elastic properties) \n')
    inputFile.write('# Bedrock shear wave velocity (m/s) \n')
    inputFile.write('set rockVS $Vs(0) \n')
    inputFile.write('# Bedrock mass density (Mg/m^3) \n')
    inputFile.write('set rockDen $rho(0) \n')
    inputFile.write('\n')

    inputFile.write('# Relative density (%) \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('# Layer %s: Dr = %s \n' % (ele, Dr[ele]))
    inputFile.write('\n')

    inputFile.write('#--- MODEL PARAMETERS \n')
    inputFile.write('\n')

    inputFile.write('# Consitutive model \n')
    inputFile.write('# constModelFlag:  0 for Linear | 1 for PDMY02 (sands) | 2 for PIDMY (clays) \n')
    for ele in reversed(np.arange(1, len(layerThick))):
        inputFile.write('set constModelFlag(%s) %s \n' % (ele, constModelFlag[ele]))
    inputFile.write('\n')

    inputFile.write('# Soil friction angle (°) \n')
    for ele in reversed(np.arange(1, len(layerThick))):
        inputFile.write('set phi(%s) %s \n' % (ele, phi[ele]))
    inputFile.write('\n')

    inputFile.write('# kNot (kPa) \n')
    for ele in reversed(np.arange(1, len(layerThick))):
        inputFile.write('set kNot(%s) %s \n' % (ele, K0[ele]))
    inputFile.write('\n')

    inputFile.write('# Peak shear strain \n')
    inputFile.write('set gammaPeak 0.1 \n')
    inputFile.write('\n')

    inputFile.write('# Phase transformation angle (not for layer 0) (°) \n')
    for ele in reversed(np.arange(1, len(layerThick))):
        if constModelFlag[ele] == 1:
            inputFile.write('set phaseAng(%s) %s \n' % (ele, PDMY02model.PTAng()[ele]))
    inputFile.write('\n')

    inputFile.write('# Contraction coefficients (not for layer 0) \n')
    for ele in reversed(np.arange(1, len(layerThick))):
        if constModelFlag[ele] == 1:
            inputFile.write('set contract1(%s) %s \n' % (ele, PDMY02model.contrac1()[ele]))
    for ele in reversed(np.arange(1, len(layerThick))):
        if constModelFlag[ele] == 1:
            inputFile.write('set contract2(%s) %s \n' % (ele, PDMY02model.contrac2()[ele]))
    for ele in reversed(np.arange(1, len(layerThick))):
        if constModelFlag[ele] == 1:
            inputFile.write('set contract3(%s) %s \n' % (ele, PDMY02model.contrac3()[ele]))
    inputFile.write('\n')

    inputFile.write('# Dilation coefficients (not for layer 0) \n')
    for ele in reversed(np.arange(1, len(layerThick))):
        if constModelFlag[ele] == 1:
            inputFile.write('set dilate1(%s) %s \n' % (ele, PDMY02model.dilat1()[ele]))
    for ele in reversed(np.arange(1, len(layerThick))):
        if constModelFlag[ele] == 1:
            inputFile.write('set dilate2(%s) %s \n' % (ele, PDMY02model.dilat2()[ele]))
    for ele in reversed(np.arange(1, len(layerThick))):
        if constModelFlag[ele] == 1:
            inputFile.write('set dilate3(%s) %s \n' % (ele, PDMY02model.dilat3()[ele]))
    inputFile.write('\n')

    inputFile.write('# Cohesion (not for layer 0) (kPa) \n')
    for ele in reversed(np.arange(1, len(layerThick))):
        if constModelFlag[ele] == 2:
            inputFile.write('set cohesion(%s) %s \n' % (ele, cohesion[ele]))
    inputFile.write('\n')

    inputFile.write('# Void ratio (needed for layer 0 for element definition) \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('set voidR(%s) %s \n' % (ele, PDMY02model.e()[ele]))
    inputFile.write('\n')

    inputFile.write('# Modulus reduction curves \n')
    for ele in reversed(np.arange(1, len(layerThick))):
        curve = MRcurves[ele]
        #inputFile.write('set MR_list_%s %s \n' % (ele, curve))
        inputFile.write('set MR_list_%s ' % (ele))
        inputFile.write('{')
        for i, currCurve in enumerate(curve):
            if i % 2 == 0:
                inputFile.write('%.0f ' % (currCurve))
            else:
                inputFile.write('%.10f ' % (currCurve))
        inputFile.write('}')
        inputFile.write('\n')

    inputFile.write('\n')

    inputFile.write('# Flags for water table and Vs inversion \n')
    inputFile.write('# set VsInvTopLayer to "Yes" if there is a velocity inversion immediately below upper layer (else "No") \n')
    inputFile.write('# set waterTopLayer to "Yes" if water table within upper most layer and layer was split in two (else "No") \n')
    inputFile.write('# if waterTopLayer == "Yes", should set refDepth(numLayers) = 1.0 and refDepth(numLayers-1) = 0.0 \n')
    inputFile.write('set VsInvTopLayer %s \n' % VsInvTopLayer)
    inputFile.write('set waterTopLayer %s \n' % waterTopLayer)
    inputFile.write('\n')

    inputFile.write('# Rayleigh damping parameters) \n')
    inputFile.write('set damp %s \n' % damp)
    inputFile.write('set omega1 %s \n' % omega1)
    inputFile.write('set omega2 %s \n' % omega2)
    inputFile.write('\n')

    return

# =====================================================================================================================
# FUNCTION create_modelOpenSees
# =====================================================================================================================

### To run site-response analysis in HPC Workflow

def create_modelOpenSees(modelOption, rootDir, templateDir, siteID, pathOutput, layerThick, GWL, gamma, Vs, constModelFlag, phi, Dr, cohesion, MRcurves, refDepth, pressCoeff, K0, VsInvTopLayer, waterTopLayer, LR, damp, omega1, omega2):

    # Model option
    # modelOption = 1: Default MR curves
    # modelOption = 2: User-defined MR curves

    # Define reference Vs
    Vs_ref = Vs[-1]

    # Open the .tcl file and write
    initialFile = open(os.path.join(pathOutput, '%s_TotalStress_LR%.1f_%.0f.tcl' % (siteID, LR, Vs_ref)), 'w')

    initialFile.write('################################################################## \n')
    initialFile.write('#                                                                # \n')
    initialFile.write('# Total-Stress Site-Response Analysis                            # \n')
    initialFile.write('# [Pore-Pressure Generation Allowed: No]                         # \n')
    initialFile.write('# Constitutive Models: Linear, PDMY02 (sands), and PIMY (clays)  # \n')
    initialFile.write('# Site: %s                                                     # \n' % siteID)
    initialFile.write('#                                                                # \n')
    initialFile.write('################################################################## \n')
    initialFile.write('\n')

    initialFile.write('# Extract exterior inputs and give them a variable name \n')
    initialFile.write('set site [lindex $argv 0] \n')
    initialFile.write('set modelID [lindex $argv 1] \n')
    initialFile.write('set gMotionPath [lindex $argv 2] \n')
    initialFile.write('set gMotionName [lindex $argv 3] \n')
    initialFile.write('set npts [lindex $argv 4] \n')
    initialFile.write('set dt [lindex $argv 5] \n')
    initialFile.write('set saveDir [lindex $argv 6] \n')
    initialFile.write('\n')

    initialFile.write('set dash "_" \n')
    initialFile.write('set analID $site$dash$modelID$dash$gMotionName\n')
    initialFile.write('set siteModelID $site$dash$modelID\n')
    initialFile.write('\n')

    # Close input file
    initialFile.close()

    # Create input parameters file, which can be used to run/test site-response analysis in personal computer
    create_inputOpenSees(modelOption, siteID, pathOutput, layerThick, GWL, gamma, Vs, constModelFlag, phi, Dr, cohesion, MRcurves, refDepth, pressCoeff, K0, VsInvTopLayer, waterTopLayer, damp, omega1, omega2)

    # Open and copy the lines written in the input model file
    inputParametersFile = open(os.path.join(rootDir, 'Output', 'inputForOpenSees', '%s_input.tcl' % siteID))
    addLines1 = inputParametersFile.read()
    inputParametersFile.close()

    # Open and copy the lines written in the OpenSees model template
    if modelOption == 1:
        OpenSeesTemplate = open(os.path.join(templateDir, 'totalStress_1_HPC.tcl'))
    elif modelOption == 2:
        OpenSeesTemplate = open(os.path.join(templateDir, 'totalStress_2_HPC.tcl'))
    addLines2 = OpenSeesTemplate.read()
    OpenSeesTemplate.close()

    # Path to the OpenSees file
    OpenSeesFilePath = os.path.join(pathOutput, '%s_TotalStress_LR%.1f_%.0f.tcl' % (siteID, LR, Vs_ref))

    # Add the copied lines to the input file
    with open(OpenSeesFilePath, 'a') as OpenSeesFile:
        OpenSeesFile.write(addLines1)
        OpenSeesFile.write(addLines2)

    # Close final file
    OpenSeesFile.close()

    return

