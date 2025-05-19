################################################################## 
#                                                                # 
# Total-Stress Site-Response Analysis                            # 
# [Pore-Pressure Generation Allowed: No]                         # 
# Constitutive Models: Linear, PDMY02 (sands), and PIMY (clays)  # 
# Site: PRPC                                                     # 
#                                                                # 
################################################################## 

# Extract exterior inputs and give them a variable name 
set site [lindex $argv 0] 
set modelID [lindex $argv 1] 
set gMotionPath [lindex $argv 2] 
set gMotionName [lindex $argv 3] 
set npts [lindex $argv 4] 
set dt [lindex $argv 5] 
set saveDir [lindex $argv 6] 

set dash "_" 
set analID $site$dash$modelID$dash$gMotionName
set siteModelID $site$dash$modelID

# Define gravity and pi 
set g 9.80665 
set pi [expr atan(1)*4] 

# ----------------------------------------------------------------------------------------- 
#  1. DEFINE SOIL GEOMETERY AND MATERIAL PARAMETERS 
# ----------------------------------------------------------------------------------------- 

#---SOIL GEOMETRY 

# Thicknesses of profile (m) (counting layer 0, which is the elastic half-space) 
set soilThick 28.1 

# Number of soil layers (not counting layer 0, which is the elastic half-space) 
set numLayers 12 

# Layer thicknesses 
set layerThick(12) 0.7 
set layerThick(11) 1.5 
set layerThick(10) 0.7 
set layerThick(9) 1.1 
set layerThick(8) 4.0 
set layerThick(7) 4.0 
set layerThick(6) 4.0 
set layerThick(5) 4.0 
set layerThick(4) 2.0 
set layerThick(3) 3.0 
set layerThick(2) 2.0 
set layerThick(1) 1.0 
set layerThick(0) 0.1 

# Depth of water table (create a new layer at WT) 
# If water not present set waterTable anywhere below depth of model 
set waterTable 2.2 

#---BASIC MATERIAL PROPERTIES 

# Shear-wave velocity, Vs (m/s) 
set Vs(12) 121 
set Vs(11) 200 
set Vs(10) 140 
set Vs(9) 140 
set Vs(8) 170 
set Vs(7) 170 
set Vs(6) 240 
set Vs(5) 240 
set Vs(4) 160 
set Vs(3) 270 
set Vs(2) 170 
set Vs(1) 170 
set Vs(0) 400 

# Reference depth, refDepth (0 is ToL, 1 is BoL, 0.5 is in the middle) 
set refDepth(12) 0.5 
set refDepth(11) 1.0 
set refDepth(10) 0.5 
set refDepth(9) 0.5 
set refDepth(8) 1.0 
set refDepth(7) 0.5 
set refDepth(6) 1.0 
set refDepth(5) 0.5 
set refDepth(4) 0.5 
set refDepth(3) 0.5 
set refDepth(2) 0.5 
set refDepth(1) 0.5 
set refDepth(0) 0.5 

# Pressure dependency coefficient 
set pressCoeff(12) 0.25 
set pressCoeff(11) 0.25 
set pressCoeff(10) 0.0 
set pressCoeff(9) 0.0 
set pressCoeff(8) 0.5 
set pressCoeff(7) 0.0 
set pressCoeff(6) 0.5 
set pressCoeff(5) 0.0 
set pressCoeff(4) 0.0 
set pressCoeff(3) 0.0 
set pressCoeff(2) 0.0 
set pressCoeff(1) 0.0 
set pressCoeff(0) 0.0 

# Mass density, rho (Mg/m^3) 
set rho(12) 2.0285423037716614 
set rho(11) 1.9266055045871557 
set rho(10) 1.7533129459734962 
set rho(9) 2.038735983690112 
set rho(8) 2.059123343527013 
set rho(7) 2.0489296636085625 
set rho(6) 2.0489296636085625 
set rho(5) 2.059123343527013 
set rho(4) 1.763506625891947 
set rho(3) 2.0998980632008157 
set rho(2) 1.8960244648318043 
set rho(1) 1.763506625891947 
set rho(0) 2.038735983690112 
set rhoWater 1.0 

# Poisson ratio of the soil 
set nu 0.25 

# Rock elastic properties) 
# Bedrock shear wave velocity (m/s) 
set rockVS $Vs(0) 
# Bedrock mass density (Mg/m^3) 
set rockDen $rho(0) 

# Relative density (%) 
# Layer 12: Dr = 80 
# Layer 11: Dr = 80 
# Layer 10: Dr = 53 
# Layer 9: Dr = 92 
# Layer 8: Dr = 85 
# Layer 7: Dr = 76 
# Layer 6: Dr = 70 
# Layer 5: Dr = 67 
# Layer 4: Dr = 49 
# Layer 3: Dr = 78 
# Layer 2: Dr = 31 
# Layer 1: Dr = 15 
# Layer 0: Dr = 85 

#--- MODEL PARAMETERS 

# Consitutive model 
# constModelFlag:  0 for Linear | 1 for PDMY02 (sands) | 2 for PIDMY (clays) 
set constModelFlag(12) 1 
set constModelFlag(11) 1 
set constModelFlag(10) 1 
set constModelFlag(9) 1 
set constModelFlag(8) 1 
set constModelFlag(7) 1 
set constModelFlag(6) 1 
set constModelFlag(5) 1 
set constModelFlag(4) 2 
set constModelFlag(3) 1 
set constModelFlag(2) 1 
set constModelFlag(1) 2 

# Soil friction angle (°) 
set phi(12) 33.510854702707384 
set phi(11) 38.23617831705739 
set phi(10) 35.801777866632825 
set phi(9) 38.23617831705739 
set phi(8) 38.23617831705739 
set phi(7) 38.23617831705739 
set phi(6) 36.999370912308606 
set phi(5) 35.801777866632825 
set phi(4) 25.376933525152303 
set phi(3) 36.999370912308606 
set phi(2) 29.270831676723198 
set phi(1) 25.376933525152303 

# kNot (kPa) 
set kNot(12) 0.3843385246743418 
set kNot(11) 0.33086939364114176 
set kNot(10) 0.35721239031346075 
set kNot(9) 0.33086939364114176 
set kNot(8) 0.33086939364114176 
set kNot(7) 0.33086939364114176 
set kNot(6) 0.34394097100949284 
set kNot(5) 0.35721239031346075 
set kNot(4) 0.5 
set kNot(3) 0.34394097100949284 
set kNot(2) 0.4408070965292531 
set kNot(1) 0.5 

# Peak shear strain 
set gammaPeak 0.1 

# Phase transformation angle (not for layer 0) (°) 
set phaseAng(12) 26.2724 
set phaseAng(11) 26.2724 
set phaseAng(10) 25.7538 
set phaseAng(9) 26.5 
set phaseAng(8) 26.3859 
set phaseAng(7) 26.1816 
set phaseAng(6) 26.0454 
set phaseAng(5) 26.02 
set phaseAng(3) 26.227 
set phaseAng(2) 30.9 

# Contraction coefficients (not for layer 0) 
set contract1(12) 0.0176 
set contract1(11) 0.0176 
set contract1(10) 0.0476 
set contract1(9) 0.016 
set contract1(8) 0.0166 
set contract1(7) 0.0184 
set contract1(6) 0.0196 
set contract1(5) 0.024 
set contract1(3) 0.018000000000000002 
set contract1(2) 0.08499999999999999 
set contract2(12) 1.4724 
set contract2(11) 1.4724 
set contract2(10) 4.1538 
set contract2(9) 1.4 
set contract2(8) 1.4609 
set contract2(7) 1.4816 
set contract2(6) 1.4954 
set contract2(5) 1.7 
set contract2(3) 1.477 
set contract2(2) 4.95 
set contract3(12) 0.144 
set contract3(11) 0.144 
set contract3(10) 0.2386 
set contract3(9) 0.14 
set contract3(8) 0.1415 
set contract3(7) 0.146 
set contract3(6) 0.149 
set contract3(5) 0.16 
set contract3(3) 0.145 
set contract3(2) 0.297 

# Dilation coefficients (not for layer 0) 
set dilate1(12) 0.204 
set dilate1(11) 0.204 
set dilate1(10) 0.0624 
set dilate1(9) 0.25 
set dilate1(8) 0.22649999999999998 
set dilate1(7) 0.186 
set dilate1(6) 0.159 
set dilate1(5) 0.134 
set dilate1(3) 0.195 
set dilate1(2) 0.011 
set dilate2(12) 3.0 
set dilate2(11) 3.0 
set dilate2(10) 3.0 
set dilate2(9) 3.0 
set dilate2(8) 3.0 
set dilate2(7) 3.0 
set dilate2(6) 3.0 
set dilate2(5) 3.0 
set dilate2(3) 3.0 
set dilate2(2) 3.0 
set dilate3(12) 0.0 
set dilate3(11) 0.0 
set dilate3(10) 0.0 
set dilate3(9) 0.0 
set dilate3(8) 0.0 
set dilate3(7) 0.0 
set dilate3(6) 0.0 
set dilate3(5) 0.0 
set dilate3(3) 0.0 
set dilate3(2) 0.0 

# Cohesion (not for layer 0) (kPa) 
set cohesion(4) 85 
set cohesion(1) 100 

# Void ratio (needed for layer 0 for element definition) 
set voidR(12) 0.6116 
set voidR(11) 0.6116 
set voidR(10) 0.6907 
set voidR(9) 0.58 
set voidR(8) 0.5956 
set voidR(7) 0.6244000000000001 
set voidR(6) 0.6436000000000001 
set voidR(5) 0.652 
set voidR(4) 0.703 
set voidR(3) 0.618 
set voidR(2) 0.757 
set voidR(1) 0.76 
set voidR(0) 0.5956 

# Modulus reduction curves 
set MR_list_12 {1 0.0000001000 2 0.9984664121 3 0.0000001610 4 0.9976259989 5 0.0000002593 6 0.9963267299 7 0.0000004175 8 0.9943204310 9 0.0000006723 10 0.9912279623 11 0.0000010826 12 0.9864745796 13 0.0000017433 14 0.9791994899 15 0.0000028072 16 0.9681376469 17 0.0000045204 18 0.9514845133 19 0.0000072790 20 0.9267856986 21 0.0000117210 22 0.8909539837 23 0.0000188739 24 0.8406019961 25 0.0000303920 26 0.7729256368 27 0.0000489390 28 0.6872075579 29 0.0000788046 30 0.5864448512 31 0.0001268961 32 0.4778841875 33 0.0002043360 34 0.3713733159 35 0.0003290345 36 0.2760503558 37 0.0005298317 38 0.1975549445 39 0.0008531679 40 0.1392969418 41 0.0013738238 42 0.0970553822 43 0.0022122163 44 0.0662726590 45 0.0035622479 46 0.0442390768 47 0.0057361525 48 0.0289284389 49 0.0092367086 50 0.0186087615 51 0.0148735211 52 0.0118280915 53 0.0239502662 54 0.0074565404 55 0.0385662042 56 0.0046751052 57 0.0621016942 58 0.0029208757 59 0.1000000000 60 0.0018207844 }
set MR_list_11 {1 0.0000001000 2 0.9989984313 3 0.0000001610 4 0.9984491143 5 0.0000002593 6 0.9975992447 7 0.0000004175 8 0.9962853878 9 0.0000006723 10 0.9942566386 11 0.0000010826 12 0.9911297454 13 0.0000017433 14 0.9863238752 15 0.0000028072 16 0.9789694525 17 0.0000045204 18 0.9677892983 19 0.0000072790 20 0.9509633191 21 0.0000117210 22 0.9260197994 23 0.0000188739 24 0.8898577937 25 0.0000303920 26 0.8390910555 27 0.0000489390 28 0.7709479905 29 0.0000788046 30 0.6847878234 31 0.0001268961 32 0.5837178093 33 0.0002043360 34 0.4750820406 35 0.0003290345 36 0.3687546162 37 0.0005298317 38 0.2739093204 39 0.0008531679 40 0.2007898166 41 0.0013738238 42 0.1460772752 43 0.0022122163 44 0.1040766249 45 0.0035622479 46 0.0721132109 47 0.0057361525 48 0.0485748395 49 0.0092367086 50 0.0319364061 51 0.0148735211 52 0.0206100883 53 0.0239502662 54 0.0131254276 55 0.0385662042 56 0.0082839891 57 0.0621016942 58 0.0051975550 59 0.1000000000 60 0.0032486851 }
set MR_list_10 {1 0.0000001000 2 0.9993477895 3 0.0000001610 4 0.9989898866 5 0.0000002593 6 0.9984358904 7 0.0000004175 8 0.9975787919 9 0.0000006723 10 0.9962537838 11 0.0000010826 12 0.9942078743 13 0.0000017433 14 0.9910546709 15 0.0000028072 16 0.9862086916 17 0.0000045204 18 0.9787936599 19 0.0000072790 20 0.9675231521 21 0.0000117210 22 0.9505652470 23 0.0000188739 24 0.9254351166 25 0.0000303920 26 0.8890215654 27 0.0000489390 28 0.8379395880 29 0.0000788046 30 0.7694428769 31 0.0001268961 32 0.6829493799 33 0.0002043360 34 0.5816500101 35 0.0003290345 36 0.4729618431 37 0.0005298317 38 0.3667773755 39 0.0008531679 40 0.2721234243 41 0.0013738238 42 0.1973993882 43 0.0022122163 44 0.1473297416 45 0.0035622479 46 0.1123619008 47 0.0057361525 48 0.0860933053 49 0.0092367086 50 0.0650692024 51 0.0148735211 52 0.0478201538 53 0.0239502662 54 0.0339351439 55 0.0385662042 56 0.0232623511 57 0.0621016942 58 0.0154832215 59 0.1000000000 60 0.0100753259 }
set MR_list_9 {1 0.0000001000 2 0.9991826324 3 0.0000001610 4 0.9987342137 5 0.0000002593 6 0.9980402690 7 0.0000004175 8 0.9969670369 9 0.0000006723 10 0.9953088199 11 0.0000010826 12 0.9927505954 13 0.0000017433 14 0.9888129836 15 0.0000028072 16 0.9827737170 17 0.0000045204 18 0.9735613486 19 0.0000072790 20 0.9596247391 21 0.0000117210 22 0.9388035062 23 0.0000188739 24 0.9082712767 25 0.0000303920 26 0.8647012494 27 0.0000489390 28 0.8048818155 29 0.0000788046 30 0.7269655751 31 0.0001268961 32 0.6321548860 33 0.0002043360 34 0.5258927437 35 0.0003290345 36 0.4172324098 37 0.0005298317 38 0.3160566066 39 0.0008531679 40 0.2297430607 41 0.0013738238 42 0.1638884034 43 0.0022122163 44 0.1209734845 45 0.0035622479 46 0.0922743925 47 0.0057361525 48 0.0717878885 49 0.0092367086 50 0.0559266205 51 0.0148735211 52 0.0428188119 53 0.0239502662 54 0.0317691032 55 0.0385662042 56 0.0226971921 57 0.0621016942 58 0.0156309275 59 0.1000000000 60 0.0104360002 }
set MR_list_8 {1 0.0000001000 2 0.9992816490 3 0.0000001610 4 0.9988874917 5 0.0000002593 6 0.9982774340 7 0.0000004175 8 0.9973337365 9 0.0000006723 10 0.9958751762 11 0.0000010826 12 0.9936238222 13 0.0000017433 14 0.9901558121 15 0.0000028072 16 0.9848303449 17 0.0000045204 18 0.9766917465 19 0.0000072790 20 0.9643448307 21 0.0000117210 22 0.9458202868 23 0.0000188739 24 0.9184849064 25 0.0000303920 26 0.8791205755 27 0.0000489390 28 0.8243815521 29 0.0000788046 30 0.7518514067 31 0.0001268961 32 0.6616603775 33 0.0002043360 34 0.5579616137 35 0.0003290345 36 0.4489496791 37 0.0005298317 38 0.3446310153 39 0.0008531679 40 0.2534051922 41 0.0013738238 42 0.1825099932 43 0.0022122163 44 0.1357783456 45 0.0035622479 46 0.1039388088 47 0.0057361525 48 0.0806378936 49 0.0092367086 50 0.0622144401 51 0.0148735211 52 0.0469134725 53 0.0239502662 54 0.0341902257 55 0.0385662042 56 0.0240038380 57 0.0621016942 58 0.0162861182 59 0.1000000000 60 0.0107490545 }
set MR_list_7 {1 0.0000001000 2 0.9993717906 3 0.0000001610 4 0.9990270455 5 0.0000002593 6 0.9984933985 7 0.0000004175 8 0.9976677394 9 0.0000006723 10 0.9963912320 11 0.0000010826 12 0.9944199661 13 0.0000017433 14 0.9913812235 15 0.0000028072 16 0.9867097762 17 0.0000045204 18 0.9795585733 19 0.0000072790 20 0.9686815854 21 0.0000117210 22 0.9522987385 23 0.0000188739 24 0.9279830644 25 0.0000303920 26 0.8926694864 27 0.0000489390 28 0.8429700191 29 0.0000788046 30 0.7760311843 31 0.0001268961 32 0.6910167752 33 0.0002043360 34 0.5907503984 35 0.0003290345 36 0.4823222668 37 0.0005298317 38 0.3755336923 39 0.0008531679 40 0.2796178573 41 0.0013738238 42 0.2037651346 43 0.0022122163 44 0.1540429109 45 0.0035622479 46 0.1205262108 47 0.0057361525 48 0.0962050880 49 0.0092367086 50 0.0768470273 51 0.0148735211 52 0.0602457562 53 0.0239502662 54 0.0456709835 55 0.0385662042 56 0.0332359467 57 0.0621016942 58 0.0232280364 59 0.1000000000 60 0.0156797152 }
set MR_list_6 {1 0.0000001000 2 0.9994308895 3 0.0000001610 4 0.9991185477 5 0.0000002593 6 0.9986350193 7 0.0000004175 8 0.9978868081 9 0.0000006723 10 0.9967298096 11 0.0000010826 12 0.9949425509 13 0.0000017433 14 0.9921861582 15 0.0000028072 16 0.9879456874 17 0.0000045204 18 0.9814469992 19 0.0000072790 20 0.9715456637 21 0.0000117210 22 0.9565938986 23 0.0000188739 24 0.9343166644 25 0.0000303920 26 0.9017800601 27 0.0000489390 28 0.8556169117 29 0.0000788046 30 0.7927437165 31 0.0001268961 32 0.7117164636 33 0.0002043360 34 0.6144191245 35 0.0003290345 36 0.5070289315 37 0.0005298317 38 0.3989861709 39 0.0008531679 40 0.2999575623 41 0.0013738238 42 0.2200768143 43 0.0022122163 44 0.1658204090 45 0.0035622479 46 0.1271939818 47 0.0057361525 48 0.0975574562 49 0.0092367086 50 0.0734903647 51 0.0148735211 52 0.0536804036 53 0.0239502662 54 0.0378271601 55 0.0385662042 56 0.0257625144 57 0.0621016942 58 0.0170578166 59 0.1000000000 60 0.0110572897 }
set MR_list_5 {1 0.0000001000 2 0.9994734636 3 0.0000001610 4 0.9991844685 5 0.0000002593 6 0.9987370559 7 0.0000004175 8 0.9980446663 9 0.0000006723 10 0.9969738350 11 0.0000010826 12 0.9953193173 13 0.0000017433 14 0.9927667757 15 0.0000028072 16 0.9888378536 17 0.0000045204 18 0.9828117796 19 0.0000072790 20 0.9736192201 21 0.0000117210 22 0.9597118541 23 0.0000188739 24 0.9389326869 25 0.0000303920 26 0.9084586239 27 0.0000489390 28 0.8649643542 29 0.0000788046 30 0.8052350452 31 0.0001268961 32 0.7274120887 33 0.0002043360 34 0.6326781062 35 0.0003290345 36 0.5264538869 37 0.0005298317 38 0.4177797780 39 0.0008531679 40 0.3165433388 41 0.0013738238 42 0.2339742646 43 0.0022122163 44 0.1780071838 45 0.0035622479 46 0.1382449046 47 0.0057361525 48 0.1076259673 49 0.0092367086 50 0.0824249341 51 0.0148735211 52 0.0612131117 53 0.0239502662 54 0.0437873329 55 0.0385662042 56 0.0301913661 57 0.0621016942 58 0.0201773615 59 0.1000000000 60 0.0131662655 }
set MR_list_4 {1 0.0000001000 2 0.9997271016 3 0.0000001610 4 0.9995772595 5 0.0000002593 6 0.9993451965 7 0.0000004175 8 0.9989858721 9 0.0000006723 10 0.9984296777 11 0.0000010826 12 0.9975691829 13 0.0000017433 14 0.9962389362 15 0.0000028072 16 0.9941849653 17 0.0000045204 18 0.9910194030 19 0.0000072790 20 0.9861545850 21 0.0000117210 22 0.9787110902 23 0.0000188739 24 0.9673981609 25 0.0000303920 26 0.9503783384 27 0.0000489390 28 0.9251606739 29 0.0000788046 30 0.8886292282 31 0.0001268961 32 0.8373996921 33 0.0002043360 34 0.7687377666 35 0.0003290345 36 0.6820890402 37 0.0005298317 38 0.5806835597 39 0.0008531679 40 0.4719722451 41 0.0013738238 42 0.3713924962 43 0.0022122163 44 0.2948013448 45 0.0035622479 46 0.2323905283 47 0.0057361525 48 0.1786918524 49 0.0092367086 50 0.1324540976 51 0.0148735211 52 0.0942878822 53 0.0239502662 54 0.0646545298 55 0.0385662042 56 0.0430010555 57 0.0621016942 58 0.0279550978 59 0.1000000000 60 0.0178857173 }
set MR_list_3 {1 0.0000001000 2 0.9995079811 3 0.0000001610 4 0.9992379169 5 0.0000002593 6 0.9988197921 7 0.0000004175 8 0.9981726784 9 0.0000006723 10 0.9971717538 11 0.0000010826 12 0.9956249709 13 0.0000017433 14 0.9932379805 15 0.0000028072 16 0.9895623152 17 0.0000045204 18 0.9839209943 19 0.0000072790 20 0.9753067573 21 0.0000117210 22 0.9622545341 23 0.0000188739 24 0.9427085019 25 0.0000303920 26 0.9139459710 27 0.0000489390 28 0.8726935564 29 0.0000788046 30 0.8156544156 31 0.0001268961 32 0.7406538685 33 0.0002043360 34 0.6482968168 35 0.0003290345 36 0.5433291180 37 0.0005298317 38 0.4343667985 39 0.0008531679 40 0.3313988088 41 0.0013738238 42 0.2464137849 43 0.0022122163 44 0.1882273393 45 0.0035622479 46 0.1463076462 47 0.0057361525 48 0.1135829314 49 0.0092367086 50 0.0864704998 51 0.0148735211 52 0.0637214815 53 0.0239502662 54 0.0452189234 55 0.0385662042 56 0.0309597379 57 0.0621016942 58 0.0205764694 59 0.1000000000 60 0.0133726624 }
set MR_list_2 {1 0.0000001000 2 0.9995917501 3 0.0000001610 4 0.9993676367 5 0.0000002593 6 0.9990206143 7 0.0000004175 8 0.9984834452 9 0.0000006723 10 0.9976523444 11 0.0000010826 12 0.9963674415 13 0.0000017433 14 0.9943832535 15 0.0000028072 16 0.9913246926 17 0.0000045204 18 0.9866230185 19 0.0000072790 20 0.9794261064 21 0.0000117210 22 0.9684809003 23 0.0000188739 24 0.9519982759 25 0.0000303920 26 0.9275410951 27 0.0000489390 28 0.8920360107 29 0.0000788046 30 0.8420950996 31 0.0001268961 32 0.7748829038 33 0.0002043360 34 0.6896069642 35 0.0003290345 36 0.5891550988 37 0.0005298317 38 0.4810103203 39 0.0008531679 40 0.3957932678 41 0.0013738238 42 0.3365679650 43 0.0022122163 44 0.2908712486 45 0.0035622479 46 0.2506959965 47 0.0057361525 48 0.2115241779 49 0.0092367086 50 0.1720082897 51 0.0148735211 52 0.1334506535 53 0.0239502662 54 0.0984909523 55 0.0385662042 56 0.0694086585 57 0.0621016942 58 0.0470818579 59 0.1000000000 60 0.0310272221 }
set MR_list_1 {1 0.0000001000 2 0.9997468267 3 0.0000001610 4 0.9996078109 5 0.0000002593 6 0.9993925089 7 0.0000004175 8 0.9990591228 9 0.0000006723 10 0.9985430439 11 0.0000010826 12 0.9977445305 13 0.0000017433 14 0.9965099058 15 0.0000028072 16 0.9946031136 17 0.0000045204 18 0.9916632722 19 0.0000072790 20 0.9871427150 21 0.0000117210 22 0.9802197979 23 0.0000188739 24 0.9696837619 25 0.0000303920 26 0.9538001562 27 0.0000489390 28 0.9301937264 29 0.0000788046 30 0.8958424797 31 0.0001268961 32 0.8473610386 33 0.0002043360 34 0.7818095206 35 0.0003290345 36 0.6981352801 37 0.0005298317 38 0.5988378522 39 0.0008531679 40 0.4907051694 41 0.0013738238 42 0.3890753123 43 0.0022122163 44 0.3103805722 45 0.0035622479 46 0.2451873020 47 0.0057361525 48 0.1884818634 49 0.0092367086 50 0.1394679665 51 0.0148735211 52 0.0990520274 53 0.0239502662 54 0.0677707466 55 0.0385662042 56 0.0449918901 57 0.0621016942 58 0.0292101042 59 0.1000000000 60 0.0186712695 }

# Flags for water table and Vs inversion 
# set VsInvTopLayer to "Yes" if there is a velocity inversion immediately below upper layer (else "No") 
# set waterTopLayer to "Yes" if water table within upper most layer and layer was split in two (else "No") 
# if waterTopLayer == "Yes", should set refDepth(numLayers) = 1.0 and refDepth(numLayers-1) = 0.0 
set VsInvTopLayer No 
set waterTopLayer No 

# Rayleigh damping parameters) 
set damp 0.02674533368657426 
set omega1 10.628383251856668 
set omega2 53.14191625928334 

# Allow excess pore pressure generation? Yes or No
# If No, permeability is automatically set very high for dynamic analysis
set allowPWP       No

# Define layer boundaries
set layerBound(0) $layerThick(0)
puts "layer boundary 0 = $layerBound(0)"
for {set i 1} {$i <= $numLayers} {incr i 1} {
    set layerBound($i) [expr $layerBound([expr $i-1])+$layerThick($i)]
    puts "layer boundary $i = $layerBound($i)"
}

# Reference pressure
# computed as mean confining pressure at refDepth for each layer (0 is ToL, 1 is BoL)

## for top layer
if {$layerThick($numLayers) > $waterTable} {
    set vertStress($numLayers) [expr ($rho($numLayers) - $rhoWater) * $g * $layerThick($numLayers) * $refDepth($numLayers)]
}   else {
    set vertStress($numLayers) [expr $rho($numLayers) * $g * $layerThick($numLayers) * $refDepth($numLayers)]
}
#set kNot($numLayers) [expr (1.0 - sin($phi($numLayers) * 2.0 * $pi / 360.0))]
set meanStress [expr $vertStress($numLayers) * (1.0 + 2.0 * $kNot($numLayers)) / 3.0]
set refPress($numLayers) $meanStress
## for other layers (not necessary for layer 0 bc it's lin elastic bedrock)
set bottomDepth($numLayers) $layerThick($numLayers)

for {set k [expr $numLayers - 1]} {$k > 0 && $k <= [expr $numLayers - 1]} {incr k -1} {
    set bottomDepth($k) [expr $bottomDepth([expr $k + 1]) + $layerThick($k)]
    if {$bottomDepth($k) > $waterTable && $bottomDepth([expr $k + 1]) > $waterTable} {
        set vertStress($k) [expr $vertStress([expr $k + 1]) + ($rho([expr $k + 1]) - $rhoWater) * $g * $layerThick([expr $k + 1]) * (1.0 - $refDepth([expr $k + 1])) + ($rho($k) - $rhoWater) * $g * $layerThick($k) * $refDepth($k)]
} elseif {$bottomDepth($k) > $waterTable} {
    set vertStress($k) [expr $vertStress([expr $k + 1]) + $rho([expr $k + 1]) * $g * $layerThick([expr $k + 1]) * (1.0 - $refDepth([expr $k + 1])) + ($rho($k) - $rhoWater) * $g * $layerThick($k) * $refDepth($k)]
} else {
    set vertStress($k) [expr $vertStress([expr $k + 1]) + $rho([expr $k + 1]) * $g * $layerThick([expr $k + 1]) * (1.0 - $refDepth([expr $k + 1])) + $rho($k) * $g * $layerThick($k) * $refDepth($k)]
}
#set kNot($k) [expr 1.0 - sin($phi($k) * 2.0 * $pi / 360.0)]
set meanStress [expr $vertStress($k) * (1.0 + 2.0 * $kNot($k)) / 3.0]
set refPress($k) $meanStress
}

# Compute Vs_not, the constant required for an exponential function of Vs to have
# equal travel time to a layer of constant Vs
set Vs_not [expr $Vs($numLayers)*pow($layerThick($numLayers), -$pressCoeff($numLayers) / 2.0) / (1.0 - $pressCoeff($numLayers) / 2.0)]
puts "Vs_not = $Vs_not"

# Soil shear modulus for each layer (kPa)

# for top layer
if {$VsInvTopLayer == "Yes" || $waterTopLayer == "Yes"} {
    set G($numLayers) [expr $rho($numLayers)*$Vs($numLayers)*$Vs($numLayers)]
} else {
    set G($numLayers) [expr $rho($numLayers) * pow($Vs_not, 2.0 ) * pow($rho($numLayers) * $g  * (1.0 + 2.0 * $kNot($numLayers)) / (3.0 * $refPress($numLayers)), -$pressCoeff($numLayers)) ]
}
# for all other layers    
for {set k 0} {$k < $numLayers} {incr k 1} {
    set G($k) [expr $rho($k)*$Vs($k)*$Vs($k)]    
}

# Soil elastic modulus for each layer (kPa)
for {set k 0} {$k <= $numLayers} {incr k 1} {
    set E($k)       [expr 2*$G($k)*(1+$nu)]
}
# Soil bulk modulus for each layer (kPa)
for {set k 0} {$k <= $numLayers} {incr k 1} {
    set bulk($k)    [expr $E($k)/(3*(1-2*$nu))]
}

# Print the value of Shear Modulus to output
if {[file exists "$saveDir/$siteModelID.nodeInfo.txt"] == "0"} {
    parray G
} 

#-----------------------------------------------------------------------------------------
#  2. MESH GEOMETRY
#-----------------------------------------------------------------------------------------

#---MESH GEOMETRY
# highest frequency desired to be well resolved (Hz)
set fMax    25.0
# Wavelength of highest resolved frequency for each layer
for {set k 0} {$k <= $numLayers} {incr k 1} {
    set wave($k) [expr $Vs($k)/$fMax]
}
# Number of elements per one wavelength
set nEle     8

# Determine number of elements in column
set nElemT 0
for {set k 0} {$k <= $numLayers} {incr k 1} {

    # maximum element height
    set hEleMax($k) [expr $wave($k)/$nEle]

    # interger number of elements
    set nElemY($k) [expr int(floor($layerThick($k)/$hEleMax($k))+1)]
    puts "number of vertical elements in layer $k: $nElemY($k)"

    set nElemT [expr $nElemT + $nElemY($k)]

    # actual element height
    set sElemY($k) [expr {$layerThick($k)/$nElemY($k)}] 

    puts "vertical size of elements in layer $k: $sElemY($k)"
}
puts "total number of vertical elements: $nElemT"

# Number of nodes in vertical direction
set nNodeY  [expr 2*$nElemT+1]

# Number of elements in horizontal direction
set nElemX  1
# Number of nodes in horizontal direction
set nNodeX  [expr 2*$nElemX+1]
# Horizontal element size (m)
set sElemX  2.0

# Total number of nodes
set nNodeT  [expr $nNodeX*$nNodeY]

#-----------------------------------------------------------------------------------------
#  3. CREATE PORE PRESSURE NODES AND FIXITIES
#-----------------------------------------------------------------------------------------

model BasicBuilder -ndm 2 -ndf 3

set ppNodesInfo [open $saveDir/$siteModelID.ppNodesInfo.dat w]
set LeftPPNodesInfo [open $saveDir/$siteModelID.LeftPPNodesInfo.dat w]
set dryNodeCount 1
set PPNodeCount 1
set layerNodeCount 0
set yCoordCount 0
# loop over soil layers
for {set k 0} {$k <= $numLayers} {incr k 1} {
    # loop in horizontal direction
    for {set i 1} {$i <= $nNodeX} {incr i 2} {
        # loop in vertical direction
        if {$k == 0} {
            set bump 1
    } else {
        set bump 0
    }
    for {set j 1} {$j <= [expr 2*$nElemY($k)+$bump]} {incr j 2} {

        set xCoord       [expr ($i-1)*$sElemX/2.0]
        set yctr    [expr $j + $layerNodeCount]
        if {$k == 0} {
            set yCoord       [expr ($j-1)*$sElemY($k)/2.0]
      } else {
          set yCoord       [expr $layerBound([expr $k - 1]) + $sElemY($k) + ($j - 1)*$sElemY($k)/2.0]
      }
      set nodeNum      [expr $i + ($yctr-1)*$nNodeX]

      node $nodeNum  $xCoord  $yCoord

      # puts "yctr = $yctr"
      # puts "xCoord = $xCoord"
      # puts "yCoord = $yCoord"
      # puts "nodeNum = $nodeNum"

      set PPNode($PPNodeCount) $nodeNum
      set PPNodeCount [expr $PPNodeCount + 1]

      # output nodal information to data file
      puts $ppNodesInfo "$nodeNum  $xCoord  $yCoord"

      if {$xCoord == 0} {
          puts $LeftPPNodesInfo "$nodeNum  $xCoord  $yCoord"
          #puts "this is a PP Node"
      }

      # designate nodes above water table
      set waterHeight [expr $soilThick-$waterTable]
      if {$yCoord>=$waterHeight} {
          set dryNode($dryNodeCount) $nodeNum
          set dryNodeCount [expr $dryNodeCount+1]
        }
    }
}
set layerNodeCount [expr $yctr + 1]
}
close $ppNodesInfo
close $LeftPPNodesInfo
puts "Finished creating all -ndf 3 nodes..."

# define fixities for pore pressure nodes above water table
for {set i 1} {$i < $dryNodeCount} {incr i 1} {
    fix $dryNode($i)  0 0 1
}

# define fixities for pore pressure nodes at base of soil column
fix 1  0 1 0
fix 3  0 1 0
puts "Finished creating all -ndf 3 boundary conditions..."


# define equal degrees of freedom for pore pressure nodes
for {set i 1} {$i <= [expr 3*$nNodeY-2]} {incr i 6} {
    equalDOF $i [expr $i+2]  1 2
}
puts "Finished creating equalDOF for pore pressure nodes..."

if {[file exists "$saveDir/$siteModelID.PPnodeInfo.txt"] == "0"} {
    print $saveDir/$siteModelID.PPnodeInfo.txt -node
} 
#-----------------------------------------------------------------------------------------
#  4. CREATE INTERIOR NODES AND FIXITIES
#-----------------------------------------------------------------------------------------

model BasicBuilder -ndm 2 -ndf 2

# central column of nodes
set xCoord  [expr $sElemX/2]
# loop over soil layers
set layerNodeCount 0
for {set k 0} {$k <= $numLayers} {incr k 1} {
    # loop in vertical direction
    if {$k == 0} {
        set bump 1
} else {
    set bump 0
}
for {set j 1} {$j <= [expr 2 * $nElemY($k) + $bump]} {incr j 1} {

    set yctr    [expr $j + $layerNodeCount]
    if {$k == 0} {
        set yCoord  [expr ($j - 1) * $sElemY($k) / 2.0]
  } else {
      set yCoord  [expr $layerBound([expr $k - 1]) + $sElemY($k) / 2.0 +  ($j - 1) * $sElemY($k) / 2.0]
  }      
  set nodeNum [expr 3 * $yctr - 1] 

  node  $nodeNum  $xCoord  $yCoord 
}
set layerNodeCount $yctr
}

# Interior nodes on the element edges
# loop over layers
set layerNodeCount 0
for {set k 0} {$k <= $numLayers} {incr k 1} {
    # loop in vertical direction
    for {set j 1} {$j <= $nElemY($k)} {incr j 1} {

        set yctr [expr $j + $layerNodeCount]
        if {$k == 0} {
            set yCoord   [expr $sElemY($k) * ($j - 0.5)]
  } else {
      set yCoord   [expr $layerBound([expr $k - 1]) + ($j - 0.5) * $sElemY($k)]
  }
  set nodeNumL [expr 6*$yctr - 2]
  set nodeNumR [expr $nodeNumL + 2]

  node  $nodeNumL  0.0  $yCoord
  node  $nodeNumR  $sElemX  $yCoord
}
set layerNodeCount $yctr
}
puts "Finished creating all -ndf 2 nodes..."

# Define fixities for interior node at base of soil column
fix 2  0 1
puts "Finished creating all -ndf 2 boundary conditions..."

# Define equal degrees of freedom which have not yet been defined
for {set i 1} {$i <= [expr 3*$nNodeY-6]} {incr i 6} {
    equalDOF $i          [expr $i+1]  1 2
    equalDOF [expr $i+3] [expr $i+4]  1 2
    equalDOF [expr $i+3] [expr $i+5]  1 2
}
equalDOF [expr $nNodeT-2] [expr $nNodeT-1]  1 2
puts "Finished creating equalDOF constraints..."

# Print all node info to text file (if statement so it only does it once when batching)
puts "ground motion is $gMotionName"
if {[file exists "$saveDir/$siteModelID.nodeInfo.txt"] == "0"} {
    print $saveDir/$siteModelID.nodeInfo.txt -node
} 

#-----------------------------------------------------------------------------------------
#  5. CREATE SOIL MATERIALS
#-----------------------------------------------------------------------------------------

# Define grade of slope (%)
set grade 0.0
set slope [expr atan($grade/100.0)]
set g -9.81
set bulkWater 2.2e6
set materialProperties [open $saveDir/$siteModelID.materialProperties.dat w]

# Define nonlinear material for soil
for {set i 1} {$i <= $numLayers} {incr i 1} {

    if {$constModelFlag($i) == 0} {
	
		# Linear model
		nDMaterial ElasticIsotropic $i $E($i) $nu $rho($i)

		set thick($i) 1.0
		set xWgt($i)  [expr $g*sin($slope)]
		set yWgt($i)  [expr $g*cos($slope)]
		set porosity($i) [expr $voidR($i) / (1 + $voidR($i))]
		set uBulk($i) [expr $bulkWater/$porosity($i)]
		set hPerm($i) 1.0e-4
		set vPerm($i) 1.0e-4

        # Output material information to data file
        puts $materialProperties "
        Material ID: $i
		E: $E($i) 
		nu: $nu
        rho: $rho($i) 
        "

} elseif {$constModelFlag($i) == 1} {

		array set MR {}
		array set MR [set MR_list_$i]
    
		# PDMY02 Model
        nDMaterial PressureDependMultiYield02 $i 2 $rho($i) $G($i) $bulk($i) $phi($i) $gammaPeak \
        $refPress($i) $pressCoeff($i) $phaseAng($i) \
        $contract1($i) $contract3($i) $dilate1($i) $dilate3($i) \
        -30 $MR(1) $MR(2) $MR(3) $MR(4) $MR(5) $MR(6) $MR(7) $MR(8) $MR(9) $MR(10) \
		$MR(11) $MR(12) $MR(13) $MR(14) $MR(15) $MR(16) $MR(17) $MR(18) $MR(19) $MR(20) \
		$MR(21) $MR(22) $MR(23) $MR(24) $MR(25) $MR(26) $MR(27) $MR(28) $MR(29) $MR(30) \
		$MR(31) $MR(32) $MR(33) $MR(34) $MR(35) $MR(36) $MR(37) $MR(38) $MR(39) $MR(40) \
		$MR(41) $MR(42) $MR(43) $MR(44) $MR(45) $MR(46) $MR(47) $MR(48) $MR(49) $MR(50) \
		$MR(51) $MR(52) $MR(53) $MR(54) $MR(55) $MR(56) $MR(57) $MR(58) $MR(59) $MR(60) \ 
		$contract2($i) $dilate2($i) 1 0 $voidR($i) 0.9 0.02 0.7 101.0  

        set thick($i) 1.0
        set xWgt($i)  [expr $g*sin($slope)]
        set yWgt($i)  [expr $g*cos($slope)]
        set porosity($i) [expr $voidR($i) / (1 + $voidR($i))]
        set uBulk($i) [expr $bulkWater/$porosity($i)]
        set hPerm($i) 1.0e-4
        set vPerm($i) 1.0e-4

        # Output material information to data file
        puts $materialProperties "
        Material ID: $i
        rho: $rho($i) 
        G: $G($i)
        K: $bulk($i)
        Phi: $phi($i)
        GammaPeak: $gammaPeak
        refPress: $refPress($i)
        pressCoeff: $pressCoeff($i)
        phaseAng: $phaseAng($i)
        contract1: $contract1($i)
        contract3: $contract3($i)
        dilate1: $dilate1($i)
        dilate3: $dilate3($i)
        "

  }  elseif {$constModelFlag($i) == 2} {
  
		array set MR {}
		array set MR [set MR_list_$i]
    
		# PIMY Model
		nDMaterial PressureIndependMultiYield $i 2 $rho($i) $G($i) $bulk($i) $cohesion($i) $gammaPeak 0.0 \
		$refPress($i) $pressCoeff($i) -30 $MR(1) $MR(2) $MR(3) $MR(4) $MR(5) $MR(6) $MR(7) $MR(8) $MR(9) $MR(10) \
		$MR(11) $MR(12) $MR(13) $MR(14) $MR(15) $MR(16) $MR(17) $MR(18) $MR(19) $MR(20) \
		$MR(21) $MR(22) $MR(23) $MR(24) $MR(25) $MR(26) $MR(27) $MR(28) $MR(29) $MR(30) \
		$MR(31) $MR(32) $MR(33) $MR(34) $MR(35) $MR(36) $MR(37) $MR(38) $MR(39) $MR(40) \
		$MR(41) $MR(42) $MR(43) $MR(44) $MR(45) $MR(46) $MR(47) $MR(48) $MR(49) $MR(50) \
		$MR(51) $MR(52) $MR(53) $MR(54) $MR(55) $MR(56) $MR(57) $MR(58) $MR(59) $MR(60) \ 

		set thick($i) 1.0
		set xWgt($i)  [expr $g*sin($slope)]
		set yWgt($i)  [expr $g*cos($slope)]
		set porosity($i) [expr $voidR($i) / (1 + $voidR($i))]
		set uBulk($i) [expr $bulkWater/$porosity($i)]
		set hPerm($i) 1.0e-4
		set vPerm($i) 1.0e-4

		# Output material information to data file
		puts $materialProperties "
		Material ID: $i
		rho: $rho($i) 
		G: $G($i)
		K: $bulk($i)
		coheison: $cohesion($i)
		Phi: 0.0
		GammaPeak: $gammaPeak
		refPress: $refPress($i)
		pressCoeff: $pressCoeff($i) 
		"

  }
  
}

close $materialProperties

# Define linear elastic material for "bedrock"
nDMaterial ElasticIsotropic 0 $E(0) 0.3 $rho(0)

set thick(0) 1.0
set xWgt(0)  [expr $g*sin($slope)]
set yWgt(0)  [expr $g*cos($slope)]
set porosity(0) [expr $voidR(0) / (1 + $voidR(0))]
set uBulk(0) [expr $bulkWater/$porosity(0)]
set hPerm(0) 1.0e-4
set vPerm(0) 1.0e-4

puts "Finished creating all soil materials..."

#-----------------------------------------------------------------------------------------
#  6. CREATE SOIL ELEMENTS
#-----------------------------------------------------------------------------------------

# Define the top element number of each layer
set elemCount 0.0
set layerTopEleNum(-1) 0
for {set i 0} {$i <= $numLayers} {incr i 1} {
    set layerTopEleNum($i) [expr $elemCount + $nElemY($i)]
    set elemCount [expr $elemCount + $nElemY($i)]
}

for {set j 1} {$j <= $nElemT} {incr j 1} {

    set nI  [expr 6*$j - 5]
    set nJ  [expr $nI + 2]
    set nK  [expr $nI + 8]
    set nL  [expr $nI + 6]
    set nM  [expr $nI + 1]
    set nN  [expr $nI + 5]
    set nP  [expr $nI + 7]
    set nQ  [expr $nI + 3]
    set nR  [expr $nI + 4]

    set lowerBound 0.0
    for {set i 0} {$i <= $numLayers} {incr i 1} {

        if {$j <= $layerTopEleNum($i) && $j > $layerTopEleNum([expr $i - 1])} {

            # permeabilities are initially set at 10.0 m/s for gravity analysis, values are updated post-gravity
            element 9_4_QuadUP $j $nI $nJ $nK $nL $nM $nN $nP $nQ $nR \
            $thick($i) $i $uBulk($i) 1.0 10.0 10.0 $xWgt($i) $yWgt($i)

            puts "element 9_4_QuadUP $j $nI $nJ $nK $nL $nM $nN $nP $nQ $nR $thick($i) $i $uBulk($i) 1.0 1.0 1.0 $xWgt($i) $yWgt($i)"

    }
    set lowerBound $layerBound($i)
}
}
puts "Finished creating all soil elements..."

# print all node info to text file (if statement so it only does it once when batching)
if {[file exists "$saveDir/$siteModelID.eleInfo.txt"] == "0"} {
     print $saveDir/$siteModelID.eleInfo.txt -ele
}
#-----------------------------------------------------------------------------------------
#  7. LYSMER DASHPOT
#-----------------------------------------------------------------------------------------

# Define dashpot nodes
set dashF [expr $nNodeT+1]
set dashS [expr $nNodeT+2]
#    nodeTag x   y
node $dashF  0.0 0.0
node $dashS  0.0 0.0
# Define fixities for dashpot nodes
#   nodeTag dof1 dof2
fix $dashF  1    1
fix $dashS  0    1

# Define equal DOF for dashpot and base soil node
#        rNodeTag cNodeTag dof1
equalDOF 1        $dashS   1
puts "Finished creating dashpot nodes and boundary conditions..."

# Define dashpot material
set colArea       [expr $sElemX*$thick(0)]
set dashpotCoeff  [expr $rockVS*$rockDen]
#                         matTag              C: damping coefficient       alpha (1=linear damping)
uniaxialMaterial Viscous [expr $numLayers+1] [expr $dashpotCoeff*$colArea] 1

# Define dashpot element
#                   eleTag           iNode  jNode        matTag                   direction (1=X)
element zeroLength [expr $nElemT+1]  $dashF $dashS -mat [expr $numLayers+1]  -dir 1
puts "Finished creating dashpot material and element..."

#-----------------------------------------------------------------------------------------
#  8. CREATE GRAVITY RECORDERS
#-----------------------------------------------------------------------------------------

# Create list for pore pressure nodes
set nodeList3 {}
set channel [open "$saveDir/$siteModelID.ppNodesInfo.dat" r]
set count 0;
foreach line [split [read -nonewline $channel] \n] {
    set count [expr $count+1];
    set lineData($count) $line
    set nodeNumber [lindex $lineData($count) 0]
    lappend nodeList3 $nodeNumber
}
set nodeList3 [lsort -integer $nodeList3]
puts "pore pressure nodes are: $nodeList3"
close $channel

# record nodal displacment, acceleration, and porepressure
#eval "recorder Node -file Gdisplacement.out -time -node $nodeList3 -dof 1 2  disp"
#eval "recorder Node -file Gacceleration.out -time -node $nodeList3 -dof 1 2  accel"
#eval "recorder Node -file GporePressure.out -time -node $nodeList3 -dof 3 vel"
# record elemental stress and strain (files are names to reflect GiD gp numbering)
#recorder Element -file Gstress1.out   -time  -eleRange 1 $nElemT  material 1 stress
#recorder Element -file Gstress2.out   -time  -eleRange 1 $nElemT  material 2 stress
#recorder Element -file Gstress3.out   -time  -eleRange 1 $nElemT  material 3 stress
#recorder Element -file Gstress4.out   -time  -eleRange 1 $nElemT  material 4 stress
#recorder Element -file Gstress9.out   -time  -eleRange 1 $nElemT  material 9 stress
#recorder Element -file Gstrain1.out   -time  -eleRange 1 $nElemT  material 1 strain
#recorder Element -file Gstrain2.out   -time  -eleRange 1 $nElemT  material 2 strain
#recorder Element -file Gstrain3.out   -time  -eleRange 1 $nElemT  material 3 strain
#recorder Element -file Gstrain4.out   -time  -eleRange 1 $nElemT  material 4 strain
#recorder Element -file Gstrain9.out   -time  -eleRange 1 $nElemT  material 9 strain
puts "Finished creating gravity recorders..."

#-----------------------------------------------------------------------------------------
#  9. DEFINE ANALYSIS PARAMETERS
#-----------------------------------------------------------------------------------------

#---GROUND MOTION PARAMETERS

# Number of steps in ground motion record. 
# Retrieved from GM file name (exterior input)
set nSteps $npts
puts "duration = [expr $nSteps*$dt] sec"

#---RAYLEIGH DAMPING PARAMETERS
set pi      3.141592654
# Damping ratio
#set damp    0.05
# Lower frequency
#set omega1  [expr 2*$pi*f_site]
#set omega1  [expr 2*$pi*0.2]
# Upper frequency
#set omega2  [expr 2*$pi*5*f_site]
#set omega2  [expr 2*$pi*20]
# Damping coefficients
set a0      [expr 2*$damp*$omega1*$omega2/($omega1 + $omega2)]
set a1      [expr 2*$damp/($omega1 + $omega2)]
puts "damping coefficients: a_0 = $a0;  a_1 = $a1"

#---DETERMINE STABLE ANALYSIS TIME STEP USING CFL CONDITION
#STABILITY CHECK NOT REQUIRED FOR IMPLICIT (NEWMARK) ANALYSIS
# maximum shear wave velocity (m/s)
#set vsMax $Vs(0)
#for {set i 1} {$i <= $numLayers} {incr i 1} {
#    if {$Vs($i) > $vsMax} {
#        set vsMax $Vs($i)
#}
#}
# duration of ground motion (s)
#set duration    [expr $motionDT*$motionSteps]
# minimum element size
#set minSize $sElemY(0)
#for {set i 1} {$i <= $numLayers} {incr i 1} {
#    if {$sElemY($i) < $minSize} {
#        set minSize $sElemY($i)
#}
#}
#
# trial analysis time step
#set kTrial      [expr $minSize/(pow($vsMax,0.5))]
# define time step and number of steps for analysis
#if { $motionDT <= $kTrial } {
#    set nSteps  $motionSteps
#    set dt      $motionDT
#} else {
#    set nSteps  [expr int(floor($duration/$kTrial)+1)]
#    set dt      [expr $duration/$nSteps]
#}
puts "number of steps in analysis: $nSteps"
puts "analysis time step: $dt"

#---ANALYSIS PARAMETERS

# Newmark parameters
# Average acceleration method:
# Unconditionally stable, without numerical damping
set gamma           0.50
set beta            0.25

#-----------------------------------------------------------------------------------------
#  10. GRAVITY ANALYSIS
#-----------------------------------------------------------------------------------------

# Update materials to ensure elastic behavior
for {set k 1} {$k <= $numLayers} {incr k} {
    updateMaterialStage -material $k -stage 0
}

# Create analysis
constraints Penalty 1.e14 1.e14
test        NormDispIncr 1e-6 35 1
algorithm   KrylovNewton
numberer    RCM
system      ProfileSPD
integrator  Newmark $gamma $beta
analysis    Transient

set startT  [clock seconds]

# Elastic gravity loading
# Large time steps are used to damp out the waves generated by loading
#       numincr dt
analyze 10      5.0e2
puts "Finished with elastic gravity analysis..."

# Update materials to consider plastic behavior
for {set k 1} {$k <= $numLayers} {incr k} {
    updateMaterialStage -material $k -stage 1
}

# Plastic gravity loading
#       numincr dt
analyze 40      5.0e-2
puts "Finished with plastic gravity analysis..."

#-----------------------------------------------------------------------------------------
#  11. UPDATE ELEMENT PERMEABILITY VALUES FOR POST-GRAVITY ANALYSIS
#-----------------------------------------------------------------------------------------

# If excess pore pressure generation is not allowed (i.e., allowPWP = No),
# permeabilites are left high for dynamic analysis

if {$allowPWP == Yes} {

    # Choose base number for parameter IDs which is higer than other tags used in analysis
    set ctr 10000.0
    # Loop over elements to define parameter IDs 
    for {set i 1} {$i<=$nElemT} {incr i 1} {
        parameter [expr int($ctr+1.0)] element $i vPerm
        parameter [expr int($ctr+2.0)] element $i hPerm
        set ctr [expr $ctr+2.0]
}

# Update permeability parameters for each element using parameter IDs
set ctr 10000.0
for {set j 1} {$j <= $nElemT} {incr j 1} {

    set lowerBound 0.0
    for {set i 0} {$i <= $numLayers} {incr i 1} {

        if {$j <= $layerTopEleNum($i) && $j > $layerTopEleNum([expr $i - 1])} {
            updateParameter [expr int($ctr+1.0)] $vPerm($i)
            updateParameter [expr int($ctr+2.0)] $hPerm($i)
            puts "$j  updateParameter [expr int($ctr+1.0)] $vPerm($i)"

        }
        set lowerBound $layerBound($i)
    }
    set ctr [expr $ctr+2.0]
}
puts "Finished updating permeabilities for dynamic analysis..."
}

#-----------------------------------------------------------------------------------------
#  12. CREATE POST-GRAVITY RECORDERS
#-----------------------------------------------------------------------------------------

# Create list for pore pressure nodes
set leftPPNodeList {}
set channel1 [open "$saveDir/$siteModelID.LeftPPNodesInfo.dat" r]
set count1 0;
foreach line [split [read -nonewline $channel1] \n] {
    set count1 [expr $count1+1];
    set lineData($count1) $line
    set nodeNumber [lindex $lineData($count1) 0]
    lappend leftPPNodeList $nodeNumber
    puts "node: $nodeNumber"
}
set leftPPNodeList [lsort -integer $leftPPNodeList]
puts "leftPPNodes are: $leftPPNodeList"
close $channel1

# Reset time and analysis
setTime 0.0
wipeAnalysis
remove recorders

# Recorder time step
set recDT  [expr 10*$dt]

# Record nodal displacment, acceleration, and porepressure
#eval "recorder Node -file $saveDir/$analID.disp.out -time -dt $dt -node $leftPPNodeList -dof 1 2  disp"
#eval "recorder Node -file $saveDir/$analID.vel.out -time -dt $dt -node $leftPPNodeList -dof 1 2  vel"
#eval "recorder Node -file $saveDir/$analID${dash}acc.out -time -node $leftPPNodeList -dof 1 2 accel"
#recorder Node -file $saveDir/$analID${dash}accHorSur.out -time -node $nNodeT -dof 1 accel
# record horizontal acceleration at the top node
eval "recorder Node -file $saveDir/$analID${dash}accHorSur.out -time -node $nNodeT -dof 1 accel"

# record elemental stress and strain
#recorder Element -file $saveDir/$analID.stress1.out   -time  -eleRange 1 $nElemT  material 1 stress 3
#recorder Element -file $saveDir/$analID.stress2.out   -time  -eleRange 1 $nElemT  material 2 stress 3
#recorder Element -file $saveDir/$analID.stress3.out   -time  -eleRange 1 $nElemT  material 3 stress 3
#recorder Element -file $saveDir/$analID.stress4.out   -time  -eleRange 1 $nElemT  material 4 stress 3
#recorder Element -file $saveDir/$analID.stress5.out   -time  -eleRange 1 $nElemT  material 5 stress 3
#recorder Element -file $saveDir/$analID.stress6.out   -time  -eleRange 1 $nElemT  material 6 stress 3
#recorder Element -file $saveDir/$analID.stress7.out   -time  -eleRange 1 $nElemT  material 7 stress 3
#recorder Element -file $saveDir/$analID.stress8.out   -time  -eleRange 1 $nElemT  material 8 stress 3
#recorder Element -file $saveDir/$analID${dash}stress9.out   -time  -eleRange 1 $nElemT  material 9 stress 3
#recorder Element -file $saveDir/$analID.strain1.out   -time  -eleRange 1 $nElemT  material 1 strain
#recorder Element -file $saveDir/$analID.strain2.out   -time  -eleRange 1 $nElemT  material 2 strain
#recorder Element -file $saveDir/$analID.strain3.out   -time  -eleRange 1 $nElemT  material 3 strain
#recorder Element -file $saveDir/$analID.strain4.out   -time  -eleRange 1 $nElemT  material 4 strain
#recorder Element -file $saveDir/$analID.strain5.out   -time  -eleRange 1 $nElemT  material 5 strain
#recorder Element -file $saveDir/$analID.strain6.out   -time  -eleRange 1 $nElemT  material 6 strain
#recorder Element -file $saveDir/$analID.strain7.out   -time  -eleRange 1 $nElemT  material 7 strain
#recorder Element -file $saveDir/$analID.strain8.out   -time  -eleRange 1 $nElemT  material 8 strain
#recorder Element -file $saveDir/$analID${dash}strain9.out   -time  -eleRange 1 $nElemT  material 9 strain
puts "Finished creating all recorders..."


#-----------------------------------------------------------------------------------------
#  13. DYNAMIC ANALYSIS
#-----------------------------------------------------------------------------------------

# Create a 2D model with 3 degrees of freedom
model BasicBuilder -ndm 2 -ndf 3

# Define constant scaling factor for applied velocity
set cFactor [expr $colArea*$dashpotCoeff]

# Define velocity time history file
set velocityFile $gMotionPath

# Timeseries object for force history
set mSeries "Path -dt $dt -filePath $velocityFile -factor $cFactor"

# Loading object
#             patternTag timeSeries
pattern Plain 13         $mSeries {
#        nodeTag dof1 dof2 dof3
    load 1       1.0  0.0  0.0
}
puts "Dynamic loading created..."

# Create analysis
constraints Penalty 1.e16 1.e16
test        NormDispIncr 1.0e-5 35 1
algorithm   KrylovNewton
numberer    RCM
system      ProfileSPD
integrator  Newmark $gamma $beta
rayleigh    $a0 $a1 0.0 0.0
analysis    Transient

# Perform analysis with timestep reduction loop
#               numincr	dt
set ok [analyze $nSteps $dt]

# If analysis fails, reduce timestep and continue with analysis
if {$ok != 0} {
    puts "did not converge, reducing time step"
    set curTime  [getTime]
    set mTime $curTime
    puts "curTime: $curTime"
    set curStep  [expr $curTime/$dt]
    puts "curStep: $curStep"
    set rStep  [expr ($nSteps-$curStep)*2.0]
    set remStep  [expr int(($nSteps-$curStep)*2.0)]
    puts "remStep: $remStep"
    set dt       [expr $dt/2.0]
    puts "dt: $dt"
    puts "Ground Motion is: $gMotionName"

    set ok [analyze  $remStep  $dt]

    # If analysis fails again, reduce timestep and continue with analysis
    if {$ok != 0} {
        puts "did not converge, reducing time step"
        set curTime  [getTime]
        puts "curTime: $curTime"
        set curStep  [expr ($curTime-$mTime)/$dt]
        puts "curStep: $curStep"
        set remStep  [expr int(($rStep-$curStep)*2.0)]
        puts "remStep: $remStep"
        set dt       [expr $dt/2.0]
        puts "dt: $dt"
        puts "Ground Motion is: $gMotionName"

        analyze  $remStep  $dt
}
}
set endT    [clock seconds]
puts "Finished with dynamic analysis..."
puts "Analysis execution time: [expr $endT-$startT] seconds"

wipe
