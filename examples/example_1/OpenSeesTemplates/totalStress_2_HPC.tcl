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
