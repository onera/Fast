import Converter.PyTree as CP
import Converter.Internal as CI
import numpy
import os
import sys
import math
from collections import OrderedDict

#from userData import *

def addInletData(data,famName):
   inletDict = {}
   inletDict['type']                = 'inj1'
   inletDict['stagnation_pressure'] = Pio
   inletDict['stagnation_enthalpy'] = Hio
   inletDict['txv']                 = math.cos(alphaAm*pi/180.)
   inletDict['tyv']                 = math.sin(alphaAm*pi/180.)
   inletDict['tzv']                 = 0.0
   inletDict['inj_tur1'] = k1_ini
   if turbmod not in ['SA']: inletDict['inj_tur2'] = k2_ini
   nodeInlet = CI.getNodeFromName(data,famName)
   childInletL = []
   for k in inletDict.keys(): childInletL.append(CI.createNode(k, 'DataArray_t', value=inletDict[k], children=[]))
   childinlet = CI.createChild(nodeInlet,'.Solver#BC','UserDefinedData_t', value=None, children=childInletL)

def addOutletData(data,famName):
   outletDict = {}
   outletDict['type']     = 'outpres'
   outletDict['pressure'] = Ps_aval
   nodeOutlet = CI.getNodeFromName(data,famName)
   childOutletL = []
   for k in outletDict.keys(): childOutletL.append(CI.createNode(k, 'DataArray_t', value=outletDict[k], children=[]))
   childOutlet = CI.createChild(nodeOutlet,'.Solver#BC','UserDefinedData_t', value=None, children=childOutletL)

def addOmega(data,famName,value):
   nodeOmg  = CI.getNodeFromName(data,famName)
   childOmg = CI.createNode('omega', 'DataArray_t', value=value, children=[])
   childOmg = CI.createChild(nodeOmg,'.Solver#BC','UserDefinedData_t', value=None, children=[childOmg])

def addOutradeq(data,famName,prespiv):
   nodeOutradeq = CI.getNodeFromName(data,famName)
   outradeqDict = {}
   outradeqDict['type'    ]= 'outradeq'
   outradeqDict['indpiv'  ]= 1
   outradeqDict['dirorder']= -1
   outradeqDict['prespiv' ]= prespiv
   childOutradeqL = []
   for k in outradeqDict.keys(): childOutradeqL.append(CI.createNode(k, 'DataArray_t', value=outradeqDict[k], children=[]))
   childOutradeq = CI.createChild(nodeOutradeq,'.Solver#BC','UserDefinedData_t', value=None, children=childOutradeqL)

def addMxpl(data,famName1av, famName2am):
   spDict = {}
   spDict['type']            = 'stage_mxpl'
   spDict['mxpl_avermean']   = 'riemann'
   spDict['mxpl_avertur']    = 'conservative'
   spDict['mxpl_num']        = 'characteristic'
   spDict['jtype']           = 'nomatch_rad_line'
   spDict['nomatch_special'] = 'none'
   spDict['globborder']      = famName1av
   spDict['globborderdonor'] = famName2am
   spDict['glob_dir_i']      =  2
   spDict['glob_dir_j']      =  3
   data = addStageBC(data,famName1av, famName2am, spDict)
   return data

def addRNA(data,famName1av, famName2am, TrefAdim):
   spDict = {}
   spDict['type']            = 'stage_red'
   spDict['jtype']           = 'nomatch_rad_line'
   spDict['stage_red_type']  = 'half_sum'
   spDict['stage_ref_time']  = TrefAdim
   spDict['nomatch_special'] = 'none'
   spDict['globborder']      = famName1av
   spDict['globborderdonor'] = famName2am
   spDict['glob_dir_i']      =  2
   spDict['glob_dir_j']      =  3
   data = addStageBC(data,famName1av, famName2am, spDict)
   return data

def addStageBC(data,famName1av, famName2am, spDict):
   # - recuperer noeud dont famille == famName1av
   (roue1av,c) = CI.getParentOfNode(data,CI.getNodesFromValue(data,famName1av)[0])
   #(roue1av,c) = CI.getParentOfNode(data,CI.getNodeFromName(CI.getNodesFromValue(data,famName1av), 'FamilyName')) # correction Alan
   roue1avPtRange = CI.getNodeFromName(roue1av,'PointRange')
   # - recuperer noeud dont famille == famName2am
   (roue2am,c) = CI.getParentOfNode(data,CI.getNodesFromValue(data,famName2am)[0])
   roue2amPtRange = CI.getNodeFromName(roue2am,'PointRange')
   # - s'echanger PointRange => PointRangeDonor
   roue1avPRdonor = CI.copyNode(roue2amPtRange)
   CI.setName(roue1avPRdonor,'PointRangeDonor')
   CI.addChild(roue1av, roue1avPRdonor)
   roue2amPRdonor = CI.copyNode(roue1avPtRange)
   CI.setName(roue2amPRdonor,'PointRangeDonor')
   CI.addChild(roue2am, roue2amPRdonor)
   # - ajouter noeud Transform (valeur = array(1,2,3))
   transform = CI.createNode('Transform', 'DataArray_t', value=numpy.array([1, 2, 3], numpy.int32), children=[])
   CI.addChild(roue1av, transform)
   CI.addChild(roue2am, transform)
   # - ajouter noeud .Solver#Property et infos mxpl sous .SolverProperty
   #   + roue1av
   spChildL = []
   for k in spDict.keys(): spChildL.append(CI.createNode(k, 'DataArray_t', value=spDict[k], children=[]))
   roue1avSp = CI.createChild(roue1av,'.Solver#Property','UserDefinedData_t', value=None, children=spChildL)
   #   + roue2am
   spChildL = []
   tmp = spDict['globborder']
   spDict['globborder']      = spDict['globborderdonor']
   spDict['globborderdonor'] = tmp
   for k in spDict.keys(): spChildL.append(CI.createNode(k, 'DataArray_t', value=spDict[k], children=[]))
   roue2amSp = CI.createChild(roue2am,'.Solver#Property','UserDefinedData_t', value=None, children=spChildL)
   # - bouger noeud vers ZoneGridConnectivity + changer type voire valeur
   #   + recuperation nom zone_t roue1av
   roue1avPath = CI.getPath(data,roue1av)
   roue1avPL   = roue1avPath.split('/')
   roue1avBlk  = getParentNodeWithType(data, 'Zone_t', roue1avPath)
   roue1avName = CI.getName(roue1avBlk)
   #   + recuperation nom zone_t roue2am
   roue2amPath = CI.getPath(data,roue2am)
   roue2amPL   = roue2amPath.split('/')
   roue2amBlk  = getParentNodeWithType(data, 'Zone_t', roue2amPath)
   roue2amName = CI.getName(roue2amBlk)
   #   + modification valeur
   CI.setValue(roue1av,roue2amName)
   CI.setValue(roue2am,roue1avName)
   #   + ajout dans liste child de GridConnectivity
   roue1avGC = CI.getChildren(CI.getByType(roue1avBlk, 'ZoneGridConnectivity_t'))[0]
   CI.addChild(roue1avGC, roue1av)
   roue2amGC = CI.getChildren(CI.getByType(roue2amBlk, 'ZoneGridConnectivity_t'))[0]
   CI.addChild(roue2amGC, roue2am)
   #   + modification du type
   CI.setType(roue1av, 'GridConnectivity_t')
   CI.setType(roue2am, 'GridConnectivity_t')
   #   + suppression dans noeud BC
   data = CI.rmNodeByPath(data,roue1avPath)
   data = CI.rmNodeByPath(data,roue2amPath)
   return data

def getParentNodeWithType(data, typeName, path):
   pathL = path.split('/')
   test  = False
   ind   = -1
   while not test:
      node = CI.getNodeFromName(data,pathL[ind])
      test = CI.isType(node, 'Zone_t')
      #print test,ind, CI.getName(node)
      ind-=1
   return node

def addTriggerAndOutput(data, fam):
   stDict = {}
   stDict['file']           = 'script_coupling.py'
   stDict['next_iteration'] = 1
   stDict['next_state']     = 16
   stChildL = []
   for k in stDict.keys(): stChildL.append(CI.createNode(k, 'UserDefinedData_t', value=stDict[k], children=[]))
   soDict = {}
   #soDict['var']           = 'PressureStagnation EnthalpyStagnation'
   soDict['var']           = 'PressureStagnation'
   soDict['period']        = 1
   soDict['writingmode']   = 2
   soDict['loc']           = 4
   soChildL = []
   for k in soDict.keys(): soChildL.append(CI.createNode(k, 'DataArray_t', value=soDict[k], children=[]))
   # ajout Trigger dans famille
   for inlet in CI.getNodesFromName(data,fam):
     childInletTrigger = CI.createChild(inlet,'.Solver#Trigger','UserDefinedData_t', value=None, children=stChildL)
     childInletOutput = CI.createChild(inlet,'.Solver#Output','UserDefinedData_t', value=None, children=soChildL)
   return data

# =======================
# Add-ons FastS
# =======================
def addFlowEquationSet2FastS(data):
   bases=CI.getByType(data,'CGNSBase_t')[2]

   for base in bases:
      Prop=CI.createChild(base,'FlowEquationSet','FlowEquationSet_t')
      CI.createChild(Prop,'EquationDimension','"int"',value=Equation_Dimension)
      CI.createChild(Prop,'GoverningEquations','GoverningEquations_t',value=Governing_Equation)
      CI.createChild(Prop,'TurbulenceModel','TurbulenceModel_t',value=Turbulence_Model)
      
def addFlowSolution2FastS(data, nbBlock, rootDataFileName):
	
	data2setList   = ['Density', 'MomentumX', 'MomentumY', 'MomentumZ', 'EnergyStagnationDensity', 'TurbulentEnergyKineticDensity']
	nbData2set = len(data2setList)

#~ ==============Functions=============
	
	def getMeshBlocks(tree):
		blockMesh = []
		base   = CI.getNodesFromType(tree,'CGNSBase_t')
		baseN  = CI.getName(base)
		
		zonesList = []
		zones = CI.getChildren(baseN)
		for zone in zones:             #we get a list of all the children nodes' name of the base
			name = CI.getName(zone)
			zonesList.append(name)
			
		nbZones = len(zonesList)
		for p in range (nbZones):      # we only keep the mesh blocks' zone
			zone = zonesList[p]
			for q in range(1, nbBlock+1):
				if str(q) in zone:
					blockMesh.append(zone)
		return blockMesh	
	
	def getDataFiles(nbBlock, rootFileName):
		dataFiles = []
		
		for i in range(nbBlock):
			fileName  = rootFileName + str(i+1)
			dataFiles.append(fileName)
		
		return dataFiles

	def getIndexOfBlock(tree, block):
		node  = CI.getNodeFromName(tree, block)
		index = CI.getValue(node)
		
		im = index[0][0]
		jm = index[1][0]
		km = index[2][0]
		
		return im, jm, km
		
	def getArrayAndReshapeF(tree, data2extract):                       
			node  = CI.getNodeFromName(tree, data2extract )                #we get the node
			array = CI.getValue(node)                                   #we get the array associated 
			
			#~ size = array.size
			#~ array = numpy.reshape(array, size, order='F')               #we reshape the array with a Fortran like order		
			return array
		
	def getFlowSolNode(tree, block):
		base        = CI.getNodesFromType(tree,'CGNSBase_t')
		baseN       = CI.getName(base)
		nodePath    =  '/' + CI.getName(baseN) + '/' + block + '/FlowSolution#Centers/'
		
		flowSolNode = CI.getNodeFromPath(tree, nodePath)
		
		return flowSolNode
			
	def getFlowSolChildrenDict(tree, flowSolNode):

		childrenNodeList = CI.getChildren(flowSolNode) #we get all the node of FlowSolution#Centers to keep what is needed and replace with the cartos 
		nbChildren       = len(childrenNodeList)
		childrenDict     = OrderedDict()               #we create a dictionary with key: the name of the node, value: the node
		
		for p in range(nbChildren):        
			childNode                = childrenNodeList[p]
			childnName               = CI.getName(childNode)
			childrenDict[childnName] = childNode			

		return childrenDict		
		
#~ ================Scripts================
	blockMeshList = getMeshBlocks(data)   						#we get a list of all the blocks in the input file
	dataFilesList = getDataFiles(nbBlock,rootDataFileName) 		#we get a list of all the files used to set the FlowSolution node (restart files)
	
	
	for p in range(nbBlock):                                     #for each block (block by block)
		blockMesh = blockMeshList[p]
		dataFile = dataFilesList[p]                              #we get the right file
		tdata2set = CP.convertFile2PyTree(dataFile,'bin_tp')     #that we transform into a tree
			
		flowSolNode = getFlowSolNode(data , blockMesh)           #we get the corresponding FlowSolution node
		chilDict    = getFlowSolChildrenDict(data, flowSolNode)  #we create a dict with = key : the children nodes' name/ value = the corresponding node

		
		for q in range(nbData2set):								 #for data to set
			data2set = data2setList[q]
			data2setArray = getArrayAndReshapeF(tdata2set, data2set) #we extract the corresponding array
			print 'shape array:', data2set, "value ====", data2setArray.shape
			
			if data2set in chilDict :                            #if the data already exists in the flowsolution's node: we change it's value by the new one
				CI.setValue(chilDict[data2set], data2setArray)
					
			else:                                                #else: create a new child with the right value
				newNode = CI.createChild(flowSolNode, data2set, 'DataArray_t', value=data2setArray)		
					      


def convertRelative2Absolute(data, nbBlock, omega, gamma ):
	#~ ==============Functions=============
	
	def getMeshBlocks(tree):
		blockMesh = []
		base   = CI.getNodesFromType(tree,'CGNSBase_t')
		baseN  = CI.getName(base)
		
		zonesList = []
		zones = CI.getChildren(baseN)
		for zone in zones:             #we get a list of all the children nodes' name of the base
			name = CI.getName(zone)
			zonesList.append(name)
			
		nbZones = len(zonesList)
		for p in range (nbZones):      # we only keep the mesh blocks' zone
			zone = zonesList[p]
			for q in range(1, nbBlock+1):
				if str(q) in zone:
					blockMesh.append(zone)
		return blockMesh	
	
	def getFlowSolNode(tree, block):
		base        = CI.getNodesFromType(tree,'CGNSBase_t')
		baseN       = CI.getName(base)
		nodePath    =  '/' + CI.getName(baseN) + '/' + block + '/FlowSolution#Centers/'
		
		flowSolNode = CI.getNodeFromPath(tree, nodePath)
		
		return flowSolNode


# Definition of the functions needed to convert our data
	def norm3D(x,y,z)                                 : return math.sqrt(x**2 + y**2 + z**2)
	def atan2(x,y)                                    : return math.atan2(x,y)
	
	def Ps(gamma, Density, E, Urel, Vrel, Wrel)       : return ((gamma-1)*Density)*(E-0.5*norm3D(Urel, Vrel, Wrel))
	def roE_abs(Ps, gamma, Density, Uabs, Vabs, Wabs) : return (Ps/(gamma-1))+0.5*Density*norm3D(Uabs, Vabs, Wabs)
	
	def Vt_rel(Vrel, Wrel, teta)                      : return Vrel*math.cos(teta) - Wrel*math.sin(teta)
	def Vr_rel(Vrel, Wrel, teta)                      : return Vrel*math.sin(teta)+ Wrel*math.cos(teta)
	
	def Vt_abs(Vt_rel, omg, CoordinateY, CoordinateZ)   : return Vt_rel - omg*norm3D(CoordinateY, CoordinateZ, 0)
	
	def Vabs(Vr_abs, Vt_abs, teta)                    : return Vr_abs*math.sin(teta) + Vt_abs*math.sin(teta)
	def Wabs(Vr_abs, Vt_abs, teta)                    : return Vr_abs*math.cos(teta) - Vt_abs*math.cos(teta)
	
#~ ================Scripts================
	blockMeshList = getMeshBlocks(data)   						#we get a list of all the blocks in the input file

	for p in range(nbBlock):                                     #for each block (block by block)
		blockMesh = blockMeshList[p]

		zoneInit = CI.getNodeFromName(data, blockMesh)
		zone = CI.getNodeFromName(data, blockMesh)     #we extract the mesh zone
		zone = CP.node2Center(zone)                  #we change the coordinates (ONLY) to have them at the center of the cell
		
		#we add what is needed to convert our data
		zone = CP.initVars(zone, 'omg', omega)       
		zone = CP.initVars(zone, 'gamma', gamma)
		zone = CP.initVars(zone, 'teta', atan2, ['CoordinateY', 'CoordinateZ'])
		
		#relative velocity
		zone = CP.initVars(zone, '{Urel}={MomentumX}/{Density}')   
		zone = CP.initVars(zone, '{Vrel}={MomentumY}/{Density}')	
		zone = CP.initVars(zone, '{Wrel}={MomentumZ}/{Density}')
		
	
		zone = CP.initVars(zone, '{E}={EnergyStagnationDensity}/{Density}')
		zone = CP.initVars(zone, 'Ps', Ps, ['gamma', 'Density', 'E', 'Urel', 'Vrel', 'Wrel'])
		
		#relative velocity components
		zone = CP.initVars(zone, 'Vt_rel', Vt_rel, ['Vrel', 'Wrel', 'teta']) 
		zone = CP.initVars(zone, 'Vr_rel', Vr_rel, ['Vrel', 'Wrel', 'teta'])
		
		#absolute velocity components
		zone = CP.initVars(zone, 'Vt_abs', Vt_abs , ['Vt_rel', 'omg', 'CoordinateY', 'CoordinateZ']) 
		zone = CP.initVars(zone, '{Vr_abs}={Vr_rel}')
		
		#absolute velocity
		zone = CP.initVars(zone, '{Uabs}={Urel}')
		zone = CP.initVars(zone, 'Vabs', Vabs, ['Vr_abs', 'Vt_abs', 'teta']) 
		zone = CP.initVars(zone, 'Wabs', Wabs, ['Vr_abs', 'Vt_abs', 'teta'])
		
		#reconstruction of absolute conservative variables
		zone = CP.initVars(zone, '{MomentumX}={Uabs}*{Density}') 
		zone = CP.initVars(zone, '{MomentumY}={Vabs}*{Density}') 
		zone = CP.initVars(zone, '{MomentumZ}={Wabs}*{Density}')
		zone = CP.initVars(zone, 'EnergyStagnationDensity', roE_abs, ['Ps', 'gamma', 'Density', 'Uabs', 'Vabs', 'Wabs'])
		
		#we get the FlowSolution node to set the data
		zone = CI.renameNode(zone, 'FlowSolution', 'FlowSolution#Centers') 
		flowSol = getFlowSolNode(data, blockMesh)
		
		print 'Processing relative to absolute velocity ===', blockMesh
		ConsVariables    = ['MomentumX', 'MomentumY', 'MomentumZ', 'EnergyStagnationDensity']
		
		for i in range(len(ConsVariables)): 
			ConsVariable = ConsVariables[i]
			
			#we delete the previous nodes
			relativeNode   = CI.getNodeFromName(flowSol, ConsVariable)
			CI.rmNode(flowSol, relativeNode) 
			#we add what we want on the node
			absoluteNode  = CI.getNodeFromName(zone, ConsVariable)
			CI.addChild(flowSol, absoluteNode) 


def convertRelative2Absolute2(data, nbBlock, omega, gamma ):
	#~ ==============Functions=============
	
	def getMeshBlocks(tree):
		blockMesh = []
		base   = CI.getNodesFromType(tree,'CGNSBase_t')
		baseN  = CI.getName(base)
		
		zonesList = []
		zones = CI.getChildren(baseN)
		for zone in zones:             #we get a list of all the children nodes' name of the base
			name = CI.getName(zone)
			zonesList.append(name)
			
		nbZones = len(zonesList)
		for p in range (nbZones):      # we only keep the mesh blocks' zone
			zone = zonesList[p]
			for q in range(1, nbBlock+1):
				if str(q) in zone:
					blockMesh.append(zone)
		return blockMesh	
	
	def getFlowSolNode(tree, block):
		base        = CI.getNodesFromType(tree,'CGNSBase_t')
		baseN       = CI.getName(base)
		nodePath    =  '/' + CI.getName(baseN) + '/' + block + '/FlowSolution#Centers/'
		
		flowSolNode = CI.getNodeFromPath(tree, nodePath)
		
		return flowSolNode


# Definition of the functions needed to convert our data
	def norm3D(x,y,z)                                 : return math.sqrt(x**2 + y**2 + z**2)
	def atan2(x,y)                                    : return math.atan2(x,y)
	
	def Ps(gamma, Density, E, Urel, Vrel, Wrel)       : return ((gamma-1)*Density)*(E-0.5*norm3D(Urel, Vrel, Wrel))
	def roE_abs(Ps, gamma, Density, Uabs, Vabs, Wabs) : return (Ps/(gamma-1))+0.5*Density*norm3D(Uabs, Vabs, Wabs)
	
	def Vt_rel(Vrel, Wrel, teta)                      : return Vrel*math.cos(teta) - Wrel*math.sin(teta)
	def Vr_rel(Vrel, Wrel, teta)                      : return Vrel*math.sin(teta)+ Wrel*math.cos(teta)
	
	def Vt_abs(Vt_rel, omg, CoordinateY, CoordinateZ)   : return Vt_rel - omg*norm3D(CoordinateY, CoordinateZ, 0)
	
	def Vabs(Vr_abs, Vt_abs, teta)                    : return Vr_abs*math.sin(teta) + Vt_abs*math.sin(teta)
	def Wabs(Vr_abs, Vt_abs, teta)                    : return Vr_abs*math.cos(teta) - Vt_abs*math.cos(teta)
	
	def Ts(Ps, density)                               : return Ps/(density*287.04)
	
#~ ================Scripts================
	blockMeshList = getMeshBlocks(data)   						#we get a list of all the blocks in the input file

	for p in range(nbBlock):                                     #for each block (block by block)
		blockMesh = blockMeshList[p]

		zoneInit = CI.getNodeFromName(data, blockMesh)
		zone = CI.getNodeFromName(data, blockMesh)     #we extract the mesh zone
		zone = CP.node2Center(zone)                  #we change the coordinates (ONLY) to have them at the center of the cell
		
		#we add what is needed to convert our data
		zone = CP.initVars(zone, 'omg', omega)       
		zone = CP.initVars(zone, 'R', 287.04)       
		zone = CP.initVars(zone, 'gamma', gamma)
		zone = CP.initVars(zone, 'teta', atan2, ['CoordinateY', 'CoordinateZ'])
		
		#relative velocity
		zone = CP.initVars(zone, '{Urel}={MomentumX}/{Density}')   
		zone = CP.initVars(zone, '{Vrel}={MomentumY}/{Density}')	
		zone = CP.initVars(zone, '{Wrel}={MomentumZ}/{Density}')
		
	
		zone = CP.initVars(zone, '{E}={EnergyStagnationDensity}/{Density}')
		zone = CP.initVars(zone, 'Ps', Ps, ['gamma', 'Density', 'E', 'Urel', 'Vrel', 'Wrel'])
		
		#relative velocity components
		zone = CP.initVars(zone, 'Vt_rel', Vt_rel, ['Vrel', 'Wrel', 'teta']) 
		zone = CP.initVars(zone, 'Vr_rel', Vr_rel, ['Vrel', 'Wrel', 'teta'])
		
		#absolute velocity components
		zone = CP.initVars(zone, 'Vt_abs', Vt_abs , ['Vt_rel', 'omg', 'CoordinateY', 'CoordinateZ']) 
		zone = CP.initVars(zone, '{Vr_abs}={Vr_rel}')
		
		#absolute velocity
		zone = CP.initVars(zone, '{VelocityX}={Urel}')
		zone = CP.initVars(zone, 'VelocityY', Vabs, ['Vr_abs', 'Vt_abs', 'teta']) 
		zone = CP.initVars(zone, 'VelocityZ', Wabs, ['Vr_abs', 'Vt_abs', 'teta'])
		
		#reconstruction of absolute conservative variables
		zone = CP.initVars(zone, 'EnergyStagnationDensity', roE_abs, ['Ps', 'gamma', 'Density', 'VelocityX', 'VelocityY', 'VelocityZ'])
		zone = CP.initVars(zone, '{EnergyStagnation} = {EnergyStagnationDensity}/{Density}')
		
		zone = CP.initVars(zone, '{TurbulentSANuTilde} = {TurbulentSANuTildeDensity}/{Density}')
		
		zone = CP.initVars(zone, '{TurbulentEnergyKinetic} = {TurbulentEnergyKineticDensity}/{Density}')
		
		
		zone = CP.initVars(zone, 'Temperature', Ts, ['Ps', 'Density'])
		
		#~ ================================
		#~ zone = CP.initVars(zone, '{VelocityX_M1}={VelocityX}')
		#~ zone = CP.initVars(zone, '{VelocityY_M1}={VelocityY}')
		#~ zone = CP.initVars(zone, '{VelocityZ_M1}={VelocityZ}')
		#~ 
		#~ zone = CP.initVars(zone, '{Density_M1}={Density}')
		#~ zone = CP.initVars(zone, '{Temperature_M1}={Temperature}')
		#~ zone = CP.initVars(zone, '{TurbulentSANuTilde_M1}={TurbulentSANuTilde}')
	        #~ 
		#~ 
		#~ zone = CP.initVars(zone, '{VelocityX_P1}={VelocityX_M1}')
		#~ zone = CP.initVars(zone, '{VelocityY_P1}={VelocityY_M1}')
		#~ zone = CP.initVars(zone, '{VelocityZ_P1}={VelocityZ_M1}')
		#~ 
		#~ zone = CP.initVars(zone, '{Density_P1}={Density_M1}')
		#~ zone = CP.initVars(zone, '{Temperature_P1}={Temperature_M1}')
		#~ zone = CP.initVars(zone, '{TurbulentSANuTilde_P1}={TurbulentSANuTilde_M1}')
		
		#~ ================================
		
		#we get the FlowSolution node to set the data
		zone = CI.renameNode(zone, 'FlowSolution', 'FlowSolution#Centers') 
		flowSol = getFlowSolNode(data, blockMesh)
		
		print 'Processing===', blockMesh
		ConsVariables    = ['MomentumX', 'MomentumY', 'MomentumZ', 'EnergyStagnationDensity', 'TurbulentSANuTildeDensity', 'TurbulentEnergyKineticDensity']
		
		PrimVariables_0  =  ['VelocityX', 'VelocityY', 'VelocityZ', 'Temperature', 'TurbulentSANuTilde', 'TurbulentEnergyKinetic']
		#~ PrimVariables_M1 = ['Density_M1','VelocityX_M1', 'VelocityY_M1', 'VelocityZ_M1', 'Temperature_M1', 'TurbulentSANuTilde_M1']
		#~ PrimVariables_P1 = ['Density_P1', 'VelocityX_P1', 'VelocityY_P1', 'VelocityZ_P1', 'Temperature_P1', 'TurbulentSANuTilde_P1']
		#~ PrimVariables    = PrimVariables_0 + PrimVariables_M1 + PrimVariables_P1
		PrimVariables    = PrimVariables_0 

		
		for i in range(len(ConsVariables)): 
			#we delete the previous nodes
			ConsVariable = ConsVariables[i]
			relativeNode   = CI.getNodeFromName(flowSol, ConsVariable)
			CI.rmNode(flowSol, relativeNode)
			
		for j in range(len(PrimVariables)):
			#we add what we want on the node
			PrimVariable = PrimVariables[j]
			absoluteNode  = CI.getNodeFromName(zone, PrimVariable)
			CI.addChild(flowSol, absoluteNode)
						    
def addReferenceState2FastS(data):
    VelocityX                 = rou_ini/ro_ini
    VelocityY                 = rov_ini/ro_ini
    VelocityZ                 = row_ini/ro_ini
    Density                   = ro_ini
    MomentumX                 = rou_ini
    MomentumY                 = rov_ini
    MomentumZ                 = row_ini
    EnergyStagnationDensity   = roe_ini
    Pressure                  = ro_ini*R_gaz*Tso
    Temperature               = Tso
    Cv                        = cv
    Gamma                     = Gam
    Rok                       = k1_ini
    RoOmega                   = k2_ini
    TurbulentSANuTildeDensity = k1_ini
    Mus                       = musuth
    Cs                        = csuth
    Ts                        = Tsuth
    Pr                        = Prandtl
    Reynolds                  = 4.*VelocityX*Density/musuth
    
    dict_reference_state = OrderedDict()
    dict_reference_state['Density']=                    Density,
    dict_reference_state['MomentumX']=                  MomentumX,
    dict_reference_state['MomentumY']=                  MomentumY,
    dict_reference_state['MomentumZ']=                  MomentumZ,
    dict_reference_state['EnergyStagnationDensity']=    EnergyStagnationDensity,  
    dict_reference_state['Pressure']=                   Pressure,
    dict_reference_state['Temperature']=                Temperature,
    dict_reference_state['Cv']=                         Cv,
    dict_reference_state['Mach']=                       Mo,
    dict_reference_state['Reynolds']=                   Reynolds,
    dict_reference_state['TurbulentSANuTildeDensity']=  TurbulentSANuTildeDensity,
    dict_reference_state['Gamma']=                      Gamma,
    dict_reference_state['VelocityX']=                  VelocityX,
    dict_reference_state['VelocityY']=                  VelocityY,
    dict_reference_state['VelocityZ']=                  VelocityZ,
    dict_reference_state['Mus']=                        Mus,
    dict_reference_state['Cs']=                         Cs,
    dict_reference_state['Ts']=                         Ts,
    dict_reference_state['Pr']=                         Pr,
    dict_reference_state['Rok']=                        Rok,
    dict_reference_state['RoOmega']=                    RoOmega
    
    bases=CI.getByType(data,'CGNSBase_t')[2]
    for base in bases:
       Prop=CI.createChild(base,'ReferenceState','ReferenceState_t')
       for d in dict_reference_state: CI.createChild(Prop,d,'DataArray_t',value=dict_reference_state[d],pos=-1)
      
def addMobileCoef2FastS(data, famNameHub, famNameCasing):
    dict_wall_hub    = {'mobile_coef':1.0}
    dict_wall_casing = {'mobile_coef':0.0}
    
    for bc in CI.getNodesFromName(data,famNameHub):
	Prop=CI.createChild(bc,'.Solver#Property','UserDefinedData_t')
	for d in dict_wall_hub: CI.createChild(Prop,d,'DataArray_t',value=dict_wall_hub[d])
	
    for bc in CI.getNodesFromName(data,famNameCasing):
	Prop=CI.createChild(bc,'.Solver#Property','UserDefinedData_t')
	for d in dict_wall_casing: CI.createChild(Prop,d,'DataArray_t',value=dict_wall_casing[d])
	
def addOutletData2FastS(data,Name):
    dict_press = OrderedDict()
    dict_press['pressure'] = Ps_aval
    dict_press['k_piv']    = kpiv
    
    for bc in CI.getNodesFromName(data,Name):
       CI.setValue(bc,'BCOutpres')
       Prop=CI.createChild(bc,'.Solver#Property','UserDefinedData_t')
       for d in dict_press: CI.createChild(Prop,d,'DataArray_t',value=dict_press[d])
    
def addUniformOutletData2FastS(data,Name,Pressure):
    dict_pres = OrderedDict()
    
    for bc in CI.getNodesFromName(data,Name):
       CI.setValue(bc,'BCOutpres')
       pr = CI.getNodesFromName(bc,'PointRange')
       index = CI.getValue(pr[0])
       i1 = index[0][0]
       im = index[0][1]
       j1 = index[1][0]
       jm = index[1][1]
       k1 = index[2][0]
       km = index[2][1]
       
       pres = []
       
       # Important : Ordre des boucles = ordre FastS
       if (i1 == im):
	   for k in range(km-1) :     
	     for j in range(jm-1) :
	       pres.append(Pressure)
       if (j1 == jm):
	   for k in range(km-1) :     
	     for i in range(im-1) :
	       pres.append(Pressure)
       if (k1 == km):
	   for j in range(jm-1) :     
	     for i in range(im-1) :
	       pres.append(Pressure)
	   
       dict_pres['pressure'] = numpy.array(pres)
       dict_pres['k_piv']  = kpiv
    
       Prop=CI.createChild(bc,'.Solver#Property','UserDefinedData_t')
       for d in dict_pres: CI.createChild(Prop,d,'DataArray_t',value=dict_pres[d])

def addNonUniformOutletData2FastS(data,Name,DataFile):
   kpiv = 1
   data_carto = CP.convertFile2PyTree(DataFile,'fmt_tp')
   data2extract = ['p']
   
   for bc in CI.getNodesFromName(data,Name):

       CI.setValue(bc,'BCOutpres')
       pr = CI.getNodesFromName(bc,'PointRange')
       index = CI.getValue(pr[0])
       
       i1 = index[0][0]
       im = index[0][1]
       j1 = index[1][0]
       jm = index[1][1]
       k1 = index[2][0]
       km = index[2][1]
   
       Prop=CI.createChild(bc,'.Solver#Property','UserDefinedData_t')

   for data in data2extract: #we sort each carto on a specific array
      
       data_out = 'out_' + data
       data_out = []

       data_node2extract  = CI.getNodeFromName(data_carto, data)      #we get the node
       data_array2extract = CI.getValue(data_node2extract)            #we extract the array associated to the specific node

       sp                 = data_array2extract.shape
       size               = sp[0]*sp[1]                               #dimension of the array

       data_array2extract = numpy.reshape(data_array2extract, size, order='F') #we reshape the array with a Fortran like order
       data_array2extract = list(data_array2extract)
       p = 0
     #~ # Important : Ordre des boucles = ordre FastS
       if (i1 == im):
		   for k in range(km-1) :
			   for j in range(jm-1) :
				   data_extracted = data_array2extract[p]
				   data_out.append(data_extracted)
				   p += 1

	
       if (j1 == jm):
		   for k in range(km-1) :     
			   for i in range(im-1) :
				   data_extracted = data_array2extract[p]
				   data_out.append(data_extracted)
				   p += 1		     


	
       if (k1 == km):
		   for j in range(jm-1) :     
			   for i in range(im-1) :
				   data_extracted = data_array2extract[p]
				   data_out.append(data_extracted)
				   p += 1
	     
       data_out = numpy.array(data_out)
       CI.createChild(Prop,data,'DataArray_t',value=data_out) #we create a child node for each para
       

   CI.createChild(Prop,'k_piv','DataArray_t',value=kpiv)
           
def addInletData2FastS(data,Name):
    inletDict = OrderedDict()
    inletDict['txv']                 = math.cos(alphaAm*pi/180.)
    inletDict['tyv']                 = math.sin(alphaAm*pi/180.)
    inletDict['tzv']                 = 0.0
    inletDict['stagnation_pressure'] = Pio
    inletDict['stagnation_enthalpy'] = Hio
    inletDict['inj_tur1']            = k1_ini
    
    if turbmod not in ['SA']: inletDict['inj_tur2'] = k2_ini

    for bc in CI.getNodesFromName(data,Name):
       #~ CI.setValue(bc,'BCInflow')
       CI.setValue(bc,'BCInj1')
       Prop=CI.createChild(bc,'.Solver#Property','UserDefinedData_t')
       for d in inletDict: CI.createChild(Prop,d,'DataArray_t',value=inletDict[d])

def addUniformInletData2FastS2(data,Name,d0x,d0y,d0z,pa,ha,turb1):    
   inletDict = OrderedDict()
   
   for bc in CI.getNodesFromName(data,Name):
       #~ CI.setValue(bc,'BCInflow')
       CI.setValue(bc,'BCInj1')
       pr = CI.getNodesFromName(bc,'PointRange')
       index = CI.getValue(pr[0])
       i1 = index[0][0]
       im = index[0][1]
       j1 = index[1][0]
       jm = index[1][1]
       k1 = index[2][0]
       km = index[2][1]
       
       inj_d0x = []
       inj_d0y = []
       inj_d0z = []
       inj_pa  = []
       inj_ha  = []
       inj_turb1 = []
       
       # Important : Ordre des boucles = ordre FastS
       if (i1 == im):
	   for k in range(km-1) :     
	     for j in range(jm-1) :
	       inj_d0x.append(d0x)
	       inj_d0y.append(d0y)
	       inj_d0z.append(d0z)
	       inj_pa.append(pa)
	       inj_ha.append(ha)
	       inj_turb1.append(turb1)
	       
       if (j1 == jm):
	   for k in range(km-1) :     
	     for i in range(im-1) :
	       inj_d0x.append(d0x)
	       inj_d0y.append(d0y)
	       inj_d0z.append(d0z)
	       inj_pa.append(pa)
	       inj_ha.append(ha)
	       inj_turb1.append(turb1)
	       
       if (k1 == km):
	   for j in range(jm-1) :     
	     for i in range(im-1) :
	       inj_d0x.append(d0x)
	       inj_d0y.append(d0y)
	       inj_d0z.append(d0z)
	       inj_pa.append(pa)
	       inj_ha.append(ha)
	       inj_turb1.append(turb1)
	       
       inletDict['txv']                 = numpy.array(inj_d0x)
       inletDict['tyv']                 = numpy.array(inj_d0y)
       inletDict['tzv']                 = numpy.array(inj_d0z)
       inletDict['stagnation_pressure'] = numpy.array(inj_pa)
       inletDict['stagnation_enthalpy'] = numpy.array(inj_ha)
       inletDict['inj_tur1']            = numpy.array(inj_turb1)
	    
       for bc in CI.getNodesFromName(data,Name):
         Prop=CI.createChild(bc,'.Solver#Property','UserDefinedData_t')
         for d in inletDict: CI.createChild(Prop,d,'DataArray_t',value=inletDict[d])

def addUniformInletData2FastS(data,Name,d0x,d0y,d0z,pa,ha,turb1):    
   inletDict = OrderedDict()
   
   for bc in CI.getNodesFromName(data,Name):
       #~ CI.setValue(bc,'BCInflow')
       CI.setValue(bc,'BCInj1')
       pr = CI.getNodesFromName(bc,'PointRange')
       index = CI.getValue(pr[0])
       i1 = index[0][0]
       im = index[0][1]
       j1 = index[1][0]
       jm = index[1][1]
       k1 = index[2][0]
       km = index[2][1]
       
       inj_d0x = []
       inj_d0y = []
       inj_d0z = []
       inj_pa  = []
       inj_ha  = []
       inj_turb1 = []
       
       # Important : Ordre des boucles = ordre FastS
       if (i1 == im):
	   for k in range(km-1) :     
	     for j in range(jm-1) :
	       inj_d0x.append(d0x)
	       inj_d0y.append(d0y)
	       inj_d0z.append(d0z)
	       inj_pa.append(pa)
	       inj_ha.append(ha)
	       inj_turb1.append(turb1)
	       
       if (j1 == jm):
	   for k in range(km-1) :     
	     for i in range(im-1) :
	       inj_d0x.append(d0x)
	       inj_d0y.append(d0y)
	       inj_d0z.append(d0z)
	       inj_pa.append(pa)
	       inj_ha.append(ha)
	       inj_turb1.append(turb1)
	       
       if (k1 == km):
	   for j in range(jm-1) :     
	     for i in range(im-1) :
	       inj_d0x.append(d0x)
	       inj_d0y.append(d0y)
	       inj_d0z.append(d0z)
	       inj_pa.append(pa)
	       inj_ha.append(ha)
	       inj_turb1.append(turb1)
	       
       inletDict['txv']                 = numpy.array(inj_d0x)
       inletDict['tyv']                 = numpy.array(inj_d0y)
       inletDict['tzv']                 = numpy.array(inj_d0z)
       inletDict['stagnation_pressure'] = numpy.array(inj_pa)
       inletDict['stagnation_enthalpy'] = numpy.array(inj_ha)
       inletDict['inj_tur1']            = numpy.array(inj_turb1)
	    
       for bc in CI.getNodesFromName(data,Name):
         Prop=CI.createChild(bc,'.Solver#Property','UserDefinedData_t')
         for d in inletDict: CI.createChild(Prop,d,'DataArray_t',value=inletDict[d])
	 
def addNonUniformInletData2FastS(data, Name, DataFile):    
   data_carto = CP.convertFile2PyTree(DataFile,'fmt_v3d')
   data2extract = ['d0x', 'd0y', 'd0z', 'pa', 'ha', 'tur1']
   
   for bc in CI.getNodesFromName(data,Name):

       CI.setValue(bc,'BCInj1')
       pr = CI.getNodesFromName(bc,'PointRange')
       index = CI.getValue(pr[0])
       
       i1 = index[0][0]
       im = index[0][1]
       j1 = index[1][0]
       jm = index[1][1]
       k1 = index[2][0]
       km = index[2][1]
   
       Prop=CI.createChild(bc,'.Solver#Property','UserDefinedData_t')

   for data in data2extract: #we sort each carto on a specific array
      
       data_inj = 'inj_' + data
       data_inj = []

       data_node2extract  = CI.getNodeFromName(data_carto, data)      #we get the node
       data_array2extract = CI.getValue(data_node2extract)            #we extract the array associated to the specific node

       sp                 = data_array2extract.shape
       size               = sp[0]*sp[1]                               #dimension of the array

       data_array2extract = numpy.reshape(data_array2extract, size, order='F') #we reshape the array with a Fortran like order
       data_array2extract = list(data_array2extract)
       p = 0
     #~ # Important : Ordre des boucles = ordre FastS
       if (i1 == im):
		   for k in range(km-1) :
			   for j in range(jm-1) :
				   data_extracted = data_array2extract[p]
				   data_inj.append(data_extracted)
				   p += 1

	
       if (j1 == jm):
		   for k in range(km-1) :     
			   for i in range(im-1) :
				   data_extracted = data_array2extract[p]
				   data_inj.append(data_extracted)
				   p += 1		     


	
       if (k1 == km):
		   for j in range(jm-1) :     
			   for i in range(im-1) :
				   data_extracted = data_array2extract[p]
				   data_inj.append(data_extracted)
				   p += 1
	     
       data_inj = numpy.array(data_inj)
       CI.createChild(Prop,data,'DataArray_t',value=data_inj) #we create a child node for each para
       	 
def addWallData2FastS(data,Name):
   for bc in CI.getNodesFromName(data,Name):
       if Governing_Equation == 'Euler': 
	   CI.setValue(bc,'BCWallInviscid')
       else:
	   CI.setValue(bc,'BCWallViscous')

def addMobileCoef2FastS(data, famNameHub, famNameCasing):
    dict_wall_hub    = {'mobile_coef':1.0}
    dict_wall_casing = {'mobile_coef':0.0}
    
    for bc in CI.getNodesFromName(data, 'FamilyName'):
		if CI.getValue(bc) == famNameHub:
			(parent, c) = CI.getParentOfNode(data, bc) 
			CI.setValue(parent,'BCWallViscous')
			Prop=CI.createChild(parent,'.Solver#Property','UserDefinedData_t')
			for d in dict_wall_hub: CI.createChild(Prop,d,'DataArray_t',value=dict_wall_hub[d])
			
		if CI.getValue(bc) == famNameCasing:
			(parent, c) = CI.getParentOfNode(data, bc) 
			CI.setValue(parent,'BCWallViscous')
			Prop=CI.createChild(parent,'.Solver#Property','UserDefinedData_t')
			for d in dict_wall_casing: CI.createChild(Prop,d,'DataArray_t',value=dict_wall_casing[d])
			
def addPeriodicData2FastS(data):
    bases=CI.getByType(data,'CGNSBase_t')[2]
    for base in bases:
	zones=CI.getByType(base,'Zone_t')[2]
	j = 1
	for zone in zones:	
	    zonebc = CI.getNodeFromName(zone, 'ZoneBC')
	    connectivity1to1_t = CI.getNodesFromType(zone, 'GridConnectivity1to1_t')
	    i = 0
	    # Recherche des noeuds de raccord periodique
	    for node in connectivity1to1_t :
	      children1to1 = CI.getChildren(node)
	      if children1to1 != [] : 
	       connect_name = CI.getNodesFromName(children1to1, 'GridConnectivityProperty')
	       if connect_name != [] :
		 periodic = CI.getNodesFromName(children1to1, 'Periodic')
		
		 if periodic != [] :
		    i += 1
		    # Recuperation des donnees des noeuds de raccord periodique
		    pointRange       = CI.getNodesFromName(children1to1, 'PointRange')
		    pointRangeN      = CI.getName(pointRange[0])
		    pointRangeV      = CI.getValue(pointRange[0])
		    
		    pointRangeDonor  = CI.getNodesFromName(children1to1, 'PointRangeDonor')
		    pointRangeDonorN = CI.getName(pointRangeDonor[0])
		    pointRangeDonorV = CI.getValue(pointRangeDonor[0])
		    
		    transform        = CI.getNodesFromName(children1to1, 'Transform')
		    transformN       = CI.getName(transform[0])
		    transformV       = CI.getValue(transform[0])
		    
		    rotation_center  = CI.getNodesFromName(periodic, 'RotationCenter')
		    rotation_centerN = CI.getName(rotation_center[0])
		    rotation_centerV = CI.getValue(rotation_center[0])
		    
		    rotation_angle   = CI.getNodesFromName(periodic, 'RotationAngle')
		    rotation_angleN  = CI.getName(rotation_angle[0])
		    rotation_angleV  = CI.getValue(rotation_angle[0])
		    
		    translation      = CI.getNodesFromName(periodic, 'Translation')
		    translationN     = CI.getName(translation[0])
		    translationV     = CI.getValue(translation[0])
		    
		    # Creation de la condition aux limites BCPeriodic. Si le noeud zoneBC n existe pas, il faut le creer.
		    if zonebc == None : 
			CI.newZoneBC(parent=zone)
			zonebc = CI.getNodeFromName(zone, 'ZoneBC')
			
		    bcper = CI.newBC(name='bndper%i%i' % (i,j), pointRange=pointRangeV, btype='UserDefined',parent=zonebc)
		    CI.setValue(bcper,'BCPeriodic')
		    
		    dict_per = OrderedDict()
		    dict_per[rotation_angleN]  = rotation_angleV
		    dict_per[translationN]     = translationV
		    dict_per[rotation_centerN] = rotation_centerV
		    
		    # Creation du noeud .Solver#Property
		    Prop=CI.createChild(bcper,'.Solver#Property','UserDefinedData_t')
		    
		    CI.createChild(Prop,pointRangeDonorN,'IndexRange_t',value=pointRangeDonorV)
		    CI.createChild(Prop,transformN,'"int[IndexDimension]"',value=transformV)
		    for d in dict_per: CI.createChild(Prop,d,'DataArray_t',value=dict_per[d])
		    
		    #~ # Correction du noeud RotationAngle pour rotation sens trigonometrique et non direct!
		    #~ PropAngle = CI.getNodesFromName(Prop, 'RotationAngle')
		    #~ PropAngle_Corrected_Value = -1*CI.getValue(PropAngle[0])
		    #~ CI.setValue(PropAngle[0],PropAngle_Corrected_Value)
		    
	    j += 1

def cleanPeriodicData2FastS(data):
    bases=CI.getByType(data,'CGNSBase_t')[2]
    for base in bases:
	zones=CI.getByType(base,'Zone_t')[2]
	for zone in zones:	
	    zonebc = CI.getNodeFromName(zone, 'ZoneBC')
	    connectivity1to1_t = CI.getNodesFromType(zone, 'GridConnectivity1to1_t')
	    # Recherche des noeuds de raccord periodique
	    for node in connectivity1to1_t :
	      data = CI._rmNodesByName(node,'.Solver#Property')
	      children1to1 = CI.getChildren(node)
	      if children1to1 != [] : 
	       connect_name = CI.getNodesFromName(children1to1, 'GridConnectivityProperty')
	       if connect_name != [] :
		 periodic = CI.getNodesFromName(children1to1, 'Periodic')
		
		 if periodic != [] :
		    # Modification du noeud rotationAngle pour rotation d axe Ox
		    rotation_angle   = CI.getNodesFromName(periodic, 'RotationAngle')
		    rotation_angleV  = CI.getValue(rotation_angle[0])
		    
		    # rotation_angleV[0] = Nb d aubes
		    # rotation_angleV[1] = Nb de canaux
		    # Pour FastS :
		    rotation_angleV[0] = -1*360.*rotation_angleV[1]/rotation_angleV[0]
		    rotation_angleV[1] = 0.
		    rotation_angleV[2] = 0.
		    
		    CI.setValue(rotation_angle[0],numpy.array(rotation_angleV))

def prepareDebit2Fast(t_debit):
    zones      = CI.getNodesFromType2(t_debit, 'Zone_t') 
    tmp        = CI.getNodesFromName1(zones,'FlowSolution#Centers')
    node_debit = CI.getNodesFromName1(tmp,'Density')
    
    return node_debit
	
def extractDebit2Fast(node_debit):
    debit =0.
    for dens in node_debit:
	    density  = CI.getValue(dens)
	    debit   +=  numpy.sum(density)
	    
    return debit
		    	   
def cleanUpByType2FastS(data,Type):
    data = CI._rmNodesByType(data,Type)

def cleanUpByName2FastS(data,Name):
    data = CI._rmNodesByName(data,Name)
    
def cleanUpByValue2FastS(data,Value):
    data = CI._rmNodesByValue(data,Value)
