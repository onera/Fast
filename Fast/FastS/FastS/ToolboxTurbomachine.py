import Converter.PyTree as CP
import Converter.Internal as CI
import numpy
import math
from collections import OrderedDict
from itertools import takewhile, dropwhile
#import matplotlib.pyplot as plt
import Converter   as CV

#from scipy import interpolate
#~ ====================
#~ from etc.post  import AzimutalAverage2D, AzimutalAverage3D
#~ from   mpi4py import MPI
#~ import etc.toolbox.internal  as tgi
#~ import Converter.Array3D as CA
#~ ====================

from userData2 import *

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
        ind -= 1
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

def getDataFiles(nbBlock, rootFileName):
    dataFiles = []

    for i in range(nbBlock):
        fileName = rootFileName + str(i+1)
        dataFiles.append(fileName)

    return dataFiles

def getArrayAndReshapeF(tree, data2extract):
    node  = CI.getNodeFromName(tree, data2extract)               #we get the node
    array = CI.getValue(node)                                   #we get the array associated
    #~ size = array.size
    #~ array = numpy.reshape(array, size, order='F')               #we reshape the array with a Fortran like order
    return array

def cleanFlowSolution2FastS(data):
    var = ['Vt_abs','rhoVt_abs', 'theta', 'Radius' ]
    for v in var: CI.rmNodesByName(data, v)

def createFlowSolution(data):
    bases = CI.getByType(data,'Zone_t')[2]
    for base in bases:
        Prop=CI.createChild(base,'FlowSolution#Centers','FlowSolution_t')

# XXX COMPLETE
def addFlowSolution2FastS(data, nbBlock, rootDataFileName):
    #This function creates a flowSolution node with concervative variable and set the data (at the center of the cells) with restart files

    #DATA
    # t                 : input          [CGNS tree] FastS
    # nbBlock           : input userData [int] total number of blocks
    # rootDataFileName  : input userData [str] root name of the restart data files

    bases=CI.getByType(data,'Zone_t')[2]
    for base in bases:
        Prop=CI.createChild(base,'FlowSolution#Centers','FlowSolution_t')

    data2setList   = ['Density', 'MomentumX', 'MomentumY', 'MomentumZ', 'EnergyStagnationDensity', 'TurbulentEnergyKineticDensity']

    dataFilesList  = getDataFiles(nbBlock,rootDataFileName)       #we get a list of all the files used to set the FlowSolution node (restart files)
    dict_flowSol = OrderedDict()

    p = 0
    bases=CI.getByType(data,'Zone_t')[2]                          #list of the blocks
    for base in bases:                                            #for each block
        dataFile = dataFilesList[p]                               #we get the right file
        tdata2set = CP.convertFile2PyTree(dataFile,'bin_tp')      #that we transform into a tree

        Prop=CI.getNodeFromName(base, 'FlowSolution#Centers')     #we get the flowSol node
        for d in data2setList:                                                                    #we set the value
            value2set = getArrayAndReshapeF(tdata2set, d)
            CI.createChild(Prop,d,'DataArray_t',value=value2set,pos=-1)

        p += 1

def addUniformFlowSolution2FastS(t, Model = 'NSTurbulent'):

    ret = CP.isNamePresent(t, 'centers:Density')
    if ret != 1: # Density not present
        state = CI.getNodeFromType(t, 'ReferenceState_t')
        if state is None:
            raise ValueError('Reference state is missing in input cgns.')
        vars = ['Density', 'MomentumX', 'MomentumY', 'MomentumZ',
                'EnergyStagnationDensity']
        for v in vars:
            node = CI.getNodeFromName(state, v)
            if node is not None:
                val = float(node[1][0])
                CP._initVars(t, 'centers:'+v, val)
            else:
                raise ValueError(v + ' is missing in ReferenceState.')
        if Model == 'NSTurbulent':
            vars = ['TurbulentSANuTildeDensity']
            for v in vars:
                node = CI.getNodeFromName(state, v)
                if node is not None:
                    val = float(node[1][0])
                    CP._initVars(t, 'centers:'+v, val)

def getArray(tree, data2extract):
    node  = CI.getNodeFromName(tree, data2extract )             #we get the node
    array = CI.getValue(node)                                   #we get the array associated

    return array

def getMeshBlocks(tree, nbBlock):
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

def getIndexOfBlock(tree, block):
    node  = CI.getNodeFromName(tree, block)
    index = CI.getValue(node)

    im = index[0][0]
    jm = index[1][0]
    km = index[2][0]

    return im, jm, km

def closestSearch(myList, nbDef):
    g = lambda x: x <= nbDef

    try:
        nbLTnbDef = list(takewhile( g, myList)) #we store all the numbers that are lower than the user defined number
        nbGTnbDef = list(dropwhile( g, myList)) #we store all the numbers that are larger than the user defined number

        nbMax     = max(nbLTnbDef) #the closest from the left
        nbMin     = min(nbGTnbDef) #the closest from the right

        difMax = abs(nbDef - nbMax)
        difMin = abs(nbDef - nbMin)


        if difMax < difMin : return nbMax, list(myList).index(nbMax)   #return the value and its position
        else : return nbMin, list(myList).index(nbMin)

    except:
        if nbDef not in myList:
            print('==============ERROR ==============')
            print('The value selected is out of range')
            print('==================================')

def pivSearch(radiusList, rad_def, zone='inside'):
    #There are 3 options for zone ; 'inside', 'shaft', 'casing'
    #depending on the data given by the user

    try :

        #The data are chosen inside the radius distribution
        if zone == 'inside':
            rpiv, indpiv = closestSearch(radiusList, rad_def)
            return rpiv, indpiv

        #The data are choosen at the shaft
        if zone == 'shaft':
            indpiv = 0
            rpiv   = radiusList[0]
            return rpiv, indpiv


        #The data are choosen at the casing
        if zone == 'casing':
            indpiv = len(radiusList)-1
            rpiv   = radiusList[indpiv]
            return rpiv, indpiv

    except :
        print('===============================ERROR ===============================')
        print('Option not understood: select "casing", "shaft", or default:"inside"')
        print('====================================================================')

def addFlowEquationSet2FastS(data):
    bases=CI.getByType(data,'CGNSBase_t')[2]

    for base in bases:
        Prop=CI.createChild(base,'FlowEquationSet','FlowEquationSet_t')
        CI.createChild(Prop,'EquationDimension','"int"',value=Equation_Dimension)
        CI.createChild(Prop,'GoverningEquations','GoverningEquations_t',value=Governing_Equation)
        CI.createChild(Prop,'TurbulenceModel','TurbulenceModel_t',value=Turbulence_Model)

def cleanGhostCells(data, nbBlock, nbGhostRank):

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

    def getNode(tree, block, nodeName):
        base        = CI.getNodesFromType(tree,'CGNSBase_t')
        baseN       = CI.getName(base)
        nodePath    =  '/' + CI.getName(baseN) + '/' + block + '/' + nodeName + '/'

        flowSolNode = CI.getNodeFromPath(tree, nodePath)

        return flowSolNode

    def getArray(tree, data2extract):
        node  = CI.getNodeFromName(tree, data2extract )             #we get the node
        array = CI.getValue(node)                                   #we get the array associated

        return array

    def getChildrenDict(tree, nodeName):

        childrenNodeList = CI.getChildren(nodeName) #we get all the children nodes of the nodes
        nbChildren       = len(childrenNodeList)
        childrenDict     = OrderedDict()            #we create a dictionary with key: the name of the node, value: the node

        for p in range(nbChildren):
            childNode                = childrenNodeList[p]
            childnName               = CI.getName(childNode)
            childrenDict[childnName] = childNode

        return childrenDict

    def getIndexOfBlock(tree, block):
        node  = CI.getNodeFromName(tree, block)
        index = CI.getValue(node)

        im = index[0][1]
        jm = index[1][1]
        km = index[2][1]

        return [im, jm, km]


    Data2setCoordinatesList = ['CoordinateX','CoordinateY', 'CoordinateZ']
    Data2setListFlowSol_0   = ['Density', 'TurbulentDistance', 'TurbulentEnergyKineticDensity', 'ViscosityEddy', 'VelocityX', 'VelocityY', 'VelocityZ', 'Temperature', 'TurbulentSANuTilde']
    Data2setListFlowSol_M1  = ['Density_M1', 'VelocityX_M1', 'VelocityY_M1', 'VelocityZ_M1', 'Temperature_M1', 'TurbulentSANuTilde_M1']
    Data2setListFlowSol_P1  = ['Density_P1', 'VelocityX_P1', 'VelocityY_P1', 'VelocityZ_P1', 'Temperature_P1', 'TurbulentSANuTilde_P1']
    Data2setListFlowSolList = Data2setListFlowSol_0 + Data2setListFlowSol_M1 + Data2setListFlowSol_P1

    nb2setCoordinates = len(Data2setCoordinatesList) #number of arrays where the ghost cells will be removed
    nb2setFlowSol     = len(Data2setListFlowSolList)


#~ ================Scripts================

    blockMeshList = getMeshBlocks(data)

    for p in range(nbBlock):                                     #for each block (block by block)
        blockMesh = blockMeshList[p]

        coordinateNode = getNode(data , blockMesh, 'GridCoordinates')                #we get the corresponding coordinate node =====> Nom a changer
        coordinateNode2 = getNode(data , blockMesh, 'GridCoordinates#Init')                #we get the corresponding coordinate node =====> Nom a changer
        flowSolNode    = getNode(data , blockMesh, 'FlowSolution#Centers')           #we get the corresponding FlowSolution node =====> Nom a changer

        indexList = getIndexOfBlock(data, blockMesh)                                #we get the index of the mesh including the ghost cells

        chilDictCoordinates = getChildrenDict(data, coordinateNode)  #we create a dict with = key : the children nodes' name/ value = the corresponding node
        chilDictCoordinates2 = getChildrenDict(data, coordinateNode2)  #we create a dict with = key : the children nodes' name/ value = the corresponding node
        chilDictFlowSol     = getChildrenDict(data, flowSolNode)

        i = 0
        j = 0

#~ For the Flow solution nodes
        for q in range(nb2setFlowSol):
            FlowSolNode  = Data2setListFlowSolList[q]
            FlowSolArray = getArray(data, FlowSolNode)     #we extract the corresponding array for each children node

            for axis in range(3):
                indexMax = indexList[axis]    #the max index in the "axis" direction (for each axis)

                FlowSolArray = numpy.delete(FlowSolArray, 0,  axis = axis) #we delete the 2 first ranks
                FlowSolArray = numpy.delete(FlowSolArray, 0,  axis = axis)

                FlowSolArray = numpy.delete(FlowSolArray, indexMax-3, axis = axis) # and the last 2 ranks (we deleted 2 ranks so, there are (indexMax-1)-2 ranks left)
                FlowSolArray = numpy.delete(FlowSolArray, indexMax-4, axis = axis)

                CI.setValue(chilDictFlowSol[FlowSolNode], FlowSolArray) #we replace with the new array without the ghost cells




#~ For the coordinate nodes
        for q in range(nb2setCoordinates):
            CoordinateNode  = Data2setCoordinatesList[q]
            CoordinateArray = getArray(data, CoordinateNode)     #we extract the corresponding array for each children node

            for axis in range(3):
                indexMax = indexList[axis]    #the max index in the "axis" direction (for each axis)

                CoordinateArray = numpy.delete(CoordinateArray, 0,  axis = axis) #we delete the 2 first ranks
                CoordinateArray = numpy.delete(CoordinateArray, 0,  axis = axis)

                CoordinateArray = numpy.delete(CoordinateArray, indexMax-2, axis = axis) # and the last 2 ranks (we deleted 2 ranks so, there are (indexMax)-2 ranks left)
                CoordinateArray = numpy.delete(CoordinateArray, indexMax-3, axis = axis) #it's indexMax and not indexMax-1 because the coordinate nodes are note centered (+1)

                CI.setValue(chilDictCoordinates[CoordinateNode], CoordinateArray)
                CI.setValue(chilDictCoordinates2[CoordinateNode], CoordinateArray)



        return data

def getFlowSolNode(tree, block):
    base        = CI.getNodesFromType(tree,'CGNSBase_t')
    baseN       = CI.getName(base)
    nodePath    =  '/' + CI.getName(baseN) + '/' + block + '/FlowSolution#Centers/'

    flowSolNode = CI.getNodeFromPath(tree, nodePath)

    return flowSolNode

# XXX
def convertRelative2Absolute(data, nbBlock, omega, gamma ):
    # This function convert the data from



    # Definition of the functions needed to convert our data
    def norm3D(x,y,z)                                 : return math.sqrt(x**2 + y**2 + z**2)
    def atan2(x,y)                                    : return math.atan2(x,y)

    def Ps(gamma, Density, E, Urel, Vrel, Wrel)       : return ((gamma-1)*Density)*(E-0.5*norm3D(Urel, Vrel, Wrel)**2)
    def roE_abs(Ps, gamma, Density, Uabs, Vabs, Wabs) : return (Ps/(gamma-1))+0.5*Density*norm3D(Uabs, Vabs, Wabs)**2

    def Vabs(V_rel, omega, CoordinateZ)               : return V_rel - omega*CoordinateZ
    def Wabs(W_rel, omega, CoordinateY)               : return W_rel + omega*CoordinateY

#~ ================Scripts================
    blockMeshList = getMeshBlocks(data, nbBlock)                                            #we get a list of all the blocks in the input file

    for p in range(nbBlock):                                     #for each block (block by block)
        blockMesh = blockMeshList[p]
        print('Processing conversion relative to absolute ===', blockMesh)

        zoneInit = CI.getNodeFromName(data, blockMesh)
        zone = CI.getNodeFromName(data, blockMesh)     #we extract the mesh zone
        zone = CP.node2Center(zone)                    #we change the coordinates (ONLY) to have them at the center of the cell

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


        #absolute velocity
        zone = CP.initVars(zone, '{Uabs}={Urel}')
        zone = CP.initVars(zone, 'Vabs', Vabs, ['Vrel', 'omg', 'CoordinateZ'])
        zone = CP.initVars(zone, 'Wabs', Wabs, ['Wrel', 'omg', 'CoordinateY'])

        #reconstruction of absolute conservative variables
        zone = CP.initVars(zone, '{MomentumX}={Uabs}*{Density}')
        zone = CP.initVars(zone, '{MomentumY}={Vabs}*{Density}')
        zone = CP.initVars(zone, '{MomentumZ}={Wabs}*{Density}')
        zone = CP.initVars(zone, 'EnergyStagnationDensity', roE_abs, ['Ps', 'gamma', 'Density', 'Uabs', 'Vabs', 'Wabs'])

        #we get the FlowSolution node to set the data
        zone = CI.renameNode(zone, 'FlowSolution', 'FlowSolution#Centers')
        flowSol = getFlowSolNode(data, blockMesh)


        ConsVariables    = ['MomentumX', 'MomentumY', 'MomentumZ', 'EnergyStagnationDensity']

        for i in range(len(ConsVariables)):
            ConsVariable = ConsVariables[i]

            #we delete the previous nodes
            relativeNode   = CI.getNodeFromName(flowSol, ConsVariable)
            CI.rmNode(flowSol, relativeNode)
            #we add what we want on the node
            absoluteNode  = CI.getNodeFromName(zone, ConsVariable)
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

# XXX
def addOutletData2FastS(data,Name):
    dict_press = OrderedDict()
    dict_press['pressure'] = Ps_aval
    dict_press['k_piv']    = kpiv

    for bc in CI.getNodesFromName(data,Name):
        CI.setValue(bc,'BCOutpres')
        Prop=CI.createChild(bc,'.Solver#Property','UserDefinedData_t')
        for d in dict_press: CI.createChild(Prop,d,'DataArray_t',value=dict_press[d])

# XXX
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

# XXX
def addNonUniformOutletData2FastS(data,Name,DataFile, kpiv):

    data_carto = CP.convertFile2PyTree(DataFile)
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

# XXX
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

# XXX
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

def addNonUniformInletData2FastS(t, Name, DataFile):
    data_carto = CP.convertFile2PyTree(DataFile)
    data2extract = ['d0x', 'd0y', 'd0z', 'pa', 'ha', 'tur1']

    zones = CI.getZones(t)
    for z in zones:
        bcs =  CI.getNodesFromType(z,'BC_t')

        for bc in bcs:
            if Name in bc[0]:
                CI.setValue(bc,'BCInj1')
                pr = CI.getNodesFromName(bc,'PointRange')

                Prop=CI.createUniqueChild(bc,'.Solver#Property','UserDefinedData_t')

                zone_carto = CI.getNodeFromName(data_carto, z[0])
                for data in data2extract: #we sort each carto on a specific array

                    data_node2extract  = CI.getNodeFromName(zone_carto, data)      #we get the node
                    data_array2extract = CI.getValue(data_node2extract)            #we extract the array associated to the specific node

                    sp       = data_array2extract.shape
                    size     = sp[0]*sp[1]                               #dimension of the array

                    data_array2extract = numpy.reshape(data_array2extract, size, order='F') #we reshape the array with a Fortran like order

                    CI.createUniqueChild(Prop,data,'DataArray_t',value=data_array2extract) #we create a child node for each para

    return None

###
###
###
def addNonUniformBCData(t, BCName, BCValue, DataFile):

    data_carto   = CP.convertFile2PyTree(DataFile)

    zones = CI.getZones(t)
    for z in zones:
        bcs =  CI.getNodesFromType(z,'BC_t')

        for bc in bcs:
            if BCName in bc[0]:
                CI.setValue(bc, BCValue)

                Prop       = CI.createUniqueChild(bc,'.Solver#Property','UserDefinedData_t')
                zone_carto = CI.getNodeFromName(data_carto, z[0])
                sol_carto  = CI.getNodeFromName(zone_carto, 'FlowSolution')
                vars       = CI.getNodesFromType(sol_carto, 'DataArray_t')
                for v in vars:
                    varBC = numpy.reshape( v[1], v[1].size, order='F')                 #we reshape the array with a Fortran like order

                    CI.createUniqueChild( Prop, v[0],'DataArray_t',value= varBC)          #we create a child node for each para

    return None


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

#~ =====================================================================
#~ RADIAL EQUILIBRIUM
#~ =====================================================================
#~ Initialization of radial equilibrium

# XXX
def printStep(step):

    print('#################################################')
    print(step)
    print('#################################################')

# XXX COMPLETE
def rhoVtheta(t, outBlock):
    #This function creates the nodes and calculates the value of 'Vt_abs', 'rhoVt_abs', 'theta', 'Radius'
    # at the cell center (FlowSolution#Centers) in the outflow block, and set the value at 0 for the other blocks

    #DATA
    # t             : input          [CGNS tree] FastS
    # outBlock      : input userData [int] name of the upstream block

    def atan2(x,y)                                     : return math.atan2(x,y)
    #def Vt_abs(Vabs, Wabs, theta)                      : return Vabs*math.cos(theta) - Wabs*math.sin(theta)

    zones = CI.getZones(t)                                                  #we get a list of all the blocks in the input file
    for z in zones :

        if z[0] == outBlock:                                            #for each block (block by block)
            print('init rhoVtheta Processing=== for zone:', outBlock)

            zc = CP.node2Center(z)                                                  #we change the coordinates (ONLY) to have them at the center of the cell

            #we add what is needed to convert our data
            zc = CP.initVars(zc, 'theta', atan2, ['CoordinateY', 'CoordinateZ'])
            zc = CP.initVars(zc,'{Radius}=numpy.sqrt({CoordinateY}**2+{CoordinateZ}**2)')
            zc = CP.initVars(zc,'Vt_abs'   , 0.)

            #zone = CP.initVars(zone, 'omg', omega)
            #zone = CP.initVars(zone, '{Vabs}={MomentumY}/{Density}')
            #zone = CP.initVars(zone, '{Wabs}={MomentumZ}/{Density}')
            #zone = CP.initVars(zone, 'Vt_abs', Vt_abs, ['Vabs', 'Wabs', 'theta'])
            #absolute velocity theta component
            #z = CP.initVars(z, '{rhoVt_abs}', Vt_abs, ['MomentumY', 'MomentumZ', 'theta'])
            #z = CP.initVars(z, '{Vt_abs}={Vt_abs}/{Density}')

            #we get the FlowSolution node to set the data
            zc = CI.renameNode(zc, 'FlowSolution', 'FlowSolution#Centers')
            flowSol = CI.getNodeFromName(z,'FlowSolution#Centers')
            if flowSol == None:
                node =['FlowSolution#Centers', None, [], 'FlowSolution_t' ]
                CI.addChild(z, node)
                flowSol = CI.getNodeFromName(z,'FlowSolution#Centers')


            #var = ['Vt_abs', 'rhoVt_abs', 'theta', 'Radius']
            var = ['Vt_abs', 'theta', 'Radius']
            for v in var:
                node = CI.getNodeFromName(zc, v)
                CI.addChild(flowSol, node)

# XXX ALMOST, c
'''
def azimAverage2D(t, nbband, c):
#This function calculates the 2D density, rhoVt(absolute), Vt(absolute) and Radius azimutale average at the exit of the outflow block

#DATA
# t             : input          [CGNS tree] FastS
# nbband        : input userData [int] number of slices of azimutal average 
# c             : input userData [int] ?????
# radius_star   : input          [1D numpy array] radial repartition of the azimutale averaged radius distribution

# tree          : return         [CGNS tree] new tree with the azimutale average of "Density", "rhoVt_abs", "Vt_abs", "Radius"
    
#~ ====================
        from etc.post  import AzimutalAverage2D, AzimutalAverage3D
        from   mpi4py import MPI 
        import etc.toolbox.internal  as tgi 
        import Converter.Array3D as CA
#~ ====================
        
        variables     = ["Density", "rhoVt_abs", "Vt_abs", "Radius"]            
        t = CP.initVars(t,'{Radius}=numpy.sqrt({CoordinateY}**2+{CoordinateZ}**2)')
        t = CP.node2Center(t,'Radius')
        
#~ -------------------------------------------------------------------------------
        # Azimean 2D extraction from 3D to 2D only the outblow block
        CellCenter   = CI.__FlowSolutionCenters__ = 'FlowSolution#Centers'
        bc           = tgi.renameBC(CP.extractBCOfName(t, 'FamilySpecified:Aval'))  
        bc           = tgi.keepNodesByNameAndType(bc, [CellCenter], 'FlowSolution_t')
#~ -------------------------------------------------------------------------------      
        
        azm          = AzimutalAverage2D(bc)
        wbsum, sbsum = azm.compute(variables, fsname=CellCenter, nbband=nbband, c=c)
        
        
        # Creation of the tree 
        tree    = CI.newCGNSTree()
        base    = CI.newCGNSBase('Base', 3, 3, parent=tree)
        zone    = CI.newZone('Zone',[[wbsum.shape[0]],[0],[0]],  'Structured', parent=base)
        flowSol = CI.newFlowSolution(name='FlowSolution',gridLocation='CellCenter', parent=zone)
        
        # Extraction and setting of the data
        for v in variables:
           node = CI.createNode(v,'DataArray_t',value=wbsum[:,variables.index(v)])
           CI.addChild(flowSol,node)
           
    return tree 
# XXX COMPLETE
def interpolationAzim(t, tinit, interpolation): 
#This function interpolates the azimutale averaged derivative of the pressure from an azimutale averaged radius repartition 
#to the real derivative of the pressure with respect to the real radius distribution

#DATA
# t             : input          [CGNS tree] tree with the azimutale average data (dpdr*)
# tinit         : input          [CGNS tree] FastS initial tree
# interpolation : input userData [str] option for interpolation between computed azimutal radius and real radius of the outBlock mesh

# bc_1d         : return         [CGNS tree] tree with 1D distribution of the pressure derivative with respect to the radius and radius distribution

#~ ====================
        from etc.post  import AzimutalAverage2D, AzimutalAverage3D
        from   mpi4py import MPI 
        import etc.toolbox.internal  as tgi 
        import Converter.Array3D as CA
#~ ====================
        
        #~ ------------------these lines are used to extract the radius at the cell's center
        CellCenter  = CI.__FlowSolutionCenters__ = 'FlowSolution#Centers'
        bc          = tgi.renameBC(CP.extractBCOfName(tinit, 'FamilySpecified:Aval')) 
        bc          = tgi.keepNodesByNameAndType(bc, [CellCenter], 'FlowSolution_t')


        # Radius 'r' at the cell center
        bc_c  = CP.node2Center(bc) 
        CI.rmNode(bc_c,CI.getNodeFromName(bc_c,'FlowSolution'))
        bc_c  = CP.initVars(bc_c,'{radius_c}=numpy.sqrt({CoordinateY}**2+{CoordinateZ}**2)')
        bc_c  = CP.initVars(bc_c,'{dpdr_c}=0.0') #2D 
        bc_1d = T.subzone(bc_c,(1,1,1),(1,-1,1)) # attention on suppose que axe 'rayon' == axe maillage 'j', a retravailler pour trouver cet axe 'rayon'
        ray_c = CI.getValue(CI.getNodeFromName(bc_1d,'radius_c'))
        
        # Radius 'r*' after azimutal average
        ray_star = CI.getValue(CI.getNodeFromName(t,'Radius'))
        dpdr_star = CI.getValue(CI.getNodeFromName(t,'dpdr_star')) 

        #we can't interpolate if the data range of r is not include in the r* one => 
        #if rmin <r*min we replace rmin by r*min
        n=0
        for i in range(ray_c.size): 
                if ray_c[i] < ray_star[0]:
                        print 'we change at n = ', n, 'old radius', ray_c[i], 'new radius', ray_star[i]
                        ray_c[i] = ray_star[i]
                        n+=1
        
        # Interpolation 'r*' -> 'r'
        print 'interpolate not available'
        #f_dpdr = interpolate.interp1d(ray_star, dpdr_star, kind=interpolation) #interpolation function generator
        #dpdr   = f_dpdr(ray_c)
        dpdr   = 0
        dpdr_c = CI.getNodeFromName(bc_1d,'dpdr_c')  
        CI.setValue(dpdr_c,dpdr)
        
        return bc_1d
'''

# XXX COMPLETE
def RadEq(t, ppiv_def, rpiv_def, zone = 'inside' ):
    #This function calculates the 1D radial pressure distribution according to radial equilibrium (recurrence)

    #DATA
    # t             : input          [CGNS tree] tree with 1D distribution of the pressure derivative with respect to the radius and radius distribution
    # ppiv_def      : input userData [int] pivot pressure
    # rpiv_def      : input userData [int] pivot radius
    # zone          : input userData [str] casing/shaft/inside for the starting point of the reccurence algorithm for the static pressure

    # tree          : return          [CGNS tree] tree with 1D distribution of the pressure derivative with respect to the radius and radius distribution

    radius   = getArray(t, 'radius_c')
    dpsdr    = getArray(t, 'dpdr_c')
    nb_rad   = len(list(radius))


    #============================CASING===========================
    if zone == 'casing':
        p = [ppiv_def]
        rpiv, indPiv = pivSearch(radius, rpiv_def, zone='casing')

        for i in list(reversed(range(nb_rad-1))): #(-1) avoids border effect and goes from the casing to the shaft (reversed)
            p_before = p[0] - 0.5*(dpsdr[i] + dpsdr[i+1])*(radius[i+1]-radius[i])
            p.insert(0, p_before) #the last pressure is put at the begining of the list



    #=============================SHAFT===========================
    if zone == 'shaft':
        rpiv, indPiv = pivSearch(radius, rpiv_def, zone='shaft')
        p = [ppiv_def]

        for i in range(nb_rad-1): #avoid border effect
            p_next = p[i] + 0.5*(dpsdr[i] + dpsdr[i+1])*(radius[i+1]-radius[i])
            p.append(p_next)

    #============================INSIDE===========================
    if zone == 'inside':
        rpiv, indPiv = pivSearch(radius, rpiv_def, zone='inside')
        p = [ppiv_def]

        for i in range(nb_rad-indPiv-1):  #from piv to the casing

            ind = indPiv + i
            p_next = p[i] + 0.5*(dpsdr[ind] + dpsdr[ind+1])*(radius[ind+1]-radius[ind])
            p.append(p_next)


        for i in list(reversed(range(indPiv))): #from piv to the shaft
            p_before = p[0] - 0.5*(dpsdr[i] + dpsdr[i+1])*(radius[i+1]-radius[i])
            p.insert(0, p_before) #the last pressure is put at the begining of the list
    #=============================================================

    pArray = numpy.array(p)
    #PLOT-------------------------------------------------
    #plt.figure(facecolor="white")
    #plt.title("Pressure evolution along the span",  fontsize=20)
    #plt.xlabel('Radius (m)', fontsize=18)
    #plt.ylabel('Pressure (Pa)', fontsize=18)

    #plt.plot(radius, pArray ,marker="*", label=zone)
    #plt.legend(loc = 'best', prop={'size':20})
    #plt.show()
    #~ #------------------------------------------------------

    # Data saving
    tree    = CI.newCGNSTree()
    base    = CI.newCGNSBase('Base', 3, 3, parent=tree)
    zone    = CI.newZone('Zone', [[0],[0],[88]],  'Structured', parent=base)
    flowSol = CI.newFlowSolution(name='FlowSolution',gridLocation='CellCenter', parent=zone)

    node1 = CI.createNode('radius','DataArray_t',radius)
    CI.addChild(flowSol,node1)
    node2 = CI.createNode('p','DataArray_t',  pArray )
    CI.addChild(flowSol,node2)

    return tree

# XXX COMPLETE
def outradeqExtension(t, tdata, outBlock):
    #This function transforms p(r) in p(r,theta) , from 1D to 2D

    #DATA
    # t             : input          [CGNS tree] initial tree with the outflow bloc
    # tdata         : input          [CGNS tree] tree with 1D distribution of the pressure derivative with respect to the radius and radius distribution
    # outBlock      : input userData [int] name of the upstream block

    # data          : return         [2D numpy arrays] ['r', 'theta', 'p']

    #~ ====================
    import Converter.Array3D as CA
#~ ====================

    im, jm, km   = getIndexOfBlock(t, outBlock)

    # ####################################################################
    #1D data distribution to extend in 2D
    # ####################################################################

    ra = CI.getValue(CI.getNodeFromName(tdata, 'radius'))                                           #array of radius 1D
    pa = CI.getValue(CI.getNodeFromName(tdata, 'p'))                                                        #array of p 1D

    nr = km-1                                                                                                                                       # nb of points in  radial direction (-1 because cell centered)
    nt = jm-1
    # nb of points in  azimutal direction (-1 because cell centered)
    # ####################################################################
    #1D angular distribution
    # ####################################################################

    #~ def atan2(x,y)  : return math.atan2(x,y)*180/math.pi
    zone = CI.getNodeFromName(t, outBlock)
    zone = CP.node2Center(zone)                                                                             #from vertex to cell centered
    zone = CP.initVars(zone, 'theta', lambda x,y: math.atan2(x,y)*180/math.pi , ['CoordinateY', 'CoordinateZ'])
    thetaNode = CI.getNodeFromName(zone, 'theta')                               #3D
    theta = CI.getValue(thetaNode)
    theta2D = theta[im-2, :, :].transpose()                                                                         #extraction 2D slice : at imax = outflow

    # ####################################################################
    # Transforms p(r) in p(r,theta)
    # ####################################################################


    rt  = numpy.multiply.outer(numpy.multiply.outer(ra,numpy.ones(nt)),numpy.ones(1))
    tt  = numpy.multiply.outer(theta2D,numpy.ones(1))
    prt = numpy.multiply.outer(numpy.multiply.outer(pa,numpy.ones(nt)),numpy.ones(1))                       # p(r,theta)

    rt = numpy.transpose(rt)
    tt = numpy.transpose(tt)
    prt = numpy.transpose(prt)

    var = []
    var.append(rt)
    var.append(tt)
    var.append(prt)
    varname = ['r', 'theta', 'p']


    data = CA.convertArrays3D2Arrays([[varname,var]])
    return data


# XXX ALMOST, c
def addOutradeq2FastS(t, nbBlock,  outBlock, omega,  nbband, c, ppiv_def, rpiv_def = 0, location ='inside', interpolation = 'linear' ):
    #This function calculates the 2D pressure field of the outflow block (.Solver#Property) according to the radial equilibrium
    #it creates the nodes needed to update the solution while running computation. It saves the pressure fields on 'outradeq.dat'

    #DATA
    #~ There are 3 options : casing (requires only the pressure at the casing, the radius is calculated)
    #~                                           shaft(requires only the pressure at the shaft, the radius is calculated)
    #~                                           inside(requires the pressure ans the radius)  in any case, the value of radius is set at 0 if not given

    # t             : input          [CGNS tree] FastS
    # nbBlock       : input userData [int] total number of blocks
    # outBlock      : input userData [int] name of the upstream block
    # omega         : input userData [int] rotational speed in rad/s
    # nbband        : input userData [int] number of slices of azimutal average
    # c             : input userData [int] ?????
    # ppiv_def      : input userData [int] pivot pressure
    # rpiv_def      : input userData [int] pivot radius
    # location      : input userData [str] casing/shaft/inside for the starting point of the reccurence algorithm for the static pressure
    # interpolation : input userData [str] option for interpolation between computed azimutal radius and real radius of the outBlock mesh

    #~ ====================
    #from etc.post  import AzimutalAverage2D, AzimutalAverage3D
    #~ ====================

    tinit = t

# ====================================================================================================
#  STEP1 = rhoVtheta calculus
# ====================================================================================================
    rhoVtheta(t, nbBlock, outBlock, omega)                                                          #3D calculates rhoVtheta and creates the node
    printStep('STEP1 = rhoVtheta calculus DONE' )
# =====================================================================================================
#  STEP2 = azimutal average calculus = *
# =====================================================================================================
    #t2average = azimAverage2D(t, nbband, c)                                                                #2D azim average

    def dpdr(rhoVtheta, density, radius):   return (rhoVtheta**2)/(radius*density)
    t2average = CP.initVars(t2average,'dpdr_star', dpdr, ['rhoVt_abs', 'Density', 'Radius'])

    printStep('STEP2 = azimutal average calculus = * DONE' )
# ======================================================================================================
#  STEP3 =  interpolation from azim average *
# =======================================================================================================
    #t3 = interpolationAzim(t2average, tinit, interpolation)
    printStep('STEP3 =  interpolation from azim average * DONE' )
# ========================================================================================================
#  STEP4 =  Radial equilibrium
# ========================================================================================================
    t4 = RadEq(t3, ppiv_def, rpiv_def, zone = location)                             #1D
    printStep('STEP4 =  Radial equilibrium DONE' )
# ========================================================================================================
#  STEP5 =  azimutal extension/carto creation
# ========================================================================================================
    data = outradeqExtension(tinit, t4, outBlock)                                           #2D
    CV.convertArrays2File(data,'outradeq.plt','bin_tp')
    CV.convertArrays2File(data,'outradeq.dat')
    printStep('STEP5 =  azimutal extension DONE' )

# ====================================================================================================
#~ RADIAL EQUILIBRIUM
# ====================================================================================================
#~ Update of radial equilibrium

# XXX COMPLETE
def realRadius1D(t, outBlock):
    #This function calculates the real radius distribution along the span of the outflow block exit

    #DATA
    # t         : input          [CGNS tree] FastS
    # outBlock  : input userData [int] name of the upstream block

    # radius_1D : return         [1D numpy array] real radial distribution along the span of the outflow block exit

    # Radius at the cell center
    blk_aval = CI.getNodeFromName(t,outBlock)
    im,jm,km = getIndexOfBlock(t, outBlock)

    radius      = CI.getValue(CI.getNodeFromName(blk_aval,"Radius"))
    #radius_1D   = radius[im-2,jm-2,:] # attention on suppose que axe 'rayon' == axe maillage 'j', a retravailler pour trouver cet axe 'rayon'
    radius_1D   = radius[im-3,jm-3,:] # attention on suppose que axe 'rayon' == axe maillage 'j', a retravailler pour trouver cet axe 'rayon'

    return radius_1D

# XXX ALMOST c
'''
def avgRadius1D(t, outBlock, nbband, c):
#This function calculates the azimutale averaged radius distribution along the span of the outflow block exit

#DATA
# t         : input          [CGNS tree] FastS  
# outBlock  : input userData [int] name of the upstream block
# nbband    : input userData [int] number of slices of azimutal average 
# c         : input userData [int] ?????

# radius_1D : return         [1D numpy array] azimutale averaged radial distribution along the span of the outflow block exit
#~ ====================
        from etc.post  import AzimutalAverage2D, AzimutalAverage3D
        from   mpi4py import MPI 
        import etc.toolbox.internal  as tgi 
        import Converter.Array3D as CA
#~ ====================
        
        variables     = ["Radius"]              
        
#~ -------------------------------------------------------------------------------
        # Azimean 2D extraction from 3D to 2D only the outblow block
        CellCenter   = CI.__FlowSolutionCenters__ = 'FlowSolution#Centers'
        bc           = tgi.renameBC(CP.extractBCOfName(t, 'FamilySpecified:Aval'))  
        bc           = tgi.keepNodesByNameAndType(bc, [CellCenter], 'FlowSolution_t')
        print 'bc',CI.printTree(bc)
#~ -------------------------------------------------------------------------------      
        
        azm          = AzimutalAverage2D(bc)
        print 'azm',CI.printTree(azm)
        
        nbband_loc = nbband
        c_loc      = c
        wbsum, sbsum = azm.compute(variables, fsname=CellCenter, nbband=nbband_loc, c=c_loc)    
        print 'azm.radius',azm.radius
        print 'wbsum',wbsum

        radius_1D =  wbsum[:,variables.index('Radius')]
        
        return radius_1D
'''
# XXX COMPLETE
def adaptRadius(radius_real, radius_star):
    #we can't interpolate if the data range of r is not include in the r* one =>
    #if rmin <r*min we replace rmin by r*0, #if rmax > r*max we replace rmin by r*max

    #DATA
    # radius_real : input  [1D numpy array] real radial distribution along the span of the outflow block exit
    # radius_star : input  [1D numpy array] azimutale averaged radial distribution along the span of the outflow block exit

    # radius_real : return  [1D numpy array] real radial distribution included in  azimutale averaged radial distribution

    nbRad = len(radius_real)

    for i in range(nbRad):
        if radius_real[i] < radius_star[0]:
            radius_real[i] = radius_star[0]
            print('LOWER PART ==> we change at i = ', i, 'old radius=', radius_real[i], 'new radius=', radius_star[i])

        if radius_real[i] > radius_star[nbRad-1]:
            radius_real[i] = radius_star[nbRad-1]
            print('UPPER PART ==> we change at i = ', i, 'old radius=', radius_real[i], 'new radius=', radius_star[nbRad-1])

    return radius_real


#~ COMPLETE
def _rhoVthetaUpdate(t, outBlock):
    #This function updates the value of Vt_abs and rhoVt_abs (FlowSolution#Centers) in the outflow block

    #DATA
    # t             : input [CGNS tree] FastS
    # nbBlock       : input userData [int] total number of blocks
    # outBlock      : input userData [int] name of the upstream block

    zones = CI.getZones(t)
    for zone in zones :

        if zone[0] == outBlock:                                         #for the  outflow block
            sol = CI.getNodeFromName(zone, 'FlowSolution#Centers')                                  #we extract the mesh zone
            vteta = CI.getNodeFromName(sol, 'Vt_abs')[1]
            vy    = CI.getNodeFromName(sol, 'VelocityY')[1]
            vz    = CI.getNodeFromName(sol, 'VelocityZ')[1]
            theta = CI.getNodeFromName(sol, 'theta')[1]

            #vteta[0:]= vy[0:]*numpy.cos(theta[0:]) - vz[0:]*numpy.sin(theta[0:])
            vteta[0:]= vz[0:]*numpy.cos(theta[0:]) - vy[0:]*numpy.sin(theta[0:])

            #for k in range(vteta.shape[2]):
            #       for j in range(vteta.shape[1]):
            #               for i in range(vteta.shape[0]):
            #print 'shape',vteta.shape
            #for k in range(vteta.shape[2]):
            #  for j in range(0,44):
            #print 'vteta=', vteta[82, j, k], j,k

# XXX ALMOST c
def azimAverage2DUpdate_(t, varlist, dir_avg, outBlock, bctype):

    zones = CI.getZones(t)
    for zone in zones :

        if zone[0] == outBlock:                                          #for the  outflow block

            o         = CI.getNodeFromName1( zone, '.Solver#ownData')
            bcs       = CI.getNodesFromType2(zone, 'BC_t')
            param_int = CI.getNodeFromName1( o, 'Parameter_int')[1]
            param_real= CI.getNodeFromName1( o, 'Parameter_real')[1]

            varout=[]

            for bc in bcs:
                name    = CI.getValue(bc)
                if name == bctype:
                    sol     = CI.getNodeFromName(zone, 'FlowSolution#Centers')           #we extract the mesh zone
                    Ptrange = CI.getNodeFromType1( bc , 'IndexRange_t')
                    indrange= CI.getValue(Ptrange)
                    ind_bc  = numpy.zeros(6, numpy.int32)
                    ind_bc[0] = indrange[0][0]-1
                    ind_bc[1] = indrange[0][1]-1
                    ind_bc[2] = indrange[1][0]-1
                    ind_bc[3] = indrange[1][1]-1
                    ind_bc[4] = indrange[2][0]-1
                    ind_bc[5] = indrange[2][1]-1
                    if  ind_bc[1]==  ind_bc[0]: #imin/imax
                        if    ind_bc[1] == 0: ind_bc[1]=           2
                        else:                 ind_bc[1]= ind_bc[1]-3
                        ind_bc[0]= ind_bc[1]
                        idir = 1
                    elif ind_bc[3] == ind_bc[2]: #jmin/jmax
                        if    ind_bc[3] == 0: ind_bc[3]=           2
                        else:                 ind_bc[3]= ind_bc[3]-3
                        ind_bc[2]= ind_bc[3]
                        idir = 2
                    else:
                        if    ind_bc[5] == 0: ind_bc[5]=           2
                        else:                 ind_bc[5]= ind_bc[3]-3
                        ind_bc[4]= ind_bc[5]
                        idir = 3

                    if   dir_avg==1:
                        size = max(1,(ind_bc[3]-ind_bc[2]))*max(1,(ind_bc[5]-ind_bc[4]))

                    elif dir_avg==2:
                        size = max(1,(ind_bc[1]-ind_bc[0]))*max(1,(ind_bc[5]-ind_bc[4]))
                    else           :
                        size = max(1,(ind_bc[1]-ind_bc[0]))*max(1,(ind_bc[3]-ind_bc[2]))

                    if   idir ==1: ind_bc[1] = ind_bc[1] + 1
                    elif idir ==2: ind_bc[3] = ind_bc[3] + 1
                    else         : ind_bc[5] = ind_bc[5] + 1

                    for v in varlist:
                        #print 'name',name, size, v
                        var_avg   = numpy.zeros(size, numpy.float64)
                        var = CI.getNodeFromName(sol, v)[1]


                        if   dir_avg==1:
                            c = 0
                            for k in range(ind_bc[4],ind_bc[5]):
                                for j in range(ind_bc[2],ind_bc[3]):
                                    nb= 0
                                    for i in range(ind_bc[0],ind_bc[1]):
                                        var_avg[c]= var_avg[c] +  var[i, j, k]
                                        nb+=1
                                    var_avg[c]=var_avg[c]/float(nb)
                                    #print 'var_avg',var_avg[c]
                                    c+=1
                        elif dir_avg==2:

                            c = 0
                            for k in range(ind_bc[4],ind_bc[5]):
                                for i in range(ind_bc[0],ind_bc[1]):
                                    nb= 0
                                    for j in range(ind_bc[2],ind_bc[3]):
                                        var_avg[c]= var_avg[c] +  var[i, j, k]
                                        nb+=1
                                    var_avg[c]=var_avg[c]/float(nb)
                                    #print 'var_avg',var_avg[c],c,k,i,nb
                                    c+=1
                        else:

                            c = 0
                            for j in range(ind_bc[2],ind_bc[3]+1):
                                for i in range(ind_bc[0],ind_bc[1]+1):
                                    nb= 0
                                    for k in range(ind_bc[4],ind_bc[5]+1):
                                        var_avg[c]= var_avg[c] +  var[i, j, k]
                                        nb+=1
                                    var_avg[c]=var_avg[c]/float(nb)
                                    #print 'var_avg',var_avg[c]
                                    c+=1
                        varout.append([outBlock, v, var_avg, nb])

    return varout

# XXX ALMOST c
def azimAverage2DUpdate(t, nbband, c, radius_star):
    #This function calculates the 2D density and rhoVt absolute azimutale average at the exit of the outflow block

    #DATA
    # t             : input          [CGNS tree] FastS
    # nbband        : input userData [int] number of slices of azimutal average
    # c             : input userData [int] ?????
    # radius_star   : input          [1D numpy array] radial repartition of the azimutale averaged radius distribution

    # tree          : return         [CGNS tree] new tree with the azimutale average of rho and rhoVt absolute

    #~ ====================
    from etc.post  import AzimutalAverage2D, AzimutalAverage3D
    from   mpi4py import MPI
    import etc.toolbox.internal  as tgi
    import Converter.Array3D as CA
#~ ====================

    variables    = ["Density", "rhoVt_abs"]

#~ -------------------------------------------------------------------------------
    # Azimean 2D extraction from 3D to 2D only the outblow block
    CellCenter   = CI.__FlowSolutionCenters__ = 'FlowSolution#Centers'
    bc           = tgi.renameBC(CP.extractBCOfName(t, 'FamilySpecified:Aval'))
    bc           = tgi.keepNodesByNameAndType(bc, [CellCenter], 'FlowSolution_t')
#~ -------------------------------------------------------------------------------

    azm          = AzimutalAverage2D(bc)
    wbsum, sbsum = azm.compute(variables, fsname=CellCenter, nbband=nbband, c=c)


    # Creation of the tree
    tree    = CI.newCGNSTree()
    base    = CI.newCGNSBase('Base', 3, 3, parent=tree)
    zone    = CI.newZone('Zone',[[wbsum.shape[0]],[0],[0]],  'Structured', parent=base)
    flowSol = CI.newFlowSolution(name='FlowSolution',gridLocation='CellCenter', parent=zone)

    # Extraction and setting of the data

    for v in variables:
        node = CI.createNode(v,'DataArray_t',value=wbsum[:,variables.index(v)])
        CI.addChild(flowSol,node)

    radius_starNode = CI.createNode('radius_star','DataArray_t',value=radius_star)
    CI.addChild(flowSol,radius_starNode)

    return tree

# XXX COMPLETE
def dpdr_starCalculus(tAverage):
    # This function calculates the 1D azimutale averaged derivative of the pressure with respect to the averaged radius

    #DATA
    # tAverag    : input  [CGNS tree] tree with the azimutale average of rho and rhoVt absolute
    # dpdr_star  : return [1D numpy array] azimutale averaged derivative of the pressure with respect to the averaged radius

    def dpdr(rhoVtheta, density, radius):   return (rhoVtheta**2)/(radius*density)

    tAverage  = CP._initVars(tAverage,'dpdr_star', dpdr, ['rhoVt_abs', 'Density', 'radius_star'])
    dpdr_star = CI.getValue(CI.getNodeFromName(tAverage, 'dpdr_star'))

    return dpdr_star

# XXX COMPLETE
def interpolationAzimUpdate(radius_star, dpdr_star, radius_real, interpolation='linear'):
    #This function interpolates the azimutale averaged derivative of the pressure from an azimutale averaged radius repartition
    #to the real derivative of the pressure with respect to the real radius distribution

    #DATA
    # radius_star   : input          [1D numpy array] azimutale averaged radial distribution along the span of the outflow block exit
    # dpdr_star     : input          [1D numpy array] azimutale averaged derivative of the pressure with respect to the averaged radius
    # radius_real   : input          [1D numpy array] real radial distribution along the span of the outflow block exit
    # interpolation : input userData [str] option for interpolation between computed azimutal radius and real radius of the outBlock mesh

    # dpdr          : return         [1D numpy array] real derivative of the pressure with respect to the real radius distribution

    print('interpolate not available')
    #f_dpdr = interpolate.interp1d(radius_star, dpdr_star, kind=interpolation) #interpolation function generator
    #dpdr   = f_dpdr(radius_real)
    dpdr   = 0

    return dpdr

# XXX COMPLETE
def RadEqUpdateUpdate(dpdr1D, radius_real1D, ppiv_def, rpiv_def, zone = 'inside' ):
    #This function calculates the radial pressure distribution according to radial equilibrium (recurrence)

    #DATA
    # dpdr1D        : input          [1D numpy array] real derivative of the pressure with respect to the real radius distribution
    # radius_real   : input          [1D numpy array] real radial distribution
    # ppiv_def      : input userData [int] pivot pressure
    # rpiv_def      : input userData [int] pivot radius
    # zone          : input userData [str] casing/shaft/inside for the starting point of the reccurence algorithm for the static pressure

    # pArray        : return         [1D numpy array] presure distribution along the span

    radius   = radius_real1D
    dpsdr    = dpdr1D
    nb_rad   = len(list(radius))


    #============================CASING===========================
    if zone == 'casing':
        p = [ppiv_def]
        rpiv, indPiv = pivSearch(radius, rpiv_def, zone='casing')

        for i in list(reversed(range(nb_rad-1))): #(-1) avoids border effect and goes from the casing to the shaft (reversed)
            p_before = p[0] - 0.5*(dpsdr[i] + dpsdr[i+1])*(radius[i+1]-radius[i])
            p.insert(0, p_before) #the last pressure is put at the begining of the list



    #=============================SHAFT===========================
    if zone == 'shaft':
        rpiv, indPiv = pivSearch(radius, rpiv_def, zone='shaft')
        p = [ppiv_def]

        for i in range(nb_rad-1): #avoid border effect
            p_next = p[i] + 0.5*(dpsdr[i] + dpsdr[i+1])*(radius[i+1]-radius[i])
            p.append(p_next)

    #============================INSIDE===========================
    if zone == 'inside':
        rpiv, indPiv = pivSearch(radius, rpiv_def, zone='inside')
        p = [ppiv_def]

        for i in range(nb_rad-indPiv-1):  #from piv to the casing

            ind = indPiv + i
            p_next = p[i] + 0.5*(dpsdr[ind] + dpsdr[ind+1])*(radius[ind+1]-radius[ind])
            p.append(p_next)


        for i in list(reversed(range(indPiv))): #from piv to the shaft
            p_before = p[0] - 0.5*(dpsdr[i] + dpsdr[i+1])*(radius[i+1]-radius[i])
            p.insert(0, p_before) #the last pressure is put at the begining of the list
    #=============================================================

    pArray = numpy.array(p)

    return pArray

def outradeqExtensionUpdate_(t, outBlock, press1D, dir, nb_avg, iteration, display = 100000000):

    size    = press1D.size*nb_avg
    press2d = numpy.empty(size, numpy.float64)

    #a generaliser (permutation jk pour matcher ordre eciture ficher dataArray
    if dir == 1:
        c=0
        for j in range(press1D.size):
            for i in range(nb_avg):
                press2d[c] =  press1D[j]
                c+=1

    elif dir == 2:

        c=0
        for j in range(press1D.size):
            for i in range(nb_avg):
                press2d[c] =  press1D[j]
                c+=1

    return press2d

# XXX COMPLETE
def outradeqExtensionUpdate(tinit, outBlock,  radius1D, press1D, iteration, display = 100000000):
    #This function transforms p(r) in p(r,theta) , from 1D to 2D and save the results every "display"

    #DATA
    # tinit         : input          [CGNS tree] initial tree with the outflow block
    # outBlock      : input userData [int] name of the upstream block
    # radius1D      : input          [1D numpy array] real radial distribution
    # press1D       : input          [1D numpy array] 1D presure distribution along the span
    # iteration     : input          [int] number of iteration
    # display       : input userData [int] data save parameter, save the 2D pressure field every "display"

    # prt           : return         [2D numpy array] 2D presure distribution along the span


    import Converter.Array3D as CA
    im, jm, km   = getIndexOfBlock(tinit, outBlock)

    # ####################################################################
    #1D data distribution to extend in 2D
    # ####################################################################

    ra = radius1D
    pa = press1D

    nr = km-5                                                                                               # nb of points in  radial direction (-1 because cell centered)
    nt = jm-5                                                                                               # nb of points in  azimutal direction (-1 because cell centered)


    # ####################################################################
    #1D angular distribution
    # ####################################################################


    zone      = CI.getNodeFromName(tinit, outBlock)
    thetaNode = CI.getNodeFromName(zone, 'theta')            #3D

    theta = CI.getValue(thetaNode)
    theta2D = theta[im-2, :, :].transpose()                          #extraction 2D slice : at imax = outflow

    # ####################################################################
    # Transforms p(r) in p(r,theta)
    # ####################################################################

    rt  = numpy.multiply.outer(numpy.multiply.outer(ra,numpy.ones(nt)),numpy.ones(1))
    tt  = numpy.multiply.outer(theta2D,numpy.ones(1))
    prt = numpy.multiply.outer(numpy.multiply.outer(pa,numpy.ones(nt)),numpy.ones(1))       # p(r,theta)

    rt  = numpy.transpose(rt)
    tt  = numpy.transpose(tt)
    prt = numpy.transpose(prt)

    var = []
    var.append(rt)
    var.append(tt)
    var.append(prt)
    varname = ['r', 'theta', 'p']

    prt2 = prt[0, :, :]
    press2NoGhost =  prt[0 ,2:-2 ,2:-2]
    press2NoGhost = numpy.reshape(press2NoGhost, press2NoGhost.size, order='F')

    text = CA.convertArrays3D2Arrays([[varname,var]])

    if iteration % display == 0:
        fileName = 'outradeq' + str(iteration)
        CV.convertArrays2File(text, fileName + '.dat')

    return press2NoGhost

#~ XXX ALMOST c
def _updateOutradeq2FastS(iteration, t, rad_real1D, rad_star1D,  nbBlock,  outBlock, omega,   nbband, c, ppiv_def, rpiv_def, location ='inside', interpolation = 'linear', display=100000000 ):
    #This function updates the 2D pressure field of the outflow block (.Solver#Property) according to the radial equilibrium

    #DATA
    #~ There are 3 options : casing (requires only the pressure at the casing, the radius is calculated)
    #~                                           shaft(requires only the pressure at the shaft, the radius is calculated)
    #~                                           inside(requires the pressure ans the radius)  in any case, the value of radius is set at 0 if not given

    # iteration     : input userData [int] ????
    # t             : input          [CGNS tree] FastS
    # rad_real1D    : input              [1D numpy array] real radial distribution along the span of the outflow block exit
    # radius_star   : input          [1D numpy array] azimutale averaged radial distribution along the span of the outflow block exit
    # nbBlock       : input userData [int] total number of blocks
    # outBlock      : input userData [int] name of the upstream block
    # omega         : input userData [int] rotational speed in rad/s
    # nbband        : input userData [int] number of slices of azimutal average
    # c             : input userData [int] ?????
    # ppiv_def      : input userData [int] pivot pressure
    # rpiv_def      : input userData [int] pivot radius
    # location      : input userData [str] casing/shaft/inside for the starting point of the reccurence algorithm for the static pressure
    # interpolation : input userData [str] option for interpolation between computed azimutal radius and real radius of the outBlock mesh
    # display       : input userData [int] data save parameter, save the 2D pressure field every "display"

    # t             : return [CGNS tree] the starting tree with the updated pressure field at the exit

    #tinit = t


    # ====================================================================================================
    #  STEP1 = rhoVtheta calculus
    # ====================================================================================================
    #printStep('STEP1 = rhoVtheta calculus Processing' )
    _rhoVthetaUpdate(t, outBlock)                                                                                      #3D
    #printStep('STEP1 = rhoVtheta calculus DONE' )
# =====================================================================================================
#  STEP2 = azimutal average calculus = *
# =====================================================================================================
    #tAverage = azimAverage2DUpdate(t, nbband, c, rad_star1D)       #2D azim average
    dir = 2    # on moyenne suivant l'indice j
    bctype = 'BCOutpres'
    varlist=['Vt_abs','Density','Radius']

    var_avg = azimAverage2DUpdate_(t, varlist, dir, outBlock, bctype )

    #print 'Variable moyennee',var_avg
    #CI.printTree(tAverage)

    #printStep('STEP2 = azimutal average calculus = * DONE' )
# ======================================================================================================
#  STEP3 =  interpolation from azim average * to real radius distribution
# ======================================================================================================
    #dpdr_star1D = dpdr_starCalculus(tAverage)                                                        #1D
    #dpdr_real1D = interpolationAzimUpdate(rad_star1D, dpdr_star1D, rad_real1D, interpolation=interpolation)

    vt_abs = var_avg[0][2]
    density= var_avg[1][2]
    radius = var_avg[2][2]

    dpdr_real1D = numpy.empty(vt_abs.size, numpy.float64)
    dpdr_real1D[0:]= vt_abs[0:]*vt_abs[0:]*density[0:]/radius[0:]

    #printStep('STEP3 =  interpolation from azim average * DONE' )
# ======================================================================================================
#  STEP4 = radial equilibrium
# ======================================================================================================
    #press1D = RadEqUpdateUpdate(dpdr_real1D, rad_real1D, ppiv_def, rpiv_def, zone = location )      #1D
    press1D = RadEqUpdateUpdate(dpdr_real1D, radius, ppiv_def, rpiv_def, zone = location )      #1D
    #print('press_carter=',press1D[87])
    #printStep('STEP4 =  Radial equilibrium DONE' )
# ======================================================================================================
#  STEP5 =  azimutal extension/carto creation
# ======================================================================================================
    #press2D = outradeqExtensionUpdate(tinit, outBlock,  rad_real1D, press1D, iteration)                                      #2D
    press2D = outradeqExtensionUpdate_(t, outBlock, press1D, dir, var_avg[0][3], iteration)                                   #2D
    #print press2D.size
    #printStep('STEP5 =  azimutal extension DONE' )
# ======================================================================================================
#  STEP6 =  2D pressure setting
# ======================================================================================================
    zones = CI.getZones(t)
    for zone in zones:
        if zone[0] == outBlock:                                      #for the  outflow block
            bcs       = CI.getNodesFromType2(zone, 'BC_t')
            for bc in bcs:
                name    = CI.getValue(bc)
                if name == bctype:
                    CI.getNodeFromName(t, 'p')[1][:]=press2D[:]
                    #pressure = CI.getNodeFromName(bc, 'p')[1]
                    #pressure[:]= press2D[:]
                    #for i in range (pressure.size):
                    #   pressure[i]= press2D[i]

    #CI._setValue(CI.getNodeFromName(t, 'p'), numpy.copy(press2D))

    #print press2D
    #printStep('STEP6 =  2D pressure setting DONE' )

    return None
