"""FastC for IBM preprocessing"""
import Converter.PyTree as C
import Converter.Distributed as Distributed
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import Converter
import Connector.Mpi as Xmpi
import Connector.OversetData as XOD
import Connector.PyTree as X
import Connector.IBM as C_IBM
import Connector
import Generator.PyTree as G
import Generator.IBMmodelHeight as G_IBM_Height
import Initiator.PyTree as I
import Geom.IBM as D_IBM
import Transform.PyTree as T
import Compressor.PyTree as Compressor
import Dist2Walls.PyTree as DTW
import Distributor2.PyTree as D2
import Post.Mpi as Pmpi
import Post.PyTree as P
import KCore.test as test
import Generator
import Transform
import KCore
import numpy
import math

from . import fastc

varsn       = ['gradxTurbulentDistance','gradyTurbulentDistance','gradzTurbulentDistance']
varsnDouble = ['gradxTurbulentDistanceDouble','gradyTurbulentDistanceDouble','gradzTurbulentDistanceDouble']
TOLDIST     = 1.e-14
SHIFTF      = 1.e-10
EPSCART     = 1.e-6
TOLCELLN    = 0.01

LBM_IBC_NUM = 113 #from param_solver.h

TypesOfIBC  = XOD.TypesOfIBC

def prepareIBMData(t_case, t_out, tc_out, t_in=None, to=None, tbox=None, tinit=None, tbCurvi=None,
                   snears=0.01, snearsf=None, dfars=10., dfarDir=0, vmin=21, depth=2, frontType=1, octreeMode=0,
                   IBCType=1, verbose=True, expand=3, ext=-1, order=2, extrap=1, dTarget=1000,
                   check=False, twoFronts=False, cartesian=True, cleanCellN=True, frontFile=None,
                   yplus=100., Lref=1., correctionMultiCorpsF42=False, blankingF42=False, wallAdaptF42=None, heightMaxF42=-1.):

    import Generator.IBM as G_IBM
    import time as python_time

    optimized=-1  #sinon  faut appeler connector legacy
    nature   = 0
    if ext == -1:
        ext = depth+1
        if optimized==-1: ext=depth

    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = Internal.copyTree(t_case)

    ## Note: cartesian = True is left as an input argument to avoid regressing  during the non-regression test.
    ##       In the near future the ref. values for the non-regression tests will be updated with cartesian=True.
    ##       At this point, cartesian=True input argument can be deleted.
    ## Note: when cartesian = True is deleted as an input argument the line below must be uncommented and the cartesian in the if statement but be deleted.
    #cartesian = True
    if t_in and cartesian:
        cartesian = G_IBM.checkCartesian(t_in, nghost=2)
        if cartesian:
            RED  = "\033[1;31;40m"
            END  = "\033[0m"
            print("===========================================")
            print("Note: t_in is a " + RED + "CARTESIAN " + END + "grid")
            print("===========================================")
        else:
            RED  = "\033[1;31;40m"
            END  = "\033[0m"
            print("===========================================")
            print("Note: t_in is " + RED + "NOT" + END + " a " + RED + "CARTESIAN " + END + "grid")
            print("===========================================")

    refstate = Internal.getNodeFromName(tb, 'ReferenceState')
    flowEqn  = Internal.getNodeFromName(tb, 'FlowEquationSet')

    Reynolds = Internal.getNodeFromName(tb, 'Reynolds')
    if Reynolds is not None:
        Reynolds = Internal.getValue(Reynolds)
        if Reynolds < 1.e5: frontType = 1
    else: Reynolds = 1.e6

    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    if dimPb is None: raise ValueError('prepareIBMData: EquationDimension is missing in input geometry tree.')
    dimPb = Internal.getValue(dimPb)
    if dimPb == 2: C._initVars(tb, 'CoordinateZ', 0.)
    #else: cleanCellN = True

    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('prepareIBMData: GoverningEquations is missing in input geometry tree.')
    model = Internal.getValue(model)

    ibctypes = Internal.getNodesFromName(tb, 'ibctype')
    if ibctypes is None: raise ValueError('prepareIBMData: ibc type is missing in input geometry tree.')
    ibctypes = list(set(Internal.getValue(ibc) for ibc in ibctypes))

    if model == 'Euler':
        if any(ibc in ['Musker', 'MuskerMob', 'Mafzal', 'Log', 'TBLE', 'TBLE_FULL'] for ibc in ibctypes):
            raise ValueError("prepareIBMData: governing equations (Euler) not consistent with ibc types %s"%(ibctypes))

    if frontType == 42 and tbox is None:
        print("Info: prepareIBMData: frontType 42 is used, but no tbox has been provided to ensure that the near-wall resolution is sufficiently propagated. Forcing expand=4.")
        expand = 4

    #===================
    # STEP 0 : GET FILAMENT BODIES
    #===================
    tb, tbFilament              = D_IBM.determineClosedSolidFilament__(tb)
    isFilamentOnly, isWireModel = D_IBM.localWMMFlags__(tb, tbFilament)


    #===================
    # STEP 1 : GENERATE MESH
    #===================
    if isFilamentOnly: tbLocal=tbFilament
    elif tbFilament: tbLocal=Internal.merge([tb,tbFilament])
    else: tbLocal=tb


    if t_in is None:
        if verbose: pt0 = python_time.time(); C_IBM.printTimeAndMemory__('generate Cartesian mesh', time=-1)
        test.printMem("Info: prepareIBMData: generate Cartesian mesh [start]")
        t = G_IBM.generateIBMMesh(tbLocal, vmin=vmin, snears=snears, dimPb=dimPb, dfars=dfars, tbox=tbox,
                                  snearsf=snearsf, check=check, to=to, ext=ext, optimized=optimized,
                                  expand=expand, dfarDir=dfarDir, octreeMode=octreeMode)
        Internal._rmNodesFromName(tb,"SYM")
        test.printMem("Info: prepareIBMData: generate Cartesian mesh [end]")


        #C_IBM._redispatch__(t=t)
        if verbose: C_IBM.printTimeAndMemory__('generate Cartesian mesh', time=python_time.time()-pt0)

    else:
        t = t_in

    for b in Internal.getBases(t):
        Internal.addChild(b, refstate, pos=0)
        Internal.addChild(b, flowEqn , pos=0)

    #===================
    # STEP 2 : DIST2WALL
    #===================
    if verbose: pt0 = python_time.time(); C_IBM.printTimeAndMemory__('compute wall distance', time=-1)
    C_IBM._dist2wallIBM(t, tb, dimPb=dimPb, frontType=frontType, Reynolds=Reynolds, yplus=yplus, Lref=Lref,
                        correctionMultiCorpsF42=correctionMultiCorpsF42, heightMaxF42=heightMaxF42, dTarget=dTarget,
                        tbFilament=tbFilament, cleanCellN=cleanCellN)
    #tbFilament=tbFilament, cleanCellN=cleanCellN, verbose=verbose)
    if verbose: C_IBM.printTimeAndMemory__('compute wall distance', time=python_time.time()-pt0)

    #===================
    # STEP 3 : BLANKING IBM
    #===================
    if verbose: pt0 = python_time.time(); C_IBM.printTimeAndMemory__('blank by IBC bodies', time=-1)
    C_IBM._blankingIBM(t, tb, dimPb=dimPb, frontType=frontType, IBCType=IBCType, depth=depth,
                       Reynolds=Reynolds, yplus=yplus, Lref=Lref, twoFronts=twoFronts,
                       heightMaxF42=heightMaxF42, correctionMultiCorpsF42=correctionMultiCorpsF42,
                       wallAdaptF42=wallAdaptF42, blankingF42=blankingF42,
                       tbFilament=tbFilament, cleanCellN=cleanCellN)


    Cmpi.barrier()
    #C_IBM._redispatch__(t=t)
    if verbose: C_IBM.printTimeAndMemory__('blank by IBC bodies', time=python_time.time()-pt0)
    #===================
    # STEP 4 : INTERP DATA CHIM
    #===================
    #mise a zero cellN dans Ghost BC
    for z in Internal.getZones(t):
        bcs   = Internal.getNodesFromType2(z, 'BC_t')
        sol   = Internal.getNodeFromName1(z, "FlowSolution#Centers")
        cellN = Internal.getNodeFromName1(sol,"cellNChim")[1]
        for bc in bcs:
            btype = Internal.getValue(bc)
            if btype != "BCOverlap":
                ptrange = Internal.getNodesFromType1(bc, 'IndexRange_t')
                rg      = ptrange[0][1]
                if rg[0,1]==rg[0,0]:
                    if rg[0,1]==1: cellN[0:2,:,:]=1
                    else         : cellN[-2:,:,:]=1
                elif rg[1,1]==rg[1,0]:
                    if rg[1,1]==1: cellN[:, 0:2,:]=1
                    else         : cellN[:, -2:,:]=1
                elif rg[2,1]==rg[2,0] and dimPb==3:
                    if rg[2,1]==1: cellN[:,:,0:2]=1
                    else         : cellN[:,:,2: ]=1

    if verbose: pt0 = python_time.time(); C_IBM.printTimeAndMemory__('compute interpolation data (Abutting & Chimera)', time=-1)
    tc = C.node2Center(t)

    if Internal.getNodeFromType(t, "GridConnectivity1to1_t") is not None:
        Xmpi._setInterpData(t, tc, nature=1, loc='centers', storage='inverse', sameName=1, dim=dimPb, itype='abutting', order=2, cartesian=cartesian)


    setInterpDataAndSetInterpTransfer__(t,tc, nature=nature, loc='centers', storage='inverse', sameName=1, sameBase=1, dim=dimPb, order=order, extrap=extrap, cartesian=cartesian, corner=True)

    if verbose: C_IBM.printTimeAndMemory__('compute interpolation data (Abutting & Chimera)', time=python_time.time()-pt0)
    #===================
    # STEP 4 : BUILD FRONT
    #===================
    if verbose: pt0 = python_time.time(); C_IBM.printTimeAndMemory__('build IBM front', time=-1)

    t, tc, front, front2, frontWMM = buildFrontIBM(t, tc, tb=tb, dimPb=dimPb, frontType=frontType,
                                                   cartesian=cartesian, twoFronts=twoFronts, check=check,
                                                   tbFilament=tbFilament, optimized=optimized)

    if frontFile is not None: front = C.convertFile2PyTree(frontFile)
    '''
    for z in Internal.getZones(t):
         fastc._updateNatureForIBMGhost(z,
                                        Internal.__GridCoordinates__,
                                        Internal.__FlowSolutionNodes__,
                                        Internal.__FlowSolutionCenters__)
    '''
    if verbose: C_IBM.printTimeAndMemory__('build IBM front', time=python_time.time()-pt0)

    #on recalcule les interp chimere en supprimant les coin maintenant que le calcul du front est ok
    if verbose: pt0 = python_time.time(); C_IBM.printTimeAndMemory__('compute interpolation data 2nd pass corner (Abutting & Chimera)', time=-1)

    setInterpDataAndSetInterpTransfer__(t, tc, nature=nature, loc='centers', storage='inverse', sameName=1, sameBase=1, dim=dimPb,\
                                        extrap=extrap, order=order, cartesian=cartesian, corner=False)

    if verbose: C_IBM.printTimeAndMemory__('compute interpolation data 2nd pass corner (Abutting & Chimera)', time=python_time.time()-pt0)

    #C.convertPyTree2File(t, 'verif_cellNChim.cgns')
    #===================
    # STEP 5 : INTERP DATA IBM
    #===================
    if verbose: pt0 = python_time.time(); C_IBM.printTimeAndMemory__('compute interpolation data (IBM)', time=-1)
    nature_loc = 1
    extrap_loc = 1
    val        = 1 #ghost donneuse autorisee

    _setInterpDataIBM(t, tc, tb, front, front2=front2, dimPb=dimPb, frontType=frontType, IBCType=IBCType, depth=depth,
                      Reynolds=Reynolds, yplus=yplus, Lref=Lref,
                      cartesian=cartesian, twoFronts=twoFronts, check=check, optimized=optimized, nature=nature_loc, penalty=1, extrap=extrap_loc, val=val,
                      tbFilament=tbFilament, frontWMM=frontWMM)
    if verbose: C_IBM.printTimeAndMemory__('compute interpolation data (IBM)', time=python_time.time()-pt0)
    #===================
    # STEP 6 : INIT IBM
    #===================
    if verbose: pt0 = python_time.time(); C_IBM.printTimeAndMemory__('initialize and clean', time=-1)


    t, tc, tc2 = C_IBM.initializeIBM(t, tc, tb, tinit=tinit, tbCurvi=tbCurvi, dimPb=dimPb, twoFronts=twoFronts,
                                     tbFilament=tbFilament, cleanCellN=cleanCellN)

    #C_IBM._redispatch__(t=t, tc=tc, tc2=tc2)

    C_IBM._setInjOutlet__(tc, tb)

    if isinstance(tc_out, str):
        if cartesian: tcp = Compressor.compressCartesian(tc)
        else: tcp = tc
        Cmpi.convertPyTree2File(tcp, tc_out, ignoreProcNodes=True)

        if tc2:
            if cartesian: tcp2 = Compressor.compressCartesian(tc2)
            else: tcp2 = tc2
            tc2_out = tc_out.replace('tc', 'tc2') if 'tc' in tc_out else 'tc2.cgns'
            Cmpi.convertPyTree2File(tcp2, tc2_out, ignoreProcNodes=True)

    if isinstance(t_out, str):
        if cartesian: tp = Compressor.compressCartesian(t)
        else: tp = t
        Cmpi.convertPyTree2File(tp, t_out, ignoreProcNodes=True)

    C_IBM._computeMeshInfo(t)

    Cmpi.barrier()
    if verbose: C_IBM.printTimeAndMemory__('initialize and clean', time=python_time.time()-pt0)

    if tc2 is not None: return t, tc, tc2
    else: return t, tc

def prepareIBMDataExtrude(t_case, t_out, tc_out, t, to=None,
                          depth=2, frontType=1, octreeMode=0, IBCType=1, nature=1, order=2, optimized=-1, extrap=1,
                          verbose=True, check=False, twoFronts=False, cartesian=True, frontFile=None, cleanCellN=False,
                          yplus=100., Lref=1., correctionMultiCorpsF42=False, blankingF42=False, wallAdaptF42=None, heightMaxF42=-1.,
                          tbox=None, extrusion='cart'):
    import Generator.IBM as G_IBM
    import time as python_time

    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = Internal.copyTree(t_case)

    refstate = Internal.getNodeFromName(tb, 'ReferenceState')
    flowEqn  = Internal.getNodeFromName(tb, 'FlowEquationSet')

    Reynolds = Internal.getNodeFromName(tb, 'Reynolds')
    if Reynolds is not None:
        Reynolds = Internal.getValue(Reynolds)
        if Reynolds < 1.e5: frontType = 1
    else: Reynolds = 1.e6

    expand = 3 if frontType != 42 else 4

    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    if dimPb is None: raise ValueError('prepareIBMDataPara: EquationDimension is missing in input geometry tree.')
    dimPb = Internal.getValue(dimPb)
    if dimPb == 2: C._initVars(tb, 'CoordinateZ', 0.)

    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('prepareIBMDataPara: GoverningEquations is missing in input geometry tree.')
    model = Internal.getValue(model)

    ibctypes = Internal.getNodesFromName(tb, 'ibctype')
    if ibctypes is None: raise ValueError('prepareIBMDataPara: ibc type is missing in input geometry tree.')
    ibctypes = list(set(Internal.getValue(ibc) for ibc in ibctypes))

    if model == 'Euler':
        if any(ibc in ['Musker', 'MuskerMob', 'Mafzal', 'Log', 'TBLE', 'TBLE_FULL'] for ibc in ibctypes):
            raise ValueError("prepareIBMDataPara: governing equations (Euler) not consistent with ibc types %s"%(ibctypes))

    #===================
    # STEP 0 : GET FILAMENT BODIES
    #===================
    tb, tbFilament              = D_IBM.determineClosedSolidFilament__(tb)
    isFilamentOnly, isWireModel = D_IBM.localWMMFlags__(tb, tbFilament)

    if extrusion=='cyl': cartesian = False
    cellN = Internal.getNodeFromName(t, 'cellN')
    if cellN is not None: C._initVars(t, '{centers:cellN}=minimum(1, {centers:cellN})') # modification needed for extrude

    #===================
    # STEP 1 : GENERATE MESH
    #===================
    ##SKIPPED - mesh is provided as an input

    #===================
    # STEP 2 : DIST2WALL
    #===================
    ##SKIPPED - mesh is provided as an input & has TurbulentDistance already

    #===================
    # STEP 3 : BLANKING IBM
    #===================
    if verbose: pt0 = python_time.time(); C_IBM.printTimeAndMemory__('blank by IBC bodies', time=-1, functionName='prepareIBMDataExtrude')
    C_IBM._blankingIBM(t, tb, dimPb=dimPb, frontType=frontType, IBCType=IBCType, depth=depth,
                       Reynolds=Reynolds, yplus=yplus, Lref=Lref, twoFronts=twoFronts, cleanCellN=cleanCellN,
                       heightMaxF42=heightMaxF42, correctionMultiCorpsF42=correctionMultiCorpsF42,
                       wallAdaptF42=wallAdaptF42, blankingF42=blankingF42,
                       tbFilament=tbFilament)

    ##set the kmin et kmax Ghost cells are potential donors                                          #__
    listvars_local =['cellNChim','cellNIBC']                                                         #  |
    for z in Internal.getZones(t):                                                                   #  |
        sol            = Internal.getNodeFromName(z,'FlowSolution#Centers')                          #  |
        for var in listvars_local:                                                                   #  |
            cellN          = Internal.getNodeFromName(sol,var)[1]                                    #  | Modification needed for extrude &
            sh             = numpy.shape(cellN)                                                      #  | to replicate previous behavior
            for k in [0,1, sh[2]-2, sh[2]-1]:                                                        #  |
                for j in range(sh[1]):                                                               #  |
                    for i in range(sh[0]):                                                           #  |
                        if  cellN[i,j,k] != 0:  cellN[i,j,k] =1                                      #  |

    C._initVars(t,'{centers:cellN}=maximum(0.,{centers:cellNChim})')# vaut -3, 0, 1, 2 initialement  #__

    #C.convertPyTree2File(t,'verifBlank.cgns')


    Cmpi.barrier()
    #C_IBM._redispatch__(t=t)
    if verbose: C_IBM.printTimeAndMemory__('blank by IBC bodies', time=python_time.time()-pt0, functionName='prepareIBMDataExtrude')
    #===================
    # STEP 4 : INTERP DATA CHIM
    #===================
    ## REQUIREMENT:: cellN mush be correct here      --> _setInterpData uses cellN
    ##               if cellN* is correct henceforth --> correct values at the end of prepareIBMDataExtrude
    #mise a zero cellN dans Ghost BC
    for z in Internal.getZones(t):
        bcs   = Internal.getNodesFromType2(z, 'BC_t')
        sol   = Internal.getNodeFromName1(z, "FlowSolution#Centers")
        cellN = Internal.getNodeFromName1(sol,"cellNChim")[1]
        for bc in bcs:
            btype = Internal.getValue(bc)
            if btype != "BCOverlap":
                ptrange = Internal.getNodesFromType1(bc, 'IndexRange_t')
                rg      = ptrange[0][1]
                if rg[0,1]==rg[0,0]:
                    if rg[0,1]==1: cellN[0:2,:,:]=1
                    else         : cellN[-2:,:,:]=1
                elif rg[1,1]==rg[1,0]:
                    if rg[1,1]==1: cellN[:, 0:2,:]=1
                    else         : cellN[:, -2:,:]=1
                elif rg[2,1]==rg[2,0] and dimPb==3:
                    if rg[2,1]==1: cellN[:,:,0:2]=1
                    else         : cellN[:,:,2: ]=1

    if verbose: pt0 = python_time.time(); C_IBM.printTimeAndMemory__('compute interpolation data (Abutting & Chimera)', time=-1, functionName='prepareIBMDataExtrude')
    tc = C.node2Center(t)

    if Internal.getNodeFromType(t, "GridConnectivity1to1_t") is not None:
        Xmpi._setInterpData(t, tc, nature=1, loc='centers', storage='inverse', sameName=1, dim=dimPb, itype='abutting', order=2, cartesian=cartesian)

    setInterpDataAndSetInterpTransfer__(t,tc, nature=nature, loc='centers', storage='inverse', sameName=1, sameBase=1, dim=dimPb,\
                                        order=order, extrap=extrap, cartesian=cartesian, corner=True)

    if verbose: C_IBM.printTimeAndMemory__('compute interpolation data (Abutting & Chimera)', time=python_time.time()-pt0, functionName='prepareIBMDataExtrude')
    #===================
    # STEP 5 : BUILD FRONT
    #===================
    if verbose: pt0 = python_time.time(); C_IBM.printTimeAndMemory__('build IBM front Extrude', time=-1, functionName='prepareIBMDataExtrude')



    #C.convertPyTree2File(t, 'AvtFront.cgns')

    t, tc, front, front2, frontWMM = buildFrontIBM(t, tc, tb=tb, dimPb=dimPb, frontType=frontType,
                                                   cartesian=cartesian, twoFronts=twoFronts, check=check,
                                                   tbFilament=tbFilament, optimized=optimized)

    #C.convertPyTree2File(t, 'AprFront.cgns')
    if frontFile is not None: front = C.convertFile2PyTree(frontFile)

    '''
    for z in Internal.getZones(t):
        fastc._updateNatureForIBMGhost(z, Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__, Internal.__FlowSolutionCenters__)
    '''

    if verbose: C_IBM.printTimeAndMemory__('build IBM front Extrude', time=python_time.time()-pt0)

    # on recalcule les interp chimere en supprimant les coin maintenant que le calcul du front est ok
    setInterpDataAndSetInterpTransfer__(t, tc, nature=nature, loc='centers', storage='inverse', sameName=1, sameBase=1, dim=dimPb,\
                                        extrap=extrap, order=order, cartesian=cartesian, corner=False)

    if verbose: C_IBM.printTimeAndMemory__('build IBM front Extrude', time=python_time.time()-pt0, functionName='prepareIBMDataExtrude')
    #===================
    # STEP 6 : INTERP DATA IBM
    #===================
    if verbose: pt0 = python_time.time(); C_IBM.printTimeAndMemory__('compute interpolation data (IBM)', time=-1, functionName='prepareIBMDataExtrude')
    nature_loc = 1
    extrap_loc = 1
    val        = 1

    ##set the kmin et kmax Ghost cells are potential donor
    '''
    '''
    #on empeche la creation de point IBc en k=-1 et 0  et on rend ces cellules non donneuse
    todo =[['cellNIBC',t,'FlowSolution#Centers', 1.], ['cellN',tc,'FlowSolution', 1.25]]
    for task in todo:
        val=task[3]
        for z in Internal.getZones(task[1]):
            sol            = Internal.getNodeFromName(z,task[2])
            var = task[0]
            cellN          = Internal.getNodeFromName(sol,var)[1]
            sh             = numpy.shape(cellN)
            for k in [0,1, sh[2]-2, sh[2]-1]:
                for j in range(sh[1]):
                    for i in range(sh[0]):
                        if  cellN[i,j,k] != 0:  cellN[i,j,k] =val

    #C.convertPyTree2File(t, 'AvtIBM.cgns')
    _setInterpDataIBM(t, tc, tb, front, front2=front2, dimPb=dimPb, frontType=frontType, IBCType=IBCType, depth=depth,
                      Reynolds=Reynolds, yplus=yplus, Lref=Lref,
                      cartesian=cartesian, twoFronts=twoFronts, check=check,optimized=optimized, nature=nature_loc, penalty=1, extrap=extrap_loc, val=val,
                      tbFilament=tbFilament, frontWMM=frontWMM)
    #C.convertPyTree2File(t, 'AprIBM.cgns')

    if verbose: C_IBM.printTimeAndMemory__('compute interpolation data (IBM)', time=python_time.time()-pt0, functionName='prepareIBMDataExtrude')
    #===================
    # STEP 7 : INIT IBM
    #===================
    if verbose: pt0 = python_time.time(); C_IBM.printTimeAndMemory__('initialize and clean', time=-1, functionName='prepareIBMDataExtrude')
    tsave = Internal.copyTree(t)    # Modification needed to by pass the initialization of t in the macro function initializeIBM
    t     = None                    #
    t, tc, tc2 = C_IBM.initializeIBM(t, tc, tb, dimPb=dimPb, twoFronts=twoFronts, tbFilament=tbFilament)
    t = Internal.copyTree(tsave)    # Modification needed to by pass the initialization of t in the macro function initializeIBM
    #C_IBM._redispatch__(t=t, tc=tc, tc2=tc2)

    if extrusion == 'cyl':                                                                              #__
        T._cyl2Cart(t, (0,0,0),(1,0,0))                                                                 #  |
        T._cyl2Cart(tc,(0,0,0),(1,0,0))                                                                 #  |
        # modif info mesh in zonesubregion_t                                                            #  |
        for z in Internal.getZones(tc):                                                                 #  |
            for zsr in Internal.getNodesFromType(z, "ZoneSubRegion_t"):                                 #  |
                zsrname = Internal.getName(zsr)                                                         #  |
                zsrname = zsrname.split('_')                                                            #  |Modification needed for extrude in cylindrical coordinates
                if zsrname[0]=='IBCD':                                                                  #  |
                    for var in ['C','W','I']:                                                           #  |
                        r     = Internal.getNodeFromName(zsr,'CoordinateY_P'+var)[1]                    #  |
                        theta = Internal.getNodeFromName(zsr,'CoordinateZ_P'+var)[1]                    #  |
                        for l in range(numpy.size(r)):                                                  #  |
                            yy  = r[l]*numpy.cos( theta[l] )                                            #  |
                            zz  = r[l]*numpy.sin( theta[l] )                                            #  |
                            r[l]= yy; theta[l] = zz                                                     #  |
                            #__

    if isinstance(tc_out, str):
        tcp = Compressor.compressCartesian(tc)
        Cmpi.convertPyTree2File(tcp, tc_out, ignoreProcNodes=True)

        if twoFronts:
            tcp2 = Compressor.compressCartesian(tc2)
            tc2_out = tc_out.replace('tc', 'tc2') if 'tc' in tc_out else 'tc2.cgns'
            Cmpi.convertPyTree2File(tcp2, tc2_out, ignoreProcNodes=True)

    if isinstance(t_out, str):
        tp = Compressor.compressCartesian(t)
        Cmpi.convertPyTree2File(tp, t_out, ignoreProcNodes=True)

    C_IBM._computeMeshInfo(t)

    if Cmpi.size > 1: Cmpi.barrier()
    if verbose: C_IBM.printTimeAndMemory__('initialize and clean', time=python_time.time()-pt0, functionName='prepareIBMDataExtrude')

    if tc2 is not None: return t, tc, tc2
    else: return t, tc

def buildFrontIBM(t, tc, tb=None, dimPb=3, frontType=1, cartesian=True, twoFronts=False, check=False,
                  tbFilament=None, optimized=1):
    """Build the IBM front for IBM pre-processing."""

    isFilamentOnly, isWireModel = D_IBM.localWMMFlags__(tb, tbFilament)

    tbbc = Cmpi.createBBoxTree(tc)

    interpDataType = 0 if cartesian else 1

    C._initVars(t,'{centers:cellNIBCDnr}=minimum(2.,abs({centers:cellNIBC}))')
    #C._initVars(t,'{centers:cellNIBC}=maximum(0.,{centers:cellNIBC})')# vaut -3, 0, 1, 2, 3 initialement
    if optimized==-1:
        npass=2
        #met les valeur 3 (ghost non masquee) a 1.5 (donneuse dans algo opt=-1). la valeur 1.5 est modifiee a 1 avant interpIBM
        for z in Internal.getZones(t):
            tmp = Internal.getNodeFromName(z,'cellNIBC')[1]
            tmp[tmp >= 2.5] = 1.5
            #tmp[tmp >= 2.5] = 1.05
            #tmp[tmp >= 2.5] = 1.02
            tmp[tmp <=-2.5] = 2

        C._cpVars(t,'centers:cellNIBC',t,'centers:cellN')
        Internal._rmNodesByName(tc,'cellN')
        C._initVars(tc,'{cellN}=1.')
        for z in Internal.getZones(tc):
            cellN = Internal.getNodeFromName(z,'cellN')[1]
            z1    = Internal.getNodeFromName(t,z[0])
            sol   = Internal.getNodeFromName(z1,'FlowSolution#Centers')
            cellNc= Internal.getNodeFromName(sol,'cellN')[1]
            sh = numpy.shape(cellNc)
            sh2 = numpy.shape(cellN)
            if   len(sh)==3 and sh[2]!=1: cellN[:,:,:]=cellNc[:,:,:]
            elif len(sh)==3 and sh[2]==1: cellN[:,:]  =cellNc[:,:,0]
            else: cellN[:,:]=cellNc[:,:]

    else:
        C._initVars(t,'{centers:cellNIBC}=maximum(0.,{centers:cellNIBC})')# vaut -3, 0, 1, 2, 3 initialement
        npass = 1
        C._initVars(t,'{centers:cellNIBC}={centers:cellNIBC}*({centers:cellNIBC}<2.5)') #met les valeur 3 (ghost non masquee) a zero (non donneuse dans algo legacy)
        C._cpVars(t,'centers:cellNIBC',t,'centers:cellN')
        C._cpVars(t,'centers:cellN',tc,'cellN')

    for l in range(npass):
        # Transfert du cellNFront
        C._cpVars(t,'centers:cellNFront',tc,'cellNFront')
        # propager cellNVariable='cellNFront'
        Xmpi._setInterpTransfers(t, tc, variables=['cellNFront'], cellNVariable='cellNFront', compact=0)

    #C.convertPyTree2File(t,'t_FrontAprsTrans.cgns')

    if twoFronts:
        C._cpVars(t,'centers:cellNFront_2',tc,'cellNFront_2')
        Xmpi._setInterpTransfers(t, tc, variables=['cellNFront_2'], cellNVariable='cellNFront_2', compact=0)

    if frontType == 2: C_IBM._pushBackImageFront2__(t, tc, tbbc, cartesian=cartesian)

    C._rmVars(t,['centers:cellNFront'])
    if twoFronts: C._rmVars(t,['centers:cellNFront_2', 'centers:cellNIBC_2'])

    C._cpVars(t,'centers:TurbulentDistance',tc,'TurbulentDistance')

    print('Minimum distance: %f.'%C.getMinValue(t,'centers:TurbulentDistance'))
    P._computeGrad2(t, 'centers:TurbulentDistance', ghostCells=True, withCellN=False)

    front = C_IBM.getIBMFront(tc, 'cellNFront', dim=dimPb, frontType=frontType)
    front = C_IBM.gatherFront(front)

    if twoFronts:
        front2 = getIBMFront(tc, 'cellNFront_2', dim=dimPb, frontType=frontType)
        front2 = gatherFront(front2)
    else:
        front2 = None

    if isWireModel:
        frontWMM = C_IBM.getIBMFront(tc, 'cellNFrontFilWMM', dim=dimPb, frontType=frontType)
        frontWMM = C_IBM.gatherFront(frontWMM)
    else:
        frontWMM = None


    #C.convertPyTree2File(t, 'Infront.cgns')


    if optimized==-1:
        for l in range(npass):
            # Transfert du cellNIBC (= cellN a cet endroit)
            C._cpVars(t,'centers:cellN',tc,'cellN')
            #  propager cellNle='cellNFront'
            #Xmpi._setInterpTransfers(t, tc, variables=['cellN'], cellNVariable='cellN', compact=0)  # cellNVarriable ne permet pas de transmettre des valeur=2, qui se transforme en 0
            Xmpi._setInterpTransfers(t, tc, variables=['cellN'], compact=0)  # contournement provisoire car risque cellN= 1.652...  avec interp sur 8 pts
        for z in Internal.getZones(t):
            tmp = Internal.getNodeFromName(z,'cellN')[1]
            sh = numpy.shape(tmp)
            hole= Internal.getNodeFromName(z,'cellNIBC_hole')
            if hole is not None:
                hole=hole[1]
                tmp[ hole == 2. ] = 2.
            #tmp[ tmp<=1.9 and tmp >=2.1 ] = 2.
            for k in range(sh[2]):
                for j in range(sh[1]):
                    for i in range(sh[0]):
                        #if hole[i,j,k] == 2:
                        if tmp[i,j,k] >= 1.99 and tmp[i,j,k] <= 2.01: tmp[i,j,k]=2.
                        #if tmp[i,j,k] >= 1.03 and tmp[i,j,k] <= 2.01: tmp[i,j,k]=2.
            '''
            if z[0]=='Cart.356X0':  C.convertPyTree2File(z, '356Avtpdate.cgns')
            '''

        for z in Internal.getZones(t):
            fastc._updateNatureForIBMGhost(z,
                                           Internal.__GridCoordinates__,
                                           Internal.__FlowSolutionNodes__,
                                           Internal.__FlowSolutionCenters__)

            #if z[0]=='Cart.356X0':  C.convertPyTree2File(z, '356AprUpdate.cgns')

    if check and Cmpi.rank == 0:
        C.convertPyTree2File(front, 'front.cgns')
        if twoFronts:
            C.convertPyTree2File(front2, 'front2.cgns')
        if isWireModel:
            C.convertPyTree2File(frontWMM, 'frontWMM.cgns')

    return t, tc, front, front2, frontWMM

#=========================================================================
# Compute the transfer coefficients and data for IBM pre-processing.
# IN: t (tree): computational tree
# IN: tc (tree): connectivity tree
# IN: tb (tree): geometry tree (IBM bodies)
# IN: front (tree): front of image points
# IN: front2 (tree, optional): front of second image points
# IN: dimPb (2 or 3): problem dimension
# IN: frontType (0,1,2 or 42): type of IBM front
# IN: IBCType (-1 or 1): type of IBM, -1: IB target points are located inside the solid, 1: IB target points are located in the fluid
# IN: depth (int): depth of overlapping regions
# IN: Reynolds (float): Reynolds number (F42)
# IN: yplus (float): estimated yplus at the first computed cells (F42)
# IN: Lref (float): characteristic length of the geometry (F42)
# IN: cartesian (boolean): if True, activates optimized algorithms for Cartesian meshes
# IN: twoFronts (boolean): if True, performs the IBM pre-processing for an additional image point positioned farther away
# OUT: IBCD* zones inside tc
# OUT: (optional) 2_IBCD* zones inside tc
#=========================================================================
def _setInterpDataIBM(t, tc, tb, front, front2=None, dimPb=3, frontType=1, IBCType=1, depth=2, Reynolds=1.e6,
                      yplus=100, Lref=1., cartesian=True, twoFronts=False, check=False, optimized=1, nature=1, penalty=1, extrap=1, val=0,
                      tbFilament=None, frontWMM=None):
    """Compute the transfer coefficients and data for IBM pre-processing."""

    isFilamentOnly, isWireModel = D_IBM.localWMMFlags__(tb, tbFilament)

    for zc in Internal.getZones(tc):
        proc = Cmpi.getProc(zc)
        if proc == -1: Cmpi._setProc(zc, 0)

    tbbc = Cmpi.createBBoxTree(tc)

    interpDataType = 0 if cartesian else 1

    zonesRIBC = []
    for zrcv in Internal.getZones(t):
        if C.getMaxValue(zrcv, 'centers:cellNIBC')==2.:
            zonesRIBC.append(zrcv)

    nbZonesIBC = len(zonesRIBC)
    if nbZonesIBC == 0:
        res = [{},{},{}]
        if twoFronts or isWireModel: res2 = [{},{},{}]
    else:
        tb_local = tb
        if tbFilament:
            if isFilamentOnly:tb_local= tbFilament
            else:             tb_local= Internal.merge([tb,tbFilament])
        res = C_IBM.getAllIBMPoints(zonesRIBC, loc='centers',tb=tb_local, tfront=front, frontType=frontType,
                                    cellNName='cellNIBC', depth=depth, IBCType=IBCType, Reynolds=Reynolds, yplus=yplus, Lref=Lref,
                                    isOrthoFirst=isFilamentOnly, check=check)
        if twoFronts:
            res2 = C_IBM.getAllIBMPoints(zonesRIBC, loc='centers',tb=tb, tfront=front2, frontType=frontType,
                                         cellNName='cellNIBC', depth=depth, IBCType=IBCType, Reynolds=Reynolds, yplus=yplus, Lref=Lref, check=check,
                                         twoFronts=twoFronts)
    # cleaning
    C._rmVars(tc,['cellNChim','cellNIBC','TurbulentDistance','cellNFront'])
    # dans t, il faut cellNChim et cellNIBCDnr pour recalculer le cellN a la fin
    varsRM = ['centers:gradxTurbulentDistance','centers:gradyTurbulentDistance','centers:gradzTurbulentDistance','centers:cellNFront','centers:cellNIBC']
    C._rmVars(t, varsRM)
    front = None
    if twoFronts: front2 = None

    graph = {}; datas = {}
    procDict = Cmpi.getProcDict(tc)

    # graph d'intersection des pts images de ce proc et des zones de tbbc
    zones  = Internal.getZones(tbbc)
    allBBs = []
    dictOfCorrectedPtsByIBCType = res[0]
    dictOfWallPtsByIBCType      = res[1]
    dictOfInterpPtsByIBCType    = res[2]
    interDictIBM={}

    if twoFronts or isWireModel:
        dictOfCorrectedPtsByIBCType2 = res2[0]
        dictOfWallPtsByIBCType2      = res2[1]
        dictOfInterpPtsByIBCType2    = res2[2]
    else:
        dictOfCorrectedPtsByIBCType2={}
        dictOfWallPtsByIBCType2     ={}
        dictOfInterpPtsByIBCType2   ={}
    interDictIBM2={}

    if dictOfCorrectedPtsByIBCType!={}:
        for ibcTypeL in dictOfCorrectedPtsByIBCType:
            allCorrectedPts = dictOfCorrectedPtsByIBCType[ibcTypeL]
            allWallPts      = dictOfWallPtsByIBCType[ibcTypeL]
            allInterpPts    = dictOfInterpPtsByIBCType[ibcTypeL]
            for nozr in range(nbZonesIBC):
                if allCorrectedPts[nozr] != []:
                    ##[AJ] Keep temporarily
                    ##if ibcTypeL=="140#filament":
                    ##    nlen  = numpy.shape(allInterpPts[nozr][1])[1]
                    ##    save2file = numpy.zeros((nlen,2),dtype=float)
                    ##    save2file[:,0]=allInterpPts[nozr][1][0][:]
                    ##    save2file[:,1]=allInterpPts[nozr][1][1][:]
                    ##    numpy.savetxt('allInterpPts.txt', save2file, delimiter=',')   # X is an array
                    zrname = zonesRIBC[nozr][0]
                    interpPtsBB = Generator.BB(allInterpPts[nozr])
                    for z in zones:
                        bba = C.getFields('GridCoordinates', z, api=1)[0]
                        if Generator.bboxIntersection(interpPtsBB, bba, isBB=True):
                            zname = z[0]
                            popp  = Cmpi.getProc(z)
                            if Cmpi.size > 1: Distributed.updateGraph__(graph, popp, Cmpi.rank, zname)

                            if zrname not in interDictIBM: interDictIBM[zrname]=[zname]
                            else:
                                if zname not in interDictIBM[zrname]: interDictIBM[zrname].append(zname)
        if twoFronts or isWireModel:
            for ibcTypeL in dictOfCorrectedPtsByIBCType2:
                allCorrectedPts2 = dictOfCorrectedPtsByIBCType2[ibcTypeL]
                allWallPts2      = dictOfWallPtsByIBCType2[ibcTypeL]
                allInterpPts2    = dictOfInterpPtsByIBCType2[ibcTypeL]
                for nozr in range(nbZonesIBC):
                    if allCorrectedPts2[nozr] != []:
                        ##[AJ] Keep temporarily
                        ##if ibcTypeL=="140#filament":
                        ##    nlen  = numpy.shape(allInterpPts2[nozr][1])[1]
                        ##    save2file = numpy.zeros((nlen,2),dtype=float)
                        ##    save2file[:,0]=allInterpPts2[nozr][1][0][:]
                        ##    save2file[:,1]=allInterpPts2[nozr][1][1][:]
                        ##    numpy.savetxt('allInterpPts2.txt', save2file, delimiter=',')   # X is an array
                        zrname = zonesRIBC[nozr][0]
                        interpPtsBB2 = Generator.BB(allInterpPts2[nozr])
                        for z in zones:
                            bba = C.getFields('GridCoordinates', z, api=1)[0]
                            if Generator.bboxIntersection(interpPtsBB2,bba,isBB=True):
                                zname = z[0]
                                popp  = Cmpi.getProc(z)
                                if Cmpi.size > 1: Distributed.updateGraph__(graph, popp, Cmpi.rank, zname)
                                if zrname not in interDictIBM2: interDictIBM2[zrname]=[zname]
                                else:
                                    if zname not in interDictIBM2[zrname]: interDictIBM2[zrname].append(zname)
    else: graph={}
    if Cmpi.KCOMM is not None: allGraph = Cmpi.KCOMM.allgather(graph)
    else: allGraph = [graph]

    graph = {}
    for i in allGraph:
        for k in i:
            if not k in graph: graph[k] = {}
            for j in i[k]:
                if not j in graph[k]: graph[k][j] = []
                graph[k][j] += i[k][j]
                graph[k][j] = list(set(graph[k][j])) # pas utile?

    # keyword subr=False to avoid memory overflow
    Cmpi._addXZones(tc, graph, variables=['cellN'], cartesian=cartesian, subr=False)

    ReferenceState = Internal.getNodeFromType2(t, 'ReferenceState_t')
    model = Internal.getNodeFromName(t, 'GoverningEquations')
    if model is not None: model = Internal.getValue(model)
    else:                 model = "Euler"

    for i in range(Cmpi.size): datas[i] = [] # force

    if dictOfCorrectedPtsByIBCType!={}:
        for ibcTypeL in dictOfCorrectedPtsByIBCType:
            if '#' in ibcTypeL: ibcNameL = 'IBCD_'+'_'.join(ibcTypeL.split('#')) #ibctype with familyname
            else: ibcNameL = 'IBCD_'+ibcTypeL #regular ibctype
            allCorrectedPts = dictOfCorrectedPtsByIBCType[ibcTypeL]
            allWallPts = dictOfWallPtsByIBCType[ibcTypeL]
            allInterpPts = dictOfInterpPtsByIBCType[ibcTypeL]
            for nozr in range(nbZonesIBC):
                if allCorrectedPts[nozr] != []:
                    zrcv = zonesRIBC[nozr]
                    zrname = zrcv[0]
                    dnrZones = []
                    for zdname in interDictIBM[zrname]:
                        zd = Internal.copyRef(Internal.getNodeFromName2(tc, zdname))
                        if zd is None: print('!!!Zone None', zrname, zdname)
                        else: dnrZones.append(zd)

                    if optimized ==-1 : # correction ghost non masque pour les rendre donneuse: val=1
                        valRef=1.5
                        for z in dnrZones:
                            cellN = Internal.getNodeFromName(z,'cellN')[1]
                            #print("shape", z[0], z[1], numpy.shape(cellN))
                            C._initVars(z,'{tmp}={cellN}')
                            sh_R    = numpy.shape(cellN)
                            ni =sh_R[0]; nj =sh_R[1];nk =sh_R[2];
                            for k in range(nk):
                                for j in range(nj):
                                    for i in range(2):
                                        if abs(cellN[i     , j, k]-valRef) < 0.01: cellN[i,j,k       ]= val
                                        if abs(cellN[ni-1-i, j, k]-valRef) < 0.01: cellN[ni-1-i, j, k]= val
                            for k in range(nk):
                                for j in range(2):
                                    for i in range(ni):
                                        if abs(cellN[i     , j     , k     ]-valRef) < 0.01: cellN[i,j     , k]= val
                                        if abs(cellN[i     , nj-1-j, k     ]-valRef) < 0.01: cellN[i,nj-1-j, k]= val
                            if dimPb ==3:
                                for k in range(2):
                                    for j in range(nj):
                                        for i in range(ni):
                                            if abs(cellN[i, j, k     ]-valRef) < 0.01: cellN[i, j , k     ]= val
                                            if abs(cellN[i, j, nk-1-k]-valRef) < 0.01: cellN[i, j , nk-1-k]= val

                    '''
                    if zrname in ['Cart.359X0','Cart.396X0'] :
                      C.convertPyTree2File(zrcv,'RecIBM_'+zrname+'.cgns')
                      C.convertPyTree2File(dnrZones,'DnrIBM_'+zrname+'.cgns')
                      import pickle
                      with open ("Interp_"+zrname+".lst", "wb" ) as inter:
                         pickle.dump ( allInterpPts[nozr]  , inter )
                      with open ("wall_"+zrname+".lst", "wb" ) as wall:
                         pickle.dump ( allWallPts[nozr]  , wall )
                      with open ("Corrected_"+zrname+".lst", "wb" ) as corrected:
                         pickle.dump ( allCorrectedPts[nozr]  , corrected )
                    '''
                    XOD._setIBCDataForZone__(zrcv, dnrZones, allCorrectedPts[nozr], allWallPts[nozr], allInterpPts[nozr],
                                             nature=nature, penalty=penalty, extrap=extrap, loc='centers', storage='inverse', dim=dimPb,
                                             interpDataType=interpDataType, ReferenceState=ReferenceState, bcType=ibcTypeL,model=model)

                    if optimized==-1 :
                        C._initVars(dnrZones,'{cellN}={tmp}')
                        Internal.rmNodesByName(dnrZones,'tmp')

                    nozr += 1
                    for zd in dnrZones:
                        zdname = zd[0]
                        destProc = procDict[zdname]

                        IDs = []
                        for i in zd[2]:
                            if i[0][0:4] == 'IBCD':
                                #where i[0] == 'IBCD_#ibctype_#familyname_#zname' or 'IBCD_#ibctype_#zname'
                                if ibcNameL == '_'.join(i[0].split('_')[:-1]) and Internal.getValue(i)==zrname:
                                    IDs.append(i) #add only ibcd zones related to local ibcType to avoid doublons

                        if IDs != []:
                            if destProc == Cmpi.rank:
                                zD = Internal.getNodeFromName2(tc, zdname)
                                zD[2] += IDs
                            else:
                                if destProc not in datas: datas[destProc]=[[zdname,IDs]]
                                else: datas[destProc].append([zdname,IDs])
                        else:
                            if destProc not in datas: datas[destProc] = []

    if dictOfCorrectedPtsByIBCType2!={}:
        for ibcTypeL in dictOfCorrectedPtsByIBCType2:
            if '#' in ibcTypeL: ibcNameL = '2_IBCD_'+'_'.join(ibcTypeL.split('#'))
            else: ibcNameL = '2_IBCD_'+ibcTypeL
            allCorrectedPts2 = dictOfCorrectedPtsByIBCType2[ibcTypeL]
            allWallPts2      = dictOfWallPtsByIBCType2[ibcTypeL]
            allInterpPts2    = dictOfInterpPtsByIBCType2[ibcTypeL]
            for nozr in range(nbZonesIBC):
                if allCorrectedPts2[nozr] != []:
                    zrcv     = zonesRIBC[nozr]
                    zrname   = zrcv[0]
                    dnrZones = []
                    for zdname in interDictIBM2[zrname]:
                        zd = Internal.copyRef(Internal.getNodeFromName2(tc, zdname))
                        if zd is None: print('!!!Zone None', zrname, zdname)
                        else: dnrZones.append(zd)
                    XOD._setIBCDataForZone2__(zrcv, dnrZones, allCorrectedPts2[nozr], allWallPts2[nozr], None, allInterpPts2[nozr],
                                              nature=1, penalty=1, loc='centers', storage='inverse', dim=dimPb,
                                              interpDataType=interpDataType, ReferenceState=ReferenceState, bcType=ibcTypeL)

                    nozr += 1
                    for zd in dnrZones:
                        zdname = zd[0]
                        destProc = procDict[zdname]

                        IDs = []
                        for i in zd[2]:
                            if i[0][0:6] == '2_IBCD':
                                if ibcNameL == '_'.join(i[0].split('_')[:-1]) and Internal.getValue(i)==zrname:
                                    IDs.append(i) #add only ibcd zones related to local ibcType to avoid doublons

                        if IDs != []:
                            if destProc == Cmpi.rank:
                                zD = Internal.getNodeFromName2(tc, zdname)
                                zD[2] += IDs
                            else:
                                if destProc not in datas: datas[destProc]=[[zdname,IDs]]
                                else: datas[destProc].append([zdname,IDs])
                        else:
                            if destProc not in datas: datas[destProc] = []

    Cmpi._rmXZones(tc)
    dictOfCorrectedPtsByIBCType = None
    dictOfWallPtsByIBCType      = None
    dictOfInterpPtsByIBCType    = None
    interDictIBM = None
    if twoFronts or isWireModel:
        dictOfCorrectedPtsByIBCType2 = None
        dictOfWallPtsByIBCType2      = None
        dictOfInterpPtsByIBCType2    = None
        interDictIBM2 = None

    Internal._rmNodesByName(tc, Internal.__FlowSolutionNodes__)
    destDatas = Cmpi.sendRecv(datas, graph)
    for i in destDatas:
        for n in destDatas[i]:
            zname = n[0]
            IBCDs = n[1]
            if IBCDs != []:
                zD = Internal.getNodeFromName2(tc, zname)
                zD[2] += IBCDs

    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is not None: model = Internal.getValue(model)
    else: model = "NSTurbulent"

    C._initVars(t,'{centers:cellN}=minimum({centers:cellNChim}*{centers:cellNIBCDnr},2.)')
    varsRM = ['centers:cellNChim', 'centers:cellNIBCDnr']
    if model == 'Euler': varsRM += ['centers:TurbulentDistance']
    C._rmVars(t, varsRM)

    if check:
        C_IBM.extractIBMInfo(tc, IBCNames="IBCD_*", fileout='IBMInfo.cgns')

        if twoFronts:
            C_IBM.extractIBMInfo(tc, IBCNames="2_IBCD_*", fileout='IBMInfo2.cgns')

    return None

#=============================================================================
# Version local explicit [GJ]
#=============================================================================
def doInterp2(t, tc, tbb, tb=None, typeI='ID', dim=3, dictOfADT=None, front=None, frontType=0, depth=2, IBCType=1, interpDataType=1, Reynolds=6.e6, yplus=100., Lref=1., check=False):
    ReferenceState = Internal.getNodeFromType2(t, 'ReferenceState_t')
    model = Internal.getNodeFromName(t, 'GoverningEquations')
    if model is not None: model = Internal.getValue(model)
    else:                 model = "Euler"

    bases  = Internal.getNodesFromType1(t     , 'CGNSBase_t')
    dimmm  = Internal.getNodeFromName2(bases[0], 'EquationDimension')
    dimPb   = Internal.getValue(dimmm)
    dxmax = 0.0


    zones = Internal.getZones(t)

    dico_dx = {}
    dico_dy = {}
    dico_dz = {}

    for z in zones:
        nodes = Internal.getNodesFromName(z, 'GridCoordinates')
        coordx = nodes[0][2][0][1]
        coordy = nodes[0][2][1][1]
        coordz = nodes[0][2][2][1]

        dxx  = abs(coordx[1,0,0]   - coordx[0,0,0])
        dyy  = abs(coordy[0,1,0]   - coordy[0,0,0])
        dzz  = abs(coordz[0,0,1]   - coordz[0,0,0])

        dico_dx[z[0]] = dxx
        dico_dy[z[0]] = dyy
        if dimPb == 2:dico_dz[z[0]] = 1
        else : dico_dz[z[0]] = dzz

        if dimPb == 2:dzz=max(dxx,dyy)

        dx = min(dxx,dyy,dzz)
        if dx > dxmax:dxmax=dx

    niveaux_temps = {}
    cx = {}

    for z in zones:
        nodes = Internal.getNodesFromName(z, 'GridCoordinates')
        coordx = nodes[0][2][0][1]
        coordy = nodes[0][2][1][1]
        coordz = nodes[0][2][2][1]

        dxx  = abs(coordx[1,0,0]   - coordx[0,0,0])
        dyy  = abs(coordy[0,1,0]   - coordy[0,0,0])
        dzz  = abs(coordz[0,0,1]   - coordz[0,0,0])

        if dimPb == 2:dzz=max(dxx,dyy)


        dx = min(dxx,dyy,dzz)

        #cx[z[0]]= coordx[1,0,0]

        N = math.log(dxmax/dx)/math.log(2.0)
        N = round(N) - 2
        if N < 0:
            niveaux_temps[z[0]] = 2**0
        else :
            niveaux_temps[z[0]] = 2**N
        ##if (N < 6):N=0
        ##else:N=1
        #if (cx[z[0]] < 0.98): niveaux_temps[z[0]] = 2**N
        #else:  niveaux_temps[z[0]] = 1
        #niveaux_temps['cart.264'] = 2
        #niveaux_temps['cart.294'] = 2

        #niveaux_temps[z[0]] = 2**N

        print(niveaux_temps[z[0]])
        #print(round(dxmax/dx))


    if typeI == 'ID':
        # toutes les zones sont interpolables en Chimere
        intersectionsDict = X.getIntersectingDomains(tbb, method='AABB', taabb=tbb)
        rcvZones = []
        for zrcv in Internal.getZones(t):
            if C.getMaxValue(zrcv,'centers:cellN')==2.:
                zrcvname = zrcv[0]; rcvZones.append(zrcv)
        nozr = 0
        nbZonesChim = len(rcvZones)
        for nozr in range(nbZonesChim):
            zrcv = rcvZones[nozr]
            zrcvname = zrcv[0]
            nozr += 1; hook0 = []
            nobOfDnrBases = []; nobOfDnrZones=[]; dnrZones=[]
            for nobd in range(len(tc[2])):
                if tc[2][nobd][3] == 'CGNSBase_t':
                    for nozd in range(len(tc[2][nobd][2])):
                        zdnr = tc[2][nobd][2][nozd]
                        if zdnr[3] == 'Zone_t':
                            zdnrname = zdnr[0]
                            if zdnrname in intersectionsDict[zrcvname]:
                                nobOfDnrBases.append(nobd)
                                nobOfDnrZones.append(nozd)
                                dnrZones.append(zdnr)
                                if interpDataType==1 and dictOfADT is not None:
                                    hook0.append(dictOfADT[zdnrname])
            if interpDataType == 0: hook0 = None

            X._setInterpData(zrcv, dnrZones, nature=1,penalty=1,loc='centers',storage='inverse',sameName=1,\
                             interpDataType=interpDataType, itype='chimera')


            levelrcv = niveaux_temps[zrcv[0]]


            for nod in range(len(dnrZones)):
                dim__ = Internal.getZoneDim(dnrZones[nod])
                prange = numpy.zeros(6,dtype=Internal.E_NpyInt)
                prangedonor = numpy.zeros(6,dtype=Internal.E_NpyInt)
                profondeur=numpy.zeros(1,dtype=Internal.E_NpyInt)
                dirD=numpy.zeros(1,dtype=Internal.E_NpyInt)
                dirR=numpy.zeros(1,dtype=Internal.E_NpyInt)

                plist = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][2][1]
                plistdnr = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][3][1]
                coeff = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][4][1]
                typ = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][5][1]

                leveldnr = niveaux_temps[dnrZones[nod][0]]

                nobd = nobOfDnrBases[nod]
                nozd = nobOfDnrZones[nod]
                tc[2][nobd][2][nozd] = dnrZones[nod]

                prangebis=numpy.reshape(prange,6)

                info = dnrZones[nod][2][len(dnrZones[nod][2])-1]
                info[2].append(['PointRange', prangebis , [], 'IndexArray_t'])

                transfo=numpy.zeros(3,dtype=Internal.E_NpyInt)#XOD.getTransfo(dnrZones[nod],zrcv)

                Connector.connector.indiceToCoord2(plist,prangedonor,transfo,profondeur,dirD,typ,dirR,plist.size,dim__[1],dim__[2],dim__[3])


                #connector.correctCoeffList(plist, coeff, typ, plist.size , dim__[1] , dim__[2] , dim__[3])

                NMratio = numpy.zeros(3,dtype=Internal.E_NpyInt)
                NMratio[0]=1
                NMratio[1]=1
                NMratio[2]=1

                info[2].append(['PointRangeDonor', prangedonor , [], 'IndexArray_t'])
                info[2].append(['DirReceveur', dirR , [], 'IndexArray_t'])
                info[2].append(['DirDonneur', dirD , [], 'IndexArray_t'])
                info[2].append(['Transform', transfo , [], 'IndexArray_t'])
                info[2].append(['Profondeur', profondeur , [], 'IndexArray_t'])
                info[2].append(['PointPivot', transfo , [], 'IndexArray_t'])
                info[2].append(['NMratio', NMratio , [], 'IndexArray_t'])
                info[2].append(['LevelZRcv', levelrcv , [], 'IndexArray_t'])
                info[2].append(['LevelZDnr', leveldnr , [], 'IndexArray_t'])



    elif typeI == 'IBCD':
        # detection des zones IBC
        zonesRIBC = []
        for zrcv in Internal.getZones(t):
            if C.getMaxValue(zrcv,'centers:cellNIBC')==2.:
                zrcvname = zrcv[0]; zonesRIBC.append(zrcv)

        if zonesRIBC == []: return tc

        res = C_IBM.getAllIBMPoints(zonesRIBC, loc='centers',tb=tb, tfront=front, frontType=frontType, \
                                    cellNName='cellNIBC', depth=depth, IBCType=IBCType, Reynolds=Reynolds, yplus=yplus, Lref=Lref, check=check)
        nbZonesIBC = len(zonesRIBC)
        dictOfADT = {}
        dictOfCorrectedPtsByIBCType = res[0]
        dictOfWallPtsByIBCType = res[1]
        dictOfInterpPtsByIBCType = res[2]
        for ibcTypeL in  dictOfCorrectedPtsByIBCType:
            allCorrectedPts = dictOfCorrectedPtsByIBCType[ibcTypeL]
            allWallPts = dictOfWallPtsByIBCType[ibcTypeL]
            allInterpPts = dictOfInterpPtsByIBCType[ibcTypeL]
            for nozr in range(nbZonesIBC):
                if allCorrectedPts[nozr] != []:
                    interpPtsBB=Generator.BB(allInterpPts[nozr])
                    zrcv = zonesRIBC[nozr]
                    zrcvname = zrcv[0]
                    nobOfDnrBases = []; nobOfDnrZones=[]; dnrZones=[]
                    if interpDataType == 1: hook0 = []
                    else: hook0 = None
                    for nobd in range(len(tc[2])):
                        if tc[2][nobd][3] == 'CGNSBase_t':
                            for nozd in range(len(tc[2][nobd][2])):
                                zdnr = tc[2][nobd][2][nozd]
                                if zdnr[3] == 'Zone_t':
                                    zdnrname = zdnr[0]
                                    zbb = tbb[2][nobd][2][nozd]
                                    bba = C.getFields(Internal.__GridCoordinates__, zbb, api=1)[0]
                                    if Generator.bboxIntersection(interpPtsBB,bba,isBB=True) == 1:
                                        if interpDataType == 1:
                                            if zdnrname not in dictOfADT:
                                                HOOKADT = C.createHook(zdnr, 'adt')
                                                dictOfADT[zdnrname] = HOOKADT
                                            hook0.append(dictOfADT[zdnrname])

                                        dnrZones.append(zdnr)
                                        nobOfDnrBases.append(nobd)
                                        nobOfDnrZones.append(nozd)

                    XOD._setIBCDataForZone__(zrcv, dnrZones, allCorrectedPts[nozr], allWallPts[nozr], allInterpPts[nozr], \
                                             nature=1, penalty=1, loc='centers', storage='inverse',
                                             interpDataType=interpDataType, hook=hook0, dim=dim, \
                                             ReferenceState=ReferenceState, bcType=ibcTypeL,model=model)

                    nozr += 1

                    levelrcv = niveaux_temps[zrcv[0]]

                    for nod in range(len(dnrZones)):

                        dim__ = Internal.getZoneDim(dnrZones[nod])
                        prange = numpy.zeros(6,dtype=Internal.E_NpyInt)
                        prangedonor = numpy.zeros(6,dtype=Internal.E_NpyInt)
                        profondeur=numpy.zeros(1,dtype=Internal.E_NpyInt)
                        dirD=numpy.zeros(1,dtype=Internal.E_NpyInt)
                        dirR=numpy.zeros(1,dtype=Internal.E_NpyInt)

                        plist = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][2][1]
                        plistdnr = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][3][1]
                        coeff = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][4][1]
                        typ = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][5][1]

                        leveldnr = niveaux_temps[dnrZones[nod][0]]

                        nobd = nobOfDnrBases[nod]
                        nozd = nobOfDnrZones[nod]

                        tc[2][nobd][2][nozd] = dnrZones[nod]


                        prangebis=numpy.reshape(prange,6)

                        info = dnrZones[nod][2][len(dnrZones[nod][2])-1]
                        info[2].append(['PointRange', prangebis , [], 'IndexArray_t'])

                        transfo=numpy.zeros(3,dtype=Internal.E_NpyInt)#XOD.getTransfo(dnrZones[nod],zrcv)

                        Connector.connector.indiceToCoord2(plist,prangedonor,transfo,profondeur,dirD,typ,dirR,plist.size,dim__[1],dim__[2],dim__[3])

                        #connector.correctCoeffList(plist, coeff, typ, plist.size , dim__[1] , dim__[2] , dim__[3])

                        NMratio = numpy.zeros(3,dtype=Internal.E_NpyInt)
                        NMratio[0]=1
                        NMratio[1]=1
                        NMratio[2]=1

                        info[2].append(['PointRangeDonor', prangedonor , [], 'IndexArray_t'])
                        info[2].append(['DirReceveur', dirR , [], 'IndexArray_t'])
                        info[2].append(['DirDonneur', dirD , [], 'IndexArray_t'])
                        info[2].append(['Transform', transfo , [], 'IndexArray_t'])
                        info[2].append(['Profondeur', profondeur , [], 'IndexArray_t'])
                        info[2].append(['PointPivot', transfo , [], 'IndexArray_t'])
                        info[2].append(['NMratio', NMratio , [], 'IndexArray_t'])
                        info[2].append(['LevelZRcv', levelrcv , [], 'IndexArray_t'])
                        info[2].append(['LevelZDnr', leveldnr , [], 'IndexArray_t'])


        if dictOfADT is not None:
            for dnrname in dictOfADT: C.freeHook(dictOfADT[dnrname])

    return tc

def prepareIBMData2(t, tbody, DEPTH=2, loc='centers', frontType=1, inv=False, interpDataType=1):
    tb =  Internal.copyRef(tbody)

    # tb: fournit model et dimension
    dimPb = Internal.getNodeFromName(tb,'EquationDimension')
    if dimPb is None: raise ValueError('prepareIBMData: EquationDimension is missing in input body tree.')
    dimPb = Internal.getValue(dimPb)

    # type de traitement paroi: pts interieurs ou externes
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('prepareIBMData: GoverningEquations is missing in input body tree.')
    # model: Euler, NSLaminar, NSTurbulent
    model = Internal.getValue(model)

    if model == 'Euler': IBCType =-1
    else: IBCType = 1 # Points cibles externes
    if loc == 'nodes':
        raise NotImplemented("prepareIBMData: prepareIBMData at nodes not yet implemented.")

    #------------------------
    # Ghost cells (overlaps)
    #------------------------
    X._applyBCOverlaps(t, depth=DEPTH,loc='centers',val=2, cellNName='cellN')
    C._initVars(t,'{centers:cellNChim}={centers:cellN}')

    #------------------------
    # Blanking IBM
    #------------------------
    C._initVars(t,'centers:cellN',1.)
    if dimPb == 2:
        z0 = Internal.getNodeFromType2(t, 'Zone_t')
        dims = Internal.getZoneDim(z0)
        npts = dims[1]*dims[2]*dims[3]
        zmin = C.getValue(z0,'CoordinateZ',0)
        zmax = C.getValue(z0,'CoordinateZ',npts-1)
        dz = zmax-zmin
        # Creation du corps 2D pour le preprocessing IBC
        T._addkplane(tb)
        T._contract(tb, (0,0,0), (1,0,0), (0,1,0), dz)

    t = blankByIBCBodies(t,tb,'centers',dimPb)
    if not inv: C._initVars(t,'{centers:cellNIBC}={centers:cellN}')
    if inv: C._initVars(t,'{centers:cellNIBC}=1-{centers:cellN}') # ecoulement interne



    #-----------------------------------------
    # calcul de la normale et distance signee
    #-----------------------------------------
    COMPDIST = False # distance deja calculee ou non
    if Internal.getNodeFromName(t, 'TurbulentDistance') is None: COMPDIST=True
    if COMPDIST:
        print('Computing distance field...')
        DTW._distance2Walls(t,tb,loc='centers',type='ortho',signed=0)
    else: pass
    _signDistance(t)

    #-----------------------------------------
    # Pts IBC
    #-----------------------------------------
    C._initVars(t,'{centers:cellN}={centers:cellNIBC}')
    # determination des pts IBC
    if IBCType == -1: X._setHoleInterpolatedPoints(t,depth=-DEPTH,dir=0,loc='centers',cellNName='cellN',addGC=False)
    elif IBCType == 1:
        X._setHoleInterpolatedPoints(t,depth=1,dir=1,loc='centers',cellNName='cellN',addGC=False) # pour les gradients
        if frontType < 2:
            X._setHoleInterpolatedPoints(t,depth=DEPTH,dir=0,loc='centers',cellNName='cellN',addGC=False)
        else:
            DEPTHL=DEPTH+1
            X._setHoleInterpolatedPoints(t,depth=DEPTHL,dir=0, loc='centers',cellNName='cellN',addGC=False)
            #cree des pts extrapoles supplementaires
            # _blankClosestTargetCells(t,cellNName='cellN', depth=DEPTHL)
    else:
        raise ValueError('prepareIBMData: not valid IBCType. Check model.')
    _removeBlankedGrids(t, loc='centers')
    print('Nb of Cartesian grids=%d.'%len(Internal.getZones(t)))
    npts = 0
    for i in Internal.getZones(t):
        dims = Internal.getZoneDim(i)
        npts += dims[1]*dims[2]*dims[3]
    print('Final number of points=%5.4f millions.'%(npts/1000000.))

    C._initVars(t,'{centers:cellNIBC}={centers:cellN}')

    #------------------------------------------------------------------------
    # Nature des points en fonction de leur nature Chimere et leur nature IBC
    #------------------------------------------------------------------------
    # -3 : agit comme un point masque - non donneur pour le type de point
    #  3  : agit comme donneur uniquement
    # updateNatureForIBM: modifie cellNChim, cellNFront, cellNIBM
    # cellNChim=-3, si cellNIBC=0 (masque)
    if IBCType == 1: # Points corriges IBM externes
        C._initVars(t,'{centers:cellNFront}=logical_and({centers:cellNIBC}>0.5, {centers:cellNIBC}<1.5)')
        for z in Internal.getZones(t):
            Connector.connector._updateNatureForIBM(z, IBCType,
                                                    Internal.__GridCoordinates__,
                                                    Internal.__FlowSolutionNodes__,
                                                    Internal.__FlowSolutionCenters__)

    else: # EN 2 PARTIES : NECESSITE LE TRANSFERT DU FRONT PAR INTERPOLATION, QUI EST CALCULEE APRES
        print('Euler: on repousse le front un peu plus loin.')
        C._initVars(t,'{centers:dummy}={centers:cellN}') # sauvegarde
        C._initVars(t,'{centers:cellN}=({centers:cellNIBC}>0.5)*({centers:cellNIBC}<1.5)')
        X._setHoleInterpolatedPoints(t,depth=1,dir=1,loc='centers',cellNName='cellN',addGC=False)
        C._initVars(t,'{centers:cellNFront}=logical_and({centers:cellN}>0.5, {centers:cellN}<1.5)')
        C._cpVars(t,'centers:dummy',t,'centers:cellN')
        C._rmVars(t, ['centers:dummy'])
        for z in Internal.getZones(t):
            Connector.connector._updateNatureForIBM(z, IBCType,
                                                    Internal.__GridCoordinates__,
                                                    Internal.__FlowSolutionNodes__,
                                                    Internal.__FlowSolutionCenters__)
    #------------------------------------------------------------------------
    # setInterpData - Chimere
    C._initVars(t,'{centers:cellN}=maximum(0.,{centers:cellNChim})')# vaut -3, 0, 1, 2 initialement

    # maillage donneur: on MET les pts IBC comme donneurs
    tc = C.node2Center(t)
    FSN = Internal.getNodesFromName(tc, Internal.__FlowSolutionNodes__)
    Internal._rmNodesByName(FSN,'cellNFront')
    Internal._rmNodesByName(FSN,'cellNIBC')
    Internal._rmNodesByName(FSN, "TurbulentDistance")

    tbb = G.BB(tc)

    # Creation du dictionnaire des ADT pour les raccords
    if interpDataType == 1:
        dictOfADT = {}
        for zdnr in Internal.getZones(tc):
            zdnrname = zdnr[0]
            if zdnrname not in dictOfADT:
                HOOKADT = C.createHook(zdnr, 'adt')
                dictOfADT[zdnrname] = HOOKADT
    else: dictOfADT = None
    print('Interpolations Chimere.')
    tc = doInterp2(t, tc, tbb, tb=None, typeI='ID', dim=dimPb,
                   interpDataType=interpDataType, dictOfADT=dictOfADT)
    if dictOfADT is not None:
        for dnrname in dictOfADT: C.freeHook(dictOfADT[dnrname])

    # setIBCData - IBC
    C._initVars(t,'{centers:cellNIBCDnr}=minimum(2.,abs({centers:cellNIBC}))')
    C._initVars(t,'{centers:cellNIBC}=maximum(0.,{centers:cellNIBC})')# vaut -3, 0, 1, 2, 3 initialement
    C._initVars(t,'{centers:cellNIBC}={centers:cellNIBC}*({centers:cellNIBC}<2.5)')
    C._cpVars(t,'centers:cellNIBC',t,'centers:cellN')
    C._cpVars(t,'centers:cellN',tc,'cellN')

    #-----------------------------------------------
    # Transfert du cellNFront
    C._cpVars(t,'centers:cellNFront',tc,'cellNFront')

    for zc in Internal.getZones(tc):
        cellNFront = Internal.getNodeFromName2(zc,'cellNFront')
        if cellNFront != []:
            cellNFront = cellNFront[1]
            sizeTot = cellNFront.shape[0]*cellNFront.shape[1]*cellNFront.shape[2]
            sizeOne =  int(numpy.sum(cellNFront))
            if sizeOne < sizeTot:
                XOD._setInterpTransfers(t,zc,variables=['cellNFront'],cellNVariable='cellNFront',compact=0)

    if frontType==2 or frontType==3: _pushBackImageFront2(t, tc, tbb, interpDataType=interpDataType)

    ## Fin traitement specifique, vaut 0 ou 1 apres la ligne suivante
    C._cpVars(t,'centers:cellNFront',tc,'cellNFront')
    C._rmVars(t,['centers:cellNFront'])
    C._cpVars(t,'centers:TurbulentDistance',tc,'TurbulentDistance')

    print('Minimum distance: %f.'%C.getMinValue(t,'centers:TurbulentDistance'))
    P._computeGrad2(t, 'centers:TurbulentDistance', withCellN=False)
    print('Building the IBM front.')
    front = getIBMFront(tc, 'cellNFront', dimPb, frontType)
    print('Interpolations IBM')
    tc = doInterp2(t,tc,tbb, tb=tb,typeI='IBCD',dim=dimPb, dictOfADT=None, front=front, frontType=frontType, depth=DEPTH, IBCType=IBCType, interpDataType=interpDataType)

    # cleaning...
    Internal._rmNodesByName(tc, Internal.__FlowSolutionNodes__)
    Internal._rmNodesByName(tc, Internal.__GridCoordinates__)
    C._initVars(t,'{centers:cellN}=minimum({centers:cellNChim}*{centers:cellNIBCDnr},2.)')
    varsRM = ['centers:gradxTurbulentDistance','centers:gradyTurbulentDistance','centers:gradzTurbulentDistance','centers:cellNFront','centers:cellNIBCDnr']
    varsRM += ['centers:cellNChim','centers:cellNIBC']
    C._rmVars(t, varsRM)
    C._rmVars(tc,['cellNChim','cellNIBC','TurbulentDistance'])
    #----------
    # SORTIE
    #----------
    return t, tc

#=============================================================================
# never called ?
#=============================================================================
def doInterp3(t, tc, tbb, tb=None, typeI='ID', dim=3, dictOfADT=None, frontType=0, depth=2, IBCType=1, interpDataType=1, Reynolds=6.e6, yplus=100., Lref=1., check=False):


    ReferenceState = Internal.getNodeFromType2(t, 'ReferenceState_t')
    model = Internal.getNodeFromName(t, 'GoverningEquations')
    if model is not None: model = Internal.getValue(model)
    else:                 model = "Euler"

    bases  = Internal.getNodesFromType1(t     , 'CGNSBase_t')
    dimmm  = Internal.getNodeFromName2(bases[0], 'EquationDimension')
    dimPb   = Internal.getValue(dimmm)
    dxmax = 0.0


    zones = Internal.getZones(t)

    dico_dx = {}
    dico_dy = {}
    dico_dz = {}

    for z in zones:
        nodes = Internal.getNodesFromName(z, 'GridCoordinates')
        coordx = nodes[0][2][0][1]
        coordy = nodes[0][2][1][1]
        coordz = nodes[0][2][2][1]

        dxx  = abs(coordx[1,0,0]   - coordx[0,0,0])
        dyy  = abs(coordy[0,1,0]   - coordy[0,0,0])
        dzz  = abs(coordz[0,0,1]   - coordz[0,0,0])

        dico_dx[z[0]] = dxx
        dico_dy[z[0]] = dyy
        if dimPb == 2: dico_dz[z[0]] = 1
        else: dico_dz[z[0]] = dzz

        if dimPb == 2:dzz=max(dxx,dyy)

        dx = min(dxx,dyy,dzz)
        if dx > dxmax:dxmax=dx

    niveaux_temps = {}
    cx = {}

    for z in zones:
        nodes = Internal.getNodesFromName(z, 'GridCoordinates')
        coordx = nodes[0][2][0][1]
        coordy = nodes[0][2][1][1]
        coordz = nodes[0][2][2][1]

        dxx  = abs(coordx[1,0,0]   - coordx[0,0,0])
        dyy  = abs(coordy[0,1,0]   - coordy[0,0,0])
        dzz  = abs(coordz[0,0,1]   - coordz[0,0,0])

        if dimPb == 2:dzz=max(dxx,dyy)

        dx = min(dxx,dyy,dzz)

        #cx[z[0]]= coordx[1,0,0]

        N = math.log(dxmax/dx)/math.log(2.0)
        N = round(N)
        if N < 6: N=0
        else: N=1
        #if (cx[z[0]] < 0.98): niveaux_temps[z[0]] = 2**N
        #else:  niveaux_temps[z[0]] = 1
        #niveaux_temps['cart.264'] = 2
        #niveaux_temps['cart.294'] = 2

        niveaux_temps[z[0]] = 2**N


        print(niveaux_temps[z[0]])
        #print(round(dxmax/dx))


    if typeI == 'ID':
        # toutes les zones sont interpolables en Chimere
        intersectionsDict = X.getIntersectingDomains(tbb, method='AABB', taabb=tbb)


        rcvZones = []
        for zrcv in Internal.getZones(t):
            if C.getMaxValue(zrcv,'centers:cellN')==2.:
                zrcvname = zrcv[0]; rcvZones.append(zrcv)

        #dico={}
        #for zrcv in Internal.getZones(t):
                # listofjoins = Internal.getNodesFromType2(zrcv, 'GridConnectivity_t')
                # if listofjoins is not None:
                #    prange_list=[]
                #    dir_list=[]
                #    for join in listofjoins:
                #        prange_ = Internal.getNodeFromName1(join,'PointRange')[1]
                #dirR = CGC.getDirBorderStruct__(prange_,dimPb)
                #dir_list.append(dirR)
                #print 'prange_= ', prange_
                #       for i in range(3):
                #           if prange_[i,1] == prange_[i,0] and prange_[i,1] != 1:
                #               prange_[i,1] =  prange_[i,1]-1
                #               prange_[i,0] =  prange_[i,0]-1
                #           elif prange_[i,1] != prange_[i,0] and prange_[i,1] != 1 :
                #              prange_[i,1] =  prange_[i,1]-1
                #      prange_=numpy.reshape(prange_,6)
                #      prange_list.append(prange_)
                #  dico[zrcv[0]]=prange_list
                #dico[zrcv[0]]=dir_list
                # print prange_, zrcv[0]


        nozr = 0
        nbZonesChim = len(rcvZones)
        for nozr in range(nbZonesChim):
            zrcv = rcvZones[nozr]
            dim_ = Internal.getZoneDim(zrcv)
            zrcvname = zrcv[0]
            nozr += 1; hook0 = []
            nobOfDnrBases = []; nobOfDnrZones=[]; dnrZones=[]
            for nobd in range(len(tc[2])):
                if tc[2][nobd][3] == 'CGNSBase_t':
                    for nozd in range(len(tc[2][nobd][2])):
                        zdnr = tc[2][nobd][2][nozd]
                        if zdnr[3] == 'Zone_t':
                            zdnrname = zdnr[0]
                            if zdnrname in intersectionsDict[zrcvname]:
                                nobOfDnrBases.append(nobd)
                                nobOfDnrZones.append(nozd)
                                dnrZones.append(zdnr)
                                hook0.append(dictOfADT[zdnrname])

            dnrZones = X.setInterpData(zrcv,dnrZones,nature=1,penalty=1,loc='centers',storage='inverse',sameName=1,\
                                       hook=hook0, itype='chimera')


            levelrcv = niveaux_temps[zrcv[0]]

            for nod in range(len(dnrZones)):

                dim__ = Internal.getZoneDim(dnrZones[nod])
                prange = numpy.zeros(6,dtype=Internal.E_NpyInt)
                prangedonor = numpy.zeros(6,dtype=Internal.E_NpyInt)
                profondeur=numpy.zeros(1,dtype=Internal.E_NpyInt)
                dirD=numpy.zeros(1,dtype=Internal.E_NpyInt)
                dirR=numpy.zeros(1,dtype=Internal.E_NpyInt)

                plist = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][2][1]
                plistdnr = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][3][1]
                coeff = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][4][1]
                typ = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][5][1]

                leveldnr = niveaux_temps[dnrZones[nod][0]]

                nobd = nobOfDnrBases[nod]
                nozd = nobOfDnrZones[nod]

                tc[2][nobd][2][nozd] = dnrZones[nod]


                prangebis=numpy.reshape(prange,6)

                info = dnrZones[nod][2][len(dnrZones[nod][2])-1]
                info[2].append(['PointRange', prangebis , [], 'IndexArray_t'])

                transfo = XOD.getTransfo(dnrZones[nod],zrcv)

                Connector.connector.indiceToCoord2(plist,prangedonor,transfo,profondeur,dirD,typ,dirR,plist.size,dim__[1],dim__[2],dim__[3])


                #connector.correctCoeffList(plist, coeff, typ, plist.size , dim__[1] , dim__[2] , dim__[3])

                NMratio = numpy.zeros(3,dtype=Internal.E_NpyInt)
                NMratio[0]=1
                NMratio[1]=1
                NMratio[2]=1

                info[2].append(['PointRangeDonor', prangedonor , [], 'IndexArray_t'])
                info[2].append(['DirReceveur', dirR , [], 'IndexArray_t'])
                info[2].append(['DirDonneur', dirD , [], 'IndexArray_t'])
                info[2].append(['Transform', transfo , [], 'IndexArray_t'])
                info[2].append(['Profondeur', profondeur , [], 'IndexArray_t'])
                info[2].append(['PointPivot', transfo , [], 'IndexArray_t'])
                info[2].append(['NMratio', NMratio , [], 'IndexArray_t'])
                info[2].append(['LevelZRcv', levelrcv , [], 'IndexArray_t'])
                info[2].append(['LevelZDnr', leveldnr , [], 'IndexArray_t'])



    elif typeI == 'IBCD':
        # detection des zones IBC
        zonesRIBC = []
        for zrcv in Internal.getZones(t):
            if C.getMaxValue(zrcv,'centers:cellNIBC')==2.:
                zrcvname = zrcv[0]; zonesRIBC.append(zrcv)

        if zonesRIBC == []: return tc

        print('Building the IBM front.')

        front = getIBMFront(tc, 'cellNFront', dim, frontType)
        # Sortie du front pour debug
        #C.convertPyTree2File(front, 'front.cgns')

        res = C_IBM.getAllIBMPoints(zonesRIBC, loc='centers',tb=tb, tfront=front, frontType=frontType, \
                                    cellNName='cellNIBC', depth=depth, IBCType=IBCType, Reynolds=Reynolds, yplus=yplus, Lref=Lref, check=check)
        nbZonesIBC = len(zonesRIBC)
        dictOfADT = {}
        dictOfCorrectedPtsByIBCType = res[0]
        dictOfWallPtsByIBCType = res[1]
        dictOfInterpPtsByIBCType = res[2]
        for ibcTypeL in dictOfCorrectedPtsByIBCType:
            allCorrectedPts = dictOfCorrectedPtsByIBCType[ibcTypeL]
            allWallPts = dictOfWallPtsByIBCType[ibcTypeL]
            allInterpPts = dictOfInterpPtsByIBCType[ibcTypeL]
            for nozr in range(nbZonesIBC):
                if allCorrectedPts[nozr] != []:
                    interpPtsBB=Generator.BB(allInterpPts[nozr])

                    zrcv = zonesRIBC[nozr]
                    zrcvname = zrcv[0]
                    nobOfDnrBases = []; nobOfDnrZones=[]; dnrZones=[]; hook0 = []
                    for nobd in range(len(tc[2])):
                        if tc[2][nobd][3] == 'CGNSBase_t':
                            for nozd in range(len(tc[2][nobd][2])):
                                zdnr = tc[2][nobd][2][nozd]
                                if zdnr[3] == 'Zone_t':
                                    zdnrname = zdnr[0]
                                    zbb = tbb[2][nobd][2][nozd]
                                    bba = C.getFields(Internal.__GridCoordinates__, zbb, api=1)[0]
                                    if Generator.bboxIntersection(interpPtsBB,bba,isBB=True) == 1:
                                        if zdnrname not in dictOfADT:
                                            HOOKADT = C.createHook(zdnr, 'adt')
                                            dictOfADT[zdnrname] = HOOKADT
                                        dnrZones.append(zdnr)
                                        hook0.append(dictOfADT[zdnrname])
                                        nobOfDnrBases.append(nobd)
                                        nobOfDnrZones.append(nozd)

                    XOD._setIBCDataForZone__(zrcv, dnrZones, allCorrectedPts[nozr], allWallPts[nozr], allInterpPts[nozr], \
                                             loc='centers', storage='inverse',  hook=hook0, dim=dim, ReferenceState=ReferenceState, bcType=ibcTypeL,model=model)
                    nozr += 1

                    levelrcv = niveaux_temps[zrcv[0]]

                    for nod in range(len(dnrZones)):

                        dim__ = Internal.getZoneDim(dnrZones[nod])
                        prange = numpy.zeros(6,dtype=Internal.E_NpyInt)
                        prangedonor = numpy.zeros(6,dtype=Internal.E_NpyInt)
                        profondeur=numpy.zeros(1,dtype=Internal.E_NpyInt)
                        dirD=numpy.zeros(1,dtype=Internal.E_NpyInt)
                        dirR=numpy.zeros(1,dtype=Internal.E_NpyInt)

                        plist = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][2][1]
                        plistdnr = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][3][1]
                        coeff = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][4][1]
                        typ = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][5][1]

                        leveldnr = niveaux_temps[dnrZones[nod][0]]

                        nobd = nobOfDnrBases[nod]
                        nozd = nobOfDnrZones[nod]

                        tc[2][nobd][2][nozd] = dnrZones[nod]


                        prangebis=numpy.reshape(prange,6)

                        info = dnrZones[nod][2][len(dnrZones[nod][2])-1]
                        info[2].append(['PointRange', prangebis , [], 'IndexArray_t'])

                        transfo=XOD.getTransfo(dnrZones[nod],zrcv)

                        Connector.connector.indiceToCoord2(plist,prangedonor,transfo,profondeur,dirD,typ,dirR,plist.size,dim__[1],dim__[2],dim__[3])

                        #connector.correctCoeffList(plist, coeff, typ, plist.size , dim__[1], dim__[2], dim__[3])

                        NMratio = numpy.zeros(3,dtype=Internal.E_NpyInt)
                        NMratio[0]=1
                        NMratio[1]=1
                        NMratio[2]=1

                        info[2].append(['PointRangeDonor', prangedonor , [], 'IndexArray_t'])
                        info[2].append(['DirReceveur', dirR , [], 'IndexArray_t'])
                        info[2].append(['DirDonneur', dirD , [], 'IndexArray_t'])
                        info[2].append(['Transform', transfo , [], 'IndexArray_t'])
                        info[2].append(['Profondeur', profondeur , [], 'IndexArray_t'])
                        info[2].append(['PointPivot', transfo , [], 'IndexArray_t'])
                        info[2].append(['NMratio', NMratio , [], 'IndexArray_t'])
                        info[2].append(['LevelZRcv', levelrcv , [], 'IndexArray_t'])
                        info[2].append(['LevelZDnr', leveldnr , [], 'IndexArray_t'])

                        print('LEVELS= ', levelrcv, leveldnr)


        for dnrname in dictOfADT.keys(): C.freeHook(dictOfADT[dnrname])

        for dnrname in dictOfADT: C.freeHook(dictOfADT[dnrname])

    return tc


##########################################################################
#Interpolation chimere en 7 etapes pour prepIBM si optimized=-1
#limite le recouvrement
##  etape1 : interp ordre2 nature 1 entre grille de meme niveau
##  etape2 : modif cellN=0 pour les points interpole dans etape 1
##  etape3 : interp ordre2 nature 1 entre grille de niveau N (Receveur) et N-1 (Donneur fin)
##  etape4 : modif cellN=0 pour les points interpole dans etape 3
##  etape5 : interp ordre=order nature 1 entre grille de niveau N (Receveur) et  N+1 (Donneur grossier)
##  etape6 : modif cellN=0 pour les points interpole dans etape 5 si nature =0
##  etape7 : interp ordre=order si  nature=0 entre grille de niveau N (Receveur) et  N+1 (Donneur grossier)
##########################################################################
def setInterpDataAndSetInterpTransfer__(t, tc, nature=1, loc='centers', storage='inverse', sameName=1, sameBase=1, dim=3, order=2, extrap=1, cartesian=True, corner=True ):

    # create dictOfModels for adaptRANSLES in OversetData
    dictOfModels = {}
    for b in Internal.getBases(t):
        model_b = Internal.getNodeFromName2(b, 'GoverningEquations')
        if model_b is not None: model_b = Internal.getValue(model_b)
        else: model_b = 'None'
        for z in Internal.getZones(b):
            model = Internal.getNodeFromName2(z, 'GoverningEquations')
            if model is None: model = model_b
            else: model = Internal.getValue(model)
            dictOfModels[z[0]] = [model]

    dictOfModels = Cmpi.allgatherDict(dictOfModels)
    dictOfModels = {key:value[0] for key,value in dictOfModels.items()}

    #
    ## determine dx=dy for each zone & store per zone & levelzone
    #
    levelZone={}
    hmin_loc = 1.e30
    for z in Internal.getZones(t):
        h = abs(C.getValue(z,'CoordinateX',0)-C.getValue(z,'CoordinateX',1))
        levelZone[z[0]] = h
        if h < hmin_loc : hmin_loc = h
    hmin_loc = Cmpi.allgather(hmin_loc)
    if Cmpi.size > 1:
        hmin = 1.e30
        for h in hmin_loc:
            if h < hmin: hmin = h
    else:
        hmin = hmin_loc

    Nlevels = 1
    for i in levelZone:
        levelZone[i]= math.log( int(levelZone[i]/hmin + 0.00000001), 2)
        if levelZone[i] +1  > Nlevels: Nlevels = int(levelZone[i]) +1

    ## partage des info level en mpi
    if Cmpi.size > 1:
        levelZone = Cmpi.allgatherDict2(levelZone)

    #
    #on vire les ghost de type coin pour alleger les calcul NS, puisque Fast fait sa propre sauce
    #
    if not corner:

        print("ON optimise les coins a interpoler")
        #Internal._rmNodesByName(tc,'ID_*')
        #Sauvegarde cellN passe 2
        C._initVars(t, '{centers:cellNSave2}= {centers:cellN}') #cellN pointe sur cellNIBC
        C._initVars(tc, '{cellNSave2}= {cellN}')
        C._initVars(t, '{centers:cellN}= {centers:cellNSave1}') #modifie cellNIBC

        _correctCellNCorner(t, tc, dim=dim, verbose=0)

        C._cpVars(t,'centers:cellN',tc,'cellN')
        C._initVars(t, '{centers:cellNSave3}= {centers:cellN}')

    else:
        #Sauvegarde cellN passe 1
        C._initVars(t, '{centers:cellNSave1}= {centers:cellN}')

    tbbc = Cmpi.createBBoxTree(tc)
    interDict = X.getIntersectingDomains(tbbc, taabb=tbbc)
    procDict = Cmpi.getProcDict(tc)

    # Get baseName for each zone
    baseNames = {}
    for b in Internal.getBases(tbbc):
        for z in Internal.getZones(b): baseNames[z[0]] = b[0]

    # on ne conserve que les intersections inter bases
    if sameBase == 0:
        for i in interDict:
            bi = baseNames[i]
            out = []
            for z in interDict[i]:
                if bi != baseNames[z]: out.append(z)
            interDict[i] = out

    # Perform addXZones on tc
    graph = Cmpi.computeGraph(tbbc, type='bbox', intersectionsDict=interDict, reduction=False)
    Cmpi._addXZones(tc, graph, variables=['cellN'], noCoordinates=False, cartesian=cartesian, zoneGC=False, keepOldNodes=False)
    # serialisation eventuelle
    #graphs = Cmpi.splitGraph(graph)
    #for g in graphs:
    #    Cmpi._addXZones(tc, g, variables=['centers:cellN'], noCoordinates=False,
    #                    cartesian=False, zoneGC=False, keepOldNodes=False)

    # Build hook on local tc zones
    hooks = {};
    for b in Internal.getBases(tc):
        if b[0] == 'CARTESIAN':
            for z in Internal.getZones(b):
                hooks[z[0]] = None # must be None for Cartesian
        else:
            for z in Internal.getZones(b):
                hooks[z[0]] = C.createHook(z, 'adt')

    datas = {}

    ## parcours des zones par niveau en vue traitememnt particulier
    #
    #etape1: calcul interpolation si grille donneuse meme niveau ou plus fine que receuveuse (ordre2, nature1, pas d'extrap)
    #
    cellNModif={}
    for level in range(Nlevels):
        #filtre les zones par niveau de resolution
        zones=[]
        for z in Internal.getZones(t):
            if levelZone[z[0]]==level:
                zones.append(z)
                #print("level=", level,"zone In", z[0])

        print("Interp level=",level, "etape 1")
        for zr in zones:
            zrname = zr[0]
            baseNameRcv = baseNames[zrname]

            dnrZones = []
            for zdname in interDict[zrname]:
                zd = Internal.getNodeFromName2(tc, zdname)
                baseNameDnr = baseNames[zd[0]]
                if levelZone[zd[0]] == level and (sameBase ==1 or baseNameDnr == baseNameRcv): dnrZones.append(zd)

            hookL = []; interpDataTypeL = []
            for z in dnrZones:
                h = hooks[z[0]]
                hookL.append(h)
                if h is None: interpDataTypeL.append(0)
                else: interpDataTypeL.append(1)

            #print("zd=", zd[0],  "donors", dnrZones)
            if dnrZones != []:
                X._setInterpData(zr, dnrZones, nature=1, penalty=1, extrap=0, loc='centers', storage='inverse',
                                 interpDataType=interpDataTypeL, hook=hookL, sameName=1, order=2, itype='chimera', verbose=0)
                fix='Nocorner'
                if corner: fix='corner'
                '''
                if zr[0]=='Cart.553X0':
                      C.convertPyTree2File(zr,'t553_E1_'+fix+'.cgns')
                      C.convertPyTree2File(dnrZones,'tc553_E1_'+fix+'.cgns')
                '''

            modCellN = []
            for zd in dnrZones:
                zdname = zd[0]
                destProc = procDict[zdname]

                IDs = []
                for i in zd[2]:
                    modif=0
                    if i[0][0:2] == 'ID':
                        if Internal.getValue(i)==zrname:
                            IDs.append(i)
                            modCellN.append(i)

                if IDs != []:
                    if destProc == Cmpi.rank:
                        zD = Internal.getNodeFromName2(tc, zdname)
                        #zD[2] += IDs
                    else:
                        if destProc not in datas: datas[destProc] = [[zdname,IDs]]
                        else: datas[destProc].append([zdname,IDs])
                else:
                    if destProc not in datas: datas[destProc] = []

            cellNModif[zrname]=modCellN  # a terminer: attention surememnt indice "bandelette" en mpi a cette position

    Cmpi._rmXZones(tc)
    #for h in hooks: C.freeHook(hooks[h])
    #test.printMem(">>> Interpdata [after rmXZones]")
    destDatas = Cmpi.sendRecv(datas, graph)
    for i in destDatas:
        for n in destDatas[i]:
            zname = n[0]
            IDs = n[1]
            if IDs != []:
                zD = Internal.getNodeFromName2(tc, zname)
                zD[2] += IDs


    datas = {}; destDatas = None;# graph={}

    #C.convertPyTree2File(t,'tE1.cgns')
    #C.convertPyTree2File(tc,'tcE1.cgns')

    #
    #etape 2 : modif CellN=0 pour les points interpoles a l'etape 1  et on vire info orphelin
    #
    #
    # Ne fonctionne pas en mpi: faudrait tranferer le ptlist receveur pour modifer cellN
    #
    #
    val=0  #on force cellN a zero car etape 1 et 3 = nature 1
    #if nature==0: val=1 #Ghost donneuse acceptee
    for zd in Internal.getZones(tc):
        #subRegions =  Internal.getNodesFromType1(zd, 'ZoneSubRegion_t')
        subRegions =  cellNModif[ zd[0] ]
        for s in subRegions:
            zRname = Internal.getValue(s)
            #faire test sur dest proc

            zr = Internal.getNodeFromName2(t,zRname)
            sol= Internal.getNodeFromName2(zr, 'FlowSolution#Centers')
            cellN= Internal.getNodeFromName1(sol, 'cellN')[1]


            sh   = numpy.shape(cellN)
            nxny = sh[0]*sh[1]
            nx   = sh[0]
            pointlistD =  Internal.getNodeFromName1(s, 'PointListDonor')[1]
            if len(sh)==2:
                for l in range(numpy.size(pointlistD)):
                    j  = pointlistD[l]//nx
                    i  = pointlistD[l] -j*nx
                    cellN[i,j]=val #si nature=0 , 0 sinon
            else:
                for l in range(numpy.size(pointlistD)):
                    k    = pointlistD[l]//nxny
                    rest = pointlistD[l]-k*nxny
                    j    = rest//nx
                    i    = rest -j*nx

                    cellN[i,j,k]=val #si nature=0 , 0 sinon

    '''
    '''
    fix='Nocorner'
    if corner: fix='corner'
    #zr = Internal.getNodeFromName2(t,'Cart.553X0')
    #if zr[0]=='Cart.553X0':
    #    C.convertPyTree2File(zr,'t553_E2_'+fix+'.cgns')

    C._cpVars(t,'centers:cellN',tc,'cellN')
    Cmpi._addXZones(tc, graph, variables=['cellN'], noCoordinates=False, cartesian=cartesian, zoneGC=False, keepOldNodes=False)


    #
    #etape3: calcul interpolation si grille donneuse plus fine que receuveuse (ordre2, nature1, pas d'extrap)
    #
    cellNModif={}
    for level in range(Nlevels):
        #filtre les zones par niveau de resolution
        zones=[]
        for z in Internal.getZones(t):
            if levelZone[z[0]]==level:
                zones.append(z)
                #print("level=", level,"zone In", z[0])

        print("Interp level=",level, "etape 3")
        for zr in zones:
            zrname = zr[0]
            baseNameRcv = baseNames[zrname]

            dnrZones = []
            for zdname in interDict[zrname]:
                zd = Internal.getNodeFromName2(tc, zdname)
                baseNameDnr = baseNames[zd[0]]
                if levelZone[zd[0]] <  level and (sameBase ==1 or baseNameDnr == baseNameRcv): dnrZones.append(zd)

            hookL = []; interpDataTypeL = []
            for z in dnrZones:
                h = hooks[z[0]]
                hookL.append(h)
                if h is None: interpDataTypeL.append(0)
                else: interpDataTypeL.append(1)

            if dnrZones != []:
                X._setInterpData(zr, dnrZones, nature=1, penalty=1, extrap=0, loc='centers', storage='inverse',
                                 interpDataType=interpDataTypeL, hook=hookL, sameName=1, order=2, itype='chimera', verbose=0)

                fix='Nocorner'
                if corner: fix='corner'
                '''
                if zr[0]=='Cart.553X0':
                    C.convertPyTree2File(zr,'t553_E3_'+fix+'.cgns')
                    C.convertPyTree2File(dnrZones,'tc553_E3_'+fix+'.cgns')
                '''

            modCellN = []
            for zd in dnrZones:
                zdname = zd[0]
                destProc = procDict[zdname]

                IDs = []
                for i in zd[2]:
                    modif=0
                    if i[0][0:2] == 'ID':
                        if Internal.getValue(i)==zrname:
                            IDs.append(i)
                            modCellN.append(i)

                if IDs != []:
                    if destProc == Cmpi.rank:
                        zD = Internal.getNodeFromName2(tc, zdname)
                        #zD[2] += IDs
                    else:
                        if destProc not in datas: datas[destProc] = [[zdname,IDs]]
                        else: datas[destProc].append([zdname,IDs])
                else:
                    if destProc not in datas: datas[destProc] = []

            cellNModif[zrname]=modCellN  # a terminer: attention surememnt indice "bandelette" en mpi a cette position

    Cmpi._rmXZones(tc)
    #for h in hooks: C.freeHook(hooks[h])
    #test.printMem(">>> Interpdata [after rmXZones]")
    destDatas = Cmpi.sendRecv(datas, graph)
    for i in destDatas:
        for n in destDatas[i]:
            zname = n[0]
            IDs = n[1]
            if IDs != []:
                zD = Internal.getNodeFromName2(tc, zname)
                zD[2] += IDs


    datas = {}; destDatas = None;# graph={}

    #
    #etape 4 : modif CellN=0 pour les points interpoles a l'etape 1  et on vire info orphelin
    #
    #
    # Ne fonctionne pas en mpi: faudrait tranferer le ptlist receveur pour modifer cellN
    #
    #
    val=0
    #if nature==0: val=1 #modif suite etape 6 et 7
    for zd in Internal.getZones(tc):
        subRegions =  Internal.getNodesFromType1(zd, 'ZoneSubRegion_t')
        #subRegions =  cellNModif[ zd[0] ]
        for s in subRegions:
            zRname = Internal.getValue(s)
            #faire test sur dest proc

            zr = Internal.getNodeFromName2(t,zRname)
            sol= Internal.getNodeFromName2(zr, 'FlowSolution#Centers')
            cellN= Internal.getNodeFromName1(sol, 'cellN')[1]

            sh   = numpy.shape(cellN)
            nxny = sh[0]*sh[1]
            nx   = sh[0]
            pointlistD =  Internal.getNodeFromName1(s, 'PointListDonor')[1]
            if len(sh)==2:
                for l in range(numpy.size(pointlistD)):
                    j  = pointlistD[l]//nx
                    i  = pointlistD[l] -j*nx
                    cellN[i,j]=val #si nature=0 , 0 sinon
            else:
                for l in range(numpy.size(pointlistD)):
                    k    = pointlistD[l]//nxny
                    rest = pointlistD[l]-k*nxny
                    j    = rest//nx
                    i    = rest -j*nx

                    cellN[i,j,k]=val #si nature=0 , 0 sinon

    fix='Nocorner'
    if corner: fix='corner'
    #zr = Internal.getNodeFromName2(t,'Cart.553X0')
    #if zr[0]=='Cart.553X0':
    #    C.convertPyTree2File(zr,'t553_E4_'+fix+'.cgns')

    C._cpVars(t,'centers:cellN',tc,'cellN')

    Cmpi._addXZones(tc, graph, variables=['cellN'], noCoordinates=False, cartesian=cartesian, zoneGC=False, keepOldNodes=False)


    Internal._rmNodesByName(tc,'*Orphan*')
    Internal._rmNodesByName(t,'*Orphan*')
    #
    #etape 5: calcul interpolation entre grille de niveau N (Receveur) et N+1 (donneurs)
    #
    for level in range(Nlevels):
        print("Interp level=",level, "etape 5")
        #filtre les zones par niveau de resolution
        zones=[]
        for z in Internal.getZones(t):
            if levelZone[z[0]]==level: zones.append(z)

        for zr in zones:
            dnrZones = []
            zrname = zr[0]
            baseNameRcv = baseNames[zrname]
            for zdname in interDict[zrname]:
                zd = Internal.getNodeFromName2(tc, zdname)
                baseNameDnr = baseNames[zd[0]]
                if levelZone[zd[0]]==level+1 and (sameBase ==1 or baseNameDnr == baseNameRcv):
                    dnrZones.append(zd)
                #print("level=", level,"zone In", z[0])

            hookL = []; interpDataTypeL = []
            for z in dnrZones:
                h = hooks[z[0]]
                hookL.append(h)
                if h is None: interpDataTypeL.append(0)
                else: interpDataTypeL.append(1)

            if dnrZones != []:
                #X._setInterpData(zr, dnrZones, nature=nature, penalty=1, loc='centers', storage='inverse', extrap=extrap, verbose=3,
                #                 sameName=1, interpDataType=interpDataTypeL, hook=hookL, order=order, itype='chimera')
                X._setInterpData(zr, dnrZones, nature=1, penalty=1, loc='centers', storage='inverse', extrap=0, verbose=0,
                                 sameName=1, interpDataType=interpDataTypeL, hook=hookL, order=order, itype='chimera')
                fix='Nocorner'
                if corner: fix='corner'
                '''
                if zr[0]=='Cart.553X0':
                      C.convertPyTree2File(zr,'t553_E5_'+fix+'.cgns')
                      C.convertPyTree2File(dnrZones,'tc553_E5_'+fix+'.cgns')
                '''
            for zd in dnrZones:
                zdname = zd[0]
                destProc = procDict[zdname]
                IDs = []
                for i in zd[2]:
                    if i[0][0:2] == 'ID':
                        if Internal.getValue(i)==zrname:
                            IDs.append(i)
                if IDs != []:
                    if destProc == Cmpi.rank:
                        zD = Internal.getNodeFromName2(tc, zdname)
                        #zD[2] += IDs
                    else:
                        if destProc not in datas: datas[destProc] = [[zdname,IDs]]
                        else: datas[destProc].append([zdname,IDs])
                else:
                    if destProc not in datas: datas[destProc] = []

    Cmpi._rmXZones(tc)
    for h in hooks: C.freeHook(hooks[h])
    #test.printMem(">>> Interpdata [after rmXZones]")
    destDatas = Cmpi.sendRecv(datas, graph)
    for i in destDatas:
        for n in destDatas[i]:
            zname = n[0]
            IDs = n[1]
            if IDs != []:
                zD = Internal.getNodeFromName2(tc, zname)
                zD[2] += IDs
    datas = {}; destDatas = None; graph={}

    if nature==0:
        #
        #etape 6 : modif CellN=0 pour les points interpoles a l'etape 5  et on vire info orphelin
        #
        #
        # Ne fonctionne pas en mpi: faudrait tranferer le ptlist receveur pour modifer cellN
        #
        #
        val=1 #Ghost donneuse acceptee, mais pas les coin car cellN=0 si etape1 et not corner
        for zd in Internal.getZones(tc):
            subRegions =  Internal.getNodesFromType1(zd, 'ZoneSubRegion_t')
            #subRegions =  cellNModif[ zd[0] ]
            for s in subRegions:
                zRname = Internal.getValue(s)
                #faire test sur dest proc

                zr = Internal.getNodeFromName2(t,zRname)
                sol= Internal.getNodeFromName2(zr, 'FlowSolution#Centers')
                cellN= Internal.getNodeFromName1(sol, 'cellN')[1]

                sh   = numpy.shape(cellN)
                nxny = sh[0]*sh[1]
                nx   = sh[0]
                pointlistD =  Internal.getNodeFromName1(s, 'PointListDonor')[1]
                if len(sh)==2:
                    for l in range(numpy.size(pointlistD)):
                        j  = pointlistD[l]//nx
                        i  = pointlistD[l] -j*nx
                        cellN[i,j]=val #si nature=0 , 0 sinon
                else:
                    for l in range(numpy.size(pointlistD)):
                        k    = pointlistD[l]//nxny
                        rest = pointlistD[l]-k*nxny
                        j    = rest//nx
                        i    = rest -j*nx

                        cellN[i,j,k]=val #si nature=0 , 0 sinon

        fix='Nocorner'
        if corner: fix='corner'
        #zr = Internal.getNodeFromName2(t,'Cart.553X0')
        #if zr[0]=='Cart.553X0':
        #   C.convertPyTree2File(zr,'t553_E6_'+fix+'.cgns')

        C._cpVars(t,'centers:cellN',tc,'cellN')

        Cmpi._addXZones(tc, graph, variables=['cellN'], noCoordinates=False, cartesian=cartesian, zoneGC=False, keepOldNodes=False)

        Internal._rmNodesByName(tc,'*Orphan*')
        Internal._rmNodesByName(t,'*Orphan*')

        #
        #etape 7: calcul interpolation entre grille de niveau N (Receveur) et N+1 (donneurs)
        #
        for level in range(Nlevels):
            print("Interp level=",level, "etape 7")
            #filtre les zones par niveau de resolution
            zones=[]
            for z in Internal.getZones(t):
                if levelZone[z[0]]==level: zones.append(z)

            for zr in zones:
                dnrZones = []
                zrname = zr[0]
                baseNameRcv = baseNames[zrname]
                for zdname in interDict[zrname]:
                    zd = Internal.getNodeFromName2(tc, zdname)
                    baseNameDnr = baseNames[zd[0]]
                    if levelZone[zd[0]]==level+1 and (sameBase ==1 or baseNameDnr == baseNameRcv):
                        dnrZones.append(zd)
                    #print("level=", level,"zone In", z[0])

                hookL = []; interpDataTypeL = []
                for z in dnrZones:
                    h = hooks[z[0]]
                    hookL.append(h)
                    if h is None: interpDataTypeL.append(0)
                    else: interpDataTypeL.append(1)

                if dnrZones != []:
                    X._setInterpData(zr, dnrZones, nature=nature, penalty=1, loc='centers', storage='inverse', extrap=extrap, verbose=3,
                                     sameName=1, interpDataType=interpDataTypeL, hook=hookL, order=order, itype='chimera')

                    fix='Nocorner'
                    if corner: fix='corner'
                    '''
                    if zr[0]=='Cart.553X0':
                       C.convertPyTree2File(zr,'t553_E7_'+fix+'.cgns')
                       C.convertPyTree2File(dnrZones,'tc553_E7_'+fix+'.cgns')
                    '''

                for zd in dnrZones:
                    zdname = zd[0]
                    destProc = procDict[zdname]
                    IDs = []
                    for i in zd[2]:
                        if i[0][0:2] == 'ID':
                            if Internal.getValue(i)==zrname: IDs.append(i)
                    if IDs != []:
                        if destProc == Cmpi.rank:
                            zD = Internal.getNodeFromName2(tc, zdname)
                            #zD[2] += IDs
                        else:
                            if destProc not in datas: datas[destProc] = [[zdname,IDs]]
                            else: datas[destProc].append([zdname,IDs])
                    else:
                        if destProc not in datas: datas[destProc] = []

        Cmpi._rmXZones(tc)
        for h in hooks: C.freeHook(hooks[h])
        #test.printMem(">>> Interpdata [after rmXZones]")
        destDatas = Cmpi.sendRecv(datas, graph)
        for i in destDatas:
            for n in destDatas[i]:
                zname = n[0]
                IDs = n[1]
                if IDs != []:
                    zD = Internal.getNodeFromName2(tc, zname)
                    zD[2] += IDs
        datas = {}; destDatas = None; graph={}

    #Recuperation du vrai cellN
    if corner:
        C._initVars(t, '{centers:cellN}= {centers:cellNSave1}')
    else:
        C._initVars(t, '{centers:cellN}= {centers:cellNSave2}')
        C._cpVars(t,'centers:cellN',tc,'cellN')

        Internal._rmNodesByName(t,'cellNSave1')
        Internal._rmNodesByName(t,'cellNSave2')
        Internal._rmNodesByName(t,'cellNSave3')

        Internal._rmNodesByName(tc,'cellNSave2')
        '''
        '''

    Internal._rmNodesByName(tc,'*Orphan*')

    return None

##
## modif cellN coin pour interpoler uniquemement les cellules "coin" donneuse en nature=0
##
def _correctCellNCorner(t, tc, dim=3, verbose=0):

    dimZones={}
    dicZones={}
    #on flag les coins a 1, saud 1ere rangee
    C._initVars(tc, '{cellNCornerD}= 0')
    C._initVars(tc, '{cellNCornerR}= 0')
    for z in Internal.getZones(tc):
        dimR          = Internal.getZoneDim(z)
        dicZones[z[0]]= z
        dimZones[z[0]]= dimR
        cellN = Internal.getNodeFromName(z, "cellNCornerR")[1]
        sh_R    = numpy.shape(cellN)

        val=0  #on force cellN a zero car etape 1 et 3 = nature 1
        if dim==2:
            #ghost 2eme couronne a 1
            ni =sh_R[0]; nj =sh_R[1]
            for i in range(2):
                cellN[i     , 0      ]=1
                cellN[ni-1-i, 0      ]=1
                cellN[i     , nj-1   ]=1
                cellN[ni-1-i, nj-1   ]=1
            cellN[0     , 1      ]=1
            cellN[ni-1  , 1      ]=1
            cellN[0     , nj-2   ]=1
            cellN[ni-1  , nj-2   ]=1

            #ghost classique et cell reelle a 3
            cellN[ :     , 2:nj-2]=3
            cellN[2:ni-2, 0:2    ]=3
            cellN[2:ni-2, nj-2:nj]=3

        else:
            ni =sh_R[0]; nj =sh_R[1];nk =sh_R[2]
            for i in range(2):
                cellN[i     , 0      , :]=1
                cellN[ni-1-i, 0      , :]=1
                cellN[i     , nj-1   , :]=1
                cellN[ni-1-i, nj-1   , :]=1
            cellN[0     , 1      , :]=1
            cellN[ni-1  , 1      , :]=1
            cellN[0     , nj-2   , :]=1
            cellN[ni-1  , nj-2   , :]=1

            for j in range(2):
                cellN[2:ni-2, j      , 0     ]=1
                cellN[2:ni-2, nj-1-j , 0     ]=1
                cellN[2:ni-2, j      , nk-1  ]=1
                cellN[2:ni-2, nj-1-j , nk-1  ]=1
            cellN[2:ni-2, 0      ,    1  ]=1
            cellN[2:ni-2, nj-1   ,    1  ]=1
            cellN[2:ni-2, 0      ,  nk-2 ]=1
            cellN[2:ni-2, nj-1   ,  nk-2 ]=1

            for i in range(2):
                cellN[i     ,2:nj-2, 0   ]=1
                cellN[ni-1-i,2:nj-2, 0   ]=1
                cellN[i     ,2:nj-2, nk-1]=1
                cellN[ni-1-i,2:nj-2, nk-1]=1
            cellN[0     , 2:nj-2 , 1   ]=1
            cellN[ni-1  , 2:nj-2 , 1   ]=1
            cellN[0     , 2:nj-2 , nk-2]=1
            cellN[ni-1  , 2:nj-2 , nk-2]=1

            #flag point reel a 3
            cellN[:     , 2:nj-2  , 2:nk-2 ]=3
            cellN[2:ni-2, 0:2     , 2:nk-2 ]=3
            cellN[2:ni-2, nj-2:nj , 2:nk-2 ]=3
            cellN[2:ni-2, 2:nj-2 ,  0:2    ]=3
            cellN[2:ni-2, 2:nj-2 ,  nk-2:nk]=3


    for z in Internal.getZones(tc):
        #zd = Internal.getNodeFromName(t, z[0])
        cellNCornerD_D = Internal.getNodeFromName(z, "cellNCornerD")[1]
        cellNCornerR_D = Internal.getNodeFromName(z, "cellNCornerR")[1]
        subRegions  =  Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        for s in subRegions:
            count =0
            zRname = Internal.getValue(s)

            dimD = Internal.getZoneDim(z)
            dimR = dimZones[zRname]

            zr = dicZones[zRname]
            flowSol        = Internal.getNodeFromName1(zr, "FlowSolution")
            cellNCornerR_R = Internal.getNodeFromName1(flowSol, "cellNCornerR")[1]

            pointlist  =  Internal.getNodeFromName1(s, 'PointList')[1]
            pointlistD =  Internal.getNodeFromName1(s, 'PointListDonor')[1]
            Interptype =  Internal.getNodeFromName1(s, 'InterpolantsType')[1]
            coeff      =  Internal.getNodeFromName1(s, 'InterpolantsDonor')[1]

            nxnyD = dimD[1]*dimD[2]
            nxD   = dimD[1]
            nxnyR = dimR[1]*dimR[2]
            nxR   = dimR[1]

            #Flag donor
            for l in range(numpy.size(pointlist)):

                kD    = pointlist[l]//nxnyD
                rest = pointlist[l]-kD*nxnyD
                jD    = rest//nxD
                iD    = rest -jD*nxD

                kR    = pointlistD[l]//nxnyR
                rest = pointlistD[l]-kR*nxnyR
                jR    = rest//nxR
                iR    = rest -jR*nxR

                #if zRname=='Cart.2X0': print('ijkD', iD, jD,kD, 'zD:', z[0], 'R_R', cellNCornerR_R[iR,jR,kR], 'ijkR', iR,jR,kR, 'type', Interptype[l],'R_D', cellNCornerR_D[iD,jD,kD] )

                if dim==3:
                    if cellNCornerR_R[iR,jR,kR]==3:

                        if Interptype[l]==1 and cellNCornerR_D[iD,jD,kD]<=1:
                            cellNCornerD_D[iD,jD,kD]=1
                            count +=1

                        elif Interptype[l]==2:
                            count_loc=0
                            for kk in range(2):
                                for jj in range(2):
                                    for ii in range(2):
                                        if coeff[count +count_loc] > 1.e-8 and cellNCornerR_D[iD +ii ,jD +jj, kD+ kk]<=1: cellNCornerD_D[iD +ii ,jD +jj, kD+ kk]=1
                                        count_loc+=1
                            count+=8

                        elif Interptype[l]==44:

                            for kk in range(4):
                                for jj in range(4):
                                    for ii in range(4):
                                        val = coeff[count +ii ] * coeff[count +jj + 4 ] * coeff[count +kk +8 ]
                                        #if zd[0]=='Cart.79X0': print('type 44: ijkD', iD, jD,kD, 'zR:', zRname, val, 'iijjkk', ii,jj,kk )
                                        if abs(val) > 1.e-11 and cellNCornerR_D[iD +ii, jD +jj, kD +kk]<=1: cellNCornerD_D[iD +ii, jD +jj, kD +kk ]=1
                            count+=12

                        elif Interptype[l]==22:

                            count_loc=0
                            for jj in range(2):
                                for ii in range(2):
                                    if coeff[count +count_loc] > 1.e-8 and cellNCornerR_D[iD +ii ,jD +jj, kD]<=1: cellNCornerD_D[iD +ii ,jD +jj, kD]=1
                                    count_loc+=1
                            count+=4
                else: #dim=2
                    if cellNCornerR_R[iR,jR]==3:

                        if Interptype[l]==1 and cellNCornerR_D[iD,jD]<=1:
                            cellNCornerD_D[iD,jDD]=1
                            count +=1

                        elif Interptype[l]==22:

                            count_loc=0
                            for jj in range(2):
                                for ii in range(2):
                                    if coeff[count +count_loc] > 1.e-8 and cellNCornerR_D[iD +ii ,jD +jj]<=1: cellNCornerD_D[iD +ii ,jD +jj]=1
                                    count_loc+=1
                            count+=4

    Internal._rmNodesByName(tc,'ID_*')

    for z in Internal.getZones(t):
        cellN       = Internal.getNodeFromName(z, "cellN")[1]
        zc          = dicZones[z[0]]
        cellNCorner = Internal.getNodeFromName(zc, "cellNCornerD")[1]
        sh_R    = numpy.shape(cellN)

        val=0  #on force cellN a zero car etape 1 et 3 = nature 1
        if dim==2:
            ni =sh_R[0]; nj =sh_R[1]
            for j in range(2):
                for i in range(2):
                    if cellNCorner[i      , j     ] < 0.1: cellN[i      , j     ]= min(val, cellN[i      , j ] )
                    if cellNCorner[ni-2+i , j     ] < 0.1: cellN[ni-2+i , j     ]= min(val, cellN[ni-2+i , j ] )
                    if cellNCorner[ni-2+i , nj-2+j] < 0.1: cellN[ni-2+i ,nj-2+j ]= min(val, cellN[ni-2+i , nj-2+j ] )
                    if cellNCorner[i      , nj-2+j] < 0.1: cellN[ i     ,nj-2+j ]= min(val, cellN[ i     , nj-2+j ] )
        else:
            ni =sh_R[0]; nj =sh_R[1];nk =sh_R[2];
            for k in range(nk):
                for j in range(2):
                    for i in range(2):
                        if cellNCorner[i     , j     , k] < 0.1: cellN[i     , j      , k]=min(val, cellN[i     , j     , k] )
                        if cellNCorner[i     , nj-2+j, k] < 0.1: cellN[i     , nj-2+j , k]=min(val, cellN[i     , nj-2+j, k] )
                        if cellNCorner[ni-2+i, j     , k] < 0.1: cellN[ni-2+i, j      , k]=min(val, cellN[ni-2+i, j     , k] )
                        if cellNCorner[ni-2+i, nj-2+j, k] < 0.1: cellN[ni-2+i, nj-2+j , k]=min(val, cellN[ni-2+i, nj-2+j, k] )
            for k in range(2):
                for j in range(2):
                    for i in range(2,ni-2):
                        if cellNCorner[i     , j      , k     ] < 0.1: cellN[i     , j      , k     ]= min(val, cellN[i     , j     , k] )
                        if cellNCorner[i     , j      , nk-2+k] < 0.1: cellN[i     , j      , nk-2+k]= min(val, cellN[i     , j     , nk-2+k] )
                        if cellNCorner[i     , nj-2+j , k     ] < 0.1: cellN[i     , nj-2+j , k     ]= min(val, cellN[i     , nj-2+j, k     ] )
                        if cellNCorner[i     , nj-2+j , nk-2+k] < 0.1: cellN[i     , nj-2+j , nk-2+k]= min(val, cellN[i     , nj-2+j, nk-2+k] )
            for k in range(2):
                for j in range(2,nj-2):
                    for i in range(2):
                        if cellNCorner[i     , j      , k     ] < 0.1: cellN[i     , j      , k     ]= min(val, cellN[i     , j     , k     ] )
                        if cellNCorner[i     , j      , nk-2+k] < 0.1: cellN[i     , j      , nk-2+k]= min(val, cellN[i     , j     , nk-2+k] )
                        if cellNCorner[ni-2+i, j      , k     ] < 0.1: cellN[ni-2+i, j      , k     ]= min(val, cellN[ni-2+i, j     , k     ] )
                        if cellNCorner[ni-2+i, j      , nk-2+k] < 0.1: cellN[ni-2+i, j      , nk-2+k]= min(val, cellN[ni-2+i, j     , nk-2+k] )

    Internal._rmNodesByName(tc,'cellNCorner*')
