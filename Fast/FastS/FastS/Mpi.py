# FastS + MPI
from . import PyTree
from . import fasts

from .PyTree import display_temporal_criteria, _displayTemporalCriteria, displayTemporalCriteria, _createConvergenceHistory, createConvergenceHistory, extractConvergenceHistory, _extractConvergenceHistory, createStressNodes, createStatNodes, _computeStats, initStats, _computeEnstrophy, _computeVariables, _computeGrad, _applyBC, _init_metric, allocate_metric, _movegrid, _computeVelocityAle, copy_velocity_ale,  checkBalance, itt, allocate_ssor, setIBCData_zero, display_cpu_efficiency, _postStats, _stretch, _computePhaseStats, _ConservativeWallIbm, create_add_t_converg_hist, calc_global_convergence

import timeit
import time as Time
import numpy
import sys

import FastC.fastc

try:
    import Converter.PyTree as C
    import Converter.Mpi    as Cmpi
    import Distributor2.PyTree as D2
    import Converter.Internal as Internal
    import Connector.Mpi as Xmpi
    import Connector.PyTree as X
    import Connector
    import FastC.PyTree as FastC
    import RigidMotion.PyTree as R
    import os
    import math
except:
    raise ImportError("FastS: requires Converter, Connector and Distributor2 modules.")

try: range = xrange
except: pass

#==============================================================================
def _compute(t, metrics, nitrun, tc=None, graph=None, tc2=None, graph2=None, layer="c", NIT=1, ucData=None, vtune=False):
    gradP      =False
    TBLE       =False
    isWireModel=False
    if tc is not None:
        base       = Internal.getBases(tc)[0]
        solverIBC  = Internal.getNodeFromName(base ,'.Solver#IBCdefine')
        if solverIBC is not None:
            gradP      = eval(Internal.getValue(Internal.getNodeFromName(solverIBC, 'isgradP')))
            TBLE       = eval(Internal.getValue(Internal.getNodeFromName(solverIBC, 'isTBLE')))
            isWireModel= eval(Internal.getValue(Internal.getNodeFromName(solverIBC, 'isWireModel')))

    if isWireModel and tc2 is not None:
        print("Wire model does not currently support tc2 options...exiting...")
        exit()
	
    if isinstance(graph, list):
        #test pour savoir si graph est une liste de dictionnaires (explicite local)
        #ou juste un dictionnaire (explicite global, implicite)
        grapheliste=True
    else:
        grapheliste=False
    
    if graph is not None and not grapheliste:
        procDict  = graph['procDict']
        graphID   = graph['graphID']
        graphIBCD = graph['graphIBCD']
        graphInvIBCD_WM  = None
        if tc2 is not None:
          graphIBCD2 = graph2['graphIBCD']
          graphInvIBCD = Cmpi.computeGraph(tc, type='INV_IBCD', procDict=procDict)
          graphInvIBCD_WM  = Cmpi.computeGraph(tc, type='INV_IBCD', procDict=procDict)
    else: 
        procDict=None; graphID=None; graphIBCD=None; graphInvIBCD_WM=None

    own  = Internal.getNodeFromName1(t, '.Solver#ownData')  
    dtloc= Internal.getNodeFromName1(own , '.Solver#dtloc')

    zones= Internal.getZones(t)

    node = Internal.getNodeFromName1(t, '.Solver#define')
    omp_node = Internal.getNodeFromName1(node, 'omp_mode')
    ompmode  = FastC.OMP_MODE
    if  omp_node is not None: ompmode = Internal.getValue(omp_node)

    dtloc   = Internal.getValue(dtloc) # tab numpy
    nitmax  = int(dtloc[0])

    #### a blinder...
    param_int_firstZone = Internal.getNodeFromName2( zones[0], 'Parameter_int' )[1]
    itypcp = param_int_firstZone[29]
    rk     = param_int_firstZone[52]
    exploc = param_int_firstZone[54]
    #### a blinder...
    if nitrun == 1: print('Info: using layer trans=%s (ompmode=%d)'%(layer, ompmode))

    #determine layer pour transfert et sous-iteration
    if          "Python" == layer: layer_mode= 0
    elif "Python_ucdata" == layer: layer_mode= 2
    else                         : layer_mode= 1

    if layer_mode>=1: FastC.HOOK["mpi"] = 1

    timelevel_target = int(dtloc[4])
    tps_cp = Time.time(); tps_cp = tps_cp-tps_cp     
    tps_tr = Time.time(); tps_tr = tps_tr-tps_tr
    tps_tr1= Time.time(); tps_tr1= tps_tr1-tps_tr1
    if layer_mode==0 or layer_mode==2:

      if exploc==1 and tc is not None and layer_mode==0:
         tc_compact = Internal.getNodeFromName1(tc, 'Parameter_real')
         if tc_compact is not None:
                param_real_tc= tc_compact[1]
                param_int_tc = Internal.getNodeFromName1(tc, 'Parameter_int' )[1]

                zones_tc    = Internal.getZones(tc)
                nbcomIBC    = param_int_tc[1]
                shift_graph = nbcomIBC + param_int_tc[2+nbcomIBC] + 2
                comm_P2P    = param_int_tc[0]
                pt_ech      = param_int_tc[comm_P2P + shift_graph]
                dest        = param_int_tc[pt_ech]

      hookTransfer = []
      hook1        = FastC.HOOK.copy()
      for nstep in range(1, nitmax+1): # pas RK ou ssiterations
         skip      = 0
         if layer_mode==0:
           hook1.update(  FastC.fastc.souszones_list(zones, metrics, FastC.HOOK, nitrun, nstep, 0) )
           if hook1["lssiter_verif"] == 0 and nstep == nitmax and itypcp ==1: skip = 1

         # calcul Navier Stokes + appli CL
         if skip == 0:
            nstep_deb = nstep
            nstep_fin = nstep
            nit_c     = 1
            tic = Time.time()

            tps_trIt  = fasts._computePT(zones, metrics, nitrun, nstep_deb, nstep_fin, layer_mode, nit_c, hook1)
            tps_cp += Time.time()-tic - tps_trIt  
            tps_tr += tps_trIt  

            # dtloc GJeanmasson
            if exploc==1 and tc is not None and layer_mode==0:

               no_transfert = 1#comm_P2P
               process = Cmpi.rank

               zonesD    = Internal.getZones(tc)
               
               fasts.dtlocal2para_mpi(zones,zonesD,param_int_tc,param_real_tc,hook1,0,nstep,ompmode,no_transfert,process)

               if    nstep%2 == 0 and itypcp == 2: vars = ['Density'  ] 
               elif  nstep%2 == 1 and itypcp == 2: vars = ['Density_P1'] 

               _applyBC(zones,metrics, hook1, nstep, var=vars[0])    

               FastC.switchPointers2__(zones,nitmax,nstep)
               # Ghostcell

               if    nstep%2 == 0 and itypcp == 2: vars = ['Density'  ]  # Choix du tableau pour application transfer et BC
               elif  nstep%2 == 1 and itypcp == 2: vars = ['Density_P1']

               if graph is not None and grapheliste == True:
                       procDict  = graph[nstep-1]['procDict']
                       graphID   = graph[nstep-1]['graphID']
                       graphIBCD = graph[nstep-1]['graphIBCD']
               else: 
                       procDict=None; graphID=None; graphIBCD=None          
               _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, ompmode, hook1, graphID, graphIBCD, procDict, nitmax, rk, exploc, isWireModel=isWireModel)

               no_transfert = 1#comm_P2P      

               fasts.recup3para_mpi(zones,zones_tc, param_int_tc, param_real_tc, hook1, 0, nstep, ompmode, 1) 

               if nstep%2 == 0:
                   vars = ['Density']
                   if graph is not None and grapheliste == True:
                       procDict  = graph[nitmax+nstep-1]['procDict']
                       graphID   = graph[nitmax+nstep-1]['graphID']
                       graphIBCD = graph[nitmax+nstep-1]['graphIBCD']
                   else: 
                       procDict=None; graphID=None; graphIBCD=None   
                   _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, nitmax, hook1, graphID, graphIBCD, procDict, nitmax, rk, exploc, num_passage=2, isWireModel=isWireModel)

               if    nstep%2 == 0 and itypcp == 2: vars = ['Density'  ] 
               elif  nstep%2 == 1 and itypcp == 2: vars = ['Density_P1'] 
               _applyBC(zones,metrics, hook1, nstep, var=vars[0])


            else: ### Autres schemas
               # Ghostcell
               vars = FastC.varsP
               if nstep%2 == 0 and itypcp == 2: vars = FastC.varsN  # Choix du tableau pour application transfer et BC

               tic=Time.time()

               if not tc2:
                  if layer_mode==0:
                     _fillGhostcells(zones, tc, metrics, timelevel_target , vars, nstep, ompmode, hook1, graphID, graphIBCD, procDict, gradP=gradP, TBLE=TBLE, isWireModel=isWireModel, graphInvIBCD_WM=graphInvIBCD_WM)
               else:
                _fillGhostcells2(zones, tc, tc2, metrics, timelevel_target, vars, nstep, ompmode, hook1, graphID, graphIBCD, procDict, gradP=gradP, TBLE=TBLE, graphInvIBCD=graphInvIBCD, graphIBCD2=graphIBCD2)

               # add unsteady Chimera transfers (motion) here
               if ucData is not None:
                   (graphX, intersectionDict, dictOfADT, 
                   dictOfNobOfRcvZones, dictOfNozOfRcvZones,
                   dictOfNobOfDnrZones, dictOfNozOfDnrZones, 
                   dictOfNobOfRcvZonesC, dictOfNozOfRcvZonesC,
                   time, procDict, interpInDnrFrame, varType, tfreq, order, verbose) = ucData

                   if nstep%2 == 0 and itypcp == 2: 
                       VARS = ['Density', 'VelocityX', 'VelocityY', 'VelocityZ', 'Temperature']
                       if varType == 1: VARS += ['TurbulentSANuTilde']
                   else:
                       VARS = ['Density_P1', 'VelocityX_P1', 'VelocityY_P1', 'VelocityZ_P1', 'Temperature_P1']
                       if varType == 1: VARS += ['TurbulentSANuTilde_P1']
                   for v in VARS: C._cpVars(t, 'centers:'+v, tc, v)
                   C._cpVars(t, "centers:cellN", tc, "cellN")

                   if nstep == nitmax or nstep%tfreq == 0:
                     Xmpi._transfer2(t, tc, VARS, graphX, intersectionDict, dictOfADT, 
                                     dictOfNobOfRcvZones, dictOfNozOfRcvZones,
                                     dictOfNobOfDnrZones, dictOfNozOfDnrZones, 
                                     dictOfNobOfRcvZonesC, dictOfNozOfRcvZonesC, 
                                     time=time, absFrame=True,
                                     procDict=procDict, cellNName='cellN#Motion', 
                                     interpInDnrFrame=interpInDnrFrame, order=order, 
                                     hook=hookTransfer, verbose=verbose)

                     _applyBC(t,metrics, hook1, nstep, var=VARS[0])

               tps_tr1+= Time.time()-tic  

    else: 
      nstep_deb = 1
      nstep_fin = nitmax
      layer_mode= 1
      nit_c     = NIT
      tic = Time.time()
      tps_tr  = fasts._computePT(zones, metrics, nitrun, nstep_deb, nstep_fin, layer_mode, nit_c, FastC.HOOK)
      tps_cp += Time.time()-tic - tps_tr  
      tps_tr1 =0
    #switch pointer a la fin du pas de temps
    if exploc==1 and tc is not None:
         if layer_mode == 0: FastC.switchPointers__(zones, 1, 3)
         else: FastC.switchPointers3__(zones, nitmax)
    else:
         case = NIT%3
         if case != 0 and itypcp < 2: FastC.switchPointers__(zones, case)
         if case != 0 and itypcp ==2: FastC.switchPointers__(zones, case, nitmax%2+2)

    # flag pour derivee temporelle 1er pas de temps implicit
    FastC.HOOK["FIRST_IT"]  = 1
    FastC.FIRST_IT          = 1

    if vtune is None:
      return None
    else:
      return tps_cp,tps_tr,tps_tr1

#==============================================================================
def _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, omp_mode, hook1, graphID, graphIBCD, procDict, nitmax=1, rk=1, exploc=0, num_passage=1, gradP=False, TBLE=False, isWireModel=False, graphInvIBCD_WM=None):

   #timecount = numpy.zeros(4, dtype=numpy.float64)
   timecount = []

   if hook1['lexit_lu'] ==0:

       #transfert
       if tc is not None:
           tc_compact = Internal.getNodeFromName1(tc, 'Parameter_real')
           if tc_compact is not None:
              param_real= tc_compact[1]              
              param_int = Internal.getNodeFromName1( tc, 'Parameter_int' )[1]
              zonesD    = Internal.getZones(tc)

              if hook1['neq_max'] == 5: varType = 2
              else                    : varType = 21

              dtloc = hook1['dtloc']

              for v in vars: C._cpVars(zones, 'centers:'+v, zonesD, v)

              if gradP: 
                    varsGrad = ['Density', 'gradxDensity']
                    varType = 22
                    # varType 22 means that there is 6 variables and 6 gradVariables to transfers
                    # to avoid interfering with the original ___setInterpTransfers, we instead use ___setInterpTransfers4GradP
                    for v in varsGrad: C._cpVars(zones, 'centers:'+v, zonesD,  v)
                    type_transfert = 1
                    Xmpi.__setInterpTransfers4GradP(zones , zonesD, varsGrad, param_int, param_real, type_transfert, timelevel_target,
                                          nstep, nitmax, rk, exploc, num_passage, varType=varType, compact=1, graph=graphIBCD, procDict=procDict)
                    varType = 21

              #if Cmpi.rank == 0: print('fillGC: timeleveltarget= ', timelevel_target)
              if isWireModel:                  
                  ##The approach is the following
                  ##1) interpolate the values at the additional points for WM
                  ##2) copy these interpolated values into the tc
                  ##3) perfom the ibc & id interpolation as normal, excluding the additional points for WM

                  
                  ##1) Interpolate the values at the additional interpolated points for WM
                  C._cpVars(zones, 'centers:Density_WM', zonesD, 'Density_WM')
                  nvars_local     = hook1['neq_max']
                  type_transfert  = 1  

                  Xmpi.__setInterpTransfers(zones, zonesD, vars, dtloc, param_int, param_real, type_transfert, timelevel_target,#timecount,
                                            nstep, nitmax, rk, exploc, num_passage, varType=varType, compact=1, graph=graphIBCD, procDict=procDict,
                                            isWireModel_int=1)

                  ##2) Transfers info at target points from flow field to the tc 
                  Xmpi.__setInterpTransfers_WireModel(zones, zonesD, vars, dtloc, param_int, param_real, type_transfert, timelevel_target,
                                                      nstep, nitmax, rk, exploc, num_passage, varType=varType, compact=1, graph=graphIBCD, procDict=procDict,
                                                      graphIBCD=graphIBCD, graphInvIBCD_WM=graphInvIBCD_WM, nvars=nvars_local)
                  
                  C._rmVars(zonesD, 'Density_WM')
                  
                  ##3)
                  isWireModel_int = -1  # is a local flag for Xmpi.__setInterpTransfers
                  type_transfert  =  1  # 0= ID uniquememnt, 1= IBC uniquememnt, 2= All 
                  Xmpi.__setInterpTransfers(zones, zonesD, vars, dtloc, param_int, param_real, type_transfert, timelevel_target,#timecount,
                                            nstep, nitmax, rk, exploc, num_passage, varType=varType, compact=1, graph=graphID, procDict=procDict,
                                            isWireModel_int=isWireModel_int)
                  
                  isWireModel_int = 0
                  type_transfert  = 0  # 0= ID uniquememnt, 1= IBC uniquememnt, 2= All 
                  Xmpi.__setInterpTransfers(zones, zonesD, vars, dtloc, param_int, param_real, type_transfert, timelevel_target,#timecount,
                                            nstep, nitmax, rk, exploc, num_passage, varType=varType, compact=1, graph=graphID, procDict=procDict,
                                            isWireModel_int=isWireModel_int)
              else:
                  #recuperation Nb pas instationnaire dans tc
                  type_transfert = 1  # 0= ID uniquememnt, 1= IBC uniquememnt, 2= All 
                  Xmpi.__setInterpTransfers(zones, zonesD, vars, dtloc, param_int, param_real, type_transfert, timelevel_target,#timecount,
                                            nstep, nitmax, rk, exploc, num_passage, varType=varType, compact=1, graph=graphIBCD, procDict=procDict,
                                            isWireModel_int=int(isWireModel))
                  type_transfert = 0  # 0= ID uniquememnt, 1= IBC uniquememnt, 2= All 
                  Xmpi.__setInterpTransfers(zones, zonesD, vars, dtloc, param_int, param_real, type_transfert, timelevel_target,#timecount,
                                            nstep, nitmax, rk, exploc, num_passage, varType=varType, compact=1, graph=graphID, procDict=procDict,
                                            isWireModel_int=int(isWireModel))

              #toc = Time.time() - tic
       # if Cmpi.rank == 0:
       #     print("Time in MPI send_buffer, irecv ","%.6f"%timecount[0])
       #     print("Time InterpTransfert (Inter)  ","%.6f"%timecount[1])
       #     print("Time InterpTransfert (Intra)  ","%.6f"%timecount[2])
       #     print("Time in getTransfersInter ","%.6f"%timecount[3])
       # if Cmpi.rank == 0: t0=timeit.default_timer()
       #apply BC
       #tic = Time.time()
       if exploc != 1:
           _applyBC(zones, metrics, hook1, nstep, var=vars[0])
       #toc = Time.time() - tic
       # if Cmpi.rank == 0:
       #     t1=timeit.default_timer()
       #     print("time/it (s) BC only (=t_RK X 1.0): ",(t1-t0))
   #return toc
   return None

#==============================================================================
# modified _fillGhostCells for two image points (tc + tc2) for gradP 
#==============================================================================
def _fillGhostcells2(zones, tc, tc2, metrics, timelevel_target, vars, nstep, omp_mode, hook1, graphID, graphIBCD, procDict, nitmax=1, rk=1, exploc=0, num_passage=1, gradP=False, TBLE=False, graphInvIBCD=None, graphIBCD2=None): 
   
   #timecount = numpy.zeros(4, dtype=numpy.float64)
   timecount = []

   if hook1['lexit_lu'] ==0:

       #transfert
       if tc is not None:
           tc_compact = Internal.getNodeFromName1(tc, 'Parameter_real')
           if tc_compact is not None:
              param_real= tc_compact[1]              
              param_int = Internal.getNodeFromName1( tc, 'Parameter_int' )[1]
              zonesD    = Internal.getZones(tc)
              zonesD2   = Internal.getZones(tc2)

              if hook1['neq_max'] == 5: varType = 2
              else                    : varType = 21

              dtloc = hook1['dtloc']

              for v in vars: C._cpVars(zones, 'centers:'+v, zonesD, v)

              if gradP: 
                  tc2_compact = Internal.getNodeFromName1( tc2, 'Parameter_real')
                  param_real2= tc2_compact[1]
                  param_int2 = Internal.getNodeFromName1(tc2, 'Parameter_int' )[1]
                  varsGrad = ['Density', 'gradxDensity']
                  for v in varsGrad: 
                    C._cpVars(zones, 'centers:'+v, zonesD,  v)
                    C._cpVars(zones, 'centers:'+v, zonesD2,  v)
                  #tc2 -> RCV ZONES IBCD
                  type_transfert = 1; varType = 23
                  Xmpi.__setInterpTransfers4GradP(zones , zonesD2, varsGrad, param_int2, param_real2, type_transfert, timelevel_target,
                                        nstep, nitmax, rk, exploc, num_passage, varType=varType, compact=1, graph=graphIBCD2, procDict=procDict)
                  #RCV ZONES -> tc when RCV ZONES are on the same proc as IBCD ZONES
                  type_transfert = 1; varType = 24
                  Xmpi.__setInterpTransfers4GradP(zones , zonesD, varsGrad, param_int, param_real, type_transfert, timelevel_target,
                                       nstep, nitmax, rk, exploc, num_passage, varType=varType, compact=1, graph=graphIBCD, procDict=procDict)
                  #RCV ZONES -> tc when RCV ZONES are NOT on the same proc as IBCD ZONES
                  X._setIBCTransfers4GradP3(zones, zonesD, graphInvIBCD, graphIBCD, procDict, variablesIBC=varsGrad, alpha=Internal.getNodeFromName(zones, 'Parameter_real')[1][57])
                  varType = 21

              # #recuperation Nb pas instationnaire dans tc
              type_transfert = 1  # 0= ID uniquememnt, 1= IBC uniquememnt, 2= All 
              Xmpi.__setInterpTransfers(zones , zonesD, vars, dtloc, param_int, param_real, type_transfert, timelevel_target,#timecount,
                                        nstep, nitmax, rk, exploc, num_passage, varType=varType, compact=1,
                                        graph=graphIBCD, procDict=procDict, isWireModel=False)
              type_transfert = 0  # 0= ID uniquememnt, 1= IBC uniquememnt, 2= All 
              Xmpi.__setInterpTransfers(zones , zonesD, vars, dtloc, param_int, param_real, type_transfert, timelevel_target,#timecount,
                                        nstep, nitmax, rk, exploc, num_passage, varType=varType, compact=1, graph=graphID, procDict=procDict,
                                        isWireModel=False)

              #toc = Time.time() - tic
       # if Cmpi.rank == 0:
       #     print("Time in MPI send_buffer, irecv ","%.6f"%timecount[0])
       #     print("Time InterpTransfert (Inter)  ","%.6f"%timecount[1])
       #     print("Time InterpTransfert (Intra)  ","%.6f"%timecount[2])
       #     print("Time in getTransfersInter ","%.6f"%timecount[3])
       # if Cmpi.rank == 0: t0=timeit.default_timer()
       #apply BC
       #tic = Time.time()
       if exploc != 1:
           _applyBC(zones, metrics, hook1, nstep, var=vars[0])
       #toc = Time.time() - tic
       # if Cmpi.rank == 0:
       #     t1=timeit.default_timer()
       #     print("time/it (s) BC only (=t_RK X 1.0): ",(t1-t0))
   #return toc
   return None

#==============================================================================
# Compute linelets info
# Adaptive h0
#==============================================================================
def computeLineletsInfo(tc, Re=6.e6, Cf_law='ANSYS', Lref=1., q=1.2):
    from mpi4py import MPI

    if Cf_law == 'ANSYS':
        def compute_Cf(Re):
            return 0.058*Re**(-0.2)
    elif Cf_law == 'PW':
        def compute_Cf(Re):
            return 0.026*Re**(-1/7.)
    elif Cf_law == 'PipeDiameter':
        def compute_Cf(Re):
            return 0.079*Re**(-0.25)
    elif Cf_law == 'Laminar':
        def compute_Cf(Re):
            return 1.46*Re**(-0.5)

    hmod = 0.

    zones_tc = Internal.getZones(tc)

    for z in zones_tc:
        subRegions =  Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        if subRegions is not None:
              for s in subRegions:
                zsrname = s[0]
                sname  = zsrname[0:2]
                if sname == 'IB':
                  zI = Internal.getNodeFromName1(s,'CoordinateZ_PI')[1]
                  yI = Internal.getNodeFromName1(s,'CoordinateY_PI')[1]
                  xI = Internal.getNodeFromName1(s,'CoordinateX_PI')[1]
                  zW = Internal.getNodeFromName1(s,'CoordinateZ_PW')[1]
                  yW = Internal.getNodeFromName1(s,'CoordinateY_PW')[1]
                  xW = Internal.getNodeFromName1(s,'CoordinateX_PW')[1]
                  pointlistD    = Internal.getNodeFromName1(s, 'PointListDonor')
                  Nbpts_D       = numpy.shape(pointlistD[1])[0]

                  for i in range(Nbpts_D):
                    hmod = max(hmod, math.sqrt( (xI[i]-xW[i])**2 + (yI[i]-yW[i])**2 + (zI[i]-zW[i])**2 ))

    if Cmpi.KCOMM is not None:
      hmod = numpy.array([hmod])
      hmod_max = numpy.zeros(1)
      Cmpi.KCOMM.Allreduce(hmod, hmod_max, MPI.MAX)
      hmod = hmod_max[0]

    h0 = (Lref*math.sqrt(2))/(Re*math.sqrt(compute_Cf(Re)))
    # hmod = sum(0,n) h0*q**n
    n = int(math.ceil(math.log(1-((hmod/h0)*(1-q)))/(math.log(q)))) - 1
    # (n+1) intervalles / (n+2) points

    return (h0, h0*q**n, n+2)   


#==============================================================================
# Compute linelets info
# Fixed h0 1e-6
#==============================================================================
def computeLineletsInfo2(tc, q=1.2):
    from mpi4py import MPI
    hmod = 0.

    zones_tc = Internal.getZones(tc)

    for z in zones_tc:
        subRegions =  Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        if subRegions is not None:
              for s in subRegions:
                zsrname = s[0]
                sname  = zsrname[0:2]
                if sname == 'IB':
                  zI = Internal.getNodeFromName1(s,'CoordinateZ_PI')[1]
                  yI = Internal.getNodeFromName1(s,'CoordinateY_PI')[1]
                  xI = Internal.getNodeFromName1(s,'CoordinateX_PI')[1]
                  zW = Internal.getNodeFromName1(s,'CoordinateZ_PW')[1]
                  yW = Internal.getNodeFromName1(s,'CoordinateY_PW')[1]
                  xW = Internal.getNodeFromName1(s,'CoordinateX_PW')[1]
                  pointlistD    = Internal.getNodeFromName1(s, 'PointListDonor')
                  Nbpts_D       = numpy.shape(pointlistD[1])[0]

                  for i in range(Nbpts_D):
                    hmod = max(hmod, math.sqrt( (xI[i]-xW[i])**2 + (yI[i]-yW[i])**2 + (zI[i]-zW[i])**2 ))

    if Cmpi.KCOMM is not None:
      hmod = numpy.array([hmod])
      hmod_max = numpy.zeros(1)
      Cmpi.KCOMM.Allreduce(hmod, hmod_max, MPI.MAX)
      hmod = hmod_max[0]

    h0 = 1.e-6
    n = int(math.ceil(math.log(1-((hmod/h0)*(1-q)))/(math.log(q)))) - 1

    return (h0, h0*q**n, n+2)   


#==============================================================================
# Prepare IBC ODE (Spalart 1D)
# Info stored and compacted in IBCD ZONES
# Para friendly
#==============================================================================
def _createTBLESA2(t, tc, h0, hn, nbpts_linelets=45):
      nfields_lin  = 6 # u,nutild,y,matm,mat,matp
      count_racIBC = 0
      addr_IBC = 0
      addrIBC  = []
      size_IBC = 0
      isODEDataPresent = 0
      distIW = 0.

      zones_tc = Internal.getZones(tc)

      for z in zones_tc:
          subRegions =  Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
          if subRegions is not None:
                for s in subRegions:
                  zsrname = s[0]
                  sname  = zsrname[0:2]
                  if sname == 'IB':
                     zsrname = zsrname.split('_')
                     if len(zsrname)<3:
                        print('Warning: createTBLESA: non consistent with the version of IBM preprocessing.')
                     else:
                       if zsrname[1] == '16':
                        print('Info: using MuskerSA2 wall model')
                        pointlistD    = Internal.getNodeFromName1(s, 'PointListDonor')
                        Nbpts_D       = numpy.shape(pointlistD[1])[0]
                        addrIBC.append(size_IBC)
                        size_IBC      = size_IBC + Nbpts_D*(nbpts_linelets*nfields_lin + 1)
                        count_racIBC  = count_racIBC + 1
                        u_ODE = Internal.getNodeFromName1(s, 'VelocityT_ODE')
                        if u_ODE is not None: isODEDataPresent = 1
                       elif zsrname[1] == '17':
                        print('Info: using FULL_TBLE_SA wall model')
                        pointlistD    = Internal.getNodeFromName1(s, 'PointListDonor')
                        Nbpts_D       = numpy.shape(pointlistD[1])[0]
                        addrIBC.append(size_IBC)
                        size_IBC      = size_IBC + Nbpts_D*(nbpts_linelets*nfields_lin + 1)
                        count_racIBC  = count_racIBC + 1
                        u_ODE = Internal.getNodeFromName1(s, 'VelocityT_ODE')
                        if u_ODE is not None: isODEDataPresent = 1
                       elif zsrname[1] == '18':
                        print('Info: using FULL_TBLE_Prandtl wall model')
                        pointlistD    = Internal.getNodeFromName1(s, 'PointListDonor')
                        Nbpts_D       = numpy.shape(pointlistD[1])[0]
                        addrIBC.append(size_IBC)
                        size_IBC      = size_IBC + Nbpts_D*(nbpts_linelets*nfields_lin + 1)
                        count_racIBC  = count_racIBC + 1
                        u_ODE = Internal.getNodeFromName1(s, 'VelocityT_ODE')
                        if u_ODE is not None: isODEDataPresent = 1
                       elif zsrname[1] == '19':
                        print('Info: using MafzalSA wall model')
                        pointlistD    = Internal.getNodeFromName1(s, 'PointListDonor')
                        Nbpts_D       = numpy.shape(pointlistD[1])[0]
                        addrIBC.append(size_IBC)
                        size_IBC      = size_IBC + Nbpts_D*(nbpts_linelets*nfields_lin + 1)
                        count_racIBC  = count_racIBC + 1
                        u_ODE = Internal.getNodeFromName1(s, 'VelocityT_ODE')
                        if u_ODE is not None: isODEDataPresent = 1

      if size_IBC > 0:
        if not isODEDataPresent:   # No ODE Data in tc
         print('ODE data missing in tc, they will be built...')

         if Internal.getNodeFromName1(t, 'nbpts_linelets') is None: 
          Internal.createUniqueChild(t, 'nbpts_linelets', 'DataArray_t', value=nbpts_linelets)

         # Nbpts_Dtot = size_IBC/(nbpts_linelets*nfields_lin)
         Nbpts_Dtot = size_IBC//(nbpts_linelets*nfields_lin+1)
         _addrIBC = numpy.asarray([nbpts_linelets,count_racIBC-1,Nbpts_Dtot] + addrIBC + [0]*Nbpts_Dtot, dtype=Internal.E_NpyInt)

         linelets = numpy.zeros(size_IBC + Nbpts_Dtot, dtype=numpy.float64)

         shift = 0
         for z in zones_tc:
           subRegions =  Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
           if subRegions is not None:
              for s in subRegions:
                  zsrname = s[0]
                  sname  = zsrname[0:2]
                  if sname == 'IB':
                     zsrname = zsrname.split('_')
                     if len(zsrname)<3:
                       print('Warning: createTBLESA2: non consistent with the version of IBM preprocessing.')
                     else:
                       if zsrname[1] == '16' or zsrname[1] == '17' or zsrname[1] == '18' or zsrname[1] == '19':
                         pointlistD    = Internal.getNodeFromName1(s, 'PointListDonor')
                         Nbpts_D       = numpy.shape(pointlistD[1])[0]

                         gradxP = Internal.getNodeFromName1(s, 'gradxPressure')
                         if gradxP is None:
                            Internal.createUniqueChild(s, 'gradxPressure', 'DataArray_t',
                                                    numpy.zeros(Nbpts_D, numpy.float64))
                            Internal.createUniqueChild(s, 'gradyPressure', 'DataArray_t',
                                                    numpy.zeros(Nbpts_D, numpy.float64))
                            Internal.createUniqueChild(s, 'gradzPressure', 'DataArray_t',
                                                    numpy.zeros(Nbpts_D, numpy.float64))


                         zI = Internal.getNodeFromName1(s,'CoordinateZ_PI')
                         valzI = Internal.getValue(zI)
                         yI = Internal.getNodeFromName1(s,'CoordinateY_PI')
                         valyI = Internal.getValue(yI)
                         xI = Internal.getNodeFromName1(s,'CoordinateX_PI')
                         valxI = Internal.getValue(xI)
                         zW = Internal.getNodeFromName1(s,'CoordinateZ_PW')
                         valzW = Internal.getValue(zW)
                         yW = Internal.getNodeFromName1(s,'CoordinateY_PW')
                         valyW = Internal.getValue(yW)
                         xW = Internal.getNodeFromName1(s,'CoordinateX_PW')
                         valxW = Internal.getValue(xW)

                         if Nbpts_D > 1:

                             for noind in range(0, Nbpts_D):

                                 b0 = valxI[noind]-valxW[noind]
                                 b1 = valyI[noind]-valyW[noind]
                                 b2 = valzI[noind]-valzW[noind]
                                 normb = numpy.sqrt(b0*b0+b1*b1+b2*b2)

                                 if hn < 0:
                                   hn = 0.6*normb

                                 # _stretch(linelets[shift + noind*nbpts_linelets: shift + (noind+1)*nbpts_linelets],nbpts_linelets,0.0,normb, 1.0e-6, 0.6*normb, 2)
                                 _stretch(linelets[shift + noind*nbpts_linelets: shift + (noind+1)*nbpts_linelets], nbpts_linelets, 0.0, normb, h0, hn, 2)
                         else:

                                 b0 = valxI-valxW
                                 b1 = valyI-valyW
                                 b2 = valzI-valzW
                                 normb = numpy.sqrt(b0*b0+b1*b1+b2*b2)

                                 if hn < 0: hn = 0.6*normb

                                 # _stretch(linelets[shift : shift + nbpts_linelets],nbpts_linelets,0.0,normb, 1.0e-6, 0.6*normb, 2)
                                 _stretch(linelets[shift : shift + nbpts_linelets], nbpts_linelets, 0.0, normb, h0, hn, 2)


                         Internal.createUniqueChild(s, 'CoordinateN_ODE', 'DataArray_t',
                                                    linelets[shift : shift + Nbpts_D*nbpts_linelets])
                         Internal.createUniqueChild(s, 'VelocityT_ODE', 'DataArray_t',
                                                    linelets[shift + 1*Nbpts_D*nbpts_linelets :  shift + 2*Nbpts_D*nbpts_linelets ])
                         Internal.createUniqueChild(s, 'TurbulentSANuTilde_ODE', 'DataArray_t',
                                                    linelets[shift + 2*Nbpts_D*nbpts_linelets :  shift + 3*Nbpts_D*nbpts_linelets ])
                         Internal.createUniqueChild(s, 'Psi_ODE', 'DataArray_t',
                                                    numpy.zeros(Nbpts_D*nbpts_linelets, numpy.float64))
                         Internal.createUniqueChild(s, 'Matm_ODE', 'DataArray_t',
                                                    linelets[shift + 3*Nbpts_D*nbpts_linelets :  shift + 4*Nbpts_D*nbpts_linelets ])
                         Internal.createUniqueChild(s, 'Mat_ODE', 'DataArray_t',
                                                    linelets[shift + 4*Nbpts_D*nbpts_linelets :  shift + 5*Nbpts_D*nbpts_linelets ])
                         Internal.createUniqueChild(s, 'Matp_ODE', 'DataArray_t',
                                                    linelets[shift + 5*Nbpts_D*nbpts_linelets :  shift + 6*Nbpts_D*nbpts_linelets ])

                         Internal.createUniqueChild(s, 'alphasbeta_ODE', 'DataArray_t',
                                                    numpy.zeros(Nbpts_D, numpy.float64))
                         Internal.createUniqueChild(s, 'index_ODE', 'DataArray_t',
                                                    numpy.zeros(Nbpts_D, numpy.float64))

                         # Internal.createUniqueChild(s, 'Matm_ODE', 'DataArray_t',
                         #                            numpy.ones(Nbpts_D*45, numpy.float64)*2)

                         # Internal.createUniqueChild(s, 'Mat_ODE', 'DataArray_t',
                         #                            numpy.ones(Nbpts_D*45, numpy.float64)*3)
                         
                         # Internal.createUniqueChild(s, 'Matp_ODE', 'DataArray_t',
                         #                            numpy.ones(Nbpts_D*45, numpy.float64)*4)

                         # Internal.createUniqueChild(s, 'alphasbeta_ODE', 'DataArray_t',
                         #                            numpy.ones(Nbpts_D, numpy.float64)*5)

                         # Internal.createUniqueChild(s, 'index_ODE', 'DataArray_t',
                         #                            numpy.ones(Nbpts_D, numpy.float64)*6)


                         shift = shift + Nbpts_D*(nbpts_linelets*nfields_lin + 1)

      return None

#==============================================================================
def _updateGradPInfo(t, tc, metrics, type='IBCD'): 

  varGrad = ['Density', 'Temperature']
  variablesIBC = ["Density", "Temperature", "gradxDensity", "gradyDensity", "gradzDensity", "gradxTemperature", "gradyTemperature", "gradzTemperature"]

  for var in varGrad:
   C._rmVars(t, 'gradx{}'.format(var))
   C._rmVars(t, 'grady{}'.format(var))
   C._rmVars(t, 'gradz{}'.format(var))

   _computeGrad(t, metrics, varGrad)
   for v in variablesIBC: C._cpVars(t, 'centers:'+v, tc, v)
   Xmpi._setInterpTransfers(t, tc, variables=variablesIBC, variablesIBC=[], compact=0, type='ID')
   Xmpi._setInterpTransfers(t, tc, variables=[], variablesIBC=variablesIBC, compact=0, type=type)
   C._rmVars(tc, variablesIBC)

  return None
#==============================================================================
def warmup(t, tc, graph=None, infos_ale=None, Adjoint=False, tmy=None, list_graph=None, Padding=None, verbose=0):
    Re  =-1
    Lref= 1.
    gradP      =False
    isWireModel=False
    if tc is not None:
        base       = Internal.getBases(tc)[0]
        solverIBC  = Internal.getNodeFromName(base ,'.Solver#IBCdefine')
        if solverIBC is not None:
            gradP      = eval(Internal.getValue(Internal.getNodeFromName(solverIBC, 'isgradP')))
            isWireModel= eval(Internal.getValue(Internal.getNodeFromName(solverIBC, 'isWireModel')))

            Re  = Internal.getValue(Internal.getNodeFromName(solverIBC, 'Reref'))
            Lref= Internal.getValue(Internal.getNodeFromName(solverIBC, 'Lref'))

    graphInvIBCD_WM = None
    if isinstance(graph, list):
        #test pour savoir si graph est une liste de dictionnaires (explicite local)
        #ou juste un dictionnaire (explicite global, implicite)
        grapheliste=True
    else:
        grapheliste=False
    
    if graph is not None and not grapheliste:
        procDict  = graph['procDict']
        graphID   = graph['graphID']
        graphIBCD = graph['graphIBCD']
        if isWireModel: graphInvIBCD_WM  = Cmpi.computeGraph(tc, type='INV_IBCD', procDict=procDict)
    elif graph is not None and grapheliste:
        procDict  = graph[0]['procDict']
        graphID   = graph[0]['graphID']
        graphIBCD = graph[0]['graphIBCD']
        if isWireModel: graphInvIBCD_WM = Cmpi.computeGraph(tc, type='INV_IBCD', procDict=procDict)
    else: 
        procDict=None; graphID=None; graphIBCD=None; graphInvIBCD_WM=None
        
    if isWireModel and graphInvIBCD_WM is None:
        print("Wire Model REQUIRES graphInvIBCD...exiting")
        exit()

    # compute info linelets
    nbpts_linelets = 0
    if tc is not None:
      if Re > 0: #Adaptive method
        h0, hn, nbpts_linelets = computeLineletsInfo(tc, Re=Re, Lref=Lref, q=1.1)
      elif Re == 0: 
        h0, hn, nbpts_linelets = computeLineletsInfo2(tc, q=1.1)
      else: #Alferez' og method
        h0, nbpts_linelets = 1.e-6, 45

    first = Internal.getNodeFromName1(t, 'NbptsLinelets')
    if first is None: Internal.createUniqueChild(t, 'NbptsLinelets', 'DataArray_t', value=nbpts_linelets)

    # Get omp_mode
    ompmode = FastC.OMP_MODE
    node = Internal.getNodeFromName1(t, '.Solver#define')
    if node is not None:
        node = Internal.getNodeFromName1(node, 'omp_mode')
        if  node is not None: ompmode = Internal.getValue(node)

    # Reordone les zones pour garantir meme ordre entre t et tc
    FastC._reorder(t, tc)

    # Construction param_int et param_real des zones
    FastC._buildOwnData(t, Padding)

    t = Internal.rmNodesByName(t, 'NbptsLinelets')
    
    # determination taille des zones a integrer (implicit ou explicit local)
    #evite probleme si boucle en temps ne commence pas a it=0 ou it=1. ex: range(22,1000)

    tmp   = Internal.getNodeFromName1(t  , '.Solver#ownData')
    dtlocPy = Internal.getNodeFromName1(tmp, '.Solver#dtloc')  # noeud
    dtloc   = Internal.getValue(dtlocPy)# tab numpy
    nssiter = dtloc[0]   

    if isinstance(graph, list):
        #test pour savoir si graph est une liste de dictionnaires (explicite local)
        #ou juste un dictionnaire (explicite global, implicite)
        grapheliste=True
    else:
        grapheliste=False
    
    if graph is not None and not grapheliste:
        procDict  = graph['procDict']
        graphID   = graph['graphID']
        graphIBCD = graph['graphIBCD']
        if isWireModel: graphInvIBCD_WM = Cmpi.computeGraph(tc, type='INV_IBCD', procDict=procDict)
    elif graph is not None and grapheliste:  ### Dans warmup tous les transferts doivent etre faits
        procDict  = graph[nssiter-2]['procDict']   ### On va chercher le graphe a nssiter-2 car a cette ssite tous les transferts
        graphID   = graph[nssiter-2]['graphID']    ### sont faits pour le schema a pas de temps local
        graphIBCD = graph[nssiter-2]['graphIBCD']
        if isWireModel: graphInvIBCD_WM  = Cmpi.computeGraph(tc, type='INV_IBCD', procDict=procDict)
    else: 
        procDict=None; graphID=None; graphIBCD=None; graphInvIBCD_WM = None
    zones = Internal.getZones(t)
    f_it = FastC.FIRST_IT
    if FastC.HOOK is None: FastC.HOOK = FastC.createWorkArrays__(zones, dtloc, f_it ); FastC.FIRST_IT = f_it

    # allocation d espace dans param_int pour stockage info openmp
    FastC._build_omp(t) 

    # alloue metric: tijk, ventijk, ssiter_loc
    # init         : ssiter_loc
    metrics = allocate_metric(t)

    # Contruction BC_int et BC_real pour CL
    FastC._BCcompact(t) 

    # Contruction Conservative_int
    FastC._Fluxcompact(t) 
 
    #determination taille des zones a integrer (implicit ou explicit local)
    #evite probleme si boucle en temps ne commence pas a it=0 ou it=1. ex: range(22,1000)
    for nstep in range(1, int(dtloc[0])+1):
        hook1 = FastC.HOOK.copy()
        hook1.update(  FastC.fastc.souszones_list(zones, metrics, FastC.HOOK, 1, nstep, verbose) )

    _init_metric(t, metrics, hook1)

    ssors = allocate_ssor(t, metrics, hook1, ompmode)

    # compact + align + init numa
    rmConsVars=True
    adjoint   =Adjoint

    t, FIRST_IT, zones2compact = FastC.createPrimVars(t, ompmode, rmConsVars, adjoint, gradP, isWireModel)
    FastC.HOOK['FIRST_IT']= FIRST_IT

    #compactage des champs en fonction option de calcul  
    count = -1
    if ompmode == 1: count = 0          
    for data in zones2compact:
        if ompmode == 1: count += 1
        zone    = data[0]
        varnames= data[1]
        for fields in varnames:
            FastC._compact(zone, fields=fields, mode=count, dtloc=dtlocPy)

    #corection pointeur ventijk si ale=0: pointeur Ro perdu par compact.
    zones = Internal.getZones(t) # car create primvar rend zones caduc
    c   = 0
    ale = 0
    for z in zones:
        motion = 'none'
        b = Internal.getNodeFromName2(z, 'motion')
        if b is not None: motion = Internal.getValue(b)
        if motion == 'none':
            sol = Internal.getNodeFromName1(z, 'FlowSolution#Centers')
            ro = Internal.getNodeFromName1(sol, 'Density')
            metrics[c][2] = ro[1]
        elif motion == 'deformation': ale = 2
        else: ale = 1
        c += 1

    #
    # mise a jour vitesse entrainememnt
    #
    if ale == 1 and infos_ale is not None:
        print("ale actif. Teta et tetap=", infos_ale)
        teta = infos_ale[0];  tetap = infos_ale[1]
        FastC._motionlaw(t, teta, tetap)
        _computeVelocityAle(t,metrics)
    elif ale == 2:
        first = Internal.getNodeFromName1(t, 'Time')
        if first is not None: time = Internal.getValue(first)
        else: time = 0.
        R._evalPosition(t, time)
        R._evalGridSpeed(t, time)
        copy_velocity_ale(t, metrics)
    #
    # Compactage arbre transfert
    #
    if tc is not None:      
      nbpts_linelets = 0
      if Re > 0: #Adaptive method
        # h0, hn, nbpts_linelets = computeLineletsInfo(tc, Re=Re, Lref=Lref, q=1.1)
        # print(h0, hn, nbpts_linelets)
        _createTBLESA2(t, tc, h0=h0, hn=hn, nbpts_linelets=nbpts_linelets)
      elif Re == 0: 
        # h0, hn, nbpts_linelets = computeLineletsInfo2(tc, q=1.1)
        # print(h0, hn, nbpts_linelets)
        _createTBLESA2(t, tc, h0=h0, hn=hn, nbpts_linelets=nbpts_linelets)
      else: #Alferez' og method
        # h0, nbpts_linelets = 1.e-6, 45
        # print(h0, -1, nbpts_linelets)
        _createTBLESA2(t, tc, h0=h0, hn=-1, nbpts_linelets=nbpts_linelets)

      X.miseAPlatDonorTree__(zones, tc, graph=graph,list_graph=list_graph, nbpts_linelets=nbpts_linelets)

      FastC.HOOK['param_int_tc'] = Internal.getNodeFromName1( tc, 'Parameter_int')[1]
      param_real_tc               = Internal.getNodeFromName1( tc, 'Parameter_real')
      if param_real_tc is not None: FastC.HOOK['param_real_tc']= param_real_tc[1]
    else:
        FastC.HOOK['param_real_tc'] = None
        FastC.HOOK['param_int_tc']  = None 


    if ssors is not []:
        FastC.HOOK['ssors'] = ssors
    else:
        FastC.HOOK['ssors'] = None

    # Compactage arbre moyennes stat
    #
    if tmy is not None:
        sol = Internal.getNodesFromName3(tmy, 'FlowSolution#Centers')
        var = Internal.getNodesFromType1(sol[0] , 'DataArray_t')
        varmy=[]
        for v in var: varmy.append('centers:'+v[0])
        FastC._compact(tmy, fields=varmy)

    #
    # remplissage ghostcells
    #
    hook1['lexit_lu'] = 0
    nstep             = 1
    nitrun            = 0
    if infos_ale is not None and len(infos_ale) == 3: nitrun = infos_ale[2]
    timelevel_target = int(dtloc[4])

    _fillGhostcells(zones, tc, metrics, timelevel_target, ['Density'], nstep, ompmode, hook1,graphID, graphIBCD, procDict, isWireModel=isWireModel, graphInvIBCD_WM=graphInvIBCD_WM)

    #
    # initialisation Mut
    #
    if infos_ale is not None and len(infos_ale) == 3: nitrun = infos_ale[2]
    fasts._computePT_mut(zones, metrics, hook1)

    return (t, tc, metrics)

#==============================================================================
# computeStress avec reduction 
#==============================================================================
def _computeStress(t, teff, metrics, xyz_ref=(0.,0.,0.)):
    """Compute efforts in teff.""" 
    ret = PyTree._computeStress(t, teff, metrics, xyz_ref)
    ret = numpy.array(ret, dtype=numpy.float64)
    #ret1 = numpy.empty(ret.shape, dtype=numpy.float64)
    ret1 = numpy.zeros(ret.shape, dtype=numpy.float64)
    #Cmpi.Allreduce(ret, ret1, op=Cmpi.SUM)
    Cmpi.KCOMM.Allreduce(ret, ret1, op=Cmpi.SUM)
    return ret1.tolist()

#==============================================================================
# display CPU efficiency diagnostic
#==============================================================================
def display_cpu_efficiency(t, mask_cpu=0.08, mask_cell=0.01, diag='compact', FILEOUT='listZonesSlow.dat', FILEOUT1='diagCPU.dat', RECORD=None):

    rank = Cmpi.rank
    names   = FILEOUT.split('.')
    fileout = names[0]+str(rank)+'.'+names[1]
    names1   = FILEOUT1.split('.')
    fileout1 = names1[0]+str(rank)+'.'+names1[1]

    if RECORD is not None: 
      tape = PyTree.display_cpu_efficiency(t, mask_cpu=mask_cpu, mask_cell=mask_cell, diag=diag, FILEOUT=fileout, FILEOUT1=fileout1, RECORD=RECORD)
      return tape
    else: 
      PyTree.display_cpu_efficiency(t, mask_cpu=mask_cpu, mask_cell=mask_cell, diag=diag, FILEOUT=fileout, FILEOUT1=fileout1, RECORD=RECORD)
      return None

# For periodic unsteady chimera join, parameter must be updated peridicaly 
#==============================================================================
def _UpdateUnsteadyJoinParam(t, tc, tc_skel, graph, omega, timelevelInfos, split='single',root_steady='tc_steady', root_unsteady='tc_', dir_steady='.', dir_unsteady='.', init=False, layer='Python', Rotors=['Base02','Base04','Base06'], Stators=['Base01','Base03','Base05','Base07']):

    #on cree/initialise le dico infos graph
    if graph =={} or graph is None:
       graph= {'graphID':None, 'graphIBCD':None, 'procDict':None, 'procList':None}

    #on cree les noeud infos insta pour chimere perio s'il n'existe pas 
    TimeLevelOpts=['TimeLevelMotion','TimeLevelTarget','Iteration'] 
    for opt in TimeLevelOpts:
       tmp = Internal.getNodeFromName1(t, opt)
       if tmp is None: Internal.createUniqueChild(t, opt, 'DataArray_t', value=0)

    #mise a jour des info temporelle
    if not init:
      Internal.getNodeFromName1(t, 'TimeLevelMotion')[1][0] +=1
      Internal.getNodeFromName1(t, 'TimeLevelTarget')[1][0] +=1
      Internal.getNodeFromName1(t, 'Iteration'      )[1][0] +=1
 
    timelevel_motion = Internal.getNodeFromName1(t, 'TimeLevelMotion')[1][0]
    timelevel_target = Internal.getNodeFromName1(t, 'TimeLevelTarget')[1][0]
    iteration        = Internal.getNodeFromName1(t, 'Iteration')[1][0]

    timelevel_period = timelevelInfos["TimeLevelPeriod"]
    timelevel_360    = timelevelInfos["TimeLevel360"]
    timelevel_perfile= timelevelInfos["TimeLevelPerFile"]
    timelevel_axeRot = timelevelInfos["TimeLevelRotationAxis"]

    No_period = timelevel_motion//timelevel_period 
    #
    #target in no more in tc; need need data in a new file
    #

    if timelevel_target == timelevel_perfile or tc is None or timelevel_motion%timelevel_period == 0: 

       graph= {'graphID':None, 'graphIBCD':None, 'procDict':None, 'procList':None}
       tc =None; tc_skel =None

       rank = Cmpi.rank
       tmp  = No_period*timelevel_period
       root = timelevel_perfile + ( (timelevel_motion - tmp)//timelevel_perfile)*timelevel_perfile
       if root > timelevel_period : root=timelevel_period ### Comme 8000 pas multiple de 60 on force le load de tc_8000
    
       if split == 'single':

          # steady part
          FILE = dir_steady +'/'+ root_steady + '.cgns'
          if os.access(FILE, os.F_OK): 
             t0 = timeit.default_timer()
             tc = Cmpi.convertFile2SkeletonTree(FILE)
             tc_skel= Internal.copyRef(tc)
             tc = Cmpi.readZones(tc, FILE, rank=rank)
             tc = Cmpi.convert2PartialTree(tc, rank=rank)
             t1 = timeit.default_timer()
             #print( "cout steady= ", t1-t0, 'rank=', rank)
             #sys.stdout.flush
          # unsteady part
          FILE = dir_unsteady +'/'+root_unsteady+str(root)+'.cgns'
          if rank==0: print('chargement tc: timelevel_target=',timelevel_target, 'file=',FILE, root)
          if os.access(FILE, os.F_OK): 
             t0=timeit.default_timer()
             tc_inst = Cmpi.convertFile2SkeletonTree(FILE)
             tc_inst_skel= Internal.copyRef(tc_inst)
             tc_inst = Cmpi.readZones(tc_inst, FILE, rank=rank)
             tc_inst = Cmpi.convert2PartialTree(tc_inst, rank=rank)
             t1=timeit.default_timer()
             #print( "cout unsteady= ", t1-t0, 'rank=', rank)
             #sys.stdout.flush

       else:
          # steady part
          t0=timeit.default_timer()
          FILE_Skeleton = dir_steady +'/'+ root_steady + 'skeleton.cgns'
          if os.access(FILE_Skeleton, os.F_OK): 
             tc_skel = C.convertFile2PyTree(FILE_Skeleton)

          t1=timeit.default_timer()
          cpu_skelS = t1 -t0
          t0=t1

          FILE = dir_steady +'/'+ root_steady + str(rank) +'.cgns'
          if os.access(FILE, os.F_OK): 
             tc = C.convertFile2PyTree(FILE)
          t1=timeit.default_timer()
          cpu_S = t1 -t0
          t0=t1
          #print( "cout steady= ", t1-t0, 'rank=', rank)
          #sys.stdout.flush

          # unsteady part: toto_Nit_skeleton.cgns
          t0=timeit.default_timer()
          FILE_Skeleton = dir_unsteady +'_'+str(root) +'/'+ root_unsteady +str(root) + '_skeleton.cgns'
          if os.access(FILE_Skeleton, os.F_OK): 
             tc_inst_skel = C.convertFile2PyTree(FILE_Skeleton)

          t1=timeit.default_timer()
          cpu_skelU = t1 -t0
          t0=t1
          # toto_Nit_proc.cgns
          FILE = dir_unsteady+'_'+str(root) +'/'+root_unsteady+str(root)+'_'+str(rank)+'.cgns'
          if os.access(FILE, os.F_OK): 
             tc_inst = C.convertFile2PyTree(FILE)

          t1=timeit.default_timer()
          cpu_U = t1 -t0
          t0=t1
          #print( "cout unsteady= ", t1-t0, 'rank=', rank)
          #sys.stdout.flush

       #mise a zero timelevel_target en fin de fichier ou fin de periode azymutale
       if timelevel_target == timelevel_perfile or timelevel_motion%timelevel_period == 0:
          timelevel_target = 0
          Internal.getNodeFromName1(t, 'TimeLevelTarget')[1][0] = timelevel_target

       if rank==0: print("timelevel_motion= ", timelevel_motion,"timelevel_target= ", timelevel_target)

       #
       #timelevel_motion larger than calculated peridicity; need to modify angle of rotation for azymuth periodicity
       #
       if timelevel_motion >= timelevel_period:

           No_period     = timelevel_motion//timelevel_period
           iteration_loc = timelevel_motion - No_period*timelevel_period 

           bases  = Internal.getNodesFromType1(tc_inst , 'CGNSBase_t')       # noeud

           sign =-1
           if omega > 0: sign = 1
           for base in bases:
             if   base[0] in Rotors  : teta = -2*math.pi*timelevel_period/timelevel_360*No_period*sign
             elif base[0] in Stators : teta =  2*math.pi*timelevel_period/timelevel_360*No_period*sign
             zones  = Internal.getNodesFromType1(base , 'Zone_t')       # noeud
             for z in zones:
               angles = Internal.getNodesFromName2(z, 'RotationAngle')
               for angle in angles: angle[1][:]= angle[1][:] + teta*timelevel_axeRot[:]
       else:
            iteration_loc = timelevel_motion

       tc      = Internal.merge( [tc     , tc_inst     ] )
       tc_skel = Internal.merge( [tc_skel, tc_inst_skel] )

       #on reordone les bases par alphabet  au cas ou tc_inst foireux
       for tmp in [tc, tc_skel]:
         Internal._sortByName(tmp,recursive=False)
         bases  = Internal.getBases(tmp)
         cgns   = Internal.getNodeFromName1(tmp,'CGNSLibraryVersion')
         tmp[2] = [cgns]+bases
         #on reordone les rac instat au cas ou tc_inst foireux
         for z in Internal.getZones(tmp):
            _sortByName(z)

       t1=timeit.default_timer()
       cpu_merge = t1 -t0
       t0=t1

       graph['procDict'] = D2.getProcDict(tc_skel)
       graph['procList'] = D2.getProcList(tc_skel, sort=True)
       graph['graphIBCD']= Cmpi.computeGraph(tc_skel, type='IBCD', reduction=False, procDict=graph['procDict'])

       graph['graphID_Steady'],graph['graphID_Unsteady'] = Cmpi.computeGraph(tc_skel, type='ID_Unsteady', reduction=False, procDict=graph['procDict'])

       graph['graphID'] = Cmpi.mergeGraph( graph['graphID_Steady'], graph['graphID_Unsteady'][iteration_loc] )

       t1=timeit.default_timer()
       cpu_graph = t1 -t0
       t0=t1

       if rank==0: print("timelevel_motion= ", timelevel_motion,"timelevel_target= ", timelevel_target,"itGraph=",iteration_loc,"cout Graph=", t1-t0)

       #print("steady",graph['graphID'] )

       # 
       # Reordone les zones pour garantir meme ordre entre t et tc
       FastC._reorder(t, tc)

       # Compactage arbre transfert. Si init, compactage dans warmup plus tard
       #
       if not init: 

          zones = Internal.getZones(t)

          t0=timeit.default_timer()
          X.miseAPlatDonorTree__(zones, tc, graph=graph)
          t1=timeit.default_timer()
          cpu_plat = t1 -t0
          t0=t1
          
          #print("cout skelS=",cpu_skelS, 'skelU=', cpu_skelU, "S=",cpu_S, 'U=', cpu_U, 'merge', cpu_merge, 'graph', cpu_graph, 'plat', cpu_plat, 'rank=', rank)
          #sys.stdout.flush

          FastC.HOOK['param_int_tc'] = Internal.getNodeFromName1( tc, 'Parameter_int')[1]
          param_real_tc              = Internal.getNodeFromName1( tc, 'Parameter_real')
          if param_real_tc is not None: FastC.HOOK['param_real_tc']= param_real_tc[1]

    else:
      if layer =='Python':

        if timelevel_motion >= timelevel_period:
            No_period = timelevel_motion//timelevel_period
            iteration_loc = timelevel_motion - No_period*timelevel_period 
        else:
            iteration_loc = timelevel_motion
    
        rank = Cmpi.rank

        t0=timeit.default_timer()
        graph['graphID'] = Cmpi.mergeGraph( graph['graphID_Steady'], graph['graphID_Unsteady'][iteration_loc] )
        #graph['graphID']  = Cmpi.computeGraph(tc_skel, type='ID', reduction=False, procDict=graph['procDict'], it=iteration_loc)
        t1=timeit.default_timer()
        if rank==0:
          print('calcul du graph it', iteration_loc, "cout graphID= ", t1-t0)
          sys.stdout.flush

    if timelevel_motion > timelevel_360:
       #print('remise a ZERO dans updateUnsteady')
       timelevel_motion = 1
       Internal.getNodeFromName1(t, 'TimeLevelMotion')[1][0] = timelevel_motion

    own  = Internal.getNodeFromName1(t, '.Solver#ownData')
    if own is not None: 
       dtloc= Internal.getNodeFromName1(own , '.Solver#dtloc')[1]
       dtloc[3]=timelevel_motion
       dtloc[4]=timelevel_target

    return tc, tc_skel, graph

#====================================================================================
# reordonne  les raccords instationnaires par No croissant  pour chaque zone
#====================================================================================
def _sortByName(t):
    names = []
    for i in t[2]:
        toto = i[0].split('#')
        if len(toto) == 2:
            if toto[1][0].isdigit():
                search =True
                c=1
                while search:
                   if toto[1][c].isdigit():
                      c+=1
                   else: search = False
                names.append(int(toto[1][0:c]))
                #print i[0], toto[1][0:c]
            else:
                names.append(-1)

        else: names.append(-1)
    nodes = t[2]
    zipped = zip(names, nodes)
    zipped = sorted(zipped, key=lambda x: x[0])
    ret = []
    for i in zip(*zipped): ret = i
    t[2] = list(ret)

#====================================================================================
# Calcule la CFL en chaque maille et la met dans l'arbre t pour dtloc instationnaire
#====================================================================================
def computeCFL_dtlocal(t):

    import math
    import Transform.PyTree as T        
    from mpi4py import MPI  
 
    t = C.initVars(t,'centers:CFL', 0.)

    dimPb = 3
    node = Internal.getNodeFromName(t, 'EquationDimension')
    if node is not None: dimPb = Internal.getValue(node)

    liste_BCPeriodiques = []
    
    """
    BCMatch = Internal.getNodesFromType(t, 'GridConnectivity1to1_t')
    for match in BCMatch:
        perio = Internal.getNodesFromName(match, 'GridConnectivityProperty')
        if perio != [] :
            rotCenter  = Internal.getNodeFromName(perio, 'RotationCenter')[1]
            rotAngle   = Internal.getNodeFromName(perio, 'RotationAngle')[1]
            translation = Internal.getNodeFromName(perio, 'Translation')[1]
         
            list1 = [abs(translation).tolist(),abs(rotAngle).tolist(),rotCenter.tolist()]
 
            if (list1 not in liste_BCPeriodiques) :            
                 liste_BCPeriodiques.append(list1)
  
    """

               
    print(liste_BCPeriodiques)
                 

    (t,tc,metrics) = warmup(t,tc=None)

    zones = Internal.getZones(t)

    dim = Internal.getZoneDim(zones[0])
    
    print('dimPb= ', dimPb)
    
    fasts.prep_cfl(zones, metrics,1,1,1)

    #t = Internal.rmGhostCells(t, t, 2)
    """
    C._rmVars(t, 'Temperature_M1')
    C._rmVars(t, 'Temperature_P1')
    C._rmVars(t, 'Density_M1')
    C._rmVars(t, 'Density_P1')
    C._rmVars(t, 'VelocityX_M1')
    C._rmVars(t, 'VelocityX_P1')
    C._rmVars(t, 'VelocityY_M1')
    C._rmVars(t, 'VelocityY_P1')
    C._rmVars(t, 'VelocityZ_M1')
    C._rmVars(t, 'VelocityZ_P1')
    """

    cflmax_glob = 0.0
    cflmin_glob = 1.0e15 

    for z in zones:
        dim = Internal.getZoneDim(z)
        zp = T.subzone(z, (3,3,3), (dim[1]-2,dim[2]-2,dim[3]-2))
        cflmax = C.getMaxValue(zp, 'centers:CFL')
        cflmin = C.getMinValue(zp, 'centers:CFL')
        if cflmax > cflmax_glob:cflmax_glob = cflmax
        if cflmin < cflmin_glob:cflmin_glob = cflmin

    dtmin = 1.0/cflmax
    dtmax = 1.0/cflmin

    dts = numpy.array([dtmin,dtmax])
    dts_rcv = numpy.empty(6, dtype='d')

    MPI.COMM_WORLD.Allgather(dts,dts_rcv)

    dtmax = max(dts_rcv)
    dtmin = min(dts_rcv)
    

    timeLevel_max = math.floor(math.log((dtmax/dtmin))/math.log(2))

    print('dtmin, dtmax, Maxlevel = ', dtmin,',',dtmax,',',timeLevel_max)

    return (t,dtmin,timeLevel_max)

#==============================================================================
# decoupe maillage pour dtloc instationnaire
#==============================================================================
def _decoupe2(t, niveauMax, dtmin):

    import Transform.PyTree as T

    fes = Internal.getNodeFromName2(t,'FlowEquationSet')
    rs = Internal.getNodeFromName2(t,'ReferenceState')

    dtmin = 100000.
    taille_bloc = 25
    bases  = Internal.getNodesFromType1(t     , 'CGNSBase_t') 
    zones = []
    zones_decoupe = []

    dicoTps = {}
    numZone = []

    """
    liste_BCPeriodiques = []
    
    BCMatch = Internal.getNodesFromType(t, 'GridConnectivity1to1_t')
    for match in BCMatch:
        perio = Internal.getNodesFromName(match, 'GridConnectivityProperty')
        if perio != [] :
            rotCenter  = Internal.getNodeFromName(perio, 'RotationCenter')[1]
            rotAngle   = Internal.getNodeFromName(perio, 'RotationAngle')[1]
            translation = Internal.getNodeFromName(perio, 'Translation')[1]
         
            list1 = [abs(translation).tolist(),abs(rotAngle).tolist(),rotCenter.tolist()]
 
            if (list1 not in liste_BCPeriodiques) :            
                 liste_BCPeriodiques.append(list1)
         
    print(liste_BCPeriodiques)
    """     

    for b in bases:
        zones += Internal.getNodesFromType1(b, 'Zone_t')

    somme = 0
    numZone.append(somme)
    
    for z in zones:

        C._addVars(z,'centers:niveaux_temps')

        dim = Internal.getZoneDim(z)
        ni = dim[1]-5
        nj = dim[2]-5
        nk = dim[3]-5

        dimPb = dim[4]
        
        if dimPb==2:nk=1

        nbbloc_i = ni/taille_bloc
        nbbloc_i = int(nbbloc_i)
        if (nbbloc_i*taille_bloc < ni) : nbbloc_i = nbbloc_i + 1
        nbbloc_j = nj/taille_bloc
        nbbloc_j = int(nbbloc_j)
        if (nbbloc_j*taille_bloc < nj) : nbbloc_j = nbbloc_j + 1
        nbbloc_k = nk/taille_bloc
        nbbloc_k = int(nbbloc_k)
        if (nbbloc_k*taille_bloc < nk) : nbbloc_k = nbbloc_k + 1

        somme += nbbloc_i*nbbloc_j*nbbloc_k
        
        numZone.append(somme)

        print('decoupage de la zone', z[0])

        print(nbbloc_i,nbbloc_j,nbbloc_k)
        
         
        for k in range(0,nbbloc_k):
            for j in range(0,nbbloc_j):
                for i in range(0,nbbloc_i):

                    imin = i*taille_bloc + 3
                    jmin = j*taille_bloc + 3
                    kmin = k*taille_bloc + 3

                    imax = min((i+1)*taille_bloc+1,ni)
                    jmax = min((j+1)*taille_bloc+1,nj)
                    kmax = min((k+1)*taille_bloc+1,nk)

                    if dimPb==2: kmax=1
                   
                    zp = T.subzone(z,(imin,jmin,kmin),(imax,jmax,kmax))

                    cflmax_loc = C.getMaxValue(zp, 'centers:CFL')
                    dtmin_loc = 1./cflmax_loc
                    

                    niveauTps = math.log(((2**niveauMax)*dtmin)/dtmin_loc)/math.log(2)
     
 
                    fs = Internal.getNodeFromName(t,'FlowSolution#Centers')
                    nt = Internal.getNodeFromName(fs,'niveaux_temps')[1]
                    for i_ in range(imin, imax):
                        for j_ in range(jmin, jmax):
                            for k_ in range(kmin, kmax):
                                 nt[i,j,k] = float(niveauTps)


    return t
