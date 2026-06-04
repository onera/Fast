# FastS + MPI
from . import PyTree
from . import fastc

import numpy

try:
    import Converter.PyTree as C
    import Converter.Mpi as Cmpi
    import FastC.PyTree as FastC
    import os
    import math
except:
    raise ImportError("FastS: requires Converter, Connector and Distributor2 modules.")

try: range = xrange
except: pass

#===============================================================================
# __setInterpTransfers - version optimisee de _setInterpTransfers: arbre t et tc compact,
# moins de python + de C
#
# Warning: inverse storage!
# IN: zones: list zones receveurs
# IN: zoneD: list zones donneurs
# IN: type: ID: interpolation, IBCD: IBCs, ALLD: interp+IBCs
# IN: bcType  0: glissement
#             1: adherence
#             2: loi de paroi log
#             3: loi de paroi Musker
# IN: varType=1,2,3: variablesIBC define (ro,rou,rov,row,roE(,ronutilde)),(ro,u,v,w,t(,nutilde)),(ro,u,v,w,p(,nutilde))
# Adim: KCore.adim1 for Minf=0.1
#===============================================================================
def __setInterpTransfers(zones, zonesD, vars, dtloc, param_int, param_real, it_target,
                         nstep, nitmax, rk, exploc, num_passage, varType=2,
                         graph=None, procDict=None):

    ##for moving IBMs
    isIbmMoving_int  = 0
    # Transferts locaux/globaux
    # Calcul des solutions interpolees par arbre donneur
    # On envoie aussi les indices receveurs pour l'instant
    datas = {}
    nbcomID = param_int[2]
    shift_graph = nbcomID + 2

    for comm_P2P in range(1,param_int[1]+1):
        pt_ech = param_int[comm_P2P + shift_graph]
        dest   = param_int[pt_ech]


        rank  = Cmpi.rank
        #print('dest:', dest, 'rank', rank, flush=True)
        no_transfert = comm_P2P
        if dest == Cmpi.rank: # transfert intra_processus


            FastC.fastc.___setInterpTransfers( zones, zonesD, vars, dtloc, param_int, param_real, it_target, varType, no_transfert, nstep,
                                               nitmax, rk, exploc, num_passage)
            #print('apres ___setInterpTransfers',  flush=True)
        else:
            rank  = Cmpi.rank
            type_transfert =2 #inutile
            infos = FastC.fastc.__setInterpTransfersD(zones, zonesD, vars, dtloc, param_int, param_real, it_target, varType,
                                                      type_transfert, no_transfert, nstep, nitmax, rk, exploc, num_passage, rank,
                                                      isIbmMoving_int)

            #infos = []
            print("ATTENTION fillghost: transferD Reactive")
            if infos != []:
                for n in infos:
                    rcvNode = dest
                    if rcvNode not in datas: datas[rcvNode] = [n]
                    else: datas[rcvNode] += [n]

    # Envoie des numpys suivant le graph
    if graph is not None:
        print("graph is not None dans fillghost", flush=True)
        rcvDatas = Cmpi.sendRecvC(datas, graph)
        #rcvDatas = Cmpi.sendRecv(datas, graph)
    else: rcvDatas = {}

    # Remise des champs interpoles dans l'arbre receveur
    for i in rcvDatas:
        #if Cmpi.rank==0: print(Cmpi.rank, 'recoit de',i, '->', len(rcvDatas[i]), 'nstep=',nstep,flush=True)
        #print(Cmpi.rank, 'recoit de',i, '->', len(rcvDatas[i]), 'nstep=',nstep,flush=True)
        for n in rcvDatas[i]:
            rcvName = n[0]
            field = n[1]

            isSetPartialFields = True
            #if isSetPartialFieldsCheck==1 and field != []:
            #    minfld = numpy.ndarray.min(field[1][0])
            #    maxfld = numpy.ndarray.max(field[1][0])
            #    if maxfld == minfld and maxfld < -1.e05: isSetPartialFields=False

            if isSetPartialFields:
                listIndices = n[2]
                z = zones[rcvName]
                C._setPartialFields(z, [field], [listIndices], loc='centers')

    return None
