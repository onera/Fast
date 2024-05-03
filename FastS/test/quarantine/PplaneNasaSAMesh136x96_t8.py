# - compute (pyTree) -
# - Monobloc laminar flat plate -
import Converter.PyTree as C
import Fast.PyTree as Fast
import FastS.PyTree as FastS
import Connector.PyTree as X
import Converter.Internal as Internal
import KCore.test as test
import sys

t = C.convertFile2PyTree('Pplane_136_96.cgns')
t = C.addState(t, MInf=0.2, ReInf=25.e6, MutSMuInf=15)
#t = C.addState(t, 'GoverningEquations', 'NSLaminar')

mycfl = 10.000

numb = { 'temporal_scheme':'implicit', 'ss_iteration':1, 'modulo_verif':1, 'omp_mode':1 }
numz = {'ssdom_IJK': [600, 30, 20000], 
        'io_thread': -4, 
        'time_step':0.0000152280, 
        'time_step_nature': 'local', 
        'cfl': mycfl, 
        'scheme':'roe', 
        'slope':'o3', 
        'epsi_newton': 0.5}

numz["implicit_solver"] = "gmres"
numz["nb_krylov"] = 50
numz["nb_restart"] =  1
numz["epsi_linear"] = 0.01
numz["nb_relax"] = 2
numz["lu_match"] = 1

Fast._setNum2Zones(t, numz) ; Fast._setNum2Base(t, numb)

(t, tc, metrics) = FastS.warmup(t, None)

zones = Internal.getZones(t)

nit = 30; time = 0.
for it in range(nit):
    FastS._compute(t, metrics, it)
    if it ==0: residu_0 = FastS.display_temporal_criteria(t, metrics, it, gmres=True)
    if it%1 == 0:
        print('it=', it)
        residu = FastS.display_temporal_criteria(t, metrics, it, format='double', gmres=True)
        ratio2 = residu['L2']/residu_0['L2']; ratioo = residu['Loo']/residu_0['Loo']; r= max(ratio2,ratioo)
        for z in zones:
            num = Internal.getNodeFromName1(z, '.Solver#ownData')
            cfl = Internal.getNodeFromName1(num,'Parameter_real')[1][15]
            cflmax = min( 1./r*cfl , 500000000000.)
            Internal.getNodeFromName1(num,'Parameter_real')[1][15]=max(mycfl, cflmax)
            print('cfl', max(mycfl, cflmax))

Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
#C.convertPyTree2File(t, 'gmres.cgns')
test.testT(t, 1)
