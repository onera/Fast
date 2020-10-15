# - compute (pyTree) -
# - Lamb vortex [Euler/implicit] -
import Generator.PyTree as G
import Converter.PyTree as C
import Initiator.PyTree as I
import Fast.PyTree as Fast
import FastS.PyTree as FastS
import KCore.Adim as Adim
import KCore.test as test
import Converter.Internal as Internal

mach = 0.5
a = G.cart((0,0,0), (0.25,0.25,0.25), (400,200,2))
a = I.initLamb(a, position=(25.,25.), Gamma=2., MInf=mach, loc='centers')
t = C.newPyTree(['Base']) ; t[2][1][2] += [a]
t = C.addState(t, 'GoverningEquations', 'Euler')
t = C.addState(t, MInf=mach)
# Numerics
numb = {}
numb["temporal_scheme"]    = "implicit"
numb["ss_iteration"]       = 4
numb["modulo_verif"]       = 10 
numz = {}
numz["time_step"]           = 0.01
numz["scheme"]             = "ausmpred"
Fast._setNum2Zones(t, numz) ; Fast._setNum2Base(t, numb)

(t, tc, metrics) = FastS.warmup(t, None)

#iptvar=[]
#zones  = Internal.getNodesFromType(t, 'Zone_t')
#iptro.append( Internal.getNodeFromName2( zones[0], var )[1] )
#for var in ['Temperature','Temperature_M1','Temperature_P1','Density','Density_M1','Density_P1']:
#        print var
#        iptvar.append( Internal.getNodeFromName2( zones[0], var )[1] )

#iptvar[3][:,:,:] = iptvar[3][:,:,:] - 0.103957
#iptvar[4][:,:,:] = iptvar[4][:,:,:] - 0.103957
#iptvar[5][:,:,:] = iptvar[5][:,:,:] - 0.103957
#iptvar[0][:,:,:] = 1.
#iptvar[1][:,:,:] = 1.
#iptvar[2][:,:,:] = 1.
#iptvar[0][:,:,:] = iptvar[0][:,:,:]*iptvar[3][:,:,:]*0.714285
#iptvar[1][:,:,:] = 1./iptvar[3][:,:,:]
#iptvar[2][:,:,:] = 1./iptvar[3][:,:,:]

#C.convertPyTree2File(t, 'funk.cgns')

nit = 102; time = 0.
timeStep = numz['time_step']
for it in range(1,nit):
    #print('it=', it)
    FastS._compute(t, metrics, it)
    if it%10 == 0: FastS.display_temporal_criteria(t, metrics, it)
    time += timeStep
    if it%100 == 0: print(time)

Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
test.testT(t, 1)
