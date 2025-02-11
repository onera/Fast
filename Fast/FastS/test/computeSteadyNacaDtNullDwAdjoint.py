# - compute (pyTree) -
#=====================================================================================
#
import Generator.PyTree as G
import Converter.PyTree as C
import FastC.PyTree as FastC
import FastS.PyTree as FastS
import Post.PyTree as P
import Connector.PyTree as X
import Converter.Internal as Internal
import os,math, sys

#=====================================================================================
DOSSIER= './NOMINAL/' #'./VARS_CONS/POINT_100_2/'
SUFFIXE= ''    #'_P100_2_W1_1em5'

NIT = 5
NP = 0               # NP: 0: Seq, NP> 0: running with mpi, distributed on NP procs
rank=0               # sequential 8 cores
it0 = 0

t_m129,tc_m129,ts,graph=FastC.load(DOSSIER+'restart_129'+SUFFIXE+'.cgns',DOSSIER+'tc_m129.cgns', 'ts.cgns', split='single', restart=False)

modulo_verif = 1
mycfl = 0.00000001

numb = {}
numb["temporal_scheme"]    = "implicit"
numb["ss_iteration"]       = 0         # 0 au lieu de 1 pour eviter un "calcul double" ou "double affichage des valeurs"
numb["modulo_verif"]       = modulo_verif
numz = {}
numz["time_step"]          = 0.005
numz["time_step_nature"]   = "local"
numz["cfl"]                = mycfl
numz["scheme"]             = "ausmpred"
numz['cache_blocking_J']   = 1000
numz['cache_blocking_K']   = 1000
numz["io_thread"]          =    1

FastC._setNum2Zones(t_m129, numz); FastC._setNum2Base(t_m129, numb)

#(t_m129, tc_m129, metrics) = FastS.warmup(t_m129,tc_m129, graph=graph)

#====================================================================================

nrec = NIT/modulo_verif
[RoInf, RouInf, RovInf, RowInf, RoeInf, PInf, TInf, cvInf, MInf,
 ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
 Mus, Cs, Ts, Pr] = C.getState(t_m129)

print(" RoInf ",RoInf)
print(" RouInf ",RouInf)
print(" RovInf ",RovInf)
print(" RowInf ",RowInf)
print(" RoeInf ",RoeInf)

QInf   = math.sqrt(RouInf**2+RovInf**2+RowInf**2)
cosAoA = RouInf/QInf
sinAoA = RovInf/QInf
print(" cosAoA ",cosAoA)
print(" sinAoA ",sinAoA)


#====================================================================================
infos_ale=None
(t_m129, tc_m129, metrics) = FastS.warmup(t_m129,tc_m129,graph, infos_ale, Adjoint=True)

FastS.createConvergenceHistory(t_m129, nrec )
t_eff  = FastS.createStressNodes(t_m129,['BCWall'])

zmin1  = C.getMinValue(t_m129,"CoordinateZ")
zmax1  = C.getMaxValue(t_m129,"CoordinateZ")
span   = zmax1-zmin1
print(" pseudo-span ",span)
corde = 1.
nit=NIT

time_step = Internal.getNodeFromName(t_m129, 'time_step')
time_step = Internal.getValue(time_step)

#=======================================================================================

titre='ZONE T="NACA0012 M=0,8 AoA=1.25 I=%f - AUSMP cfl=%d" F=POINT\n'%(nrec,mycfl)
f=open(DOSSIER+'effort_m129_adj'+SUFFIXE+'.dat','w')
f.write('variables="it","CLp","CDp"\n')
f.write(titre)

#=======================================================================================

for it in range(NIT):
    FastS._compute(t_m129, metrics, it, tc_m129, graph)
    #FastS._compute(t_m129, metrics, it, 0, tc_m129, graph)

    if it%modulo_verif == 0:
        print('- %d / %d '%(it+it0, NIT+it0))
        surfinv = 1/(corde*span)
        FastS.display_temporal_criteria(t_m129, metrics, it, format='double')
        eff      = FastS._computeStress(t_m129, t_eff, metrics)
#        eff      = FastS._computeStress(t_m129, t_eff, metrics,cosAoA,sinAoA,surfinv,0)
        drag = (eff[0]*cosAoA+eff[1]*sinAoA)/(corde*span)
        lift = (eff[1]*cosAoA-eff[0]*sinAoA)/(corde*span)
        print('it, lift, drag:', it+it0, lift, drag)
        line="{0} {1} {2} \n".format(it+it0,lift,drag)
        f.write(line)
        f.flush()
#eff2 = FastS._computeStress(t_m129, t_eff, metrics,1)
#eff2 = FastS._computeStress(t_m129, t_eff, metrics,cosAoA,sinAoA,surfinv,1)
eff2 = FastS._compute_dpJ_dpW(t_m129, t_eff, metrics, cosAoA, sinAoA, surfinv)
#eff2 = FastS._computeAdjoint(t_m129, t_eff, metrics, cosAoA, sinAoA, surfinv, 1)

#FastS._compute(t_m129, metrics, it, 1, tc_m129, graph)

#computeAdjoint(t, teff, metrics, nitrun, tc=None, graph=None, cosAoA, sinAoA, surfinv, indFunc):
#computeStress(t, teff, metrics, cosAoA, sinAoA, surfinv, cmpjac)
#===================================================================================

[RoInf, RouInf, RovInf, RowInf, RoeInf, PInf, TInf, cvInf, MInf,
 ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
 Mus, Cs, Ts, Pr] = C.getState(t_m129)
Q2Inf = (RouInf**2+RovInf**2+RowInf**2)/(RoInf**2)
Href = Gamma/(Gamma-1)*PInf/RoInf + .5 * Q2Inf
print('Q2Inf = %f , Href = %f'%(Q2Inf,Href))

#%f(C.isNamePresent(t_eff,'centers:CoefPressure')):
#        print C.getValue( t_eff, 'centers:CoefPressure', 0 )
#      	 C._initVars(t_eff,'{centers:errh}=abs({centers:CoefPressure}-%f)/%f'%(Href,Href))
#	 errtot = P.integ(t_eff,'centers:errh')
#	 print 'erreur relative sur enthalpie totale :', errtot[0]/span
#	 Internal.createUniqueChild(t_eff, 'ErreurEnthalpie', 'DataArray_t', value=errtot[0])

#===================================================================================

#FastC.save(t_m129, DOSSIER+'restart_129'+'.cgns', split='single')
C.convertPyTree2File(t_eff, DOSSIER+'effort_m129_adj'+SUFFIXE+'.cgns')
FastC.save(t_m129, DOSSIER+'restart_129_bis'+'.cgns', split='single')

#===================================================================================
