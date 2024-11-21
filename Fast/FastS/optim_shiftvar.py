# - compute (pyTree) -
# - Lamb vortex / Const -
import Fast.PyTree as Fast
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import Initiator.PyTree as I
import KCore.Adim as Adim
import numpy as ny
import sys
import timeit
import os
import FastS.PyTree

fileInput = 'listZonesSlow0.dat'

Fast.FastI.MX_SYNCHRO= 2001
Fast.FastI.MX_OMP_SIZE_INT= 4000
# NP: 0: Seq, NP> 0: running with mpi, distributed on NP procs
NP = 1
if NP > 0:
    import Converter.Mpi as Cmpi
    import FastS.Mpi as FastS
    rank = Cmpi.rank
else:
    import FastS.PyTree as FastS
    rank=0

runs         =  5   
NIT_internal = 10

CacheI       =128
CacheJ       =  2
CacheK       =  7


FIn = open(fileInput, 'r')
lines_FIn = FIn.readlines()

sizeIJK=[]
for l in lines_FIn:
    mots = l.split(',')
    print(mots[0],  mots[1], mots[2], mots[3])
    sizeIJK.append([ int( mots[1] ),  int( mots[2] ),  int( mots[3]) ] )
FIn.close
typezone = int(mots[4])
iflow    = int(mots[5])
itypcp   = int(mots[6])
kfludom  = int(mots[7])

Nbthread=os.environ['OMP_NUM_THREADS']
platform = 'BDW'

if iflow==1:
    eqs  = 'Euler'
elif iflow==2:
    eqs  = 'NSLaminar'
else:
    eqs  = 'NSTurbulent'

if kfludom ==1:
    scheme  = 'ausmpred'
elif kfludom ==2:
    scheme  = 'senseur'
else:
    scheme  = 'roe_min'

if itypcp ==2:
    temporal ='explicit'
else:
    temporal ='implicit'


shifttab=[]
for k in range(0,71,1):
    shifttab.append(k)

shifttab=shifttab + [128,256,384,512,640,768,1442]

sample = len(shifttab)

for size in sizeIJK:

    sizeI= size[0]+5 
    sizeJ= size[1]+5 
    sizeK= size[2]+5 
    if size[2] == 1:
        sizeK= 2

    f = open('Shift', 'a')
    f.write("#########  Database of the shift parameter used in FastS ############"+'\n')
    f.write("#The first line recall the Cache Block parameters used during the scan: I,J,K"+'\n')
    f.write("#The format is : Nbthread, SizeI-5, SizeJ-5, SizeK-5, TimeCPU, Shiftdata"+'\n')
    f.write(str(CacheI)+","+str(CacheJ)+","+str(CacheK)+'\n')
    f.close


    timeshift=ny.zeros(len(shifttab), ny.float64) 


    print("#########################"+'\n')
    print("Sizes: ",sizeI,sizeJ,sizeK,'\n')
    print("#########################"+'\n')

    mach = 0.7
    adim = Adim.adim1(MInf=mach)
    itshift=0

    a = G.cart((0,0,0), (0.5,0.5,0.25), (sizeI,sizeJ,sizeK))
    a = I.initLamb(a, position=(25.,25.), Gamma=2., MInf=mach, loc='centers')
    #modifie la grille
    zones = Internal.getZones(a)
    for z in zones:
        x = Internal.getNodeFromName2(z,'CoordinateX')[1]
        shape = x.shape
        if typezone==0:
            for k in range(shape[2]):
                x[4,4,k]+=0.00001
            x[4,4,2]+=0.00002
        elif typezone==1:
            for k in range(shape[2]):
                x[4,4,k]+=0.00001

    time20   = ny.empty(sample*runs, ny.float64 )    
    for run in range(runs):
        itshift=0
        for shift in shifttab:

            #print("#########################"+'\n')
            #print("Shift is : ",shift,        '\n')
            print("#########################"+'\n')

            #print(shift,(sizeI-5)*(sizeJ-5))

            t = C.newPyTree(['Base']) ; t[2][1][2] += [a]
            C._addState(t, MInf=mach)
            C._addState(t, 'GoverningEquations', eqs)
            if eqs == 'NSTurbulent': C._initVars(t, 'centers:TurbulentDistance', 1.)
            # Numerics
            ss_iteration   = 1
            numb = {}
            numb["temporal_scheme"]    = temporal
            numb["ss_iteration"]       = ss_iteration
            numb["modulo_verif"]       = 200

            numz = {}
            numz["io_thread"]          =  -1
            numz["time_step"]          = 0.0001
            numz["scheme"]             = scheme
            numz["time_step_nature"]   = "local"
            numz["cache_blocking_I"]   = CacheI
            numz["cache_blocking_J"]   = CacheJ
            numz["cache_blocking_K"]   = CacheK
            numz["shiftvar"]           = shift
            Fast._setNum2Zones(t, numz) ; Fast._setNum2Base(t, numb)

            #Initialisation parametre calcul: calcul metric + var primitive + compactage + alignement + placement DRAM
            (t,tc, metrics)  = FastS.warmup(t,tc=None)

            scale=1.0/(sizeK-5)/(sizeI-5)/(sizeJ-5)/ss_iteration*int(Nbthread)
            # Compute 
            t0=timeit.default_timer()
            FastS._compute(t, metrics,  0, tc=None, NIT=NIT_internal)
            t1=timeit.default_timer()

            #cout = (t1-t0)*scale/NIT_internal
            #print("time per cell",cout,  'run=', run,' shift=',shift)

            time20[ itshift+run*sample ] = FastS.display_cpu_efficiency(t, mask_cpu=-0.45, mask_cell=0.0008, FILEOUT='listZonesSlow_test'+str(rank)+'.dat', RECORD=True)

            print("time per cell",time20[ itshift+run*sample ] ,  'run=', run,' shift=',shift)

            for hook in FastS.PyTree.HOOK: del(hook)
            FastS.PyTree.FIRST_IT = 0
            FastS.PyTree.HOOK = None
            FastS.PyTree.HOOKIBC = None

            timeshift[itshift]+=time20[ itshift + run*sample       ]/float(runs)

            itshift=itshift+1            

    filename='scanPadding_'+str(sizeI-5)+'_'+str(sizeJ-5)+'_'+str(sizeK-5)+'.dat'
    f = open(filename, 'w')
    itshift=0
    for shift in shifttab:
        line = str(shift)
        for r in range(runs):
            line =line +'  '+ str( time20[itshift+sample*r ] )
        line =line + '  '+ str( timeshift[itshift] )
        f.write(line+'\n')
        itshift+=1            
    f.close()

    f = open('Padding.dat', 'a')
    f.write(str(Nbthread)+","+str(sizeI-5)+","+str(sizeJ-5)+","+str(sizeK-5)+","+str(min(timeshift))+","+str(shifttab[ny.argmin(timeshift)])+'\n')
    f.close()

print('option numerique=', temporal, scheme, eqs, typezone)
