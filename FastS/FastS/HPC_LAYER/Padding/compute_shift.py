# - compute (pyTree) -
# - Lamb vortex / Const -
import Generator.PyTree as G
import Converter.PyTree as C
import Initiator.PyTree as I
import Fast.PyTree as Fast
import FastS.PyTree as FastS
import KCore.Adim as Adim
import numpy as ny
import sys
import timeit
import os
import re


# Cree le fichier Shift contenant les infos tailles Pb, tmps cpu et shift. Ce fichier peut etre utiliser pour visualiser (en 2D ou 1D) le tps CPU en fonction des tailles en I,J a l'aide du script plotData.py 

Imin = 119
Imax = 200

Jmin = 6
Jmax = 200

Nbthread=os.environ['OMP_NUM_THREADS']
platform = 'BDW'
CacheI=1000
CacheJ=3
CacheK=3

f = open('Shift', 'a')
f.write("#########  Database of the shift parameter used in FastS ############"+'\n')
f.write("#The first line recall the Cache Block parameters used during the scan: I,J,K"+'\n')
f.write("#The format is : Nbthread, SizeI-5, SizeJ-5, SizeK-5, TimeCPU, Shiftdata"+'\n')
f.write(str(CacheI)+","+str(CacheJ)+","+str(CacheK)+'\n')
f.close

sizeK=80
shifttab=[0,64,128,256,384,512,640,768]
timeshift=ny.empty(len(shifttab), dtype=float)    
timeshift[:] = 1000

for sizeI in xrange(Imin,Imax,1):
    for sizeJ in xrange(Jmin,Jmax,1):

        print "#########################"+'\n'
        print "Sizes: ",sizeI,sizeJ,sizeK,'\n'
        print "#########################"+'\n'

        scale=1.0/(sizeK-5)/(sizeI-5)/(sizeJ-5)/3.0*int(Nbthread)

        mach = 0.7
        adim = Adim.adim1(MInf=mach)
        itshift=0

        for shift in shifttab:

            print "#########################"+'\n'
            print "Shift is : ",shift,        '\n'
            print "#########################"+'\n'

            print shift,(sizeI-5)*(sizeJ-5)
            if shift >= (sizeI-5)*(sizeJ-5): 
                itshift=itshift+1
                continue
            if ( (sizeI <= 40) or (sizeJ <= 40)):
                if shift >= 100:
                    itshift=itshift+1
                    continue
            a = G.cart((0,0,0), (0.5,0.5,0.25), (sizeI,sizeJ,sizeK))
            a = I.initLamb(a, position=(25.,25.), Gamma=2., MInf=mach, loc='centers')
            t = C.newPyTree(['Base']) ; t[2][1][2] += [a]
            del(a)
            t = C.addState(t, 'GoverningEquations', 'NSLaminar')
            t = C.addState(t, MInf=mach)
        #  
# Numerics
            numb = {}
            numb["temporal_scheme"]    = "explicit"
#            numb["ss_iteration"]       = 3
#            numb["modulo_verif"]       = 1

            numz = {}
            numz["io_thread"]          =   1
            numz["time_step"]          = 0.0001
            numz["scheme"]             = "ausmpred"
            numz["cache_blocking_I"]   = CacheI
            numz["cache_blocking_J"]   = CacheJ
            numz["cache_blocking_K"]   = CacheK
            numz["shiftvar"]           = shift
            Fast._setNum2Zones(t, numz) ; Fast._setNum2Base(t, numb)

#Initialisation parametre calcul: calcul metric + var primitive + compactage + alignement + placement DRAM
            (t,tc, metrics)  = FastS.warmup(t,tc=None)

            nit =     1; time = 0.
            time20=ny.empty(nit, dtype=float)    
            timeStep = numz['time_step']
        
            t00=timeit.default_timer()

# Compute 
            for it in xrange(nit):
                t0=timeit.default_timer()
                FastS._compute(t, metrics, 1, tc=None, NIT=40)
                t1=timeit.default_timer()
                print "Compute time",(t1-t0)*scale/40
                time20[it]=(t1-t0)*scale/40
#                if (it%10 == 0):
#                    FastS.display_temporal_criteria(t, metrics, it)                        
#                    print '- %d - %f'%(it, time)
#                    time += timeStep
   #         FastS.itt('pause')
            del(t,metrics,numz,numb)
            t11=timeit.default_timer()
            print "Compute time all",(t11-t00)/40
           
            for a in FastS.HOOK: del(a)
            FastS.FIRST_IT = 0
            FastS.HOOK = None
            FastS.HOOKIBC = None

            timeshift[itshift]=time20[0]
            itshift=itshift+1            
        f = open('Shift', 'a')
        f.write(str(Nbthread)+","+str(sizeI-5)+","+str(sizeJ-5)+","+str(sizeK-5)+","+str(min(timeshift))+","+str(shifttab[ny.argmin(timeshift)])+'\n')
        f.close()
sys.exit()
