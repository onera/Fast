# - compute (pyTree) -
# - TGVa
#
#mpirun -x OMP_NUM_THREADS -np 8 python tgv_mpi.py
#
import Generator.PyTree as G
import FastC.PyTree as FastC
import FastS.PyTree as FastS
import FastS.Mpi as FastSmpi
import Converter.PyTree as C
import Distributor2.PyTree as D2
import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import sys
import math, numpy

# Grandeurs
pi = math.pi
L0 = 1
U0 = 1. ; RO0 = 1. ; P0 = 71.345295521 ; T0 = 1.
L1 = 2*pi*L0
r  = P0 /(RO0*T0)
K0 = RO0*U0*U0/16.
O0 = r*T0

# Grille
N = 33

rank = Cmpi.rank; size = Cmpi.size
NP =size

mpi_i = 2
mpi_j = 2
mpi_k = 2

ipos = int(rank%mpi_i)
jpos = int((rank-ipos)/mpi_i%mpi_j)
kpos = int((rank-ipos-jpos*mpi_i)/(mpi_i*mpi_j))

#print "rank", rank, "pos", ipos,jpos,kpos
# Bloc avec ghost cells
hx = L1/(N-1.)/mpi_i
hy = L1/(N-1.)/mpi_j
hz = L1/(N-1.)/mpi_k

xdeb = -L1*0.5 + ipos*L1/mpi_i  -2*hx
ydeb = -L1*0.5 + jpos*L1/mpi_j  -2*hy
zdeb = -L1*0.5 + kpos*L1/mpi_k  -2*hz

a = G.cart( (xdeb, ydeb, zdeb) , (hx,hy,hz), (N+4,N+4,N+4))
a = D2.addProcNode(a, rank)

# initialisation
C._initVars(a, '{centers:Density}=(%f+%f*(cos(2.*{centers:CoordinateX}/%f)+cos(2.*{centers:CoordinateY}/%f))*(cos(2.*{centers:CoordinateZ}/%f)+2.))/%f'%(P0,K0,L0,L0,L0,O0))
C._initVars(a, '{centers:VelocityX}= %f*sin({centers:CoordinateX}/%f)*cos({centers:CoordinateY}/%f)*cos({centers:CoordinateZ}/%f)'%(U0,L0,L0,L0))
C._initVars(a, '{centers:VelocityY}= -%f*cos({centers:CoordinateX}/%f)*sin({centers:CoordinateY}/%f)*cos({centers:CoordinateZ}/%f)'%(U0,L0,L0,L0))
C._initVars(a, '{centers:VelocityZ} = 0.')
C._initVars(a, '{centers:Temperature}=%f'%T0)

t = C.newPyTree(['Base']) ; t[2][1][2] += [a]
t[2][1][2][0][0] ='cart'+str(rank)
t = C.addState(t, 'GoverningEquations', 'NSLaminar')
t = C.addState(t, MInf=0.1, ReInf=1600, adim='adim2funk')

#periodicite par BC si pas de decoupe MPI dans une direction
if mpi_k ==1 :
   t[2][1][2][0] = C.addBC2Zone(t[2][1][2][0], 'period kmin', 'BCautoperiod', 'kmin')
   t[2][1][2][0] = C.addBC2Zone(t[2][1][2][0], 'period kmax', 'BCautoperiod', 'kmax')
if mpi_j ==1 :
   t[2][1][2][0] = C.addBC2Zone(t[2][1][2][0], 'period jmin', 'BCautoperiod', 'jmin')
   t[2][1][2][0] = C.addBC2Zone(t[2][1][2][0], 'period jmax', 'BCautoperiod', 'jmax')
if mpi_i ==1 :
   t[2][1][2][0] = C.addBC2Zone(t[2][1][2][0], 'period imin', 'BCautoperiod', 'imin')
   t[2][1][2][0] = C.addBC2Zone(t[2][1][2][0], 'period imax', 'BCautoperiod', 'imax')

tc = C.node2Center(t)
z = tc[2][1][2][0]
C._rmVars(z, 'FlowSolution')

racs =[]
if mpi_i != 1:
   #Imin
   mpi_rac = rank -1
   if ipos ==0: mpi_rac = rank + mpi_i-1
   racs.append(['imin', mpi_rac ])
   #Imax
   mpi_rac = rank +1
   if ipos == mpi_i-1: mpi_rac = rank - mpi_i +1
   racs.append(['imax', mpi_rac ])

if mpi_j != 1:
   #Jmin
   mpi_rac = rank - mpi_i
   if jpos ==0: mpi_rac = rank + (mpi_j-1)*mpi_i
   racs.append(['jmin', mpi_rac ])
   #Jmax
   mpi_rac = rank + mpi_i
   if jpos == mpi_j-1: mpi_rac = rank - (mpi_j -1)*mpi_i
   racs.append(['jmax', mpi_rac ])

if mpi_k != 1:
   #Kmin
   mpi_rac = rank - mpi_i*mpi_j
   if kpos ==0: mpi_rac = rank + (mpi_k-1)*mpi_i*mpi_j
   racs.append(['kmin', mpi_rac ])
   #Kmax
   mpi_rac = rank + mpi_i*mpi_j
   if kpos == mpi_k-1: mpi_rac = rank - (mpi_k -1)*mpi_i*mpi_j
   racs.append(['kmax', mpi_rac ])



#if rank == 0: print "racs", racs, ipos, rank
proc ={}
for rac in racs:
   if str(rac[1]) in proc.keys():
       tmp = proc[ str(rac[1]) ]
       proc[ str(rac[1]) ]= tmp+ [rac[0]]
       tmp = proc[ str(rac[1]) ]
   else:
       proc[ str(rac[1]) ]= [rac[0]]

#if rank == 0:  print "proc", proc
#size fenetre raccor sous hypothese domaine cubique
size_fen = 2*((N-1)**2)
c= 0
for key in proc.keys():
   dirs= proc[key]
   nfen= len(dirs)
   #print "dirs", dirs,nfen, rank
   if rank == 0:  print("dirs", dirs,nfen, rank)
   o = Internal.createUniqueChild( z , 'ID_cart'+key, 'ZoneSubRegion_t', 'cart'+key  )
   subRegion = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
   o = Internal.createChild( subRegion[c] , 'ZoneRole', 'DataArray_t', 'Donor')
   o = Internal.createChild( subRegion[c] , 'GridLocation', 'GridLocation_t', 'CellCenter')
   datap = numpy.ones( size_fen*nfen , numpy.float64)
   Internal.createUniqueChild(subRegion[c], 'InterpolantsDonor', 'DataArray_t', datap)
   datap = numpy.ones(size_fen*nfen, dtype=Internal.E_NpyInt)
   Internal.createUniqueChild(subRegion[c], 'InterpolantsType', 'DataArray_t', datap)

   #pointlistreceveur
   datar = numpy.zeros(size_fen*nfen, dtype=Internal.E_NpyInt)
   datad = numpy.zeros(size_fen*nfen, dtype=Internal.E_NpyInt)
   l = 0
   ni   = (N-1+4)
   ninj = ni*ni
   for  dir in dirs:

     if dir == 'imin':
      for k in range(1,N):
       for j in range(1,N):
         for i in range(-1,1):
               datar[l]= 0 + (i+1) + (j+1)*ni + (k+1)*ninj
               datad[l]= 0 + (i+N) + (j+1)*ni + (k+1)*ninj
               l +=1
     elif dir == 'imax':
       for k in range(1,N):
          for j in range(1,N):
            for i in range(N,N+2):
               datar[l]= 0 + (i+1) + (j+1)*ni + (k+1)*ninj
               datad[l]= 0 + (i-N+2) + (j+1)*ni + (k+1)*ninj
               l +=1


     elif dir == 'jmin':
        for j in range(-1,1):
          for k in range(1,N):
            for i in range(1,N):
               datar[l]=  0 + (i+1) + (j  +1)*ni + (k+1)*ninj
               datad[l]=  0 + (i+1) + (j+N  )*ni + (k+1)*ninj
               l +=1
     elif dir == 'jmax':
        for j in range(N,N+2):
          for k in range(1,N):
            for i in range(1,N):
               datar[l]=  0 + (i+1) + (j  +1)*ni + (k+1)*ninj
               datad[l]=  0 + (i+1) + (j-N+2)*ni + (k+1)*ninj
               l +=1
     elif dir == 'kmin':
       for k in range(-1,1):
          for j in range(1,N):
            for i in range(1,N):
               datar[l]=  0 + (i+1) + (j+1)*ni + (k  +1)*ninj
               datad[l]=  0 + (i+1) + (j+1)*ni + (k+N  )*ninj
               l +=1
     elif dir == 'kmax':
        for k in range(N,N+2):
          for j in range(1,N):
            for i in range(1,N):
               datar[l]=  0 + (i+1) + (j+1)*ni + (k  +1)*ninj
               datad[l]=  0 + (i+1) + (j+1)*ni + (k-N+2)*ninj
               l +=1

   Internal.createUniqueChild(subRegion[c], 'PointListDonor', 'DataArray_t', datar)
   Internal.createUniqueChild(subRegion[c], 'PointList', 'DataArray_t', datad)
   c +=1

#graph
procList=[]
procDict={}
graphID ={}
procs={}
for i in range(size):
  proc ={}
  procDict['cart'+str(i)]= i
  procList.append(['cart'+str(i)])
  rank =i
  
  ipos = int(rank%mpi_i)
  jpos = int((rank-ipos)/mpi_i%mpi_j)
  kpos = int((rank-ipos-jpos*mpi_i)/(mpi_i*mpi_j))
  
  racs =[]
  if mpi_i != 1:
     #Imin
     mpi_rac = rank -1
     if ipos ==0: mpi_rac = rank + mpi_i-1
     racs.append(['imin', mpi_rac ])
     #Imax
     mpi_rac = rank +1
     if ipos == mpi_i-1: mpi_rac = rank - mpi_i +1
     racs.append(['imax', mpi_rac ])

  if mpi_j != 1:
     #Jmin
     mpi_rac = rank - mpi_i
     if jpos ==0: mpi_rac = rank + (mpi_j-1)*mpi_i
     racs.append(['jmin', mpi_rac ])
     #Jmax
     mpi_rac = rank + mpi_i
     if jpos == mpi_j-1: mpi_rac = rank - (mpi_j -1)*mpi_i
     racs.append(['jmax', mpi_rac ])

  if mpi_k != 1:
     #Kmin
     mpi_rac = rank - mpi_i*mpi_j
     if kpos ==0: mpi_rac = rank + (mpi_k-1)*mpi_i*mpi_j
     racs.append(['kmin', mpi_rac ])
     #Kmax
     mpi_rac = rank + mpi_i*mpi_j
     if kpos == mpi_k-1: mpi_rac = rank - (mpi_k -1)*mpi_i*mpi_j
     racs.append(['kmax', mpi_rac ])

  for rac in racs:
     if rac[1] not in proc.keys():
       proc[ rac[1] ]= [ 'cart'+str(rac[1]) ]


  procs[i]= proc
 
 
graph={}
graph['graphID']=procs
graph['graphIBCD']={}
graph['procDict']=procDict
graph['procList']=procList

rank = Cmpi.rank

if rank ==0: print("proc", graph)
#for s in subRegions:
#   print "rac",s[0], rank



# Numerics
numb = {}
numb["temporal_scheme"] = "explicit"
numb["modulo_verif"]    = 50
numz = {}
numz["time_step"]          = 0.003
numz["scheme"]             = "senseur"
FastC._setNum2Zones(t, numz); FastC._setNum2Base(t, numb)

#Initialisation parametre calcul: calcul metric + var primitive + compactage + alignement + placement DRAM
(t, tc, metrics) = FastSmpi.warmup(t, tc, graph)

C.convertPyTree2File(t, 't_'+str(rank)+'.cgns')
C.convertPyTree2File(tc, 'tc_'+str(rank)+'.cgns')

nit = 5000; time = 0.
for it in range(nit):
    FastSmpi._compute(t, metrics, it, tc, graph)
    #FastSmpi._compute(t, metrics, it)
    if rank == 0 and it%50 == 0:
        print('- %d - %f -'%(it, time)); sys.stdout.flush()
        FastS.display_temporal_criteria(t, metrics, it)
    #time += numz['time_step']
t1 = C.node2Center(t)
Cmpi.convertPyTree2File(t1, 't.cgns')
