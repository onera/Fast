# Entree:
# - un maillage 2D (un seul plan) multibloc
# - ou un maillage 3D multibloc
# avec BCs + raccords + (solution initiale et/ou Reference State)
# Sortie: t.cgns et tc.cgns pour compute.py

import Converter.PyTree as C
import Converter.Mpi    as Cmpi
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.Internal as Internal
import Converter.GhostCells as Ghost
import Connector.PyTree as X
import Fast.PyTree as Fast
import Fast.Utils as Utils
import FastS.ToolboxChimera as TBX
import copy, numpy
import sys
import math

NP = 0 # NB DE PROCS

FILE = 'case_light_perio.cgns' # FICHIER DE MAILLAGE MB
#FILE = 'case_withElsaCL_lisse.cgns' # FICHIER DE MAILLAGE MB

# ROTATION : CENTRE, AXE, ANGLE
(XC0,YC0,ZC0)      =(0 , 0, 0)
(AXISX,AXISY,AXISZ)=(1., 0, 0)

#-----------------------------------------------------------------------------
# PARAMETRES CALCULS A REGLER
#-----------------------------------------------------------------------------
THETADEG       = -11.25 # periodicite azymuthal en degres
TimeLevelPerio = 1000.  # Nombre de dt pour parcourir la periode azym
NGhostCells    = 2
# INPUT: maillage + BC + GC
t = Cmpi.convertFile2SkeletonTree(FILE)

ZoneTarget=['interface','intarface']
ZoneNames = []
for b in Internal.getBases(t):
    for z in Internal.getZones(b):
        for target in ZoneTarget:
            if target in z[0]:
                ZoneNames.append(b[0]+'/'+z[0])
Cmpi._readZones(t, FILE, zoneNames=ZoneNames)
t = Cmpi.convert2PartialTree(t)
#sys.exit()

t = Internal.addGhostCells(t, t, NGhostCells, adaptBCs=1, fillCorner=0)
dim =3
if dim == 2:
    t = T.addkplane(t)
    t = T.contract(t, (0,0,0), (1,0,0), (0,1,0), 0.01)
    t = T.makeDirect(t)

# Preprocessing Chimere periodique
C._initVars(t,"centers:cellN",1.)
t = TBX.modifyBCOverlapsForGhostMesh(t,NGhostCells)

#
#ROTOR/Stator plan de raccord chimere instationnaire
#
etages   = ['Stator','Rotor']
tc_etage = []
tc_RS    ={}
tcyl_RS  ={}
info_PtlistRebuild={}


for etage in etages:
    tmp     = Internal.newCGNSBase(etage)
    base    = Internal.getNodeFromName(t, etage)
    zones   = TBX.ZonePrecond( base , NGhostCells, info_PtlistRebuild )
    tmp[2] += zones
    T._cart2Cyl(tmp,(XC0,YC0,ZC0),(AXISX,AXISY,AXISZ))
    tcyl_RS[etage]= tmp
    tc            = C.node2Center(tmp)
    tc_RS[etage]  = tc

del t

IT_DEB=int(sys.argv[1]);  IT_FIN=int(sys.argv[2])
print('----------------')
print('Plan chimere instationnaire', etages)
print('Nbre de pas de temps a calculer: ', TimeLevelPerio)
print('Calcul des raccord', etages, ' entre les pas de temps ', IT_DEB,' et ', IT_FIN)

THETARAD= THETADEG*math.pi/180.
DTHETA  = THETARAD/float(TimeLevelPerio)

tc_out= TBX.setInterpDataRS( tcyl_RS , tc_RS, THETARAD, DTHETA, IT_DEB, IT_FIN, info_PtlistRebuild,  (XC0,YC0,ZC0), (AXISX,AXISY,AXISZ) )


#ajout du noeud State pour RANS/LES
U0= 140.001
P0= 89280.81
L0=1.
R0=1.111711
T0=279.15
Model = 'NSTurbulent'
tc_out = C.addState(tc_out, 'GoverningEquations', Model)
tc_out = C.addState(tc_out, UInf=U0, RoInf=R0, PInf=P0, LInf=L0, alphaZ=0., adim='dim3')


# Copie la distribution de t dans tc
if NP > 0: D2._copyDistribution(tc_out,t)
C.convertPyTree2File(tc_out, 'tc_'+str(IT_FIN)+'.cgns')
