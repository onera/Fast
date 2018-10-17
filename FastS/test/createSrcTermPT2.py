import Converter.PyTree as C
import Generator.PyTree as G
import Initiator.PyTree as I
import CPlot.PyTree as CPlot
import Fast.PyTree as Fast
import FastS.PyTree as FastS
import Converter.Internal as Internal

mach = 0.7
a = G.cart((0,0,0), (0.5,0.5,0.25), (200,100,2))
I._initLamb(a, position=(25.,25.), Gamma=2., MInf=mach, loc='centers')
t = C.newPyTree(['Base', a])
C._addState(t, 'GoverningEquations', 'Euler')
C._addState(t, MInf=mach)

# Create source term
coefa = 2.71128
coefb = 2.4
x0    = 65.
y0    = 30.
z0    = 0.
amp   = 0.01
per   = 0.014
phi   = 0.

alpha  = log(coefa)/(coefb*coefb)
omega  = 2.*pi/per
omt    = omega*temps+phi
eps    = amp*sin(omt)
seuil = exp(-alpha*30.)


C._initVars(t,'{centers:Density_src}=0.')
C._initVars(t,'{centers:MomentumX_src}=0.')
C._initVars(t,'{centers:MomentumY_src}=0.')
C._initVars(t,'{centers:MomentumZ_src}=0.')
C._initVars(t,'{centers:EnergyStagnationDensity_src}=0.')

# Numerics
numb = {}
numb["temporal_scheme"]    = "explicit"
numb["ss_iteration"]       = 10
numb["modulo_verif"]       = 5
numz = {}
numz["time_step"]          = 0.05
numz["scheme"]             = "roe"
numz["source"]		   = 1
Fast._setNum2Zones(t, numz); Fast._setNum2Base(t, numb)

(t, tc, metrics) = FastS.warmup(t, None)

zones      = Internal.getZones(t)

nit = 100 ; time = 0.
timeStep = numz['time_step']
for it in xrange(nit):
    FastS._compute(t, metrics, it)
    if (it%1 == 0):
        print '- %d - %g'%(it, time)
        CPlot.display(t, dim=2, mode=3, isoEdges=1)
    time += timeStep

Fast.save(t,'restart.cgns')

