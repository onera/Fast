# - FastS: Cylindre IBC Octree -

import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Geom.PyTree as D
import Post.PyTree as P
import numpy
import Dist2Walls.PyTree as DTW
import Transform.PyTree as T
import Initiator.PyTree as I
import Fast.PyTree as Fast
import FastS.PyTree as FastS
import CPlot.PyTree as CPlot
import Converter.Internal as Internal
import KCore.Adim as Adim
import sys

restart = False

if restart == False:
    MInf = 0.5 ; alpha = 0.
    adim = Adim.adim1(MInf=MInf)

    # Solver settings
    numb = {'temporal_scheme':'explicit', 'ss_iteration':20} 
    numz = {'time_step':0.01, 'scheme':'ausmpred'}

    # Cylindre
    s = G.cylinder((0,0,0), 1., 2.0, 360., 0., 5*0.01, (160,1,1))

    # Grilles cartesiennes
    NI = 501 ; dh = 40./(NI-1); vmin = 31
    snears = [dh*(vmin-1)]
    o = G.octree([s], snears,dfar = 100., balancing=1)
    C.convertPyTree2File([o], "octree.cgns")

    res = G.octree2Struct(o,vmin=vmin,ext=3,optimized=0,merged=1)
    res = C.fillEmptyBCWith(res, 'far', 'BCFarfield', dim=2) 
    res = T.addkplane(res)

    t = C.newPyTree(['Cart']); t[2][1][2] = res
    t = C.fillMissingVariables(t)

    # Interpolations pour Octree
    t = X.applyBCOverlaps(t, depth=2)
    tc = C.node2Center(t)
    tc = X.setInterpData(t, tc, nature=1, loc='centers', storage='inverse', 
                         sameName=1, method='leastsquares',dim=2)
    t = C.initVars(t, 'centers:cellN', 1.) # init pour les IBCs
    # Blanking
    s2 = T.addkplane(s)
    bodies = [[s2]]
    BM = numpy.array([[1]],numpy.int32)
    t = X.blankCells(t, bodies, BM, blankingType='center_in',dim=2)
    t = X.setHoleInterpolatedPoints(t, depth=-1)
    t = X.setHoleInterpolatedPoints(t, depth=+2)
    tc = C.cpVars(t, 'centers:cellN', tc, 'cellN')

    # Dist2Walls
    sc = C.node2Center(s2)
    C.convertPyTree2File(sc,"body.cgns")

    t = DTW.distance2Walls(t, [sc], type='ortho', loc='centers', signed=1,dim=2)
    t = C.center2Node(t, 'centers:TurbulentDistance')
    # Gradient de distance localise en centres => normales
    t = P.computeGrad(t, 'TurbulentDistance')
    tc = X.setIBCData(t, tc, loc='centers', nature=1, storage='inverse', hi=dh, he=dh,
                      method='leastsquares',dim=2)
    t = C.rmVars(t,['centers:gradxTurbulentDistance','centers:gradyTurbulentDistance','centers:gradzTurbulentDistance',\
                        'TurbulentDistance','centers:TurbulentDistance'])
    tc = C.rmVars(tc, 'cellN') # tres important pour l'instant

    # Init
    t = C.addState(t, 'GoverningEquations', 'Euler')
    t = C.addState(t, MInf=MInf, alphaZ=alpha)
    t = I.initConst(t, MInf=MInf, alphaZ=alpha, loc='centers')
    Fast._setNum2Zones(t, numz) ; Fast._setNum2Base(t, numb)
    #C.convertPyTree2File(t, 'out.cgns') ; sys.exit()
    C.convertPyTree2File(t ,"mesh.cgns")
    C.convertPyTree2File(tc, "tc.cgns")
else:
    t = C.convertFile2PyTree("mesh.cgns")
    tc = C.convertFile2PyTree("tc.cgns")

(t, tc, metrics) = FastS.warmup(t, tc)

nit = 1001
for it in xrange(nit):
    FastS._compute(t, metrics, nit, tc)
    if it%100 == 0:
        print('- %d -'%it); sys.stdout.flush()
        CPlot.display(t, dim=2, mode=3, scalarField=1)
C.convertPyTree2File(t, "out.cgns")
