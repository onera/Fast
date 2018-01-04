import Converter.PyTree as C
import Converter.Internal as Internal
import Connector.PyTree as X
import Transform.PyTree as T
import Generator.PyTree as G
import numpy
import math

#==============================================================================
# masquage par des surfaces definies dans un arbre tb
# gridType = single or composite - composite means that an off body grid exists
#==============================================================================
def blankByBodies(t, tb, loc, dim, gridType='single'):
    blankalgo = 'tri'
    DIM = dim
    if DIM == 2: blankalgo = 'xray'
    bodies = []
    for b in Internal.getBases(tb):
        wallsl = Internal.getNodesFromType1(b, 'Zone_t')
        if wallsl != []: 
            wallsl = C.convertArray2Tetra(wallsl)
            wallsl = T.join(wallsl)
            wallsl = G.close(wallsl)
            if DIM ==3:
                try: P.exteriorFaces(wallsl) 
                except: bodies.append([wallsl])
            else: bodies.append([wallsl])

    nbodies = len(bodies)
    print 'Blanking mesh by %d bodies'%nbodies
    if loc == 'centers': typeb = 'cell_intersect'
    else: typeb = 'node_in'
    nbases = len(Internal.getBases(t))
    BM = numpy.ones((nbases,nbodies),dtype=numpy.int32)
    if gridType=='composite':
        if nbodies < nbases:
           for nob in xrange(nbodies): BM[nob,nob]=0
    if blankalgo == 'xray' or DIM == 2:
        XRAYDIM1 = 2000; XRAYDIM2 = XRAYDIM1
        if DIM == 2:  XRAYDIM2 = 2
        t = X.blankCells(t, bodies,BM,blankingType=typeb,delta=1.e-10,XRaydim1=XRAYDIM1,XRaydim2=XRAYDIM2,dim=DIM)
    else:
        t = X.blankCellsTri(t,bodies,BM,blankingType=typeb)
    return t

#----------------------------------------------------------------------------------------------------------
# Mark cellN to 0 for NGhostCells layers of cells near window of range win=[imin,imax,jmin,jmax,kmin,kmax]
#----------------------------------------------------------------------------------------------------------
def _blankCellsInRange(z, win, NGhostCells):
    [i1,i2,j1,j2,k1,k2] = win # fenetre en noeuds
    dims = Internal.getZoneDim(z)
    
    ni = dims[1]; nj = dims[2]; nk = dims[3]; ninj = ni*nj
    nic=ni-1; njc=nj-1; nicnjc=(ni-1)*(nj-1)
    cellN = C.getField("centers:cellN",z)[0]
    il1=max(1,i1-1); il2=max(1,i2-1); jl1=max(1,j1-1); jl2=max(1,j2-1); kl1=max(1,k1-1); kl2=max(1,k2-1)
    dir1 = il2-il1; dir2 = jl2-jl1; 
    dir3 = kl2-kl1
    il3 = il1+NGhostCells-1; il4 = il2-NGhostCells
    jl3 = jl1+NGhostCells-1; jl4 = jl2-NGhostCells
    kl3 = kl1+NGhostCells-1; kl4 = kl2-NGhostCells

    if dir1==0: 
        if il1 == 1: il1 = 0; il2 = il1+NGhostCells; il4 = il3+NGhostCells
        else: il1 = il1-NGhostCells; il3 = il1-NGhostCells
        
    elif dir2==0:
        if j1 == 1: jl1 = 0; jl2 = jl1+NGhostCells; jl4 = jl3+NGhostCells
        else: jl1 = jl1-NGhostCells; jl3 = jl1-NGhostCells

    elif dir3==0:
        if k1 == 1: kl1 = 0; kl2 = kl1+NGhostCells; kl4 = kl3+NGhostCells
        else: kl1 = kl1-NGhostCells; kl3 = kl1-NGhostCells
    
    cellNT = cellN[1][0,:]
    if nk == 2:
        for j in xrange(jl1,jl2):
            for i in xrange(il1,il2):
                ind = i + j*nic
                cellNT[ind]=0.

        for j in xrange(jl3,jl4):
            for i in xrange(il3,il4):
                ind = i + j*nic
                cellNT[ind]=2.
    else:
        for k in xrange(kl1,kl2):
            for j in xrange(jl1,jl2):
                for i in xrange(il1,il2):
                    ind = i + j*nic+k*nicnjc
                    cellNT[ind]=0.
        for k in xrange(kl3,kl4):
            for j in xrange(jl3,jl4):
                for i in xrange(il3,il4):
                    ind = i + j*nic+k*nicnjc
                    cellNT[ind]=2.
    C.setFields([cellN],z,loc='centers')
    return None


#------------------------------------------------------------------
# Mark cellN to 0 for the NGhostCells first layers near a BCOverlap
#------------------------------------------------------------------
def modifyBCOverlapsForGhostMesh(t,NGhostCells):
    for z in Internal.getZones(t):
        #parent,noz = Internal.getParentOfNode(t,z)
        GCS = Internal.getNodesFromType(z,'GridConnectivity_t')
        for GC in GCS:
            gctype=Internal.getNodeFromType1(GC,'GridConnectivityType_t')
            PR = Internal.getNodeFromName1(GC,'PointRange')
            if gctype is not None:
                if Internal.getValue(gctype)=='Overset':
                    win = Internal.range2Window(PR[1])
                    _blankCellsInRange(z,win,NGhostCells)
    return t

#-----------------------------------------------------------------------------
# INTERPOLATIONS CHIMERE PERIODIQUE ENTRE LE ROTOR ET LE STATOR
# Les donnees d interpolation sont calculees dans le repere cylindrique
#-----------------------------------------------------------------------------
def setInterpDataRS(tcyl,tc,THETA,DTHETA, IT_DEB, IT_FIN, (XC0,YC0,ZC0), (AXISX,AXISY,AXISZ), check=False):

    basenames  = ['Stator', 'Rotor' ]
    donorBases = {}
    theta_meanDonor = {}
    theta_meanRecept= {}
    theta_perioDonor= {}

    theta_abs  = abs(THETA)
    dtheta_abs = abs(DTHETA)
    for name in basenames:

       base   = Internal.getNodeFromName(tc[name], name)
       zones  = Internal.getZones(base)


       if AXISX > 0. or AXISZ > 0.:
          thetameanD = C.getMeanValue(zones, 'CoordinateY')
          zonesPerioD= T.translate(zones,(0, theta_abs, 0))
          zonesPerioG= T.translate(zones,(0,-theta_abs, 0))
       elif AXISY > 0.: 
          thetameanD = C.getMeanValue(zones, 'CoordinateX')
          zonesPerioD= T.translate(zones,( theta_abs, 0, 0))
          zonesPerioG= T.translate(zones,(-theta_abs, 0, 0))

       for c in range( len(zones)):
          zonesPerioD[c][0]=  zones[c][0]+'_PeriodicD'
          zonesPerioG[c][0]=  zones[c][0]+'_PeriodicG'

       newbase = Internal.newCGNSBase( name             ); newbase[2]  += zones
       newbaseD= Internal.newCGNSBase( name+'_PeriodicD'); newbaseD[2] += zonesPerioD
       newbaseG= Internal.newCGNSBase( name+'_PeriodicG'); newbaseG[2] += zonesPerioG

       tmp = [newbaseG,newbase,newbaseD]
       tmp = Internal.rmNodesByName( tmp,'ID_*')

       theta_perioDonor[name]= [      -theta_abs      , 0.        ,       theta_abs        ]
       theta_meanDonor[name] = [ thetameanD- theta_abs, thetameanD, thetameanD+ theta_abs, ]
       donorBases[name] = tmp
       print name, ' is Donor.  <angle>= ', theta_meanDonor[name][1]/math.pi*180,' <angle> donnorGD=', theta_meanDonor[name][0]/math.pi*180, theta_meanDonor[name][2]/math.pi*180


       #calcul theta_meanRecepteur
       base   = Internal.getNodeFromName(tcyl[name], name)
       zones  = Internal.getZones(base)
       thetameanR= C.getMeanValue(zones, 'CoordinateY')
       theta_meanRecept[name] = [thetameanR] 
       print  name, ' is Recept. <angle>= ', thetameanR/math.pi*180
      
    #rotation parameter
    THETADEG  = THETA/math.pi*180
    RotCenter = numpy.zeros((3), numpy.float64)
    RotCenter[0] = XC0 
    RotCenter[1] = YC0
    RotCenter[2] = ZC0
    
    #tcyl['Rotor']       = C.initVars(tcyl['Rotor']      ,"CoordinateX0={CoordinateX}")
    #tcyl['Rotor']       = C.initVars(tcyl['Rotor']      ,"CoordinateY0={CoordinateY}")
    #donorBases['Rotor'] = C.initVars(donorBases['Rotor'],"CoordinateX0={CoordinateX}")
    #donorBases['Rotor'] = C.initVars(donorBases['Rotor'],"CoordinateY0={CoordinateY}")

    if AXISX > 0. or AXISZ > 0.:
      T._translate( tcyl['Rotor']       ,  (0, (IT_DEB-1)*DTHETA, 0))
      T._translate( donorBases['Rotor'] ,  (0, (IT_DEB-1)*DTHETA, 0))
    elif AXISY > 0.:
      T._translate( tcyl['Rotor']       ,  ( (IT_DEB-1)*DTHETA, 0, 0))
      T._translate( donorBases['Rotor'] ,  ( (IT_DEB-1)*DTHETA, 0, 0))

    ReceptCyl={}
    for it in xrange(IT_DEB, IT_FIN):
        print '-----------'
        print 'theta(radians,degree) = ', it*DTHETA,' ,', it*DTHETA/math.pi*180
        print '-----------'
        RotAngleDG= [ numpy.zeros( (3), numpy.float64 ) , numpy.zeros( (3), numpy.float64), numpy.zeros( (3), numpy.float64) ]

        for rcpt in basenames:
        #for rcpt in ['Rotor']:
           print '  '

           ReceptCyl[ rcpt ] = Internal.getNodeFromName(tcyl[ rcpt ], rcpt )

           if rcpt == 'Rotor':
                zones2translateR =  ReceptCyl['Rotor']
                donor            = 'Stator'
                #angle moyen receveur
                theta_meanR= theta_meanRecept[rcpt][0]+it*DTHETA
                #angle moyen donneur
                theta_meanDonorD= theta_meanDonor[donor][2]
                theta_meanDonorG= theta_meanDonor[donor][0]
                print 'Rcpt  rotor.  <angle> receveur=', theta_meanR*180/math.pi, '<angle> donnorGD=',  theta_meanDonorG*180/math.pi, theta_meanDonorD*180/math.pi
           else:
                zones2translateR =  donorBases['Rotor']
                donor            = 'Rotor'
                #angle moyen receveur
                theta_meanR= theta_meanRecept[rcpt][0]
                #angle moyen donneur
                theta_meanDonorD= theta_meanDonor[donor][2] +it*DTHETA
                theta_meanDonorG= theta_meanDonor[donor][0] +it*DTHETA
                print 'Rcpt stator. <angle> receveur=', theta_meanR*180/math.pi, '<angle> donnorGD=',  theta_meanDonorG*180/math.pi, theta_meanDonorD*180/math.pi

           Internal._rmNodesByName( donorBases[donor],'ID_*')

           if AXISX > 0. or AXISZ > 0.:

                #C._initVars( zones2translateR ,'CoordinateY={CoordinateY0}')
                #T._translate(zones2translateR ,(0, it*DTHETA, 0))
                T._translate(zones2translateR ,(0, DTHETA, 0))

                coef_secure = 1.3
                if    theta_meanR +theta_abs*coef_secure <  theta_meanDonorD:

                      T._translate(donorBases[donor][2],(0, -3*theta_abs, 0))

                      #mise a jour base Gauche, centre, droite
                      tmp1 = donorBases[donor][0]
                      donorBases[donor][0] =   donorBases[donor][2]
                      donorBases[donor][2] =   donorBases[donor][1]
                      donorBases[donor][1] =   tmp1
                      #mise a jour angle moyen base Gauche, centre, droite
                      tmp2 = theta_meanDonor[donor][0]
                      theta_meanDonor[donor][0] =   theta_meanDonor[donor][2] -3*theta_abs
                      theta_meanDonor[donor][2] =   theta_meanDonor[donor][1]
                      theta_meanDonor[donor][1] =   tmp2
                     
                      #mise a jour peridicite base Gauche, centre, droite
                      tmp3 = theta_perioDonor[donor][0]
                      theta_perioDonor[donor][0] =   -2*theta_abs
                      theta_perioDonor[donor][2] =   theta_perioDonor[donor][1]
                      theta_perioDonor[donor][1] =   tmp3

                      print 'Donnor_D du ', rcpt, ' moved by', -3*theta_abs*180/3.1415,'. <Angle> receveur + perio=', (theta_meanR+theta_abs)*180/math.pi,'. <Angle> donor D=',theta_meanDonorD*180/math.pi

                elif theta_meanR -theta_abs*coef_secure >  theta_meanDonorG:

                      T._translate(donorBases[donor][0],(0, 3*theta_abs, 0))

                      #mise a jour base Gauche, centre, droite
                      tmp1 = donorBases[donor][2]
                      donorBases[donor][2] =   donorBases[donor][0]
                      donorBases[donor][0] =   donorBases[donor][1]
                      donorBases[donor][1] =   tmp1

                      #mise a jour angle moyen base Gauche, centre, droite
                      tmp2 = theta_meanDonor[donor][2]
                      theta_meanDonor[donor][2] =   theta_meanDonor[donor][0]
                      theta_meanDonor[donor][0] =   theta_meanDonor[donor][1]
                      theta_meanDonor[donor][1] =   tmp2

                      #mise a jour peridicite base Gauche, centre, droite
                      tmp3 = theta_perioDonor[donor][2]
                      theta_perioDonor[donor][2] =   2*theta_abs
                      theta_perioDonor[donor][0] =   theta_perioDonor[donor][1]
                      theta_perioDonor[donor][1] =   tmp3

                      print 'Donnor_G du ', rcpt, ' moved by', 3*theta_abs*180/3.1415,'. <Angle> receveur + perio=', (theta_meanR+theta_abs)*180/math.pi,'. <Angle> donor G=',theta_meanDonorG*180/math.pi
                      
           elif AXISY > 0.: 

                #C._initVars(  zones2translateR , 'CoordinateX={CoordinateX0}')
                #T._translate( zones2translateR , (it*DTHETA, 0, 0))
                T._translate(zones2translateR ,( DTHETA, 0, 0))

                coef_secure = 1.3

                if    theta_meanR +theta_abs*coef_secure <  theta_meanDonorD:

                      T._translate(donorBases[donor][2],( -3*theta_abs, 0, 0))

                      #mise a jour base Gauche, centre, droite
                      tmp1 = donorBases[donor][0]
                      donorBases[donor][0] =   donorBases[donor][2]
                      donorBases[donor][2] =   donorBases[donor][1]
                      donorBases[donor][1] =   tmp1
                      #mise a jour angle moyen base Gauche, centre, droite
                      tmp2 = theta_meanDonor[donor][0]
                      theta_meanDonor[donor][0] =   theta_meanDonor[donor][2] -3*theta_abs
                      theta_meanDonor[donor][2] =   theta_meanDonor[donor][1]
                      theta_meanDonor[donor][1] =   tmp2
                     
                      #mise a jour peridicite base Gauche, centre, droite
                      tmp3 = theta_perioDonor[donor][0]
                      theta_perioDonor[donor][0] =   -2*theta_abs
                      theta_perioDonor[donor][2] =   theta_perioDonor[donor][1]
                      theta_perioDonor[donor][1] =   tmp3

                      print 'Donnor_D du ', rcpt, ' moved by', -3*theta_abs*180/3.1415,'. <Angle> receveur + perio=', (theta_meanR+theta_abs)*180/math.pi,'. <Angle> donor D=',theta_meanDonorD*180/math.pi

                elif theta_meanR -theta_abs*coef_secure >  theta_meanDonorG:

                      T._translate(donorBases[donor][0],( 3*theta_abs, 0, 0))

                      #mise a jour base Gauche, centre, droite
                      tmp1 = donorBases[donor][2]
                      donorBases[donor][2] =   donorBases[donor][0]
                      donorBases[donor][0] =   donorBases[donor][1]
                      donorBases[donor][1] =   tmp1

                      #mise a jour angle moyen base Gauche, centre, droite
                      tmp2 = theta_meanDonor[donor][2]
                      theta_meanDonor[donor][2] =   theta_meanDonor[donor][0]
                      theta_meanDonor[donor][0] =   theta_meanDonor[donor][1]
                      theta_meanDonor[donor][1] =   tmp2

                      #mise a jour peridicite base Gauche, centre, droite
                      tmp3 = theta_perioDonor[donor][2]
                      theta_perioDonor[donor][2] =   2*theta_abs
                      theta_perioDonor[donor][0] =   theta_perioDonor[donor][1]
                      theta_perioDonor[donor][1] =   tmp3

                      print 'Donnor_G du ', rcpt, ' moved by', 3*theta_abs*180/3.1415,'. <Angle> receveur + perio=', (theta_meanR+theta_abs)*180/math.pi,'. <Angle> donor G=',theta_meanDonorG*180/math.pi

           RD = [ ReceptCyl[rcpt] , donorBases[donor]] 

           #C._rmVars( RD[0], 'FlowSolution')
           #C._rmVars( RD[1], 'FlowSolution')
           #C.convertPyTree2File( RD[0],rcpt+'recept_'+str(it)+'.cgns')
           #C.convertPyTree2File( RD[1],rcpt+'donors_'+str(it)+'.cgns')

           RD[1] = X.setInterpData( RD[0], RD[1], storage='inverse',loc='centers',penalty=1,nature=1,itype='chimera')
        
           # bloc donneur periodiques -> remise dans le bloc original
           # bloc original     =RD[1][0] 
           # bloc perio Droite =RD[1][1] 
           # bloc perio Gauche =RD[1][2]
           #bloc_target= RD[1][0]  #bloc original
           c = 0
           for perio in theta_perioDonor[donor]:
             if perio == 0. : blk_target = c 
             c+=1

           bloc_target= RD[1][blk_target]  #bloc original
           c = 0
           for bloc_Perio in RD[1][0:]:
          

             if c != blk_target:
               RotAngleDG[c][0] =AXISX*theta_perioDonor[donor][c]/math.pi*180
               RotAngleDG[c][1] =AXISY*theta_perioDonor[donor][c]/math.pi*180
               RotAngleDG[c][2] =AXISZ*theta_perioDonor[donor][c]/math.pi*180



               for nozd in xrange(len(bloc_Perio[2])):
                 zdperio = bloc_Perio[2][nozd]
                 if zdperio[3] == "Zone_t":
                    zname  = zdperio[0].split('_Perio')[0]
                    zdorig = Internal.getNodeFromName(tc[donor],zname)

                    for zsr in Internal.getNodesFromName(zdperio,"ID_*"):
                      print rcpt,' est le receveur.  teta_perio=', theta_perioDonor[donor][c], 'present sur bloc', c

                      srname = zsr[0].split('_')
                      srname = srname[1]
                      zsr[0] = 'IDPER#%d_%s'%(it,srname)
                      # ajout des infos de periodicite 
                      Internal.createChild(zsr,'RotationAngle' ,'DataArray_t',value=RotAngleDG[c])
                      Internal.createChild(zsr,'RotationCenter','DataArray_t',value=RotCenter)   
                      zdorig[2].append(zsr)

                    zd = bloc_target[2][nozd]
                    for zsr in Internal.getNodesFromName(zd,"ID_*"):                
                      srname = zsr[0].split('_'); srname = srname[1]
                      zsr[0] = 'ID#%d_%s'%(it,srname)
                      zdorig[2].append(zsr)
             c+=1

        #--------------CHECKS----------------------------
        if check:
            #baseRotorCylc = C.node2Center(baseRotorCyl)
            C.convertPyTree2File(ReceptCyl['Rotor'],"rotor.cgns")
            C.convertPyTree2File(donorBases['Stator'],"stator.cgns")
    
        Internal._rmNodesByName(tc[rcpt],'GridCoordinates')

    connec_tree=[]
    for name in basenames:  connec_tree.append( tc[name] )

    return connec_tree
