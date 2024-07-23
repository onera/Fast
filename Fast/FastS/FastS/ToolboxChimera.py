import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.GhostCells as Ghost
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
    print('Blanking mesh by %d bodies'%nbodies)
    if loc == 'centers': typeb = 'cell_intersect'
    else: typeb = 'node_in'
    nbases = len(Internal.getBases(t))
    BM = numpy.ones((nbases,nbodies), dtype=Internal.E_NpyInt)
    if gridType=='composite':
        if nbodies < nbases:
           for nob in range(nbodies): BM[nob,nob]=0
    if blankalgo == 'xray' or DIM == 2:
        XRAYDIM1 = 2000; XRAYDIM2 = XRAYDIM1
        if DIM == 2:  XRAYDIM2 = 2
        t = X.blankCells(t, bodies, BM,blankingType=typeb, delta=1.e-10, XRaydim1=XRAYDIM1, XRaydim2=XRAYDIM2, dim=DIM)
    else:
        t = X.blankCellsTri(t, bodies, BM, blankingType=typeb)
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
        for j in range(jl1,jl2):
            for i in range(il1,il2):
                ind = i + j*nic
                cellNT[ind]=0.

        for j in range(jl3,jl4):
            for i in range(il3,il4):
                ind = i + j*nic
                cellNT[ind]=2.
    else:
        for k in range(kl1,kl2):
            for j in range(jl1,jl2):
                for i in range(il1,il2):
                    ind = i + j*nic+k*nicnjc
                    cellNT[ind]=0.
        for k in range(kl3,kl4):
            for j in range(jl3,jl4):
                for i in range(il3,il4):
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
def setInterpDataRS(tcyl,tc,THETA,DTHETA, IT_DEB, IT_FIN, infos_PtlistRebuild, tCoord0, tAxis, etage_names=None , check=False):

    import timeit

    (XC0,YC0,ZC0) = tCoord0
    (AXISX,AXISY,AXISZ) = tAxis

    if etage_names is None :
        basenames  = ['Stator', 'Rotor' ]
    else :
        basenames = etage_names


    baseRole = []
    name = basenames[0] 
    if int(name[-1])%2 == 0:
             baseRole=['Rotor','Stator']
    else:
             baseRole=['Stator','Rotor']
    #print ('baseRole',baseRole)

    donorBases = {}
    theta_meanDonor = {}
    theta_meanRecept= {}
    theta_perioDonor= {}

    theta_abs  = abs(THETA)
    dtheta_abs = abs(DTHETA)

   
    print('tAxis= ', tAxis)
    print('theta_abs= ', theta_abs)

    bbox={}
    for name in basenames:
       
       base   = Internal.getNodeFromName(tc[name], name)
      
       zones  = Internal.getZones(base)
       for z in zones:
         bbox[ z[0] ] = G.bbox(z)

       
       if AXISX > 0. or AXISZ > 0.:
          print("NAME=",name)
          thetameanD = C.getMeanValue(zones, 'CoordinateZ')
          zonesPerioD= T.translate(zones,(0,0, theta_abs))
          zonesPerioG= T.translate(zones,(0,0,-theta_abs))
          for c in range( len(zones)):
            zonesPerioD[c][0]=  zones[c][0]+'_PeriodicD'
            zonesPerioG[c][0]=  zones[c][0]+'_PeriodicG'
            tmp=[]
            for i in  bbox[ zones[c][0] ]: 
               tmp.append(i)
            bbox[ zonesPerioD[c][0] ] = tmp 
            bbox[ zonesPerioD[c][0] ][2]= bbox[ zonesPerioD[c][0] ][2] + theta_abs
            bbox[ zonesPerioD[c][0] ][5]= bbox[ zonesPerioD[c][0] ][5] + theta_abs

            tmp1=[] 
            for i in  bbox[ zones[c][0] ]: 
              tmp1.append(i)
            bbox[ zonesPerioG[c][0] ] = tmp1  
            bbox[ zonesPerioG[c][0] ][2]= bbox[ zonesPerioG[c][0] ][2] - theta_abs
            bbox[ zonesPerioG[c][0] ][5]= bbox[ zonesPerioG[c][0] ][5] - theta_abs


       elif AXISY > 0.: 
          thetameanD = C.getMeanValue(zones, 'CoordinateX')
          zonesPerioD= T.translate(zones,( theta_abs, 0, 0))
          zonesPerioG= T.translate(zones,(-theta_abs, 0, 0))
          for c in range( len(zones)):
            zonesPerioD[c][0]=  zones[c][0]+'_PeriodicD'
            zonesPerioG[c][0]=  zones[c][0]+'_PeriodicG'
            tmp=[]
            for i in  bbox[ zones[c][0] ]: 
               tmp.append(i)
            bbox[ zonesPerioD[c][0] ] = tmp 
            bbox[ zonesPerioD[c][0] ][0]= bbox[ zonesPerioD[c][0] ][0] + theta_abs
            bbox[ zonesPerioD[c][0] ][3]= bbox[ zonesPerioD[c][0] ][3] + theta_abs
            tmp1=[] 
            for i in  bbox[ zones[c][0] ]: 
              tmp1.append(i)
            bbox[ zonesPerioG[c][0] ] = tmp1  
            bbox[ zonesPerioG[c][0] ][0]= bbox[ zonesPerioG[c][0] ][0] - theta_abs
            bbox[ zonesPerioG[c][0] ][3]= bbox[ zonesPerioG[c][0] ][3] - theta_abs

       newbase = Internal.newCGNSBase( name             ); newbase[2]  += zones
       newbaseD= Internal.newCGNSBase( name+'_PeriodicD'); newbaseD[2] += zonesPerioD
       newbaseG= Internal.newCGNSBase( name+'_PeriodicG'); newbaseG[2] += zonesPerioG

       tmp = [newbaseG,newbase,newbaseD]
       tmp = Internal.rmNodesByName( tmp,'ID_*')



       theta_perioDonor[name]= [      -theta_abs      , 0.        ,       theta_abs        ]
       theta_meanDonor[name] = [ thetameanD- theta_abs, thetameanD, thetameanD+ theta_abs, ]
       donorBases[name] = tmp

       #C.convertPyTree2File(donorBases[name],name+'donor_.cgns')

       
       print(name, ' is Donor.  <angle>= ', theta_meanDonor[name][1]/math.pi*180,' <angle> donorGD=', theta_meanDonor[name][0]/math.pi*180, theta_meanDonor[name][2]/math.pi*180)


       #calcul theta_meanRecepteur
       base   = Internal.getNodeFromName(tcyl[name], name)
       zones  = Internal.getZones(base)
       thetameanR= C.getMeanValue(zones, 'CoordinateZ')
       theta_meanRecept[name] = [thetameanR] 
       print(name, ' is Recept. <angle>= ', thetameanR/math.pi*180)
      
    #rotation parameter
    THETADEG  = THETA/math.pi*180
    RotCenter = numpy.zeros((3), numpy.float64)
    RotCenter[0] = XC0 
    RotCenter[1] = YC0
    RotCenter[2] = ZC0
    
    for name in basenames:
        i = basenames.index(name)
        if (baseRole[i]=='Rotor'): 
            if AXISX > 0. or AXISZ > 0.:
                T._translate( tcyl[etage_names[i]]       ,  (0, 0, (IT_DEB-1)*DTHETA))
                T._translate( donorBases[etage_names[i]] ,  (0, 0, (IT_DEB-1)*DTHETA))  
                for z in Internal.getZones(donorBases[etage_names[i]]):
                    bbox[ z[0] ][2]= bbox[ z[0] ][2] + (IT_DEB-1)*DTHETA
                    bbox[ z[0] ][5]= bbox[ z[0] ][5] + (IT_DEB-1)*DTHETA

            elif AXISY > 0.:
                T._translate( tcyl[etage_names[i]]       ,  ( (IT_DEB-1)*DTHETA, 0, 0))
                T._translate( donorBases[etage_names[i]] ,  ( (IT_DEB-1)*DTHETA, 0, 0))
                for z in Internal.getZones(donorBases[etage_names[i]]):
                    bbox[ z[0] ][0]= bbox[ z[0] ][0] + (IT_DEB-1)*DTHETA
                    bbox[ z[0] ][3]= bbox[ z[0] ][3] + (IT_DEB-1)*DTHETA


    bboxR={}
    for key in bbox:
      bboxR[key]= bbox[key][:]
    #Hook optimisation recherche

    coef_secure = 1.3
    ReceptCyl={}
    for rcpt in basenames:
        ReceptCyl[ rcpt ] = Internal.getNodeFromName(tcyl[ rcpt ], rcpt )

    for it in range(IT_DEB, IT_FIN):
        print('-----------')
        print('theta(radians,degree,it) = ', it*DTHETA,' ,', it*DTHETA/math.pi*180,' ,', it)
        print('-----------')

        for rcpt in basenames:
           print('  ')
 
           t0=timeit.default_timer()

           RotAngleDG= [ numpy.zeros( (3), numpy.float64 ) , numpy.zeros( (3), numpy.float64), numpy.zeros( (3), numpy.float64) ]

           #ReceptCyl[ rcpt ] = Internal.getNodeFromName(tcyl[ rcpt ], rcpt )

           #if rcpt == 'Rotor':
           i = basenames.index(rcpt)
           if (baseRole[i]=='Rotor'): 
              
                #zones2translateR =  ReceptCyl['Rotor']
                #donor            = 'Stator'
                zones2translateR =  ReceptCyl[etage_names[i]]
                donor            = etage_names[baseRole.index('Stator')]
                choice_tree      = 0
                #angle moyen receveur
                theta_meanR= theta_meanRecept[rcpt][0]+it*DTHETA
                #angle moyen donneur
                theta_meanDonorD= theta_meanDonor[donor][2]
                theta_meanDonorG= theta_meanDonor[donor][0]
                print('Rcpt  rotor.  <angle> receveur=', theta_meanR*180/math.pi, '<angle> donorGD=',  theta_meanDonorG*180/math.pi, theta_meanDonorD*180/math.pi)
           else:
                #zones2translateR =  donorBases['Rotor']
                #donor            = 'Rotor'
                zones2translateR =  donorBases[etage_names[baseRole.index('Rotor')]]
                donor            = etage_names[baseRole.index('Rotor')]
                choice_tree      = 1
                #angle moyen receveur
                theta_meanR= theta_meanRecept[rcpt][0]
                #angle moyen donneur
                theta_meanDonorD= theta_meanDonor[donor][2] +it*DTHETA
                theta_meanDonorG= theta_meanDonor[donor][0] +it*DTHETA
                print('Rcpt stator. <angle> receveur=', theta_meanR*180/math.pi, '<angle> donorGD=',  theta_meanDonorG*180/math.pi, theta_meanDonorD*180/math.pi)

           Internal._rmNodesByName( donorBases[donor],'ID_*')


           if AXISX > 0. or AXISZ > 0.:

                T._translate(zones2translateR ,(0, 0, DTHETA))
                #mjr bbox
                for z in Internal.getZones(zones2translateR):
                    if choice_tree==0:
                      bboxR[ z[0] ][2]= bboxR[ z[0] ][2] + DTHETA 
                      bboxR[ z[0] ][5]= bboxR[ z[0] ][5] + DTHETA 
                    else:
                      bbox[ z[0] ][2]= bbox[ z[0] ][2] + DTHETA 
                      bbox[ z[0] ][5]= bbox[ z[0] ][5] + DTHETA 

                print('theta_meanR= ',theta_meanR*180/math.pi , 'theta_meanDonorD= ', theta_meanDonorD*180/math.pi)

                if    theta_meanR +theta_abs*coef_secure <  theta_meanDonorD:

                      T._translate(donorBases[donor][2],(0,0,-3*theta_abs))
                      for z in Internal.getZones(donorBases[donor][2]):
                         bbox[ z[0] ][2]= bbox[ z[0] ][2] -3*theta_abs
                         bbox[ z[0] ][5]= bbox[ z[0] ][5] -3*theta_abs
                    
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


                      
                      print('Donor_D du ', rcpt, ' moved by', -3*theta_abs*180/3.1415,'. <Angle> receveur + perio=', (theta_meanR+theta_abs)*180/math.pi,'. <Angle> donor D=',theta_meanDonor[donor][2]*180/math.pi)

                elif theta_meanR -theta_abs*coef_secure >  theta_meanDonorG:

                      T._translate(donorBases[donor][0],(0,0,3*theta_abs))
                      for z in Internal.getZones(donorBases[donor][0]):
                         bbox[ z[0] ][2]= bbox[ z[0] ][2] +3*theta_abs
                         bbox[ z[0] ][5]= bbox[ z[0] ][5] +3*theta_abs

                      #mise a jour base Gauche, centre, droite
                      tmp1 = donorBases[donor][2]
                      donorBases[donor][2] =   donorBases[donor][0]
                      donorBases[donor][0] =   donorBases[donor][1]
                      donorBases[donor][1] =   tmp1

                      #mise a jour angle moyen base Gauche, centre, droite
                      tmp2 = theta_meanDonor[donor][2]
                      theta_meanDonor[donor][2] =   theta_meanDonor[donor][0]+3*theta_abs #!!!
                      theta_meanDonor[donor][0] =   theta_meanDonor[donor][1]
                      theta_meanDonor[donor][1] =   tmp2

                      #mise a jour peridicite base Gauche, centre, droite
                      tmp3 = theta_perioDonor[donor][2]
                      theta_perioDonor[donor][2] =   2*theta_abs
                      theta_perioDonor[donor][0] =   theta_perioDonor[donor][1]
                      theta_perioDonor[donor][1] =   tmp3


                      print('Donor_G du ', rcpt, ' moved by', 3*theta_abs*180/3.1415,'. <Angle> receveur + perio=', (theta_meanR+theta_abs)*180/math.pi,'. <Angle> donor G=',theta_meanDonor[donor][0]*180/math.pi)
                      
           elif AXISY > 0.: 

                #C._initVars(  zones2translateR , 'CoordinateX={CoordinateX0}')
                #T._translate( zones2translateR , (it*DTHETA, 0, 0))
                T._translate(zones2translateR ,( DTHETA, 0, 0))
                for z in Internal.getZones(zones2translateR):
                    if choice_tree==0:
                      bboxR[ z[0] ][0]= bboxR[ z[0] ][0] + DTHETA
                      bboxR[ z[0] ][3]= bboxR[ z[0] ][3] + DTHETA
                    else:
                      bbox[ z[0] ][0]= bbox[ z[0] ][0] + DTHETA
                      bbox[ z[0] ][3]= bbox[ z[0] ][3] + DTHETA

                if    theta_meanR +theta_abs*coef_secure <  theta_meanDonorD:

                      T._translate(donorBases[donor][2],( -3*theta_abs, 0, 0))
                      for z in Internal.getZones(donorBases[donor][2]):
                         bbox[ z[0] ][0]= bbox[ z[0] ][0] -3*theta_abs
                         bbox[ z[0] ][3]= bbox[ z[0] ][3] -3*theta_abs

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

                      print('Donor_D du ', rcpt, ' moved by', -3*theta_abs*180/3.1415,'. <Angle> receveur + perio=', (theta_meanR+theta_abs)*180/math.pi,'. <Angle> donor D=',theta_meanDonorD*180/math.pi)

                elif theta_meanR -theta_abs*coef_secure >  theta_meanDonorG:

                      T._translate(donorBases[donor][0],( 3*theta_abs, 0, 0))
                      for z in Internal.getZones(donorBases[donor][0]):
                         bbox[ z[0] ][0]= bbox[ z[0] ][0] +3*theta_abs
                         bbox[ z[0] ][3]= bbox[ z[0] ][3] +3*theta_abs

                      #mise a jour base Gauche, centre, droite
                      tmp1 = donorBases[donor][2]
                      donorBases[donor][2] =   donorBases[donor][0]
                      donorBases[donor][0] =   donorBases[donor][1]
                      donorBases[donor][1] =   tmp1

                      #mise a jour angle moyen base Gauche, centre, droite
                      tmp2 = theta_meanDonor[donor][2]
                      theta_meanDonor[donor][2] =   theta_meanDonor[donor][0] +3*theta_abs
                      theta_meanDonor[donor][0] =   theta_meanDonor[donor][1]
                      theta_meanDonor[donor][1] =   tmp2

                      #mise a jour peridicite base Gauche, centre, droite
                      tmp3 = theta_perioDonor[donor][2]
                      theta_perioDonor[donor][2] =   2*theta_abs
                      theta_perioDonor[donor][0] =   theta_perioDonor[donor][1]
                      theta_perioDonor[donor][1] =   tmp3

                      print('Donor_G du ', rcpt, ' moved by', 3*theta_abs*180/3.1415,'. <Angle> receveur + perio=', (theta_meanR+theta_abs)*180/math.pi,'. <Angle> donor G=',theta_meanDonorG*180/math.pi)

           t1=timeit.default_timer()
           print( "cout rotation= ", t1-t0, 'it= ', it,rcpt)
           t0=timeit.default_timer()

           #tentative optim
           #OPTIM_INTERP = False
           OPTIM_INTERP = True
           if OPTIM_INTERP == True:
              baseR = Internal.getBases( ReceptCyl[rcpt] )
              tcout = Internal.copyRef(donorBases[donor] )
              baseD = Internal.getBases( tcout )

              namebaseD=[]
              for base in baseD:
                namebaseD.append(base[0])

              for zr in Internal.getZones(ReceptCyl[rcpt]):

                 rminR= bboxR[zr[0]][1]; rmaxR= bboxR[zr[0]][4]; tminR= bboxR[zr[0]][2]; tmaxR= bboxR[zr[0]][5]
          
                 #print('receuse', zr[0],  rminR , bboxR[ zr[0] ][1],  rmaxR, bboxR[ zr[0] ][4], 'tmin=', tminR, bboxR[ zr[0] ][2], 'tmax=',tmaxR, bboxR[ zr[0] ][5], DTHETA)

                 donor_filtre  = C.newPyTree(namebaseD)
                 cpt=1
                 for base in baseD:
                   zDout=[] 
                   for zd in Internal.getZones(base ):
                      rminD= bbox[ zd[0] ][1]; rmaxD= bbox[ zd[0] ][4]; tminD= bbox[ zd[0] ][2]; tmaxD= bbox[ zd[0] ][5]
                      #bbox[ zd[0] ] = G.bbox(zd)
                      #print('doneuseCC', zd[0],  rminD,rmaxD, tminD, tmaxD)
                      if rmaxD >= rminR and rminD <= rmaxR and tmaxD >= tminR and tminD <= tmaxR:
                         #print('doneuse', zd[0],  rminD - bbox[ zd[0] ][1],  rmaxD- bbox[ zd[0] ][4], tminD- bbox[ zd[0] ][2], tmaxD- bbox[ zd[0] ][5])
                         zDout.append(zd)
                      
                   donor_filtre[2][cpt][2]+= zDout
                   cpt +=1

                 RD = [  zr,  donor_filtre]
                 #C.convertPyTree2File( RD[1] ,'RD1'+rcpt+str(it)+zr[0]+'.cgns')
                 #C.convertPyTree2File( RD[0] ,'RD0'+rcpt+str(it)+zr[0]+'.cgns')
  
                 RD[1] = X.setInterpData( RD[0], RD[1], storage='inverse',loc='centers',penalty=1,nature=1,itype='chimera')

                 for z in Internal.getZones(RD[1]):
                   racs = Internal.getNodesFromName(z,'ID_*')
                   ztg  = Internal.getNodeFromName(donorBases[donor], z[0])
                   for rac in racs:
                      check  = Internal.getNodeFromName(ztg, rac[0])
                      if check is not None: rac[0] = rac[0]+"bis"
                      ztg[2].append(rac)

           RD = [ ReceptCyl[rcpt] , donorBases[donor]]

           if OPTIM_INTERP == False:
              RD[1] = X.setInterpData( RD[0], RD[1], storage='inverse',loc='centers',penalty=1,nature=1,itype='chimera')

           t1=timeit.default_timer()
           print( "cout interp= ", t1-t0, 'it= ', it, rcpt)
           t0=timeit.default_timer()

           #C.convertPyTree2File( RD[1],'donors'+str(it)+'.cgns')

           ##C._rmVars( RD[0], 'FlowSolution')
           ##C._rmVars( RD[1], 'FlowSolution')

           #if it == 10 :
           #C.convertPyTree2File( RD[0],'recept'+str(it)+'.cgns')
           #C.convertPyTree2File( RD[1],'donors'+str(it)+'.cgns')

           # bloc perio Droite =RD[1][2] 
           # bloc original     =RD[1][1] 
           # bloc perio Gauche =RD[1][0]
           # bloc donneur periodiques -> remise dans le bloc original +  ajout des infos de periodicite 
           c = 0
           for bloc_Perio in RD[1][0:]:
          
             RotAngleDG[c][0] =AXISX*theta_perioDonor[donor][c]
             RotAngleDG[c][1] =AXISY*theta_perioDonor[donor][c]
             RotAngleDG[c][2] =AXISZ*theta_perioDonor[donor][c]

             if   c==0: source='M'
             elif c==1: source='C'
             else     : source='P'

             for nozd in range(len(bloc_Perio[2])):
               zdperio = bloc_Perio[2][nozd]
               if zdperio[3] == "Zone_t":
                  zname  = zdperio[0].split('_Perio')[0]
                  zdorig = Internal.getNodeFromName1(tc[donor],zname)
                  #for zsr in Internal.getNodesFromName(zdperio,"ID_*"):
                  for zsr in Internal.getNodesFromType1(zdperio,'ZoneSubRegion_t'):
                    srname = Internal.getValue(zsr)
                    zsr[0] = 'IDPER'+source+'#%d_%s'%(it,srname)
                    #print('zoneD=',zdperio[0],'zoneR=',srname,'subRegionName=', zsr[0],'. teta_perio=', theta_perioDonor[donor][c], 'present sur bloc', c)
                    #modif pointlist pour calcul sur zone reduite
                    _adaptRange(zname, zsr, infos_PtlistRebuild ) 

                    Internal.createChild(zsr,'RotationAngle' ,'DataArray_t',value=RotAngleDG[c])
                    Internal.createChild(zsr,'RotationCenter','DataArray_t',value=RotCenter)   
                    zdorig[2].append(zsr)

                  ### efface les noeuds sinon perte HPC
                  Internal._rmNodesByType1(zdperio,'ZoneSubRegion_t')
                  #Internal._rmNodesByName(zdperio,'IDPER*')
             c+=1
           '''
           '''

           t1=timeit.default_timer()
           print( "cout adaptcoef= ", t1-t0, 'it= ', it, rcpt)

        #--------------CHECKS----------------------------
        if check:
            #baseRotorCylc = C.node2Center(baseRotorCyl)
            #C.convertPyTree2File(ReceptCyl['Rotor'],"rotor.cgns")
            #C.convertPyTree2File(donorBases['Stator'],"stator.cgns")

            C.convertPyTree2File(ReceptCyl[etage_names[1]], etage_names[1]+'.cgns'  )
            C.convertPyTree2File(donorBases[etage_names[0]], etage_names[0]+'.cgns'  )
    
        Internal._rmNodesByName(tc[rcpt],'GridCoordinates')

    tc_out =  C.newPyTree(['Base'])
    for name in basenames:  tc_out[2].append( tc[name] ) 
    Internal._rmNodesByName(tc_out,'Base')

    vars= ['FlowSolution','CoordinateX','CoordinateY','CoordinateZ']
    for v in vars: C._rmVars(tc_out, v)

    #mise a jour info taille de zone donneuse
    zones = Internal.getZones(tc_out)
    for z in zones:
       nijk    = infos_PtlistRebuild[ z[0]  ] [2]
       z[1][0:4,0] = nijk[0:4]
       z[1][0:4,1] = nijk[0:4]-1 

    return tc_out

#------------------------------------------------------------------
# modify Pointlist and PointlistDonor in order to compute unsteady chimera interpolation on smaller grid
#  z       = zone donneuse
#  s       = subregion
#  dim_old = dictionnary of (ni,nj,nk) of reduced zone
#  dim_new = dictionnary of (ni,nj,nk) of true    zone
#  shift   = dictionnary of (ni,nj,nk) of true    zone
#------------------------------------------------------------------
def _adaptRange(zname, s, infos_Ptlist):

     # zonename du receveur
     recpt = Internal.getValue(s)

     ptlist  = Internal.getNodeFromName( s, "PointList" )[1]
     ptlistD = Internal.getNodeFromName( s, "PointListDonor" )[1]

     shift   = infos_Ptlist[ zname ] [0]
     nijkOpt = infos_Ptlist[ zname ] [1]
     nijk    = infos_Ptlist[ zname ] [2]

     shiftD  = infos_Ptlist[ recpt ] [0]
     nijkOptD= infos_Ptlist[ recpt ] [1]
     nijkD   = infos_Ptlist[ recpt ] [2]

     ninj   = nijkOpt[0]*nijkOpt[1]
     ni     = nijkOpt[0]
     ninjnew= nijk[0]*nijk[1]
     ninew  = nijk[0]

     ninjD   = nijkOptD[0]*nijkOptD[1]
     niD     = nijkOptD[0]
     ninjDnew= nijkD[0]*nijkD[1]
     niDnew  = nijkD[0]

     k     = numpy.empty( ptlist.size, dtype=numpy.int32)
     j     = numpy.empty( ptlist.size, dtype=numpy.int32)
     i     = numpy.empty( ptlist.size, dtype=numpy.int32)
     kD    = numpy.empty( ptlist.size, dtype=numpy.int32)
     jD    = numpy.empty( ptlist.size, dtype=numpy.int32)
     iD    = numpy.empty( ptlist.size, dtype=numpy.int32)



     k[:] =  ptlist[:]//ninj
     j[:] = (ptlist[:] -k[:]*ninj)//ni
     i[:] =  ptlist[:] -k[:]*ninj -j*ni

     kD[:]= ptlistD[:]//ninjD
     jD[:]= (ptlistD[:] -kD[:]*ninjD)//niD
     iD[:]=  ptlistD[:] -kD[:]*ninjD -jD*niD

     #print iD[:]
     '''
     if z[0] == 'D_Stator_10.11' and  recpt == 'D_Rotor_69.1':
       print("shift ,nijkOpt , nijk =", shift, nijkOpt, nijk)
       print("shiftD,nijkOptD, nijkD=", shiftD, nijkOptD, nijkD)
       for l in range(ptlist.size):
           print('ptlist',  ptlist[l], k[l],j[l],i[l] )
     '''

     
     ptlistD[:] = iD[:]+shiftD[0] + (jD[:]+shiftD[1])*niDnew + (kD[:]+shiftD[2])*ninjDnew
     ptlist[:]  =  i[:]+ shift[0] + ( j[:]+ shift[1])*ninew  + ( k[:]+ shift[2])*ninjnew


     return None

#------------------------------------------------------------------
# selection des subzones pour optimiser calcul chimere instationnaire
#------------------------------------------------------------------
def ZonePrecond( base, NGhostCells, info_PtlistRebuild, etages, idir_tg, depth_tg, dim=3):

     zones=[]
     for z in Internal.getZones(base):
        #print('ZONE selectionnee %s'%z[0])
        #on cherche si interface =imin, imax, ...
        connect =Internal.getNodeFromType(z ,'ZoneGridConnectivity_t')
        conns   =Internal.getNodesFromType1(connect ,'GridConnectivity_t')
        for conn in conns:
          contype =Internal.getNodeFromType(conn ,'GridConnectivityType_t')
          if Internal.getValue(contype) == 'Overset':        
             ptrange =Internal.getNodeFromType(conn ,'IndexRange_t')
             #dir=0,...5
             idir = Ghost.getDirection__(dim, [ptrange])

             if idir_tg[ z[0] ]== "imin": dir_tg = 0
             if idir_tg[ z[0] ]== "imax": dir_tg = 1
             if idir_tg[ z[0] ]== "jmin": dir_tg = 2
             if idir_tg[ z[0] ]== "jmax": dir_tg = 3
             if idir_tg[ z[0] ]== "kmin": dir_tg = 4
             if idir_tg[ z[0] ]== "kmax": dir_tg = 5

             depth = depth_tg[ z[0] ]
             #print("Verif",base[0],etages[0],etages[1],idir) 
             
             if idir== dir_tg:
                rg = Internal.getValue(ptrange)
                nijk     = numpy.empty(3, dtype=numpy.int32)
                nijkOpt  = numpy.empty(3, dtype=numpy.int32)
                shift_ijk= numpy.zeros(3, dtype=numpy.int32)
                dimzone  = Internal.getZoneDim(z)
                nijk[0:3]= dimzone[1:4]
                nijk[ :] -=1  #vertex2center

                #print("indice", rg[0,0],rg[0,1], rg[1,0],rg[1,1], rg[2,0],rg[2,1] )          
                if idir == 0:
                   zp = T.subzone(z, (rg[0,0]+NGhostCells,rg[1,0],rg[2,0]), (rg[0,0]+depth,rg[1,1],rg[2,1]))
                elif idir == 1:
                   zp = T.subzone(z, (rg[0,1]-depth,rg[1,0],rg[2,0]), (rg[0,1]-NGhostCells ,rg[1,1],rg[2,1]))
                elif idir == 2:
                   zp = T.subzone(z, (rg[0,0], rg[1,0]+NGhostCells, rg[2,0]),   (rg[0,1],rg[1,0]+depth,rg[2,1]))
                elif idir == 3:
                   zp = T.subzone(z, (rg[0,0], rg[1,1]-depth, rg[2,0]), (rg[0,1] ,rg[1,1]-NGhostCells,rg[2,1]))
                elif idir == 4:
                   zp = T.subzone(z, (rg[0,0], rg[1,0], rg[2,0]+NGhostCells),   (rg[0,1],rg[1,1],rg[2,0]+depth))
                elif idir == 5:
                   zp = T.subzone(z, (rg[0,0], rg[1,0], rg[2,1]-depth), (rg[0,1] , rg[1,1], rg[2,1]-NGhostCells))

                zp[0] = z[0]
                dimzone      = Internal.getZoneDim(zp)
                nijkOpt[0:3] = dimzone[1:4]
                nijkOpt[ :] -=1

                if idir == 0 or idir == 2 or idir == 4:
                   shift_ijk[idir//2] = NGhostCells
                elif idir == 1:
                   shift_ijk[0] = nijk[ 0] - nijkOpt[0] -NGhostCells
                elif idir == 3:
                   shift_ijk[1] = nijk[ 1] - nijkOpt[1] -NGhostCells
                elif idir == 5:
                   shift_ijk[2] = nijk[ 2] - nijkOpt[2] -NGhostCells

                info_PtlistRebuild[zp[0]]  = [ shift_ijk , nijkOpt, nijk ]

                zones.append(zp)
     return zones
