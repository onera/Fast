# Pour Vince
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.Filter as Filter
import Transform.PyTree as T
import Connector.PyTree as X
import numpy
import Distributor2.PyTree as D2

# Nb de tranches
NbTranche = 2
# Largeur d'une tranche
Lz = 0.02

prefixFile = 't_%d.cgns'
prefixFileC = 'tc_%d.cgns'

FILE_DATA   = 'restart_domaine1_6.cgns'
FILE_CONNECT= 'tc_ivan3D.cgns'

numberofproc = 6

# TRANSLATION DE L'ARBRE SOLUTION
t0 = C.convertFile2PyTree(FILE_DATA)

'''
c=0
for z in Internal.getZones(t0):
   D2._addProcNode(z, c)
   procz = Internal.getNodeFromName2(z,'proc')
   procz = Internal.getValue(procz)
   numberofproc = max(procz+1,numberofproc)
   c+=1
'''

print 'Nb de coeurs par tranche = ',  numberofproc

print 'TRANSLATION ET COPIE de RESTART EN TRANCHES'
for tranche in xrange(NbTranche):

   dz =Lz
   if tranche==0: dz =0
   print 'Tranche', tranche,' translatee de ', dz*tranche
   t = T.translate(t0,(0,0,tranche*dz))

   ##################
   # CHANGEMENT DE NOM DES ZONES
   print 'CHANGEMENT DE NOM DES ZONES DES TRANCHES'
   zones = Internal.getZones(t)
   for z in zones:
        z[0]='tr'+str(tranche)+'_'+z[0]
       
        #modif des nom de zone des connectivite
        connect = Internal.getNodeFromType1(z,'ZoneGridConnectivity_t')
        matchs  = Internal.getNodesFromType1(connect,'GridConnectivity1to1_t')
        matchs += Internal.getNodesFromType1(connect,'GridConnectivity_t')
        for match in matchs:
           name = Internal.getValue(match)
           perio= Internal.getNodeFromName1(match,'GridConnectivityProperty')
           name2= str(tranche)
            
           if perio is not None:
              rg = Internal.getNodeFromName1(match,'PointRange')[1]
              # kmin
              if rg[2,0] == 1:
                 name2 = str(tranche-1)
                 if tranche == 0: 
                    name2 = str(NbTranche-1)
                    perio[2][0][2][2][1][2] = perio[2][0][2][2][1][2]*NbTranche
                    print 'PERIO', perio[2][0][2][2][1][2]
                 else:  
                    Internal._rmNodesFromName(match,'GridConnectivityProperty')    
              else: # kmax
                 name2 = str(tranche+1)
                 if tranche == NbTranche-1: 
                    name2 = str('0')
                    perio[2][0][2][2][1][2] = perio[2][0][2][2][1][2]*NbTranche
                    print 'PERIO', perio[2][0][2][2][1][2]
                 else: Internal._rmNodesFromName(match,'GridConnectivityProperty')    

           name ='tr'+name2+'_'+name
           Internal.setValue(match,name)

        # Change le num des procs
        proc = Internal.getNodeFromName2(z,'proc')
	if proc is not None and tranche != 0:
           proc[1][0] = proc[1][0]+ numberofproc

   for pr in xrange(numberofproc*tranche,numberofproc*(tranche+1)):
       tp = C.newPyTree(['Base'])
       zout=[]
       
       for z in zones:
          proc = Internal.getNodeFromName2(z,'proc')
	  if proc is not None:
             if pr==proc[1][0]: 
               zout.append(z)
          else:
             zout.append(z)
       floweq = Internal.getNodeFromName2(t,'FlowEquationSet')
       if floweq is not None: zout.append(floweq)
       refstat= Internal.getNodeFromName2(t,'ReferenceState')
       if floweq is not None: zout.append(refstat)
       tp[2][1][2] += zout
      
       C.convertPyTree2File(tp, './2TRANCHE/'+prefixFile%(pr))

# TRANSLATION DE L'ARBRE CONNECTIVITE
t2 = C.convertFile2PyTree(FILE_CONNECT)

'''
c=0
for z in Internal.getZones(t2):
   D2._addProcNode(z, c)
'''

tp = C.newPyTree(['Base'])
#Internal._rmNodesByName(t,'IDPER*')
print 'TRANSLATION ET COPIE de RESTART EN TRANCHES'
for tranche in xrange(0,NbTranche):

   t = Internal.copyTree(t2)

   dz =Lz
   print 'Tranche', tranche,' translatee de ', dz*tranche
   T._translate(t,(0,0,dz*tranche))

   ##################
   # CHANGEMENT DE NOM DES ZONES
   print 'CHANGEMENT DE NOM DES ZONES DES TRANCHES'
   zones = Internal.getZones(t)
   for z in zones:
        nameSave = z[0]
        z[0]='tr'+str(tranche)+'_'+z[0]

        #modif des nom de zone des connectivite
        matchs  = Internal.getNodesFromType1(z,'ZoneSubRegion_t')
        print "ZNAME",z[0]
        for match in matchs:
           name = Internal.getValue(match)

           #raccord perio
           if name == nameSave:
             #split du raccord en 2
             if NbTranche > 2:

                 ptlist    = Internal.getNodeFromName1(match,'PointList')
                 ptlistD   = Internal.getNodeFromName1(match,'PointListDonor')
                 interpD   = Internal.getNodeFromName1(match,'InterpolantsDonor')
                 interpType= Internal.getNodeFromName1(match,'InterpolantsType')

                 new_sz = numpy.size(ptlist[1])/2

                 ptlist_G     = numpy.empty( new_sz, numpy.int32)
                 ptlist_D     = numpy.empty( new_sz, numpy.int32)
                 ptlistD_G    = numpy.empty( new_sz, numpy.int32)
                 ptlistD_D    = numpy.empty( new_sz, numpy.int32)
                 interpType_G = numpy.ones(  new_sz, numpy.int32)
                 interpType_D = numpy.ones(  new_sz, numpy.int32)
                 interpD_G    = numpy.ones( new_sz, numpy.float64)
                 interpD_D    = numpy.ones( new_sz, numpy.float64)

                 ptlist_G[:]     = ptlist[1][    0 :  new_sz]
                 ptlist_D[:]     = ptlist[1][new_sz:2*new_sz]
                 ptlistD_G[:]    = ptlistD[1][    0 :  new_sz]
                 ptlistD_D[:]    = ptlistD[1][new_sz:2*new_sz]

                 p0     = ptlist[1][0]
                 p1     = ptlist[1][new_sz]
                 
                 ptlist[1] = ptlist_G
                 ptlistD[1]= ptlistD_G
                 interpD[1]= interpD_G
                 interpType[1]= interpType_G

                 match1 = Internal.copyTree(match)
                 ptlist    = Internal.getNodeFromName1(match1,'PointList')
                 ptlistD   = Internal.getNodeFromName1(match1,'PointListDonor')
                 ptlist[1] = ptlist_D
                 ptlistD[1]= ptlistD_D

                 if p0 < p1:
                     trG   = tranche-1
                     trD   = tranche+1
                     if trG==-1: trG = NbTranche-1
                     if trD== NbTranche: trD = 0
                 else:
                     trG   = tranche+1
                     trD   = tranche-1
                     if trD==-1: trD = NbTranche-1
                     if trG== NbTranche: trG = 0

                 name2 = str(trG)
                 nameG='tr'+name2+'_'+name
                 match[0] =  'ID_'+nameG
                 Internal.setValue(match,nameG)
                
                 name2 = str(trD)
                 nameD ='tr'+name2+'_'+name
                 match1[0] =  'ID_'+nameD
                 Internal.setValue(match1,nameD)
                 z[2].append(match1)


             else:
               if tranche==0            : name2 = str(NbTranche-1)
               elif tranche==NbTranche-1: name2 = str(0)
               name ='tr'+name2+'_'+name
               match[0] =  'ID_'+name
               Internal.setValue(match,name)

           else:
             name2= str(tranche)
             name ='tr'+name2+'_'+name
             match[0] =  'ID_'+name
             Internal.setValue(match,name)
        
        # Change le num des procs
        proc = Internal.getNodeFromName2(z,'proc')
	if proc is not None and tranche != 0:
           proc[1][0] = proc[1][0]+ numberofproc

   zout=zones
   tp[2][1][2] += zout

floweq = Internal.getNodeFromName2(t,'FlowEquationSet')
if floweq is not None: tp[2][1][2].append(floweq)
refstat= Internal.getNodeFromName2(t,'ReferenceState')
if floweq is not None: tp[2][1][2].append(refstat)
C.convertPyTree2File(tp,'./2TRANCHE/tc_Nbtranche'+str(NbTranche)+'.cgns')
