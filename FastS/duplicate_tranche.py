# Pour Vince
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.Filter as Filter
import Transform.PyTree as T
import Connector.PyTree as X
import numpy

# Nb de tranches
NbTranche = 2
# Nb de coeurs/tranche
numberofproc=4
# Largeur d'une tranche
Lz = 0.0008

FILE_DATA   = 't_NSLam.cgns'
FILE_CONNECT= 'tc_NSLam.cgns'

# TRANSLATION DE L'ARBRE SOLUTION
t = C.convertFile2PyTree(FILE_DATA)

print 'TRANSLATION ET COPIE de RESTART EN TRANCHES'
for tranche in xrange(NbTranche):

   dz =Lz
   if tranche==0: dz =0
   print 'Tranche', tranche,' translatee de ', dz*tranche
   T._translate(t,(0,0,dz))

   ##################
   # CHANGEMENT DE NOM DES ZONES
   print 'CHANGEMENT DE NOM DES ZONES DES TRANCHES'
   zones = Internal.getNodesFromType(t, 'Zone_t')
   for z in zones:
        if tranche==0:
           z[0]='tr'+str(tranche)+'_'+z[0]
        else:
           tmp  = z[0].split('_')
           name = 'tr'+str(tranche)
           for i in range(1, len(tmp)):
              name+='_'+tmp[i]
           z[0]=name

        #modif des nom de zone des connectivite
        connect = Internal.getNodeFromType1(z,'ZoneGridConnectivity_t')
        matchs  = Internal.getNodesFromType1(connect,'GridConnectivity1to1_t')
        matchs  = matchs +  Internal.getNodesFromType1(connect,'GridConnectivity_t')
        for match in matchs:
           name = Internal.getValue(match)
           perio= Internal.getNodeFromName1(match,'GridConnectivityProperty')
           name2= str(tranche)
           if perio is not None:
              rg = Internal.getNodeFromName1(match,'PointRange')[1]
              # kmin
              if rg[2,0] == 1:
                 name2 = str(tranche-1)
                 if tranche == 0: name2 = str(NbTranche-1)
              # kmax
              else:
                 name2 = str(tranche+1)
                 if tranche == NbTranche-1: name2 = str('0')

           if tranche==0:
              name ='tr'+name2+'_'+name
           else:
              tmp  = name.split('_')
              name = 'tr'+name2
              for i in range(1, len(tmp)):
                 name+='_'+tmp[i]
           Internal.setValue(match,name)
           
        # Change le num des procs
        proc = Internal.getNodeFromName2(z,'proc')
	if proc is not None:
           proc[1][0,0] = proc[1][0,0]+ numberofproc*tranche

   for pr in range(numberofproc*tranche,numberofproc*(tranche+1)):
       tp = C.newPyTree(['Base'])
       zout=[]
       for z in zones:
          proc = Internal.getNodeFromName2(z,'proc')
	  if proc is not None:
             if pr==proc[1][0,0]: 
               zout.append(z)
          else:
            zout.append(z)
       floweq = Internal.getNodeFromName2(t,'FlowEquationSet')
       if floweq is not None: zout.append(floweq)
       refstat= Internal.getNodeFromName2(t,'ReferenceState')
       if floweq is not None: zout.append(refstat)
       tp[2][1][2] += zout
            
       C.convertPyTree2File(tp, 't'+str(pr)+'.cgns')


# TRANSLATION DE L'ARBRE CONNECTIVITE
t = C.convertFile2PyTree(FILE_CONNECT)

print 'TRANSLATION ET COPIE de RESTART EN TRANCHES'
for tranche in xrange(NbTranche):

   dz =Lz
   if tranche==0: dz =0
   print 'Tranche', tranche,' translatee de ', dz*tranche
   T._translate(t,(0,0,dz))


   ##################
   # CHANGEMENT DE NOM DES ZONES
   print 'CHANGEMENT DE NOM DES ZONES DES TRANCHES'
   zones = Internal.getNodesFromType(t, 'Zone_t')
   for z in zones:
        if tranche==0:
           z[0]='tr'+str(tranche)+'_'+z[0]
        else:
           tmp  = z[0].split('_')
           name = 'tr'+str(tranche)
           for i in range(1, len(tmp)):
              name+='_'+tmp[i]
           z[0]=name

        #modif des nom de zone des connectivite
        matchs  = Internal.getNodesFromType1(z,'ZoneSubRegion_t')
        for match in matchs:
           name = Internal.getValue(match)
           name2= str(tranche)

           if tranche != 0:
              tmp  = name.split('_')
              name =''
              for i in range(1, len(tmp)):
                    if i ==  len(tmp) -1:
                      name+= tmp[i]
                    else:
                      name+= tmp[i]+'_'

           if match[0][0:5] == 'IDPER' and tranche>=0:
              donor = Internal.getNodeFromName1(match,'PointListDonor')[1]
              kmin = donor[0]
              kmax = donor[donor.size -1]
              matchNew = Internal.copyTree(match)
              #print 'period',matchNew[0],  matchNew[2]
              # kmin
              if kmin < kmax:
                 trM =  'tr'+str(tranche-1)+'_'+name 
                 trP =  'tr'+str(tranche+1)+'_'+name 
                 if tranche == 0: 
                     trM = 'tr'+str(NbTranche-1)+'_'+name

                 elif tranche == NbTranche-1: 
                     trP = 'tr0_'+name

                 match[0]   = 'IDPERM' + trM
                 matchNew[0]= 'IDPERP' + trP
                 Internal.setValue(match   ,trM)
                 Internal.setValue(matchNew,trP)
              # kmax
              else:
                 trM =  'tr'+str(tranche+1)+'_'+name 
                 trP =  'tr'+str(tranche-1)+'_'+name 
                 if tranche == 0: 
                     trP = 'tr'+str(NbTranche-1)+'_'+name

                 elif tranche == NbTranche-1: 
                     trM = 'tr0_'+name

                 match[0]   = 'IDPERP' + trM
                 matchNew[0]= 'IDPERM' + trP
                 Internal.setValue(match   ,trM)
                 Internal.setValue(matchNew,trP)

              Internal.createChild(z, matchNew[0], 'ZoneSubRegion_t', value=trP, children=matchNew[2])

              #resize des numpy: 1 numpy(N)=(N/2)
              np_list=['PointList','PointListDonor','InterpolantsDonor','InterpolantsType']
              for var in np_list:
                 tmp1  = Internal.getNodeFromName1(match,var)
                 tmp2  = Internal.getNodeFromName1(matchNew,var)
                 size  = tmp1[1].size/2
                 if var[0:7]== 'Interp':
                    newtab1 = numpy.empty(size, dtype=numpy.float64)
                    newtab2 = numpy.empty(size, dtype=numpy.float64)
                 else:
                    newtab1 = numpy.empty(size, dtype=numpy.int32)
                    newtab2 = numpy.empty(size, dtype=numpy.int32)

                 if kmin < kmax:
                    newtab2[0:size]    = tmp1[1][0:size]
                    newtab1[0:size]    = tmp1[1][size:2*size]
                 else:
                    newtab2[0:size]    = tmp1[1][size:2*size]
                    newtab1[0:size]    = tmp1[1][0:size]

                 tmp1[1]= newtab1
                 tmp2[1]= newtab2
                
           else:
              name ='tr'+name2+'_'+name
              match[0] =  'ID_'+name
              Internal.setValue(match,name)
           
        # Change le num des procs
        proc = Internal.getNodeFromName2(z,'proc')
	if proc is not None:
           proc[1][0,0] = proc[1][0,0]+ numberofproc*tranche

   if tranche == 0:
      C.convertPyTree2File(t, 'tc_'+str(NbTranche)+'Tranches.cgns')
   else:
      Filter.writeNodesFromPaths('tc_'+str(NbTranche)+'Tranches.cgns', ['CGNSTree/Base']*len(zones), zones, format='bin_cgns')


