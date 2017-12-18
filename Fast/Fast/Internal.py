"""Services for FAST solvers.
"""
import numpy
import os
try:
    import Converter.PyTree as C
    import Converter.Internal as Internal
    import Post.PyTree as P
except:
    raise ImportError("Fast.Internal: requires Converter and Post module.")


try:
    OMP_NUM_THREADS = os.environ['OMP_NUM_THREADS']
    OMP_NUM_THREADS = int(OMP_NUM_THREADS)
except: OMP_NUM_THREADS = 1

#==============================================================================
# Met un dictionnaire de numerics dans une/des zones
#==============================================================================
def setNum2Zones(t, num):
    """Set the numeric dictionary to zones."""
    tp = Internal.copyRef(t)
    _setNum2Zones(tp, num)
    return tp

def _setNum2Zones(a, num):
    zones = Internal.getNodesFromType2(a, 'Zone_t')
    for z in zones:
        cont = Internal.createUniqueChild(z, '.Solver#define', 
                                          'UserDefinedData_t')
        for k in num.keys(): # some checks?
            Internal.createUniqueChild(cont, k, 'DataArray_t', num[k])
    return None

#==============================================================================
def setNum2Base(t, num):
    """Set the numeric dictionary to bases."""
    tp = Internal.copyRef(t)
    _setNum2Base(tp, num)
    return tp

def _setNum2Base(a, num):
    bases = Internal.getNodesFromType2(a, 'CGNSBase_t')
    for b in bases:
        cont = Internal.createUniqueChild(b, '.Solver#define', 
                                          'UserDefinedData_t')
        for k in num.keys(): # some checks?
            Internal.createUniqueChild(cont, k, 'DataArray_t', num[k])
    return None

#==============================================================================
# Reorder zone pour // openmp mode=1
#==============================================================================
def _reorder(t, tc=None, omp_mode=0):

    # reordone les zones pour garantir meme ordre entre t et tc
    if tc is not None: 
       bases_tc = Internal.getNodesFromType1(tc,'CGNSBase_t')
       bases    = Internal.getNodesFromType1(t, 'CGNSBase_t')
       zones_tc = Internal.getNodesFromType2(tc,'Zone_t')
       zones    = Internal.getNodesFromType2(t, 'Zone_t')
 
       for base_tc in bases_tc: 

           zones_tc = Internal.getNodesFromType1(base_tc,'Zone_t')
 
           #on cherche le matching entre base t et tc
           c  =0
           ctg=0
           for base in bases: 
              if base[0] == base_tc[0]: ctg=c
              c+=1  

           zones    = Internal.getNodesFromType1(bases[ctg], 'Zone_t')

           zname_tc = []
           for zc in zones_tc: zname_tc.append( zc[0] )

           new_zones_tc = []
           for z in zones:
              pos = zname_tc.index( z[0] )
              new_zones_tc.append( zones_tc[pos] )

           l = base_tc[2]
           orig = []
           for i in l:
	        if i[3] != 'Zone_t': orig.append(i)
                base_tc[2] = new_zones_tc + orig

    # reorder_zone pour equilibrage omp
    if omp_mode != 0:
       zones = Internal.getNodesFromType2(t, 'Zone_t')
       if tc is not None: zones_tc = Internal.getNodesFromType2(tc, 'Zone_t')

       nzone =len(zones)
       size_zone =[]
       for z in zones:
            dim = Internal.getZoneDim(z)
            if dim[0]=='Structured':
               if dim[3] == 2: kfic = 1
               else          : kfic = 5
               ndimdx = (dim[1]-5)*(dim[2]-5)*(dim[3]-kfic) 
            else: ndimdx = dim[2]

            size_zone.append(ndimdx)
          
       target = sum(size_zone)/OMP_NUM_THREADS + 1

       new_zones       =[]
       new_zones_tc    =[]
       equi            = numpy.zeros((OMP_NUM_THREADS), numpy.int32)
       zone_par_thread = nzone/OMP_NUM_THREADS
       for z in xrange(nzone):
             th = z%OMP_NUM_THREADS
             vmax = max( size_zone )
             vmin = min( size_zone)
             pos_max = size_zone.index( vmax )
             pos_min = size_zone.index( vmin )
             if equi[th]+vmax > target*1.00:
                 equi[th] =  equi[th] + vmin
                 #print "vmin", vmin, equi[th], th
                 new_zones.append( zones[pos_min])
                 if tc is not None: 
                    new_zones_tc.append( zones_tc[pos_min])
                    del zones_tc[pos_min]
                 del size_zone[pos_min]
                 del zones[pos_min]
             else:
                 #print "vmax", vmax, equi[th], th
                 equi[th] =  equi[th] + vmax
                 new_zones.append( zones[pos_max])
                 if tc is not None:
                    new_zones_tc.append( zones_tc[pos_max])
                    del zones_tc[pos_max]
                 del size_zone[pos_max]
                 del zones[pos_max]

       l = t[2][1][2]
       orig = []
       for i in l:
	  if i[3] != 'Zone_t': orig.append(i)

       t[2][1][2]=  orig+new_zones

       if tc is not None: 
         l = tc[2][1][2]
         orig = []
         for i in l:
	    if i[3] != 'Zone_t': orig.append(i)
         tc[2][1][2]= new_zones_tc + orig

       for th in range(OMP_NUM_THREADS): print "count th=",equi[th] , th+1


#==============================================================================
# Init/create primitive Variable 
#==============================================================================
def _createPrimVars(base, zone, FIRST_IT, omp_mode, rmConsVars=True, adjoint=False):
    # Get model
    model = "Euler"
    a = Internal.getNodeFromName2(zone, 'GoverningEquations')
    if a is not None: model = Internal.getValue(a)
    else: 
        a = Internal.getNodeFromName2(base, 'GoverningEquations')
        if a is not None: model = Internal.getValue(a)

    # Get adim
    adim = None
    a = Internal.getNodeFromName2(zone, 'ReferenceState')
    if a is not None: adim = C.getState(a)
    else:
        a = Internal.getNodeFromName2(base, 'ReferenceState')
        if a is not None: adim = C.getState(a)
    if adim is None: raise ValueError("createPrimVars: ReferenceState is missing in t.")

    if (not (model == 'Euler' or model == 'euler')):
        if C.isNamePresent(zone, 'centers:ViscosityEddy') != 1: 
            C._initVars(zone, 'centers:ViscosityEddy', 0.)

    if (not (model == 'NSTurbulent' or model == 'nsspalart')): sa = False
    else: 
        sa = True
        if C.isNamePresent(zone, 'centers:TurbulentDistance') != 1: 
            raise ValueError("createPrimVars: TurbulentDistance field required at cell centers for RANS computations.")
    
    rgp = (adim[11]-1.)*adim[7]
    if C.isNamePresent(zone, 'centers:VelocityX'  ) != 1: P._computeVariables(zone, ['centers:VelocityX'  ], rgp=rgp)
    if C.isNamePresent(zone, 'centers:VelocityY'  ) != 1: P._computeVariables(zone, ['centers:VelocityY'  ], rgp=rgp)
    if C.isNamePresent(zone, 'centers:VelocityZ'  ) != 1: P._computeVariables(zone, ['centers:VelocityZ'  ], rgp=rgp)
    if C.isNamePresent(zone, 'centers:Temperature') != 1: P._computeVariables(zone, ['centers:Temperature'], rgp=rgp)

    if (sa and C.isNamePresent(zone, 'centers:TurbulentSANuTilde') != 1): 
        if C.isNamePresent(zone, 'centers:TurbulentSANuTildeDensity') != 1: C._initVars(zone, 'centers:TurbulentSANuTilde= %20.16g/{centers:Density}'%adim[14])
        else: C._initVars(zone, 'centers:TurbulentSANuTilde= {centers:TurbulentSANuTildeDensity}/{centers:Density}')

    if rmConsVars:
        C._rmVars(zone, 'centers:MomentumX')
        C._rmVars(zone, 'centers:MomentumY')
        C._rmVars(zone, 'centers:MomentumZ')
        C._rmVars(zone, 'centers:EnergyStagnationDensity')
        C._rmVars(zone, 'centers:TurbulentSANuTildeDensity')

    #on test s'il existe 2 niveau en temps dans l'arbre pour appliquer la bonne formule de derivee temporelle a la premere iteration
    FIRST_IT = 1
    if C.isNamePresent(zone, 'centers:Density_M1')   != 1: C._cpVars(zone, 'centers:Density'  , zone, 'centers:Density_M1')  ; FIRST_IT=0
    if C.isNamePresent(zone, 'centers:VelocityX_M1') != 1: C._cpVars(zone, 'centers:VelocityX', zone, 'centers:VelocityX_M1'); FIRST_IT=0
    if C.isNamePresent(zone, 'centers:VelocityY_M1') != 1: C._cpVars(zone, 'centers:VelocityY', zone, 'centers:VelocityY_M1'); FIRST_IT=0
    if C.isNamePresent(zone, 'centers:VelocityZ_M1') != 1: C._cpVars(zone, 'centers:VelocityZ', zone, 'centers:VelocityZ_M1'); FIRST_IT=0
    if C.isNamePresent(zone, 'centers:Temperature_M1') != 1: C._cpVars(zone, 'centers:Temperature', zone, 'centers:Temperature_M1'); FIRST_IT=0
    if (sa and  C.isNamePresent(zone, 'centers:TurbulentSANuTilde_M1') != 1): C._cpVars(zone, 'centers:TurbulentSANuTilde', zone, 'centers:TurbulentSANuTilde_M1'); FIRST_IT=0

    if C.isNamePresent(zone, 'centers:Density_P1')   != 1: C._cpVars(zone, 'centers:Density'  , zone, 'centers:Density_P1')  
    if C.isNamePresent(zone, 'centers:VelocityX_P1') != 1: C._cpVars(zone, 'centers:VelocityX', zone, 'centers:VelocityX_P1')
    if C.isNamePresent(zone, 'centers:VelocityY_P1') != 1: C._cpVars(zone, 'centers:VelocityY', zone, 'centers:VelocityY_P1')
    if C.isNamePresent(zone, 'centers:VelocityZ_P1') != 1: C._cpVars(zone, 'centers:VelocityZ', zone, 'centers:VelocityZ_P1')
    if C.isNamePresent(zone, 'centers:Temperature_P1') != 1: C._cpVars(zone, 'centers:Temperature', zone, 'centers:Temperature_P1')
    if (sa and  C.isNamePresent(zone, 'centers:TurbulentSANuTilde_P1') != 1): C._cpVars(zone, 'centers:TurbulentSANuTilde', zone, 'centers:TurbulentSANuTilde_P1')

    sfd = 0
    a = Internal.getNodeFromName2(zone, 'sfd')
    if a is not None: sfd = Internal.getValue(a)
    if sfd == 1:
       if C.isNamePresent(zone, 'centers:Density_f')   != 1:   C._initVars(zone, 'centers:Density_f', 0.)
       if C.isNamePresent(zone, 'centers:VelocityX_f') != 1:   C._initVars(zone, 'centers:VelocityX_f', 0.)
       if C.isNamePresent(zone, 'centers:VelocityY_f') != 1:   C._initVars(zone, 'centers:VelocityY_f', 0.)
       if C.isNamePresent(zone, 'centers:VelocityZ_f') != 1:   C._initVars(zone, 'centers:VelocityZ_f', 0.)
       if C.isNamePresent(zone, 'centers:Temperature_f') != 1: C._initVars(zone, 'centers:Temperature_f', 0.)
       if (sa and  C.isNamePresent(zone, 'centers:TurbulentSANuTilde_f') != 1): C._initVars(zone, 'centers:TurbulentSANuTilde_f', 0.)


   # TEST PYTHON CONTEXTE ADJOINT ? (dpCLp, dpCDp)  -----------------------------
    if (adjoint==True): 
       if C.isNamePresent(zone, 'centers:dpCLp_dpDensity') != 1: C._initVars(zone, 'centers:dpCLp_dpDensity', 0.)
       if C.isNamePresent(zone, 'centers:dpCLp_dpMomentumX') != 1: C._initVars(zone, 'centers:dpCLp_dpMomentumX', 0.)
       if C.isNamePresent(zone, 'centers:dpCLp_dpMomentumY') != 1: C._initVars(zone, 'centers:dpCLp_dpMomentumY', 0.)
       if C.isNamePresent(zone, 'centers:dpCLp_dpMomentumZ') != 1: C._initVars(zone, 'centers:dpCLp_dpMomentumZ', 0.)
       if C.isNamePresent(zone, 'centers:dpCLp_dpEnergyStagDens') != 1: C._initVars(zone, 'centers:dpCLp_dpEnergyStagDens', 0.)

       if C.isNamePresent(zone, 'centers:dpCDp_dpDensity') != 1: C._initVars(zone, 'centers:dpCDp_dpDensity', 0.)
       if C.isNamePresent(zone, 'centers:dpCDp_dpMomentumX') != 1: C._initVars(zone, 'centers:dpCDp_dpMomentumX', 0.)
       if C.isNamePresent(zone, 'centers:dpCDp_dpMomentumY') != 1: C._initVars(zone, 'centers:dpCDp_dpMomentumY', 0.)
       if C.isNamePresent(zone, 'centers:dpCDp_dpMomentumZ') != 1: C._initVars(zone, 'centers:dpCDp_dpMomentumZ', 0.)
       if C.isNamePresent(zone, 'centers:dpCDp_dpEnergyStagDens') != 1: C._initVars(zone, 'centers:dpCDp_dpEnergyStagDens', 0.)
       # 
       # On a fait passer rhs en variable locale (inutile ici)
       if C.isNamePresent(zone, 'centers:rhsIterAdjCLp_RDensity') != 1: C._initVars(zone, 'centers:rhsAdjCLp_RDensity', 0.)
       if C.isNamePresent(zone, 'centers:rhsIterAdjCLp_RMomentumX') != 1: C._initVars(zone, 'centers:rhsAdjCLp_RMomentumX', 0.)
       if C.isNamePresent(zone, 'centers:rhsIterAdjCLp_RMomentumY') != 1: C._initVars(zone, 'centers:rhsAdjCLp_RMomentumY', 0.)
       if C.isNamePresent(zone, 'centers:rhsIterAdjCLp_RMomentumZ') != 1: C._initVars(zone, 'centers:rhsAdjCLp_RMomentumZ', 0.)
       if C.isNamePresent(zone, 'centers:rhsIterAdjCLp_REnergyStagDens') != 1: C._initVars(zone, 'centers:rhsAdjCLp_REnergyStagDens', 0.)
       
       if C.isNamePresent(zone, 'centers:rhsIterAdjCDp_RDensity') != 1: C._initVars(zone, 'centers:rhsAdjCDp_RDensity', 0.)
       if C.isNamePresent(zone, 'centers:rhsIterAdjCDp_RMomentumX') != 1: C._initVars(zone, 'centers:rhsAdjCDp_RMomentumX', 0.)
       if C.isNamePresent(zone, 'centers:rhsIterAdjCDp_RMomentumY') != 1: C._initVars(zone, 'centers:rhsAdjCDp_RMomentumY', 0.)
       if C.isNamePresent(zone, 'centers:rhsIterAdjCDp_RMomentumZ') != 1: C._initVars(zone, 'centers:rhsAdjCDp_RMomentumZ', 0.)
       if C.isNamePresent(zone, 'centers:rhsIterAdjCDp_REnergyStagDens') != 1: C._initVars(zone, 'centers:rhsAdjCDp_REnergyStagDens', 0.)
       #

       #
       if C.isNamePresent(zone, 'centers:AdjCLp_RDensity') != 1: C._initVars(zone, 'centers:AdjCLp_RDensity', 0.)
       if C.isNamePresent(zone, 'centers:AdjCLp_RMomentumX') != 1: C._initVars(zone, 'centers:AdjCLp_RMomentumX', 0.)
       if C.isNamePresent(zone, 'centers:AdjCLp_RMomentumY') != 1: C._initVars(zone, 'centers:AdjCLp_RMomentumY', 0.)
       if C.isNamePresent(zone, 'centers:AdjCLp_RMomentumZ') != 1: C._initVars(zone, 'centers:AdjCLp_RMomentumZ', 0.)
       if C.isNamePresent(zone, 'centers:AdjCLp_REnergyStagDens') != 1: C._initVars(zone, 'centers:AdjCLp_REnergyStagDens', 0.)

       if C.isNamePresent(zone, 'centers:AdjCDp_RDensity') != 1: C._initVars(zone, 'centers:AdjCDp_RDensity', 0.)
       if C.isNamePresent(zone, 'centers:AdjCDp_RMomentumX') != 1: C._initVars(zone, 'centers:AdjCDp_RMomentumX', 0.)
       if C.isNamePresent(zone, 'centers:AdjCDp_RMomentumY') != 1: C._initVars(zone, 'centers:AdjCDp_RMomentumY', 0.)
       if C.isNamePresent(zone, 'centers:AdjCDp_RMomentumZ') != 1: C._initVars(zone, 'centers:AdjCDp_RMomentumZ', 0.)
       if C.isNamePresent(zone, 'centers:AdjCDp_REnergyStagDens') != 1: C._initVars(zone, 'centers:AdjCDp_REnergyStagDens', 0.)
       #
       if C.isNamePresent(zone, 'centers:incAdj_RDensity') != 1: C._initVars(zone, 'centers:incAdj_RDensity', 0.)
       if C.isNamePresent(zone, 'centers:incAdj_RMomentumX') != 1: C._initVars(zone, 'centers:incAdj_RMomentumX', 0.)
       if C.isNamePresent(zone, 'centers:incAdj_RMomentumY') != 1: C._initVars(zone, 'centers:incAdj_RMomentumY', 0.)
       if C.isNamePresent(zone, 'centers:incAdj_RMomentumZ') != 1: C._initVars(zone, 'centers:incAdj_RMomentumZ', 0.)
       if C.isNamePresent(zone, 'centers:incAdj_REnergyStagDens') != 1: C._initVars(zone, 'centers:incAdj_REnergyStagDens', 0.)
       #
       if C.isNamePresent(zone, 'dpCLp_dpX') != 1: C._initVars(zone, 'dpCLp_dpX', 0.)
       if C.isNamePresent(zone, 'dpCLp_dpY') != 1: C._initVars(zone, 'dpCLp_dpY', 0.)
       if C.isNamePresent(zone, 'dpCLp_dpZ') != 1: C._initVars(zone, 'dpCLp_dpZ', 0.)

       if C.isNamePresent(zone, 'dpCDp_dpX') != 1: C._initVars(zone, 'dpCDp_dpX', 0.)
       if C.isNamePresent(zone, 'dpCDp_dpY') != 1: C._initVars(zone, 'dpCDp_dpY', 0.)
       if C.isNamePresent(zone, 'dpCDp_dpZ') != 1: C._initVars(zone, 'dpCDp_dpZ', 0.)

    # fin du if adjoint ------------------------------------
    return sa, FIRST_IT

#==============================================================================
# create work arrays
#==============================================================================
def createWorkArrays__(zones, dtloc, FIRST_IT):
    ndimt = 0 ; ndimcoe = 0 ; ndimwig = 0 ; c = 0;
    #global  CACHELINE
    # on check sur 1ere zone si implicite: mixte impli/expli impossible pour le moment
    scheme = "implicit"
    a = Internal.getNodeFromName2(zones[0], 'temporal_scheme')
    if a is not None: scheme = Internal.getValue(a)

    #neq = 5
    #if (model == 'nsspalart' or model =='NSTurbulent'): neq = 6

    if scheme == 'implicit': lssiter_loc = 0
    elif scheme == 'implicit_local':  lssiter_loc = 1  # sous-iteration local
    else: lssiter_loc = 0

    neq_max =0
    for z in zones:
        model = "Euler"
        a = Internal.getNodeFromName2(z, 'model')
        model  = Internal.getValue(a)
        neq = 5
        if (model == 'nsspalart' or model =='NSTurbulent'): neq = 6
        
        if scheme == 'implicit' or scheme == 'implicit_local':   neq_coe = neq     
        else:                                                    neq_coe = 1    # explicit

        neq_max = max(neq, neq_max)

        # dim from getZoneDim
        dims = Internal.getZoneDim(z)
        if dims[0]=='Structured' : nijk = (dims[1]-1)*(dims[2]-1)*(dims[3]-1)
        else                     : nijk = dims[2]
        ndimt   +=     neq*nijk       # surdimensionne ici
        ndimcoe += neq_coe*nijk       # surdimensionne ici
        ndimwig +=        3*nijk      # surdimensionne ici
        c       += 1

    mx_sszone   = 20
    mx_synchro  = 2000 
    mx_thread   = OMP_NUM_THREADS       # surdimensionne : doit etre = a OMP_NUM_THREADS
    verrou      = mx_sszone*c*mx_synchro*mx_thread

#    wig   = KCore.empty(ndimwig, CACHELINE)
#    coe   = KCore.empty(ndimcoe, CACHELINE)
#    drodm = KCore.empty(ndimt  , CACHELINE)
    wig   = numpy.empty(ndimwig, dtype=numpy.float64)
    coe   = numpy.empty(ndimcoe, dtype=numpy.float64)
    drodm = numpy.empty(ndimt  , dtype=numpy.float64)


    #drodm1= numpy.empty(ndimt  , dtype=numpy.float64)
    #drodm1[rodm = numpy.empty(ndimt  , dtype=numpy.float640:ndimt] = 1.

#    for z in zones:
#         param_int = Internal.getNodeFromName2(z, 'Parameter_int')  # noeud
#    
#    ptr = drodm1.reshape((nijk*5), order='Fortran') # no copy I hope
#    for c in xrange(5):
#        print'c=',c
#        #print 'a0',a[0]
#        # Copy elements
#        fasts.initNuma( ptr, drodm, param_int, c )
#        # Replace numpys with slice
#        #drodm1 = drodm[c*nijk:(c+1)*nijk]
#        #a[1] = a[1].reshape(sh, order='Fortran')
#
#    drodm1 = drodm[0:ndimt]
#    print "verif", ndimt, (132**3)*5, drodm1.shape
#    #drodm1 = drodm1.reshape((132,132,132,5), order='Fortran')
#    drodm1 = drodm1.reshape((132*132*132*5), order='Fortran')
#    print drodm1.shape

    lok      = numpy.empty(verrou     , dtype=numpy.int32  )
    iskip_lu = numpy.empty(dtloc[0]*2 , dtype=numpy.int32  )   # (dtloc[0] = nitmax

    
    hook = [wig, coe, drodm, lok, iskip_lu, dtloc, lssiter_loc, mx_sszone, mx_synchro, FIRST_IT, neq_max]
    
    return hook

#==============================================================================
# Retourne un tag suivant bcname
#==============================================================================
def tagBC(bcname):
  if   bcname == "BCExtrapolate"    :       tag = 0;
  elif bcname == "BCFarfield":              tag = 1;
  elif bcname == "BCInflowSupersonic":      tag = 2;
  elif bcname == "BCWallInviscid":          tag = 3;
  elif bcname == "BCSymmetryPlane":         tag = 3;
  elif bcname == "BCWall":                  tag = 4;
  elif bcname == "FamilySpecified":         tag = 5;
  elif bcname == "BCWallViscous":           tag = 6;
  elif bcname == "BCWallViscous_isot_fich": tag = 7;
  elif bcname == "BC_fich":                 tag = 8;
  elif bcname == "Nearmatch":               tag = 9;
  elif bcname == "BCOutflow":               tag =10;
  elif bcname == "BCautoperiod":            tag =11;
  elif bcname == "BCWallViscous_transition":tag =12;
  elif bcname == "BCInflow":                tag =13;
  elif bcname == "BCExtrapolateRANS":       tag =14;
  elif bcname  == "BCPeriodic":             tag =15;
  elif bcname == "BCOutpres":               tag =16;
  elif bcname == "BCInj1":                  tag =17;
  elif bcname == "BCDegenerateLine":        tag = 0;
  else: print "Warning: unknown BC type %s."%bcname
  return tag

#==============================================================================
# Retourne la liste des parametres necessaires aux transferts IBC 
# OUT: [bcType, Gamma, Cv, MuS, Cs, Ts]
# Pour l'instant, on recupere ces variables du premier reference state
#==============================================================================
def getIBCInfo__(t):
    model = "Euler"
    zones = Internal.getZones(t)
    a = Internal.getNodeFromName2(zones[0], 'model')
    model = Internal.getValue(a)
    if model == "Euler": bcType = 0
    elif model == "NSLaminar": bcType = 3
    else: bcType = 3 # 1: lineaire, 2: loi log, 3: Musker

    state = Internal.getNodeFromName(t, 'ReferenceState')
    if state is None: raise ValueError("getIBCInfo: ReferenceState is missing in t.")
    
    Gamma = Internal.getNodeFromName1(state,'Gamma')
    Gamma = Internal.getValue(Gamma)
    cvInf = Internal.getNodeFromName1(state,'Cv')
    cvInf = Internal.getValue(cvInf)
    Mus   = Internal.getNodeFromName1(state,'Mus')
    Mus   = Internal.getValue(Mus)
    Cs    = Internal.getNodeFromName1(state,'Cs')
    Cs    = Internal.getValue(Cs)
    Ts    = Internal.getNodeFromName1(state,'Ts')
    Ts    = Internal.getValue(Ts)
    return [bcType, Gamma, cvInf, Mus, Cs, Ts]

#==============================================================================
# echange M1 <- current, current <- P1, P1 <- M1
#==============================================================================
def switchPointers__(zones, order=3):
    
    if order == 1 or order == 3:
        for z in zones:
            sol =  Internal.getNodeFromName1(z, 'FlowSolution#Centers')
            own =  Internal.getNodeFromName1(z, '.Solver#ownData')
            ca = Internal.getNodeFromName1(sol, 'Density')
            cb = Internal.getNodeFromName1(sol, 'VelocityX')
            cc = Internal.getNodeFromName1(sol, 'VelocityY')
            cd = Internal.getNodeFromName1(sol, 'VelocityZ')
            ce = Internal.getNodeFromName1(sol, 'Temperature')

            caM1 = Internal.getNodeFromName1(sol, 'Density_M1')
            cbM1 = Internal.getNodeFromName1(sol, 'VelocityX_M1')
            ccM1 = Internal.getNodeFromName1(sol, 'VelocityY_M1')
            cdM1 = Internal.getNodeFromName1(sol, 'VelocityZ_M1')
            ceM1 = Internal.getNodeFromName1(sol, 'Temperature_M1')

            caP1 = Internal.getNodeFromName1(sol, 'Density_P1')
            cbP1 = Internal.getNodeFromName1(sol, 'VelocityX_P1')
            ccP1 = Internal.getNodeFromName1(sol, 'VelocityY_P1')
            cdP1 = Internal.getNodeFromName1(sol, 'VelocityZ_P1')
            ceP1 = Internal.getNodeFromName1(sol, 'Temperature_P1')

            # sauvegarde M1
            ta = caM1[1]; tb = cbM1[1]; tc = ccM1[1]; td = cdM1[1]; te = ceM1[1]

            # M1 <- current
            caM1[1] = ca[1]; cbM1[1] = cb[1]; ccM1[1] = cc[1]; cdM1[1] = cd[1]; ceM1[1] = ce[1]

            # current <- P1
            ca[1] = caP1[1]; cb[1] = cbP1[1]; cc[1] = ccP1[1]; cd[1] = cdP1[1]; ce[1] = ceP1[1]

            # P1 <- temp
            caP1[1] = ta; cbP1[1] = tb; ccP1[1] = tc; cdP1[1] = td; ceP1[1] = te

            model  = Internal.getNodeFromName1(own, 'model')
            model  = Internal.getValue(model)

            if model == 'nsspalart' or model=='NSTurbulent': 
                cf   =  Internal.getNodeFromName1(sol, 'TurbulentSANuTilde')
                cfM1 =  Internal.getNodeFromName1(sol, 'TurbulentSANuTilde_M1')
                cfP1 =  Internal.getNodeFromName1(sol, 'TurbulentSANuTilde_P1')
                tf   = cfM1[1]; cfM1[1] = cf[1]; cf[1] = cfP1[1]; cfP1[1] = tf


    elif order == 2:
        for z in zones:
            caP1 = Internal.getNodeFromName2(z, 'Density_P1')
            cbP1 = Internal.getNodeFromName2(z, 'VelocityX_P1')
            ccP1 = Internal.getNodeFromName2(z, 'VelocityY_P1')
            cdP1 = Internal.getNodeFromName2(z, 'VelocityZ_P1')
            ceP1 = Internal.getNodeFromName2(z, 'Temperature_P1')
        
            caM1 = Internal.getNodeFromName2(z, 'Density_M1')
            cbM1 = Internal.getNodeFromName2(z, 'VelocityX_M1')
            ccM1 = Internal.getNodeFromName2(z, 'VelocityY_M1')
            cdM1 = Internal.getNodeFromName2(z, 'VelocityZ_M1')
            ceM1 = Internal.getNodeFromName2(z, 'Temperature_M1')

            caP1[1] = caM1[1]; cbP1[1] = cbM1[1]; ccP1[1] = ccM1[1]; cdP1[1] = cdM1[1]; ceP1[1] = ceM1[1]

            model  = Internal.getNodeFromName2(z, 'model')
            model  = Internal.getValue(model)

            if (model == 'nsspalart' or model=='NSTurbulent'): 
                cfM1 =  Internal.getNodeFromName2(z, 'TurbulentSANuTilde_M1')
                cfP1 =  Internal.getNodeFromName2(z, 'TurbulentSANuTilde_P1')
                cfP1[1] = cfM1[1]

#==============================================================================
# Compact vars in Container or defined in fields
#==============================================================================
def getField2Compact__(z, field):
    field = field.split(':')
    if len(field) == 1: 
        container = Internal.__FlowSolutionNodes__ ; name = field[0]
    elif field[0] == 'nodes': 
        container = Internal.__FlowSolutionNodes__ ; name = field[1]
    else: 
        container = Internal.__FlowSolutionCenters__ ; name = field[1]
    conts = Internal.getNodesFromName(z, container)
    for c in conts:
        ar = Internal.getNodeFromName1(c, name)
        if ar is not None: return ar
    return None
    
def getFields2Compact__(z, containers, fields):
    out = []
    if fields is not None:
        for n in fields:
            v = getField2Compact__(z, n)
            if v is not None: out.append(v)
    else: # containeurs
        conts = []
        for c in containers:
            conts += Internal.getNodesFromName(z, c)
        for c in conts:
            ars = Internal.getNodesFromType(c, 'DataArray_t')
            out += ars
    return out

#==============================================================================
# Retourne un dictionnaire num a partir des donnees solveur d'une zone
#==============================================================================
def getNumFromTag(t):
    num = {}
    zones = Internal.getNodesFromType(t, 'Zone_t')
    if zones == []: return num
    z = zones[0]
    node = Internal.getNodesFromName1(z, '.Solver#define')
    if node == []: return num
    node = node[0]
    for n in node[2]: num[n[0]] = Internal.getValue(n)
    return num

#==============================================================================
# Retourne la valeur d'un champ de tag
# IN: name: champ voulu (solver...)
#==============================================================================
def getValueFromTag(t, name):
    zones = Internal.getNodesFromType(t, 'Zone_t')
    if zones == []: return None
    z = zones[0]
    node = Internal.getNodesFromName1(z, '.Solver#define')
    if node == []: return None
    node = node[0]
    node = Internal.getNodesFromName1(node, name)
    if node == []: return None
    node = node[0]
    val = Internal.getValue(node)
    return val

#==============================================================================
# Retourne le ref state si il existe dans t
#==============================================================================
def getRefState__(t):
    state = None
    ref = Internal.getNodesFromType2(t, 'ReferenceState_t')
    if len(ref) > 0:
        ref = ref[0]
        Density = Internal.getNodesFromName1(ref, 'Density')
        MomentumX = Internal.getNodesFromName1(ref, 'MomentumX')
        MomentumY = Internal.getNodesFromName1(ref, 'MomentumY')
        MomentumZ = Internal.getNodesFromName1(ref, 'MomentumZ')
        Energy = Internal.getNodesFromName1(ref, 'EnergyStagnationDensity')
        if len(Density) > 0: Density = Internal.getValue(Density[0])
        else: return None
        if len(MomentumX) > 0: MomentumX = Internal.getValue(MomentumX[0])
        else: return None
        if len(MomentumY) > 0: MomentumY = Internal.getValue(MomentumY[0])
        else: return None
        if len(MomentumZ) > 0: MomentumZ = Internal.getValue(MomentumZ[0])
        else: return None
        if len(Energy) > 0: Energy = Internal.getValue(Energy[0])
        else: return None
        state = [Density, MomentumX, MomentumY, MomentumZ, Energy]
        return state
    else: return None

#==============================================================================
# Retourne les data du BCDataSet
#==============================================================================
def getBCData__(node):
    out = None
    data = Internal.getNodesFromType(node, 'BCData_t')
    for d in data:
        array = Internal.getNodesFromType(d, 'DataArray_t')
        if array != []: out = array[0][1] ; break
    return out

#==============================================================================
# Retourne le type de la BC, sa fenetre et ses donnees eventuelles
# On regarde 1. Le BCDataSet, 2. le ref state de la zone, 3. le ref state
# de l'arbre
#==============================================================================
def getBC__(node, state):
    range = Internal.getNodesFromName(node, 'PointRange')
    win = Internal.range2Window(range[0][1])
    btype = Internal.getValue(node) # attention aux familles
    data = None
    statel = getBCData__(node)
    if statel is not None: data = statel
    else:
        statel = getRefState__(node)
        if statel is not None: data = statel
        else: data = state
    return (btype, win, data)

#==============================================================================
# IN: d: container
# IN: keys: les cles possibles
#==============================================================================
def checkKeys(d, keys):
    fail = False
    for i in d[2]:
      if keys.has_key(i[0]):
        k = keys[i[0]]
        val = Internal.getValue(i)
        if k == 0: # must be int
          if not isinstance(val, (int, long)):
            fail = True; print 'Error: Fast: keyword %s requires an int.'%i[0]
        elif k == 1: # must be float but we accept int
          if not isinstance(val, (float, int, long)):
            fail = True; print 'Error: Fast: keyword %s requires a float.'%i[0]
        elif k == 2: # must be any string
          if not isinstance(val, string):
            fail = True; print 'Error: Fast: keyword %s requires a string.'%i[0]
        elif k == 3: # Array/list of ints
          if not isinstance(val, numpy.ndarray):
            fail = True; print 'Error: Fast: keyword %s requires an array of ints.'%i[0]
        elif k == 4: # Array/list of floats
          if not isinstance(val, numpy.ndarray):
            fail = True; print 'Error: Fast: keyword %s requires an array of ints.'%i[0]
        elif isinstance(k, list): # list of given strings
          if not val in k:
            fail = True; print 'Error: Fast: keyword %s requires a string in %s.'%(i[0], str(k))
      else:
        print 'Error: Fast: keyword %s is INVALID.'%i[0]
        import difflib
        words = keys.keys()
        m = difflib.get_close_matches(i[0], words)
        print 'Error: Fast: do you mean %s?'%str(m)
        fail = True
    if fail:
      import sys; sys.exit(0) # much stronger os._exit(0)?
