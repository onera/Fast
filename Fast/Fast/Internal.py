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

MX_SYNCHRO = 1000
MX_SSZONE  = 10
MX_OMP_SIZE_INT = 250*OMP_NUM_THREADS

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

    if tc is not None: 

       #reordone les zones de tc par taille decroissante dans chaque base pour optim openmp
       bases_tc = Internal.getNodesFromType1(tc,'CGNSBase_t')
       for base_tc in bases_tc: 
          zones = Internal.getNodesFromType1(base_tc,'Zone_t')
          #calcul taille de la zone
          size_zone =[]
          for z in zones:
             dim = Internal.getZoneDim(z)
             if dim[0]=='Structured':
                if dim[3] == 1: kfic = 0
                else          : kfic = 2
                ndimdx = (dim[1]-4)*(dim[2]-4)*(dim[3]-kfic) 
             else: ndimdx = dim[2]
             size_zone.append(ndimdx)

          #Tri les zone par taille decroissante
          new_zones =[]
          for z in xrange(len(size_zone)):
              vmax    = max( size_zone )
              pos_max = size_zone.index( vmax )
              new_zones.append( zones[pos_max])
              size_zone.pop( pos_max) 
              zones.pop(     pos_max) 

          l = base_tc[2]
          orig = []
          for i in l:
	     if i[3] != 'Zone_t': orig.append(i)
          base_tc[2]=  orig+new_zones

       # reordone les zones de t pour garantir meme ordre entre t et tc
       bases    = Internal.getNodesFromType1(t, 'CGNSBase_t')
       for base in bases: 

           zones = Internal.getNodesFromType1(base,'Zone_t')
 
           #on cherche le matching entre base t et tc
           c  =0
           ctg=0
           for base_tc in bases_tc: 
              if base[0] == base_tc[0]: ctg=c
              c+=1  

           zones_tc  = Internal.getNodesFromType1(bases_tc[ctg], 'Zone_t')

           zname = []
           for zc in zones: zname.append( zc[0] )

           new_zones = []
           for z in zones_tc:
              pos = zname.index( z[0] )
              new_zones.append( zones[pos] )

           l = base[2]
           orig = []
           for i in l:
	        if i[3] != 'Zone_t': orig.append(i)
           base[2] = new_zones + orig

#==============================================================================
# Init/create primitive Variable 
#==============================================================================
def _createPrimVars(base, zone, omp_mode, rmConsVars=True, adjoint=False):
    import timeit
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
           C._initVars(zone, 'centers:ViscosityEddy', 1e-15)

    if (not (model == 'NSTurbulent' or model == 'nsspalart')): sa = False
    else: 
        sa = True
        if C.isNamePresent(zone, 'centers:TurbulentDistance') != 1: 
            raise ValueError("createPrimVars: TurbulentDistance field required at cell centers for RANS computations.")
    
    rgp = (adim[11]-1.)*adim[7]
    t0=timeit.default_timer()
    if C.isNamePresent(zone, 'centers:VelocityX'  ) != 1: P._computeVariables(zone, ['centers:VelocityX'  ], rgp=rgp)
    ''' 
    t1=timeit.default_timer()
    print "cout calcul VX= ", t1-t0, zone[0]
    t0=timeit.default_timer()
    '''
    if C.isNamePresent(zone, 'centers:VelocityY'  ) != 1: P._computeVariables(zone, ['centers:VelocityY'  ], rgp=rgp)
    '''
    t1=timeit.default_timer()
    print "cout calcul VY= ", t1-t0, zone[0]
    t0=timeit.default_timer()
    '''
    if C.isNamePresent(zone, 'centers:VelocityZ'  ) != 1: P._computeVariables(zone, ['centers:VelocityZ'  ], rgp=rgp)
    '''
    t1=timeit.default_timer()
    print "cout calcul VZ= ", t1-t0, zone[0]
    t0=timeit.default_timer()
    '''
    if C.isNamePresent(zone, 'centers:Temperature') != 1: P._computeVariables(zone, ['centers:Temperature'], rgp=rgp)
    '''
    t1=timeit.default_timer()
    print "cout calcul T = ", t1-t0, zone[0]
    t0=timeit.default_timer()
    '''

    if (sa and C.isNamePresent(zone, 'centers:TurbulentSANuTilde') != 1): 
        if C.isNamePresent(zone, 'centers:TurbulentSANuTildeDensity') != 1: C._initVars(zone, '{centers:TurbulentSANuTilde}= %20.16g/{centers:Density}'%adim[14])
        else: C._initVars(zone, '{centers:TurbulentSANuTilde}= {centers:TurbulentSANuTildeDensity}/{centers:Density}')
    '''
    t1=timeit.default_timer()
    print "cout calcul Nu = ", t1-t0, zone[0]
    t0=timeit.default_timer()
    '''

    if rmConsVars:
        C._rmVars(zone, 'centers:MomentumX')
        C._rmVars(zone, 'centers:MomentumY')
        C._rmVars(zone, 'centers:MomentumZ')
        C._rmVars(zone, 'centers:EnergyStagnationDensity')
        C._rmVars(zone, 'centers:TurbulentSANuTildeDensity')

    '''
    t1=timeit.default_timer()
    print "cout rmvars    = ", t1-t0, zone[0]
    t0=timeit.default_timer()
    '''

    #on test s'il existe 2 niveau en temps dans l'arbre pour appliquer la bonne formule de derivee temporelle a la premere iteration
    FIRST_IT = 1
    if C.isNamePresent(zone, 'centers:Density_M1')   != 1: C._cpVars(zone, 'centers:Density'  , zone, 'centers:Density_M1')  ; FIRST_IT=0
    '''
    t1=timeit.default_timer()
    print "cout calcul RoM1= ", t1-t0, zone[0]
    t0=timeit.default_timer()
    '''
    if C.isNamePresent(zone, 'centers:VelocityX_M1') != 1: C._cpVars(zone, 'centers:VelocityX', zone, 'centers:VelocityX_M1'); FIRST_IT=0
    '''
    t1=timeit.default_timer()
    print "cout calcul VXM1= ", t1-t0, zone[0]
    t0=timeit.default_timer()
    '''
    if C.isNamePresent(zone, 'centers:VelocityY_M1') != 1: C._cpVars(zone, 'centers:VelocityY', zone, 'centers:VelocityY_M1'); FIRST_IT=0
    '''
    t1=timeit.default_timer()
    print "cout calcul VYM1= ", t1-t0, zone[0]
    t0=timeit.default_timer()
    '''
    if C.isNamePresent(zone, 'centers:VelocityZ_M1') != 1: C._cpVars(zone, 'centers:VelocityZ', zone, 'centers:VelocityZ_M1'); FIRST_IT=0
    '''
    t1=timeit.default_timer()
    print "cout calcul VZM1= ", t1-t0, zone[0]
    t0=timeit.default_timer()
    '''
    if C.isNamePresent(zone, 'centers:Temperature_M1') != 1: C._cpVars(zone, 'centers:Temperature', zone, 'centers:Temperature_M1'); FIRST_IT=0
    '''
    t1=timeit.default_timer()
    print "cout calcul TM1 = ", t1-t0, zone[0]
    t0=timeit.default_timer()
    '''
    if (sa and  C.isNamePresent(zone, 'centers:TurbulentSANuTilde_M1') != 1): C._cpVars(zone, 'centers:TurbulentSANuTilde', zone, 'centers:TurbulentSANuTilde_M1'); FIRST_IT=0
    '''
    t1=timeit.default_timer()
    print "cout calcul NuM1= ", t1-t0, zone[0]
    t0=timeit.default_timer()
    '''

    if C.isNamePresent(zone, 'centers:Density_P1')   != 1: C._cpVars(zone, 'centers:Density'  , zone, 'centers:Density_P1')  
    '''
    t1=timeit.default_timer()
    print "cout calcul RoP1= ", t1-t0, zone[0]
    t0=timeit.default_timer()
    '''
    if C.isNamePresent(zone, 'centers:VelocityX_P1') != 1: C._cpVars(zone, 'centers:VelocityX', zone, 'centers:VelocityX_P1')
    '''
    t1=timeit.default_timer()
    print "cout calcul VXP1= ", t1-t0, zone[0]
    t0=timeit.default_timer()
    '''
    if C.isNamePresent(zone, 'centers:VelocityY_P1') != 1: C._cpVars(zone, 'centers:VelocityY', zone, 'centers:VelocityY_P1')
    '''
    t1=timeit.default_timer()
    print "cout calcul VYP1= ", t1-t0, zone[0]
    t0=timeit.default_timer()
    '''
    if C.isNamePresent(zone, 'centers:VelocityZ_P1') != 1: C._cpVars(zone, 'centers:VelocityZ', zone, 'centers:VelocityZ_P1')
    '''
    t1=timeit.default_timer()
    print "cout calcul VZP1= ", t1-t0, zone[0]
    t0=timeit.default_timer()
    '''
    if C.isNamePresent(zone, 'centers:Temperature_P1') != 1: C._cpVars(zone, 'centers:Temperature', zone, 'centers:Temperature_P1')
    '''
    t1=timeit.default_timer()
    print "cout calcul TP1 = ", t1-t0, zone[0]
    t0=timeit.default_timer()
    '''
    if (sa and  C.isNamePresent(zone, 'centers:TurbulentSANuTilde_P1') != 1): C._cpVars(zone, 'centers:TurbulentSANuTilde', zone, 'centers:TurbulentSANuTilde_P1')
    '''
    t1=timeit.default_timer()
    print "cout calcul NuP1= ", t1-t0, zone[0]
    t0=timeit.default_timer()
    '''

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

    DES_debug = "none"
    a = Internal.getNodeFromName2(zone, 'DES_debug')
    if a is not None: DES_debug = Internal.getValue(a)
    des = "none"
    a = Internal.getNodeFromName2(zone, 'DES')
    if a is not None: des = Internal.getValue(a)
    else:
        a = Internal.getNodeFromName2(zone, 'des')
        if a is not None: des = Internal.getValue(a)
        else:
            a = Internal.getNodeFromName2(base, 'DES')
            if a is not None: des = Internal.getValue(a)
            else:
                a = Internal.getNodeFromName2(base, 'des')
                if a is not None: des = Internal.getValue(a)
    if des != 'none':
       if   des == "ZDES1"   or des == "zdes1":   ides = 2
       elif des == "ZDES1_w" or des == "zdes1_w": ides = 3
       elif des == "ZDES2"   or des == "zdes2":   ides = 4
       elif des == "ZDES2_w" or des == "zdes2_w": ides = 5
       elif des == "ZDES3"   or des == "zdes3":   ides = 6
       elif des == "ZDES3_w" or des == "zdes3_w": ides = 7
       if ides >= 2 and DES_debug != "none": 
          if C.isNamePresent(zone, 'centers:delta')!= 1: C._initVars(zone, 'centers:delta', 0.)
       if (ides == 4 or ides == 5) and DES_debug != "none":  
          if C.isNamePresent(zone, 'centers:fd')   != 1: C._initVars(zone, 'centers:fd', 0.)

       if (ides == 6 or ides == 7) and C.isNamePresent(z, 'centers:zgris') != 1:
           raise ValueError("createPrimVars: 'In ZDES mode 3, zgris field required at cell centers.")


    extract_res = 0
    a = Internal.getNodeFromName2(zone, 'extract_res')
    if a is not None: extract_res = Internal.getValue(a)
    if extract_res == 1:
       if C.isNamePresent(zone, 'centers:Res_Density') != 1:   C._initVars(zone, 'centers:Res_Density', 0.)
       if C.isNamePresent(zone, 'centers:Res_MomentumX') != 1:   C._initVars(zone, 'centers:Res_MomentumX', 0.)
       if C.isNamePresent(zone, 'centers:Res_MomentumY') != 1:   C._initVars(zone, 'centers:Res_MomentumY', 0.)
       if C.isNamePresent(zone, 'centers:Res_MomentumZ') != 1:   C._initVars(zone, 'centers:Res_MomentumZ', 0.)
       if C.isNamePresent(zone, 'centers:Res_EnergyStagnationDensity') != 1: C._initVars(zone, 'centers:Res_EnergyStagnationDensity', 0.)
       if (sa and  C.isNamePresent(zone, 'centers:Res_TurbulentSANuTildeDensity') != 1): C._initVars(zone, 'centers:Res_TurbulentSANuTildeDensity', 0.)
    if extract_res == 2:
       C._initVars(zone, 'centers:drodm0_1', 1e-15)
       C._initVars(zone, 'centers:drodm0_2', 1e-15)
       C._initVars(zone, 'centers:drodm0_3', 1e-15)
       C._initVars(zone, 'centers:drodm0_4', 1e-15)
       C._initVars(zone, 'centers:drodm0_5', 1e-15)
       C._initVars(zone, 'centers:drodmN_1', 1e-15)
       C._initVars(zone, 'centers:drodmN_2', 1e-15)
       C._initVars(zone, 'centers:drodmN_3', 1e-15)
       C._initVars(zone, 'centers:drodmN_4', 1e-15)
       C._initVars(zone, 'centers:drodmN_5', 1e-15)
       if sa:
          C._initVars(zone, 'centers:drodm0_6', 1e-15)
          C._initVars(zone, 'centers:drodmN_6', 1e-15)

    #KRYLOV SUBSPACE VECTOR ALLOCATION
    implicit_solver = 'none'
    a = Internal.getNodeFromName2(zone, 'implicit_solver')
    if a is not None: implicit_solver = Internal.getValue(a)
    if (implicit_solver == 'gmres'):
        #on nettoie l'arbre des vecteur de krylov residuel.Evite la phase de compactage
        flowsol = Internal.getNodeFromName1(zone, 'FlowSolution#Centers')
        if flowsol is not None:
            vars    = Internal.getNodesFromType1(flowsol, 'DataArray_t')
            for var in vars:
                if 'Kry' in var[0]:
                     C._rmVars(t, 'centers:'+var[0])

        #dimesionnement tableau krylov
        nbr_krylov = 20
        a = Internal.getNodeFromName2(zone, 'nb_krylov')
        if a is not None: nbr_krylov = Internal.getValue(a)
        dens = Internal.getNodeFromName2(zone, 'Density')
        ndimdx = dens[1].size
        size   = nbr_krylov*ndimdx*5
        if sa: size = nbr_krylov*ndimdx*6
        kry = numpy.empty( size , dtype=numpy.float64)

        sol  = Internal.getNodeFromName2(zone, 'FlowSolution#Centers')
        vars = ['Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity']
        if sa : vars.append('TurbulentSANuTildeDensity')
        shift =0
        VarNew = Internal.copyTree(dens)
        for Vector in range(nbr_krylov):
            for var in vars:
               name='Kry_' + str(Vector) + '_'+var
               tmp = kry[ shift + 0: shift + ndimdx ]
               Internal.createChild(sol , name, 'DataArray_t', value=tmp, children=VarNew[2])
               shift += ndimdx

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

    # init termes sources
    source = 0
    a = Internal.getNodeFromName2(zone, 'source')
    if a is not None: source = Internal.getValue(a)
    if source == 1:
       if C.isNamePresent(zone, 'centers:Density_src')   != 1:   C._initVars(zone, 'centers:Density_src', 0.)
       if C.isNamePresent(zone, 'centers:MomentumX_src') != 1:   C._initVars(zone, 'centers:MomentumX_src', 0.)
       if C.isNamePresent(zone, 'centers:MomentumY_src') != 1:   C._initVars(zone, 'centers:MomentumY_src', 0.)
       if C.isNamePresent(zone, 'centers:MomentumZ_src') != 1:   C._initVars(zone, 'centers:MomentumZ_src', 0.)
       if C.isNamePresent(zone, 'centers:EnergyStagnationDensity_src') != 1: C._initVars(zone, 'centers:EnergyStagnationDensity_src', 0.)
       if (sa and C.isNamePresent(zone, 'centers:TurbulentSANuTildeDensity_src') != 1): C._initVars(zone, 'centers:TurbulentSANuTildeDensity_src', 0.)



    return sa, FIRST_IT

#==============================================================================
# create work arrays
#==============================================================================
def createWorkArrays__(zones, dtloc, FIRST_IT):
    ndimt = 0 ; ndimcoe = 0 ; ndimwig = 0 ; c = 0;
    global  MX_OMP_SIZE_INT
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

        param_int = Internal.getNodeFromName2(z, 'Parameter_int')
        if param_int is not None:
           shiftvar  = param_int[1][66]
        else:
           shiftvar  = 100
           print 'shiftVar=', shiftvar
           print 'create workarray'
           print 'Danger sur optimisation cache shiftVar'      
           print 'Danger '      
           print 'Danger '      

        # dim from getZoneDim
        dims = Internal.getZoneDim(z)
        if dims[0]=='Structured' : nijk = (dims[1]-1)*(dims[2]-1)*(dims[3]-1)+shiftvar
        else                     : nijk = dims[2]+shiftvar
        ndimt   +=     neq*nijk       # surdimensionne ici
        ndimcoe += neq_coe*nijk       # surdimensionne ici
        ndimwig +=        3*nijk      # surdimensionne ici
        c       += 1

    mx_thread   = OMP_NUM_THREADS       # surdimensionne : doit etre = a OMP_NUM_THREADS
    verrou      = MX_SSZONE*c*MX_SYNCHRO*mx_thread

    timerOmp = numpy.zeros(  mx_thread*2*dtloc[0] + mx_thread*len(zones)+1, dtype=numpy.float64)

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

    lok      = numpy.zeros(verrou     , dtype=numpy.int32  )
    iskip_lu = numpy.empty(dtloc[0]*2 , dtype=numpy.int32  )   # (dtloc[0] = nitmax

    #tableau pour ibc from C function
    param_int_ibc  = numpy.empty((2), numpy.int32)         
    param_real_ibc = numpy.empty((5), numpy.float64)


    hook = {}
    hook['wiggle']         = wig
    hook['coe']            = coe
    hook['rhs']            = drodm
    hook['verrou_omp']     = lok
    hook['skip_lu']        = iskip_lu
    hook['dtloc']          = dtloc
    hook['lssiter_loc']    = lssiter_loc
    hook['MX_SSZONE']      = MX_SSZONE
    hook['MX_SYNCHRO']     = MX_SYNCHRO
    hook['MX_OMP_SIZE_INT']= MX_OMP_SIZE_INT
    hook['TIMER_OMP']      = timerOmp
    hook['FIRST_IT']       = FIRST_IT
    hook['neq_max']        = neq_max
    hook['param_int_ibc']  = param_int_ibc
    hook['param_real_ibc'] = param_real_ibc
    hook['param_int_tc']   = None
    hook['param_real_tc']  = None
    hook['mpi']            = 0     
    hook['linelets_int']   = None
    hook['linelets_real']  = None
    return hook

#==============================================================================
# Retourne un tag suivant bcname
#==============================================================================
def tagBC(bcname):
  if   bcname == "BCExtrapolate"    :       tag = 0;
  elif bcname == "BCDegenerateLine":        tag = 0;
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
  elif bcname == "BCPeriodic":              tag =15;
  elif bcname == "BCOutpres":               tag =16;
  elif bcname == "BCInj1":                  tag =17;
  elif bcname == "BCRacinf":                tag =18;
  else:
    tag = -1
    print "Warning: unknown BC type %s."%bcname
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
def switchPointers__(zones, case, order=3):
    
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

            model  = Internal.getNodeFromName1(own, 'model')
            model  = Internal.getValue(model)

            if case ==1:
              # sauvegarde M1
              ta = caM1[1]; tb = cbM1[1]; tc = ccM1[1]; td = cdM1[1]; te = ceM1[1]

              # M1 <- current
              caM1[1] = ca[1]; cbM1[1] = cb[1]; ccM1[1] = cc[1]; cdM1[1] = cd[1]; ceM1[1] = ce[1]

              # current <- P1
              ca[1] = caP1[1]; cb[1] = cbP1[1]; cc[1] = ccP1[1]; cd[1] = cdP1[1]; ce[1] = ceP1[1]

              # P1 <- temp
              caP1[1] = ta; cbP1[1] = tb; ccP1[1] = tc; cdP1[1] = td; ceP1[1] = te

              if model == 'nsspalart' or model=='NSTurbulent': 
                cf   =  Internal.getNodeFromName1(sol, 'TurbulentSANuTilde')
                cfM1 =  Internal.getNodeFromName1(sol, 'TurbulentSANuTilde_M1')
                cfP1 =  Internal.getNodeFromName1(sol, 'TurbulentSANuTilde_P1')
                tf   = cfM1[1]; cfM1[1] = cf[1]; cf[1] = cfP1[1]; cfP1[1] = tf

            elif case ==2:
              # sauvegarde  P1
              ta = caP1[1]; tb = cbP1[1]; tc = ccP1[1]; td = cdP1[1]; te = ceP1[1]

              # P1 <-  current    
              caP1[1] = ca[1]; cbP1[1] = cb[1]; ccP1[1] = cc[1]; cdP1[1] = cd[1]; ceP1[1] = ce[1]

              # current <- M1
              ca[1] = caM1[1]; cb[1] = cbM1[1]; cc[1] = ccM1[1]; cd[1] = cdM1[1]; ce[1] = ceM1[1]

              # M1 <- temp
              caM1[1] = ta; cbM1[1] = tb; ccM1[1] = tc; cdM1[1] = td; ceM1[1] = te

              if model == 'nsspalart' or model=='NSTurbulent': 
                cf   =  Internal.getNodeFromName1(sol, 'TurbulentSANuTilde')
                cfM1 =  Internal.getNodeFromName1(sol, 'TurbulentSANuTilde_M1')
                cfP1 =  Internal.getNodeFromName1(sol, 'TurbulentSANuTilde_P1')
                tf   = cfP1[1]; cfP1[1] = cf[1];  cf[1]= cfM1[1]; cfM1[1] = tf

    elif order == 2:
        print 'adapter switch poinyer'
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

            if model == 'nsspalart' or model == 'NSTurbulent': 
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
            fail = True; print 'Error: Fast: keyword %s requires an array of floats.'%i[0]
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
