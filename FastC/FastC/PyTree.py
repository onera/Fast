"""Services for FAST solvers.
"""
import numpy
import os
from . import fastc

FIRST_IT = 0
HOOK     = None

# transfered variables
varsN    = ['Density']
varsP    = ['Density_P1']
varsM    = ['Density_M1']
varsPLBM = ['Q1']
varsMLBM = ['Q1_M1']
varsMacro= ['Density']

try:
    import Converter.PyTree as C
    import Converter.Internal as Internal
    import Post.PyTree as P
    import math
except:
    raise ImportError("Fast.Internal: requires Converter and Post module.")

try: range = xrange
except: pass

try:
    OMP_NUM_THREADS = os.environ['OMP_NUM_THREADS']
    OMP_NUM_THREADS = int(OMP_NUM_THREADS)
except: OMP_NUM_THREADS = 1

import Converter.Mpi as Cmpi

MX_SYNCHRO = 1000
MX_SSZONE  = 10
MX_OMP_SIZE_INT = 250*OMP_NUM_THREADS

#==============================================================================
# Met un dictionnaire de numerics dans une/des zones
#==============================================================================
def _setNum2Zones(a, num):
    zones = Internal.getNodesFromType2(a, 'Zone_t')
    for z in zones:
        cont = Internal.createUniqueChild(z, '.Solver#define', 
                                          'UserDefinedData_t')
        for k in num: # some checks?
            Internal.createUniqueChild(cont, k, 'DataArray_t', num[k])
    return None

#==============================================================================
def _setNum2Base(a, num):
    bases = Internal.getNodesFromType1(a, 'CGNSBase_t')
    for b in bases:
        cont = Internal.createUniqueChild(b, '.Solver#define', 
                                          'UserDefinedData_t')
        for k in num: # some checks?
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
          for z in range(len(size_zone)):
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
# Converti les variables conservatives de l'arbre en variables primitives
# compactees
# IN: t: arbre devant contenir les variables conservatives aux centres
# IN: adim: etat de reference, on utilise uniquement cvInf pour la temperature.
# IN: rmConsVars: if True, remove the conservative variables
#==============================================================================
def createPrimVars(t, omp_mode, rmConsVars=True, Adjoint=False):
    tp = Internal.copyRef(t)
    first_iter , vars_zones = _createPrimVars(tp, omp_mode, rmConsVars, Adjoint)
    return tp, first_iter, vars_zones

#==============================================================================
def _createPrimVars(t, omp_mode, rmConsVars=True, Adjoint=False):
    vars_zones=[]
    bases = Internal.getNodesFromType1(t, 'CGNSBase_t')
    for b in bases:  
        zones = Internal.getNodesFromType1(b, 'Zone_t')
        count = -1
        if omp_mode == 1: count = 0
        for z in zones:
            if omp_mode == 1: count += 1

            #Sauvegarde noeud specifique car delete par init!!! Christophe??
            FA_intext=  Internal.getNodeFromPath(z, 'NFaceElements/IntExt')
            if FA_intext is not None:
               FA_indx  =  Internal.getNodeFromPath(z, 'NFaceElements/ElementIndex')
               NG_intext=  Internal.getNodeFromPath(z, 'NGonElements/IntExt')
               NG_indx  =  Internal.getNodeFromPath(z, 'NGonElements/FaceIndex')
               NG_pe    =  Internal.getNodeFromPath(z, 'NGonElements/ParentElements')

            sa, lbmflag, neq_lbm, FIRST_IT  = _createVarsFast(b, z, omp_mode, rmConsVars,Adjoint)
            
            if FA_intext is not None:
               Internal.getNodeFromPath(z, 'NFaceElements')[2] +=[FA_indx, FA_intext] 
               Internal.getNodeFromPath(z, 'NGonElements' )[2] +=[NG_pe, NG_indx, NG_intext] 
            
            #HOOK['FIRST_IT']= FIRST_IT
            source = 0
            a = Internal.getNodeFromName2(z,'Density_src')
            if a is not None: source = 1
            sfd = 0
            a = Internal.getNodeFromName2(z, 'sfd')
            if a is not None: sfd = Internal.getValue(a)
            extract_res = 0
            a = Internal.getNodeFromName2(z, 'extract_res')
            if a is not None: extract_res = Internal.getValue(a)
            motion = None 
            a = Internal.getNodeFromName2(z, 'motion')
            if a is not None: motion = Internal.getValue(a)

            vars_p=['Density', 'VelocityX', 'VelocityY','VelocityZ', 'Temperature']
            vars_c=['Density', 'MomentumX', 'MomentumY','MomentumZ', 'EnergyStagnationDensity']
            if sa: 
               vars_p.append('TurbulentSANuTilde')
               vars_c.append('TurbulentSANuTildeDensity')

            vars=[]
            timelevel = ['', '_M1', '_P1']
            if lbmflag: 
               fields2compact =[]
               for i in range(1,neq_lbm+1):
                 fields2compact.append('centers:Q'+str(i))
               vars.append(fields2compact)
               fields2compact =[]
               for i in range(1,neq_lbm+1):
                 fields2compact.append('centers:Q'+str(i)+'_M1')
               vars.append(fields2compact)
               timelevel = ['']
            for level in timelevel:  #champs primitives
               fields2compact =[]
               for v in vars_p:
                 fields2compact.append('centers:'+v+level)
               vars.append(fields2compact)
            if source == 1:                   #terme source volumique
               fields2compact =[]
               for v in vars_c:
                 fields2compact.append('centers:'+v+'_src')
               vars.append(fields2compact)
            if motion == 'deformation':     #ale deformable
               fields2compact =[]
               for v in ['VelocityX','VelocityY','VelocityZ']:
                 fields2compact.append('Motion:'+v)
               vars.append(fields2compact)
            if sfd == 1:                       #sfd
               fields2compact =[]
               for v in vars_p:
                 fields2compact.append('centers:'+v+'_f')
               vars.append(fields2compact)
            if extract_res == 1:
               fields2compact =[]
               for v in vars_c:
                 fields2compact.append('centers:Res_'+v)
               vars.append(fields2compact)
    
            # on compacte les variables "visqueuse"
            loc            ='centers:'
            fields         = [loc+'ViscosityEddy',loc+'TurbulentDistance', loc+'zgris']
            for sufix in ['0','N']:
               for i in range(1,7):
                  fields.append(loc+'drodm'+sufix+'_'+str(i))

            fields2compact = []
            for field in fields: 
               if C.isNamePresent(z, field) == 1: fields2compact.append(field)
            if len(fields2compact) != 0: vars.append(fields2compact)

            # on compacte les variables SA debug
            fields         = [loc+'delta',loc+'fd']
            fields2compact = []
            for field in fields: 
               if C.isNamePresent(z, field) == 1: fields2compact.append(field)
            if len(fields2compact) != 0: vars.append(fields2compact)

            # on compacte CellN
            fields2compact = []
            if  C.isNamePresent(z, loc+'cellN') == 1: fields2compact.append(loc+'cellN')
            if len(fields2compact) != 0: vars.append(fields2compact)
            
            # on compacte CellN_IBC
            fields2compact = []
            if  C.isNamePresent(z, loc+'cellN_IBC') == 1: fields2compact.append(loc+'cellN_IBC')
            if len(fields2compact) != 0: vars.append(fields2compact)

            #  adjoint 
            if  C.isNamePresent(z, 'centers:dpCLp_dpDensity') == 1: 
               fields =['dpCDp_dp', 'dpCLp_dp', 'rhsIterAdjCLp_R', 'rhsIterAdjCDp_R','AdjCLp_R', 'AdjCDp_R','incAdj_R']

               for field in fields:
                  fields2compact =[]
                  for v in vars_c:
                      fields2compact.append('centers:'+field+v)
                  vars.append(fields2compact)

               vars.append(['dpCLp_dpX','dpCLp_dpY','dpCLp_dpZ','dpCDp_dpX','dpCDp_dpY','dpCDp_dpZ']) 
 
            ##on stocke zone et variable car fonction compact specific fastS ou FastP: ne peut peut etre appelee ici
            vars_zones.append( [z, vars] )
            
    return FIRST_IT, vars_zones

#==============================================================================
# Init/create primitive Variable 
#==============================================================================
def _createVarsFast(base, zone, omp_mode, rmConsVars=True, adjoint=False):
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

    if not (model == 'Euler' or model == 'euler'):
        if C.isNamePresent(zone, 'centers:ViscosityEddy') != 1: 
           C._initVars(zone, 'centers:ViscosityEddy', 1e-15)

    if not (model == 'NSTurbulent' or model == 'nsspalart'): sa = False
    else: 
        sa = True
        if C.isNamePresent(zone, 'centers:TurbulentDistance') != 1: 
            raise ValueError("createPrimVars: TurbulentDistance field required at cell centers for RANS computations.")
    
    if not (model == 'LBMLaminar' or model == 'lbmlaminar'): lbm = False
    else: lbm = True

    rgp = (adim[11]-1.)*adim[7]
    t0 = timeit.default_timer()
    if C.isNamePresent(zone, 'centers:VelocityX'  ) != 1: P._computeVariables(zone, ['centers:VelocityX'  ], rgp=rgp)

    if C.isNamePresent(zone, 'centers:VelocityY'  ) != 1: P._computeVariables(zone, ['centers:VelocityY'  ], rgp=rgp)
    if C.isNamePresent(zone, 'centers:VelocityZ'  ) != 1: P._computeVariables(zone, ['centers:VelocityZ'  ], rgp=rgp)
    if C.isNamePresent(zone, 'centers:Temperature') != 1: P._computeVariables(zone, ['centers:Temperature'], rgp=rgp)

    if (sa and C.isNamePresent(zone, 'centers:TurbulentSANuTilde') != 1): 
        if C.isNamePresent(zone, 'centers:TurbulentSANuTildeDensity') != 1: C._initVars(zone, '{centers:TurbulentSANuTilde}= %20.16g/{centers:Density}'%adim[14])
        else: C._initVars(zone, '{centers:TurbulentSANuTilde}= {centers:TurbulentSANuTildeDensity}/{centers:Density}')

    if rmConsVars:
        C._rmVars(zone, 'centers:MomentumX')
        C._rmVars(zone, 'centers:MomentumY')
        C._rmVars(zone, 'centers:MomentumZ')
        C._rmVars(zone, 'centers:EnergyStagnationDensity')
        C._rmVars(zone, 'centers:TurbulentSANuTildeDensity')


    FIRST_IT = 1
    if model != 'LBMLaminar':
       neq_lbm = 0
       #on test s'il existe 2 niveaux en temps dans l'arbre pour appliquer la bonne formule de derivee temporelle a la premere iteration
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
    else:
       neq_lbm = Internal.getNodeFromName2(zone, 'Parameter_int')[1][86]

       #print('Internal fast: neq_LBM=', neq_lbm)
       for i in range(1,neq_lbm+1):
           if C.isNamePresent(zone, 'centers:Q'+str(i)) != 1: C._initVars(zone, 'centers:Q'+str(i), 0.)
           if C.isNamePresent(zone, 'centers:Q'+str(i)+'_M1') != 1: C._initVars(zone, 'centers:Q'+str(i)+'_M1', 0.)

    # gestion ale maillage deformable
    ale=0
    b = Internal.getNodeFromName1(zone, '.Solver#define')
    if b is not None:
       a = Internal.getNodeFromName1(b, 'motion')
       if a is not None: ale = Internal.getValue(a)
       if ale =='deformation':
          Motion = Internal.getNodeFromName1(zone, 'Motion')
          vx = None
          if Motion is not None: vx = Internal.getNodeFromName1(Motion, 'VelocityX')
          else:
              Internal.createNode('Motion', 'ArbitraryGridMotion_t', value=None, children=[], parent=zone) 
              Motion   = Internal.getNodeFromType1(zone, 'ArbitraryGridMotion_t')
              Internal.createNode('ArbitraryGridMotion', 'ArbitraryGridMotionType_t', value='DeformingGrid', children=[], parent=Motion) 

          if vx is None: 
              a             = Internal.getNodeFromName1(zone, 'GridCoordinates')
              size_velocity = numpy.size( Internal.getNodeFromName1(a, 'CoordinateX')[1])
              datax = numpy.zeros(size_velocity, numpy.float64)
              datay = numpy.zeros(size_velocity, numpy.float64)
              dataz = numpy.zeros(size_velocity, numpy.float64)
              Internal.createUniqueChild(Motion, 'VelocityX', 'DataArray_t', datax)
              Internal.createUniqueChild(Motion, 'VelocityY', 'DataArray_t', datay)
              Internal.createUniqueChild(Motion, 'VelocityZ', 'DataArray_t', dataz)

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

       if (ides == 6 or ides == 7) and C.isNamePresent(zone, 'centers:zgris') != 1:
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
    if implicit_solver == 'gmres':
        #on nettoie l'arbre des vecteur de krylov residuel.Evite la phase de compactage
        flowsol = Internal.getNodeFromName1(zone, 'FlowSolution#Centers')
        if flowsol is not None:
            vars    = Internal.getNodesFromType1(flowsol, 'DataArray_t')
            for var in vars:
                if 'Kry' in var[0]:
                     C._rmVars(a, 'centers:'+var[0])

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
    if adjoint: 
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

    return sa, lbm, neq_lbm, FIRST_IT

#==============================================================================
# Construit les donnees compactees pour traiter les verrou et plage omp
#==============================================================================
def _build_omp(t):
    # Data for each Base
    bases = Internal.getNodesFromType1(t, 'CGNSBase_t')

    #dimensionnememnt tableau
    size_int       = HOOK['MX_OMP_SIZE_INT']
    for b in bases:
        zones = Internal.getZones(b)
        for z in zones:
            
            o    = Internal.getNodeFromName1( z, '.Solver#ownData')
            dims = Internal.getZoneDim(z)

            #on concatene les donnes omp dans param_int
            param_int = Internal.getNodeFromName1( o, 'Parameter_int')
            size = numpy.shape(param_int[1])
            c = 1
            for s in size: c=c*s

            size_scheduler = 0
            if dims[3] == 'NGON':
               CellScheduler = Internal.getNodeFromName1( z, 'CellScheduler')[1]
               size_scheduler += numpy.size(CellScheduler)
               Nbr_BlkIntra = Internal.getNodeFromName1( z, 'Nbr_BlkIntra')[1]
               size_scheduler += numpy.size(Nbr_BlkIntra)

            #print("SIZE int", size_int,"c=", c,"size_sched=",  size_scheduler)
            datap = numpy.zeros(size_int + c + size_scheduler, numpy.int32)

            datap[0:c]   = param_int[1][0:c]
            datap[69]    = c                   # debut tableau omp dans param_int
            datap[87]    = c + size_int        # debut cellscheduler omp dans param_int
            param_int[1] = datap

            if dims[3] == 'NGON':
               deb = param_int[1][87] 
               size = numpy.size(Nbr_BlkIntra)
               for l in range(size):
                 param_int[1][ deb +l ]=Nbr_BlkIntra [l]
               Nbr_BlkIntra = param_int[1][ deb:deb+size ]

               deb  += size
               shap          = numpy.shape(CellScheduler)
               size_scheduler= numpy.size(CellScheduler)
               c = 0
               for l in range(shap[3]):
                 for k in range(shap[2]):
                   for j in range(shap[1]):
                     for i in range(shap[0]):
                       #print("deb+c", deb+c, c, l,k,j,i,param_int[1][87] )
                       param_int[1][ deb +c ]= CellScheduler[ i, j, k, l]
                       c +=1
               CellScheduler = param_int[1][ deb:deb+size_scheduler ]

    return None
#==============================================================================
# Construit les datas possedees par FastP
#==============================================================================
def _buildOwnData(t, Padding):
    # Data for each Base
    bases = Internal.getNodesFromType1(t, 'CGNSBase_t')

    #init time et No iteration
    it0 =0; temps=0.; timelevel_motion= 0; timelevel_target= 0
    first = Internal.getNodeFromName1(t, 'Iteration')
    if first is not None: it0 = Internal.getValue(first)
    first = Internal.getNodeFromName1(t, 'Time')
    if first is not None: temps = Internal.getValue(first)
    first = Internal.getNodeFromName1(t, 'TimeLevelMotion')
    if first is not None: timelevel_motion = Internal.getValue(first)
    first = Internal.getNodeFromName1(t, 'TimeLevelTarget')
    if first is not None: timelevel_target = Internal.getValue(first)

     
    levelg=[]; leveld=[]; val=1; i=0
    veclevel = []; mod="";posg=[] 
    for b in bases:

           zones = Internal.getNodesFromType1(b, 'Zone_t')
           for z in zones:
               a = Internal.getNodeFromName2(z, 'GoverningEquations')
               if a is not None: mod = Internal.getValue(a)
               if mod == "nsspalart" or mod == "NSTurbulent": neq=6
               else: neq = 5
        
           ### Recuperation du niveau en temps de chaque zone si on est en dtloc instationnaire ###          
           d = Internal.getNodeFromName1(b, '.Solver#define')
           if d is not None:
               a = Internal.getNodeFromName1(d, 'temporal_scheme')
               if a is not None:
                   temporal_scheme = Internal.getValue(a)
                   if temporal_scheme == 'explicit_local':
                       for z in zones:
                          dim = Internal.getZoneDim(z)
                          i = dim[1]//2
                          j = dim[2]//2
                          k = dim[3]//2
                          level = C.getValue( z, 'centers:niveaux_temps', (i,j,k))
                          veclevel.append(int(level))
                   else:
                       for z in zones:
                          veclevel.append(1)                        

           else:
               for z in zones:
                   veclevel.append(1)              
               
                          
    maxlevel = max(veclevel)

    #print('veclevel= ', veclevel)

    levelg = numpy.roll(veclevel,1)
    leveld = numpy.roll(veclevel,-1)

    # Available keys for bases and zones
    # 0: requires and int, 1: requires a float, 2: requires any string, 
    # 3: requires array/list of ints, 4: requires array/list of floats,
    # []: requires given strings
    keys4Base = {
    'temporal_scheme':['explicit', 'implicit', 'implicit_local','explicit_local'],
    'ss_iteration':0,
    'rk':0, 
    'modulo_verif':0,
    'exp_local':0,
    'time_begin_ale':1,
    'omp_mode':0,
    'explicit_local_type':0    ## 0:explicit local non conservatif,  1:explicit local conservatif (correction du bilan de flux aux interface)  
    }
    keys4Zone = {
    'scheme':['ausmpred', 'senseur', 'roe_min', 'roe', 'roe_nul', 'roe_kap'],
    'implicit_solver':['lussor', 'gmres'],
    'nb_relax':0,
    'nb_krylov':0,
    'nb_restart':0,
    'motion':['none', 'rigid', 'deformation'],
    'rotation':4,
    'time_step':1,
    'io_thread':0,
    'sgsmodel': ['smsm','Miles'],
    'wallmodel': ['musker','power'],
    'wallmodel_sample': 0,
    'ransmodel': ['SA','SA_comp','SA_diff'],
    'cache_blocking_I':0,
    'cache_blocking_J':0,
    'cache_blocking_K':0,
    'shiftvar':0,
    'time_step_nature':['local', 'global'],
    'ssdom_IJK':3,
    'lu_match':1,
    'epsi_newton':1,
    'epsi_linear':1,
    'inj1_newton_tol':1,
    'inj1_newton_nit':0,
    'cfl':1, 
    'niveaux_temps':0, 
    'psiroe':1, 
    'prandtltb':1, 
    'sfd':0, 
    'sfd_chi':1, 
    'sfd_delta':1, 
    'sfd_init_iter':0, 
    'slope':["o1", "o3", "minmod"],
    'DES':["zdes1", "zdes1_w", "zdes2", "zdes2_w", "zdes3", "zdes3_w"],
    'snear': 1, # ignored
    'DES_debug':['none','active'],
    'extract_res':0,
    'IBC':3,
    'source':0,
    'Cups':4,
    'senseurType':0,
    'ratiom':1
    }


    for b in bases:
        # Base define data
        d = Internal.getNodeFromName1(b, '.Solver#define')

        # - defaults -
        temporal_scheme = "implicit"
        ss_iteration    = 30
        rk              = 3
        modulo_verif    = 200
        restart_fields  = 1
        exploc          = 0
        t_init_ale      = temps
        timelevel_prfile= 0
        exploctype      = 0

        if d is not None:
            checkKeys(d, keys4Base)
            a = Internal.getNodeFromName1(d, 'temporal_scheme')
            if a is not None: temporal_scheme = Internal.getValue(a)
            a = Internal.getNodeFromName1(d, 'ss_iteration')
            if a is not None: ss_iteration = Internal.getValue(a)
            if temporal_scheme == "implicit_local": modulo_verif = 7
            a = Internal.getNodeFromName1(d, 'modulo_verif')
            if a is not None: modulo_verif = Internal.getValue(a)
            a = Internal.getNodeFromName1(d, 'restart_fields')
            if a is not None: restart_fields = Internal.getValue(a)
            a = Internal.getNodeFromName1(d, 'rk')
            if a is not None: rk = Internal.getValue(a)
            if temporal_scheme == "implicit" or temporal_scheme =="implicit_local": rk=3
            a = Internal.getNodeFromName1(d, 'exp_local')
            if a is not None: exploc = Internal.getValue(a)
            if temporal_scheme == "implicit" or temporal_scheme =="implicit_local": exploc=0
            a = Internal.getNodeFromName1(d, 'explicit_local_type')         
            if a is not None: exploctype = Internal.getValue(a)
            a = Internal.getNodeFromName1(d, 'time_begin_ale')
            if a is not None: t_init_ale = Internal.getValue(a)


           
        # Base ownData (generated)
        o = Internal.createUniqueChild(b, '.Solver#ownData', 
                                       'UserDefinedData_t')
        if temporal_scheme == "explicit": nssiter = 3
        elif temporal_scheme == "explicit_lbm": nssiter = 1
        elif temporal_scheme == "implicit": nssiter = ss_iteration+1
        elif temporal_scheme == "implicit_local": nssiter = ss_iteration+1
        elif temporal_scheme == "explicit_local": #explicit local instationnaire ordre 3
            nssiter = 4*maxlevel 
            rk = 3
            exploc = 1
        else: print('Warning: FastS: invalid value %s for key temporal_scheme.'%temporal_scheme)
        try: ss_iteration = int(ss_iteration)
        except: print('Warning: Fast: invalid value %s for key ss_iteration.'%ss_iteration)
        try: modulo_verif = int(modulo_verif)
        except: print('Warning: Fast: invalid value %s for key modulo_verif.'%modulo_verif)
        if rk == 1 and exploc==0 and temporal_scheme == "explicit": nssiter = 1 # explicit global
        if rk == 2 and exploc==0 and temporal_scheme == "explicit": nssiter = 2 # explicit global
        if rk == 3 and exploc==0 and temporal_scheme == "explicit": nssiter = 3 # explicit global

        

        dtdim = nssiter + 7
        datap = numpy.empty((dtdim), numpy.int32) 
        datap[0] = nssiter 
        datap[1] = modulo_verif
        datap[2] = restart_fields-1
        datap[3] = timelevel_motion
        datap[4] = timelevel_target 
        datap[5] = timelevel_prfile
        datap[6:dtdim-1] = 1
        datap[dtdim-1] = rk
        Internal.createUniqueChild(o, '.Solver#dtloc', 'DataArray_t', datap)

    # Data for each zone
    #==== Check if padding file exists (potential cache associativity issue)
    pad = False
    if Padding is not None:
      try:
          f       = open(Padding,'rb')
          pad     = True
      except IOError:
          print('Warning: Fast: Padding file not found, using default values.')
          pad   = False
      if pad: 
          lines_f = f.readlines()

          sizeIJK=[]
          shiftvars=[]
          for l in lines_f:
             mots = l.split(',')
             sizeIJK.append([ int( mots[1] ),  int( mots[2] ),  int( mots[3]) ] )
             shiftvars.append( int( mots[5] ) )

          f.close()

    bases = Internal.getNodesFromType2(t, 'CGNSBase_t')
    
    i=0
    for b in bases:
        zones = Internal.getNodesFromType2(b, 'Zone_t')
        nzones=len(zones)
        for z in zones:
            shiftvar   = 0
            #=== check for a padding file  (cache associativity issue)
            dims = Internal.getZoneDim(z)
            if dims[3] == 'NGON':
               ngon=True
            else: 
               ngon = False
               iv = dims[1]-5
               jv = dims[2]-5
               kv = dims[3]-5
               if kv == -3 : kv =1
               if pad: 
                  target = [iv,jv,kv]
                  if target in sizeIJK:
                       l        = sizeIJK.index( target)
                       shiftvar =  shiftvars[l]

            # zone ownData (generated)
            o = Internal.createUniqueChild(z, '.Solver#ownData', 
                                           'UserDefinedData_t')

            # - defaults -
            model    = "Euler"
            ransmodel= 'SA'
            sgsmodel = "Miles"
            wallmodel= "power"
            wallmodel_sample= 4  # position du point (en maille )par rapport paroi pour calcul utau
            des      = "none"
            implicit_solver = "lussor"
            nbr_relax = 1
            nbr_krylov    = 20
            nbr_restart   =  1
            lu_match      =  0
            temporal_scheme = "implicit"
            scheme          = "ausmpred"
            slope           = "o3"
            if ngon : slope ='o2'
            motion          = "none"
            filtrage        = "off"
            io_th      = 0
            cacheblckI = 2048
            cacheblckJ = 2         # modifier sir LES car 3 valeur minimal dans ce cas
            cacheblckK = 7           
            dtnature   = "global"
            dtc        = -0.000001
            epsi_newton  = 0.1 
            epsi_linear  = 0.01 
            psiroe       = 0.1
            cfl          = 1.
            rotation     = [ 0.,0.,0., 0.,0.,0.,0.,0.]
            ssdom_IJK    = [240,20,900]
            sfd          = 0
            sfd_chi      = 0.
            sfd_delta    = 1.e15
            sfd_init_iter= 1
            nit_inflow   = 10
            epsi_inflow  = 1.e-5
            DES_debug    = 0
            extract_res  = 0
            ibc          = numpy.zeros(7, dtype=numpy.int32)
            source       = 0
            cups         = [1.,1.]
            ratiom       = 10000.
            meshtype     = 1  #stuctured
            senseurtype  = 1  #version celia laurent du schema senseur
            
            a = Internal.getNodeFromName1(z, 'ZoneType')
            tmp = Internal.getValue(a)
            if tmp == 'Unstructured':
              meshtype =2
              print('Warning: Fast: zone %s is unstructured.'%z[0])

            a = Internal.getNodeFromName2(z, 'GoverningEquations')
            if a is not None: model = Internal.getValue(a)
            else:
                a = Internal.getNodeFromName2(b, 'GoverningEquations')
                if a is not None: model = Internal.getValue(a)
            a = Internal.getNodeFromName2(z, 'DES')
            if a is not None: des = Internal.getValue(a)
            else:
                a = Internal.getNodeFromName2(z, 'des')
                if a is not None: des = Internal.getValue(a)
                else:
                    a = Internal.getNodeFromName2(b, 'DES')
                    if a is not None: des = Internal.getValue(a)
                    else:
                        a = Internal.getNodeFromName2(b, 'des')
                        if a is not None: des = Internal.getValue(a)
            ref = None
            a = Internal.getNodeFromName1(z, 'ReferenceState')
            if a is not None: ref = a
            else:
                a = Internal.getNodeFromName1(b, 'ReferenceState')
                if a is not None: ref = a
            adim = None
            if ref is not None: 
              adim      = C.getState(ref)
              prandtltb = adim[18]        #prandtl laminar 
            else: 
              print('Warning: Fast: can not find a refState.')
              print('Warning: Fast: Prandtl turb by default= 1.')
              prandtltb = 1.

            d = Internal.getNodeFromName1(z, '.Solver#define')
            if d is not None:
                checkKeys(d, keys4Zone)
                a = Internal.getNodeFromName2(b, 'temporal_scheme')
                if a is not None: temporal_scheme = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'implicit_solver')
                if a is not None: implicit_solver = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'nb_relax')
                if a is not None: nbr_relax = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'nb_krylov')
                if a is not None: nbr_krylov = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'nb_restart')
                if a is not None: nbr_restart = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'lu_match')
                if a is not None: lu_match = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'scheme')
                if a is not None: scheme = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'slope')
                if a is not None: slope  = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'motion')
                if a is not None: motion = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'rotation')
                if a is not None: rotation = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'filtrage')
                if a is not None: filtrage = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'io_thread')
                if a is not None: io_th = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'cache_blocking_I')
                if a is not None: cacheblckI = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'cache_blocking_J')
                if a is not None: cacheblckJ = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'cache_blocking_K')
                if a is not None: cacheblckK = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'shiftvar')
                if a is not None: shiftvar = Internal.getValue(a)                
                a = Internal.getNodeFromName1(d, 'time_step_nature')            
                if a is not None: dtnature = Internal.getValue(a)                
                a = Internal.getNodeFromName1(d, 'sgsmodel')
                if a is not None: sgsmodel = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'wallmodel')
                if a is not None: wallmodel = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'wallmodel_sample')
                if a is not None: wallmodel_sample = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'time_step')
                if a is not None: dtc = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'inj1_newton_tol')
                if a is not None: epsi_inflow = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'inj1_newton_nit')
                if a is not None: nit_inflow = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'epsi_newton')
                if a is not None: epsi_newton = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'epsi_linear')
                if a is not None: epsi_linear = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'cfl')
                if a is not None: cfl = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'ssdom_IJK')
                if a is not None: ssdom_IJK = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'prandtltb')
                if a is not None: prandtltb = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'psiroe')
                if a is not None: psiroe = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'sfd')
                if a is not None: sfd = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'sfd_chi')
                if a is not None: sfd_chi = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'sfd_delta')
                if a is not None: sfd_delta = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'sfd_init_iter')
                if a is not None: sfd_init_iter = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'ransmodel')
                if a is not None: ransmodel = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'DES_debug')
                if a is not None: 
                    tmp = Internal.getValue(a)
                    if tmp == 'active': DES_debug =1
                a = Internal.getNodeFromName1(d, 'extract_res')
                if a is not None: extract_res = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'IBC')
                if a is not None: ibc = a[1]
                a = Internal.getNodeFromName1(d, 'source')
                if a is not None: source = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'Cups')
                if a is not None: cups = Internal.getValue(a)  
                a = Internal.getNodeFromName1(d, 'ratiom')
                if a is not None: ratiom = Internal.getValue(a)  
                a = Internal.getNodeFromName1(d, 'senseurType')
                if a is not None: senseurtype = Internal.getValue(a)  
               
            iflow  = 1
            ides   = 0; idist = 1; ispax = 2; izgris = 0; iprod = 0
            azgris = 0.01; addes = 0.2

            if not ngon:
              neq_lbm=19
              if kv == 1: neq_lbm=9  # calcul 2d
            else:
              neq_lbm=0

            if   model == "Euler"     or model == "euler":     iflow = 1
            elif model == "NSLaminar" or model == "nslaminar": iflow = 2
            elif model == "NSLes"     or model == "nsles":     iflow = 2

            elif model == "nsspalart" or model == "NSTurbulent": 
                iflow = 3
                if   ransmodel == 'SA_comp':
                    ides = 1
                elif ransmodel == 'SA_diff':
                    ides = 8

                if des != 'none':
                    if   des == "ZDES1"   or des == "zdes1":   ides = 2
                    elif des == "ZDES1_w" or des == "zdes1_w": ides = 3
                    elif des == "ZDES2"   or des == "zdes2":   ides = 4
                    elif des == "ZDES2_w" or des == "zdes2_w": ides = 5
                    elif des == "ZDES3"   or des == "zdes3":   ides = 6
                    elif des == "ZDES3_w" or des == "zdes3_w": ides = 7

            elif model == "LBMLaminar" or model == "lbmlaminar": 
                iflow   = 4
            else:
                 print('Warning: Fast: model %s is invalid.'%model)

            iles = 0
            if sgsmodel == 'smsm': iles = 1
            if iles ==1:
               cacheblckI = max(cacheblckI,4)
               cacheblckJ = max(cacheblckJ,4)
               cacheblckK = max(cacheblckK,4)

            iwallmodel =1
            if wallmodel == 'musker': iwallmodel = 0

            #par defaut valeur petite, car determine le niveau de parallelisme pour calcul residu en implicit et explict
            #si valeur grande, alors calcul sequentiel de cprdus1
            size_ssdom_I = 100
            size_ssdom_J = 10
            size_ssdom_K = 10


            if temporal_scheme == "explicit" or temporal_scheme == "explicit_local": itypcp = 2
            else: itypcp = 1

            if temporal_scheme == "implicit_local": 
                itypcp = 1
                #ssdom_IJK = [80,15,80]
                a = Internal.getNodeFromName1(d, 'ssdom_IJK')
                if a is not None: 
                    ssdom_IJK = Internal.getValue(a)
                size_ssdom_I = ssdom_IJK[0]
                size_ssdom_J = ssdom_IJK[1]
                size_ssdom_K = ssdom_IJK[2]

            islope=2
            if   slope == "o1"    : islope = 1
            elif slope == "o2"    : islope = 4
            elif slope == "o3"    : islope = 2
            elif slope == "minmod": islope = 3

            kfludom = 1
            if   scheme == "ausmpred" : kfludom = 1
            elif scheme == "senseur"  : kfludom = 2
            elif scheme == "dfcentre6": kfludom = 3
            elif scheme == "roe_kap"  :
              kfludom = 5
              islope  = 2
            elif scheme == "roe_min": 
              kfludom = 5
              islope  = 3
            elif scheme == "roe_nul": 
              kfludom = 5
              islope  = 1
            elif scheme == "roe": 
              kfludom = 5
            elif scheme == "ausmpred_pattern":  
                kfludom = 7

            else: print('Warning: Fast: scheme %s is invalid.'%scheme)

            lale   = 0; size_ale =0
            if    motion == "none"       : lale = 0
            elif  motion == "rigid"      : lale = 1; size_ale = 11
            elif  motion == "deformation": lale = 2
            else: print('Warning: Fast: motion %s is invalid.'%motion)

            iflagflt = 0
            if filtrage == "on": iflagflt = 1

            dtloc = 0
            if dtnature == "global": dtloc = 0
            elif dtnature == "local": dtloc = 1
            else: print('Warning: FastS: time_step_nature %s is invalid.'%dtnature)

            ImplicitSolverNum = 0
            if implicit_solver == "lussor": ImplicitSolverNum = 0
            elif implicit_solver == "gmres": ImplicitSolverNum = 1
            else: print('Warning: FastS: implicit_solver %s is invalid.'%implicit_solver)

            # creation noeud parametre integer
            # Determination de levelg et leveld             
            levelg=0
            leveld=0

            
            datap = numpy.empty(90, numpy.int32)
            datap[0:25]= -1

            #Structure ou polyhedric
            ngon   = Internal.getNodeFromName1(z   , 'NGonElements')
            if ngon is not None:
               nface  = Internal.getNodeFromName1(z   , 'NFaceElements')
               #Nombre de vextex
               grid      = Internal.getNodeFromName1(z    , 'GridCoordinates')
               coordx    = Internal.getNodeFromName1(grid , 'CoordinateX')
               nvtx      = numpy.size( coordx[1])
               datap[0 ] = nvtx

               #Nombre de face
               ng_intext   = Internal.getNodeFromName1(ngon, 'IntExt')[1]
               NF_I0  = ng_intext[0]
               NF_RAC = ng_intext[1]
               NF_I1  = ng_intext[2]
               NF_BC1 = ng_intext[3]
               NF_BC0 = ng_intext[4]
               NF_TOT = ng_intext[6]

               #datap[1 ] = NF_I0+ NF_I1 + NF_I2 + NF_RAC + NF_BC0 + NF_BC1  #nombre de face total
               datap[1 ] = NF_TOT                                           #nombre de face total
               datap[2 ] = NF_I0                                            #Face couche zero (sans BC et rac)
               datap[3 ] = NF_I1                                            #Face couche  un  (sans BC et rac)
               datap[4 ] = NF_RAC                                           #Face couche zero (BC)
               datap[5 ] = NF_BC0                                           #Face couche zero (BC)
               datap[6 ] = NF_BC1                                           #Face couche  un  (BC)

               nf_intext = Internal.getNodeFromName1(nface, 'IntExt')[1]
               NG_0  = nf_intext[0]  #element couche zero
               NG_1  = nf_intext[1]  #element couche 1 de type raccord
               NG_2  = nf_intext[2]  #element couche 2 de type raccord
               NG_3  = nf_intext[3]  #element couche 1 et 2 de type BC
               #NG_4  = nf_intext[4]  #element couche 2 de type raccord

               #datap[7 ] = NG_0 + NG_1 + NG_2 + NF_BC0 + NF_BC1 # elements total (avec les element BC )
               datap[7 ] = NG_0 + NG_1 + NG_2 +NG_3           # elements total (avec les element BC et racc ) NELTS
               print("shiftvar NGON=",shiftvar)
               datap[41] = NG_0 + NG_1 + NG_2 +NG_3 +shiftvar # elements total (avec les element BC et racc ) NDIMDX pour partage transfert Fast
               datap[8 ] = NG_0
               datap[9 ] = NG_1
               datap[10] = NG_3
            
               #size ElementConnectivity (NGON)
               ngon_elconn = Internal.getNodeFromName1(ngon, 'ElementConnectivity')
               ngon_elconn_size = numpy.size( ngon_elconn[1])
               datap[11 ]  = ngon_elconn_size
               #size ElementConnectivity (NFACE)
               nfac_elconn = Internal.getNodeFromName1(nface, 'ElementConnectivity')
               nfac_elconn_size = numpy.size( nfac_elconn[1])
               datap[12]  = nfac_elconn_size
               #Nb cache bloc
               Nb_Cache_Bloc = numpy.shape(Internal.getNodeFromName1(z, 'CellScheduler')[1])[2]
               print("Nb_Cache_Bloc",Nb_Cache_Bloc)
               datap[13]  = Nb_Cache_Bloc

            datap[25]  = 0     # zone 3D curvi par defaut
            datap[26]  = 0     # Vent 3D par defaut
            datap[27]  = iflow
            datap[28]  = iles
            datap[29]  = itypcp
            datap[30]  = size_ssdom_I
            datap[31]  = size_ssdom_J
            datap[32]  = size_ssdom_K
            datap[33]  = kfludom
            datap[34]  = lale
            datap[35]  = iflagflt
            if ngon is not None:
               datap[36:41]= -1
               datap[42:45]= -1
            else:
               datap[36:45]= -1
            datap[45]  = io_th    
            datap[46]  = dtloc
            datap[47]  = ides  
            datap[48]  = idist 
            datap[49]  = ispax
            datap[50]  = izgris
            datap[51]  = iprod
            datap[52]  = rk
            datap[53]  = veclevel[i]
            datap[54]  = exploc
            datap[55]  = exploctype
            datap[56]  = levelg
            datap[57]  = leveld
            datap[58]  = nssiter
            datap[59]  = cacheblckI
            datap[60]  = cacheblckJ
            datap[61]  = cacheblckK
            datap[62]  = sfd
            datap[63]  = sfd_init_iter
            datap[64]  = islope
            datap[65]  = nit_inflow
            datap[66]  = shiftvar
            datap[67]  = extract_res
            datap[68]  = DES_debug     #index 69 et 70 pour PT_BC et PT_OMP
            datap[71]   = nbr_relax
            datap[72]   = nbr_restart
            datap[73]   = nbr_krylov
            datap[74]   = ImplicitSolverNum
            datap[75]   = lu_match          
            datap[76:83]= ibc[0:7]
            datap[83]   = source
            datap[84]   = meshtype
            datap[85]   = senseurtype
            datap[86]   = neq_lbm
            datap[87]   = -1
            datap[88]   = iwallmodel
            datap[89]   = wallmodel_sample

            i += 1
         
            Internal.createUniqueChild(o, 'Parameter_int', 'DataArray_t', datap)
            
            # creation noeud parametre real 
            #size_real = 23 +  size_ale   
            size_real = 41
            datap = numpy.zeros(size_real, numpy.float64)
            if dtc < 0: 
                print('Warning: Fast: time_step set to default value (0.0001).')
                dtc = 0.0001
            datap[0] = dtc
            
            if adim is not None:
                pinf    = adim[5]
                ainf    = math.sqrt(adim[11]*pinf/adim[0])
                datap[ 1]=  adim[11]    # gamma
                datap[ 2]=  adim[ 7]    # cv
                datap[ 3]=  adim[ 0]    # roinf
                datap[ 4]=  adim[ 5]    # pinf
                datap[ 5]=  math.sqrt(adim[1]**2 + adim[2]**2 + adim[3]**2) / adim[0] # uinf
                datap[ 6]=  adim[6]      # Tinf
                datap[ 7]=  adim[ 8]     # Minf
                datap[ 8]=  adim[14]     # roNutildeinf
                datap[ 9]=  adim[ 9]     # Re
                datap[10]=  adim[18]     # Pr
                datap[11]=  adim[15]     # Mus
                datap[12]=  adim[17]     # Ts (sutherland)
                datap[13]=  adim[16]     # Cs (sutherland)
                datap[19]= adim[1]/adim[0] # Vx inf
                datap[20]= adim[2]/adim[0] # Vy inf
                datap[21]= adim[3]/adim[0] # Vz inf
            else: datap[1:14] = 0. ; datap[19:22] = 0.

            datap[14]=  epsi_newton
            datap[15]=  cfl        
            datap[16]=  azgris
            datap[17]=  addes
            datap[18]=  ratiom
            datap[22]=  temps
 
            if lale ==1: 
                datap[23:31]=  rotation[0:8]
                datap[29]   = datap[29]*2.  *math.pi
                datap[30]   = datap[30]/180.*math.pi
                datap[31]   = t_init_ale
                datap[32]   = 0.     #teta
                datap[33]   = 0.     #tetap

            datap[34]=  psiroe
            datap[35]=  sfd_chi
            datap[36]=  sfd_delta
            datap[37]=  epsi_inflow
            datap[38]=  prandtltb
            datap[39]=  epsi_linear
            datap[40]=  cups[1] # on garde la valeur max, pas la moyenne

            Internal.createUniqueChild(o, 'Parameter_real', 'DataArray_t', datap)

            # More
            Internal.createUniqueChild(o, 'CFL_minmaxmoy', 'DataArray_t', [0.,0.,0.])
            Internal.createUniqueChild(o, 'type_zone'    , 'DataArray_t',  0)
            Internal.createUniqueChild(o, 'model', 'DataArray_t', model)
            Internal.createUniqueChild(o, 'temporal_scheme', 'DataArray_t', temporal_scheme)
            Internal.createUniqueChild(o, 'exp_local', 'DataArray_t', exploc)
            Internal.createUniqueChild(o, 'rk', 'DataArray_t', rk)

    return None

#==============================================================================
# create work arrays
#==============================================================================
def createWorkArrays__(zones, dtloc, FIRST_IT):
    ndimt=0; ndimcoe=0; ndimwig=0; ndimgrad=0; ndimplan=0; c = 0
    global  MX_OMP_SIZE_INT
    # on check sur 1ere zone si implicite: mixte impli/expli impossible pour le moment
    scheme = "implicit"
    a = Internal.getNodeFromName2(zones[0], 'temporal_scheme')
    if a is not None: scheme = Internal.getValue(a)
    rk=3
    rk_ = Internal.getNodeFromName2(zones[0],'rk')
    if rk_ is not None: rk = Internal.getValue(rk_)
   
    exploc=0
    exploc_ = Internal.getNodeFromName2(zones[0], 'exp_local')
    if exploc_ is not None: exploc = Internal.getValue(exploc_)

    if   scheme == 'implicit': lssiter_loc = 0
    elif scheme == 'implicit_local':  lssiter_loc = 1  # sous-iteration local
    else: lssiter_loc = 0

    neq_max =0
    for z in zones:
        model = "Euler"
        a = Internal.getNodeFromName2(z, 'model')
        model  = Internal.getValue(a)
        neq = 5
        if model == 'nsspalart' or model =='NSTurbulent': neq = 6
        param_int = Internal.getNodeFromName2(z, 'Parameter_int')[1]
        if model == 'LBMLaminar' : neq = param_int[86]
        
        if scheme == 'implicit' or scheme == 'implicit_local':   neq_coe = neq     
        else:                                                    neq_coe = 1    # explicit

        neq_max = max(neq, neq_max)
        shiftvar = param_int[66]
        if param_int is not None:
           shiftvar  = param_int[66]
        else:
           shiftvar = 100
           print('shiftVar=%d'%shiftvar)
           print('create workarray')
           print('Danger sur optimisation cache shiftVar')      
           
        # dim from getZoneDim
        dims = Internal.getZoneDim(z)
        if dims[0]=='Structured' : 
            nijk = (dims[1]-1)*(dims[2]-1)*(dims[3]-1)+shiftvar
            dimJK    = (dims[2]-1)*(dims[3]-1); dimIK = (dims[1]-1)*(dims[3]-1); dimIJ=(dims[1]-1)*(dims[2]-1)
            ndimplan = max( ndimplan, dimIJ)
            ndimplan = max( ndimplan, dimIK)
            ndimplan = max( ndimplan, dimJK)

            ndimFlu  =0
        else: 
            nijk = dims[2]+shiftvar
            print('dtloc fastP: revoir dim tab jeanmasson')
            ndimplan = 10     ### a revoir
            #ndimFlu  =2500000*neq*3
            ndimFlu  =0
            print('FastC: ndimFlu=0, a revoir pour vecto')
            ndimgrad+= neq_max*3*nijk

        #print 'zone,taille =', z[0], nijk

        ndimt   +=     neq*nijk       # surdimensionne ici
        ndimcoe += neq_coe*nijk       # surdimensionne ici
        ndimwig +=        3*nijk      # surdimensionne ici

        c       += 1
        
    #           3 depth     6 faces  5 tableau
    ndimface= neq*3*ndimplan*6*5*len(zones)
    ndimface = min(ndimface, 2000000000)
    #si pas de temps local inactif (rk3)
    #if rk!=3 or exploc !=2: ndimface=1
    if exploc !=1: ndimface=1
    else: print('taille tab dtloc=%d'%ndimface)

    mx_thread   = OMP_NUM_THREADS       # surdimensionne : doit etre = a OMP_NUM_THREADS
    verrou      = MX_SSZONE*c*MX_SYNCHRO*mx_thread

    timerOmp = numpy.zeros(  mx_thread*2*dtloc[0] + (mx_thread+1)*len(zones), dtype=numpy.float64)

#    wig   = KCore.empty(ndimwig, CACHELINE)
#    coe   = KCore.empty(ndimcoe, CACHELINE)
#    drodm = KCore.empty(ndimt  , CACHELINE)
    wig   = numpy.empty(ndimwig , dtype=numpy.float64)
    coe   = numpy.empty(ndimcoe , dtype=numpy.float64)
    flu   = numpy.empty(ndimFlu , dtype=numpy.float64)
    grad  = numpy.empty(ndimgrad, dtype=numpy.float64)

    if   (rk==4 and exploc==0): size_drodm = 4*ndimt
    elif (rk==5 and exploc==0): size_drodm = 5*ndimt
    else                      : size_drodm =   ndimt
    drodm     = numpy.empty(ndimt   , dtype=numpy.float64)        
    tab_dtloc = numpy.empty(ndimface, dtype=numpy.float64)

#    for z in zones:
#         param_int = Internal.getNodeFromName2(z, 'Parameter_int')  # noeud
#    
#    ptr = drodm1.reshape((nijk*5), order='Fortran') # no copy I hope
#    for c in range(5):
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

    hook = {}
    hook['wiggle']         = wig
    hook['coe']            = coe
    hook['rhs']            = drodm
    hook['grad']           = grad
    hook['tab_dtloc']      = tab_dtloc
    hook['flu_vecto']      = flu
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
    hook['param_int_tc']   = None
    hook['param_real_tc']  = None
    hook['mpi']            = 0     
    hook['linelets_int']   = None
    hook['linelets_real']  = None
    return hook

#==============================================================================
# loi horaire (ALE) 
#==============================================================================
def _motionlaw(t, teta, tetap):
    zones = Internal.getZones(t)
    fastc._motionlaw(zones, teta, tetap)
    return None

#==============================================================================
# Retourne un tag suivant bcname
#==============================================================================
def tagBC(bcname):
  if   bcname == "BCExtrapolate"    :       tag = 0
  elif bcname == "BCDegenerateLine":        tag = 0
  elif bcname == "BCFarfield":              tag = 1
  elif bcname == "BCInflowSupersonic":      tag = 2
  elif bcname == "BCWallInviscid":          tag = 3
  elif bcname == "BCSymmetryPlane":         tag = 3
  elif bcname == "BCWall":                  tag = 4
  elif bcname == "FamilySpecified":         tag = 5
  elif bcname == "BCWallViscous":           tag = 6
  elif bcname == "BCWallViscous_isot_fich": tag = 7
  elif bcname == "BC_fich":                 tag = 8
  elif bcname == "Nearmatch":               tag = 9
  elif bcname == "BCOutflow":               tag =10
  elif bcname == "BCautoperiod":            tag =11
  elif bcname == "BCWallViscous_transition":tag =12
  elif bcname == "BCInflow":                tag =13
  elif bcname == "BCExtrapolateRANS":       tag =14
  elif bcname == "BCPeriodic":              tag =15
  elif bcname == "BCOutpres":               tag =16
  elif bcname == "BCInj1":                  tag =17
  elif bcname == "BCRacinf":                tag =18
  elif bcname == "BCInflowLund":            tag =19
  elif bcname == "BCInjMFR":                tag =20
  elif bcname == "BCOutMFR":                tag =21
  elif bcname == "BCOverlap":               tag =22
  elif bcname == "BCWallModel":             tag =30
  elif bcname == "BCWallExchange":          tag =31
  else:
    tag = -1
    print("Warning: Fast: unknown BC type %s."%bcname)
  return tag

#==============================================================================
# echange M1 <- current, current <- P1, P1 <- M1
#==============================================================================
def switchPointers__(zones, case, order=3):
    
    if order == 2: return None
    elif order == 1 or order == 3:
      for z in zones:
          sol =  Internal.getNodeFromName1(z, 'FlowSolution#Centers')
          own =  Internal.getNodeFromName1(z, '.Solver#ownData')
          vars=['Density','VelocityX','VelocityY','VelocityZ','Temperature']
          model  = Internal.getNodeFromName1(own, 'model')
          model  = Internal.getValue(model)
          if model == 'nsspalart' or model=='NSTurbulent':
            vars.append('TurbulentSANuTilde')

          for v in vars:
            ca   = Internal.getNodeFromName1(sol, v)
            caM1 = Internal.getNodeFromName1(sol, v+'_M1')
            caP1 = Internal.getNodeFromName1(sol, v+'_P1')

            if case ==1:
               #sauvegarde M1   #M1 <- current    # current <- P1   # P1 <- temp 
               ta = caM1[1];    caM1[1] = ca[1];  ca[1] = caP1[1];  caP1[1] = ta  
            elif case ==2:
               #sauvegarde P1   #P1 <- current    # current <- M1   # M1 <- temp 
               ta = caP1[1];    caP1[1] = ca[1];  ca[1] = caM1[1];  caM1[1] = ta
#==============================================================================
# Compact vars in Container or defined in fields
#==============================================================================
def getField2Compact__(z, field):
    field = field.split(':')
    if len(field) == 1: 
        container = Internal.__FlowSolutionNodes__ ; name = field[0]
    elif field[0] == 'nodes': 
        container = Internal.__FlowSolutionNodes__ ; name = field[1]
    elif field[0] == 'Motion': 
        container = 'Motion' ; name = field[1]
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
    wrange = Internal.getNodesFromName(node, 'PointRange')
    win = Internal.range2Window(wrange[0][1])
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
      if i[0] in keys:
        k = keys[i[0]]
        val = Internal.getValue(i)
        if k == 0: # must be int
          if not isinstance(val, (int)):
            fail = True; print('Error: Fast: keyword %s requires an int.'%i[0])
        elif k == 1: # must be float but we accept int
          if not isinstance(val, (float, int)):
            fail = True; print('Error: Fast: keyword %s requires a float.'%i[0])
        elif k == 2: # must be any string
          if not isinstance(val, str):
            fail = True; print('Error: Fast: keyword %s requires a string.'%i[0])
        elif k == 3: # Array/list of ints
          if not isinstance(val, numpy.ndarray):
            fail = True; print('Error: Fast: keyword %s requires an array of ints.'%i[0])
        elif k == 4: # Array/list of floats
          if not isinstance(val, numpy.ndarray):
            fail = True; print('Error: Fast: keyword %s requires an array of floats.'%i[0])
        elif isinstance(k, list): # list of given strings
          if not val in k:
            fail = True; print('Error: Fast: keyword %s requires a string in %s.'%(i[0], str(k)))
      else:
        print('Error: Fast: keyword %s is INVALID.'%i[0])
        import difflib
        words = keys.keys()
        m = difflib.get_close_matches(i[0], words)
        print('Error: Fast: do you mean %s?'%str(m))
        fail = True
    if fail:
      import sys; sys.exit(0) # much stronger os._exit(0)?



#==============================================================================
# echange M1 <- current, current <- P1, P1 <- M1
#==============================================================================
def switchPointers2__(zones,nitmax,nstep):

    for z in zones:
        dim = Internal.getZoneDim(z)
        i = dim[1]//2
        j = dim[2]//2
        k = dim[3]//2
        niveau = C.getValue( z, 'centers:niveaux_temps', (i,j,k))
    

        cycle = nitmax/niveau
        if nstep%cycle==0 and nstep != nitmax:

            #print('switch zone= ', z[0])

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

            #print 'ca[1]= ',ca[1]

            # sauvegarde current
            ta = ca[1]; tb = cb[1]; tc = cc[1]; td = cd[1]; te = ce[1]

            # M1 <- current
            #caM1[1] = ca[1]; cbM1[1] = cb[1]; ccM1[1] = cc[1]; cdM1[1] = cd[1]; ceM1[1] = ce[1]

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


#==============================================================================
# echange M1 <- current, current <- P1, P1 <- M1
#==============================================================================
def switchPointers3__(zones,nitmax):

    for z in zones:
        dim = Internal.getZoneDim(z)
        i = dim[1]//2
        j = dim[2]//2
        k = dim[3]//2
        niveau = C.getValue( z, 'centers:niveaux_temps', (i,j,k))
    
        cycle = nitmax/niveau
        if cycle==nitmax :

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

        else:


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
            caM1[1] = caP1[1]; cbM1[1] = cbP1[1]; ccM1[1] = ccP1[1]; cdM1[1] = cdP1[1]; ceM1[1] = ceP1[1]

            # P1 <- temp
            caP1[1] = ta; cbP1[1] = tb; ccP1[1] = tc; cdP1[1] = td; ceP1[1] = te


#==============================================================================
# echange M1 <- current, current <- P1, P1 <- M1
#==============================================================================
def switchPointersLBM__(zones, neq_lbm):
    
     for z in zones:
        sol =  Internal.getNodeFromName1(z, 'FlowSolution#Centers')
        for i in range(1,neq_lbm+1):
            caM1 = Internal.getNodeFromName1(sol, 'Q'+str(i)+'_M1')
            ca   = Internal.getNodeFromName1(sol, 'Q'+str(i))

            # sauvegarde M1
            ta = caM1[1]

            # M1 <- current
            caM1[1] = ca[1]

            # current <- P1
            ca[1] = ta

#==============================================================================
# Construit les donnees compactees pour traiter les BC
#==============================================================================
def _BCcompactNG(t):
    #dimensionnememnt tableau
    size_param_int =[]
    size_param_real=[]
    zones = Internal.getZones(t)
    for z in zones:

        bcs           = Internal.getNodesFromType2(z, 'BC_t')
        FaceScheduler = Internal.getNodeFromName1(z, 'FaceScheduler')[1]
        Nb_bc = len(bcs)
        size_int  =  numpy.size(FaceScheduler)
        size_real = 0
        for bc in bcs:
            bcdata  = Internal.getNodesFromType3(bc, 'DataArray_t')
            Nb_data = len(bcdata)
            #size_int= size_int  + Nb_data 
            for data in bcdata:
               print("Data in BC not implemented in FastP")
               import sys; sys.exit()
               size = numpy.shape(data[1])
               c = 1
               for s in size: c=c*s
               size_real = size_real + c

        # zone ownData (generated)
        o = Internal.getNodeFromName1( z, '.Solver#ownData')

        #on concatene les donnes BC dans param_int et param_real
        param_int = Internal.getNodeFromName1( o, 'Parameter_int')
        size = numpy.shape(param_int[1])
        c = 1
        for s in size: c=c*s
        size_param_int.append(c)

        #print("BC size facesched", size_int, "old param int=",c)
        datap = numpy.zeros(size_int+c, numpy.int32)
        datap[0:c]   = param_int[1][0:c]
        param_int[1] = datap

        param_real = Internal.getNodeFromName1( o, 'Parameter_real')
        size = numpy.shape(param_real[1])
        c = 1
        for s in size: c=c*s
        size_param_real.append(c)

        if size_real != 0:
            datap = numpy.zeros(size_real+c, numpy.float64)
            datap[0:c]    = param_real[1][0:c]
            param_real[1] = datap

    #Init param BC
    no = 0
    zones = Internal.getZones(t)
    for z in zones:

        o             = Internal.getNodeFromName1( z, '.Solver#ownData')
        param_int     = Internal.getNodeFromName1( o, 'Parameter_int')[1]
        FaceScheduler = Internal.getNodeFromName1(z, 'FaceScheduler')[1]

        pt_bcs_int   =  size_param_int[no]
        param_int[70]=  pt_bcs_int
        no           += 1

        size = numpy.shape(FaceScheduler)
        size_scheduler = 1
        for s in size: size_scheduler=size_scheduler*s

        deb = pt_bcs_int
        c   = 0
        for l in range(size[3]):
         for k in range(size[2]):
          for j in range(size[1]):
            for i in range(size[0]):
              param_int[ deb +c ]= FaceScheduler[ i, j, k, l]
              #print("verif fsched=",FaceScheduler[ i, j, k, l],l,k,j,i, "deb=",deb,"c=",c)
              c +=1
        FaceScheduler = param_int[ deb:deb+size_scheduler ] # wrong shape?

#==============================================================================
# Construit les donnees compactees pour traiter les BC
#==============================================================================
def _BCcompact(t):
    
    #dimensionnememnt tableau
    size_param_int =[]
    size_param_real=[]
    zones = Internal.getZones(t)
    nzones=len(zones)
    familyBCTypes = C.getFamilyBCNamesDict(t)
    for z in zones:

        bcs   = Internal.getNodesFromType2(z, 'BC_t')
        Nb_bc = len(bcs)
        size_int  = 1 +  Nb_bc*2 + 9*Nb_bc
        size_real = 0
        for bc in bcs:
            btype = Internal.getValue(bc)
            ## Si la bc est une wall_viscous_transition on ajoute un data qui est un vecteur contenant des nombres random
            if btype == 'BCWallViscous_transition':
                ptrange = Internal.getNodesFromType1(bc, 'IndexRange_t')
                wrange  = ptrange[0][1]
                nb_cell = max([wrange[0][1] - wrange[0][0],1])* max([wrange[1][1] - wrange[1][0],1])*max([wrange[2][1] - wrange[2][0],1])
                rand = numpy.zeros((1,nb_cell), dtype=numpy.float64)
                Internal.createUniqueChild(bc, 'random_vec', 'DataArray_t', rand)

            if btype == 'BCWallExchange':
                sol   = Internal.getNodeFromName1(z, 'FlowSolution#Centers')
                cellN = Internal.getNodeFromName1(sol, 'cellN')
                if cellN == None:
                  C._initVars(z, 'centers:cellN', 1.)
                  cellN = Internal.getNodeFromName1(sol, 'cellN')[1]
                  ptrange = Internal.getNodesFromType1(bc, 'IndexRange_t')
                  wrange  = ptrange[0][1]
                  shift =1
                  if '2couche' in bc[0]: shift=2
                  if '3couche' in bc[0]: shift=3
                  print ('nb couche loi de paroi',shift)
                  

                  #CL I
                  if wrange[0][1] - wrange[0][0]==0:
                     if wrange[0][1]==1: ijk_tg =2;              sens = 1
                     else:               ijk_tg =wrange[0][1]-4; sens =-1

                     for k in range(wrange[2][0]-1,wrange[2][1]):
                        for j in range(wrange[1][0]-1,wrange[1][1]):
                          for s in range(shift):
                            cellN[ijk_tg+ s*sens, j , k]=2.

                  #CL J
                  if wrange[1][1] - wrange[1][0]==0:
                     if wrange[1][1]==1: ijk_tg =2;              sens = 1
                     else:               ijk_tg = wrange[1][1]-4;sens =-1

                     for k in range(wrange[2][0]-1,wrange[2][1]):
                        for i in range(wrange[0][0]-1,wrange[0][1]):
                          for s in range(shift):
                             cellN[i, ijk_tg+ s*sens  , k]=2.

                  #CL K
                  if numpy.shape(cellN)[2]!=1:
                    if wrange[2][1] - wrange[2][0]==0:
                       if wrange[2][1]==1: ijk_tg =2;               sens = 1
                       else:               ijk_tg = wrange[2][1]-4; sens =-1
                       for j in range(wrange[1][0]-1,wrange[1][1]):
                         for i in range(wrange[0][0]-1,wrange[0][1]):
                           for s in range(shift):
                              cellN[i, j, ijk_tg+ s*sens ]=2.
              
            bcdata  = Internal.getNodesFromType3(bc, 'DataArray_t')
            Nb_data = len(bcdata)

            for data in bcdata:
               size = numpy.shape(data[1])
               c = 1
               for s in size: c=c*s
               # a sortie de la boucle: voir fastP
               size_int  = size_int  + Nb_data

               size_real = size_real + c

        # zone ownData (generated)
        o = Internal.getNodeFromName1(z, '.Solver#ownData')

        #on concatene les donnes BC dans param_int et param_real
        param_int = Internal.getNodeFromName1(o, 'Parameter_int')
        size = numpy.shape(param_int[1])

        c = 1
        for s in size: c=c*s
        size_param_int.append(c)

        datap = numpy.zeros(size_int+c, numpy.int32)
        datap[0:c]   = param_int[1][0:c]
        param_int[1] = datap

        param_real = Internal.getNodeFromName1(o, 'Parameter_real')
        size = numpy.shape(param_real[1])
        c = 1
        for s in size: c=c*s
        size_param_real.append(c)

        if size_real != 0:
            datap = numpy.zeros(size_real+c, numpy.float64)
            datap[0:c]    = param_real[1][0:c]
            param_real[1] = datap


    #Init param BC
    no = 0
    zones = Internal.getZones(t)
    nzones=len(zones)
    for z in zones:

        o         = Internal.getNodeFromName1( z, '.Solver#ownData')
        bcs       = Internal.getNodesFromType2(z, 'BC_t')
        param_int = Internal.getNodeFromName1( o, 'Parameter_int')[1]
        param_real= Internal.getNodeFromName1( o, 'Parameter_real')[1]

        Nb_bc = len(bcs)

        #pt_bcs_int = 64
        #pt_bcs_real= 37
        pt_bcs_int   =  size_param_int[no]
        param_int[70]=  pt_bcs_int

        pt_bcs_real=  size_param_real[no] 
        param_int[ pt_bcs_int ] = Nb_bc
        i         = 1
        size_int  = 1 + 2*Nb_bc  # shift pour nb_bc et pointeur BC_int et BC_real
        size_real = 0
        no       += 1
        for bc in bcs:
            param_int[ pt_bcs_int +i       ] = size_int  + pt_bcs_int
            param_int[ pt_bcs_int +i +Nb_bc] = size_real + pt_bcs_real

            pt_bc                            =  param_int[ pt_bcs_int +i ] 

            btype = Internal.getValue(bc)
            if btype == 'FamilySpecified':
              n = Internal.getNodeFromType1(bc, 'FamilyName_t')
              if n is None: raise ValueError('compact: familyName missing in %s.'%bc[0])
              name = Internal.getValue(n)
              if name in familyBCTypes: btype = familyBCTypes[name]
              else: raise ValueError('compact: %s is not a defined familyName.'%name)
              #print('getting family btype=', btype)
            
            param_int[pt_bc] = tagBC( btype )

            Ptrange = Internal.getNodeFromType1(bc, 'IndexRange_t')
            indrange= Internal.getValue(Ptrange)
            ind_bc  = numpy.zeros(6, numpy.int32)
            ind_bc[0] = indrange[0][0]
            ind_bc[1] = indrange[0][1]
            ind_bc[2] = indrange[1][0]
            ind_bc[3] = indrange[1][1]
            ind_bc[4] = indrange[2][0]
            ind_bc[5] = indrange[2][1]

            fastc.PygetRange(ind_bc,  param_int, pt_bc+ 1)

            bcdata  = Internal.getNodesFromType3(bc, 'DataArray_t')
            Nb_data = len(bcdata)
            param_int[pt_bc + 8] = Nb_data
            j = 1
            ctot = 0
            for data in bcdata:
               size = numpy.shape(data[1])
               c = 1
               for s in size: c=c*s
               ctot = ctot +c
               param_int[pt_bc + 8 + j ] = ctot
               deb = pt_bcs_real+size_real
               param_real[ deb:deb+c ]= data[1].flat[0:c]
               data[1]                = param_real[ deb:deb+c ] # wrong shape?
               size_real = size_real + c
               j = j+1

            size_int  = size_int  + 9 + Nb_data
            i         = i + 1

    return None



#==============================================================================
# Load t.cgns (data) tc.cgns (connectivity file) ts.cgns (stat)
# si split='single': load t.cgns
# si split='multiple': read t.cgns ou t/t_1.cgns
# autorestart possible
# return a partial tree t, a donor partial tree tc, 
# a stat partial tree,
# the communication graph for Chimera and abutting transfers
# the communication graph for IBM transfers
#==============================================================================
def load(fileName='t', fileNameC='tc', fileNameS='tstat', split='single',
         restart=False, cartesian=False):
    """Load tree and connectivity tree."""
    import os.path
    baseName = os.path.basename(fileName)
    spl = os.path.splitext(baseName)
    baseName = spl[0] # name without extension
    ext      = spl[1] # extension
    fileName = os.path.splitext(fileName)[0] # full path without extension
    
    baseNameC = os.path.basename(fileNameC)
    spl = os.path.splitext(baseNameC)
    baseNameC = spl[0] # name without extension
    extC      = spl[1] # extension
    fileNameC = os.path.splitext(fileNameC)[0] # full path without extension

    baseNameS = os.path.basename(fileNameS)
    spl = os.path.splitext(baseNameS)
    baseNameS = spl[0] # name without extension
    extS      = spl[1] # extension
    fileNameS = os.path.splitext(fileNameS)[0] # full path without extension

    if ext == '': ext = '.cgns'
    if extC == '': extC = '.cgns'
    if extS == '': extS = '.cgns'

    graph = {'graphID':None, 'graphIBCD':None, 'procDict':None, 'procList':None}
    if Cmpi.size > 1: # mpi run
        import Distributor2.PyTree as D2
        rank = Cmpi.rank; size = Cmpi.size
        
        if split == 'single':
            # Load connect (tc)
            FILE = fileNameC+extC
            tc = Cmpi.convertFile2SkeletonTree(FILE)
            tc = Cmpi.readZones(tc, FILE, rank=rank)
            graphID = Cmpi.computeGraph(tc, type='ID')
            graphIBCD = Cmpi.computeGraph(tc, type='IBCD')
            procDict = D2.getProcDict(tc)
            procList = D2.getProcList(tc, sort=True)
            graph = {'graphID':graphID, 'graphIBCD':graphIBCD, 'procDict':procDict, 'procList':procList }
            tc = Cmpi.convert2PartialTree(tc, rank=rank)
            # Load data (t)
            FILE = fileName+ext
            if restart and os.access('restart.cgns', os.F_OK): FILE = 'restart.cgns'
            if os.access(FILE, os.F_OK):
                t = Cmpi.convertFile2SkeletonTree(FILE)
                t = Cmpi.readZones(t, FILE, rank=rank)
                t = Cmpi.convert2PartialTree(t)
            else: t = None
            # Load stat (ts)
            FILE = fileNameS+extS
            if os.access(FILE, os.F_OK):
                ts = Cmpi.convertFile2SkeletonTree(FILE)
                ts = Cmpi.readZones(ts, FILE, rank=rank)
                ts = Cmpi.convert2PartialTree(ts)
            else: ts = None

        else: # un fichier loade par proc
            # Try to load graph
            if os.access('%s/graph.pk'%fileName, os.R_OK):
                try: import cPickle as pickle
                except: import pickle
                file = open('%s/graph.pk'%fileName, 'rb')
                graph = pickle.load(file)
                file.close()
            FILE = '%s/%s_%d.cgns'%(fileNameC, baseNameC, rank)
            if os.access(FILE, os.F_OK): tc = C.convertFile2PyTree(FILE)
            else:
                FILE = fileNameC+extC
                tc = Cmpi.convertFile2SkeletonTree(FILE)
                tc = Cmpi.readZones(tc, FILE, rank=rank)        
                graphID = Cmpi.computeGraph(tc, type='ID',reduction=False)
                graphIBCD = Cmpi.computeGraph(tc, type='IBCD',reduction=False)
                procDict = D2.getProcDict(tc)
                procList = D2.getProcList(tc, sort=True)
                graph = {'graphID':graphID, 'graphIBCD':graphIBCD, 'procDict':procDict, 'procList':procList }
                tc = Cmpi.convert2PartialTree(tc)
            # Load data (t)
            FILE = '%s/%s_%d%s'%(fileName, baseName, rank, ext)
            if restart and os.access('restart/restart_%d%s'%(rank,ext), os.F_OK):
                FILE = 'restart/restart_%d%s'%(rank,ext)
            if os.access(FILE, os.F_OK): t = C.convertFile2PyTree(FILE)
            else: t = None
            # Load stat (ts)
            FILE = '%s/%s_%d%s'%(fileNameS, baseNameS, rank, ext)
            if os.access(FILE, os.F_OK): ts = C.convertFile2PyTree(FILE)
            else: ts = None
    else: # sequential run
        if split == 'single':
            # Load Connectivity (tc)
            FILE = fileNameC+extC
            if os.access(FILE, os.F_OK): tc = C.convertFile2PyTree(FILE)
            else: tc = None
            # Load Data (t)
            FILE = fileName+ext
            if restart and os.access('restart.cgns', os.F_OK): FILE = 'restart.cgns'
            if os.access(FILE, os.F_OK): t = C.convertFile2PyTree(FILE)
            else: t = None
            # Load Stat (ts)
            FILE = fileNameS+extS
            if os.access(FILE, os.F_OK): ts = C.convertFile2PyTree(FILE)
            else: ts = None
        else: # multiple
            # Load connectivity (tc)
            ret = 1; no = 0; tc = []
            while ret == 1:
                FILE = '%s/%s_%d%s'%(fileNameC, baseNameC, no, extC)
                if not os.access(FILE, os.F_OK): ret = 0
                if ret == 1: tc.append(C.convertFile2PyTree(FILE))
                no += 1
            if tc != []: tc = Internal.merge(tc)
            else: 
                FILE = fileNameC+extC
                if os.access(FILE, os.F_OK): tc = C.convertFile2PyTree(FILE)
                else: tc = None
            # Load data (t)
            ret = 1; no = 0; t = []
            while ret == 1:
                FILEIN = '%s/%s_%d%s'%(fileName, baseName, no, ext)
                if os.access(FILEIN, os.F_OK): FILE = FILEIN
                elif restart and os.access('restart/restart_%d%s'%(no, ext), os.F_OK): 
                    FILE = 'restart/restart_%d%s'%(no,ext)
                else: ret = 0
                if ret == 1: t.append(C.convertFile2PyTree(FILE))
                no += 1
            if t != []: t = Internal.merge(t)
            else: t = None
            # Load stat (ts)
            ret = 1; no = 0; ts = []
            while ret == 1:
                FILE = '%s/%s_%d%s'%(fileNameS, baseNameS, no, extS)
                if not os.access(FILEIN, os.F_OK): ret = 0 
                if ret == 1: ts.append(C.convertFile2PyTree(FILE))
                no += 1
            if ts != []: ts = Internal.merge(ts)
            else: ts = None
    if cartesian: # peut etre inutile (fait dans convert2File?)
        import Compressor.PyTree as Compressor
        Compressor._uncompressCartesian(t)
        Compressor._uncompressCartesian(tc)
    
    return t, tc, ts, graph

#==============================================================================
# save t
# IN: NP: 0 (seq run), >0 (mpi run, distributed), <0 (seq, save multiple)
# IN: split: 'single', 'multiple'
# IN: fileName: name of output file or dir
# single: write restart.cgns (full)
# multiple: write restart/restart_1.cgns, ... (partial trees)
#==============================================================================
def save(t, fileName='restart', split='single', NP=0, cartesian=False):
    """Save tree and connectivity tree."""
    # Rip file ext if any
    import os.path
    baseName = os.path.basename(fileName)
    baseName = os.path.splitext(baseName)[0] # name without extension
    fileName = os.path.splitext(fileName)[0] # full path without extension

    if split != 'single':
        if NP > 0: # mpirun
            if Cmpi.rank == 0:
                if not os.path.exists(fileName): os.makedirs(fileName)
            Cmpi.barrier()
        else:
            if not os.path.exists(fileName): os.makedirs(fileName)

    # Rip some useless data (FastS)
    t2 = C.rmVars(t, 'centers:Density_P1')
    C._rmVars(t2, 'centers:VelocityX_P1')
    C._rmVars(t2, 'centers:VelocityY_P1')
    C._rmVars(t2, 'centers:VelocityZ_P1')
    C._rmVars(t2, 'centers:Temperature_P1')
    C._rmVars(t2, 'centers:TurbulentSANuTilde_P1')
    zones = Internal.getZones(t2)
    for z in zones: Internal._rmNodeByPath(z, '.Solver#ownData')
    if cartesian: 
        import Compressor.PyTree as Compressor
        Compressor._compressCartesian(t2)
        #Compressor._compressCellN(t2)

    flowsol = Internal.getNodeFromName1(zones[0], 'FlowSolution#Centers')
    if flowsol is not None:
        vars = Internal.getNodesFromType1(flowsol, 'DataArray_t')
        for var in vars:
            if ('Kry' in var[0]) and not('Kry_0' in var[0]):
                C._rmVars(t2, 'centers:'+var[0])

    for z in zones:
        node = Internal.getNodeFromName1(z, '.Solver#define')
        if node is not None:
            integ = 'None'; nature = 'None'
            n = Internal.getNodeFromName1(node, 'temporal_scheme')
            if n is not None: integ = Internal.getValue(n)
            n = Internal.getNodeFromName1(node, 'time_step_nature')
            if n is not None: nature = Internal.getValue(n)
            if integ == 'explicit' or nature == 'local':
                C._rmVars(z, 'centers:Density_M1')
                C._rmVars(z, 'centers:VelocityX_M1')
                C._rmVars(z, 'centers:VelocityY_M1')
                C._rmVars(z, 'centers:VelocityZ_M1')
                C._rmVars(z, 'centers:Temperature_M1')
                C._rmVars(z, 'centers:TurbulentSANuTilde_M1')

    # save
    if Cmpi.size > 1: # mpi run
        if split == 'single':  # output in a single file
            Cmpi.convertPyTree2File(t2, fileName+'.cgns')
        else:
            rank = Cmpi.rank
            C.convertPyTree2File(t2, '%s/%s_%d.cgns'%(fileName,baseName,rank))
            # Rebuild graph
            # skeleton -> gather -> merge -> graph
    else: # sequential run
        if split == 'single':
            C.convertPyTree2File(t2, fileName+'.cgns')
        else:
            # Get and save graph
            import Distributor2.PyTree as D2
            graphID = Cmpi.computeGraph(t2, type='ID')
            graphIBCD = Cmpi.computeGraph(t2, type='IBCD')
            
            procDict = D2.getProcDict(t2)
            procList = D2.getProcList(t2, sort=True)
            objet = {'graphID':graphID, 'graphIBCD':graphIBCD, 'procDict':procDict, 'procList':procList}

            # Rebuild local trees
            for i in range(max(1,-NP)):
                tl = Internal.copyRef(t2)
                bases = Internal.getNodesFromType1(tl, 'CGNSBase_t')
                for b in bases:
                    zones = Internal.getNodesFromType1(b, 'Zone_t')
                    for z in zones:
                        proc = Internal.getNodeFromName2(z, 'proc')
                        if proc is not None:
                            proc = Internal.getValue(proc)
                            if proc != i: Internal._rmNode(b, z)
                C.convertPyTree2File(tl, '%s/%s_%d.cgns'%(fileName,baseName,i))

            # Write graph
            try: import cPickle as pickle
            except: import pickle
            file = open('%s/graph.pk'%fileName, 'wb')
            pickle.dump(objet, file, protocol=pickle.HIGHEST_PROTOCOL)
            file.close()

#============================================================================
# Retourne le max proc pour les zones
def getMaxProc(t):
    maxProc = 0
    zones = Internal.getZones(t)
    for z in zones:
        proc = Internal.getNodeFromName2(z, 'proc')
        if proc is not None:
            proc = Internal.getValue(proc)
            maxProc = max(maxProc, proc)
    return maxProc

#==============================================================================
# Load one file
# si split='single': load t.cgns
# si split='multiple': read t.cgns ou t/t_1.cgns
# if graph=True,
# return the communication graph for Chimera and abutting transfers
# and the communication graph for IBM transfers
#==============================================================================
def loadFile(fileName='t.cgns', split='single', graph=False, 
             mpirun=False, cartesian=False, exploc=False):
    """Load tree and connectivity tree."""
    import os.path
    baseName = os.path.basename(fileName)
    baseName = os.path.splitext(baseName)[0] # name without extension
    fileName = os.path.splitext(fileName)[0] # full path without extension

    graphN = {'graphID':None, 'graphIBCD':None, 'procDict':None, 'procList':None}
    if mpirun: # mpi run
        import Distributor2.PyTree as D2
        rank = Cmpi.rank; size = Cmpi.size
        
        if split == 'single':
            FILE = fileName+'.cgns'
            t = Cmpi.convertFile2SkeletonTree(FILE)
            mp = getMaxProc(t)
            #if mp+1 != size: #### COMMENTE PAR GUILLAUME POUR PASSER LE CREATE
            #    raise ValueError('The number of mpi proc (%d) doesn t match the tree distribution (%d).'%(size,mp+1))

            if graph and exploc == False:

                graphID   = Cmpi.computeGraph(t, type='ID'  , reduction=False)
                graphIBCD = Cmpi.computeGraph(t, type='IBCD', reduction=False)
                procDict  = D2.getProcDict(t)
                procList  = D2.getProcList(t, sort=True)
                graphN = {'graphID':graphID, 'graphIBCD':graphIBCD, 'procDict':procDict, 'procList':procList }



            if graph and exploc == True : ### dtloc instationnaire : on sort une liste de graphes (un graph par ssiter)

                list_graph=[]
                graphID   = Cmpi.computeGraph(t, type='ID'  , reduction=False, exploc=True)
                graphIBCD  = Cmpi.computeGraph(t, type='IBCD', reduction=False, exploc=True)
                procDict  = D2.getProcDict(t)
                procList  = D2.getProcList(t,  sort=True)
              
  
                i=0
                graphN={}
                if graphID is not None:
                    for g in graphID:
                        graphN={'graphID':g, 'graphIBCD':{}, 'procDict':procDict, 'procList':procList}
                        list_graph.append(graphN)
                        i += 1

                
                elif graphIBCD is not None:
                    for g in graphIBCD:
                        graphN={'graphID':{}, 'graphIBCD':g, 'procDict':procDict, 'procList':procList}
                        list_graph.append(graphN)
                        i += 1


            t = Cmpi.readZones(t, FILE, rank=rank)
            t = Cmpi.convert2PartialTree(t, rank=rank)

            
        else: # load 1 fichier par proc

            if graph and exploc ==False :

                # Try to load graph from file
                if os.access('%s/graph.pk'%fileName, os.R_OK):
                    try: import cPickle as pickle
                    except: import pickle
                    file = open('%s/graph.pk'%fileName, 'rb')
                    graphN = pickle.load(file)
                    file.close()
                # Load all skeleton proc files
                else:
                    ret = 1; no = 0; t = []
                    while ret == 1:
                       FILE = '%s_%d.cgns'%(fileName, no)
                       if not os.access(FILE, os.F_OK): ret = 0
                       if ret == 1: 
                          t.append(Cmpi.convertFile2SkeletonTree(FILE))
                       no += 1
                    if t != []:
                        t        = Internal.merge(t)
                        graphID  = Cmpi.computeGraph(t, type='ID'  , reduction=False)
                        graphIBCD= Cmpi.computeGraph(t, type='IBCD', reduction=False)
                        procDict = D2.getProcDict(t)
                        procList = D2.getProcList(t,  sort=True)
                        graphN   = {'graphID':graphID, 'graphIBCD':graphIBCD, 'procDict':procDict, 'procList':procList }
                    else: print('graph non calculable: manque de fichiers connectivite.')


            if graph and exploc : ## dtloc instationnaire

                # Try to load graph from file
                if os.access('%s/graph.pk'%fileName, os.R_OK):
                    try: import cPickle as pickle
                    except: import pickle
                    file = open('%s/graph.pk'%fileName, 'rb')
                    graphN = pickle.load(file)
                    file.close()
                # Load all skeleton proc files
                else:
                    ret = 1; no = 0; t = []
                    while ret == 1:
                       FILE = '%s_%d.cgns'%(fileName, no)
                       if not os.access(FILE, os.F_OK): ret = 0
                       if ret == 1: 
                          t.append(Cmpi.convertFile2SkeletonTree(FILE))
                       no += 1
                    if t != []:
                        t        = Internal.merge(t)
                        list_graph=[]
                        graphID   = Cmpi.computeGraph(t, type='ID'  , reduction=False, exploc=True)
                        graphIBC  = Cmpi.computeGraph(t, type='IBCD', reduction=False, exploc=True)
                        procDic   = D2.getProcDict(t)
                        procLis   = D2.getProcList(t,  sort=True)
              
  
                        i=0
                        graphN={}
                        if graphID is not None:
                            for g in graphID :
                                graphN={'graphID':g, 'graphIBCD':{}, 'procDict':procDict, 'procList':procList}
                                list_graph.append(graphN)
                                i += 1

                
                        elif graphIBCD is not None:
                            for g in graphIBCD :
                                graphN ={'graphID':{}, 'graphIBCD':g, 'procDict':procDict, 'procList':procList}
                                list_graph.append(graphN)
                                i += 1


            FILE = '%s/%s_%d.cgns'%(fileName, baseName, rank)
            if os.access(FILE, os.F_OK): t = C.convertFile2PyTree(FILE)
            else: t = None

    else: # sequential run
        if split == 'single':
            FILE = fileName+'.cgns'
            if os.access(FILE, os.F_OK): t = C.convertFile2PyTree(FILE)
            else: t = None
        else: # multiple
            ret = 1; no = 0; t = []
            while ret == 1:
                FILE = '%s/%s_%d.cgns'%(fileName, baseName, no)
                if not os.access(FILE, os.F_OK): ret = 0
                if ret == 1: t.append(C.convertFile2PyTree(FILE))
                no += 1
            if t != []: t = Internal.merge(t)
            else: t = None

    if cartesian: # peut etre inutile (fait dans convert2File?)
        import Compressor.PyTree as Compressor
        Compressor._uncompressCartesian(t)
        
    if graph and not exploc : return t, graphN
    elif graph and exploc :   return t, list_graph
    else:     return t
    
#==============================================================================
# Load one file
# si split='single': load t.cgns
# si split='multiple': read t.cgns ou t/t_1.cgns
# if graph=True,
# return the communication graph for Chimera and abutting transfers
# and the communication graph for IBM transfers
#==============================================================================
def loadFileG(fileName='t.cgns', split='single', graph=False, mpirun=False):
    """Load tree and connectivity tree."""
    import os.path
    baseName = os.path.basename(fileName)
    baseName = os.path.splitext(baseName)[0] # name without extension
    fileName = os.path.splitext(fileName)[0] # full path without extension

    graphN = {'graphID':None, 'graphIBCD':None, 'procDict':None, 'procList':None}
    list_graph=[]
    if mpirun: # mpi run
        import Distributor2.PyTree as D2
        rank = Cmpi.rank; size = Cmpi.size
        
        if split == 'single':
            FILE = fileName+'.cgns'
            t = Cmpi.convertFile2SkeletonTree(FILE)

            mp = getMaxProc(t)

            if mp+1 != size: 
                raise ValueError('The number of mpi proc (%d) doesn t match the tree distribution (%d)'%(mp+1,size)) 
            if graph:
                graphID   = Cmpi.computeGraph(t, type='ID'  , reduction=False)
                graphIBCD = Cmpi.computeGraph(t, type='IBCD', reduction=False)
                procDict  = D2.getProcDict(t)
                procList  = D2.getProcList(t,  sort=True)
                graphN = {'graphID':graphID, 'graphIBCD':graphIBCD, 'procDict':procDict, 'procList':procList }

                graphID_   = Cmpi.computeGraph(t, type='ID'  , reduction=False, exploc=True)
                graphIBCD_  = Cmpi.computeGraph(t, type='IBCD', reduction=False, exploc=True)
                procDict_  = D2.getProcDict(t)
                procList_  = D2.getProcList(t,  sort=True)
              
                #print(graphID_)
  
                i=0
                graphN_={}
                if graphID_ is not None:
                    for g in graphID_:
                        graphN_={'graphID':g, 'graphIBCD':{}, 'procDict':procDict, 'procList':procList}
                        list_graph.append(graphN_)
                        i += 1

                
                elif graphIBCD_ is not None:
                    for g in graphIBCD_:
                        graphN_={'graphID':{}, 'graphIBCD':g, 'procDict':procDict, 'procList':procList}
                        list_graph.append(graphN_)
                        i += 1
                
                #print(list_graph)
            
            t = Cmpi.readZones(t, FILE, rank=rank)
            t = Cmpi.convert2PartialTree(t, rank=rank)
            
        else: # load 1 fichier par proc
            if graph:
                # Try to load graph from file
                if os.access('%s/graph.pk'%fileName, os.R_OK):
                    try: import cPickle as pickle
                    except: import pickle
                    file = open('%s/graph.pk'%fileName, 'rb')
                    graphN = pickle.load(file)
                    file.close()
                # Load all skeleton proc files
                else:
                    ret = 1; no = 0; t = []
                    while ret == 1:
                       FILE = '%s_%d.cgns'%(fileName, no)
                       if not os.access(FILE, os.F_OK): ret = 0
                       if ret == 1: 
                          t.append(Cmpi.convertFile2SkeletonTree(FILE))
                       no += 1
                    if t != []:
                        t        = Internal.merge(t)
                        graphID  = Cmpi.computeGraph(t, type='ID'  , reduction=False)
                        graphIBCD= Cmpi.computeGraph(t, type='IBCD', reduction=False)
                        procDict = D2.getProcDict(t)
                        procList = D2.getProcList(t,  sort=True)
                        graphN   = {'graphID':graphID, 'graphIBCD':graphIBCD, 'procDict':procDict, 'procList':procList }
                    else: print('graph non calculable: manque de fichiers connectivite.')

            FILE = '%s/%s_%d.cgns'%(fileName, baseName, rank)
            if os.access(FILE, os.F_OK): t = C.convertFile2PyTree(FILE)
            else: t = None

    else: # sequential run
        if split == 'single':
            FILE = fileName+'.cgns'
            if os.access(FILE, os.F_OK): t = C.convertFile2PyTree(FILE)
            else: t = None
        else: # multiple
            ret = 1; no = 0; t = []
            while ret == 1:
                FILE = '%s/%s_%d.cgns'%(fileName, baseName, no)
                if not os.access(FILE, os.F_OK): ret = 0
                if ret == 1: t.append(C.convertFile2PyTree(FILE))
                no += 1
            if t != []: t = Internal.merge(t)
            else: t = None

    if graph: return t, list_graph
    else:     return t


#==============================================================================
# Save tree in one file.
# si split='single': load t.cgns
# si split='multiple': read t.cgns ou t/t_1.cgns
# autorestart possible
# return a partial tree t, a donor partial tree tc, 
# a stat partial tree,
# the communication graph for Chimera and abutting transfers
# the communication graph for IBM transfers
#==============================================================================
def saveFile(t, fileName='restart.cgns', split='single', graph=False, NP=0, 
    mpirun=False, cartesian=False):
    """Save tree in file."""
    import os.path
    baseName = os.path.basename(fileName)
    baseName = os.path.splitext(baseName)[0] # name without extension
    fileName = os.path.splitext(fileName)[0] # full path without extension

    # Rip/add some useless/usefull data (FastS)
    t2 = C.rmVars(t, 'centers:Density_P1')
    C._rmVars(t2, 'centers:VelocityX_P1')
    C._rmVars(t2, 'centers:VelocityY_P1')
    C._rmVars(t2, 'centers:VelocityZ_P1')
    C._rmVars(t2, 'centers:Temperature_P1')
    C._rmVars(t2, 'centers:TurbulentSANuTilde_P1')
    zones = Internal.getZones(t2)
    for z in zones: Internal._rmNodeByPath(z, '.Solver#ownData')
    if cartesian:
        import Compressor.PyTree as Compressor
        Compressor._compressCartesian(t2)
        #Compressor._compressCellN(t2)

    for z in zones:
        node = Internal.getNodeFromName1(z, '.Solver#define')
        if node is not None:
            integ = 'None'; nature = 'None'
            n = Internal.getNodeFromName1(node, 'temporal_scheme')
            if n is not None: integ = Internal.getValue(n)
            n = Internal.getNodeFromName1(node, 'time_step_nature')
            if n is not None: nature = Internal.getValue(n)
            if integ == 'explicit' or nature == 'local':
                C._rmVars(z, 'centers:Density_M1')
                C._rmVars(z, 'centers:VelocityX_M1')
                C._rmVars(z, 'centers:VelocityY_M1')
                C._rmVars(z, 'centers:VelocityZ_M1' )
                C._rmVars(z, 'centers:Temperature_M1')
                C._rmVars(z, 'centers:TurbulentSANuTilde_M1')
    dtloc = Internal.getNodeFromName3(t2, '.Solver#dtloc')
    if dtloc is not None:
        dtloc = dtloc[1]
        node =  Internal.getNodeFromName1(t, 'TimeLevelMotion')
        if node is not None: node[1][0]= dtloc[3]
        else: Internal.createUniqueChild(t, 'TimeLevelMotion', 'DataArray_t', value=dtloc[3])
        node =  Internal.getNodeFromName1(t, 'TimeLevelTarget')
        if node is not None: node[1][0]= dtloc[4]
        else: Internal.createUniqueChild(t, 'TimeLevelTarget', 'DataArray_t', value=dtloc[4])

    # save
    if mpirun: # mpi run
        if split == 'single':  # output in a single file
            FILE = fileName+'.cgns'
            Cmpi.convertPyTree2File(t2, FILE)
        else:
            rank = Cmpi.rank
            if rank == 0: 
                if not os.path.exists(fileName): os.makedirs(fileName)
            Cmpi.barrier()
            FILE = '%s/%s_%d.cgns'%(fileName, baseName, rank)
            C.convertPyTree2File(t2, FILE)

    else: # sequential run
        if split == 'single':
            FILE = fileName+'.cgns'
            C.convertPyTree2File(t2, FILE)
        else:
            if NP == 0: NP = -getMaxProc(t2)-1
            if not os.path.exists(fileName): os.makedirs(fileName)
            for i in range(max(1,-NP)):
                tl = Internal.copyRef(t2)
                bases = Internal.getNodesFromType1(tl, 'CGNSBase_t')
                for b in bases:
                    zones = Internal.getNodesFromType1(b, 'Zone_t')
                    for z in zones:
                        proc = Internal.getNodeFromName2(z, 'proc')
                        if proc is not None:
                            proc = Internal.getValue(proc)
                            if proc != i: Internal._rmNode(b, z)
                C.convertPyTree2File(tl, '%s/%s_%d.cgns'%(fileName,baseName,i))

#==============================================================================
# Load t.cgns (data) tc.cgns (connectivity file) ts.cgns (stat)
# si split='single': load t.cgns
# si split='multiple': read t.cgns ou t/t_1.cgns
# autorestart possible
# return a partial tree t, a donor partial tree tc, 
# a stat partial tree,
# the communication graph for Chimera and abutting transfers
# the communication graph for IBM transfers
# dir is the directory containing files to be read 
#==============================================================================
def loadTree(fileName='t.cgns', split='single', directory='.', graph=False, mpirun=False):
    """Load tree and connectivity tree."""
    import os.path
    import Converter.PyTree as C

    fileName = os.path.splitext(fileName)[0] # full path without extension

    graphN = {'graphID':None, 'graphIBCD':None, 'procDict':None, 'procList':None}
    if mpirun: # mpi run
        import Distributor2.PyTree as D2
        rank = Cmpi.rank; size = Cmpi.size
        
        if split == 'single':
            # Load connect (tc)
            FILE = directory+'/'+fileName+'.cgns'
            t = Cmpi.convertFile2SkeletonTree(FILE)
            if graph:
                graphID   = Cmpi.computeGraph(t, type='ID'  , reduction=False)
                graphIBCD = Cmpi.computeGraph(t, type='IBCD', reduction=False)
                procDict  = D2.getProcDict(t)
                procList  = D2.getProcList(t,  sort=True)
                graphN = {'graphID':graphID, 'graphIBCD':graphIBCD, 'procDict':procDict, 'procList':procList }
            t = Cmpi.readZones(t, FILE, rank=rank)
            zones = Internal.getZones(t)
            t = Cmpi.convert2PartialTree(t, rank=rank)
            zones = Internal.getZones(t)

        else: # load 1 fichier par proc
            if graph:
                # Try to load graph
                if os.access('%s/graph.pk'%directory, os.R_OK):
                    try: import cPickle as pickle
                    except: import pickle
                    file = open('%s/graph.pk'%directory, 'rb')
                    graphN= pickle.load(file)
                    file.close()
                # Load all sqeleton proc files
                else:
                    ret = 1; no = 0; t = []
                    while ret == 1:
                       #FILE = '%s/%s_proc%d.cgns'%(directory, fileName, no)
                       FILE = '%s/%s%d.cgns'%(directory, fileName, no)
                       if not os.access(FILE, os.F_OK): ret = 0
                       if ret == 1: 
                          t.append(Cmpi.convertFile2SkeletonTree(FILE))
                       no += 1
                    if no == size and t != []:
                        t        = Internal.merge(t)
                        graphID  = Cmpi.computeGraph(t, type='ID'  , reduction=False)
                        graphIBCD= Cmpi.computeGraph(t, type='IBCD', reduction=False)
                        procDict = D2.getProcDict(t)
                        procList = D2.getProcList(t,  sort=True)
                        graphN   = {'graphID':graphID, 'graphIBCD':graphIBCD, 'procDict':procDict, 'procList':procList }
                    else: print('graph non calculable: manque de fichiers connectivite.')

            #FILE = '%s/%s_proc%d.cgns'%(directory, fileName, rank)
            FILE = '%s/%s%d.cgns'%(directory, fileName, rank)
            if os.access(FILE, os.F_OK): t = C.convertFile2PyTree(FILE)
            else: t = None

    else: # sequential run
        if split == 'single':
            FILE = fileName+'.cgns'
            if os.access(FILE, os.F_OK): t = C.convertFile2PyTree(FILE)
            else: t = None
        else: # multiple
            # Load connectivity (tc)
            ret = 1; no = 0; tc = []
            while ret == 1:
                #FILE = '%s/%s_proc%d.cgns'%(directory, fileName, no)
                FILE = '%s/%s%d.cgns'%(directory, fileName, no)
                if not os.access(FILE, os.F_OK): ret = 0
                if ret == 1: tc.append(C.convertFile2PyTree(FILE))
                no += 1
            #if no == NP and t != []: t = Internal.merge(t)
            if  t != []: t = Internal.merge(t)
            else: t = None

    if graph: return t, graphN
    else:     return t

#==============================================================================
# Load t.cgns (data) tc.cgns (connectivity file) ts.cgns (stat)
# si split='single': load t.cgns
# si split='multiple': read t.cgns ou t/t_1.cgns
# autorestart possible
# return a partial tree t, a donor partial tree tc, 
# a stat partial tree,
# the communication graph for Chimera and abutting transfers
# the communication graph for IBM transfers
# dir is the directory containing files to be read 
#==============================================================================
def saveTree(t, fileName='restart.cgns', split='single', directory='.', graph=False, mpirun=False):
    """Load tree and connectivity tree."""
    import os.path
    baseName = os.path.basename(fileName)
    baseName = os.path.splitext(baseName)[0] # name without extension
    fileName = os.path.splitext(fileName)[0] # full path without extension
    import Converter.PyTree as C

    # Rip/add some useless/usefull data (FastS)
    import Converter.PyTree as C
    t2 = C.rmVars(t, 'centers:Density_P1')
    C._rmVars(t2, 'centers:VelocityX_P1')
    C._rmVars(t2, 'centers:VelocityY_P1')
    C._rmVars(t2, 'centers:VelocityZ_P1')
    C._rmVars(t2, 'centers:Temperature_P1')
    C._rmVars(t2, 'centers:TurbulentSANuTilde_P1')

    bases  = Internal.getNodesFromType1(t2    , 'CGNSBase_t')       # noeud
    own   = Internal.getNodeFromName1(bases[0], '.Solver#ownData')  # noeud
    dtloc = None
    if own is not None: dtloc = Internal.getNodeFromName1(own, '.Solver#dtloc')[1] # numpy

    zones = Internal.getZones(t2)
    for z in zones:
        node = Internal.getNodeFromName1(z, '.Solver#define')
        if node is not None:
            integ = 'None'; nature = 'None'
            n = Internal.getNodeFromName1(node, 'temporal_scheme')
            if n is not None: integ = Internal.getValue(n)
            n = Internal.getNodeFromName1(node, 'time_step_nature')
            if n is not None: nature = Internal.getValue(n)
            if integ == 'explicit' or nature == 'local':
                C._rmVars(z, 'centers:Density_M1')
                C._rmVars(z, 'centers:VelocityX_M1')
                C._rmVars(z, 'centers:VelocityY_M1')
                C._rmVars(z, 'centers:VelocityZ_M1' )
                C._rmVars(z, 'centers:Temperature_M1')
                C._rmVars(z, 'centers:TurbulentSANuTilde_M1')

    if dtloc is not None:
       node =  Internal.getNodeFromName1(t, 'TimeLevelMotion')
       if node is not None: node[1][0]= dtloc[3]
       else: Internal.createUniqueChild(t, 'TimeLevelMotion', 'DataArray_t', value=dtloc[3])
       node =  Internal.getNodeFromName1(t, 'TimeLevelTarget')
       if node is not None: node[1][0]= dtloc[4]
       else: Internal.createUniqueChild(t, 'TimeLevelTarget', 'DataArray_t', value=dtloc[4])

    # save
    if mpirun: # mpi run
        if split == 'single':  # output in a single file
            FILE = directory+'/'+fileName+'.cgns'
            Cmpi.convertPyTree2File(t2, FILE)
        else:
            rank = Cmpi.rank
            #FILE = '%s/%s_proc%d.cgns'%(directory, fileName, rank)
            FILE = '%s/%s%d.cgns'%(directory, fileName, rank)
            C.convertPyTree2File(t2, FILE)

    else: # sequential run
        FILE = directory+'/'+fileName+'.cgns'
        C.convertPyTree2File(t2, FILE)

def convertpointwise2fast(FILEIN):
    ## This script converst Pointwise meshes to a format that can be used for FastS
    t  = C.convertFile2PyTree(FILEIN+".cgns")
    ##A few comments: Assumes unspecified is empty!!
    for fam in Internal.getNodesFromType(t,'Family_t'):
        if(fam[0] == "Unspecified"):
            Internal.rmNode(t,fam)

    ## The code below was provided by Thomas Renaud. Thank you!_____________________                             
    dicofambc={}                                                                   #|   
    for fam in Internal.getNodesFromType(t,'Family_t'):                            #|       
        fambc = Internal.getValue(Internal.getNodeFromType(fam,'FamilyBC_t'))      #|                 
        dicofambc[fam[0]]=fambc                                                    #|         
                                                                                   #|         
                                                                                   #|           
    for bc in Internal.getNodesFromType(t,'BC_t'):                                 #|             
        if Internal.getValue(bc)=='FamilySpecified':                               #|                  
            famname = Internal.getNodeFromType(bc,'FamilyName_t')                  #|                  
            valbc = dicofambc[Internal.getValue(famname)]                          #|                    
            Internal.setValue(bc,valbc)                                            #|                     
            Internal.rmNode(bc,famname)                                            #|                     
    Internal._rmNodesByType(t,'Family_t')                                          #|               
    ##                                               _______________________________|

    ## Remove further unneccesary information
    ## output is bare bones for FastS
    vars = ['Descriptor_t','FamilyName_t','DataClass_t','DimensionalUnits_t','GridLocation_t','FamilyName']
    for v in vars:
        for fam in Internal.getNodesFromType(t,v):
            Internal.rmNode(t,fam)
    
    C.convertPyTree2File(t ,FILEIN+"_fast.cgns")
    return t
    
def change_name_BC_pointwise(t):
    import sys
    zones = Internal.getNodesFromType2(t, 'Zone_t')
    for z in zones:
        count=0
        zonebc = Internal.getNodesFromType(z, 'BC_t')
        for zbc in zonebc:
            count +=1
            Internal.setName(zbc, Internal.getValue(zbc)+str(count));
    return None


def pointwise2D_2Fast(t):
    import numpy as np
    zones = Internal.getNodesFromType2(t, 'Zone_t')
    for z in zones:
        zonebc = Internal.getNodesFromType(z, 'BC_t')
        for zbc in zonebc:
            s = np.ones((3,2), dtype=np.int32)
            for j in range(0,2):
                for i in range(0,2):
                    s[j][i]=zbc[2][0][1][j][i]
            Internal.setValue(zbc[2][0],s)

        zonegrid = Internal.getNodesFromType(z, 'GridCoordinates_t')
        for zgrid in zonegrid:            
            if Internal.getNodeFromName(zgrid, 'CoordinateZ') is None:
                coordz=Internal.copyNode(Internal.getNodeFromName(zgrid, 'CoordinateX'))
                Internal.setName(coordz, 'CoordinateZ')
                Internal.setValue(coordz,np.zeros((np.shape(coordz[1])[0],np.shape(coordz[1])[1])))
                Internal.addChild(zgrid, coordz, pos=-1)
    return None
