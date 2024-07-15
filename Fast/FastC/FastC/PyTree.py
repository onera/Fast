"""Services for FAST solvers.
"""
from fileinput import filename
import numpy
import os, fnmatch
from . import fastc

FIRST_IT = 0
HOOK     = None

OMP_MODE = 0

# transfered variables
varsN    = ['Density']
varsP    = ['Density_P1']
varsM    = ['Density_M1']

# LBM Variables
varsPLBM   = ['Q1']
#varsMLBM   = ['Qstar_1']
#varsEQLBM  = ['Qeq_1']
#varsNEQLBM = ['Qneq_1']
varsMLBM = ['Q1_M1']
varsS	 = ['Sxx']
varsPSI  = ['corr_xx']
varsMacro= ['Density']

try:
    import Converter.PyTree as C
    import Converter.Internal as Internal
    import Post.PyTree as P
    import Fast.VariablesSharePyTree as VSHARE
    import Post.ExtraVariables2 as PE
    import RigidMotion.PyTree as R
    import math
except:
    raise ImportError("FastC.PyTree: requires Converter, Post and RigidMotion module.")

try: range = xrange
except: pass

try:
    OMP_NUM_THREADS = os.environ['OMP_NUM_THREADS']
    OMP_NUM_THREADS = int(OMP_NUM_THREADS)
except: OMP_NUM_THREADS = 1

import Converter.Mpi as Cmpi
import Distributor2.PyTree as D2

MX_SYNCHRO = 1000
MX_SSZONE  = 10
MX_OMP_SIZE_INT = -1 # set in warmup

#==============================================================================
# Met un dictionnaire de numerics dans une/des zones
#==============================================================================
def setNum2Zones(a, num):
    """Set numeric data dictionary  in zones."""
    ap = Internal.copyRef(a)
    _setNum2Zones(ap, num)
    return ap

def _setNum2Zones(a, num):
    """Set numeric data dictionary  in zones."""
    zones = Internal.getNodesFromType2(a, 'Zone_t')
    for z in zones:
        cont = Internal.createUniqueChild(z, '.Solver#define', 
                                          'UserDefinedData_t')
        for k in num: # some checks?
            Internal.createUniqueChild(cont, k, 'DataArray_t', num[k])
    return None

#==============================================================================
def setNum2Base(a, num):
    """Set numeric data dictionary in bases."""
    ap = Internal.copyRef(a)
    _setNum2Base(ap, num)
    return ap

def _setNum2Base(a, num):
    """Set numeric data dictionary in bases."""
    cont = Internal.createUniqueChild(a, '.Solver#define', 'UserDefinedData_t')
    for k in num: 
         Internal.createUniqueChild(cont, k, 'DataArray_t', num[k])
    return None

#==============================================================================
# Reorder zone pour // omp 
#==============================================================================
def _reorder(t, tc=None):

    if tc is not None: 
       #reordone les bases, sinon souci potentiel en MPi si ordre base != entre proc
       #Internal._sortByName(tc,recursive=False)
       #Internal._sortByName(t, recursive=False)

       #reordone les zones de tc par taille decroissante dans chaque base pour optim openmp
       bases_tc = Internal.getNodesFromType1(tc,'CGNSBase_t')
       for base_tc in bases_tc: 
          zones = Internal.getNodesFromType1(base_tc,'Zone_t')
          #calcul taille de la zone
          size_zone =[]
          for z in zones:
             dim = Internal.getZoneDim(z)
             if dim[0] == 'Structured':
                if dim[3] == 1: kfic = 0
                else          : kfic = 2
                ndimdx = (dim[1]-4)*(dim[2]-4)*(dim[3]-kfic) 
             else: ndimdx = dim[2]
             size_zone.append(ndimdx)

          #Tri les zone par taille decroissante
          new_zones = []
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
          base_tc[2] = orig+new_zones

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
           for zc in zones: 
               zname.append( zc[0] )

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
def createPrimVars(t, omp_mode, rmConsVars=True, Adjoint=False, gradP=False, isWireModel=False, lbmAJ=True):
    """Create primitive vars from conservative vars."""
    tp = Internal.copyRef(t)
    first_iter , vars_zones = _createPrimVars(tp, omp_mode, rmConsVars, Adjoint, gradP=gradP, isWireModel=isWireModel, lbmAJ=lbmAJ)
    return tp, first_iter, vars_zones

#==============================================================================
def _createPrimVars(t, omp_mode, rmConsVars=True, Adjoint=False, gradP=False, isWireModel=False, lbmAJ=True):
    """Create primitive vars from conservative vars."""
    vars_zones=[]
    flag_NSLBM = 0
    bases = Internal.getNodesFromType1(t, 'CGNSBase_t')
    for b in bases:  
        # On regarde si b est une simulation couplee NS-LBM
        model="Unknown"
        a = Internal.getNodeFromName2(b, 'GoverningEquations')
        if a is not None:
           model = Internal.getValue(a)
           if model == 'CouplageNSLBM': flag_NSLBM = 1

        zones = Internal.getNodesFromType1(b, 'Zone_t')
        count = -1
        if omp_mode == 1: count = 0
        for z in zones:
            if omp_mode == 1: count += 1

            #Sauvegarde noeud specifique car delete par init!!! Christophe??
            FA_intext = Internal.getNodeFromPath(z, 'NFaceElements/IntExt')
            if FA_intext is not None:
               FA_indx  =  Internal.getNodeFromPath(z, 'NFaceElements/ElementIndex')
               NG_intext=  Internal.getNodeFromPath(z, 'NGonElements/IntExt')
               NG_indx  =  Internal.getNodeFromPath(z, 'NGonElements/FaceIndex')
               NG_pe    =  Internal.getNodeFromPath(z, 'NGonElements/ParentElements')

            model_loc = 'Unknown'
            a = Internal.getNodeFromName2(z, 'GoverningEquations')
            if a is not None: model_loc = Internal.getValue(a)
            if model_loc == 'Unknown':
               if model != "Unknown" and model != 'CouplageNSLBM': model_loc = model
               else: 
                  model_loc = 'Euler'; print("WARNING: definition model: nothing found in base or zone: Euler imposed.")

            sa, lbmflag, neq_lbm, FIRST_IT = _createVarsFast(b, z, omp_mode, rmConsVars, Adjoint, gradP=gradP, isWireModel=isWireModel, flag_coupNSLBM=flag_NSLBM, lbmAJ=lbmAJ)
            
            if FA_intext is not None:
               Internal.getNodeFromPath(z, 'NFaceElements')[2] +=[FA_indx, FA_intext] 
               Internal.getNodeFromPath(z, 'NGonElements' )[2] +=[NG_pe, NG_indx, NG_intext] 
            
            #recuperation option de calcul
            define = Internal.getNodeFromName1(z, '.Solver#define')
            sponge = 0
            a = Internal.getNodeFromName1(define, 'lbm_sponge')
            if a is not None: sponge = Internal.getValue(a)
            source = 0
            a = Internal.getNodeFromName1(define, 'source')
            if a is not None: source = Internal.getValue(a)
            sfd = 0
            a = Internal.getNodeFromName1(define, 'sfd')
            if a is not None: sfd = Internal.getValue(a)
            extract_res = 0
            a = Internal.getNodeFromName1(define, 'extract_res')
            if a is not None: extract_res = Internal.getValue(a)
            motion = None 
            a = Internal.getNodeFromName1(define, 'motion')
            if a is not None: motion = Internal.getValue(a)
            show_senseur = 0
            a = Internal.getNodeFromName2(z, 'show_senseur')
            if a is not None: show_senseur = Internal.getValue(a)

            vars_p=['Density', 'VelocityX', 'VelocityY','VelocityZ', 'Temperature']
            vars_c=['Density', 'MomentumX', 'MomentumY','MomentumZ', 'EnergyStagnationDensity']
            vars_s=['xx','xy','xz','yy','yz','zz']
            vars_psi=['xx','yy','zz','x','y','z']
            vars_showsenseur = ['senseur_i','senseur_j','senseur_k']

            if sa: 
               vars_p.append('TurbulentSANuTilde')
               vars_c.append('TurbulentSANuTildeDensity')

            vars=[]
            if lbmflag:
                fields2compact=[]
                for i in range(1,neq_lbm+1):
                   fields2compact.append('centers:Q'+str(i))
                vars.append(fields2compact)
                fields2compact=[]
                for i in range(1,neq_lbm+1):
                  fields2compact.append('centers:Q'+str(i)+'_M1')
                vars.append(fields2compact)
                fields2compact=[]
                for i in range(len(vars_s)):
                  fields2compact.append('centers:S'+str(vars_s[i])) #Tenseur S pour HRR
                vars.append(fields2compact)
                for i in range(len(vars_psi)):
                   fields2compact.append('centers:corr_'+str(vars_psi[i])) #Gradients PSI pour HRR
                vars.append(fields2compact)

            if lbmAJ:
                fields2compact=[]
                for i in range(1,3+1):
                  fields2compact.append('centers:cellN_IBC_LBM_'+str(i))
                vars.append(fields2compact)

                vars.append(['centers:SpongeCoef'])

            #Compactage des variables dans le cas du couplage
            if flag_NSLBM == 1:
               fields2compact=[]
               for i in range(len(vars_s)):
                  fields2compact.append('centers:S'+str(vars_s[i]))
               vars.append(fields2compact)

            if model_loc == 'Euler' or model_loc == 'NSLaminar' or model_loc == 'NSTurbulent':
              timelevel = ['', '_M1','_P1']
            elif model_loc == 'LBMLaminar':
              if lbmAJ:  timelevel = ['', '_M1']
              else: timelevel = ['']
            else: raise ValueError('createPrimVars: unknown model %s.'%model_loc)

            for level in timelevel:  #champs primitives
               fields2compact=[]
               for v in vars_p:
                 fields2compact.append('centers:'+v+level)
               vars.append(fields2compact)
            if source == 1:                   #terme source volumique
               fields2compact=[]
               for v in vars_c:
                 fields2compact.append('centers:'+v+'_src')
               vars.append(fields2compact)
            if source == 2:          #terme source volumique
               fields2compact=[]
               for v in vars_c:
                 fields2compact.append('centers:'+v+'_src')
               fields2compact.append('centers:cellN_src')
               vars.append(fields2compact)
            if motion == 'deformation' or motion == 'rigid_ext':     #ale deformable
               fields2compact=[]
               for v in ['VelocityX','VelocityY','VelocityZ']:
                 fields2compact.append('Motion:'+v)
               vars.append(fields2compact)
            if sfd == 1:                       #sfd
               fields2compact=[]
               for v in vars_p:
                 fields2compact.append('centers:'+v+'_f')
               vars.append(fields2compact)
            if extract_res == 1:
               fields2compact=[]
               for v in vars_c:
                 fields2compact.append('centers:Res_'+v)
               vars.append(fields2compact)
            if show_senseur == 1:
               fields2compact=[]
               for v in vars_showsenseur:
                 fields2compact.append('centers:'+v)
               vars.append(fields2compact)

            if gradP:
              fields2compact=[]
              for v in ['Density', 'Temperature']:
                fields2compact.append('centers:'+'gradx'+v)
                fields2compact.append('centers:'+'grady'+v)
                fields2compact.append('centers:'+'gradz'+v)
              vars.append(fields2compact)

            if isWireModel:
              fields2compact = []
              for v in vars_p:
                fields2compact.append('centers:'+v+'_WM')
              vars.append(fields2compact)
    
            # on compacte les variables "visqueuse"
            loc       ='centers:'
            fields    = []
            #for v in ['ViscosityEddy','TurbulentDistance', 'zgris', 'ViscosityEddyCorrection', 'sgsCorrection']: fields.append(loc+v)
            for v in ['ViscosityEddy','TurbulentDistance', 'zgris', 'ViscosityEddyCorrection']: fields.append(loc+v)
            for sufix in ['0','N']:
               for i in range(1,len(vars_p)+1):
                  fields.append(loc+'drodm'+sufix+'_'+str(i))

            fields2compact = []
            for field in fields: 
               if C.isNamePresent(z, field) == 1:fields2compact.append(field)
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

            # on compacte Pression
            fields2compact = []
            if  C.isNamePresent(z, loc+'Pressure') == 1: fields2compact.append(loc+'Pressure')
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
def _createVarsFast(base, zone, omp_mode, rmConsVars=True, adjoint=False, gradP=False,isWireModel=False, flag_coupNSLBM=0, lbmAJ=True):
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

    gamma = adim[11]
    rgp = (gamma-1.)*adim[7]
    
    t0 = timeit.default_timer()
    if C.isNamePresent(zone, 'centers:VelocityX'  ) != 1: P._computeVariables(zone, ['centers:VelocityX'  ], rgp=rgp)
    if C.isNamePresent(zone, 'centers:VelocityY'  ) != 1: P._computeVariables(zone, ['centers:VelocityY'  ], rgp=rgp)
    if C.isNamePresent(zone, 'centers:VelocityZ'  ) != 1: P._computeVariables(zone, ['centers:VelocityZ'  ], rgp=rgp)
    if C.isNamePresent(zone, 'centers:Temperature') != 1: P._computeVariables(zone, ['centers:Temperature'], rgp=rgp, gamma=gamma)

    if sa and C.isNamePresent(zone, 'centers:TurbulentSANuTilde') != 1: 
        if C.isNamePresent(zone, 'centers:TurbulentSANuTildeDensity') != 1: C._initVars(zone, '{centers:TurbulentSANuTilde}= %20.16g/{centers:Density}'%adim[14])
        else: C._initVars(zone, '{centers:TurbulentSANuTilde}= {centers:TurbulentSANuTildeDensity}/{centers:Density}')

    if rmConsVars:
        C._rmVars(zone, 'centers:MomentumX')
        C._rmVars(zone, 'centers:MomentumY')
        C._rmVars(zone, 'centers:MomentumZ')
        C._rmVars(zone, 'centers:EnergyStagnationDensity')
        C._rmVars(zone, 'centers:TurbulentSANuTildeDensity')


    FIRST_IT = 1
    #===========================================================================
    # ZONE NS
    #===========================================================================
    if model != 'LBMLaminar':
       neq_lbm = 0
       #on test s'il existe 2 niveaux en temps dans l'arbre pour appliquer la bonne formule de derivee temporelle a la premere iteration
       if C.isNamePresent(zone, 'centers:Density_M1')     != 1: C._cpVars(zone, 'centers:Density'    , zone, 'centers:Density_M1')  ; FIRST_IT=0
       if C.isNamePresent(zone, 'centers:VelocityX_M1')   != 1: C._cpVars(zone, 'centers:VelocityX'  , zone, 'centers:VelocityX_M1'); FIRST_IT=0
       if C.isNamePresent(zone, 'centers:VelocityY_M1')   != 1: C._cpVars(zone, 'centers:VelocityY'  , zone, 'centers:VelocityY_M1'); FIRST_IT=0
       if C.isNamePresent(zone, 'centers:VelocityZ_M1')   != 1: C._cpVars(zone, 'centers:VelocityZ'  , zone, 'centers:VelocityZ_M1'); FIRST_IT=0
       if C.isNamePresent(zone, 'centers:Temperature_M1') != 1: C._cpVars(zone, 'centers:Temperature', zone, 'centers:Temperature_M1'); FIRST_IT=0
       if sa and C.isNamePresent(zone, 'centers:TurbulentSANuTilde_M1') != 1: C._cpVars(zone, 'centers:TurbulentSANuTilde', zone, 'centers:TurbulentSANuTilde_M1'); FIRST_IT=0

       if C.isNamePresent(zone, 'centers:Density_P1')   != 1: C._cpVars(zone, 'centers:Density'  , zone, 'centers:Density_P1')  
       if C.isNamePresent(zone, 'centers:VelocityX_P1') != 1: C._cpVars(zone, 'centers:VelocityX', zone, 'centers:VelocityX_P1')
       if C.isNamePresent(zone, 'centers:VelocityY_P1') != 1: C._cpVars(zone, 'centers:VelocityY', zone, 'centers:VelocityY_P1')
       if C.isNamePresent(zone, 'centers:VelocityZ_P1') != 1: C._cpVars(zone, 'centers:VelocityZ', zone, 'centers:VelocityZ_P1')
       if C.isNamePresent(zone, 'centers:Temperature_P1') != 1: C._cpVars(zone, 'centers:Temperature', zone, 'centers:Temperature_P1')
       if sa and C.isNamePresent(zone, 'centers:TurbulentSANuTilde_P1') != 1: C._cpVars(zone, 'centers:TurbulentSANuTilde', zone, 'centers:TurbulentSANuTilde_P1')

       if gradP and C.isNamePresent(zone, 'centers:gradxDensity') != 1:
           C._initVars(zone, 'centers:gradxDensity', 1e-15)
           C._initVars(zone, 'centers:gradyDensity', 1e-15)
           C._initVars(zone, 'centers:gradzDensity', 1e-15)

           C._initVars(zone, 'centers:gradxTemperature', 1e-15)
           C._initVars(zone, 'centers:gradyTemperature', 1e-15)
           C._initVars(zone, 'centers:gradzTemperature', 1e-15)

       if isWireModel:
           if C.isNamePresent(zone, 'centers:Density_WM') != 1:
               C._initVars(zone, 'centers:Density_WM'    , 1e-15)
               C._initVars(zone, 'centers:VelocityX_WM'  , 1e-15)
               C._initVars(zone, 'centers:VelocityY_WM'  , 1e-15)
               C._initVars(zone, 'centers:VelocityZ_WM'  , 1e-15)
               C._initVars(zone, 'centers:Temperature_WM', 1e-15)
               if (sa): C._initVars(zone, 'centers:TurbulentSANuTilde_WM', 1e-15)

       #Pour le couplage NS-LBM on calcule S
       if flag_coupNSLBM:
          if C.isNamePresent(zone, 'centers:Sxx') != 1: C._initVars(zone,'centers:Sxx',0.)
          if C.isNamePresent(zone, 'centers:Sxy') != 1: C._initVars(zone,'centers:Sxy',0.)
          if C.isNamePresent(zone, 'centers:Sxz') != 1: C._initVars(zone,'centers:Sxz',0.)
          if C.isNamePresent(zone, 'centers:Syy') != 1: C._initVars(zone,'centers:Syy',0.)
          if C.isNamePresent(zone, 'centers:Syz') != 1: C._initVars(zone,'centers:Syz',0.)
          if C.isNamePresent(zone, 'centers:Szz') != 1: C._initVars(zone,'centers:Szz',0.)
    # fin du if zone NS ------------------------------------
    #===========================================================================
    # ZONE LBM
    #===========================================================================
    else:
        neq_lbm = Internal.getNodeFromName2(zone, 'Parameter_int')[1][VSHARE.NEQ_LBM]
        sponge = Internal.getNodeFromName2(zone, 'Parameter_int')[1][VSHARE.LBM_SPONGE]
        
        #On cree les fonctions de distribution
        for i in range(1,neq_lbm+1):
            if C.isNamePresent(zone, 'centers:Q'+str(i))       != 1: C._initVars(zone, 'centers:Q'+str(i), 0.)
            if C.isNamePresent(zone, 'centers:Q'+str(i)+'_M1') != 1: C._initVars(zone, 'centers:Q'+str(i)+'_M1', 0.)
            #var AJ
            #if C.isNamePresent(zone, 'centers:Qstar_'+str(i)      ) != 1: C._initVars(zone, 'centers:Qstar_'+str(i)      , 0.)  /Q1_M1
            #if C.isNamePresent(zone, 'centers:Qeq_'  +str(i)      ) != 1: C._initVars(zone, 'centers:Qeq_'  +str(i)      , 0.)  /drodm
            #if C.isNamePresent(zone, 'centers:Qneq_' +str(i)      ) != 1: C._initVars(zone, 'centers:Qneq_' +str(i)      , 0.)  /drodm2
        
        #Mettre un flag si l'on veut le HRR ? Sinon le tenseur S et psi ne sont pas necessaires.
        #Tenseur grad des vitesses pour le HRR
        COMPOSANTES_S = ['xx','xy','xz','yy','yz','zz'] #en  3D 9 composantes mais que 6 ici car symetrique (xy=yx, xz=zx et yz=zy)
        for comps in COMPOSANTES_S:
            if C.isNamePresent(zone, 'centers:S'+str(comps)) != 1: C._initVars(zone, 'centers:S'+str(comps), 0.)
        #Gradiens pour le PSI
        COMPOSANTES_PSI = ['xx','yy','zz','x','y','z']
        for comppsi in COMPOSANTES_PSI:
            if C.isNamePresent(zone, 'centers:corr_'+str(comppsi)) != 1: C._initVars(zone, 'centers:corr_'+str(comppsi), 0.)
        
        #Ajout du coef des zones eponges si necessaire
        if sponge == 1:
            if C.isNamePresent(zone, 'centers:SpongeCoef') != 1: C._initVars(zone, 'centers:SpongeCoef', 0.)

        #var AJ
        if lbmAJ:
          for i in range(1,3+1):
              if C.isNamePresent(zone, 'centers:cellN_IBC_LBM'+str(i)) != 1: C._initVars(zone, 'centers:cellN_IBC_LBM_'+str(i), 0.)
            
          if C.isNamePresent(zone, 'centers:SpongeCoef')    != 1: C._initVars(zone, 'centers:SpongeCoef'    , 0.)

          if C.isNamePresent(zone, 'centers:Density_M1'    ) != 1: C._cpVars(zone, 'centers:Density'    , zone, 'centers:Density_M1'    ); FIRST_IT=0
          if C.isNamePresent(zone, 'centers:VelocityX_M1'  ) != 1: C._cpVars(zone, 'centers:VelocityX'  , zone, 'centers:VelocityX_M1'  ); FIRST_IT=0
          if C.isNamePresent(zone, 'centers:VelocityY_M1'  ) != 1: C._cpVars(zone, 'centers:VelocityY'  , zone, 'centers:VelocityY_M1'  ); FIRST_IT=0
          if C.isNamePresent(zone, 'centers:VelocityZ_M1'  ) != 1: C._cpVars(zone, 'centers:VelocityZ'  , zone, 'centers:VelocityZ_M1'  ); FIRST_IT=0
          if C.isNamePresent(zone, 'centers:Temperature_M1') != 1: C._cpVars(zone, 'centers:Temperature', zone, 'centers:Temperature_M1'); FIRST_IT=0
    # fin du if zone LBM ------------------------------------

    # init moyenne loi de paroi
    bcs   = Internal.getNodesFromType2(zone, 'BC_t')
    for bc in bcs:
      btype = Internal.getValue(bc)
      if btype == 'BCWallExchange':
         sol     = Internal.getNodeFromName1(zone, 'FlowSolution#Centers')
         ptrange = Internal.getNodesFromType1(bc, 'IndexRange_t')
         rg      = ptrange[0][1]
         sz      = max(1, rg[0,1]-rg[0,0] ) * max(1, rg[1,1]-rg[1,0] ) * max(1, rg[2,1]-rg[2,0] )
         Prop    = Internal.getNodeFromName(bc,'.Solver#Property')
         wmles   = Internal.getNodeFromName(Prop,'WMLES_parameter')[1]
         Nechant = wmles[1]
         print('Nechant',Nechant) 
         vars=['Density','VelocityX','VelocityY','VelocityZ','Temperature']
         if Nechant==0:
            wmles[1]=1  #nb echantillon initialise
            wmles[2]=1  #position dernier echantillon calculer
            moy  = Internal.getNodeFromName1(Prop, "AvgPlane-Primitive" )[1]
            shift=0
            for var in vars:
               inst = Internal.getNodeFromName1(sol, var )[1]
               #CL I
               if rg[0][1] - rg[0][0]==0:
                 if rg[0][1]==1: ijk_tg =5;     
                 else:           ijk_tg =rg[0][1]-7
                 l = shift
                 for k in range(rg[2][0]-1,rg[2][1]):
                    for j in range(rg[1][0]-1,rg[1][1]):
                      moy[l]    = inst[ijk_tg ,j,k]   #moyenne
                      moy[l+sz] = inst[ijk_tg ,j,k]   #echant 1
                      l+=1
               #CL J
               if rg[1][1] - rg[1][0]==0:
                 if rg[1][1]==1: ijk_tg =5
                 else:           ijk_tg = rg[1][1]-7
                 l = shift
                 #print("shift",shift, sz*5, var, rg[2][0]-1,rg[2][1]-1 )
                 for k in range(rg[2][0]-1,rg[2][1]-1):
                    for i in range(rg[0][0]-1,rg[0][1]-1):
                      moy[l]      = inst[i, ijk_tg ,k]   #moyenne
                      moy[l+sz*5] = inst[i, ijk_tg ,k]   #echant 1
                      #if zone[0]=='cart.0' and var== "Temperature": print l+sz*5,i 
                      l+=1
               #CL K
               if numpy.shape(inst)[2]!=1:
                  if rg[2][1] - rg[2][0]==0:
                    if rg[2][1]==1: ijk_tg =5
                    else:           ijk_tg = rg[2][1]-7
                    l = shift
                    for j in range(rg[1][0]-1,rg[1][1]):
                      for i in range(rg[0][0]-1,rg[0][1]):
                        moy[l]    = inst[i,j, ijk_tg ]   #moyenne
                        moy[l+sz] = inst[i,j, ijk_tg ]   #echant 1
                        l+=1
               shift+=sz


    ale=0
    b = Internal.getNodeFromName1(zone, '.Solver#define')
    if b is not None:
       a = Internal.getNodeFromName1(b, 'motion')
       if a is not None: ale = Internal.getValue(a)
       if ale =='deformation' or ale == 'rigid_ext':
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
            vars = Internal.getNodesFromType1(flowsol, 'DataArray_t')
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
        kry = numpy.empty(size, dtype=numpy.float64)

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

    define = Internal.getNodeFromName1(zone, '.Solver#define')
    # init termes sources
    source = 0
    a = Internal.getNodeFromName1(define,'source')
    if a is not None: source = Internal.getValue(a)
    if source >= 1:
       if C.isNamePresent(zone, 'centers:Density_src')   != 1:   C._initVars(zone, 'centers:Density_src', 0.)
       if C.isNamePresent(zone, 'centers:MomentumX_src') != 1:   C._initVars(zone, 'centers:MomentumX_src', 0.)
       if C.isNamePresent(zone, 'centers:MomentumY_src') != 1:   C._initVars(zone, 'centers:MomentumY_src', 0.)
       if C.isNamePresent(zone, 'centers:MomentumZ_src') != 1:   C._initVars(zone, 'centers:MomentumZ_src', 0.)
       if C.isNamePresent(zone, 'centers:EnergyStagnationDensity_src') != 1: C._initVars(zone, 'centers:EnergyStagnationDensity_src', 0.)
       if (sa and C.isNamePresent(zone, 'centers:TurbulentSANuTildeDensity_src') != 1): C._initVars(zone, 'centers:TurbulentSANuTildeDensity_src', 0.)
       if source == 2 and C.isNamePresent(zone, 'centers:cellN_src') != 1:  C._initVars(zone, 'centers:cellN_src', 0.)

    # init termes zone eponge
    sponge = 0
    a = Internal.getNodeFromName1(define,'lbm_sponge')
    if a is not None: sponge = Internal.getValue(a)
    if sponge == 1:
       if C.isNamePresent(zone, 'centers:ViscosityEddyCorrection') != 1: C._initVars(zone, 'centers:ViscosityEddyCorrection', 1.)
       '''
       sgsmodel='Miles'
       a = Internal.getNodeFromName1(define,'sgsmodel')
       if a is not None: sgsmodel = Internal.getValue(a)
       if sgsmodel =='smsm' or sgsmodel=='msm':
          if C.isNamePresent(zone, 'centers:sgsCorrection') != 1: C._initVars(zone, 'centers:sgsCorrection', 1.)
       '''

    return sa, lbm, neq_lbm, FIRST_IT

#==============================================================================
# Construit les donnees compactees pour traiter les verrou et plage omp
#==============================================================================
def _build_omp(t):
    # Data for each Base
    bases = Internal.getNodesFromType1(t, 'CGNSBase_t')

    #dimensionnememnt tableau
    size_int       = 0
    for b in bases:
        zones = Internal.getZones(b)
        for z in zones:
            
            o    = Internal.getNodeFromName1( z, '.Solver#ownData')
            dims = Internal.getZoneDim(z)

            #on concatene les donnes omp dans param_int
            param_int = Internal.getNodeFromName1(o, 'Parameter_int')
            size = numpy.shape(param_int[1])
            c = 1
            for s in size: c=c*s

            size_scheduler = 0
            if dims[3] == 'NGON':
               CellScheduler = Internal.getNodeFromName1(z, 'CellScheduler')[1]
               size_scheduler += numpy.size(CellScheduler)
               Nbr_BlkIntra = Internal.getNodeFromName1(z, 'Nbr_BlkIntra')[1]
               size_scheduler += numpy.size(Nbr_BlkIntra)

            #print("SIZE int", size_int,"c=", c,"size_sched=",  size_scheduler)
            datap = numpy.zeros(size_int + c + size_scheduler, Internal.E_NpyInt)

            datap[0:c]   = param_int[1][0:c]
            datap[69]    = c                   # debut tableau omp dans param_int
            datap[VSHARE.SCHEDULER] = c + size_int        # debut cellscheduler omp dans param_int
            param_int[1] = datap

            if dims[3] == 'NGON':
               deb = param_int[1][VSHARE.SCHEDULER] 
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
                       #print("deb+c", deb+c, c, l,k,j,i,param_int[1][VSHARE.SCHEDULER] )
                       param_int[1][ deb +c ]= CellScheduler[ i, j, k, l]
                       c +=1
               CellScheduler = param_int[1][ deb:deb+size_scheduler ]

    return None

#==============================================================================
# Construit les datas possedees par Fast
#==============================================================================
def _buildOwnData(t, Padding):

    global  MX_OMP_SIZE_INT
    # Data for each Base

    # Available keys for bases and zones
    # 0: requires an int, 1: requires a float, 2: requires any string, 
    # 3: requires array/list of ints, 4: requires array/list of floats,
    # []: requires given strings
    keys4Base = {
    'temporal_scheme':['explicit', 'implicit', 'implicit_local', 'explicit_local'],
    'ss_iteration':0,
    'rk':0, 
    'modulo_verif':0,
    'exp_local':0,
    'time_begin_ale':1,
    'omp_mode':0,
    'explicit_local_type':0    ## 0:explicit local non conservatif,  1:explicit local conservatif (correction du bilan de flux aux interface)  
    }
    
    keys4Zone = {
    'scheme':['ausmpred', 'senseur', 'roe_min', 'roe', 'roe_nul', 'roe_kap', 'senseur_hyper'],
    'implicit_solver':['lussor', 'gmres'],
    'nb_relax':0,
    'nb_krylov':0,
    'nb_restart':0,
    'motion':['none', 'rigid', 'rigid_ext', 'deformation'],
    'rotation':4,
    'time_step':1,
    'io_thread':0,
    'sgsmodel': ['smsm','msm','Miles'],
    'wallmodel': ['musker','power'],
    'wallmodel_sample': 0,
    'ransmodel': ['SA', 'SA_comp', 'SA_diff'],
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
    'coef_hyper':4, 
    'prandtltb':1, 
    'sfd':0, 
    'sfd_chi':1, 
    'sfd_delta':1, 
    'sfd_init_iter':0, 
    'nudging_ampli':1, 
    'nudging_vector':4, 
    'slope':["o1", "o3","o5", "minmod","o3sc","o5sc"],
    'DES':["zdes1", "zdes1_w", "zdes2", "zdes2_w", "zdes3", "zdes3_w"],
    'SA_add_LowRe': 0,
    'SA_add_RotCorr': 0,
    'snear': 1, # ignored
    'DES_debug':['none','active'],
    'extract_res':0,
    'IBC':3,
    'source':0,
    'Cups':4,
    'senseurType':0,
    'ratiom':1, 
    #=========================================================
    # LBM specific keywords
    #=========================================================
    'LBM_velocity_set':['D3Q19','D3Q27'],
    'LBM_coll_model':['BGK', 'BGK+', 'RR', 'HRR', 'TRT'], #Remarque : TRT uniquement pour FastLBM (AJ)
    'LBM_relax_time':1,
    'LBM_hrr_sigma':1,
    'LBM_compressible':0,
    'LBM_sponge':0,
    'lbm_ns':0,
    'lbm_c0':1,
    'lbm_gamma_precon':1,    
    'lbm_dif_coef':1,
    'lbm_selective_filter':0,
    'lbm_selective_filter_size':0,
    'lbm_selective_filter_sigma':1,
    'lbm_adaptive_filter_chi':1,
    'lbm_chi_spongetypeII':1,
    'lbm_sponge_size':0,
    'lbm_spng_xmin':1,
    'lbm_spng_xmax':0,
    'lbm_spng_ymin':0,
    'lbm_spng_ymax':0,
    'lbm_spng_zmin':0,
    'lbm_spng_zmax':0,
    'lbm_ibm':0,
    'lbm_isforce':0,
    'lbm_zlim':1,
    'lbm_ibm_connector':0,
    'LBM_overset':0,
    'LBM_dx':1,
    'LBM_NS':0,
    #=========================================================
    #========================================================= 
    'KWire_p':1,
    'DiameterWire_p':1,
    'CtWire_p':1    
    }

    bases = Internal.getNodesFromType1(t, 'CGNSBase_t')  # noeud

    # Checks if there is at least one LBM zone (or base) in tree
    flaglbm=False
    for base in bases:
      # First check at base level
      model    = 'NSLaminar'
      a        = Internal.getNodeFromName2(base, 'GoverningEquations')
      if a is not None: model = Internal.getValue(a)
      if model == 'LBMLaminar': 
          flaglbm=True
          break
      # Then check at zone level
      zones= Internal.getZones(base)
      for z in zones:
        model_z = None
        a = Internal.getNodeFromName2(z, 'GoverningEquations')
        if a is not None: model_z = Internal.getValue(a)
        if model_z == 'LBMLaminar': 
            flaglbm=True
            break

    #init time et No iteration
    it0 =0; temps=0.; timelevel_motion= 0; timelevel_target= 0;LBMCycleIteration=0
    first = Internal.getNodeFromName1(t, 'Iteration')
    if first is not None: it0 = Internal.getValue(first)
    first = Internal.getNodeFromName1(t, 'Time')
    if first is not None: temps = Internal.getValue(first)
    first = Internal.getNodeFromName1(t, 'TimeLevelMotion')
    if first is not None: timelevel_motion = Internal.getValue(first)
    first = Internal.getNodeFromName1(t, 'TimeLevelTarget')
    if first is not None: timelevel_target = Internal.getValue(first)
    first = Internal.getNodeFromName1(t, 'LBMCycleIteration')
    if first is not None: LBMCycleIteration = Internal.getValue(first)
    else:
       if flaglbm: Internal.createUniqueChild(t, 'LBMCycleIteration', 'DataArray_t', value=0)

    MafzalMode = 3
    first = Internal.getNodeFromName1(t, 'MafzalMode')
    if first is not None: MafzalMode = Internal.getValue(first)

    AlphaGradP = 1
    first = Internal.getNodeFromName1(t, 'AlphaGradP')
    if first is not None: AlphaGradP = Internal.getValue(first)

    NbptsLinelets = 0
    first = Internal.getNodeFromName1(t, 'NbptsLinelets')
    if first is not None: NbptsLinelets = Internal.getValue(first)

    zones = Internal.getZones(t)

    val=1; i=0
    veclevel = []; mod=""; posg=[] 
 
    # Recherche des niveaux en temps des differentes zones
    d = Internal.getNodeFromName1(t, '.Solver#define')
    a = Internal.getNodeFromName1(d, 'temporal_scheme')
    if a is not None:
       temporal_scheme = Internal.getValue(a)
       if temporal_scheme == 'explicit_local':
         for z in zones:
            dim = Internal.getZoneDim(z)
            i = dim[1]//2; j = dim[2]//2; k = dim[3]//2
            level = C.getValue(z, 'centers:niveaux_temps', (i,j,k))
            veclevel.append(int(level))
       else:
         for z in zones: 
           d = Internal.getNodeFromName1(z, '.Solver#define')
           if d is not None:
               checkKeys(d, keys4Zone)
               a = Internal.getNodeFromName1(d, 'niveaux_temps')
               if a is not None: veclevel.append( Internal.getValue(a) )
               else: veclevel.append(1)
    else:
      for z in zones: 
        d = Internal.getNodeFromName1(z, '.Solver#define')
        if d is not None:
            checkKeys(d, keys4Zone)
            a = Internal.getNodeFromName1(d, 'niveaux_temps')
            if a is not None: veclevel.append( Internal.getValue(a) )
            else: veclevel.append(1)
        else:  veclevel.append(1)

    maxlevel = max(veclevel)


    #determine les parametre globaux valable pour toue les zones et bases
    d = Internal.getNodeFromName1(t, '.Solver#define')

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
    ompmode         = 0

    if d is not None:
        checkKeys(d, keys4Base)
        a = Internal.getNodeFromName1(d, 'omp_mode')
        if a is not None: ompmode = Internal.getValue(a)
        ompmode = max(ompmode, 0); ompmode = min(ompmode, 1)
            
        a = Internal.getNodeFromName1(d, 'temporal_scheme')
        if a is not None: temporal_scheme = Internal.getValue(a)
        a = Internal.getNodeFromName1(d, 'ss_iteration')
        if a is not None: ss_iteration = Internal.getValue(a)
        if temporal_scheme == "implicit_local": modulo_verif = 7 # change default
        a = Internal.getNodeFromName1(d, 'modulo_verif')
        if a is not None: modulo_verif = Internal.getValue(a)
        a = Internal.getNodeFromName1(d, 'restart_fields')
        if a is not None: restart_fields = Internal.getValue(a)
        a = Internal.getNodeFromName1(d, 'rk')
        if a is not None: rk = Internal.getValue(a)
        if temporal_scheme == "implicit" or temporal_scheme =="implicit_local": rk=3
        a = Internal.getNodeFromName1(d, 'exp_local')
        if a is not None: exploc = Internal.getValue(a)
        if temporal_scheme == "implicit" or temporal_scheme == "implicit_local": exploc=0
        a = Internal.getNodeFromName1(d, 'explicit_local_type')         
        if a is not None: exploctype = Internal.getValue(a)
        a = Internal.getNodeFromName1(d, 'time_begin_ale')
        if a is not None: t_init_ale = Internal.getValue(a)

           
    # Base ownData (generated)
    #o = Internal.createUniqueChild(b, '.Solver#ownData', 'UserDefinedData_t')
    o = Internal.createUniqueChild(t, '.Solver#ownData', 'UserDefinedData_t')

    if temporal_scheme == "explicit": nssiter = 3
    elif temporal_scheme == "explicit_lbm": nssiter = 1
    elif temporal_scheme == "implicit": nssiter = ss_iteration+1
    elif temporal_scheme == "implicit_local": nssiter = ss_iteration+1
    elif temporal_scheme == "explicit_local": # explicit local instationnaire ordre 3
        nssiter = 4*maxlevel 
        rk = 3
        exploc = 1
    else: print('Warning: Fast: invalid value %s for key temporal_scheme.'%temporal_scheme)
    try: ss_iteration = int(ss_iteration)
    except: print('Warning: Fast: invalid value %s for key ss_iteration.'%ss_iteration)
    try: modulo_verif = int(modulo_verif)
    except: print('Warning: Fast: invalid value %s for key modulo_verif.'%modulo_verif)
    if rk == 1 and exploc == 0 and temporal_scheme == "explicit": nssiter = 1 # explicit global
    if rk == 2 and exploc == 0 and temporal_scheme == "explicit": nssiter = 2 # explicit global
    if rk == 3 and exploc == 0 and temporal_scheme == "explicit": nssiter = 3 # explicit global

    if MX_OMP_SIZE_INT == -1:
       #print("resize MX_OMP_SIZE_INT dans buildOwndata")
       ntask = len(zones)*2  # 4 sous-zones max par zone en moyenne 
       size_task =  (6 + 7*OMP_NUM_THREADS)*ntask
       size_omp  =  1 + nssiter*(2 + size_task)
       MX_OMP_SIZE_INT = size_omp
    size_omp = MX_OMP_SIZE_INT
    
    nitCyclLBM = 2**(maxlevel-1)
    dtdim = 12 + nitCyclLBM
    #print('ncycl_LBM', nitCyclLBM, maxlevel)

    datap = numpy.zeros((dtdim+ size_omp), Internal.E_NpyInt)
    datap[0] = nssiter 
    datap[1] = modulo_verif
    datap[2] = restart_fields-1
    datap[3] = timelevel_motion
    datap[4] = timelevel_target 
    datap[5] = timelevel_prfile
    datap[6] = rk
    datap[7] = it0
    datap[8] = ompmode
    datap[9] = maxlevel
    datap[10]= LBMCycleIteration # No it du cycle integration temporelle LBM
    datap[11]= dtdim             # shift pour acceder au info OMP

    # level max pour chaque it du cycle LBM
    for it in range(nitCyclLBM):
      for level in range(maxlevel,0,-1):
         it_tg = 2**(level-1) 
         if it%it_tg == 0: 
             datap[12+it] = level
             #print('Nblevel',datap[12+it],'itCycl=',it )
             break
       
    for it in range(nssiter):
      datap[dtdim  + it  ] = -1   #on initialise a -1 le nbre de "zone" a traiter pour forcer l'init dans warmup
    Internal.createUniqueChild(o, '.Solver#dtloc', 'DataArray_t', datap)

    dtloc = Internal.getNodeFromName1(o, '.Solver#dtloc')[1]
    #partage memoire entre dtloc et Noeud timeMotion,....
    first = Internal.getNodeFromName1(t, 'TimeLevelMotion')
    if first is not None: first[1] = dtloc[3:4]
    first = Internal.getNodeFromName1(t, 'TimeLevelTarget')
    if first is not None: first[1] = dtloc[4:5]
    first = Internal.getNodeFromName1(t, 'Iteration')
    if first is not None: first[1] = dtloc[7:8]
    first = Internal.getNodeFromName1(t, 'LBMCycleIteration')
    if first is not None: first[1] = dtloc[10:11]



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
    ndom=0
    for b in bases:
        zones = Internal.getNodesFromType2(b, 'Zone_t')
        nzones=len(zones)
        #print("Bases", b[0], nzones)
        for z in zones:
            #print("ZONES", b[0], z[0])
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
               #kv = dims[3]-3 #[AJ] Hard coded for LBM and 1 ghost point
               if kv == -3: kv =1
               if pad: 
                  target = [iv,jv,kv]
                  if target in sizeIJK:
                       l        = sizeIJK.index( target)
                       shiftvar =  shiftvars[l]

            # zone ownData (generated)
            o = Internal.createUniqueChild(z, '.Solver#ownData', 'UserDefinedData_t')

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
            if ngon:  slope ='o2'
            motion          = "none"
            filtrage        = "off"
            io_th           = 0
            cacheblckI      = 2048
            cacheblckJ      = 2
            cacheblckK      = 7
            dtnature        = "global"
            dtc             = -0.000001
            epsi_newton     = 0.1 
            epsi_linear     = 0.01 
            psiroe          = 0.1
            cfl             = 1.
            rotation        = [ 0.,0.,0., 0.,0.,0.,0.,0.]
            ssdom_IJK       = [240,20,900]
            sfd             = 0
            sfd_chi         = 0.
            sfd_delta       = 1.e15
            sfd_init_iter   = 1
            nudging_ampli   = 1.
            nudging_vector  = numpy.ones( 6, dtype=numpy.float64)
            nudging_vector[0]=0 
            nudging_vector[4]=0 
            nudging_vector[5]=0  #pas de rappel sur continuite, energie et SA par defaut 
            nit_inflow      = 10
            epsi_inflow     = 1.e-5
            DES_debug       = 0
            sa_dist         = 0
            extract_res     = 0
            ibc             = numpy.zeros( 7, dtype=Internal.E_NpyInt)
            source          = 0
            cups            = [1.,1.]
            ratiom          = 10000.
            meshtype        = 1  #structured
            senseurtype     = 1  #version celia laurent du schema senseur
            coef_hyper      = [0.009,0.015] # coeff schema hypersonique M Lugrin

            ##LBM
            #===================================================================
            lbm_neq                = 19               # Number of discrete velocities
            lbm_coll_model         = "BGK"            # Collision Model
            lbm_compressible       = 0                # Flag if compressible LBM
            lbm_taug               = 0.5              # Relaxation time | for BGK this is tau for TRT this is big lambda Par defaut:Viscosite nulle (Euler)
            lbm_hrr_sigma          = 0.985            # Par defaut: Valeur M2P2 pour HRR
            lbm_c0                 = 1./math.sqrt(3.) # Lattice speed of sound
            
            lbm_gamma_precon       = 1.
            lbm_dif_coef           = 0.               # nu_local
            # coll_model   = "BGK_O3" #Par defaut: Modele BGK avec eq ordre 3
            lbm_overset  = 0        #Par defaut, pas d'overset
            lbm_dx       = 1.0      #Par defaut, on fixe dx=1
            flag_nslbm   = 0        #Par defaut, zone sans couplage NS-LBM


            ##IBM
            lbm_ibm                = 21
            lbm_ibm_connector      = 1
            lbm_collision_operator = 1

            # Selective Filter
            ## S-M-C   :: streaming - macro - collision 
            ## M-C-B-S :: macro - collision - BCs - streaming
            lbm_selective_filter       = 0               # -1 M-C-B-S || 0 S-M-C   || 1 - Macro props || 2 - f's || 3 - collision operator
            lbm_selective_filter_size  = 9               # Filter size
            lbm_selective_filter_sigma = 0.              # Filter sigma
            lbm_adaptive_filter_chi    = 1.5             # chi var for adaptive filter
            lbm_chi_spongetypeII       = 0.2             # chi value for sponge type II absorbing layer Xu & Sagaut 2013

            # Non-reflecting boundary condition : Absorbing/Sponge Layer
            lbm_sponge             = 0               #  pas de zone eponge par defaut
            lbm_sponge_size        = 0               # Sponge size - Number of grid points
            lbm_sponge_prep        = 0               # Prepare Sponge
            lbm_cbc_prep           = 0               # Prepare Sponge

            lbm_spng_xmin          = 0
            lbm_spng_xmax          = 0
            lbm_spng_ymin          = 0
            lbm_spng_ymax          = 0
            lbm_spng_zmin          = 0
            lbm_spng_zmax          = 0
            lbm_ibm                = 0
            lbm_initialize         = 0
            lbm_isforce            = 0
            lbm_forcex             = [0.,0.,0.]
            lbm_zlim               = 0.

            
            sol  = Internal.getNodeFromName1(z, 'FlowSolution#Centers')
            dist = Internal.getNodeFromName1(sol, 'TurbulentDistance')
            if dist is not None: sa_dist =1

            KWire_p=0.1 
            DiameterWire_p=0.1 
            CtWire_p=0.1

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
                a = Internal.getNodeFromName2(t, 'temporal_scheme')
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
                a = Internal.getNodeFromName1(d, 'coef_hyper')
                if a is not None: coef_hyper = Internal.getValue(a)
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
                a = Internal.getNodeFromName1(d, 'nudging_ampli')
                if a is not None:  nudging_ampli = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'nudging_vector')
                if a is not None:  nudging_vector = a[1]
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

                a = Internal.getNodeFromName1(d, 'KWire_p')
                if a is not None:KWire_p = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'DiameterWire_p')
                if a is not None: DiameterWire_p = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'CtWire_p')
                if a is not None: CtWire_p = Internal.getValue(a)

                #==========================================================
                ##LBM                 
                #==========================================================
                a = Internal.getNodeFromName1(d, 'LBM_velocity_set')
                if a is not None:
                    lattice_name = Internal.getValue(a)
                    if lattice_name =='D3Q19' and kv>=2:
                        lbm_neq = 19
                        lbm_c0 = 1./numpy.sqrt(3.)
                    elif lattice_name=='D3Q27' and kv>=2:
                        lbm_neq = 27
                        lbm_c0 = 1./numpy.sqrt(3.)
                    else:
                        print(lattice_name+' lattice cannot be used on a 2D mesh.')
                        import sys; sys.exit()
                # if kv == 1: lbm_neq=9  # over-ride to default value of 9 for 2d
                else: lbm_neq = 19                  
                if ngon: lbm_neq=0
                
                a = Internal.getNodeFromName1(d, 'lbm_c0')
                if a is not None: lbm_c0 = Internal.getValue(a)

                a = Internal.getNodeFromName1(d, 'LBM_coll_model')
                if a is not None: lbm_coll_model = Internal.getValue(a)
                if   lbm_coll_model == "BGK":  lbm_collision_model = 1
                elif lbm_coll_model == "BGK+": lbm_collision_model = 2
                elif lbm_coll_model == "RR":   lbm_collision_model = 3
                elif lbm_coll_model == "HRR":  lbm_collision_model = 4
                elif lbm_coll_model == 'TRT':  lbm_collision_model = 2 #Works only for FastLBM (AJ)

                a = Internal.getNodeFromName1(d, 'LBM_compressible')
                if a is not None: lbm_compressible = Internal.getValue(a)

                a = Internal.getNodeFromName1(d, 'LBM_hrr_sigma')
                if a is not None: lbm_hrr_sigma = Internal.getValue(a)

                a = Internal.getNodeFromName1(d, 'LBM_relax_time')
                if a is not None: lbm_taug = Internal.getValue(a)

                a = Internal.getNodeFromName1(d, 'LBM_sponge')
                if a is not None: lbm_sponge = Internal.getValue(a)
                # a = Internal.getNodeFromName1(d, 'lbm_sponge')
                # if a is not None: lbm_sponge = Internal.getValue(a)

                a = Internal.getNodeFromName1(d, 'LBM_overset')
                if a is not None: lbm_overset = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'LBM_dx')
                if a is not None: lbm_dx = Internal.getValue(a)

                a = Internal.getNodeFromName1(d, 'LBM_NS')
                if a is not None: flag_nslbm = Internal.getValue(a)

                klbmoverset = 0
                if lbm_overset == 1: klbmoverset = 1



                a = Internal.getNodeFromName1(d, 'lbm_gamma_precon')
                if a is not None: lbm_gamma_precon = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'lbm_dif_coef')
                if a is not None: lbm_dif_coef = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'lbm_selective_filter')
                if a is not None: lbm_selective_filter = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'lbm_selective_filter_size')
                if a is not None: lbm_selective_filter_size = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'lbm_selective_filter_sigma')
                if a is not None: lbm_selective_filter_sigma = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'lbm_adaptive_filter_chi')
                if a is not None: lbm_adaptive_filter_chi = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'lbm_chi_spongetypeII')
                if a is not None: lbm_chi_spongetypeII = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'lbm_sponge_size')
                if a is not None: lbm_sponge_size = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'lbm_spng_xmin')
                if a is not None: lbm_spng_xmin = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'lbm_spng_xmax')
                if a is not None: lbm_spng_xmax = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'lbm_spng_ymin')
                if a is not None: lbm_spng_ymin = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'lbm_spng_ymax')
                if a is not None: lbm_spng_ymax = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'lbm_spng_zmin')
                if a is not None: lbm_spng_zmin = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'lbm_spng_zmax')
                if a is not None: lbm_spng_zmax = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'lbm_ibm')
                if a is not None: lbm_ibm = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'lbm_ibm_connector')
                if a is not None: lbm_ibm_connector = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'lbm_init_method')
                if a is not None: lbm_initialize = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'lbm_isforce')
                if a is not None: lbm_isforce = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'lbm_force')
                if a is not None: lbm_forcex = Internal.getValue(a)
                a = Internal.getNodeFromName1(d, 'lbm_zlim')
                if a is not None: lbm_zlim = Internal.getValue(a) 
                

                # Determination de level et (levelgd obsolete) pour dt local (NS jeanmasson et LBM)            
                level =1; levelg=0; leveld=0
                #a = Internal.getNodeFromName1(d, 'niveaux_temps')
                #if a is not None: level =Internal.getValue(a)
                
            iflow  = 1
            ides   = 0; idist = 1; ispax = 2; izgris = 0; iprod = 0
            azgris = 0.01; addes = 0.2

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
            if   sgsmodel == 'smsm': iles = 1
            elif sgsmodel ==  'msm': iles = 2
            if iles ==1:
               cacheblckI = max(cacheblckI,4)
               cacheblckJ = max(cacheblckJ,4)
               cacheblckK = max(cacheblckK,4)

            iwallmodel = 1
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
            elif slope == "o5"    : islope = 5
            elif slope == "minmod": islope = 3
            elif slope == "o3sc"  : islope = 6
            elif slope == "o5sc"  : islope = 7

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
            elif scheme == "senseur_hyper": 
                kfludom = 8
                islope  = 2

            else: print('Warning: Fast: scheme %s is invalid.'%scheme)

            lale = 0
            if    motion == "none"       : lale = 0
            elif  motion == "rigid"      : lale = 1
            elif  motion == "rigid_ext"  : lale = 2
            elif  motion == "deformation": lale = 3
            else: print('Warning: Fast: motion %s is invalid.'%motion)
            
            iflagflt = 0
            if filtrage == "on": iflagflt = 1

            dtloc = 0
            if dtnature == "global": dtloc = 0
            elif dtnature == "local": dtloc = 1
            else: print('Warning: Fast: time_step_nature %s is invalid.'%dtnature)

            ImplicitSolverNum = 0
            if implicit_solver == "lussor": ImplicitSolverNum = 0
            elif implicit_solver == "gmres": ImplicitSolverNum = 1
            else: print('Warning: Fast: implicit_solver %s is invalid.'%implicit_solver)

            # creation noeud parametre integer

            
            number_of_defines_param_int = 135                           # Number Param INT
            size_int                   = number_of_defines_param_int +1 # number of defines + 1


            datap      = numpy.empty(size_int, Internal.E_NpyInt)
            datap[:]   = -1000
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
               # Nb cache bloc
               cellScheduler = Internal.getNodeFromName1(z, 'CellScheduler')
               Nb_Cache_Bloc = numpy.shape(cellScheduler[1])[2]
               print("Nb_Cache_Bloc", Nb_Cache_Bloc)
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
            datap[56]  = pow(2, veclevel[i]-1)
            #datap[56]  = levelg
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
            datap[69]  = 0
            datap[70]  = 0
            datap[71]   = nbr_relax
            datap[72]   = nbr_restart
            datap[73]   = nbr_krylov
            datap[74]   = ImplicitSolverNum
            datap[75]   = lu_match          
            datap[76:83]= ibc[0:7]
            datap[83]   = source
            datap[84]   = meshtype
            datap[85]   = senseurtype
            datap[86]   = -1
            datap[87]   = iwallmodel
            datap[88]   = wallmodel_sample

            datap[VSHARE.SA_DIST]        = sa_dist
            
            ## LBM
            datap[VSHARE.NEQ_LBM]        = lbm_neq
            datap[VSHARE.LBM_COL_OP]     = lbm_collision_model
            datap[VSHARE.LBM_HLBM]       = lbm_compressible
            
            #Unused parameters at the moment
            datap[VSHARE.LBM_FILTER]     = lbm_selective_filter
            datap[VSHARE.LBM_FILTER_SZ]  = lbm_selective_filter_size
            datap[VSHARE.LBM_SPONGE]     = lbm_sponge
            datap[VSHARE.LBM_SPONGE_SIZE]= lbm_sponge_size
            datap[VSHARE.LBM_SPONGE_PREP]= lbm_sponge_prep

            datap[VSHARE.flag_streaming]           = 1
            datap[VSHARE.flag_macro]               = 1
            datap[VSHARE.flag_collision_operator]  = 1
            datap[VSHARE.flag_collision]           = 1
            datap[VSHARE.LBM_isforce]              = lbm_isforce
            datap[VSHARE.flag_BConQstar_switch]    = 0

            datap[VSHARE.LBM_CORR_TERM]     = 0
            datap[VSHARE.LBM_OVERSET]       = klbmoverset
            datap[VSHARE.LBM_NS]            = flag_nslbm

            datap[VSHARE.LBM_IBC]          = lbm_ibm
            datap[VSHARE.LBM_IBC_NUM]      = 0 
            datap[VSHARE.LBM_IBC_PREP]     = 0
            datap[VSHARE.LBM_IBC_CONNECTOR]= lbm_ibm_connector

            datap[VSHARE.LBM_spng_xmin] = lbm_spng_xmin
            datap[VSHARE.LBM_spng_xmax] = lbm_spng_xmax
            datap[VSHARE.LBM_spng_ymin] = lbm_spng_ymin
            datap[VSHARE.LBM_spng_ymax] = lbm_spng_ymax
            datap[VSHARE.LBM_spng_zmin] = lbm_spng_zmin
            datap[VSHARE.LBM_spng_zmax] = lbm_spng_zmax

            ##Setting pointers to 0 for non-regression tests
            datap[VSHARE.LBM_NQ_BC]          =0 # [non-regression redundancy]
            datap[VSHARE.LBM_BConQstar]      =0 # [non-regression redundancy]
            datap[VSHARE.PT_LBM_Cs]          =0 # [non-regression redundancy]
            datap[VSHARE.PT_LBM_Ws]          =0 # [non-regression redundancy]
            datap[VSHARE.PT_LBM_Cminus]      =0 # [non-regression redundancy]
            datap[VSHARE.PT_LBM_BC]          =0 # [non-regression redundancy]
            datap[VSHARE.PT_LBM_H2H3]        =0 # [non-regression redundancy]
            datap[VSHARE.PT_LBM_SPEC]        =0 # [non-regression redundancy]
            datap[VSHARE.PT_LBM_FILTER_WGHT] =0 # [non-regression redundancy]
            datap[VSHARE.PT_LBM_FILTER_STNCL]=0 # [non-regression redundancy]
            datap[VSHARE.PT_LBM_IBC_LIST]    =0 # [non-regression redundancy]
            datap[VSHARE.PT_LBM_IBC_DIST]    =0 # [non-regression redundancy]
            datap[VSHARE.PT_LBM_IBC_DIR]     =0 # [non-regression redundancy]

            #correction flux paroi ibm
            datap[VSHARE.IBC_PT_FLUX]   = -1

            # options SA
            datap[VSHARE.SA_LOW_RE] = 0 # active Low Reynolds correction
            datap[VSHARE.SA_ROT_CORR] = 0 # active Rotation correction
            if d is not None:
                a = Internal.getNodeFromName1(d, 'SA_add_LowRe')
                if a is not None: 
                    val = Internal.getValue(a)
                    if val == 1 or val == 'active' or val == True: datap[VSHARE.SA_LOW_RE] = 1
                a = Internal.getNodeFromName1(d, 'SA_add_RotCorr')
                if a is not None: 
                    val = Internal.getValue(a)
                    if val == 1 or val == 'active' or val == True: datap[VSHARE.SA_ROT_CORR] = 1


            datap[VSHARE.NONZ] = ndom
            i += 1
         
            Internal.createUniqueChild(o, 'Parameter_int', 'DataArray_t', datap)
            
            #=====================================================================
            # creation noeud parametre real 
            #=====================================================================
            number_of_defines_param_real = 71                                    # Number Param REAL
            size_real                    = number_of_defines_param_real+1
            datap                        = numpy.zeros(size_real, numpy.float64)
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
 
            if lale == 1: 
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

            datap[55]=  coef_hyper[0] 
            datap[56]=  coef_hyper[1] 

            # Ben's WM
            datap[57] = MafzalMode
            datap[58] = AlphaGradP
            datap[59] = NbptsLinelets

            # Wire Model
            datap[60] = numpy.sqrt(((0.25*KWire_p)**2+1))-0.25*KWire_p # DeltaVWire
            datap[61] = KWire_p
            datap[62] = DiameterWire_p
            datap[63] = CtWire_p

            ## Rotation IBM            
            datap[64]=0
            timemotion = Internal.getNodesFromName(z,'TimeMotion')
            if timemotion:
                rotlocal = Internal.getNodeFromType(timemotion, 'TimeRigidMotion_t')
                datap[64]    = Internal.getValue(Internal.getNodeFromName(rotlocal,'MotionType'))

            ##Nudging
            datap[65]    = nudging_ampli #inutile??
            datap[66:72] = nudging_vector[:]*nudging_ampli

            # LBM related stuff
            datap[VSHARE.LBM_c0]        = lbm_c0
            datap[VSHARE.LBM_taug]      = lbm_taug
            datap[VSHARE.LBM_HRR_sigma] = lbm_hrr_sigma

            datap[VSHARE.LBM_difcoef] = lbm_dif_coef
            
            datap[VSHARE.LBM_filter_sigma]        = lbm_selective_filter_sigma
            datap[VSHARE.LBM_forcex:VSHARE.LBM_forcex+3] = lbm_forcex[0:3]
            datap[VSHARE.LBM_adaptive_filter_chi] = lbm_adaptive_filter_chi
            datap[VSHARE.LBM_chi_spongetypeII]    = lbm_chi_spongetypeII
            datap[VSHARE.LBM_gamma_precon]        = lbm_gamma_precon
            datap[VSHARE.LBM_zlim]                = lbm_zlim            
            datap[VSHARE.LBM_DX]                  = lbm_dx
                                   
            Internal.createUniqueChild(o, 'Parameter_real', 'DataArray_t', datap)

            # More
            Internal.createUniqueChild(o, 'CFL_minmaxmoy', 'DataArray_t', [0.,0.,0.])
            Internal.createUniqueChild(o, 'type_zone'    , 'DataArray_t',  0)
            Internal.createUniqueChild(o, 'model', 'DataArray_t', model)
            Internal.createUniqueChild(o, 'temporal_scheme', 'DataArray_t', temporal_scheme)
            Internal.createUniqueChild(o, 'exp_local', 'DataArray_t', exploc)
            Internal.createUniqueChild(o, 'rk', 'DataArray_t', rk)

            ndom += 1

    return None

#==============================================================================
# create work arrays
#==============================================================================
def createWorkArrays__(zones, dtloc, FIRST_IT):
    ndimt=0; ndimcoe=0; ndimwig=0; ndimgrad=0; ndimplan=0; c = 0
    ndima1=0; ndima3=0; ndima4=0;
    global  MX_OMP_SIZE_INT
    # on check sur 1ere zone si implicite: mixte impli/expli impossible pour le moment
    scheme = "implicit"
    a = Internal.getNodeFromName2(zones[0], 'temporal_scheme')
    if a is not None: scheme = Internal.getValue(a)
    rk=3
    rk_ = Internal.getNodeFromName2(zones[0],'rk')
    if rk_ is not None: rk = Internal.getValue(rk_)
   
    exploc = Internal.getNodeFromName2(zones[0], 'exp_local')
    if exploc is not None: exploc = Internal.getValue(exploc)
    else: exploc = 0

    lssiter_loc = 0
    if scheme == 'implicit_local':  lssiter_loc = 1  # sous-iteration local
    if exploc == 1:                 lssiter_loc = 2  # dtloc G Jeanmasson   

    neq_max = 0
    for z in zones:
        model = "Euler"
        a = Internal.getNodeFromName2(z, 'model')
        model = Internal.getValue(a)
        neq = 5
        if model == 'nsspalart' or model =='NSTurbulent': neq = 6
        param_int = Internal.getNodeFromName2(z, 'Parameter_int')[1]
        if model == 'LBMLaminar':
           neq = param_int[VSHARE.NEQ_LBM]
        if scheme == 'implicit' or scheme == 'implicit_local':   neq_coe = neq     
        else:                                                    neq_coe = 1    # explicit

        neq_max = max(neq, neq_max)
        shiftvar = param_int[VSHARE.SHIFTVAR]
        if param_int is not None:
           shiftvar  = param_int[VSHARE.SHIFTVAR]
        else:
           shiftvar = 100
           print('shiftVar=%d'%shiftvar)
           print('create workarray')
           print('Danger sur optimisation cache shiftVar')      
           
        # dim from getZoneDim
        dims = Internal.getZoneDim(z)
        if dims[0] == 'Structured': 
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

        ndimt   +=     neq*nijk       # surdimensionne ici
        ndimcoe += neq_coe*nijk       # surdimensionne ici

        if param_int[VSHARE.KFLUDOM]==2:   #senseur wiggle
            ndimwig +=   3*nijk
        if param_int[VSHARE.SLOPE]==6 or param_int[VSHARE.SLOPE]==7 : #senseur choc
            ndimwig +=   3*nijk

        if param_int[VSHARE.IFLOW]==4:   #schema LBM
           print("workarray: affiner taille tableau lbm")
           ndima1  +=       9*nijk
           ndima3  +=       6*nijk
           ndima4  +=     neq*nijk

        c += 1
     
    #           3 depth     6 faces  5 tableau
    ndimface = neq*3*ndimplan*6*5*len(zones)
    ndimface = min(ndimface, 2000000000)
    #si pas de temps local inactif (rk3)
    #if rk!=3 or exploc !=2: ndimface=1
    if exploc != 1: ndimface=1
    #else: print('taille tab dtloc=%d'%ndimface)

    mx_thread   = OMP_NUM_THREADS       # surdimensionne : doit etre = a OMP_NUM_THREADS
    verrou      = MX_SSZONE*c*MX_SYNCHRO*mx_thread

    timerOmp = numpy.zeros(  mx_thread*2*dtloc[0] + (mx_thread*2+1)*len(zones), dtype=numpy.float64)

#    wig   = KCore.empty(ndimwig, CACHELINE)
#    coe   = KCore.empty(ndimcoe, CACHELINE)
#    drodm = KCore.empty(ndimt  , CACHELINE)
    wig   = numpy.empty(ndimwig , dtype=numpy.float64)

    coe   = numpy.empty(ndimcoe , dtype=numpy.float64)
    flu   = numpy.empty(ndimFlu , dtype=numpy.float64)
    grad  = numpy.empty(ndimgrad, dtype=numpy.float64)

    if   rk==4 and exploc==0: size_drodm = 4*ndimt
    elif rk==5 and exploc==0: size_drodm = 5*ndimt
    else                    : size_drodm =   ndimt
    #print('taille tab drodm=%d'%ndimt)
    drodm     = numpy.empty(ndimt   , dtype=numpy.float64)        
    tab_dtloc = numpy.empty(ndimface, dtype=numpy.float64)

    tab_a1pr  = numpy.empty(ndima1, dtype=numpy.float64)
    tab_a1fd  = numpy.empty(ndima1, dtype=numpy.float64)
    tab_a1hr  = numpy.empty(ndima1, dtype=numpy.float64)
    tab_a3    = numpy.empty(ndima3, dtype=numpy.float64)
    tab_psi   = numpy.empty(ndima4, dtype=numpy.float64)
    drodm2    = numpy.empty(ndima4, dtype=numpy.float64)


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

    lok      = numpy.zeros(verrou     , dtype=Internal.E_NpyInt )
    iskip_lu = numpy.empty(dtloc[0]*2 , dtype=Internal.E_NpyInt )   # (dtloc[0] = nitmax

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
    hook['hors_eq']        = drodm2
    hook['a1_pr']          = tab_a1pr
    hook['a1_fd']          = tab_a1fd
    hook['a1_hr']          = tab_a1hr
    hook['aneq_3']         = tab_a3
    hook['psi_corr']       = tab_psi

    return hook


#==============================================================================
# compactage tableau + init Numa
#==============================================================================
def _compact(t, containers=[Internal.__FlowSolutionNodes__, Internal.__FlowSolutionCenters__], dtloc=None, fields=None, mode=None, init=True):
    if  mode is not None:
      if mode == -1: thread_numa = 0
      else: thread_numa = 1
    else: thread_numa = 0

    if dtloc is None:
       own   = Internal.getNodeFromName1(t, '.Solver#ownData')  # noeud
       dtloc = Internal.getNodeFromName1(own, '.Solver#dtloc')    # noeud

    #print("compact dtloc=", dtloc)

    zones = Internal.getZones(t)
    for z in zones:
        ars = getFields2Compact__(z, containers, fields)
        sh = None ; size = None
        val = [] # valid fields
        for a in ars:
            a1 = a[1]
            if sh is None: sh = a1.shape; size = a1.size; val.append(a)
            elif a1.shape == sh: val.append(a)
        nfields = len(val)
        if nfields > 0:
            #print('ZONEEEE', z[0])
            param_int = Internal.getNodeFromName2(z, 'Parameter_int')  # noeud
            # Create an equivalent contiguous numpy [flat]
            eq = numpy.empty(nfields*(size+ param_int[1][66]), dtype=numpy.float64) # add a shift  between prim. variables (param_int[1][SHIFTVAR])
            c = 0
            if param_int is None:
                raise ValueError("_compact: Parameter_int is missing for zone %s."%z[0])
            for a in val:
                a1 = a[1]
                # Copy elements
                ptr = a1.reshape((size), order='F') # no copy I hope

                if init: fastc.initNuma( ptr, eq, param_int, dtloc, c, thread_numa)
                # Replace numpys with slice
                a[1] = eq[c*(size)+c*param_int[1][66]:(c+1)*(size)+c*param_int[1][66]]
                a[1] = a[1].reshape(sh, order='F')

                c += 1
    return None

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
  if   bcname == "BCExtrapolate":           tag = 0
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
  elif bcname == "BCReconsLBM":             tag =23#LBM 
  elif bcname == "BCdimNS":                 tag =24#LBM
  elif bcname == "BCadimcoins":             tag =25#LBM
  elif bcname == "BCOversetLBM":            tag =26#LBM
  elif bcname == "BCEquilibrium":           tag =27
  elif bcname == "BCWallModel":             tag =30
  elif bcname == "BCWallExchange":          tag =31
  elif bcname == "BCWallViscousIsothermal": tag =32
  elif bcname == "BCTDIBC":                 tag =90;#LBM - TDIBC for PR ALBATOR - DMPE/STAT R.Roncen (copy of BCEquilibrium)
  elif bcname == "BCTDIBCNSCBC":            tag =91;#LBM - TDIBC for PR ALBATOR - DMPE/STAT R.Roncen (copy of BCFarfield w/ NSCBC)
  
  elif bcname == "LBM_BCPeriodic":           tag =100;#ok
  elif bcname == "LBM_BCSymmetryPlane":      tag =101;#ok
  elif bcname == "LBM_BCSlip":               tag =101;#ok
  elif bcname == "LBM_BCWall":               tag =102;#ok
  elif bcname == "LBM_BCInflow":             tag =103;#ok
  elif bcname == "LBM_BCPressureAntiBB":     tag =105;#ok
  elif bcname == "LBM_BCPML":                tag =108;#NOT FULLY IMPLEMENTED
  elif bcname == "LBM_BCCBC_UBB":            tag =109;#work needs to be done (DO NOT USE)
  elif bcname == "LBM_BCCBC_PABB":           tag =110;#work needs to be done (BUG)
  elif bcname == "LBM_BCZerothExtrapol":     tag =111;#ok
  elif bcname == "LBM_BCLinearExtrapol":     tag =112;#ok
  elif bcname == "LBM_BCNeumannCentralDir":  tag =113;#ok
  elif bcname == "LBM_BCExtrapolJunk":       tag =114;#ok

  elif bcname == "LBM_Sponge":               tag =150;

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

            if case == 1:
               #sauvegarde M1   #M1 <- current    # current <- P1   # P1 <- temp 
               ta = caM1[1];    caM1[1] = ca[1];  ca[1] = caP1[1];  caP1[1] = ta  
            elif case == 2:
               #sauvegarde P1   #P1 <- current    # current <- M1   # M1 <- temp 
               ta = caP1[1];    caP1[1] = ca[1];  ca[1] = caM1[1];  caM1[1] = ta
    return None

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

            if model == 'nsspalart' or model == 'NSTurbulent': 
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

            if model == 'nsspalart' or model == 'NSTurbulent': 
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
def switchPointersLBM__(zones, neq_lbm, dtloc):
 
     for z in zones:
        param_int = Internal.getNodeFromName2(z, 'Parameter_int')[1]  # noeud
        level = param_int[VSHARE.LEVEL]
        maxlevel    = dtloc[ 9]
        it_cycl_lbm = dtloc[10]
        level_tg = dtloc[12 +it_cycl_lbm]

        max_it        = 2**( maxlevel-1)

        level_next_it =  maxlevel;
        if it_cycl_lbm != max_it -1 : level_next_it = dtloc[12 +it_cycl_lbm +1]
    
        #if level==1 or (level >=2 and level <= level_next_it) :
        if level <= level_next_it :
           #print("switch Pt", z[0])
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
        o = Internal.getNodeFromName1(z, '.Solver#ownData')

        # on concatene les donnes BC dans param_int et param_real
        param_int = Internal.getNodeFromName1( o, 'Parameter_int')
        size = numpy.shape(param_int[1])
        c = 1
        for s in size: c=c*s
        size_param_int.append(c)

        #print("BC size facesched", size_int, "old param int=",c)
        datap = numpy.zeros(size_int+c, Internal.E_NpyInt)
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

        o             = Internal.getNodeFromName1(z, '.Solver#ownData')
        param_int     = Internal.getNodeFromName1(o, 'Parameter_int')[1]
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

            if btype == 'BCWallViscousIsothermal':
               Prop = Internal.getNodeFromName(bc,'.Solver#Property')
               if Prop is None:
                  ptrange = Internal.getNodesFromType1(bc, 'IndexRange_t')
                  rg  = ptrange[0][1]
                  sz  = max(1, rg[0,1]-rg[0,0] ) * max(1, rg[1,1]-rg[1,0] ) * max(1, rg[2,1]-rg[2,0] )
                  print ('Error: FastC: Tableau Temperature absent pour la CL BCWallViscousIsothermal sur la zone',z[0])
                  print ('Error: FastC: besoin d un noeud .Solver#Property contenant un numpy Temperature de taille',sz,' dans le noeud BC')
                  import os; os._exit(1)

            if btype == 'BCWallExchange':
                sol   = Internal.getNodeFromName1(z, 'FlowSolution#Centers')
                cellN = Internal.getNodeFromName1(sol, 'cellN')
                if cellN is None:
                  C._initVars(z, 'centers:cellN', 1.)
                cellN = Internal.getNodeFromName1(sol, 'cellN')[1]
                ptrange = Internal.getNodesFromType1(bc, 'IndexRange_t')
                rg  = ptrange[0][1]
                shift =1
                if '2couche' in bc[0]: shift=2
                if '3couche' in bc[0]: shift=3
                Ntg= 101 #nbr echant pour calcul moyenne temporelle
                print ('nb couche loi de paroi',shift,'Nb echant pour moyenne:',Ntg)
                  
                Prop = Internal.getNodeFromName(bc,'.Solver#Property')
                if Prop is None:
                     Internal.createUniqueChild(bc,'.Solver#Property','UserDefinedData_t')
                     Prop = Internal.getNodeFromName(bc,'.Solver#Property')

                wmles = Internal.getNodeFromName(Prop, 'WMLES_parameter')
                if wmles is None:
                     sz = max(1, rg[0,1]-rg[0,0] ) * max(1, rg[1,1]-rg[1,0] ) * max(1, rg[2,1]-rg[2,0] )

                     #vars=['AvgDensity','AvgVelocityX','AvgVelocityY','AvgVelocityZ','AvgTemperature']
                     vars=['AvgPlane-Primitive']
                     c =0
                     for v in vars:
                         #tab  =  numpy.zeros(sz*Ntg, numpy.float64)
                         tab  =  numpy.zeros(sz*Ntg*5, numpy.float64)
                         Internal.createUniqueChild(Prop, v, 'DataArray_t', value=tab)
                         c += 1
                     tab   =  numpy.zeros(4, numpy.float64)
                     tab[0]= Ntg
                     tab[1]= 0
                     tab[2]= 0
                     tab[3]= 0
                     Internal.createUniqueChild( Prop, 'WMLES_parameter','DataArray_t',value= tab)
                     wmles  = Internal.getNodeFromName(Prop,'WMLES_parameter')


                #verif que coef mobile est bien le premier tableau
                datas = Internal.getNodesFromType(Prop, 'DataArray_t')
                pos = 0
                for data in datas:
                      if data[0]=='mobile_coef' and pos != 0:
                        print("zone:",z[0],"bc:",bc[0],"mobile_coef doit etre le premier noeud de .Solver#Property")
                        import sys; sys.exit()
                      pos+=1

                if numpy.shape(cellN)[2]==1:
                  kmin=0;kmax=1
                else:
                  kmin= rg[2][0]-1; kmax = rg[2][1]
                #CL I
                if rg[0][1] - rg[0][0]==0:
                     if rg[0][1]==1: ijk_tg =2;          sens = 1
                     else:           ijk_tg =rg[0][1]-4; sens =-1

                     for k in range(kmin,kmax):
                        for j in range(rg[1][0]-1,rg[1][1]):
                          for s in range(shift):
                            cellN[ijk_tg+ s*sens, j , k]=2.
                #CL J
                if rg[1][1] - rg[1][0]==0:
                     if rg[1][1]==1: ijk_tg =2;          sens = 1
                     else:           ijk_tg = rg[1][1]-4;sens =-1

                     for k in range(kmin,kmax):
                        for i in range(rg[0][0]-1,rg[0][1]):
                          for s in range(shift):
                             cellN[i, ijk_tg+ s*sens  , k]=2.
                #CL K
                if numpy.shape(cellN)[2]!=1:
                    if rg[2][1] - rg[2][0]==0:
                       if rg[2][1]==1: ijk_tg =2;           sens = 1
                       else:           ijk_tg = rg[2][1]-4; sens =-1
                       for j in range(rg[1][0]-1,rg[1][1]):
                         for i in range(rg[0][0]-1,rg[0][1]):
                           for s in range(shift):
                              cellN[i, j, ijk_tg+ s*sens ]=2.
              
            bcdata  = Internal.getNodesFromType3(bc, 'DataArray_t')
            Nb_data = len(bcdata)

            for data in bcdata:
               size = numpy.shape(data[1])
               c = 1
               for s in size: c=c*s
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

        datap = numpy.zeros(size_int+c, Internal.E_NpyInt)
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
            ind_bc  = numpy.zeros(6, Internal.E_NpyInt)
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
# Construit les donnees compactees pour traiter les interpolation temporelles LBM
#==============================================================================
def _InterpTemporelcompact(t,tc):
    
    #Test si existance data interp dans l'arbre initial
    lout = False
    zones = Internal.getZones(t)
    for z in zones:
        o = Internal.getNodeFromName1(z, 'Temporal_Interpolation')
        if o is not None:
          lout =True
          break

    # 1er cas de figure :il s'agit d'une reprise.
    # Dans ce cas, on a pas besoin de tout recreer, il suffit de mettre a jour 
    # les pointeurs d'adresse dans param_int et param_real pour prendre en compte
    # un eventuel changement de taille des tableaux pour le nouveau calcul
    if lout: 
       for z in zones:
           # zone ownData (generated)
           o = Internal.getNodeFromName1(z, 'Temporal_Interpolation')
           if o is not None:
             adresse       = Internal.getNodeFromName1(o, 'adresse')[1]
             temporal_int  = Internal.getNodeFromName1(o, 'data_int')
             size = numpy.shape(temporal_int[1])
             size_temporal_int = 1
             for s in size: size_temporal_int *=s

             temporal_real = Internal.getNodeFromName1(o, 'data_real')
             size = numpy.shape(temporal_real[1])
             size_temporal_real=1
             for s in size: size_temporal_real *=s


             #on concatene les donnes interp dans param_int et param_real
             own = Internal.getNodeFromName1(z, '.Solver#ownData')
             param_int = Internal.getNodeFromName1(own, 'Parameter_int')
             size = numpy.shape(param_int[1])
             size_param_int=1
             for s in size: size_param_int *=s

             param_real = Internal.getNodeFromName1(own, 'Parameter_real')
             size = numpy.shape(param_real[1])
             size_param_real=1
             for s in size: size_param_real *=s

             if size_temporal_real != 0:

                 #traitement real
                 new_size = size_temporal_real+size_param_real
                 datap = numpy.empty(new_size, numpy.float64)
                 datap[0:size_param_real]        =    param_real[1][0:size_param_real]
                 datap[size_param_real:new_size] = temporal_real[1][0:size_temporal_real]
                 param_real[1]    = datap
                 temporal_real[1] = datap[size_param_real:new_size]

                 #traitement int
                 new_size = size_temporal_int +size_param_int
                 datap = numpy.empty(new_size, Internal.E_NpyInt)
                 datap[0:size_param_int]        =    param_int[1][0:size_param_int]
                 datap[VSHARE.PT_INTERP]        =    size_param_int
   
                 #actualisation pointeur des raccord au cas ou les tailles du param_int/real soit different a la reprise
                 nrac = temporal_int[1][0]
                 for i in range(1,nrac+1):
                   temporal_int[1][ i      ] = temporal_int[1][ i      ] - adresse[ 0 ] + size_param_int 
                   temporal_int[1][ i +nrac] = temporal_int[1][ i +nrac] - adresse[ 1 ] + size_param_real
   
                 adresse[0] = size_param_int
                 adresse[1] = size_param_real

                 datap[size_param_int:new_size] = temporal_int[1][0:size_temporal_int]
                 param_int[1]    = datap
                 temporal_int[1] = datap[size_param_int:new_size]
   
           else: #zone sans interp

             own = Internal.getNodeFromName1(z, '.Solver#ownData')
             param_int = Internal.getNodeFromName1(own, 'Parameter_int')
             size = numpy.shape(param_int[1])
             size_param_int=1
             for s in size: size_param_int *=s

             new_size = 1 +size_param_int
             datap = numpy.empty(new_size, Internal.E_NpyInt)
             datap[0:size_param_int]        =    param_int[1][0:size_param_int]
             datap[VSHARE.PT_INTERP]        =    size_param_int
             datap[size_param_int]          =    0  #pas de raccord
             param_int[1]                   = datap

       return lout
    # Fin premier cas de figure

    # 2eme cas de figure :il s'agit d'un nouveau calcul
    # Dans ce cas, il faut creer les tableaux de stockage, les donnees pour 
    # param_int et param_real, et creer les noeuds Temporal_Interpolation
    # pour les zones concernees.

    # 2.1 : Dimensionnememnt tableau
    sizeRac={}
    zones_tc = Internal.getZones(tc)
    for z in zones_tc:

        matchs  = Internal.getNodesFromType1(z,'ZoneSubRegion_t')
        #print "ZNAME",z[0]
        for match in matchs:
           zrname     = Internal.getValue(match)
           zR         = Internal.getNodeFromName(t,zrname)
           zD         = Internal.getNodeFromName(t,z[0])
           param_intR = Internal.getNodeFromName(zR,'Parameter_int')[1]
           param_intD = Internal.getNodeFromName(zD,'Parameter_int')[1]
           nslbm      = Internal.getNodeFromName(match,'NSLBM')
           rac_nslbm  = False
           if nslbm is not None and param_intD[VSHARE.IFLOW]==4: rac_nslbm= True # on retient uniquement les raccod nslbm si zD ==lbm
           #recherche raccord entre zone de pas de temps different ou lbmns
           if (param_intD[VSHARE.LEVEL] > param_intR[VSHARE.LEVEL]) or rac_nslbm:
              print ( 'stockage actif:', match[0], 'zD=', zD[0], 'zR=',  zR[0], 'levelDR=', param_intD[VSHARE.LEVEL],param_intR[VSHARE.LEVEL],'nslbm=',rac_nslbm  )
              ptlist    = Internal.getNodeFromName1(match,'PointList')[1]
              interptyp = Internal.getNodeFromName1(match,'InterpolantsType')[1]
              imin=10000;imax=-5
              jmin=10000;jmax=-5
              kmin=10000;kmax=-5
              dims = Internal.getZoneDim(zD)
              imd = dims[1]-1; imdjmd =imd*(dims[2]-1)
              #print('dims',dims)
              #print('imd',imd,dims[2]-1, imdjmd )
              for l in range(numpy.size(ptlist)):
                 typ = interptyp[l]

                 k   = ptlist[l]//imdjmd;
                 j   = (ptlist[l]-k*imdjmd)//imd;
                 i   = (ptlist[l]-j*imd-k*imdjmd);
                 if i < imin: imin =i
                 if j < jmin: jmin =j
                 if k < kmin: kmin =k
                 molecule = typ -1 
                 if   typ == 22 or type ==0:
                    print("ERROR: interpolation leastsquare ou 2d pas prise en compte pour lbm multi dt")
                    import sys; sys.exit()
                 if i+molecule > imax: imax =i+molecule
                 if j+molecule > jmax: jmax =j+molecule
                 if k+molecule > kmax: kmax =k+molecule
             
                 #print("min",imin,imax, i,j,k,ptlist[l] ,l)

              sizeWindow = (imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)

              if   param_intR[VSHARE.IFLOW]==4 and param_intD[VSHARE.IFLOW]==4: typRac =0 #LBM-LBM
              elif param_intR[VSHARE.IFLOW]==4 and param_intD[VSHARE.IFLOW]<=3: typRac= 1 #LBM-NS
              elif param_intR[VSHARE.IFLOW]<=3 and param_intD[VSHARE.IFLOW]==4: typRac= 1 #LBM-NS
              else : typRac= 2                                                  # NS-NS

              if zD[0] in sizeRac:
                 nrac = sizeRac[ zD[0]][0]; nrac +=1
                 #Swin = sizeRac[ zD[0]][1]; Swin = Swin +  [sizeWindow, imin,imax, jmin,jmax,kmin,kmax]
                 Swin = sizeRac[ zD[0]][1]; Swin.append( [sizeWindow, imin,imax, jmin,jmax,kmin,kmax] )
                 typR = sizeRac[ zD[0]][2]; typR = typR + [typRac]
                 sizeRac[zD[0]]= [nrac , Swin ,  typR ]
              else:
                 sizeRac[zD[0]]= [1 , [[sizeWindow, imin,imax, jmin,jmax,kmin,kmax]] , [typRac]] 
     
    #print (sizeRac)


    # 2.2 : Allocation new param_int/real
    for z in Internal.getZones(t):
        # zone ownData (generated)
        o = Internal.getNodeFromName1(z, '.Solver#ownData')

        #on concatene les donnes BC dans param_int et param_real
        param_int = Internal.getNodeFromName1(o, 'Parameter_int')
        size = numpy.shape(param_int[1])
        size_param_int=1
        for s in size: size_param_int *=s
        
        if z[0] in sizeRac:
           nrac    = sizeRac[z[0]][0]
           windows = sizeRac[z[0]][1]
           typRac  = sizeRac[z[0]][2]
           Internal.createUniqueChild(z, 'Temporal_Interpolation', 'UserDefinedData_t')
        else:
           nrac    = 0
           windows =[]
           typRac  = 0

        size_int  = 1 +  nrac*(6+1+1+1+2) # 6:range, 1:type raccord, 1: niveau temps, 1:neq, 2: pointers dans param-real et param_int du  rac
        size_real = 0
        l= 0
        for win in windows:
          if typRac[l] ==0:
             nvars_additional=5
             if param_int[1][VSHARE.LBM_OVERSET]==1: nvars_additional +=12
             size_real +=3*win[0]*(param_int[1][VSHARE.NEQ_LBM]+nvars_additional) # LBM-LBM
          else:                                                                   #3 niveau temporel stockes
             size_real +=3*win[0]*(param_int[1][VSHARE.NEQ])                      # couplage NS/NS ou NS/LBM
          l+=1


        datap = numpy.zeros(size_int+size_param_int, Internal.E_NpyInt)
        datap[0:size_param_int]   = param_int[1][0:size_param_int]
        param_int[1] = datap

        param_real = Internal.getNodeFromName1(o, 'Parameter_real')
        size = numpy.shape(param_real[1])
        size_param_real=1
        for s in size: size_param_real *=s

        if size_real != 0:
            datap = numpy.zeros(size_real+size_param_real, numpy.float64)
            datap[0:size_param_real]  = param_real[1][0:size_param_real]
            param_real[1] = datap

            interp = Internal.getNodeFromName1(z,'Temporal_Interpolation')
            data_real = numpy.empty(size_real, numpy.float64)
            Internal.createUniqueChild(interp, 'data_real', 'DataArray_t', data_real)
            data_int  = numpy.empty(size_int, Internal.E_NpyInt)
            Internal.createUniqueChild(interp, 'data_int', 'DataArray_t', data_int)
            adresse   = numpy.empty(2, Internal.E_NpyInt)
            adresse[0]= size_param_int
            adresse[1]= size_param_real
            Internal.createUniqueChild(interp, 'adresse', 'DataArray_t', adresse)
            #slice pour partage avec param_real
            data_real     = Internal.getNodeFromName1(interp,'data_real')
            data_real[1]  = datap[size_param_real:size_param_real+size_real]
            data_int      = Internal.getNodeFromName1(interp,'data_int')

        ## init new param_int/real
        param_int  = param_int[1]
        param_real = param_real[1]

        pt_bcs_int   =  size_param_int
        param_int[VSHARE.PT_INTERP]=  pt_bcs_int

        pt_bcs_real=  size_param_real 
        param_int[ pt_bcs_int ] = nrac

        size_int  = 1 + 2*nrac  # shift pour nrac et pointeur BC_int et BC_real
        size_real = 0
        for i in range(1,nrac+1):
            #print('rac=', i)
            param_int[ pt_bcs_int +i      ] = size_int  + pt_bcs_int
            param_int[ pt_bcs_int +i +nrac] = size_real + pt_bcs_real

            pt_bc                           =  param_int[ pt_bcs_int +i ] 

            param_int[pt_bc  ] = typRac[i-1]
            param_int[pt_bc+1] = 0 # position stockage niveau tn+1

            #nbre variable a copier
            if typRac[i-1] ==0:
               neq = param_int[VSHARE.NEQ_LBM] + 5
               if param_int[VSHARE.LBM_OVERSET]==1: neq +=12 # LBM-LBM
            else: neq=param_int[VSHARE.NEQ]                  # couplage NS/NS ou NS/LBM
            param_int[pt_bc+2] = neq 

            #pointRange  au sens fortran :  i=1  1ere cellule calculee
            ific = param_int[VSHARE.NIJK+3]
            kfic = param_int[VSHARE.NIJK+4]
            param_int[pt_bc+3] = windows[i-1][1]-ific+1
            param_int[pt_bc+4] = windows[i-1][2]-ific+1
            param_int[pt_bc+5] = windows[i-1][3]-ific+1
            param_int[pt_bc+6] = windows[i-1][4]-ific+1
            param_int[pt_bc+7] = windows[i-1][5]-kfic+1
            param_int[pt_bc+8] = windows[i-1][6]-kfic+1

            size_real +=3*windows[i-1][0]*neq
            size_int  +=9

        #slice pour partage avec param_int
        if size_real != 0:
           data_int[1]  = param_int[size_param_int:size_param_int+size_int]
            
    return lout

#==============================================================================
# Construit les donnees compactees pour traiter des flux conservatif en IBM
#==============================================================================
def _Fluxcompact(t):
    # Data for each Base
    bases = Internal.getNodesFromType1(t, 'CGNSBase_t')

    # Recherche du nombre de famille necessitant correction de flux
    families=[]
    for b in bases:
        zones = Internal.getZones(b)
        for z in zones:
            tmp  = Internal.getNodeFromName1( z, 'Conservative_Flux')
            if tmp is not None:
              tmp1 = Internal.getNodesFromType1( tmp[2], 'UserDefinedData_t')
              for family in tmp1:
                 if family[0] not in families: families.append(family[0])

    Nfamille = len(families)
    #print("Correction de flux: Nb famille=", Nfamille, families)

    #dimensionnememnt tableau
    size_int = HOOK['MX_OMP_SIZE_INT']
    for b in bases:
        zones = Internal.getZones(b)
        for z in zones:
            
            o    = Internal.getNodeFromName1( z, '.Solver#ownData')
            tmp  = Internal.getNodeFromName1( z, 'Conservative_Flux')
            if tmp is not None:
               tmp1 = Internal.getNodesFromType1( tmp[2], 'UserDefinedData_t')
               #calcul place necessaire dans param_int
               size_int =1 + 6*Nfamille
               for family in tmp1:
                  
                  faces= Internal.getNodesFromType1( family, 'DataArray_t')
                  for f in faces:
                    size_int += numpy.size(f[1])
                    #print("compact flux", z[0],size_int)

               #mise a jour size param_int
               param_int = Internal.getNodeFromName1( o, 'Parameter_int')
               size = numpy.shape(param_int[1])
               c = 1
               for s in size: c=c*s

               #print(z[0],"SIZE int", size_int,"c=", c,"size_flux =",  size_int)
               datap = numpy.zeros(size_int + c , Internal.E_NpyInt)

               datap[0:c] = param_int[1][0:c]
               datap[VSHARE.IBC_PT_FLUX] = c                   # debut tableau flux dans param_int
               param_int[1] = datap

               deb = param_int[1][VSHARE.IBC_PT_FLUX]
               param_int[1][ deb ]= Nfamille
               #on reordone les familles pour a	voir toujour l'ordre de families
               ordre = []
               count = numpy.zeros(Nfamille, Internal.E_NpyInt)
               c = 0
               for family in tmp1:
                  i = families.index(family[0])
                  ordre.append( i )
                  count[ i ] = c
                  c += 1
                  #if len(tmp1) !=1: print("verif family",family[0])

               #on concatene les donnes flux dans param_int
               #for family in tmp1:
               deb1 = param_int[1][VSHARE.IBC_PT_FLUX] + 1 + Nfamille*6
               for l in sorted(ordre):

                  family = tmp1[ count[l] ]
                  #if len(tmp1) !=1: print("verif1 family",family[0],l)

                  i = families.index(family[0])
                  deb = param_int[1][VSHARE.IBC_PT_FLUX] + 1 + i*6
                  faces= Internal.getNodesFromType1( family, 'DataArray_t')
                  for f in faces:
                    size = numpy.size(f[1])
                    param_int[1][ deb ]= size
                    param_int[1][ deb1:deb1+size ]= f[1][:]
                    f[1]                          = param_int[1][ deb1:deb1+size ]

                    deb   += 1   
                    deb1  += size

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
         restart=False):
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
        rank = Cmpi.rank; size = Cmpi.size
        
        if split == 'single':
            # Load connect (tc)
            FILE = fileNameC+extC
            tc = Cmpi.convertFile2SkeletonTree(FILE)
            Cmpi._readZones(tc, FILE, rank=rank)
            graphID = Cmpi.computeGraph(tc, type='ID')
            graphIBCD = Cmpi.computeGraph(tc, type='IBCD')
            procDict = D2.getProcDict(tc)
            procList = D2.getProcList(tc, sort=True)
            graph = {'graphID':graphID, 'graphIBCD':graphIBCD, 'procDict':procDict, 'procList':procList }
            Cmpi._convert2PartialTree(tc, rank=rank)
            # Load data (t)
            FILE = fileName+ext
            if restart and os.access('restart.cgns', os.F_OK): FILE = 'restart.cgns'
            if os.access(FILE, os.F_OK):
                t = Cmpi.convertFile2SkeletonTree(FILE)
                Cmpi._readZones(t, FILE, rank=rank)
                Cmpi._convert2PartialTree(t)
            else: t = None
            # Load stat (ts)
            FILE = fileNameS+extS
            if os.access(FILE, os.F_OK):
                ts = Cmpi.convertFile2SkeletonTree(FILE)
                Cmpi._readZones(ts, FILE, rank=rank)
                Cmpi._convert2PartialTree(ts)
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
                Cmpi._readZones(tc, FILE, rank=rank)
                graph = prepGraphs(tc)
                Cmpi._convert2PartialTree(tc)
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

    return t, tc, ts, graph

#==============================================================================
# save t
# IN: NP: 0 (seq run), >0 (mpi run, distributed), <0 (seq, save multiple)
# IN: split: 'single', 'multiple'
# IN: fileName: name of output file or dir
# IN: compress=1 (cartesian), =2 (all)
# single: write restart.cgns (full)
# multiple: write restart/restart_1.cgns, ... (partial trees)
#==============================================================================
def save(t, fileName='restart', split='single', NP=0, compress=0):
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
    Internal._rmNodesFromName(t2, 'Displacement#0')
    Internal._rmNodesFromName(t2, 'Motion')

    Internal._rmNodeByPath(t2, '.Solver#ownData')
    zones = Internal.getZones(t2)
    for z in zones:
        Internal._rmNodeByPath(z, '.Solver#ownData')
    if compress > 0: 
        import Compressor.PyTree as Compressor
        if compress >= 1: Compressor._compressCartesian(t2)
        if compress >= 2: Compressor._compressAll(t2)

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
            objet = prepGraphs(t2)
            
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
             mpirun=False, exploc=0):
    """Load tree and connectivity tree."""
    import os.path
    baseName = os.path.basename(fileName)
    baseName = os.path.splitext(baseName)[0] # name without extension
    fileName = os.path.splitext(fileName)[0] # full path without extension

    graphN = {'graphID':None, 'graphIBCD':None, 'procDict':None, 'procList':None}
    if mpirun: # mpi run
        rank = Cmpi.rank; size = Cmpi.size
        
        if split == 'single':
            FILE = fileName+'.cgns'
            t = Cmpi.convertFile2SkeletonTree(FILE)
            mp = getMaxProc(t)
            #if mp+1 != size: #### COMMENTE PAR GUILLAUME POUR PASSER LE CREATE
            #    raise ValueError('The number of mpi proc (%d) doesn t match the tree distribution (%d).'%(size,mp+1))

            graphN = prepGraphs(t, exploc=exploc)

            t = Cmpi.readZones(t, FILE, rank=rank)
            t = Cmpi.convert2PartialTree(t, rank=rank)
            
        else: # load 1 fichier par proc

            if graph and exploc == 0:

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
                        t      = Internal.merge(t)
                        graphN = prepGraphs(t)
                    else: print('graph non calculable: manque de fichiers connectivite.')


            if graph and exploc == 1: ## dtloc instationnaire

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
                        t = Internal.merge(t)
                        list_graph = []
                        graphN = prepGraphs(t, exploc=exploc)

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
    
    if graph and not exploc: return t, graphN
    elif graph and exploc:   return t, list_graph
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
        rank = Cmpi.rank; size = Cmpi.size
        
        if split == 'single':
            FILE = fileName+'.cgns'
            t = Cmpi.convertFile2SkeletonTree(FILE)

            mp = getMaxProc(t)

            if mp+1 != size: 
                raise ValueError('The number of mpi proc (%d) doesn t match the tree distribution (%d)'%(mp+1,size)) 
            if graph:
                graphN = prepGraphs(t, exploc=0)
                graphN_= prepGraphs(t, exploc=1)
            
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
                        t      = Internal.merge(t)
                        graphN = prepGraphs(t)
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
    mpirun=False, compress=0):
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
    if compress > 0:
        import Compressor.PyTree as Compressor
        if compress >= 1: Compressor._compressCartesian(t2)
        if compress >= 2: Compressor._compressAll(t2)

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
    if len(directory) > 0 and directory[-1] == '/': directory = directory[:-1]
    if fileName[0] == '/': fileName = fileName[1:]      

    graphN = {'graphID':None, 'graphIBCD':None, 'procDict':None, 'procList':None}
    if mpirun: # mpi run
        rank = Cmpi.rank; size = Cmpi.size
        
        if split == 'single':
            # Load connect (tc)
            FILE = directory+'/'+fileName+'.cgns'
            t = Cmpi.convertFile2SkeletonTree(FILE)
            if graph: graphN = prepGraphs(t)
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
                    graphN = pickle.load(file)
                    file.close()
                # Load all skeleton proc files
                else:
                    ret = 1; no= 0; tmp= []
                    while ret == 1:
                       #FILE = '%s/%s_proc%d.cgns'%(directory, fileName, no)
                       FILE = '%s/%s%d.cgns'%(directory, fileName, no)
                       if not os.access(FILE, os.F_OK): ret = 0
                       if ret == 1: 
                          tmp.append(Cmpi.convertFile2SkeletonTree(FILE))
                          no += 1
                    if no == size and tmp != []:
                        t      = Internal.merge(tmp)
                        Internal._sortByName(t,recursive=False)
                        graphN = prepGraphs(t)
                    else: print('graph non calculable: manque de fichiers connectivite.')

            #FILE = '%s/%s_proc%d.cgns'%(directory, fileName, rank)
            FILE = '%s/%s%d.cgns'%(directory, fileName, rank)
            if os.access(FILE, os.F_OK): t = C.convertFile2PyTree(FILE)
            else: t = None

    else: # sequential run
        if split == 'single':
            FILE = directory+'/'+fileName+'.cgns'
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
# Save t.cgns (data) tc.cgns (connectivity file) ts.cgns (stat)
# si split='single': save t.cgns
# si split='multiple': save t.cgns or t/t_1.cgns
# autorestart possible
# return a partial tree t, a donor partial tree tc, 
# a stat partial tree,
# the communication graph for Chimera and abutting transfers
# the communication graph for IBM transfers
# dir is the directory containing files to be read 
#==============================================================================
def saveTree(t, fileName='restart.cgns', split='single', directory='.', graph=False, mpirun=False):
    """Save tree and connectivity tree."""
    import os.path
    baseName = os.path.basename(fileName)
    baseName = os.path.splitext(baseName)[0] # name without extension
    fileName = os.path.splitext(fileName)[0] # full path without extension

    # Rip/add some useless/usefull data (FastS)
    import Converter.PyTree as C
    t2 = C.rmVars(t, 'centers:Density_P1')
    C._rmVars(t2, 'centers:VelocityX_P1')
    C._rmVars(t2, 'centers:VelocityY_P1')
    C._rmVars(t2, 'centers:VelocityZ_P1')
    C._rmVars(t2, 'centers:Temperature_P1')
    C._rmVars(t2, 'centers:TurbulentSANuTilde_P1')


    bases = Internal.getNodesFromType1(t2    , 'CGNSBase_t')       # noeud
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

#==============================================================================
# Enhance IO (with colours) (cleaned)
#==============================================================================
def _print2screen(mssg2print,txt_colour):
    import sys
    io_fatal_warning_pass = 0
    io_fatal_warning      = 1
    io_pass               = 2
    io_warning            = 3
    io_gen_info           = 4
    io_reset              = "\033[0;0m"

    if txt_colour==io_fatal_warning:
        sys.stdout.write("\033[1;31m")
        print(mssg2print)
        sys.stdout.write(io_reset)
        exit()
    elif txt_colour==io_fatal_warning_pass:
        sys.stdout.write("\033[1;31m")
        print(mssg2print)
        sys.stdout.write(io_reset)
    elif txt_colour==io_gen_info:
        sys.stdout.write("\033[1;34m")
        print(mssg2print)
    elif txt_colour==io_pass:
        sys.stdout.write("\033[1;32m")
        print(mssg2print)
    elif txt_colour==io_warning:
        sys.stdout.write("\033[1;33m")
        print(mssg2print)
        sys.stdout.write(io_reset)
    print("---")

#==============================================================================
# Root mean squared Calculation for FastS
# Input: t is ts from FastS
# IN: mode: xyz, cylx, cylz or None
# IN: cartesian: if True and cyl, return cartesian velocities
# Output: ts w/ Mean, RMS, & Reynolds Stresses
#==============================================================================
def calc_post_stats(t, iskeeporig=False, mode=None, cartesian=True):
    if mode is None or mode == 'xyz':
        vars=['Density','MomentumX','MomentumY','MomentumZ','Pressure','Pressure^2',
              'ViscosityEddy','rou^2','rov^2','row^2','rouv','rouw','rovw']
        dict_var={'VelocityX':'rou^2','VelocityY':'rov^2','VelocityZ':'row^2'}
        list_var=[['VelocityX','VelocityY','rouv',"rho_u'v'"],
                  ['VelocityX','VelocityZ','rouw',"rho_u'w'"],
                  ['VelocityY','VelocityZ','rovw',"rho_v'w'"]]

    elif mode == 'cylx':
        vars=['Density','MomentumX','Momentum_t','Momentum_r','Pressure','Pressure^2',
          'ViscosityEddy','rou^2', 'roU_t^2','roU_r^2','rouU_t','rouU_r','roU_tU_r']
        dict_var={'Velocityt':'roU_t^2','Velocityr':'roU_r^2','VelocityX':'rou^2'}
        list_var=[['Velocityt','Velocityr','roU_tU_r',"rho_U_t'U_r'"],
                  ['Velocityt','VelocityX','rouU_t',"rho_U_t'u'"],
                  ['Velocityr','VelocityX','rouU_r',"rho_U_r'u'"]]

    elif mode == 'cylz':
        vars=['Density','Momentum_t','Momentum_r','MomentumZ','Pressure','Pressure^2',
              'ViscosityEddy','roU_t^2','roU_r^2','row^2','roU_tU_r','rowU_t','rowU_r']
        dict_var={'Velocityt':'roU_t^2','Velocityr':'roU_r^2','VelocityZ':'row^2'}
        list_var=[['Velocityt','Velocityr','roU_tU_r',"rho_U_t'U_r'"],
                  ['Velocityt','VelocityZ','rowU_t',"rho_U_t'w'"],
                  ['Velocityr','VelocityZ','rowU_r',"rho_U_r'w'"]]

    else:
        vars=['Density','MomentumX','MomentumY','MomentumZ','Pressure','Pressure^2',
              'ViscosityEddy','rou^2','rov^2','row^2','rouv','rouw','rovw']
        dict_var={'VelocityX':'rou^2','VelocityY':'rov^2','VelocityZ':'row^2'}
        list_var=[['VelocityX','VelocityY','rouv',"rho_u'v'"],
                  ['VelocityX','VelocityZ','rouw',"rho_u'w'"],
                  ['VelocityY','VelocityZ','rovw',"rho_v'w'"]]

    for z in Internal.getZones(t):
        flowsol  = Internal.getNodeFromName1(z,'FlowSolution#Centers')
        
        #Calc Velocity from Momentum & Density
        dens = Internal.getNodeFromName1(flowsol, 'Density')
        for v in vars[1:4]:
            var = Internal.getNodeFromName1(flowsol,v)
            var[0] = 'Velocity'+var[0][-1]
            var[1] = var[1]/dens[1]

        #Pressure RMS
        pres     = Internal.getNodeFromName1(flowsol,'Pressure')
        pres2    = Internal.getNodeFromName1(flowsol,'Pressure^2')
        if iskeeporig: Internal.addChild(flowsol, Internal.copyNode(pres2), pos=-1)
        pres2[0] = pres[0]+'_RMS'
        pres2[1] = numpy.sqrt(numpy.abs(pres2[1]-pres[1]**2)) # is abs necessary? reponse IM: Oui
        

        #Velocity RMS
        for i in dict_var:
            pres     = Internal.getNodeFromName1(flowsol,i)
            pres2    = Internal.getNodeFromName1(flowsol,dict_var[i])
            if iskeeporig: Internal.addChild(flowsol, Internal.copyNode(pres2), pos=-1)
            pres2[0] = pres[0]+'_RMS'
            pres2[1] = numpy.sqrt( numpy.abs(pres2[1]/dens[1] - pres[1]**2) )  # is abs necessary?

        #Reynolds Stress
        for i in list_var:
            v1       = Internal.getNodeFromName1(flowsol,i[0])
            v2       = Internal.getNodeFromName1(flowsol,i[1])
            pres2    = Internal.getNodeFromName1(flowsol,i[2])
            if iskeeporig: Internal.addChild(flowsol, Internal.copyNode(pres2), pos=-1)
            pres2[0] = i[3]
            pres2[1] = pres2[1] - v1[1]*v2[1]*dens[1]

        # Passage en xyz des vitesses si cylindrique
        if mode == "cylz" and cartesian:
            C._initVars(z,"{centers:Radius} = ( {centers:CoordinateX}**2 +{centers:CoordinateY}**2 )**0.5")
            C._initVars(z,"{centers:co} = {centers:CoordinateX}/{centers:Radius}")
            C._initVars(z,"{centers:si} = {centers:CoordinateY}/{centers:Radius}")
            C._initVars(z, '{centers:VelocityX} = {centers:Velocityr}*{centers:co} - {centers:Velocityt}*{centers:si}')
            C._initVars(z, '{centers:VelocityY} = {centers:Velocityr}*{centers:si} + {centers:Velocityt}*{centers:co}')
            C._initVars(z, "{centers:VelocityX_RMS} = {centers:Velocityr_RMS}*{centers:co}**2 + {centers:Velocityt_RMS}*{centers:si}**2 - 2*{centers:si}*{centers:co}*{centers:rho_U_t'U_r'}/{centers:Density}")
            C._initVars(z, "{centers:VelocityY_RMS} = {centers:Velocityr_RMS}*{centers:si}**2 + {centers:Velocityt_RMS}*{centers:co}**2 + 2*{centers:si}*{centers:co}*{centers:rho_U_t'U_r'}/{centers:Density}")
            # C._initVars(z,"{centers:rho_U_rp^2} = {centers:roU_r^2} -{centers:Density}*{centers:Velocityr}**2 ")
            # C._initVars(z,"{centers:rho_U_tp^2} = {centers:roU_t^2} -{centers:Density}*{centers:Velocityt}**2 ")
            # C._initVars(z,"{centers:rho_u'v'} = {centers:rho_U_rp^2} *%20.16g + {centers:rho_U_t'U_r'}*%20.16g -{centers:rho_U_tp^2}*%20.16g "%(cteta*steta , cteta**2 -steta**2,cteta*steta )) 
            C._initVars(z,"{centers:rho_u'w'} = {centers:rho_U_r'w'}*{centers:co} - {centers:rho_U_t'w'}*{centers:si}")       
            C._initVars(z,"{centers:rho_v'w'} = {centers:rho_U_r'w'}*{centers:si} + {centers:rho_U_t'w'}*{centers:co}")   
            C._rmVars(z, ['centers:Velocityr_RMS', 'centers:Velocityt_RMS' ,'centers:Velocityr', 'centers:Velocityt' ,"centers:rho_U_r'w'","centers:rho_U_t'w'","centers:rho_U_t'U_r'",'centers:co','centers:si','centers:Radius'])
        elif mode == "cylx" and cartesian:
            C._initVars(z,"{centers:Radius} = ( {centers:CoordinateY}**2 +{centers:CoordinateZ}**2 )**0.5")
            C._initVars(z,"{centers:co} = {centers:CoordinateZ}/{centers:Radius}")
            C._initVars(z,"{centers:si} = {centers:CoordinateY}/{centers:Radius}")
            C._initVars(z, '{centers:VelocityZ} = {centers:Velocityr}*{centers:co} + {centers:Velocityt}*{centers:si}')
            C._initVars(z, '{centers:VelocityY} = {centers:Velocityr}*{centers:si} - {centers:Velocityt}*{centers:co}')
            C._initVars(z, "{centers:VelocityZ_RMS} = {centers:Velocityr_RMS}*{centers:co}**2 + {centers:Velocityt_RMS}*{centers:si}**2 + 2*{centers:si}*{centers:co}*{centers:rho_U_t'U_r'}/{centers:Density}")
            C._initVars(z, "{centers:VelocityY_RMS} = {centers:Velocityr_RMS}*{centers:si}**2 + {centers:Velocityt_RMS}*{centers:co}**2 - 2*{centers:si}*{centers:co}*{centers:rho_U_t'U_r'}/{centers:Density}")
            C._initVars(z,"{centers:rho_u'w'} = {centers:rho_U_r'u'}*{centers:co} + {centers:rho_U_t'u'}*{centers:si}")
            C._initVars(z,"{centers:rho_u'v'} = {centers:rho_U_r'u'}*{centers:si} - {centers:rho_U_t'u'}*{centers:co}")
            C._rmVars(z, ['centers:Velocityr_RMS', 'centers:Velocityt_RMS' ,'centers:Velocityr', 'centers:Velocityt' ,"centers:rho_U_r'w'","centers:rho_U_t'w'","centers:rho_U_t'U_r'",'centers:co','centers:si','centers:Radius'])
    return t

#==============================================================================
# display CPU efficiency diagnostic
#==============================================================================
def display_cpu_efficiency(t, mask_cpu=0.08, mask_cell=0.01, diag='compact', FILEOUT='listZonesSlow.dat', FILEOUT1='diagCPU.dat', RECORD=None):

 own   = Internal.getNodeFromName1(t, '.Solver#ownData')  # noeud
 dtloc = Internal.getNodeFromName1(own, '.Solver#dtloc')    # noeud

 node = Internal.getNodeFromName1(t, '.Solver#define')
 node = Internal.getNodeFromName1(node, 'omp_mode')
 ompmode = OMP_MODE
 if  node is not None: ompmode = Internal.getValue(node)

 dtloc = Internal.getValue(dtloc) # tab numpy
 ss_iteration  = int(dtloc[0])

 timer_omp = HOOK["TIMER_OMP"]
 ADR = OMP_NUM_THREADS*2*(ss_iteration)
 echant    =  timer_omp[ ADR ]
 if echant == 0.:
   print('nombre iterations insuffisant pour diagnostic: nitrun * ss_iteration > 15')
   return None

 cellCount = numpy.zeros(2*OMP_NUM_THREADS, dtype=Internal.E_NpyInt)

 zones = Internal.getZones(t)
 tps_percell     =0.
 tps_percell_max =0.
 tps_percell_min =0.
 cout            =0.
 cells_tot       =0
 data            ={}
 f1 = open(FILEOUT1, 'w')
 for z in zones:
    echant    =  timer_omp[ ADR ]
    param_int = Internal.getNodeFromName2(z, 'Parameter_int')[1]
    solver_def= Internal.getNodeFromName2(z, '.Solver#define')
    if ompmode == 1:
       Ptomp       = param_int[69]
       PtrIterOmp  = param_int[Ptomp]
       PtZoneomp   = param_int[PtrIterOmp]
       NbreThreads = param_int[ PtZoneomp  + OMP_NUM_THREADS ]
    else:
       NbreThreads = OMP_NUM_THREADS

    if diag == 'compact':
      tps_zone_percell     =0.
      tps_zone_percell_max =0.
      tps_zone_percell_min =100000.
      ithread_max          = 0
      ijkv                 = param_int[20]*param_int[21]*param_int[22]
      for i in range(OMP_NUM_THREADS):
          if param_int[34]==0:#non ALE
              cellCount[2*i]+=timer_omp[ ADR + 2+i*2 ]
          else:
              cellCount[2*i+1]+=timer_omp[ ADR + 2+i*2 ]

          if ompmode == 1:
            ithread = param_int[ PtZoneomp  +  i ]
            if ithread != -2:
               #print "check", z[0],timer_omp[ ADR + 1+i*2 ]/echant, timer_omp[ ADR + 2+i*2 ], i,"echant=",echant,ijkv
               #tps_zone_percell += timer_omp[ ADR + 1+ithread ]
               tps_zone_percell += timer_omp[ ADR + 1+i*2 ]*timer_omp[ ADR + 2+i*2 ]/float(ijkv)*NbreThreads   #tps * Nb cell
               if timer_omp[ ADR + 1+i*2 ] > tps_zone_percell_max: 
                    tps_zone_percell_max = timer_omp[ ADR + 1+i*2 ]
                    ithread_max          = ithread
               if timer_omp[ ADR + 1+i*2 ] < tps_zone_percell_min: 
                    tps_zone_percell_min = timer_omp[ ADR + 1+i*2 ]
                    ithread_min          = ithread
          else:
            tps_zone_percell += timer_omp[ ADR + 1+i*2 ]*timer_omp[ ADR + 2+i*2 ]/float(ijkv)*NbreThreads
            if timer_omp[ ADR + 1+i*2 ] > tps_zone_percell_max: 
                 tps_zone_percell_max = timer_omp[ ADR + 1+i*2 ]
                 ithread_max          = i
            if timer_omp[ ADR + 1+i*2 ] < tps_zone_percell_min: 
                 tps_zone_percell_min = timer_omp[ ADR + 1+i*2 ]
                 ithread_min          = i
      
      tps_percell     += tps_zone_percell/echant/NbreThreads*ijkv
      tps_percell_max += tps_zone_percell_max/echant*ijkv
      tps_percell_min += tps_zone_percell_min/echant*ijkv
      cout_zone        = tps_zone_percell/echant/NbreThreads*ijkv
      #tps_percell     += tps_zone_percell/echant/NbreThreads/
      #tps_percell_max += tps_zone_percell_max/echant
      #cout_zone        = tps_zone_percell/echant/NbreThreads
      cout            += cout_zone
      data[ z[0] ]   = [ tps_zone_percell/echant/NbreThreads, cout_zone]

      f1.write('cpumoy/cell=  '+str(tps_zone_percell/echant/NbreThreads)+"  cpumaxmin= "+str(tps_zone_percell_max/echant)+str(tps_zone_percell_min/echant) +" th lent/rap= "+str(ithread_max)+str(ithread_min)+"  Nthtread actif="+str(NbreThreads)+" dim zone="+str(ijkv)+" dim tot="+str(cells_tot)+"  "+z[0]+'\n')

      cells_tot   += ijkv

      if RECORD is None: print('cpu/cell zone=', z[0],'moy=',tps_zone_percell/echant/NbreThreads,'maxmin=',tps_zone_percell_max/echant, tps_zone_percell_min/echant,'th maxmin=', ithread_max, ithread_min, 'Nbtread actif=', NbreThreads, ' dim zone=',ijkv,'dim tot=',cells_tot)
      if RECORD is not None:
          tape = tps_zone_percell/echant/NbreThreads

      tps_zone_percell = max(tps_zone_percell, 1.e-11)
      tps_zone_percell_max = max(tps_zone_percell_max, 1.e-11)

      perfo = numpy.empty(2, dtype=numpy.float64)
      perfo[0]= int(echant*NbreThreads/tps_zone_percell)
      perfo[1]= int(echant/tps_zone_percell_max)
      ## CUPS  : [0] MEAN VALUE [1] MIN VALUE
      ## 1/CUPS=mus/iter(time * sub-iter)/cellspercore

      Internal.createUniqueChild(solver_def, 'Cups', 'DataArray_t', value=perfo)

    else:
      for i in range(OMP_NUM_THREADS):
          if ompmode ==1:
            ithread = param_int[ PtZoneomp  +  i ]
            if ithread != -2:
               if RECORD is None: print('zone= ', z[0], 'cpu= ',timer_omp[ ADR + 1+ithread*2 ]/echant,' th=  ', ithread, 'echant= ', echant)
          else:
            ithread = i
            if RECORD is None: print('zone= ', z[0], 'cpu= ',timer_omp[ ADR + 1+ithread*2 ]/echant,' th=  ', ithread, 'echant= ', echant)

    ADR+= OMP_NUM_THREADS*2+1

 for i in range(OMP_NUM_THREADS):
   print('Nbr cellule Fixe/ALE: ',cellCount[2*i],cellCount[2*i+1] , "pour th=", i)

 tps_percell/=cells_tot
 if RECORD is None: print('cpu moyen %cell en microsec: ', tps_percell, tps_percell_max/cells_tot,tps_percell_min/cells_tot )

 f1.write('cpu moyen %cell en microsec: '+str(tps_percell)+'   '+str(tps_percell_max/cells_tot)+str(tps_percell_min/cells_tot)+ '\n')
 f1.close()

 f = open(FILEOUT, 'w')
 sizeZones=[]
 for z in zones:
    param_int = Internal.getNodeFromName2(z, 'Parameter_int')[1]
    cout_relatif   = data[z[0]][0]/tps_percell-1.
    effort_relatif = data[z[0]][1]/cout
    if cout_relatif > mask_cpu and effort_relatif > mask_cell:
       sizeZone=[param_int[20],param_int[21],param_int[22]]
       if sizeZone not in sizeZones:
          sizeZones.append(sizeZone)
          if RECORD is None: print('zone ', z[0],'(',param_int[20],param_int[21],param_int[22],'):  surcout cpu= ', cout_relatif ,' , temps necessaire a cette zone (%)=', effort_relatif)

          f.write(z[0]+','+str(param_int[20])+","+str(param_int[21])+","+str(param_int[22])+","+str(param_int[25])+","+str(param_int[27])+","+str(param_int[29])+","+str(param_int[33])+'\n')
 f.close()

 if RECORD is not None: return tape
 else: return None

#====================================================================
# Usefull functions for Chimera + motion
#====================================================================
# Return dictionary of receiver zones
def getDictOfNobNozOfRcvZones(t, intersectionDict):
    """Return the dict of Nob and Noz of receiver zones."""
    dictOfNobOfRcvZones={}
    dictOfNozOfRcvZones={}
    for nob in range(len(t[2])):
        if Internal.getType(t[2][nob]) == 'CGNSBase_t':
            for noz in range(len(t[2][nob][2])):
                z = t[2][nob][2][noz]
                if Internal.getType(z) == 'Zone_t':
                    zname = Internal.getName(z)                
                    if zname in intersectionDict and intersectionDict[zname] != []:
                        dictOfNobOfRcvZones[zname]=nob
                        dictOfNozOfRcvZones[zname]=noz
    return (dictOfNobOfRcvZones, dictOfNozOfRcvZones)

# Ajout symetrique dans un dictionnaire
def _addPair(idic, z1, z2):
    """Add z1 and z2 symetrically to dictionnary idic."""
    if z1 not in idic: idic[z1] = [z2]
    else:
        if z2 not in idic[z1]: idic[z1].append(z2)
    if z2 not in idic: idic[z2] = [z1]
    else:
        if z1 not in idic[z2]: idic[z2].append(z1)
    return None

# Retourne le nob noz des zones donneuses et remplit le dictOfAdt
# center et axis servent dans le cas d'adt cylindrique
# Filter = 'Base' or '*' or 'Base/*toto' or '*/cart*'
def getDictOfNobNozOfDnrZones(tc, intersectionDict, dictOfADT, 
                              cartFilter='CARTESIAN', cylFilter='CYLINDER*', center=(0,0,0), axis=(0,0,1), depth=2, thetaShift=0.,isIbmAle=False):
    """Fill dictOfAdt."""
    cartFilter = cartFilter.split('/')
    cartBaseFilter = cartFilter[0]
    cartZoneFilter = '*'
    if len(cartFilter) > 1: cartZoneFilter = cartFilter[1]
    cylFilter = cylFilter.split('/')
    cylBaseFilter = cylFilter[0]
    #if cylBaseFilter[-1]=='*': cylBaseFilter=cylBaseFilter[0:-1] # WildCard
    cylZoneFilter = '*'
    if len(cylFilter) > 1: cylZoneFilter = cylFilter[1]

    dnrnames=[]
    for i in intersectionDict.values(): dnrnames += i
    dnrnames = list(set(dnrnames))

    dictOfNobOfDnrZones={}; dictOfNozOfDnrZones={}
    for nob in range(len(tc[2])):
        if Internal.getType(tc[2][nob]) == 'CGNSBase_t':
            baseName       = Internal.getName(tc[2][nob])
            baseNamePref   = baseName[0:len(cylBaseFilter)]
            baseNameSelect = baseName
            if isIbmAle:baseNameSelect = baseNamePref
            for nozc in range(len(tc[2][nob][2])):
                zc = tc[2][nob][2][nozc]
                if Internal.getType(zc) == 'Zone_t':
                    zdnrname = Internal.getName(zc)
                    if zdnrname in dnrnames and zdnrname not in dictOfADT:              
                        if fnmatch.fnmatch(baseNameSelect, cartBaseFilter) and fnmatch.fnmatch(zdnrname, cartZoneFilter): 
                            print('INFO: Creating adt cart for %s.'%zdnrname)
                            adt = None
                        elif fnmatch.fnmatch(baseNameSelect, cylBaseFilter) and fnmatch.fnmatch(zdnrname, cylZoneFilter): 
                            print('INFO: Creating adt cyl for %s.'%zdnrname)
                            adt = C.createHookAdtCyl(zc, center, axis, depth=depth, thetaShift=thetaShift)
                        else:
                            print('INFO: Creating standard adt for %s.'%zdnrname) 
                            adt = C.createHook(zc, 'adt')
                        dictOfADT[zdnrname] = adt
                        dictOfNobOfDnrZones[zdnrname] = nob
                        dictOfNozOfDnrZones[zdnrname] = nozc
    return (dictOfNobOfDnrZones, dictOfNozOfDnrZones)

# In the case of deforming grids, coordinates are modified in t
# this functions update tc from t coordinates
def _pushCenters(t, tc, baseNames):
    """Push centers of t in tc."""
    for n in baseNames:
        b = Internal.getNodeFromName1(t, n)
        bc = Internal.getNodeFromName1(tc, n)
        bpc = C.node2Center(b) # centers
        for z in Internal.getZones(bpc): # push centers in tc
            zp = Internal.getNodeFromName2(bc, z[0])
            cx = Internal.getNodeFromPath(z, 'GridCoordinates/CoordinateX')
            cxp = Internal.getNodeFromPath(zp, 'GridCoordinates/CoordinateX')
            cxp[1][:] = cx[1][:]
            cy = Internal.getNodeFromPath(z, 'GridCoordinates/CoordinateY')
            cyp = Internal.getNodeFromPath(zp, 'GridCoordinates/CoordinateY')
            cyp[1][:] = cy[1][:]
            cz = Internal.getNodeFromPath(z, 'GridCoordinates/CoordinateZ')
            czp = Internal.getNodeFromPath(zp, 'GridCoordinates/CoordinateZ')
            czp[1][:] = cz[1][:]
    return None

# Scalar value ramp between it0 et itf
# retourne 0 pour it0, 1. pour itf
def ramp(it, it0, itf):
    """Return a value in [0,1] depending on it."""
    if it >= itf: return 1.
    return (it-it0)*1./(itf-it0)

# Motion ramp between it0 et itf
# Ramp motion node between 0 and correct speed
def _rampMotion(t, it, it0, itf):
    """Ramp motion node."""
    for z in Internal.getZones(t):
        m = Internal.getNodeFromName1(z, 'Motion')
        if m is not None:
            vx = Internal.getNodeFromName1(m, 'VelocityX')
            vy = Internal.getNodeFromName1(m, 'VelocityY')
            vz = Internal.getNodeFromName1(m, 'VelocityZ')
            coeff = ramp(it, it0, itf)
            vx[1][:] = coeff * vx[1][:]
            vy[1][:] = coeff * vy[1][:]
            vz[1][:] = coeff * vz[1][:]
    return None

# Modify time step in t (parameter_real + solver#define)
def _setTimeStep(t, dt):
    """Set time step."""
    zones = Internal.getZones(t)
    for z in zones:
        n = Internal.getNodeFromName2(z, 'Parameter_real')[1]
        n[0] = dt
        time_step = Internal.getNodeFromName1(z, '.Solver#define') 
        time_step = Internal.getNodeFromName1(time_step, 'time_step') 
        Internal.setValue(time_step, dt)
    return None

# Ramp time step between it0 (dt0) et itf (dtf)
def _rampTimeStep(t, it, it0, itf, dt0, dtf):
    """Ramp time step."""
    dt = dt0 + (dtf-dt0)*ramp(it, it0, itf)
    _setTimeStep(t, dt)
    return None

# Modify cfl in t (parameter_real + solver#define)
def _setCFL(t, cfl):
    """Set cfl."""
    zones = Internal.getZones(t)
    for z in zones:
        n = Internal.getNodeFromName2(z, 'Parameter_real')[1]
        n[15] = cfl
        cfl_value = Internal.getNodeFromName1(z, '.Solver#define') 
        cfl_value = Internal.getNodeFromName1(cfl_value, 'cfl') 
        Internal.setValue(cfl_value, cfl)
    return None

# Ramp CFL between it0 (cfl0) and itf (cflf)
def _rampCFL(t, it, it0, itf, cfl0, cflf):
    """Ramp cfl."""
    cfl = cfl0 + (cflf-cfl0)*ramp(it, it0, itf)
    _setCFL(t, cfl)
    return None

# Push t wall centers to teff coordinates
def _pushWalls(t, teff):
    """Push wall coordinates of t in teff."""
    walls = C.extractBCOfType(t, 'BCWall')
    # suppose ordering is the same for walls and teff
    ezones = Internal.getZones(teff)
    for c, w in enumerate(walls):
        z = ezones[c]
        contw = Internal.getNodeFromName1(w, 'GridCoordinates')
        cwx = Internal.getNodeFromName1(contw, 'CoordinateX')[1]
        cwy = Internal.getNodeFromName1(contw, 'CoordinateY')[1]
        cwz = Internal.getNodeFromName1(contw, 'CoordinateZ')[1]
        
        cont = Internal.getNodeFromName1(z, 'GridCoordinates')
        cx = Internal.getNodeFromName1(cont, 'CoordinateX')[1]
        cy = Internal.getNodeFromName1(cont, 'CoordinateY')[1]
        cz = Internal.getNodeFromName1(cont, 'CoordinateZ')[1]

        cx[:] = cwx[:]
        cy[:] = cwy[:]
        cz[:] = cwz[:]
    return None

#==============================================================================
# Stats post processing for IBM
#==============================================================================
def add2inst(tin,tout,dim_in=3,dim_out=3,direction_copy2Dto3D=3,mode=None):
    
    VARSMACRO_global =['Density','VelocityX','VelocityY','VelocityZ']
    VARSMACRO_save   =['Density','VelocityX','VelocityY','VelocityZ']
    if mode == 'cylx':
        VARSMACRO_save   =['Density','VelocityX','Velocityr','Velocityt']
    elif mode == 'cylz':
        VARSMACRO_save   =['Density','VelocityZ','Velocityr','Velocityt']

    ##Remove M1 in t tree
    VARSMACRO=VARSMACRO_global
    vars=[]
    for v in VARSMACRO:vars.append('centers:'+v+'_M1')
    vars.append('centers:Temperature_M1')
    tout = C.rmVars(tout, vars)

    ##Rename Density, Velocities, & Temperature in t
    for z in Internal.getZones(tout):
        zout = Internal.getNodeFromName(z,'FlowSolution#Centers')
        for v in VARSMACRO: Internal._renameNode(zout, v, v+"inst")
        Internal._renameNode(zout, "Temperature", "Temperatureinst")

    ## Selecting variables to copy from tstat to t
    if mode == 'cylx':
        vars=['VelocityY','VelocityZ']
    elif mode == 'cylz':
        vars=['VelocityX','VelocityY']
    tout = C.rmVars(tout, vars)
        
    VARSMACRO = VARSMACRO_save
    VARSMACRO.append('Pressure')
    
    ##Creating variables in t to which the variables in tstat will be copied to
    for v in VARSMACRO: C._initVars(tout,'{centers:'+v+'}=0')
    if dim_in == dim_out:
        for z in Internal.getZones(tout):
            zin = Internal.getNodesFromName(tin, z[0])
            for v in VARSMACRO: C._cpVars(zin,'centers:'+v,z,'centers:'+v)
    else:       
        for z in Internal.getZones(tout):
            for v in VARSMACRO:
                zout   = Internal.getNodeFromName(z,'FlowSolution#Centers')        
                zin    = Internal.getNodeFromName(Internal.getNodesFromName(tin, z[0]), 'FlowSolution#Centers')

                varout = Internal.getNodeFromName(zout, v)[1]
                varin  = Internal.getNodeFromName(zin , v)[1]
            
                sh  = numpy.shape(varout)
                
                if direction_copy2Dto3D==3:
                    for k in range( sh[direction_copy2Dto3D-1]):
                        varout[:,:,k]  = varin[:,:]
                elif direction_copy2Dto3D==2:
                    for k in range( sh[direction_copy2Dto3D-1]):
                        varout[:,k,:]  = varin[:,:]
                else:
                    for k in range( sh[direction_copy2Dto3D-1]):
                        varout[k,:,:]  = varin[:,:]
    return tout
    
def get_wall_values(t,isRANS=False,wallType='BCWall',mode=None):
    ##Removing useless Variables
    vars   =['centers:Densityinst'   ,'centers:Temperatureinst' ,
             'centers:VelocityXinst' ,'centers:VelocityYinst'   , 'centers:VelocityZinst' ,
             'centers:ViscosityEddy' ]
    if mode == 'cylx':
        vars   =['centers:Densityinst'   ,'centers:Temperatureinst' ,
                 'centers:VelocityXinst' ,'centers:Velocityrinst'   , 'centers:Velocitytinst' ,
                 'centers:ViscosityEddy' ]
    elif mode == 'cylz':
        vars   =['centers:Densityinst'   ,'centers:Temperatureinst' ,
                 'centers:Velocityrinst' ,'centers:Velocitytinst'   , 'centers:VelocityZinst' ,
                 'centers:ViscosityEddy' ]
    t = C.rmVars(t, vars)
    
    ##Getting and calculating ref. values
    [RoInf, RouInf, RovInf, RowInf, RoeInf, PInf, TInf, cvInf, MInf,
     ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
     Mus, Cs, Ts, Pr] = C.getState(t)
        
    betas    = Mus*(Ts+Cs)/(Ts**(3./2.))
    gam1cv   = (Gamma-1.)*cvInf
    RoUInf2I = 1./(RouInf*RouInf+RovInf*RovInf+RowInf*RowInf)

    ##Calculating either t avg or t stationary Temp or Press.
    if not isRANS:
        C._initVars(t,'{centers:Temperature}={centers:Pressure}/(287.053*{centers:Density})')
    else:
        C._initVars(t,'{centers:Pressure}={centers:Temperature}*(287.053*{centers:Density})')
            
    ## RM unnecessary data
    test = Internal.getNodeFromType1(t, 'FlowEquationSet_t')
    Internal._rmNode(t, test)
    
    # Version Antoine
    #t = C.center2Node(t, "FlowSolution#Centers")
    #Internal._rmNodesByName(t, "FlowSolution#Centers")
    #Internal._rmGhostCells(t, t, 2, adaptBCs=1)
    #if mode is None:
    #    C._initVars(t,'{ViscosityMolecular}=%5.12f*sqrt({Temperature})/(1.+%5.12f/{Temperature})'%(betas,Cs))
    #    t = P.computeExtraVariable(t, 'ShearStress')
    #w = C.extractBCOfType(t, wallType)

    # Version Alexis
    Internal._rmNodesFromName(t, 'Motion')
    Internal._rmNodesFromName(t, 'Densityinst')
    Internal._rmNodesFromName(t, 'VelocityXinst')
    Internal._rmNodesFromName(t, 'VelocityYinst')
    Internal._rmNodesFromName(t, 'VelocityZinst')
    Internal._rmNodesFromName(t, 'Temperatureinst')
    Internal._rmNodesFromName(t, 'TurbulentSANuTilde_M1')
    if mode is None:
      PE._extractViscosityMolecular(t)
      P._computeGrad2(t,'centers:VelocityX', ghostCells=True)
      P._computeGrad2(t,'centers:VelocityY', ghostCells=True)
      P._computeGrad2(t,'centers:VelocityZ', ghostCells=True)

    FlowSol = Internal.getNodeFromName(t, 'FlowSolution#Centers')
    t = C.center2Node(t, "FlowSolution#Centers")
    R._switchGridAndGridInit(t)
    Internal._rmNodesByName(t, "GridCoordinates#Init")
    Internal._rmGhostCells(t, t, 2, adaptBCs=1)
    t = C.center2Node(t,['centers:gradxVelocityX','centers:gradxVelocityY','centers:gradxVelocityZ','centers:gradyVelocityX','centers:gradyVelocityY','centers:gradyVelocityZ','centers:gradzVelocityX','centers:gradzVelocityY','centers:gradzVelocityZ'])
    Internal._rmNodesByName(t, "FlowSolution#Centers")
    w = C.extractBCOfType(t, 'BCWall')
    w = C.node2Center(w, 'FlowSolution')
    Internal._rmNodesByName(w, 'FlowSolution')
    if mode is None:       
      PE._extractShearStress(w)
      PE._extractFrictionVector(w)
      PE._extractFrictionMagnitude(w)

    return w

def get_skin_friction(w,RoInf,PInf,RoUInf2I):
    w = P.computeExtraVariable(w, 'SkinFriction')
    w = P.computeExtraVariable(w, 'SkinFrictionTangential')
    C._initVars(w,'{Cf}=2*sqrt({SkinFrictionTangentialX}**2+{SkinFrictionTangentialY}**2+{SkinFrictionTangentialZ}**2)*%g*%g'%(RoInf,RoUInf2I))
    C._initVars(w,'{Cp}=2*%g*({Pressure}-%g)*%g'%(RoInf,PInf,RoUInf2I))
    return w
    

def tcStat_IBC(t,tc,vartTypeIBC=2,bcTypeIB=3):
    import Connector.PyTree as X

    ##Remove flow solution
    tc   = Internal.rmNodesByName(tc, 'FlowSolution')

    ##Removing all info in tc execept for IBMs
    for z in Internal.getZones(tc):
        subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        for s in subRegions:
            sname = s[0][0:2]
            if sname == 'ID':
                Internal._rmNode(tc,s)
                
    ##Calculating either t avg or t stationary Temp or Press.
    [RoInf, RouInf, RovInf, RowInf, RoeInf, PInf, TInf, cvInf, MInf,
     ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
     Mus, Cs, Ts, Pr] = C.getState(t)
    
    C._initVars(t,'{centers:Temperature}={centers:Pressure}/(287.053*{centers:Density})')
    VARSMACRO   =['Density','VelocityX','VelocityY','VelocityZ','Temperature']
    for v in VARSMACRO: C._cpVars(t,'centers:'+v,tc,v)
    
    X._setInterpTransfers(t,tc,
                          cellNVariable="cellN",
                          variables=VARSMACRO,
                          variablesIBC=VARSMACRO,
                          varType=vartTypeIBC,
                          bcType=bcTypeIB,
                          storage=1,
                          Gamma=Gamma,Cv=cvInf,MuS=Mus,Cs=Cs,Ts=Ts,compact=0)
    
    tc = Internal.rmNodesByName(tc, 'FlowSolution')
    return tc


#==============================================================================
# Graph related functions
#==============================================================================
def prepGraphs(t, exploc=0):
    graphID   = Cmpi.computeGraph(t, type='ID'  , reduction=False, exploc=exploc)
    graphIBCD = Cmpi.computeGraph(t, type='IBCD', reduction=False, exploc=exploc)
    procDict  = D2.getProcDict(t)
    procList  = D2.getProcList(t, sort=True)
    list_graph= []
    if exploc == 0:
        graphN = {'graphID':graphID, 'graphIBCD':graphIBCD, 'procDict':procDict, 'procList':procList }
    else:
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
        graphN = list_graph
    return graphN


def printGraph(t, directory='.', exploc=0):
    graphN = prepGraphs(t, exploc)
    # Write graph
    try: import cPickle as pickle
    except: import pickle
    file = open('%s/graph.pk'%directory, 'wb')
    pickle.dump(graphN, file, protocol=pickle.HIGHEST_PROTOCOL)
    file.close()

# Add immersed boundary - method = penalization
# addTripIBCPenalization [CB:WL] ->FastS
def addTripIBCPenalization(t, listTrips):
    import FastS.PyTree as FastS
    count = 1
    for tsurf in listTrips:
        bases = Internal.getBases(tsurf)
        count = 1
        for b in bases:
            Internal.setName(b, "trip"+str(count))
            count += 1
    tsurf  = Internal.merge(listTrips)
    t      = FastS.setIBCData_zero(t, tsurf)
    return t

##Connect Near Match
def connectNearmatchAfterGhost(t, dim, depth=2):
    import Connector.PyTree as X
    t = C.fillEmptyBCWith(t, 'inactive', 'BCExtrapolate', dim=dim)
    t = C.rmBCOfType(t, 'BCNearMatch')
    t = C.fillEmptyBCWith(t, 'nm', 'BCOverlap', dim=dim)
    t = X.applyBCOverlaps(t, depth=depth)
    t = C.rmBCOfType(t, 'BCExtrapolate')
    return t


#==============================================
# Pointwise related functions
#==============================================
def convertPointwise2Fast(FILEIN):
    """This script converts Pointwise meshes to a format that can be used for FastS."""
    t = C.convertFile2PyTree(FILEIN+".cgns")
    ##A few comments: Assumes unspecified is empty!!
    for fam in Internal.getNodesFromType(t, 'Family_t'):
        if fam[0] == "Unspecified": Internal.rmNode(t,fam)

    dicofambc={}                                                  
    for fam in Internal.getNodesFromType(t, 'Family_t'):                                   
        fambc = Internal.getValue(Internal.getNodeFromType(fam, 'FamilyBC_t'))                       
        dicofambc[fam[0]] = fambc                                                             
                                                                                                                                                             
    for bc in Internal.getNodesFromType(t, 'BC_t'):                                              
        if Internal.getValue(bc) == 'FamilySpecified':                                                 
            famname = Internal.getNodeFromType(bc, 'FamilyName_t')                                    
            valbc = dicofambc[Internal.getValue(famname)]                                              
            Internal.setValue(bc, valbc)
            Internal.rmNode(bc, famname)                                                        
    Internal._rmNodesByType(t, 'Family_t')                                                         

    ## Remove further unneccesary information
    ## output is bare bones for FastS
    vars = ['Descriptor_t','FamilyName_t','DataClass_t','DimensionalUnits_t','GridLocation_t','FamilyName','DimensionalExponents_t']
    for v in vars:
        for fam in Internal.getNodesFromType(t, v):
            Internal.rmNode(t, fam)
    
    return t

def _changeBCName4Pointwise(t):
    zones = Internal.getNodesFromType2(t, 'Zone_t')
    count=0
    for z in zones:
        zonebc = Internal.getNodesFromType(z, 'BC_t')
        for zbc in zonebc:
            count += 1
            Internal.setName(zbc, Internal.getValue(zbc)+str(count));
    return None

def _pointwise2D2Fast(t):
    import numpy as np
    zones = Internal.getNodesFromType2(t, 'Zone_t')
    for z in zones:
        zonebc = Internal.getNodesFromType(z, 'BC_t')
        for zbc in zonebc:
            s = np.ones((3,2), dtype=Internal.E_NpyInt)
            for j in range(0,2):
                for i in range(0,2):
                    s[j][i]=zbc[2][0][1][j][i]
            Internal.setValue(zbc[2][0],s)
        zonebc = Internal.getNodesFromType(z, 'GridConnectivity1to1_t')
        for zbc in zonebc:
            s  = np.ones((3,2), dtype=Internal.E_NpyInt)
            s2 = np.ones((3,2), dtype=Internal.E_NpyInt)
            for j in range(0,2):
                for i in range(0,2):
                    s [j][i]=zbc[2][1][1][j][i]
                    s2[j][i]=zbc[2][2][1][j][i]
            Internal.setValue(zbc[2][1],s)
            Internal.setValue(zbc[2][2],s2)

        zonegrid = Internal.getNodesFromType(z, 'GridCoordinates_t')
        for zgrid in zonegrid:            
            if Internal.getNodeFromName(zgrid, 'CoordinateZ') is None:
                coordz=Internal.copyNode(Internal.getNodeFromName(zgrid, 'CoordinateX'))
                Internal.setName(coordz, 'CoordinateZ')
                Internal.setValue(coordz,np.zeros((np.shape(coordz[1])[0],np.shape(coordz[1])[1])))
                Internal.addChild(zgrid, coordz, pos=-1)
    return None

# cassiopee2Pointwise
def cassiopee2Pointwise(fileName):
    t = C.convertFile2PyTree(fileName)
    t = C.rmBCOfType(t, 'BC*')         
    t = C.rmBCOfType(t, 'BCMatch')     
    t = C.rmBCOfType(t, 'BCOverlap')   
    t = C.rmBCOfType(t, 'UserDefined') 
    t = Internal.rmNodesByName(t, 'FlowSolution#Init')
    baseName = os.path.basename(fileName)
    baseName = os.path.splitext(baseName)[0] # name without extension
    fileName = os.path.splitext(fileName)[0] # full path without extension
    C.convertPyTree2File(t, fileName+'_pointwise.cgns')
