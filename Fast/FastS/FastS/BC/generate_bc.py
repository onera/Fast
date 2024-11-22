#!/usr/bin/env python
import sys

#  python generate_correction.py  repertoire_flux
n = len(sys.argv)
if n != 2:
    print('Flux name folder is required as argument: python generate_flu.py  fluxFolder')
    sys.exit()

LundVar=[]
for plan in ['_in','_planlund']:
    for v in ['rho','u','v','w','t','nut']:
        LundVar.append(v+plan)
LundVar.append('lund_param')

dico= {}
dico["BCWallViscous"]           = {'name':'bvbs_wall_viscous_adia'       ,'BCname':'BCWallViscous'           ,'dataVars':[]                                                         , 'nextrank':'vide'       , 'Init':'no'}
dico["BCWallViscousIsothermal"] = {'name':'bvbs_wall_viscous_isothermal' ,'BCname':'BCWallViscousIsothermal' ,'dataVars':['tw']                                                     , 'nextrank':'vide'       , 'Init':'no'}
dico["BCWallViscous_transition"]= {'name':'bvbs_wall_viscous_transition' ,'BCname':'BCWallViscous_transition','dataVars':['random']                                                 , 'nextrank':'vide'       , 'Init':'no'}
dico["BCWallInviscid"]          = {'name':'bvbs_wall_inviscid'           ,'BCname':'BCWallInviscid'          ,'dataVars':[]                                                         , 'nextrank':'vide'       , 'Init':'no'}
dico["BCWallModel"]             = {'name':'bvbs_wallmodel'               ,'BCname':'BCWallModel'             ,'dataVars':[]                                                         , 'nextrank':'vide'       , 'Init':'yes'}
dico["BCWallModelRhs"]          = {'name':'bvbs_wallmodel_rhs'           ,'BCname':'BCWallModelRhs'          ,'dataVars':[]                                                         , 'nextrank':'vide'       , 'Init':'no'}
dico["BCInj1"]                  = {'name':'bvbs_inflow_newton'           ,'BCname':'BCInj1'                  ,'dataVars':['d0x','d0y','d0z','pa','ha','nuext']                      , 'nextrank':'BC_nextrank', 'Init':'no'}
dico["BCOutMFR"]                = {'name':'bvbs_outmfr'                  ,'BCname':'BCOutMFR'                ,'dataVars':['qp']                                                     , 'nextrank':'BC_extrap'  , 'Init':'no'}
dico["BCInjMFR"]                = {'name':'bvbs_injmfr'                  ,'BCname':'BCInjMFR'                ,'dataVars':['d0x','d0y','d0z','qp','ha','nuext']                      , 'nextrank':'BC_nextrank', 'Init':'no'}
dico["BCWallExchange"]          = {'name':'bvbs_wallexchange'            ,'BCname':'BCWallExchange'          ,'dataVars':[]                                                         , 'nextrank':'vide'       , 'Init':'no'}
dico["BCOutpres"]               = {'name':'bvbs_outpres'                 ,'BCname':'BCOutpres'               ,'dataVars':['pext']                                                   , 'nextrank':'BC_extrap'  , 'Init':'no'}
dico["BCInflowSupersonic"]      = {'name':'bvbs_inflow_supersonic'       ,'BCname':'BCInflowSupersonic'      ,'dataVars':[]                                                         , 'nextrank':'BC_nextrank', 'Init':'no'}
dico["BCInflowSupersonicFich"]  = {'name':'bvbs_inflow_supersonic_fich'  ,'BCname':'BCInflowSupersonicFich'  ,'dataVars':['rofich','rufich','rvfich','rwfich','refich','ronutfich'] , 'nextrank':'BC_nextrank', 'Init':'no'}
dico["BCInflowLund"]            = {'name':'bvbs_inflow_lund'             ,'BCname':'BCInflowLund'            ,'dataVars':LundVar                                                    , 'nextrank':'BC_nextrank', 'Init':'yes'}
dico["BCInflowFich"]            = {'name':'bvbs_inflow_fich'             ,'BCname':'BCInflowFich'            ,'dataVars':['ux','uy','uz','dens','temp','nutilde']                   , 'nextrank':'BC_nextrank', 'Init':'no'}
dico["BCInflow"]                = {'name':'bvbs_inflow'                  ,'BCname':'BCInflow'                ,'dataVars':[]                                                         , 'nextrank':'BC_nextrank', 'Init':'yes'}
dico["BCOutflow"]               = {'name':'bvbs_outflow'                 ,'BCname':'BCOutflow'               ,'dataVars':[]                                                         , 'nextrank':'BC_extrap'  , 'Init':'no'}
dico["BCFarfield"]              = {'name':'bvbs_farfield'                ,'BCname':'BCFarfield'              ,'dataVars':[]                                                         , 'nextrank':'BC_nextrank', 'Init':'yes'}
dico["BCExtrapolate"]           = {'name':'bvbs_extrapolate'             ,'BCname':'BCExtrapolate'           ,'dataVars':[]                                                         , 'nextrank':'vide'       , 'Init':'yes'}
dico["BCExtrapolateRansLes"]    = {'name':'bvbs_extrapolate_ransles'     ,'BCname':'BCExtrapolateRansLes'    ,'dataVars':[]                                                         , 'nextrank':'vide'       , 'Init':'yes'}
dico["BCPeriodic"]              = {'name':'bvbs_periodique'              ,'BCname':'BCPeriodic'              ,'dataVars':[]                                                         , 'nextrank':'vide'       , 'Init':'yes'}

bc = sys.argv[1]
if bc not in dico:
    print('BC not described in the dictionnary')
    sys.exit()

vars     = dico[ bc ]['dataVars']
name     = dico[ bc ]['name']
bcname   = dico[ bc ]['BCname']
nextrank = dico[ bc ]['nextrank']
Init     = dico[ bc ]['Init']



# ouvrir le fichier input
f     = open('template_bc.for','r')
lines = f.readlines()


for i in range( len(lines) ):
    lines[i]=lines[i].replace("bvbs_template", name)
    if 'insert0' in lines[i]:
        c_insert0=i+1
    if 'insert1' in lines[i]:
        c_insert1=i+1
    if 'insert2' in lines[i]:
        c_insert2=i+1

linsert=False
outline0='     &  ,'
c=0
for var in vars:
    #outline0+=var+' ,'
    outline0='     &                        ,'+var+' \n'
    lines.insert(c_insert0,outline0)
    c_insert0+=1
    c_insert1+=1
    c_insert2+=1
    outline1='      REAL_E '+var+'(size_data) \n'
    lines.insert(c_insert1,outline1)
    c_insert2+=1
    linsert=True
    c+=1

if linsert:
    lines.insert(c_insert0,'     &                        ,size_data, inc_bc \n')
    lines.insert(c_insert1,'      INTEGER_E size_data, indbci, inc_bc(3) \n')
else:
    for i in range( len(lines) ):
        if "      indbci(j_1,k_1) = 1 + (j"  in lines[i]:
            lines[i]=lines[i].replace("      indbci(j_1,k_1) = 1 + (j_1-inc_bc(2)) + (k_1-inc_bc(3))*inc_bc(1)","" )

if Init =='yes':
    include ='#include       "FastS/BC/'+bcname+'_init.for"  \n'
    lines.insert(c_insert2, include)

fout = dico[ bcname ]['name']+'.for'
fo = open(fout,"w")                  # ouvrir le fichier de sortie
print(bcname,' Scheme: file',fout, 'generated')

lines_del=[]
for i in range( len(lines) ):
    if "BCTarget"  in lines[i]:
        new_bc = bcname
        if 'next rank' in lines[i] and nextrank !='vide': new_bc=nextrank
        lines[i]=lines[i].replace("BCTarget", new_bc)
        if new_bc !='vide':
            lines[i]=lines[i].replace(" !next rank", '')
            #lines_del.append(i-2)
            #lines_del.append(i-1)
        else:
            lines[i]=lines[i].replace(" !next rank", '')

c=0
for i in lines_del:
    print(i), len(lines)
    del lines[i-c]
    c+=1

for l in lines: fo.write(l)
fo.close()                               # fermer le fichier output global
