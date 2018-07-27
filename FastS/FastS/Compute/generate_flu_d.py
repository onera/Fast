#!/usr/bin/env python
import sys
import subprocess
import shlex

#  python generate_flu.py  repertoire_flux
n = len(sys.argv)
if (n != 2):
    print 'Flux name folder is required as argument: python generate_flu.py  fluxFolder' 
    sys.exit()

dico= {}
dico["SENSOR_INIT"] = { 'name':'flusenseur_init', 'model':['lamin','SA','euler'], 'TypeMotion':['','ale'], 'TypeMesh':['3dfull','3dhomo','3dcart','2d'], 'TypeSlope':['o3']}
dico["SENSOR"]      = { 'name':'flusenseur'     , 'model':['lamin','SA','euler'], 'TypeMotion':['','ale'], 'TypeMesh':['3dfull','3dhomo','3dcart','2d'], 'TypeSlope':['o3']}
#dico["AUSM"]        = { 'name':'fluausm'        , 'model':['lamin','SA','euler'], 'TypeMotion':['','ale'], 'TypeMesh':['3dfull','3dhomo','3dcart','2d'], 'TypeSlope':['o3']}
dico["AUSM"]        = { 'name':'fluausm'        , 'model':['lamin','SA','euler'], 'TypeMotion':[''], 'TypeMesh':['2d'], 'TypeSlope':['o3']}
#dico["ROE"]         = { 'name':'fluroe'         , 'model':['lamin','SA','euler'], 'TypeMotion':['','ale'], 'TypeMesh':['3dfull','3dhomo','3dcart','2d'], 'TypeSlope':['minmod','o3']}
#dico["ROE"]         = { 'name':'fluroe'         , 'model':['euler'], 'TypeMotion':[''], 'TypeMesh':['2d'], 'TypeSlope':['minmod']}
dico["ROE"]         = { 'name':'fluroe'         , 'model':['lamin','SA','euler'], 'TypeMotion':[''], 'TypeMesh':['2d'], 'TypeSlope':['o3','minmod']}
rep = sys.argv[1]
if rep not in dico.keys():
    print 'Flux option not described in the dictionnary' 
    sys.exit()

rep_build='/stck1/stck7/stck7.3/imary/Cassiopee/Apps/PModules/FastS/build/x86_r8/FastS/Compute/'+rep

Model      = dico[ rep ]['model']
TypeMotion = dico[ rep ]['TypeMotion']
TypeMesh   = dico[ rep ]['TypeMesh']
TypeSlope  = dico[ rep ]['TypeSlope']
flux       = dico[ rep ]['name']
 
opt_ale = {'ale':1, '':0}
opt_flu = {'flu_ausm':1, 'flu_senseur':2, 'flu_roe':5 }
opt_mod = {'euler':1, 'lamin':2, 'SA':3 }
opt_slp = {'o1':1, 'o3':2, 'minmod':3 }
opt_mesh= {"3dfull":0, "3dhomo":1, "3dcart":2, "2d":3}


#open template file for flux selection
select= open('template_flu_select_d.for','r')
lines_select = select.readlines()

c = 0
for l in lines_select: 
   if 'ELSE' in l: c_index = c
   c+=1

lines_select_beg = lines_select[0:c_index]
lines_select_end = lines_select[c_index:c_index+8]

fselecto = open(rep+'/'+flux+'_select_d.for',"w")                  # ouvrir le fichier de sortie

#open file for sources compilation
srcs= open('../../srcs.py','r')
lines_srcs = srcs.readlines()
srcs.close() 
srcs= open('../../srcs.py','w')

for ale in TypeMotion:

   ale1 = '_'+ale+'_'
   if ale =='': ale1='_'

   for eq in Model:
      for slope in TypeSlope:
         for typezone in TypeMesh:

			option =  1000*opt_ale[ ale]  +  100*opt_slp[slope] +  10*opt_mod[eq] + opt_mesh[ typezone]

                        #determine le  fichier input pour tapenade
                        input_file  = rep_build+'/'+typezone+'/'+flux+ale1+eq+'_'+slope+'_'+typezone+'.f'
                        output_file = flux+ale1+eq+'_'+slope+'_'+typezone
                        output_rep  = rep+'/'+typezone
                        print'input  file', input_file
                        print'output file', output_file
                        print'output file', output_file

                        var_in  = 'rop'
                        var_out = 'drodm'
                        if eq == 'SA': var_in  = 'rop,xmut'

                        subprocess.call(shlex.split('./tapenade.sh '+ input_file+' "'+output_file+'('+var_out+')/('+var_in+')" '+output_file+ ' '+ output_rep ))
                        subprocess.call(shlex.split('mv '+ output_rep+'/'+output_file+'_d.f '+output_rep+'/'+output_file+'_d.for' ))
                        subprocess.call(shlex.split('rm '+ output_rep+'/'+output_file+'_d.msg' ))

                        #nettoyage fichier generer par tapenade
                        #
                        tapenade_file = output_rep+'/'+output_file+'_d.for'
                        fi= open(tapenade_file,'r')
                        lines = fi.readlines()
                        fi.close()
                        #supression operation sur drodm
                        c = 0
                        for l in lines:
                           if 'drodm(l' in l: lines = lines[:c] + lines[c+1:]; c-=1
                           c+=1
                        c = 0
                        for l in lines:
                           if 'drodmd(' in l and '= 0.0' in l: lines = lines[:c] + lines[c+1:]; c-=1
                           c+=1
                        c = 0
                        for l in lines:
                          if 'DO ii1=1,param_int' in l: 
                              c1=0
                              while ' ENDDO' not in lines[c+c1]: c1+=1
                              c1+=1
                              lines = lines[:c] + lines[c+c1:]; c-=c1
                          c+=1

                        fo = open(tapenade_file,'w')
                        for l in lines: fo.write(l)
                        fo.close()                               # fermer le fichier output

                        
                        #flux selection function
                        target = 'option.eq.'+str(option)+') THEN'
                        include = True
                        c = 0
                        for l in lines_select: 
                              if target in l: include = False
                              if 'ELSE' in l: c_index = c
                              c+=1

                        name_routine = flux+ale1+eq+'_'+slope+'_'+typezone+'_d'

                        select_out=[]
                        if include == True:
                           select_out.append('       ELSEIF (option.eq.'+str(option)+') THEN\n')
                           select_out.append('                                               \n') 
                           select_out.append('           call '+ name_routine+'(ndom, ithread,\n') 
                           select_out.append('     &                 param_int, param_real,\n') 
                           select_out.append('     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,\n') 
                           select_out.append('     &                 synchro_send_sock, synchro_send_th,\n') 
                           select_out.append('     &                 synchro_receive_sock, synchro_receive_th,\n') 
                           select_out.append('     &                 ibloc , jbloc , kbloc ,\n') 
                           select_out.append('     &                 icache, jcache, kcache,\n') 
                           select_out.append('     &                 rop, ropd, drodm, drodmd, wig,\n') 
                           select_out.append('     &                 venti, ventj, ventk,\n') 
                           if eq == 'SA':
                               select_out.append('     &                 ti, tj, tk, vol, xmut, xmutd)\n') 
                           else:
                               select_out.append('     &                 ti, tj, tk, vol, xmut)\n') 
                           select_out.append('                                               \n') 

                        lines_select_beg = lines_select_beg +   select_out 

                        #modif makefile
                        target = tapenade_file
                        include = True
                        for l in lines_srcs: 
                              if target in l: include = False
                        srcs_out=[]
                        if include == True: 
                           #recherche la fonction originelle avant passage tapenade
                           input_file  = rep+'/'+typezone+'/'+flux+ale1+eq+'_'+slope+'_'+typezone+'.for'
                           c = 0
                           for l in lines_srcs: 
                             if 'FastS/Compute/'+input_file in l: c_index = c
                             c+=1
                           c_index +=1
                           lines_srcs_beg = lines_srcs[0:c_index]
                           lines_srcs_end = lines_srcs[c_index:]
                           b = lines_srcs_beg + ["            'FastS/Compute/"+tapenade_file+"',\n"] + lines_srcs_end
                           lines_srcs  = b

c       =0
c_index =0
for i in range( len(lines_select_beg) ):
    lines_select_beg[i]=lines_select_beg[i].replace("template_flu_select_d"                ,flux+'_select_d' )
    if 'ELSEIF ' in lines_select_beg[i] and c_index ==0: 
        c_index = c
    c+=1
lines_select_beg[c_index]=lines_select_beg[c_index].replace("ELSEIF"                ,'IF ' )
 
#write of srcs.py file for makefile  
target = rep+'/'+flux+"_select_d.for',"
include = True
c = 0
for l in lines_srcs: 
    if target in l: include = False
if include == True: 
  lines_srcs_beg.append("            'FastS/Compute/"+target+"\n")
  input_file  = rep+'/'+flux+"_select.for',"
  c = 0
  for l in lines_srcs: 
    if 'FastS/Compute/'+input_file in l: c_index = c
    c+=1
  c_index +=1
  lines_srcs_beg = lines_srcs[0:c_index]
  lines_srcs_end = lines_srcs[c_index:]
  b = lines_srcs_beg + ["            'FastS/Compute/"+target+"\n"] + lines_srcs_end
  lines_srcs  = b

#write of srcs.py
for l in lines_srcs: srcs.write(l)
srcs.close()   


#write of rep/flux_select.for
for l in lines_select_beg: fselecto.write(l)
for l in lines_select_end: fselecto.write(l)
fselecto.close()                               # fermer le fichier output global

