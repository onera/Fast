#!/usr/bin/env python
import sys

#  python generate_correction.py  repertoire_flux
n = len(sys.argv)
if n != 2:
    print('Flux name folder is required as argument: python generate_flu.py  fluxFolder')
    sys.exit()

dico= {}
dico["SENSOR_INIT"] = { 'name':'flusenseur_init', 'model':['lamin','SA','euler'], 'TypeMotion':['','ale'], 'TypeMesh':['3dfull','3dhomo','3dcart','2d'], 'TypeSlope':['o3']}
dico["SENSOR"]      = { 'name':'flusenseur'     , 'model':['lamin','SA','euler'], 'TypeMotion':['','ale'], 'TypeMesh':['3dfull','3dhomo','3dcart','2d'], 'TypeSlope':['o3']}
dico["AUSM"]        = { 'name':'fluausm'        , 'model':['lamin','SA','euler'], 'TypeMotion':['','ale'], 'TypeMesh':['3dfull','3dhomo','3dcart','2d'], 'TypeSlope':['o3']}
dico["ROE"]         = { 'name':'fluroe'         , 'model':['lamin','SA','euler'], 'TypeMotion':['','ale'], 'TypeMesh':['3dfull','3dhomo','3dcart','2d'], 'TypeSlope':['minmod','o3','o1']}
dico["SENSORHYPER"] = { 'name':'flushyper'      , 'model':['lamin','SA','euler'], 'TypeMotion':['','ale'], 'TypeMesh':['3dfull','3dhomo','3dcart','2d'], 'TypeSlope':['o3']}

rep = sys.argv[1]
if rep not in dico:
    print('Flux option not described in the dictionnary')
    sys.exit()

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
select= open('template_correction_select.for','r')
lines_select = select.readlines()

#open file for sources compilation
srcs= open('../../srcs.py','r')
lines_srcs = srcs.readlines()
srcs.close() 
srcs= open('../../srcs.py','w')
c = 0
for l in lines_srcs: 
    if 'FastS/Compute/src_term.for' in l: c_index = c
    c += 1
c_index += 1
lines_srcs_beg = lines_srcs[0:c_index]
lines_srcs_end = lines_srcs[c_index:]

c = 0
for l in lines_select: 
   if 'ELSE' in l: c_index = c
   c += 1

lines_select_beg = lines_select[0:c_index]
lines_select_end = lines_select[c_index:c_index+8]

fselecto = open(rep+'/correction_'+flux+'_select.for',"w")                  # ouvrir le fichier de sortie

for ale in TypeMotion:

    ale1 = '_'+ale+'_'
    if ale =='': ale1='_'

    for eq in Model:
      for slope in TypeSlope:
         for typezone in TypeMesh:

                        option =  1000*opt_ale[ ale]  +  100*opt_slp[slope] +  10*opt_mod[eq] + opt_mesh[ typezone]

                        # ouvrir le fichier input
                        f     = open('template_correction.for','r')
                        lines = f.readlines()


                        fout = rep+'/'+typezone+'/'+'corr_'+flux+ale1+eq+'_'+slope+'_'+typezone+'.for'
                        fo = open(fout,"w")                  # ouvrir le fichier de sortie
                        print(rep,' Scheme: file',fout, 'generated')

                        for i in range( len(lines) ):
                             lines[i]=lines[i].replace("FLUX_CONV", rep)

                       # correction pour flushyper (wig dim)
                        if flux == 'flushyper':
                            for i in range( 60 ):
                                lines[i]=lines[i].replace("wig( param_int(NDIMDX)     * 3                  )", "wig( param_int(NDIMDX)     * 4                  )")

                        # suppression fluk en 2d et metrique k (pour que mode debug soit OK)
                        if typezone == '2d': 
                                c = 0
                                for l in lines:
                                    if '3D only' in l: lines = lines[:c] + lines[c+1:]; c-=1
                                    c+=1

                                #lines = lines[:189] + lines[198:250]  + lines[270:]
                                for i in range( len(lines) ):
                                      lines[i]=lines[i].replace("tcz = tk(lt)","").replace("sk      = abs (tcz)","")

                        # suppression minmod
                        if slope != 'minmod': 
                                c = 0
                                for l in lines:
                                    if 'avmin(c,r)' in l: lines = lines[:c] + lines[c+1:]; c-=1
                                    c+=1
                                for i in range( len(lines) ):
                                     lines[i]=lines[i].replace("psiroe,avmin", 'psiroe')
                          
                        # suppression Vitesse entrainement si ale=faux
                        if ale == '':
                                c = 0
                                for l in lines:
                                    if 'ALE only' in l: lines = lines[:c] + lines[c+1:]; c-=1
                                    c+=1
                                #lines = lines[:71] + lines[72:92] + lines[93:178] + lines[183:]
                                for i in range( len(lines) ):
                                      lines[i]=lines[i].replace("lven= indven( i, j, k)","")
                        eq2=''
                        # Viscous flux suppression for Euler
                        if eq == 'euler': 
                                c = 0
                                for l in lines:
                                    if 'Rans' in l: lines = lines[:c] + lines[c+1:]; c-=1
                                    c+=1
                                c = 0
                                for l in lines:
                                    if 'fluVis' in l: lines = lines[:c] + lines[c+1:]; c-=1
                                    c+=1
                        elif eq == 'lamin':
                                c = 0
                                for l in lines:
                                    if 'Rans' in l: lines = lines[:c] + lines[c+1:]; c-=1
                                    c+=1
                        # Folder modification
                        elif eq == 'SA':
                                for i in range( len(lines) ):
                                   lines[i]=lines[i].replace("fluViscRans","SA/fluViscRans").replace("assemble","SA/assemble").replace("flu_send","SA/flu_send")
                                eq2=eq+'_'

                        # creation subroutine fortran du flux
                        name_routine = 'corr_'+flux+ale1+eq+'_'+slope+'_'+typezone
                        name_fluEuler= 'fluFaceEuler' +ale1 +slope+ '_'+typezone
                        name_fluRans = 'fluFace'  +eq +ale1 +slope+ '_'+typezone
                        for i in range( len(lines) ):
                            lines[i]=lines[i].replace("!ALE only","").replace("!3D only","")

                            lines[i]=lines[i].replace("corr_lam_template",name_routine).replace("loopI_begin.for",'loopI'+ale1+'begin.for')
                            lines[i]=lines[i].replace("fluFaceEuler_i",typezone+'/'+name_fluEuler+'_i').replace("fluFaceRans_i",typezone+'/'+name_fluRans+'_i')
                            lines[i]=lines[i].replace("fluFaceEuler_j",typezone+'/'+name_fluEuler+'_j').replace("fluFaceRans_j",typezone+'/'+name_fluRans+'_j')
                            lines[i]=lines[i].replace("fluFaceEuler_k",typezone+'/'+name_fluEuler+'_k').replace("fluFaceRans_k",typezone+'/'+name_fluRans+'_k')

                            lines[i]=lines[i].replace("fluViscLaminar_i",'fluvisq_'+typezone+'_i').replace("fluViscRans_i",'fluvisq_'+eq+'_'+typezone+'_i')
                            lines[i]=lines[i].replace("fluViscLaminar_j",'fluvisq_'+typezone+'_j').replace("fluViscRans_j",'fluvisq_'+eq+'_'+typezone+'_j')
                            lines[i]=lines[i].replace("fluViscLaminar_k",'fluvisq_'+typezone+'_k').replace("fluViscRans_k",'fluvisq_'+eq+'_'+typezone+'_k')


                        for l in lines: fo.write(l)
                        fo.close()                               # fermer le fichier output global
                        
                        #modif makefile
                        target = fout
                        include = True
                        c = 0
                        for l in lines_srcs: 
                              if target in l: include = False
                        srcs_out=[]
                        if include == True: srcs_out.append("            'FastS/Compute/"+fout+"',\n")

                        lines_srcs_beg = lines_srcs_beg +   srcs_out 

                        #flux selection function
                        target = 'option.eq.'+str(option)+') THEN'
                        include = True
                        c = 0
                        for l in lines_select: 
                              if target in l: include = False
                              if 'ELSE' in l: c_index = c
                              c+=1

                        select_out=[]
                        if include == True:
                           select_out.append('       ELSEIF (option.eq.'+str(option)+') THEN\n')
                           select_out.append('                                               \n') 
                           select_out.append('         call '+ name_routine+'(ndom,\n') 
                           select_out.append('     &                 ithread, idir,\n') 
                           select_out.append('     &                 param_int, param_real,\n') 
                           select_out.append('     &                 ind_loop,\n') 
                           select_out.append('     &                 rop, drodm , wig,\n') 
                           select_out.append('     &                 venti, ventj, ventk,\n') 
                           #select_out.append('     &                 ti, tj, tk, vol)\n')  pour correction euler
                           select_out.append('     &                 ti, tj, tk, vol, xmut)\n') 
                           select_out.append('                                               \n') 

                        lines_select_beg = lines_select_beg +   select_out 

                        f.close() 

c       =0
c_index =0
for i in range( len(lines_select_beg) ):
    lines_select_beg[i]=lines_select_beg[i].replace("template_correction_select"                ,'corr_'+flux+'_select' )
    if 'ELSEIF ' in lines_select_beg[i] and c_index ==0: 
        c_index = c
    c+=1
lines_select_beg[c_index]=lines_select_beg[c_index].replace("ELSEIF"                ,'IF ' )
 
#write of srcs.py file for makefile  
target = rep+'/correction_'+flux+"_select.for',"
include = True
c = 0
for l in lines_srcs: 
    if target in l: include = False
if include == True: lines_srcs_beg.append("            'FastS/Compute/"+target+"\n")

for l in lines_srcs_beg: srcs.write(l)
for l in lines_srcs_end: srcs.write(l)
srcs.close()          

#write of rep/flux_select.for
for l in lines_select_beg: fselecto.write(l)
for l in lines_select_end: fselecto.write(l)
fselecto.close()                               # fermer le fichier output global

