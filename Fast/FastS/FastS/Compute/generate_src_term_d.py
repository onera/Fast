#!/usr/bin/env python
import sys
import subprocess
import shlex

#  python generate_correction.py  repertoire_flux
n = len(sys.argv)
'''if n != 2:
    print('Flux name folder is required as argument: python generate_flu.py  fluxFolder')
    sys.exit()

dico= {}
rep = sys.argv[1]
if rep not in dico():
    print('Flux option not described in the dictionnary')
    sys.exit()
'''
rep ='SA'
rep_build='/stck1/stck7/stck7.3/imary/Cassiopee/Apps/PModules/FastS/build/x86_r8/FastS/Compute/'+rep

models      = ['spsource_SA', 'spsource_SA_comp', 'spsource_SA_diff','vispalart']

#open file for sources compilation
srcs= open('../../srcs.py','r')
lines_srcs = srcs.readlines()
srcs.close()
srcs= open('../../srcs.py','w')

for model  in models:

    input_file  = rep_build+'/'+model+'.f'
    output_file = model
    output_rep  = rep
    print('input  file', input_file)
    print('output file', output_file)
    print('output file', output_file)

    if model == 'vispalart':
        var_in  = 'rop'
        var_out = 'xmut'
    else:
        var_in  = 'rop'
        var_out = 'drodm'

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
        if var_out+'(l' in l: lines = lines[:c] + lines[c+1:]; c-=1
        c+=1
    c = 0
    for l in lines:
        if var_out+'d(' in l and '= 0.0' in l: lines = lines[:c] + lines[c+1:]; c-=1
        c+=1
    #on supprime 4 lignes,; faure gaffe si changememnt de version tapenade
    c = 0
    for l in lines:
        if 'DO ii2=1,param_int' in l:
            c1=0
            while ' ENDDO' not in lines[c+c1]: c1+=1
            c1+=1
            lines = lines[:c] + lines[c+c1:]; c-=c1
        c+=1
    c = 0
    for l in lines:
        if 'DO ii1=1,param_int' in l:
            c1=0
            while ' ENDDO' not in lines[c+c1]: c1+=1
            c1+=1
            lines = lines[:c] + lines[c+c1:]; c-=c1
        c+=1

    c = 0
    for l in lines:
        if 'coe(l, icoe_pos) = coe(l, 5)' in l: lines = lines[:c] + lines[c+1:]; c-=1
        c+=1

    fo = open(tapenade_file,'w')
    for l in lines: fo.write(l)
    fo.close()

    #modif makefile
    target = tapenade_file
    include = True
    for l in lines_srcs:
        if target in l: include = False
    srcs_out=[]
    if include == True:
        #recherche la fonction originelle avant passage tapenade
        input_file  = 'SA/'+model+'.for'
        c_index =0
        c = 0
        for l in lines_srcs:
            if 'FastS/Compute/'+input_file in l: c_index = c
            c+=1
        c_index +=1
        #print('c_index flu', c_index, len(lines_srcs))
        lines_srcs_beg = lines_srcs[0:c_index]
        lines_srcs_end = lines_srcs[c_index:]
        b = lines_srcs_beg + ["            'FastS/Compute/"+tapenade_file+"',\n"] + lines_srcs_end
        lines_srcs  = b


#write of srcs.py
for l in lines_srcs: srcs.write(l)
srcs.close()
