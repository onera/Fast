#!/usr/bin/env python
import sys
import subprocess
import shlex

#  python generate_correction.py  repertoire_flux
n = len(sys.argv)
'''if n != 2:
    print(Flux name folder is required as argument: python generate_flu.py  fluxFolder')
    sys.exit()

dico= {}
rep = sys.argv[1]
if rep not in dico.keys():
    print('Flux option not described in the dictionnary') 
    sys.exit()
'''
rep1 = '/FastS/BC'
rep ='/stck1/stck7/stck7.3/imary/Cassiopee/Apps/PModules/FastS'+rep1
rep_build='/stck1/stck7/stck7.3/imary/Cassiopee/Apps/PModules/FastS/build/x86_r8/FastS/BC'

models  = [ 'bvbs_extrapolate',
            'bvbs_periodique',
            'bvbs_periodique_azimuthal',
            'bvbs_wall_viscous_transition',
            'bvbs_wall_viscous_adia',
            'bvbs_wall_inviscid',
            'bvbs_inflow_supersonic',
            'bvbs_farfield',
            'bvbs_outflow',
            'bvbs_outpres',
            'bvbs_inflow',
            'bvbs_inflow_fich',
            'bvbs_inflow_newton']


# 
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

    var_in  = 'rop'
    var_out = 'rop'

    subprocess.call(shlex.split('./tapenade.sh '+ input_file+' "'+output_file+'('+var_out+')/('+var_in+')" '+output_file+ ' '+ output_rep ))

    subprocess.call(shlex.split('mv '+ output_rep+'/'+output_file+'_d.f '+output_rep+'/'+output_file+'_d.for' ))
    subprocess.call(shlex.split('rm '+ output_rep+'/'+output_file+'_d.msg' ))

    #nettoyage fichier generer par tapenade
    #
    tapenade_file = rep1+'/'+output_file+'_d.for'
    tapenade_file1= '../..'+tapenade_file
    fi= open(tapenade_file1,'r')
    print('fich=', tapenade_file1)
    lines = fi.readlines()
    fi.close()
    #supression operation sur drodm
    c = 0
    for l in lines:
        if var_out+'(l' in l: 
            test_line = l.split('=')
            if var_out+'(l' in test_line[0] and '=' in l:
                c1=0
                c2 = min ( len(lines)-1, c+c1+1)
                col =min ( len(lines[c2])-1, 5 )
                while '+' == lines[c2][col]: 
                    c1+=1
                    c2 = min ( len(lines)-1, c+c1+1)
                    col =min ( len(lines[c2])-1, 5 )
                c1+=1
                lines = lines[:c] + lines[c+c1:]; c-=c1

        c+=1

    fo = open(tapenade_file1,'w')
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
        input_file  = model+'.for'
        c_index =0
        c = 0
        for l in lines_srcs: 
            if 'FastS/BC/'+input_file in l: c_index = c
            c+=1
        c_index +=1
        #print('c_index flu', c_index, len(lines_srcs))
        lines_srcs_beg = lines_srcs[0:c_index]
        lines_srcs_end = lines_srcs[c_index:]
        b = lines_srcs_beg + ["            '"+tapenade_file[1:]+"',\n"] + lines_srcs_end
        lines_srcs  = b


#write of srcs.py
for l in lines_srcs: srcs.write(l)
srcs.close() 

