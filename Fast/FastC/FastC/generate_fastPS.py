modules=['FastS', 'FastP']
#modules=['FastP']
files = []

#files.append("Compute/init_rhs.for")
#files.append("Compute/invist.for")
#files.append("Compute/cptst3.for")
#files.append("Compute/core3ark3.for")
#files.append("Compute/shape_tab_mtr.for")
files.append("Compute/core3as2.for")
files.append("Compute/core3as2_chim.for")
#files.append("Compute/src_term.for")

print(files)

for i in files:
    print('traitement fichier: %d'%i)
    f = open(i,'r')                         # ouvrir le fichier input
    lines = f.readlines()
    #loop solver
    for solver in modules:

        fout = '../../'+solver+'/'+solver+'/'+i
        print('generation fichier: %s du module %s'%(fout, solver))

        fo = open(fout,"w")                  # ouvrir le fichier de sortie

        for ligne in lines:
            if solver == 'FastP':
                if 'INTEGER_E' not in ligne:
                    ligne=ligne.replace("Fast", "FastP").replace("ind_loop(3)","1").replace("ind_loop(4)","1").replace("ind_loop(5)","1").replace("ind_loop(6)","1")
                    ligne=ligne.replace("vol( param_int(NDIMDX_MTR) )","vol( param_int(NELTS) )")
                    ligne=ligne.replace("NDIMDX_MTR" ,"NFACES")
                    ligne=ligne.replace("NDIMDX_VENT","NFACES")
                    ligne=ligne.replace("NDIMDX"     ,"NELTS")
                    fo.write(ligne)
                else:
                    #fo.write(ligne.replace("Fast", "FastP").replace("ind_loop(3)","1").replace("ind_loop(4)","1").replace("ind_loop(5)","1"))
                    ligne=ligne.replace("Fast", "FastP")
                    ligne=ligne.replace("ind_loop(3)","1")
                    ligne=ligne.replace("ind_loop(4)","1")
                    ligne=ligne.replace("ind_loop(5)","1")
                    fo.write(ligne)

            else:
                fo.write(ligne.replace("Fast", "FastS"))

        fo.close()                               # fermer le fichier output

    f.close()                               # fermer le fichier input
