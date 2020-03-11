# Utilities for solvers
# load a single file t.cgns : split='single'
#      a file per processor t_proc%d.cgns : split='multi'
# autorestart possible
# return a partial tree t, a donor partial tree tc, 
# the number of processors, 
# the communication graph for Chimera and abutting transfers
# the communication graph for IBM transfers
# dir is the directory containing files to be read 
# FILES must be t.cgns, tc.cgns and restart.cgns
# prefix='restart' : sort restart.cgns
def loadData(split='single',filedir='.', prefix='restart'):
    
    outprefix = filedir+'/'+prefix
    import os
    import Converter.PyTree as C

    procDict=None; graphID=None; graphIBCD=None
    FILEC = filedir+'/tc.cgns'
    try:
        import Converter.Mpi as Cmpi
        NP = Cmpi.size
        if NP == 1: NP = 0
    except: 
        NP = 0
    if NP > 0:
        import Converter.Mpi as Cmpi
        import Distributor2.PyTree as D2
        NP = Cmpi.size
        rank = Cmpi.rank
        tc = Cmpi.convertFile2SkeletonTree(FILEC)
        tc = Cmpi.readZones(tc,FILEC, rank=rank)
        graphID = Cmpi.computeGraph(tc, type='ID')
        graphIBCD = Cmpi.computeGraph(tc, type='IBCD')
        procDict = D2.getProcDict(tc)
        tc = Cmpi.convert2PartialTree(tc)
        if split == 'single': 
            FILEIN = outprefix+'.cgns'
            # autorestart 
            if (os.access(FILEIN, os.F_OK) == True): FILE = FILEIN
            else: FILE = filedir+'/t.cgns'
            t = Cmpi.convertFile2SkeletonTree(FILE)
            t = Cmpi.readZones(t,FILE, rank=rank)
            t = Cmpi.convert2PartialTree(t)
        else: # un fichier loade par proc
            FILEIN = outprefix+'_proc%d.cgns'%rank
            print('trying to read %s'%FILEIN)
            if os.access(FILEIN, os.F_OK): FILE = FILEIN
            else: FILE = filedir+'/t_proc%d.cgns'%rank
            t = C.convertFile2PyTree(FILE)

    else:
        rank = 0
        tc = C.convertFile2PyTree(FILEC)
        FILEIN = outprefix+'.cgns'    
        if os.access(FILEIN, os.F_OK): FILE = FILEIN
        else: FILE = filedir+'/t.cgns'
        t = C.convertFile2PyTree(FILE)

    return t, tc, NP, procDict, graphID, graphIBCD

#==============================================================================
# save files
#==============================================================================
def saveData(t, NP=0, split='single',filedir='.', prefix='restart'):
    import Converter.PyTree as C
    try:
        import Converter.Mpi as Cmpi
        NP = Cmpi.size
        if NP == 1: NP=1
    except: NP = 0

    outprefix = filedir+'/'+prefix
    if NP > 0:
        if split == 'single': 
            Cmpi.convertPyTree2File(t,outprefix+'.cgns')
        else:
            rank = Cmpi.rank
            C.convertPyTree2File(t,outprefix+'_%dproc.cgns'%rank)

    else: C.convertPyTree2File(t,outprefix+'.cgns')

#==============================================================================
# Essai de recuperer les options de la ligne de commande
# Arg1: mode = 'Full', 'Check' - mode de validation
# Arg2: gfx = 0, 1 - CPlot display ou non
# Arg3: display = 0,1 - Residual display
# Arg4: np = 0 (seq), > (nbre de process MPI)
#==============================================================================
def getArgs():
    mode = 'check'; gfx = 0; display = 0; np = 0
    try:
        import sys
        argc = len(sys.argv)
        if (argc > 1): mode = sys.argv[1]
        if (argc > 2): gfx = int(sys.argv[2])
        if (argc > 3): display = int(sys.argv[3])
        if (argc > 4): np = int(sys.argv[4])
    except: pass
    
    return (mode, gfx, display, np)

#==============================================================================
# Essai de recuperer NP de la ligne de commande
# Retourne 0 par defaut
#==============================================================================
def getArgs1():
    np = 0
    try: 
        import sys
        argc = len(sys.argv)
        if (argc > 1): np = int(sys.argv[1])
    except: pass
    return np
                             
