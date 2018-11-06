# - extractConvergenceHistory (pyTree) -
import Converter.PyTree as C
import FastS.PyTree as FastS

# lecture de l'arbre cgns
t = C.convertFile2PyTree('restart.cgns')

# extraction des residus et creation du fichier "residus.dat"
FastS.extractConvergenceHistory(t,"residus.dat")
