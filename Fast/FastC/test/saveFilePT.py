# - saveFile (pyTree) -
import Fast.PyTree as Fast
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Mpi as Cmpi

a = G.cart((0,0,0),(1,1,1),(11,11,11))
a = Cmpi.setProc(a,0)
b = G.cart((10,0,0),(1,1,1),(11,11,11))
b = Cmpi.setProc(b,1)
t = C.newPyTree(['Base',a,b])
Fast.saveFile(t, fileName='t.cgns', split='single', mpirun=False)
