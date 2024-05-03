# - saveData (pyTree) -
import Fast.Utils as FUtils
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Mpi as Cmpi

a = G.cart((0,0,0),(1,1,1),(11,11,11))
a = Cmpi.setProc(a,0)
b = G.cart((10,0,0),(1,1,1),(11,11,11))
b = Cmpi.setProc(b,1)
t = C.newPyTree(['Base',a,b])
FUtils.saveData(t, split='single', filedir='.')
FUtils.saveData(t, split='multi', filedir='.')
