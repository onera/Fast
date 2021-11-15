# - load (pyTree) -
import Fast.PyTree as Fast
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import Connector.PyTree as X
import Converter.Mpi as Cmpi

# Create case
if Cmpi.rank == 0:
    a = G.cart((0,0,0),(1,1,1),(11,11,11))
    a = Cmpi.setProc(a,0)
    b = G.cart((10,0,0),(1,1,1),(11,11,11))
    b = Cmpi.setProc(b,1)
    t = C.newPyTree(['Base',a,b])
    
    t = X.connectMatch(t, dim=3)
    t = Internal.addGhostCells(t,t,2,adaptBCs=0)
    C.convertPyTree2File(t, 't.cgns')
    tc = C.node2Center(t)
    tc = X.setInterpData(t, tc, loc='centers', storage='inverse')
    C.convertPyTree2File(tc, 'tc.cgns')
Cmpi.barrier()

t,tc,ts,graph = Fast.load('t.cgns', 'tc.cgns', split='single')

#print(t)
#print(tc)
#print(ts)
#print('procDict = ', graph['procDict'])
#print('graphID = ', graph['graphID'])
