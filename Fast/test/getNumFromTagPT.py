# - getNumFromTag (pyTree) -
import Generator.PyTree as G
import Fast.PyTree as Fast

a = G.cart((0,0,0), (1,1,1), (10,10,10))
Fast._setNum2Zones(a, {'solver':'FastIJK', 'scheme':'jameson', 'k2':1.})
num = Fast.getNumFromTag(a)
print num
