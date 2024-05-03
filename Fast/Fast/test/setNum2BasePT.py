# - setNum2Base (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Fast.PyTree as Fast

a = C.newPyTree(['Base'])

numb = { 'temporal_scheme': 'explicit',
         'ss_iteration': 30,
         'modulo_verif':20 }

Fast._setNum2Base(a, numb)
Internal.printTree(a)
