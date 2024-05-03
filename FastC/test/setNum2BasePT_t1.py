# - setNum2Base (pyTree) -
import Converter.PyTree as C
import FastC.PyTree as FastC
import KCore.test as test

a = C.newPyTree(['Base'])

numb = {'temporal_scheme': 'explicit',
        'ss_iteration': 30,
        'modulo_verif':20}

FastC._setNum2Base(a, numb)
test.testT(a, 1)
