"""Structured grid Navier-Stokes solver.
"""
__version__ = '3.2'
__author__ = "Ivan Mary, Stephanie Peron, Christophe Benoit, Thomas Renaud, Guillaume Jeanmasson"

from . import fasts

#==============================================================================
# metric
#==============================================================================
def metric(xyz, numerics):
    """Compute Metric."""
    ret = fasts.metric(xyz, numerics)
    return ret

def cart(Xo, H, N):
    """Create a cartesian mesh defined by a structured array.
    Usage: cart((xo,yo,zo), (hi,hj,hk), (ni,nj,nk))"""
    return fasts.cart(Xo, H, N)

def initVars(ro, varName, val):
    """Initialise variable sur partition omp"""
    fasts.initVars(ro, varName, val)
    return None

# def compute(array, tijk, sol, nitrun):
# #    return fasts.compute(array, tijk, sol, nitrun)
#     sol  = fasts.compute(array, tijk, sol, nitrun)
#     ro = sol[0]; ro_M1 = sol[1]; ro_P1 = sol[2]
#     t     = ro_M1[1]
#     ro_M1[1] = ro[1]
#     ro[1]    = ro_P1[1]
#     ro_P1[1] = t
#     return [ro, ro_M1, ro_P1]
    
