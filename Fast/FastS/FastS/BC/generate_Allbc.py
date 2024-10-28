#!/usr/bin/env python
import os


BC=["BCWallViscous","BCWallViscousIsothermal","BCWallViscous_transition","BCWallInviscid","BCWallModel","BCInj1","BCOutMFR"]
BC.append("BCInjMFR")
BC.append("BCWallExchange")
BC.append("BCOutpres")
BC.append("BCInflowSupersonic")
BC.append("BCInflowSupersonicFich")
BC.append("BCInflowLund")
BC.append("BCInflowFich")
BC.append("BCInflow")
BC.append("BCOutflow")
BC.append("BCFarfield")
BC.append("BCExtrapolate")
BC.append("BCExtrapolateRansLes")
BC.append("BCPeriodic")


for bc in BC:
   cmd='python generate_bc.py'+' '+bc
   print(cmd)
   os.system(cmd)

