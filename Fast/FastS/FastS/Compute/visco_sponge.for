c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 43 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine visco_sponge(ndom, param_int,ind_loop,xmut,corr_factor)
c***********************************************************************
c_U   USER : PECHIER
c
c     ACT
c_A    initialisation des tableaux de viscosite 
c_A    par la viscosite laminaire
c
c     VAL
c_V    Navier-Stokes mono-espece.
c
c     INP
c_I    ndom   : domaine courant
c_I    idom   : tableau des domaines concernes
c_I    nidom  : nombre de domaines concernes
c_I    rop     : variable conservatives
c
c     OUT
c
c     I/O
c_/    xmut   : viscosites turbulentes 
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom, param_int(0:*), ind_loop(6)

      REAL_E rop( param_int(NDIMDX)*param_int(NEQ) )
      REAL_E xmut(param_int(NDIMDX))
      REAL_E corr_factor(param_int(NDIMDX))

c Var loc
      INTEGER_E l,i,j,k,lij,v5, lt,ltij, lvo, ideb,jdeb,jfin,jp
      REAL_E temp01,cmus1,coesut,t1,r1, coeffj,lj

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      if(ind_loop(1).gt.ind_loop(2)) return 
      if(ind_loop(3).gt.ind_loop(4)) return 
      if(ind_loop(5).gt.ind_loop(6)) return 


#include   "FastC/HPC_LAYER/loop_begin.for"
             xmut(l)= xmut(l) * corr_factor(l)
#include   "FastC/HPC_LAYER/loop_end.for"

      end
