      subroutine corr_conservatif(param_int, param_real, nthreadmax,
     &                          ind_loop,  rop,  vol, celln, debit)
c***********************************************************************
c_U   USER : DANDOIS
c
c     ACT
c_A    Selective Frequency Damping
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ind_loop(6), param_int(0:*),nthreadmax
      REAL_E  debit(*), param_real(0:*)
 
      REAL_E rop( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E celln( param_int(NDIMDX) )
      REAL_E vol( param_int(NDIMDX_MTR) )

C Var loc
      INTEGER_E i,j,k,l,ltij,lij,lt,lvo 
      REAL_E amor,coef,roe

      INTEGER_E ndom,v1,v2,v3,v4,v5,v1deb,v2deb,v3deb,v4deb,v5deb,v6deb

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"


      v1 = 0
      v2 = param_int(NDIMDX)
      v3 = param_int(NDIMDX)*2
      v4 = param_int(NDIMDX)*3
      v5 = param_int(NDIMDX)*4

      v1deb =1
      v2deb =1+ nthreadmax
      v3deb =1+ nthreadmax*2
      v4deb =1+ nthreadmax*3
      v5deb =1+ nthreadmax*4
      v6deb =1+ nthreadmax*5

      if(param_int(ITYPZONE).eq.2) then

        lt  = 1
        do  k = ind_loop(5), ind_loop(6)
        do  j = ind_loop(3), ind_loop(4)
           lij  =       inddm( ind_loop(1) , j, k)
!DEC$ IVDEP
           do l = lij, lij +  ind_loop(2)- ind_loop(1)

             amor = MIN(celln(l),2.-celln(l))
             !amor = 1

             rop(l+v1) = rop(l+v1) + param_real(ROTATION)*amor

           enddo
        enddo
        enddo

      else

#include   "FastS/Compute/loop_begin.for"
             amor = MIN(celln(l),2.-celln(l))
             !amor = 1.
             rop(l+v1) = rop(l+v1) + param_real(ROTATION)*amor
#include   "FastS/Compute/loop_end.for"

      endif

      end
