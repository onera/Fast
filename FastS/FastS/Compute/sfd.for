      subroutine sfd(param_int, param_real, nitrun,
     &               ind_loop,
     &               rop, rof, coe, vol)
c***********************************************************************
c_U   USER : DANDOIS
c
c     ACT
c_A    Selective Frequency Damping
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ind_loop(6), param_int(0:*), nitrun
      REAL_E  cfl(3), param_real(0:*)
 
      REAL_E rop( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E rof( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E coe( param_int(NDIMDX) * param_int(NEQ_COE) ) ! dt/vol
      REAL_E vol( param_int(NDIMDX_MTR) )

C Var loc
      INTEGER_E i,j,k,l,ltij,lvo,vcoe,nvar,ne,lij,lt 
      REAL_E delta,chi,xsi,cd

      INTEGER_E ndom

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      delta = param_real(SFD_DELTA)
      chi   = param_real(SFD_CHI)
      cd    = chi*delta

      vcoe = 0

      if(param_int(ITYPZONE).eq.2) then

        lvo  = 1
        do  ne = 1,param_int(NEQ)
          nvar = (ne-1)*param_int(NDIMDX)
          do  k = ind_loop(5), ind_loop(6)
          do  j = ind_loop(3), ind_loop(4)
              lij  =       inddm( ind_loop(1) , j, k)
!DEC$ IVDEP
              do l = lij, lij +  ind_loop(2)- ind_loop(1)

               xsi = exp(-(chi+1./delta)*coe(l + vcoe)*vol(lvo))

               rof(l +nvar) = 1./(1.+cd)*( (1.-xsi) * rop(l +nvar)
     &                                   + (cd+xsi) * rof(l +nvar))
              enddo
          enddo
          enddo
        enddo

      else

        do  ne = 1,param_int(NEQ)
            nvar = (ne-1)*param_int(NDIMDX)
#include   "FastS/Compute/loop_begin.for"
               xsi = exp(-(chi+1./delta)*coe(l + vcoe)*vol(lvo))

               rof(l +nvar) = 1./(1.+cd)*( (1.-xsi) * rop(l +nvar)
     &                                   + (cd+xsi) * rof(l +nvar))
#include   "FastS/Compute/loop_end.for"
        enddo

      endif

      IF (nitrun.ge.param_int(SFD_INIT_IT)) THEN

       if(param_int(ITYPZONE).eq.2) then

        lvo  = 1
        do  ne = 1,param_int(NEQ)
          nvar = (ne-1)*param_int(NDIMDX)
          do  k = ind_loop(5), ind_loop(6)
          do  j = ind_loop(3), ind_loop(4)
              lij  =       inddm( ind_loop(1) , j, k)
!DEC$ IVDEP
              do l = lij, lij +  ind_loop(2)- ind_loop(1)

               xsi = exp(-(chi+1./delta)*coe(l + vcoe)*vol(lvo))

               rop(l +nvar)= 1./(xsi+cd)*( (1.+cd) * xsi * rop(l +nvar)
     &                                   + (cd*(1.-xsi)) * rof(l +nvar))
              enddo
          enddo
          enddo
        enddo

       else

        do  ne = 1,param_int(NEQ)
            nvar = (ne-1)*param_int(NDIMDX)
#include   "FastS/Compute/loop_begin.for"
               xsi = exp(-(chi+1./delta)*coe(l + vcoe)*vol(lvo))

               rop(l +nvar)= 1./(xsi+cd)*( (1.+cd) * xsi * rop(l +nvar)
     &                                   + (cd*(1.-xsi)) * rof(l +nvar))
#include   "FastS/Compute/loop_end.for"
        enddo

       endif ! typezone

      ENDIF !nitrun

      end
