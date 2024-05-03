c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine core3ark3(ndom, nitcfg,  param_int, param_real,
     &                     ind_loop,
     &                     rop_1, rop, drodm, coe)
c***********************************************************************
c_U   USER : PECHIER
c
c     ACT
c_A    Mise a jour du residu explicite
c
c     VAL
c
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      REAL_E cgamma1, cgamma2, cgamma3
      parameter( cgamma1 =  0.5             )
      parameter( cgamma2 =  0.9106836025229591 )
      parameter( cgamma3 =  0.3660254037844387 )

      include 'omp_lib.h'

      INTEGER_E ndom,nitcfg,ind_loop(6), param_int(0:*)
 
      REAL_E     rop( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E   rop_1( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E   drodm( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E     coe( param_int(NDIMDX) * param_int(NEQ_COE) )

      REAL_E param_real(0:*)
C Var loc
      INTEGER_E incmax,l,i,j,k,lij,v1,v2,v3,v4,v5,v6,b,ind
      REAL_E coefW,ro_old,u_old,v_old,w_old,t_old,nu_old,r_1,cvinv,
     &  roe_old, cv, cvinv2, tab1,tab2,tab3,tab4,tab5,tab6,
     & rho,rhou,rhov,rhoe,anulam,temp01,cmus1,coesut,t1

#include "FastS/formule_param.h"

      
      cv     = param_real( CVINF)   ! Cv
      cvinv  = 1./cv   
      cvinv2 = 0.5*cvinv

      v1 = 0
      v2 =   param_int(NDIMDX)
      v3 = 2*param_int(NDIMDX)
      v4 = 3*param_int(NDIMDX)
      v5 = 4*param_int(NDIMDX)
      v6 = 5*param_int(NDIMDX)

      b=param_int(RK)

      ind=param_int(NSSITER)/param_int(LEVEL)
    
c---------------------------------------------------------------------
c-----Determination des coefficients : Explicit global
c---------------------------------------------------------------------
      if (param_int(EXPLOC)==0) then            

         if (b==1) then ! Runge-Kutta ordre 1
            if (nitcfg.eq.1) then
               coefW = 1.
            end if
         elseif (b==2) then ! Runge-Kutta ordre 2
            if (nitcfg.eq.1) then
               coefW = 1.
            elseif (nitcfg.eq.2) then
               coefW = 0.5
            end if
         elseif (b==3) then ! Runge-Kutta ordre 3
            if (nitcfg.eq.1) then
               coefW = cgamma1
            elseif (nitcfg.eq.2) then
               coefW = cgamma2
            elseif (nitcfg.eq.3) then
               coefW = cgamma3
            end if
         
        end if
c---------------------------------------------------------------------
c-----Determination des coefficients : Explicit local (RK 2)
c---------------------------------------------------------------------

      elseif (param_int(EXPLOC)==1) then 
      
         
         if (MOD(nitcfg,ind)==0) then               
            coefW = 0.5/float(param_int(LEVEL))
            !print*, coefW
         elseif (MOD(nitcfg,ind)==1) then 
            coefW = 1./float(param_int(LEVEL))
            !print*,coefW
         elseif (MOD(nitcfg,ind)==2) then 
            coefW = 0.25/float(param_int(LEVEL))
           ! print*, coefW
         end if

      end if
c---------------------------------------------------------------------


      IF(param_int(NEQ).eq.5) THEN

         If (param_int(ITYPZONE).ne.3) then

#include  "FastS/Compute/loop_core_begin.for"
          r_1    = coefW*coe(l +v1)

          tab1 =  drodm(l +v1)*r_1
          tab2 =  drodm(l +v2)*r_1
          tab3 =  drodm(l +v3)*r_1
          tab4 =  drodm(l +v4)*r_1
          tab5 =  drodm(l +v5)*r_1
c
          ro_old      = rop(l +v1)
          u_old       = rop(l +v2)
          v_old       = rop(l +v3)
          w_old       = rop(l +v4)
          t_old       = rop(l +v5)

          rop_1(l +v1)     = rop(l +v1) + tab1

          r_1            = 1./rop_1(l +v1)

          rop_1(l +v2)     = (ro_old*u_old + tab2)*r_1

        !if (k==2 .and. j==8 .and. l-lij == 112 .and. ndom==0)
     &  !then
        !       print*, 'ro= ', rop_1(l +v1)
        ! end if

          rop_1(l +v3)     = (ro_old*v_old + tab3)*r_1
          rop_1(l +v4)     = (ro_old*w_old + tab4)*r_1
 
          roe_old        = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                             +v_old*v_old
     &                                             +w_old*w_old))

          rop_1(l +v5)     = ( roe_old + tab5 )*r_1*cvinv
     &                    - cvinv2*( rop_1(l +v2)*rop_1(l +v2)
     &                              +rop_1(l +v3)*rop_1(l +v3)
     &                              +rop_1(l +v4)*rop_1(l +v4))
#include  "FastS/Compute/loop_end.for"

         Else

#include  "FastS/Compute/loop_core_begin.for"
          r_1    = coefW*coe(l +v1)

          tab1 =  drodm(l +v1)*r_1
          tab2 =  drodm(l +v2)*r_1
          tab3 =  drodm(l +v3)*r_1
          tab5 =  drodm(l +v5)*r_1
c
          ro_old      = rop(l +v1)
          u_old       = rop(l +v2)
          v_old       = rop(l +v3)
          t_old       = rop(l +v5)

          rop_1(l +v1)     = rop(l +v1) + tab1

          r_1            = 1./rop_1(l +v1)

          rop_1(l +v2)     = (ro_old*u_old + tab2)*r_1
          rop_1(l +v3)     = (ro_old*v_old + tab3)*r_1
          rop_1(l +v4)     = 0.
 
          roe_old        = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                             +v_old*v_old))

          rop_1(l +v5)     = ( roe_old + tab5 )*r_1*cvinv
     &                    - cvinv2*( rop_1(l +v2)*rop_1(l +v2)
     &                              +rop_1(l +v3)*rop_1(l +v3))
#include  "FastS/Compute/loop_end.for"
         Endif
      
      ELSE

      cmus1  =    param_real(VISCO+4)
      temp01 = 1./param_real(VISCO+3)
      coesut =    param_real(VISCO+2) * (1.+cmus1*temp01)

         If (param_int(ITYPZONE).ne.3) then

#include  "FastS/Compute/loop_core_begin.for"
          r_1    = coefW*coe(l +v1)

          tab1 =  drodm(l +v1)*r_1
          tab2 =  drodm(l +v2)*r_1
          tab3 =  drodm(l +v3)*r_1
          tab4 =  drodm(l +v4)*r_1
          tab5 =  drodm(l +v5)*r_1
          tab6 =  drodm(l +v6)*r_1
c
          ro_old      = rop(l +v1)
          u_old       = rop(l +v2)
          v_old       = rop(l +v3)
          w_old       = rop(l +v4)
          t_old       = rop(l +v5)
          nu_old      = rop(l +v6)

          rop_1(l +v1)     = rop(l +v1) + tab1

          r_1            = 1./rop_1(l +v1)

          rop_1(l +v2)     = (ro_old*u_old + tab2)*r_1
          rop_1(l +v3)     = (ro_old*v_old + tab3)*r_1
          rop_1(l +v4)     = (ro_old*w_old + tab4)*r_1


          roe_old        = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                             +v_old*v_old
     &                                             +w_old*w_old))

          rop_1(l +v5)     = ( roe_old + tab5 )*r_1*cvinv
     &                    - cvinv2*( rop_1(l +v2)*rop_1(l +v2)
     &                              +rop_1(l +v3)*rop_1(l +v3)
     &                              +rop_1(l +v4)*rop_1(l +v4))

          anulam=coesut*sqrt(rop_1(l+v5)*temp01)/(1.+cmus1/rop_1(l +v5))
          anulam=anulam*r_1*param_real(SA_REAL+SA_RATIOM-1) 

          t_old        = (ro_old*nu_old + tab6)*r_1
          t_old        = min(t_old,anulam)
          rop_1(l +v6) = max(t_old,0.)
#include  "FastS/Compute/loop_end.for"
         Else

#include  "FastS/Compute/loop_core_begin.for"
          r_1    = coefW*coe(l +v1)

          tab1 =  drodm(l +v1)*r_1
          tab2 =  drodm(l +v2)*r_1
          tab3 =  drodm(l +v3)*r_1
          tab5 =  drodm(l +v5)*r_1
          tab6 =  drodm(l +v6)*r_1
c
          ro_old      = rop(l +v1)
          u_old       = rop(l +v2)
          v_old       = rop(l +v3)
          t_old       = rop(l +v5)
          nu_old      = rop(l +v6)

          rop_1(l +v1)     = rop(l +v1) + tab1

          r_1            = 1./rop_1(l +v1)

          rop_1(l +v2)     = (ro_old*u_old + tab2)*r_1
          rop_1(l +v3)     = (ro_old*v_old + tab3)*r_1
          rop_1(l +v4)     = 0.

          roe_old        = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                             +v_old*v_old))

          rop_1(l +v5)     = ( roe_old + tab5 )*r_1*cvinv
     &                    - cvinv2*( rop_1(l +v2)*rop_1(l +v2)
     &                              +rop_1(l +v3)*rop_1(l +v3))

          anulam=coesut*sqrt(rop_1(l+v5)*temp01)/(1.+cmus1/rop_1(l +v5))
          anulam=anulam*r_1*param_real(SA_REAL+SA_RATIOM-1) 

          t_old        = (ro_old*nu_old + tab6)*r_1
          t_old        = min(t_old,anulam)
          !rop_1(l +v6) = t_old
          rop_1(l +v6) = max(t_old,0.)
#include  "FastS/Compute/loop_end.for"

        Endif! typezone
      ENDIF! neq

      end
