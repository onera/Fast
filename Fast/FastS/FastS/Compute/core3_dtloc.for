c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine core3_dtloc(ndom, nitcfg,  param_int, param_real,
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



      REAL_E cgamma1, cgamma2, cgamma3,cgamma4
      parameter( cgamma1 =  0.5             )
      parameter( cgamma2 =  0.9106836025229591 )
      parameter( cgamma3 =  0.3660254037844387 )
      parameter( cgamma4 =  3.6427344100918355 )

      include 'omp_lib.h'

      INTEGER_E ndom,nitcfg,ind_loop(6), param_int(0:*)
 
      REAL_E     rop( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E   rop_1( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E   drodm( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E     coe( param_int(NDIMDX) * param_int(NEQ_COE) )

      REAL_E param_real(0:*)
C Var loc
      INTEGER_E incmax,l,i,j,k,lij,v1,v2,v3,v4,v5,v6,b,ind,ltij,lt,lvo
      REAL_E coefW,ro_old,u_old,v_old,w_old,t_old,nu_old,r_1,cvinv,
     &  roe_old, cv, cvinv2, tab1,tab2,tab3,tab4,tab5,tab6,
     & rho,rhou,rhov,rhoe,anulam,temp01,cmus1,coesut,t1
      REAL_E a1,a2,a3,a4,a5,a6,a7,a1_,a2_,a3_,a4_,a5_,a6_,a7_
      INTEGER_E b1,b2,b3,b4,b5,b6,b7,shift

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"


      
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
c-----Determination des coefficients : Explicit local instationnaire
c---------------------------------------------------------------------

         
         if (MOD(nitcfg,ind)==1) then
            coefW = cgamma1/float(param_int(LEVEL))
            !print*,'coeff, ndom = ',coefW, ndom
         elseif (MOD(nitcfg,ind)==ind/2) then
            coefW = cgamma2/float(param_int(LEVEL))
            !print*,'coeff, ndom = ',coefW, ndom
         elseif (MOD(nitcfg,ind)==ind-1) then
            coefW = cgamma3/float(param_int(LEVEL))
            !print*,'coeff, ndom = ',coefW, ndom
         else if (MOD(nitcfg,ind)==ind/4.and.ind.ne.4) then
            coefW = cgamma4/float(param_int(LEVEL))
            !print*,'coeff= ',coefW
                    
         end if 


c---------------------------------------------------------------------


      IF(param_int(NEQ).eq.5) THEN



         If (param_int(ITYPZONE).ne.3) then


#include  "FastC/HPC_LAYER/loop_begin.for"


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

c      if (k.eq.35.and.j.eq.35.and.l-lij.eq.24.and.ndom.eq.0)  then
c      write(*,*)'ro= ',rop_1(l +v1), tab1, drodm(l +v1),ro_old, nitcfg,l
c      end if

          rop_1(l +v3)     = (ro_old*v_old + tab3)*r_1
          rop_1(l +v4)     = (ro_old*w_old + tab4)*r_1
 
          roe_old        = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                             +v_old*v_old
     &                                             +w_old*w_old))

          rop_1(l +v5)     = ( roe_old + tab5 )*r_1*cvinv
     &                    - cvinv2*( rop_1(l +v2)*rop_1(l +v2)
     &                              +rop_1(l +v3)*rop_1(l +v3)
     &                              +rop_1(l +v4)*rop_1(l +v4))

            ! if (k==100.and.j==100.and.l-lij==50) then
            !    print*, 'coeff=',coefW, ndom
            ! end if

          !if (ndom==53.and.nitcfg==3) then
          !   print*, rop_1(l+v5), k
          !end if

          !if () then
          !   print*, rop_1(l+v4)
          !end if

#include  "FastC/HPC_LAYER/loop_end.for"
         


         Else


#include  "FastC/HPC_LAYER/loop_begin.for"


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


          !if (j==10.and. ndom==1) then          
          !   print*, 'coucou'
          !end if
         !if (j==50 .and. ndom==1 .and. nitcfg.eq.4) then
         !      print*, 'coef= ', coefW
         !end if

         !if (j==50 .and. ndom==2 .and. nitcfg.eq.4) then
         !      print*, 'coef= ', coefW
         !end if

          !if (ndom==0 .and. l-lij==100.and.j==200) then
          !   print*, tab1
          !end if


        !  print*, rop_1(l +v2)

      !if (ndom == 0.and.j==110.and.l-lij==399) then
      !   print*,'ro= ', ro_old, rop_1(l +v1), param_int(LEVEL)
      !end if

          
      !if (ndom == 0.and.j==110.and.l-lij==400) then
      !   print*,'ro= ', ro_old, rop_1(l +v1), param_int(LEVEL)
      !end if
      !if (ndom == 0.and.j==59.and.l-lij==87) then !86
      !   print*,'ro= ',rop(l+v5), drodm(l +v5), rop_1(l +v5)
      !end if


#include  "FastC/HPC_LAYER/loop_end.for"


         Endif
      
      ELSE

      cmus1  =    param_real(VISCO+4)
      temp01 = 1./param_real(VISCO+3)
      coesut =    param_real(VISCO+2) * (1.+cmus1*temp01)

         If (param_int(ITYPZONE).ne.3) then

#include  "FastC/HPC_LAYER/loop_begin.for"
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

#include  "FastC/HPC_LAYER/loop_end.for"

         Else

#include  "FastC/HPC_LAYER/loop_begin.for"
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

#include  "FastC/HPC_LAYER/loop_end.for"


        Endif! typezone
      ENDIF! neq

      end
