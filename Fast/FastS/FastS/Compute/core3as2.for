c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 39 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine core3as2(ndom,nitcfg, first_it, param_int, param_real, 
     &                    ind_loop,
     &                    rop, rop_n, rop_n1, drodm, coe)
c***********************************************************************
c_U   USER : PECHIER
c
c     ACT
c_A    Mise a jour du residu Newton
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom, nitcfg, first_it, ind_loop(6), param_int(0:*)
 
      REAL_E     rop( param_int(NDIMDX) , param_int(NEQ) )
      REAL_E   rop_n( param_int(NDIMDX) , param_int(NEQ) )
      REAL_E  rop_n1( param_int(NDIMDX) , param_int(NEQ) )
      REAL_E   drodm( param_int(NDIMDX) , param_int(NEQ) )
      REAL_E     coe( param_int(NDIMDX) , param_int(NEQ_COE) )

      REAL_E param_real(0:*)
 
C Var loc
      INTEGER_E incmax,l,i,j,k,lij,ltij,lt,lvo
      REAL_E c1,c2,c3,c4,roe,roe_n,roe_n1,cv

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      cv = param_real( CVINF )

      if(first_it.eq.0.or.param_int(DTLOC).eq.1) then !bdf1 au premier pas de temps ou dtloc
        c1=-1.
        c2= 1.
        c3= 0.
      else                    !bdf2 a partir du 2eme pas de temps
        c1=-1.5
        c2= 2.
        c3=-0.5
      endif

      !IF(iflow(ndom).eq.3.and.les(ndom).eq.0) THEN
      IF(param_int(NEQ).eq.6) THEN

       If (nitcfg.eq.1) then

         c4= c1+c2 !car ro(P+1)=ro(n) a la 1ere ss_iter

         if (param_int(ITYPZONE).ne.3) then

           if(param_int(DTLOC).eq.1) then
#include  "FastC/HPC_LAYER/loop_begin.for"
            drodm(l,1)= coe(l,1)*drodm(l,1) 
            drodm(l,2)= coe(l,1)*drodm(l,2) 
            drodm(l,3)= coe(l,1)*drodm(l,3) 
            drodm(l,4)= coe(l,1)*drodm(l,4) 
            drodm(l,5)= coe(l,1)*drodm(l,5)
            drodm(l,6)= coe(l,1)*drodm(l,6)
#include  "FastC/HPC_LAYER/loop_end.for"

           else

#include  "FastC/HPC_LAYER/loop_begin.for"
            drodm(l,1)= coe(l,1)*drodm(l,1) 
     &                     + c4*rop   (l,1)
     &                     + c3*rop_n1(l,1)

            drodm(l,2)= coe(l,1)*drodm(l,2) 
     &                     + c4*rop   (l,2)*rop   (l,1)
     &                     + c3*rop_n1(l,2)*rop_n1(l,1)

            drodm(l,3)= coe(l,1)*drodm(l,3) 
     &                     + c4*rop   (l,3)*rop   (l,1)
     &                     + c3*rop_n1(l,3)*rop_n1(l,1)

            drodm(l,4)= coe(l,1)*drodm(l,4) 
     &                     + c4*rop   (l,4)*rop   (l,1)
     &                     + c3*rop_n1(l,4)*rop_n1(l,1)

            roe =  rop(l,1)*(cv*rop(l,5) + 0.5*( rop(l,2)*rop(l,2)
     &                                          +rop(l,3)*rop(l,3)
     &                                          +rop(l,4)*rop(l,4)))

            roe_n1= rop_n1(l,1)*( cv*rop_n1(l,5) 
     &                           + 0.5*( rop_n1(l,2)*rop_n1(l,2)
     &                                  +rop_n1(l,3)*rop_n1(l,3)
     &                                  +rop_n1(l,4)*rop_n1(l,4)))

            drodm(l,5)= coe(l,1)*drodm(l,5) + c4*roe + c3*roe_n1

            drodm(l,6)= coe(l,1)*drodm(l,6) 
     &                     + c4*rop   (l,6)*rop   (l,1)
     &                     + c3*rop_n1(l,6)*rop_n1(l,1)
#include  "FastC/HPC_LAYER/loop_end.for"

           endif ! optim bdf1

        else !dom 2d

           if(param_int(DTLOC).eq.1) then
#include  "FastC/HPC_LAYER/loop_begin.for"
            drodm(l,1)= coe(l,1)*drodm(l,1) 
            drodm(l,2)= coe(l,1)*drodm(l,2) 
            drodm(l,3)= coe(l,1)*drodm(l,3) 
            drodm(l,5)= coe(l,1)*drodm(l,5)
            drodm(l,6)= coe(l,1)*drodm(l,6)
#include  "FastC/HPC_LAYER/loop_end.for"

           else

#include  "FastC/HPC_LAYER/loop_begin.for"
            drodm(l,1)= coe(l,1)*drodm(l,1) 
     &                     + c4*rop   (l,1)
     &                     + c3*rop_n1(l,1)

            drodm(l,2)= coe(l,1)*drodm(l,2) 
     &                     + c4*rop   (l,2)*rop   (l,1)
     &                     + c3*rop_n1(l,2)*rop_n1(l,1)

            drodm(l,3)= coe(l,1)*drodm(l,3) 
     &                     + c4*rop   (l,3)*rop   (l,1)
     &                     + c3*rop_n1(l,3)*rop_n1(l,1)

            roe =  rop(l,1)*(cv*rop(l,5) + 0.5*( rop(l,2)*rop(l,2)
     &                                          +rop(l,3)*rop(l,3)))


            roe_n1= rop_n1(l,1)*( cv*rop_n1(l,5) 
     &                           + 0.5*( rop_n1(l,2)*rop_n1(l,2)
     &                                  +rop_n1(l,3)*rop_n1(l,3)))

            drodm(l,5)= coe(l,1)*drodm(l,5) + c4*roe + c3*roe_n1

            drodm(l,6)= coe(l,1)*drodm(l,6) 
     &                     + c4*rop   (l,6)*rop   (l,1)
     &                     + c3*rop_n1(l,6)*rop_n1(l,1)
#include  "FastC/HPC_LAYER/loop_end.for"
           endif ! optim bdf1
        endif!2d/3d

       Else  !nitcfg > 1

         if (param_int(ITYPZONE).ne.3) then
#include  "FastC/HPC_LAYER/loop_begin.for"
            drodm(l,1)= coe(l,1)*drodm(l,1) 
     &                     + c1*rop   (l,1)
     &                     + c2*rop_n (l,1)
     &                     + c3*rop_n1(l,1)

            drodm(l,2)= coe(l,1)*drodm(l,2) 
     &                     + c1*rop   (l,2)*rop   (l,1)
     &                     + c2*rop_n (l,2)*rop_n (l,1)
     &                     + c3*rop_n1(l,2)*rop_n1(l,1)

            drodm(l,3)= coe(l,1)*drodm(l,3) 
     &                     + c1*rop   (l,3)*rop   (l,1)
     &                     + c2*rop_n (l,3)*rop_n (l,1)
     &                     + c3*rop_n1(l,3)*rop_n1(l,1)

            drodm(l,4)= coe(l,1)*drodm(l,4) 
     &                     + c1*rop   (l,4)*rop   (l,1)
     &                     + c2*rop_n (l,4)*rop_n (l,1)
     &                     + c3*rop_n1(l,4)*rop_n1(l,1)

            roe =  rop(l,1)*(cv*rop(l,5) + 0.5*( rop(l,2)*rop(l,2)
     &                                          +rop(l,3)*rop(l,3)
     &                                          +rop(l,4)*rop(l,4)))

            roe_n = rop_n(l,1)*( cv*rop_n(l,5) 
     &                          + 0.5*( rop_n(l,2)*rop_n(l,2)
     &                                 +rop_n(l,3)*rop_n(l,3)
     &                                 +rop_n(l,4)*rop_n(l,4)))

            roe_n1= rop_n1(l,1)*( cv*rop_n1(l,5) 
     &                           + 0.5*( rop_n1(l,2)*rop_n1(l,2)
     &                                  +rop_n1(l,3)*rop_n1(l,3)
     &                                  +rop_n1(l,4)*rop_n1(l,4)))

            drodm(l,5)= coe(l,1)*drodm(l,5) +c1*roe +c2*roe_n +c3*roe_n1

            drodm(l,6)= coe(l,1)*drodm(l,6) 
     &                     + c1*rop   (l,6)*rop   (l,1)
     &                     + c2*rop_n (l,6)*rop_n (l,1)
     &                     + c3*rop_n1(l,6)*rop_n1(l,1)
#include  "FastC/HPC_LAYER/loop_end.for"

        else !dom 2d

#include  "FastC/HPC_LAYER/loop_begin.for"
            drodm(l,1)= coe(l,1)*drodm(l,1) 
     &                     + c1*rop   (l,1)
     &                     + c2*rop_n (l,1)
     &                     + c3*rop_n1(l,1)

            drodm(l,2)= coe(l,1)*drodm(l,2) 
     &                     + c1*rop   (l,2)*rop   (l,1)
     &                     + c2*rop_n (l,2)*rop_n (l,1)
     &                     + c3*rop_n1(l,2)*rop_n1(l,1)

            drodm(l,3)= coe(l,1)*drodm(l,3) 
     &                     + c1*rop   (l,3)*rop   (l,1)
     &                     + c2*rop_n (l,3)*rop_n (l,1)
     &                     + c3*rop_n1(l,3)*rop_n1(l,1)

            roe =  rop(l,1)*(cv*rop(l,5) + 0.5*( rop(l,2)*rop(l,2)
     &                                          +rop(l,3)*rop(l,3)))

            roe_n = rop_n(l,1)*( cv*rop_n(l,5) 
     &                          + 0.5*( rop_n(l,2)*rop_n(l,2)
     &                                 +rop_n(l,3)*rop_n(l,3)))

            roe_n1= rop_n1(l,1)*( cv*rop_n1(l,5) 
     &                           + 0.5*( rop_n1(l,2)*rop_n1(l,2)
     &                                  +rop_n1(l,3)*rop_n1(l,3)))

            drodm(l,5)= coe(l,1)*drodm(l,5) +c1*roe +c2*roe_n +c3*roe_n1

            drodm(l,6)= coe(l,1)*drodm(l,6) 
     &                     + c1*rop   (l,6)*rop   (l,1)
     &                     + c2*rop_n (l,6)*rop_n (l,1)
     &                     + c3*rop_n1(l,6)*rop_n1(l,1)
#include  "FastC/HPC_LAYER/loop_end.for"

        endif !2d/3d
       Endif  !nitcfg=1 ou >1
      
      ELSE !neq=5

      If (nitcfg.eq.1) then

         c4= c1+c2 !car ro(P+1)=ro(n) a la 1ere ss_iter

         if (param_int(ITYPZONE).ne.3) then
#include  "FastC/HPC_LAYER/loop_begin.for"
            drodm(l,1)= coe(l,1)*drodm(l,1) 
     &                     + c4*rop   (l,1)
     &                     + c3*rop_n1(l,1)

            drodm(l,2)= coe(l,1)*drodm(l,2) 
     &                     + c4*rop   (l,2)*rop   (l,1)
     &                     + c3*rop_n1(l,2)*rop_n1(l,1)

            drodm(l,3)= coe(l,1)*drodm(l,3) 
     &                     + c4*rop   (l,3)*rop   (l,1)
     &                     + c3*rop_n1(l,3)*rop_n1(l,1)

            drodm(l,4)= coe(l,1)*drodm(l,4) 
     &                     + c4*rop   (l,4)*rop   (l,1)
     &                     + c3*rop_n1(l,4)*rop_n1(l,1)

            roe =  rop(l,1)*(cv*rop(l,5) + 0.5*( rop(l,2)*rop(l,2)
     &                                          +rop(l,3)*rop(l,3)
     &                                          +rop(l,4)*rop(l,4)))

            roe_n1= rop_n1(l,1)*( cv*rop_n1(l,5) 
     &                           + 0.5*( rop_n1(l,2)*rop_n1(l,2)
     &                                  +rop_n1(l,3)*rop_n1(l,3)
     &                                  +rop_n1(l,4)*rop_n1(l,4)))

            drodm(l,5)= coe(l,1)*drodm(l,5) + c4*roe + c3*roe_n1
#include  "FastC/HPC_LAYER/loop_end.for"

        else !dom 2d
#include  "FastC/HPC_LAYER/loop_begin.for"
            drodm(l,1)= coe(l,1)*drodm(l,1) 
     &                     + c4*rop   (l,1)
     &                     + c3*rop_n1(l,1)

            drodm(l,2)= coe(l,1)*drodm(l,2) 
     &                     + c4*rop   (l,2)*rop   (l,1)
     &                     + c3*rop_n1(l,2)*rop_n1(l,1)


            drodm(l,3)= coe(l,1)*drodm(l,3) 
     &                     + c4*rop   (l,3)*rop   (l,1)
     &                     + c3*rop_n1(l,3)*rop_n1(l,1)


            roe =  rop(l,1)*(cv*rop(l,5) + 0.5*( rop(l,2)*rop(l,2)
     &                                          +rop(l,3)*rop(l,3)))

            roe_n1= rop_n1(l,1)*( cv*rop_n1(l,5) 
     &                           + 0.5*( rop_n1(l,2)*rop_n1(l,2)
     &                                  +rop_n1(l,3)*rop_n1(l,3)))

            drodm(l,5)= coe(l,1)*drodm(l,5) + c4*roe + c3*roe_n1
#include  "FastC/HPC_LAYER/loop_end.for"

        endif!2d/3d

      ELSE!nitcfg > ou = 1

         if (param_int(ITYPZONE).ne.3) then
#include  "FastC/HPC_LAYER/loop_begin.for"
            drodm(l,1)= coe(l,1)*drodm(l,1) 
     &                     + c1*rop   (l,1)
     &                     + c2*rop_n (l,1)
     &                     + c3*rop_n1(l,1)

            drodm(l,2)= coe(l,1)*drodm(l,2) 
     &                     + c1*rop   (l,2)*rop   (l,1)
     &                     + c2*rop_n (l,2)*rop_n (l,1)
     &                     + c3*rop_n1(l,2)*rop_n1(l,1)

            drodm(l,3)= coe(l,1)*drodm(l,3) 
     &                     + c1*rop   (l,3)*rop   (l,1)
     &                     + c2*rop_n (l,3)*rop_n (l,1)
     &                     + c3*rop_n1(l,3)*rop_n1(l,1)

            drodm(l,4)= coe(l,1)*drodm(l,4) 
     &                     + c1*rop   (l,4)*rop   (l,1)
     &                     + c2*rop_n (l,4)*rop_n (l,1)
     &                     + c3*rop_n1(l,4)*rop_n1(l,1)

            roe =  rop(l,1)*(cv*rop(l,5) + 0.5*( rop(l,2)*rop(l,2)
     &                                          +rop(l,3)*rop(l,3)
     &                                          +rop(l,4)*rop(l,4)))


            roe_n = rop_n(l,1)*( cv*rop_n(l,5) 
     &                          + 0.5*( rop_n(l,2)*rop_n(l,2)
     &                                 +rop_n(l,3)*rop_n(l,3)
     &                                 +rop_n(l,4)*rop_n(l,4)))

            roe_n1= rop_n1(l,1)*( cv*rop_n1(l,5) 
     &                           + 0.5*( rop_n1(l,2)*rop_n1(l,2)
     &                                  +rop_n1(l,3)*rop_n1(l,3)
     &                                  +rop_n1(l,4)*rop_n1(l,4)))

            drodm(l,5)= coe(l,1)*drodm(l,5)+ c1*roe +c2*roe_n +c3*roe_n1
#include  "FastC/HPC_LAYER/loop_end.for"

        else !dom 2d
#include  "FastC/HPC_LAYER/loop_begin.for"
            drodm(l,1)= coe(l,1)*drodm(l,1) 
     &                     + c1*rop   (l,1)
     &                     + c2*rop_n (l,1)
     &                     + c3*rop_n1(l,1)

            drodm(l,2)= coe(l,1)*drodm(l,2) 
     &                     + c1*rop   (l,2)*rop   (l,1)
     &                     + c2*rop_n (l,2)*rop_n (l,1)
     &                     + c3*rop_n1(l,2)*rop_n1(l,1)

            drodm(l,3)= coe(l,1)*drodm(l,3) 
     &                     + c1*rop   (l,3)*rop   (l,1)
     &                     + c2*rop_n (l,3)*rop_n (l,1)
     &                     + c3*rop_n1(l,3)*rop_n1(l,1)

            roe =  rop(l,1)*(cv*rop(l,5) + 0.5*( rop(l,2)*rop(l,2)
     &                                          +rop(l,3)*rop(l,3)))

            roe_n = rop_n(l,1)*( cv*rop_n(l,5) 
     &                          + 0.5*( rop_n(l,2)*rop_n(l,2)
     &                                 +rop_n(l,3)*rop_n(l,3)))

            roe_n1= rop_n1(l,1)*( cv*rop_n1(l,5) 
     &                           + 0.5*( rop_n1(l,2)*rop_n1(l,2)
     &                                  +rop_n1(l,3)*rop_n1(l,3)))

            drodm(l,5)= coe(l,1)*drodm(l,5)+ c1*roe +c2*roe_n +c3*roe_n1
#include  "FastC/HPC_LAYER/loop_end.for"
        endif

      ENDIF!nitcfg
      ENDIF!neq

      end
