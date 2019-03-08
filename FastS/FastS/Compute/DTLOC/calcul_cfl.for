c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 57 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine calcul_cfl(ndom , param_int, param_real,
     &                  ind_loop,
     &                  cfl,
     &                  rop,xmut,venti, ti,tj,tk,vol,ipt_cfl)
c***********************************************************************
c_U   USER :  PECHIER
c
c     ACT
c_A    Calcul du CFL local max par domaine lors d'un calcul 
c_A    instationnaire.
c
c     VAL
c_V    Gaz parfait a une espece
c
c     OUT
c_O    common /infocfl/
c
c**********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom, param_int(0:*), ind_loop(6)

      REAL_E  xmut( param_int(NDIMDX) )
      REAL_E   rop( param_int(NDIMDX) , param_int(NEQ) )
      REAL_E   ipt_cfl(param_int(NDIMDX))

      REAL_E ti( param_int(NDIMDX_MTR) , param_int(NEQ_IJ) ),
     &       tj( param_int(NDIMDX_MTR) , param_int(NEQ_IJ) ),
     &       tk( param_int(NDIMDX_MTR) , param_int(NEQ_K ) ),
     &      vol( param_int(NDIMDX_MTR) ) 

      REAL_E venti( param_int(NDIMDX_VENT) * param_int(NEQ_VENT) )

      REAL_E  cfl(3), param_real(0:*)
C var local
      INTEGER_E l,ll,i,j,k,inci,incj,inck,lt,lij,ltij,
     & lvij,lv,v2ven,v3ven,inci_ven, lvo

      REAL_E xinvkt, u,v,w,r,ue,ve,we,
     &      ur,vr,wr,ur2,ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale,rgp,gam2,
     &      sp,c, xinvspc, vis, conv,comp, cfloc, u2,gam1,
     &      ddi,ddj,ddk,surf
 
#include "FastS/formule_mtr_param.h"
#include "FastS/formule_vent_param.h"
#include "FastS/formule_param.h" 

      inci = param_int(NIJK_MTR)
      incj = param_int(NIJK_MTR+1)
      inck = param_int(NIJK_MTR+2)

      rgp  = param_real(STAREF+1)*( param_real(STAREF) -1.)  !Cv(gama-1)= R (gas parfait)
      gam1 = param_real(STAREF)/param_real(VISCO +1)
      gam2 = param_real(STAREF)*rgp

c******* NS unsteady *****
      If (param_int(LALE).ge.1.and.(param_int(IFLOW).ge.2)) Then
         

         inci_ven = param_int(NIJK_VENT)

         if(param_int(NEQ_VENT).eq.2) then
            ck_vent =0.
         else
            ck_vent =1.
         endif
         v2ven =   param_int(NDIMDX_VENT)
         v3ven = 2*param_int(NDIMDX_VENT)*ck_vent


        if(param_int(ITYPZONE).eq.0) then

#include   "FastS/Compute/loop_ale_begin.for"
 
              ddi =sqrt( ti(lt+inci,1)*ti(lt+inci,1)
     &                  +ti(lt+inci,2)*ti(lt+inci,2)
     &                  +ti(lt+inci,3)*ti(lt+inci,3) )
              ddi =ddi+ sqrt( ti(lt,1)*ti(lt,1)+ti(lt,2)*ti(lt,2)
     &                       +ti(lt,3)*ti(lt,3) )

              ddj =sqrt( tj(lt+incj,1)*tj(lt+incj,1)
     &                  +tj(lt+incj,2)*tj(lt+incj,2)
     &                  +tj(lt+incj,3)*tj(lt+incj,3) )
              ddj =ddj+ sqrt( tj(lt,1)*tj(lt,1)+tj(lt,2)*tj(lt,2)
     &                       +tj(lt,3)*tj(lt,3) )

              ddk =sqrt( tk(lt+inck,1)*tk(lt+inck,1)
     &                  +tk(lt+inck,2)*tk(lt+inck,2)
     &                  +tk(lt+inck,3)*tk(lt+inck,3) )
              ddk =ddk+ sqrt( tk(lt,1)*tk(lt,1)+tk(lt,2)*tk(lt,2)
     &                       +tk(lt,3)*tk(lt,3) )

              surf= 1./max(ddi,ddj,ddk)

              u = rop(l,2)
              v = rop(l,3)
              w = rop(l,4)
              r = rop(l,1) 
              c = sqrt(gam2*rop(l,5))
 
              ue=.5*(venti(lv      )+venti(lv+inci_ven      ))
              ve=.5*(venti(lv+v2ven)+venti(lv+inci_ven+v2ven))
              we=.5*(venti(lv+v3ven)+venti(lv+inci_ven+v3ven))*ck_vent

              u  = u - ue
              v  = v - ve
              w  = w - we
              ur2 = u*u + v*v + w*w
              sp  = sqrt(ur2)
 
              xinvspc = 1./(sp+c)
              xinvkt  = gam1/xmut(l)

              vis   = 2.*r*vol(lvo)*vol(lvo)*surf*surf*xinvkt
              conv  = 2.*vol(lvo) * surf * xinvspc
              comp  = min(vis,conv)

              cfloc = 1.0/comp

              cfl(1) = max(cfl(1),cfloc)
              cfl(2) = min(cfl(2),cfloc)
              cfl(3) = cfl(3) + cfloc
              ipt_cfl(l)=cfloc
#include   "FastS/Compute/loop_end.for"

       elseif(param_int(ITYPZONE).eq.1) then

#include   "FastS/Compute/loop_ale_begin.for"
 
              ddi =sqrt( ti(lt+inci,1)*ti(lt+inci,1)
     &                  +ti(lt+inci,2)*ti(lt+inci,2) )
              ddi =ddi+ sqrt( ti(lt,1)*ti(lt,1)+ti(lt,2)*ti(lt,2) )

              ddj =sqrt( tj(lt+incj,1)*tj(lt+incj,1)
     &                  +tj(lt+incj,2)*tj(lt+incj,2) )
              ddj =ddj+ sqrt( tj(lt,1)*tj(lt,1)+tj(lt,2)*tj(lt,2) )

              ddk =sqrt( tk(lt+inck,1)*tk(lt+inck,1) )
              ddk =ddk+ sqrt( tk(lt,1)*tk(lt,1))

              surf= 1./max(ddi,ddj,ddk)

              u = rop(l,2)
              v = rop(l,3)
              w = rop(l,4)
              r = rop(l,1) 
              c = sqrt(gam2*rop(l,5))
 
              ue=.5*(venti(lv      )+venti(lv+inci_ven      ))
              ve=.5*(venti(lv+v2ven)+venti(lv+inci_ven+v2ven))
              we=.5*(venti(lv+v3ven)+venti(lv+inci_ven+v3ven))*ck_vent
 
              u  = u - ue
              v  = v - ve
              w  = w - we
              ur2 = u*u + v*v + w*w
              sp  = sqrt(ur2)
 
              xinvspc = 1./(sp+c)
              xinvkt  = gam1/xmut(l)

              vis   = 2.*r*vol(lvo)*vol(lvo)*surf*surf*xinvkt
              conv  = 2.*vol(lvo) * surf * xinvspc
              comp  = min(vis,conv)

              cfloc = 1.0/comp

              cfl(1) = max(cfl(1),cfloc)
              cfl(2) = min(cfl(2),cfloc)
              cfl(3) = cfl(3) + cfloc
              ipt_cfl(l)=cfloc
#include   "FastS/Compute/loop_end.for"

       elseif(param_int(ITYPZONE).eq.2) then

          ddi = 2.*sqrt( ti(1,1)*ti(1,1) )
          ddj = 2.*sqrt( tj(1,1)*tj(1,1) )
          ddk = 2.*sqrt( tk(1,1)*tk(1,1) )
          surf= 1./max(ddi,ddj,ddk)

          do 102 j = ind_loop(3), ind_loop(4)
             lij  =         inddm( ind_loop(1) , j, k)
             lvij =  lij-  indven( ind_loop(1) , j, k)
             lt   = 1
             do 102 l = lij,  lij + ind_loop(2) - ind_loop(1)

              lv = l - lvij
 
              u = rop(l,2)
              v = rop(l,3)
              w = rop(l,4)
              r = rop(l,1) 
              c = sqrt(gam2*rop(l,5))
 
              ue=.5*(venti(lv      )+venti(lv+inci_ven      ))
              ve=.5*(venti(lv+v2ven)+venti(lv+inci_ven+v2ven))
              we=.5*(venti(lv+v3ven)+venti(lv+inci_ven+v3ven))*ck_vent
 
              u  = u - ue
              v  = v - ve
              w  = w - we
              ur2 = u*u + v*v + w*w
              sp  = sqrt(ur2)
 
              xinvspc = 1./(sp+c)
              xinvkt  = gam1/xmut(l)

              vis   = 2.*r*vol(lvo)*vol(lvo)*surf*surf*xinvkt
              conv  = 2.*vol(lvo) * surf * xinvspc
              comp  = min(vis,conv)

              cfloc = 1.0/comp

              cfl(1) = max(cfl(1),cfloc)
              cfl(2) = min(cfl(2),cfloc)
              cfl(3) = cfl(3) + cfloc
              ipt_cfl(l)=cfloc
102   continue

       else

#include   "FastS/Compute/loop_ale_begin.for"
 
              ddi =sqrt( ti(lt+inci,1)*ti(lt+inci,1)
     &                  +ti(lt+inci,2)*ti(lt+inci,2) )
              ddi =ddi+ sqrt( ti(lt,1)*ti(lt,1)+ti(lt,2)*ti(lt,2) )

              ddj =sqrt( tj(lt+incj,1)*tj(lt+incj,1)
     &                  +tj(lt+incj,2)*tj(lt+incj,2) )
              ddj =ddj+ sqrt( tj(lt,1)*tj(lt,1)+tj(lt,2)*tj(lt,2) )

              surf= 1./max(ddi,ddj)

              u = rop(l,2)
              v = rop(l,3)
              r = rop(l,1) 
              c = sqrt(gam2*rop(l,5))
 
              ue=.5*(venti(lv      )+venti(lv+inci_ven      ))
              ve=.5*(venti(lv+v2ven)+venti(lv+inci_ven+v2ven))
 
              u  = u - ue
              v  = v - ve
              ur2 = u*u + v*v
              sp  = sqrt(ur2)
 
              xinvspc = 1./(sp+c)
              xinvkt  = gam1/xmut(l)

              vis   = 2.*r*vol(lvo)*vol(lvo)*surf*surf*xinvkt
              conv  = 2.*vol(lvo) * surf * xinvspc
              comp  = min(vis,conv)

              cfloc = 1.0/comp

              cfl(1) = max(cfl(1),cfloc)
              cfl(2) = min(cfl(2),cfloc)
              cfl(3) = cfl(3) + cfloc
              ipt_cfl(l)=cfloc
#include   "FastS/Compute/loop_end.for"

       endif

     
c******* NS steady *****

      Elseif (param_int(LALE).eq.0.and.param_int(IFLOW).ge.2) Then


        if(param_int(ITYPZONE).eq.0) then


#include   "FastS/Compute/loop_begin.for"
              ddi =sqrt( ti(lt+inci,1)*ti(lt+inci,1)
     &                  +ti(lt+inci,2)*ti(lt+inci,2)
     &                  +ti(lt+inci,3)*ti(lt+inci,3) )
              ddi =ddi+ sqrt( ti(lt,1)*ti(lt,1)+ti(lt,2)*ti(lt,2)
     &                       +ti(lt,3)*ti(lt,3) )

              ddj =sqrt( tj(lt+incj,1)*tj(lt+incj,1)
     &                  +tj(lt+incj,2)*tj(lt+incj,2)
     &                  +tj(lt+incj,3)*tj(lt+incj,3) )
              ddj =ddj+ sqrt( tj(lt,1)*tj(lt,1)+tj(lt,2)*tj(lt,2)
     &                       +tj(lt,3)*tj(lt,3) )

              ddk =sqrt( tk(lt+inck,1)*tk(lt+inck,1)
     &                  +tk(lt+inck,2)*tk(lt+inck,2)
     &                  +tk(lt+inck,3)*tk(lt+inck,3) )
              ddk =ddk+ sqrt( tk(lt,1)*tk(lt,1)+tk(lt,2)*tk(lt,2)
     &                       +tk(lt,3)*tk(lt,3) )

              surf= 1./max(ddi,ddj,ddk)

              u = rop(l,2)
              v = rop(l,3)
              w = rop(l,4)
              r = rop(l,1) 
              c = sqrt(gam2*rop(l,5))
 
              ur2 = u*u + v*v + w*w
              sp  = sqrt(ur2)
 
              xinvspc = 1./(sp+c)
              xinvkt  = gam1/xmut(l)

              vis   = 2.*r*vol(lvo)*vol(lvo)*surf*surf*xinvkt
              conv  = 2.*vol(lvo) * surf * xinvspc
              comp  = min(vis,conv)

              cfloc = 1.0/comp

              cfl(1) = max(cfl(1),cfloc)
              cfl(2) = min(cfl(2),cfloc)
              cfl(3) = cfl(3) + cfloc
              ipt_cfl(l)=cfloc
#include   "FastS/Compute/loop_end.for"


       elseif(param_int(ITYPZONE).eq.1) then



#include   "FastS/Compute/loop_begin.for"
              ddi =sqrt( ti(lt+inci,1)*ti(lt+inci,1)
     &                  +ti(lt+inci,2)*ti(lt+inci,2) )
              ddi =ddi+ sqrt( ti(lt,1)*ti(lt,1)+ti(lt,2)*ti(lt,2) )

              ddj =sqrt( tj(lt+incj,1)*tj(lt+incj,1)
     &                  +tj(lt+incj,2)*tj(lt+incj,2) )
              ddj =ddj+ sqrt( tj(lt,1)*tj(lt,1)+tj(lt,2)*tj(lt,2) )

              ddk =sqrt( tk(lt+inck,1)*tk(lt+inck,1) )
              ddk =ddk+ sqrt( tk(lt,1)*tk(lt,1))

              surf= 1./max(ddi,ddj,ddk)

              u = rop(l,2)
              v = rop(l,3)
              w = rop(l,4)
              r = rop(l,1) 
              c = sqrt(gam2*rop(l,5))
 
              ur2 = u*u + v*v + w*w
              sp  = sqrt(ur2)
 
              xinvspc = 1./(sp+c)
              xinvkt  = gam1/xmut(l)

              vis   = 2.*r*vol(lvo)*vol(lvo)*surf*surf*xinvkt
              conv  = 2.*vol(lvo) * surf * xinvspc
              comp  = min(vis,conv)

              cfloc = 1.0/comp

              cfl(1) = max(cfl(1),cfloc)
              cfl(2) = min(cfl(2),cfloc)
              cfl(3) = cfl(3) + cfloc
              ipt_cfl(l)=cfloc
#include   "FastS/Compute/loop_end.for"


       elseif(param_int(ITYPZONE).eq.2) then



          ddi = 2.*sqrt( ti(1,1)*ti(1,1) )
          ddj = 2.*sqrt( tj(1,1)*tj(1,1) )
          ddk = 2.*sqrt( tk(1,1)*tk(1,1) )
          surf= 1./max(ddi,ddj,ddk)

          do 202 k = ind_loop(5), ind_loop(6)
          do 202 j = ind_loop(3), ind_loop(4)
             lij  =        inddm( ind_loop(1) , j, k)
             lt   = 1
             lvo  = 1
             do 202 l = lij,  lij + ind_loop(2) - ind_loop(1)

              u = rop(l,2)
              v = rop(l,3)
              w = rop(l,4)
              r = rop(l,1) 
              c = sqrt(gam2*rop(l,5))
 
              ur2 = u*u + v*v + w*w
              sp  = sqrt(ur2)
 
              xinvspc = 1./(sp+c)
              xinvkt  = gam1/xmut(l)

              vis   = 2.*r*vol(lvo)*vol(lvo)*surf*surf*xinvkt
              conv  = 2.*vol(lvo) * surf * xinvspc
              comp  = min(vis,conv)

              cfloc = 1.0/comp

              cfl(1) = max(cfl(1),cfloc)
              cfl(2) = min(cfl(2),cfloc)
              cfl(3) = cfl(3) + cfloc
              ipt_cfl(l)=cfloc
202   continue

       else

 

#include   "FastS/Compute/loop_begin.for"
              ddi =sqrt( ti(lt+inci,1)*ti(lt+inci,1)
     &                  +ti(lt+inci,2)*ti(lt+inci,2) )
              ddi =ddi+ sqrt( ti(lt,1)*ti(lt,1)+ti(lt,2)*ti(lt,2) )

              ddj =sqrt( tj(lt+incj,1)*tj(lt+incj,1)
     &                  +tj(lt+incj,2)*tj(lt+incj,2) )
              ddj =ddj+ sqrt( tj(lt,1)*tj(lt,1)+tj(lt,2)*tj(lt,2) )

              surf= 1./max(ddi,ddj)

              u = rop(l,2)
              v = rop(l,3)
              r = rop(l,1) 
              c = sqrt(gam2*rop(l,5))
 
              ur2 = u*u + v*v
              sp  = sqrt(ur2)
 
              xinvspc = 1./(sp+c)
              xinvkt  = gam1/xmut(l)

              vis   = 2.*r*vol(lvo)*vol(lvo)*surf*surf*xinvkt
              conv  = 2.*vol(lvo) * surf * xinvspc
              comp  = min(vis,conv)

              cfloc = 1.0/comp

              cfl(1) = max(cfl(1),cfloc)
              cfl(2) = min(cfl(2),cfloc)
              cfl(3) = cfl(3) + cfloc
              ipt_cfl(l)=cfloc
#include   "FastS/Compute/loop_end.for"

       endif

c******* Euler unsteady *****

      Elseif (param_int(LALE).ge.1.and.param_int(IFLOW).lt.2) Then

         inci_ven = param_int(NIJK_VENT)


         if(param_int(NEQ_VENT).eq.2) then
            ck_vent =0.
         else
            ck_vent =1.
         endif
         v2ven =   param_int(NDIMDX_VENT)
         v3ven = 2*param_int(NDIMDX_VENT)*ck_vent


        if(param_int(ITYPZONE).eq.0) then



#include   "FastS/Compute/loop_ale_begin.for"
              ddi =sqrt( ti(lt+inci,1)*ti(lt+inci,1)
     &                  +ti(lt+inci,2)*ti(lt+inci,2)
     &                  +ti(lt+inci,3)*ti(lt+inci,3) )
              ddi =ddi+ sqrt( ti(lt,1)*ti(lt,1)+ti(lt,2)*ti(lt,2)
     &                       +ti(lt,3)*ti(lt,3) )

              ddj =sqrt( tj(lt+incj,1)*tj(lt+incj,1)
     &                  +tj(lt+incj,2)*tj(lt+incj,2)
     &                  +tj(lt+incj,3)*tj(lt+incj,3) )
              ddj =ddj+ sqrt( tj(lt,1)*tj(lt,1)+tj(lt,2)*tj(lt,2)
     &                       +tj(lt,3)*tj(lt,3) )

              ddk =sqrt( tk(lt+inck,1)*tk(lt+inck,1)
     &                  +tk(lt+inck,2)*tk(lt+inck,2)
     &                  +tk(lt+inck,3)*tk(lt+inck,3) )
              ddk =ddk+ sqrt( tk(lt,1)*tk(lt,1)+tk(lt,2)*tk(lt,2)
     &                       +tk(lt,3)*tk(lt,3) )

              surf= 1./max(ddi,ddj,ddk)

              u = rop(l,2)
              v = rop(l,3)
              w = rop(l,4)
              r = rop(l,1) 
              c = sqrt(gam2*rop(l,5))
 
              ue=.5*(venti(lv      )+venti(lv+inci_ven      ))
              ve=.5*(venti(lv+v2ven)+venti(lv+inci_ven+v2ven))
              we=.5*(venti(lv+v3ven)+venti(lv+inci_ven+v3ven))*ck_vent
 
              u  = u - ue
              v  = v - ve
              w  = w - we
              ur2 = u*u + v*v + w*w
              sp  = sqrt(ur2)
 
              xinvspc = 1./(sp+c)

              conv  = 2.*vol(lvo) * surf * xinvspc

              cfloc = 1.0/conv

              cfl(1) = max(cfl(1),cfloc)
              cfl(2) = min(cfl(2),cfloc)
              cfl(3) = cfl(3) + cfloc
              ipt_cfl(l)=cfloc
#include   "FastS/Compute/loop_end.for"

       elseif(param_int(ITYPZONE).eq.1) then


#include   "FastS/Compute/loop_ale_begin.for"
 
              ddi =sqrt( ti(lt+inci,1)*ti(lt+inci,1)
     &                  +ti(lt+inci,2)*ti(lt+inci,2) )
              ddi =ddi+ sqrt( ti(lt,1)*ti(lt,1)+ti(lt,2)*ti(lt,2) )

              ddj =sqrt( tj(lt+incj,1)*tj(lt+incj,1)
     &                  +tj(lt+incj,2)*tj(lt+incj,2) )
              ddj =ddj+ sqrt( tj(lt,1)*tj(lt,1)+tj(lt,2)*tj(lt,2) )

              ddk =sqrt( tk(lt+inck,1)*tk(lt+inck,1) )
              ddk =ddk+ sqrt( tk(lt,1)*tk(lt,1))

              surf= 1./max(ddi,ddj,ddk)

              u = rop(l,2)
              v = rop(l,3)
              w = rop(l,4)
              r = rop(l,1) 
              c = sqrt(gam2*rop(l,5))
 
              ue=.5*(venti(lv      )+venti(lv+inci_ven      ))
              ve=.5*(venti(lv+v2ven)+venti(lv+inci_ven+v2ven))
              we= 0.
 
              u  = u - ue
              v  = v - ve
              w  = w - we
              ur2 = u*u + v*v + w*w
              sp  = sqrt(ur2)
 
              xinvspc = 1./(sp+c)

              conv  = 2.*vol(lvo) * surf * xinvspc

              cfloc = 1.0/conv

              cfl(1) = max(cfl(1),cfloc)
              cfl(2) = min(cfl(2),cfloc)
              cfl(3) = cfl(3) + cfloc
              ipt_cfl(l)=cfloc
#include   "FastS/Compute/loop_end.for"

       elseif(param_int(ITYPZONE).eq.2) then

          ddi = 2.*sqrt( ti(1,1)*ti(1,1) )
          ddj = 2.*sqrt( tj(1,1)*tj(1,1) )
          ddk = 2.*sqrt( tk(1,1)*tk(1,1) )
          surf= 1./max(ddi,ddj,ddk)

          do 302 k = ind_loop(5), ind_loop(6)
          do 302 j = ind_loop(3), ind_loop(4)
             lij  =        inddm( ind_loop(1) , j, k)
             lt   = 1
             do 302 l = lij,  lij + ind_loop(2) - ind_loop(1)

              u = rop(l,2)
              v = rop(l,3)
              w = rop(l,4)
              r = rop(l,1) 
              c = sqrt(gam2*rop(l,5))
 
              ue=.5*(venti(lv      )+venti(lv+inci_ven      ))
              ve=.5*(venti(lv+v2ven)+venti(lv+inci_ven+v2ven))
              we=.5*(venti(lv+v3ven)+venti(lv+inci_ven+v3ven))*ck_vent
 
              u  = u - ue
              v  = v - ve
              w  = w - we
              ur2 = u*u + v*v + w*w
              sp  = sqrt(ur2)
 
              xinvspc = 1./(sp+c)

              conv  = 2.*vol(lvo) * surf * xinvspc

              cfloc = 1.0/conv

              cfl(1) = max(cfl(1),cfloc)
              cfl(2) = min(cfl(2),cfloc)
              cfl(3) = cfl(3) + cfloc
              ipt_cfl(l)=cfloc
302   continue

       else

#include   "FastS/Compute/loop_ale_begin.for"
              ddi =sqrt( ti(lt+inci,1)*ti(lt+inci,1)
     &                  +ti(lt+inci,2)*ti(lt+inci,2) )
              ddi =ddi+ sqrt( ti(lt,1)*ti(lt,1)+ti(lt,2)*ti(lt,2) )

              ddj =sqrt( tj(lt+incj,1)*tj(lt+incj,1)
     &                  +tj(lt+incj,2)*tj(lt+incj,2) )
              ddj =ddj+ sqrt( tj(lt,1)*tj(lt,1)+tj(lt,2)*tj(lt,2) )

              surf= 1./max(ddi,ddj)

              u = rop(l,2)
              v = rop(l,3)
              r = rop(l,1) 
              c = sqrt(gam2*rop(l,5))
 
              ue=.5*(venti(lv      )+venti(lv+inci_ven      ))
              ve=.5*(venti(lv+v2ven)+venti(lv+inci_ven+v2ven))
 
              u  = u - ue
              v  = v - ve
              ur2 = u*u + v*v
              sp  = sqrt(ur2)
 
              xinvspc = 1./(sp+c)

              conv  = 2.*vol(lvo) * surf * xinvspc

              cfloc = 1.0/conv

              cfl(1) = max(cfl(1),cfloc)
              cfl(2) = min(cfl(2),cfloc)
              cfl(3) = cfl(3) + cfloc
              ipt_cfl(l)=cfloc
#include   "FastS/Compute/loop_end.for"

       endif

c******* Euler steady *****

      Elseif (param_int(LALE).eq.0.and.param_int(IFLOW).lt.2) Then

          !print*, "coucou"

        if(param_int(ITYPZONE).eq.0) then

#include   "FastS/Compute/loop_begin.for"
              ddi =sqrt( ti(lt+inci,1)*ti(lt+inci,1)
     &                  +ti(lt+inci,2)*ti(lt+inci,2)
     &                  +ti(lt+inci,3)*ti(lt+inci,3) )
              ddi =ddi+ sqrt( ti(lt,1)*ti(lt,1)+ti(lt,2)*ti(lt,2)
     &                       +ti(lt,3)*ti(lt,3) )

              ddj =sqrt( tj(lt+incj,1)*tj(lt+incj,1)
     &                  +tj(lt+incj,2)*tj(lt+incj,2)
     &                  +tj(lt+incj,3)*tj(lt+incj,3) )
              ddj =ddj+ sqrt( tj(lt,1)*tj(lt,1)+tj(lt,2)*tj(lt,2)
     &                       +tj(lt,3)*tj(lt,3) )

              ddk =sqrt( tk(lt+inck,1)*tk(lt+inck,1)
     &                  +tk(lt+inck,2)*tk(lt+inck,2)
     &                  +tk(lt+inck,3)*tk(lt+inck,3) )
              ddk =ddk+ sqrt( tk(lt,1)*tk(lt,1)+tk(lt,2)*tk(lt,2)
     &                       +tk(lt,3)*tk(lt,3) )

              surf= 1./max(ddi,ddj,ddk)

              u = rop(l,2)
              v = rop(l,3)
              w = rop(l,4)
              r = rop(l,1) 
              c = sqrt(gam2*rop(l,5))
 
              ur2 = u*u + v*v + w*w
              sp  = sqrt(ur2)
 
              xinvspc = 1./(sp+c)

              conv  = 2.*vol(lvo) * surf * xinvspc

              cfloc = 1.0/conv

              cfl(1) = max(cfl(1),cfloc)
              cfl(2) = min(cfl(2),cfloc)
              cfl(3) = cfl(3) + cfloc
              ipt_cfl(l)=cfloc
#include   "FastS/Compute/loop_end.for"

       elseif(param_int(ITYPZONE).eq.1) then
#include   "FastS/Compute/loop_begin.for"
 
              ddi =sqrt( ti(lt+inci,1)*ti(lt+inci,1)
     &                  +ti(lt+inci,2)*ti(lt+inci,2) )
              ddi =ddi+ sqrt( ti(lt,1)*ti(lt,1)+ti(lt,2)*ti(lt,2) )

              ddj =sqrt( tj(lt+incj,1)*tj(lt+incj,1)
     &                  +tj(lt+incj,2)*tj(lt+incj,2) )
              ddj =ddj+ sqrt( tj(lt,1)*tj(lt,1)+tj(lt,2)*tj(lt,2) )

              ddk =sqrt( tk(lt+inck,1)*tk(lt+inck,1) )
              ddk =ddk+ sqrt( tk(lt,1)*tk(lt,1))

              surf= 1./max(ddi,ddj,ddk)

              u = rop(l,2)
              v = rop(l,3)
              w = rop(l,4)
              r = rop(l,1) 
              c = sqrt(gam2*rop(l,5))
 
              ur2 = u*u + v*v + w*w
              sp  = sqrt(ur2)
 
              xinvspc = 1./(sp+c)

              conv  = 2.*vol(lvo) * surf * xinvspc

              cfloc = 1.0/conv

              cfl(1) = max(cfl(1),cfloc)
              cfl(2) = min(cfl(2),cfloc)
              cfl(3) = cfl(3) + cfloc
#include   "FastS/Compute/loop_end.for"

       elseif(param_int(ITYPZONE).eq.2) then

          ddi = 2.*sqrt( ti(1,1)*ti(1,1) )
          ddj = 2.*sqrt( tj(1,1)*tj(1,1) )
          ddk = 2.*sqrt( tk(1,1)*tk(1,1) )
          surf= 1./max(ddi,ddj,ddk)

          do 402 k = ind_loop(5), ind_loop(6)
          do 402 j = ind_loop(3), ind_loop(4)
             lij  = inddm( ind_loop(1) , j, k)
             lt   = 1
             lvo  = 1
             do 402 l = lij,  lij + ind_loop(2) - ind_loop(1)

              u = rop(l,2)
              v = rop(l,3)
              w = rop(l,4)
              r = rop(l,1) 
              c = sqrt(gam2*rop(l,5))
 
              ur2 = u*u + v*v + w*w
              sp  = sqrt(ur2)
 
              xinvspc = 1./(sp+c)

              conv  = 2.*vol(lvo) * surf * xinvspc

              cfloc = 1.0/conv

              cfl(1) = max(cfl(1),cfloc)
              cfl(2) = min(cfl(2),cfloc)
              cfl(3) = cfl(3) + cfloc
              ipt_cfl(l)=cfloc
402   continue

       else

#include   "FastS/Compute/loop_begin.for"
              ddi =sqrt( ti(lt+inci,1)*ti(lt+inci,1)
     &                  +ti(lt+inci,2)*ti(lt+inci,2) )
              ddi =ddi+ sqrt( ti(lt,1)*ti(lt,1)+ti(lt,2)*ti(lt,2) )

              ddj =sqrt( tj(lt+incj,1)*tj(lt+incj,1)
     &                  +tj(lt+incj,2)*tj(lt+incj,2) )
              ddj =ddj+ sqrt( tj(lt,1)*tj(lt,1)+tj(lt,2)*tj(lt,2) )

              surf= 1./max(ddi,ddj)

              u = rop(l,2)
              v = rop(l,3)
              r = rop(l,1) 
              c = sqrt(gam2*rop(l,5))
 
              ur2 = u*u + v*v
              sp  = sqrt(ur2)
 
              xinvspc = 1./(sp+c)

              conv  = 2.*vol(lvo) * surf * xinvspc

              cfloc = 1.0/conv

              cfl(1) = max(cfl(1),cfloc)
              cfl(2) = min(cfl(2),cfloc)
              cfl(3) = cfl(3) + cfloc
              ipt_cfl(l)=cfloc
#include   "FastS/Compute/loop_end.for"

       endif

      Endif

      end




