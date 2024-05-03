c***********************************************************************
c     $Date: 2021-07-30 09:20:20 
c     $Revision: 1 $
c     $Author: Mathieu Lugrin $
c***********************************************************************
      subroutine sciacovelli(ndom,neq_wig, param_int, param_real,
     &                     nitrun, ind_loop, 
     &                     rop, ti, tj, tk, vol,wig)
c***********************************************************************
c
c***********************************************************************
      implicit  none

#include "FastS/param_solver.h"

      INTEGER_E ndom,neq_wig,ind_loop(6), param_int(0:*), nitrun

      REAL_E rop( param_int(NDIMDX) , param_int(NEQ) )
      REAL_E ti( param_int(NDIMDX_MTR) , param_int(NEQ_IJ) ),
     &       tj( param_int(NDIMDX_MTR) , param_int(NEQ_IJ) ),
     &       tk( param_int(NDIMDX_MTR) , param_int(NEQ_K ) )

      REAL_E wig( param_int(NDIMDX)*neq_wig)
      REAL_E vol(param_int(NDIMDX_MTR))
      REAL_E param_real(0:*)
    
c Var loc 
      INTEGER_E incmax,i,j,k,l,inci,incj,inck,inci_mtr,
     & incj_mtr,inck_mtr, l1,l2,l3,l4,l5,l6,ltij,lij,lt,
     & ndimdx,lvo,v4,wigi,wigj,wigk,lw,shift_wig,
     & icorr,jcorr,kcorr,w_tg,inc_tg, ierr

      REAL_E tci,tcj,tck,sph2,tc_tg,
     &  rotx,roty,rotz,rot,
     & tjx,tjy,tjz,tjx1,tjy1,tjz1,si,sj,sk,
     & tix,tiy,tiz,tix1,tiy1,tiz1,
     & tkx,tky,tkz,tkx1,tky1,tkz1,u1,u2,u3,u4,u5,u6,xvol,dudx,dudy,dudz,
     & dx,dy,dz,tn,tx,ty,tz,nx,ny,nz,
     & div,umag,rgp,gam,mach,coef1,coef2,
     & p0,p1,p2,p3,p4,p5,p6

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      if(ind_loop(1).gt.ind_loop(2)) return 
      if(ind_loop(3).gt.ind_loop(4)) return 
      if(ind_loop(5).gt.ind_loop(6)) return


      gam =param_real(GAMMA)
      rgp    = param_real( CVINF  )*(gam-1.)
      incmax=param_int(NIJK)*param_int(NIJK+1)*param_int(NIJK+4)


      inci = 1
      incj = param_int(NIJK)
      inck = param_int(NIJK)*param_int(NIJK+1)

      inci_mtr = param_int(NIJK_MTR)
      incj_mtr = param_int(NIJK_MTR+1)
      inck_mtr = param_int(NIJK_MTR+2)

      ndimdx   = param_int(NDIMDX)

      shift_wig=0
      if(param_int(KFLUDOM).eq.2) shift_wig = 3*param_int(NDIMDX)
      wigi = shift_wig
      wigj = shift_wig +   param_int(NDIMDX)
      wigk = shift_wig + 2*param_int(NDIMDX)

      icorr=0
      jcorr=0
      kcorr=0
      if(ind_loop(2).eq.param_int(IJKV))   icorr=1
      if(ind_loop(4).eq.param_int(IJKV+1)) jcorr=1
      if(ind_loop(6).eq.param_int(IJKV+2)) kcorr=1


      IF(param_int(ITYPZONE).eq.0)  THEN 
        DO k = ind_loop(5), ind_loop(6)
         DO j = ind_loop(3), ind_loop(4)

#include "FastS/Compute/loopI_begin.for"                  
#include     "FastS/Compute/SA/metric_3dfull.for"
#include     "FastS/Compute/SA/div_rot_3dfull.for" 
              w_tg   = wigi
              tc_tg  = tci
              inc_tg = inci
#include     "FastS/Compute/SENSOR/shock_sensor.for" 
              w_tg   = wigj
              tc_tg  = tcj
              inc_tg = incj
#include     "FastS/Compute/SENSOR/shock_sensor.for" 
              w_tg   = wigk
              tc_tg  = tck
              inc_tg = inck
#include     "FastS/Compute/SENSOR/shock_sensor.for" 
          enddo    
          if(icorr.eq.1) then !flux manquant en I
             i   = ind_loop(2) + 1
             l   = inddm(  i, j, k)
#include     "FastS/Compute/SA/metric_3dfull.for"
#include     "FastS/Compute/SA/div_rot_3dfull.for" 
              w_tg   = wigi
              tc_tg  = tci
              inc_tg = inci
#include     "FastS/Compute/SENSOR/shock_sensor.for" 
          endif !
         ENDDO
         !complement jmax
         If(jcorr.eq.1) then
           j    = ind_loop(4)+1
#include   "FastS/Compute/loopI_begin.for"
#include     "FastS/Compute/SA/metric_3dfull.for"
#include     "FastS/Compute/SA/div_rot_3dfull.for" 
              w_tg   = wigj
              tc_tg  = tcj
              inc_tg = incj
#include     "FastS/Compute/SENSOR/shock_sensor.for" 
            enddo
         Endif
        ENDDO
        !complement kmax
        If(kcorr.eq.1) then
          k    = ind_loop(6)+1               
          do j = ind_loop(3),ind_loop(4)     
#include    "FastS/Compute/loopI_begin.for"
#include     "FastS/Compute/SA/metric_3dfull.for"
#include     "FastS/Compute/SA/div_rot_3dfull.for" 
              w_tg   = wigk
              tc_tg  = tck
              inc_tg = inck
#include     "FastS/Compute/SENSOR/shock_sensor.for"  
            enddo                    
          enddo                      
        Endif


      ELSEIF(param_int(ITYPZONE).eq.1)  THEN

        DO k = ind_loop(5), ind_loop(6)
         DO j = ind_loop(3), ind_loop(4)

#include "FastS/Compute/loopI_begin.for"                  
#include     "FastS/Compute/SA/metric_3dhomo.for"
#include     "FastS/Compute/SA/div_rot_3dhomo.for" 
              w_tg   = wigi
              tc_tg  = tci
              inc_tg = inci
#include     "FastS/Compute/SENSOR/shock_sensor.for" 
              w_tg   = wigj
              tc_tg  = tcj
              inc_tg = incj
#include     "FastS/Compute/SENSOR/shock_sensor.for" 
              w_tg   = wigk
              tc_tg  = tck
              inc_tg = inck
#include     "FastS/Compute/SENSOR/shock_sensor.for" 
          enddo    
          if(icorr.eq.1) then !flux manquant en I
             i   = ind_loop(2) + 1
             l   = inddm(  i, j, k)
#include     "FastS/Compute/SA/metric_3dhomo.for"
#include     "FastS/Compute/SA/div_rot_3dhomo.for" 
              w_tg   = wigi
              tc_tg  = tci
              inc_tg = inci
#include     "FastS/Compute/SENSOR/shock_sensor.for" 
          endif !
         ENDDO
         !complement jmax
         If(jcorr.eq.1) then
           j    = ind_loop(4)+1
#include   "FastS/Compute/loopI_begin.for"
#include     "FastS/Compute/SA/metric_3dhomo.for"
#include     "FastS/Compute/SA/div_rot_3dhomo.for" 
              w_tg   = wigj
              tc_tg  = tcj
              inc_tg = incj
#include     "FastS/Compute/SENSOR/shock_sensor.for" 
            enddo
         Endif
        ENDDO
        !complement kmax
        If(kcorr.eq.1) then
          k    = ind_loop(6)+1               
          do j = ind_loop(3),ind_loop(4)     
#include    "FastS/Compute/loopI_begin.for"
#include     "FastS/Compute/SA/metric_3dhomo.for"
#include     "FastS/Compute/SA/div_rot_3dhomo.for" 
              w_tg   = wigk
              tc_tg  = tck
              inc_tg = inck
#include     "FastS/Compute/SENSOR/shock_sensor.for"  
            enddo                    
          enddo                      
        Endif


      ELSEIF(param_int(ITYPZONE).eq.2)  THEN
      !metric
      lt  = indmtr(1 , 1, 1)
      lvo = lt
      tix = ti(lt,1)
      tjy = tj(lt,1)
      tkz = tk(lt,1)
      tci     = abs (tix)
      tcj     = abs (tjy)
      tck     = abs (tkz)
      dx      = vol(lvo)/tci
      dy      = vol(lvo)/tcj
      dz      = vol(lvo)/tck
      xvol    = 0.5/vol(lvo)

      DO k = ind_loop(5), ind_loop(6)
        DO j = ind_loop(3), ind_loop(4)

#include "FastS/Compute/loopI3dcart_begin.for"
#include     "FastS/Compute/SA/div_rot_3dcart.for" 
              w_tg   = wigi
              tc_tg  = tci
              inc_tg = inci
#include     "FastS/Compute/SENSOR/shock_sensor.for" 
              w_tg   = wigj
              tc_tg  = tcj
              inc_tg = incj
#include     "FastS/Compute/SENSOR/shock_sensor.for" 
              w_tg   = wigk
              tc_tg  = tck
              inc_tg = inck
#include     "FastS/Compute/SENSOR/shock_sensor.for" 
          enddo    
          if(icorr.eq.1) then !flux manquant en I
             i   = ind_loop(2) + 1
             l   = inddm(  i, j, k)
#include     "FastS/Compute/SA/div_rot_3dcart.for" 
              w_tg   = wigi
              tc_tg  = tci
              inc_tg = inci
#include     "FastS/Compute/SENSOR/shock_sensor.for" 
          endif !
         ENDDO
         !complement jmax
         If(jcorr.eq.1) then
           j    = ind_loop(4)+1
#include   "FastS/Compute/loopI3dcart_begin.for"
#include     "FastS/Compute/SA/div_rot_3dcart.for" 
              w_tg   = wigj
              tc_tg  = tcj
              inc_tg = incj
#include     "FastS/Compute/SENSOR/shock_sensor.for" 
            enddo
         Endif
        ENDDO
        !complement kmax
        If(kcorr.eq.1) then
          k    = ind_loop(6)+1               
          do j = ind_loop(3),ind_loop(4)     
#include    "FastS/Compute/loopI3dcart_begin.for"
#include     "FastS/Compute/SA/div_rot_3dcart.for" 
              w_tg   = wigk
              tc_tg  = tck
              inc_tg = inck
#include     "FastS/Compute/SENSOR/shock_sensor.for"  
            enddo                    
          enddo                      
        Endif

      ELSE

        DO k = ind_loop(5), ind_loop(6)
         DO j = ind_loop(3), ind_loop(4)

#include "FastS/Compute/loopI_begin.for"                  
#include     "FastS/Compute/SA/metric_2d.for"
#include     "FastS/Compute/SA/div_rot_2d.for" 
              w_tg   = wigi
              tc_tg  = tci
              inc_tg = inci
#include     "FastS/Compute/SENSOR/shock_sensor.for" 
              w_tg   = wigj
              tc_tg  = tcj
              inc_tg = incj
#include     "FastS/Compute/SENSOR/shock_sensor.for" 
          enddo    
          if(icorr.eq.1) then !flux manquant en I
             i   = ind_loop(2) + 1
             l   = inddm(  i, j, k)
#include     "FastS/Compute/SA/metric_2d.for"
#include     "FastS/Compute/SA/div_rot_2d.for" 
              w_tg   = wigi
              tc_tg  = tci
              inc_tg = inci
#include     "FastS/Compute/SENSOR/shock_sensor.for" 
          endif !
         ENDDO
         !complement jmax
         If(jcorr.eq.1) then
           j    = ind_loop(4)+1
#include   "FastS/Compute/loopI_begin.for"
#include     "FastS/Compute/SA/metric_2d.for"
#include     "FastS/Compute/SA/div_rot_2d.for" 
              w_tg   = wigj
              tc_tg  = tcj
              inc_tg = incj
#include     "FastS/Compute/SENSOR/shock_sensor.for" 
            enddo
         Endif
        ENDDO

      ENDIF
      end
