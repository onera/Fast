c***********************************************************************
c     $Date: 2021-07-30 09:20:20 
c     $Revision: 1 $
c     $Author: Mathieu Lugrin $
c***********************************************************************
      subroutine ducros(ndom, param_int, param_real,
     &                     ind_loop, 
     &                     rop, ti, tj, tk, vol,wig)
c***********************************************************************
c
c***********************************************************************
      implicit  none

#include "FastS/param_solver.h"

      INTEGER_E ndom,ind_loop(6), param_int(0:*)

      REAL_E rop( param_int(NDIMDX) , param_int(NEQ) )
      REAL_E ti( param_int(NDIMDX_MTR) , param_int(NEQ_IJ) ),
     &       tj( param_int(NDIMDX_MTR) , param_int(NEQ_IJ) ),
     &       tk( param_int(NDIMDX_MTR) , param_int(NEQ_K ) )

      REAL_E wig( param_int(NDIMDX)*4),vol(param_int(NDIMDX_MTR))
      REAL_E param_real(0:*)
    
c Var loc 
      INTEGER_E incmax,i,j,k,l,inci,incj,inck,inci_mtr,
     & incj_mtr,inck_mtr, l1,l2,l3,l4,l5,l6,ltij,lij,lt,
     & ndimdx,lvo,v4,wigd

      REAL_E adtild,adtild1,ad1,adelta1,adelta2,
     &  testzg2,testzg,fw,fv1,fv2,fvv1,amulam,anulam,
     &  amutild,anutild,voldes,tci,tcj,tck,sph2,amut,xmuprov,
     &  rotx,roty,rotz,rot,chi,prod,stild,stides,r,r2,g,fwg1,fwg,
     & aseuil,destruc,dest2,anvisc,tsource,tsourceb,
     & tsourcenu,f1,f2,dist,s,t,df2dchi,dstidnu,
     & dpdnu,dfwdg,dgdr,drdnu,dfwdnu,ddesdnu,
     & temp01,cmus1,t1,t1_1,
     & auijuij,adcut,variable2,fa,testfa,ra,c1,cw1,
     & tjx,tjy,tjz,tjx1,tjy1,tjz1,si,sj,sk,
     & tix,tiy,tiz,tix1,tiy1,tiz1,
     & tkx,tky,tkz,tkx1,tky1,tkz1,u1,u2,u3,u4,u5,u6,xvol,dudx,dudy,dudz,
     & dx,dy,dz,tn,tx,ty,tz,nx,ny,nz,
     & divx,divy,divz,div,umag,rgp,gam,mach,coef1,coef2

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"
      gam = param_real(GAMMA)
      rgp = param_real( CVINF  )*(gam-1.)
      incmax=param_int(NIJK)*param_int(NIJK+1)*param_int(NIJK+4)

      inci = 1
      incj = param_int(NIJK)
      inck = param_int(NIJK)*param_int(NIJK+1)

      inci_mtr = param_int(NIJK_MTR)
      incj_mtr = param_int(NIJK_MTR+1)
      inck_mtr = param_int(NIJK_MTR+2)

      ndimdx   = param_int(NDIMDX)

      wigd = 3*param_int(NDIMDX)

      coef1 = param_real(HYPER_COEF1)
      coef2 = param_real(HYPER_COEF2)

      IF(param_int(ITYPZONE).eq.0)  THEN 
#include   "FastS/Compute/loop_begin.for"
#include       "FastS/Compute/SA/metric_3dfull.for"
#include       "FastS/Compute/SENSORHYPER/3dfull/rot_3dfull.for" 
#include   "FastS/Compute/loop_end.for"
      ELSEIF(param_int(ITYPZONE).eq.1)  THEN
#include   "FastS/Compute/loop_begin.for"
#include     "FastS/Compute/SA/metric_3dhomo.for" 
#include     "FastS/Compute/SENSORHYPER/3dhomo/rot_3dhomo.for" 
#include   "FastS/Compute/loop_end.for"
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
      sph2    = max(dx,dy,dz)
      voldes  = vol(lvo) 
      xvol    = 0.5/vol(lvo)
      do k = ind_loop(5), ind_loop(6)
       do j = ind_loop(3), ind_loop(4)
#include   "FastS/Compute/loopI3dcart_begin.for"
#include     "FastS/Compute/SENSORHYPER/3dcart/rot_3dcart.for" 
#include   "FastS/Compute/loop_end.for"

      ELSE

#include   "FastS/Compute/loop_begin.for"
#include       "FastS/Compute/SA/metric_2d.for" 
#include       "FastS/Compute/SENSORHYPER/2d/rot_2d.for" 
#include   "FastS/Compute/loop_end.for"
      ENDIF
      end
