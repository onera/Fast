c***********************************************************************
c     $Date: 2013-02-04 19:27:20 +0100 (lun. 04 févr. 2013) $
c     $Revision: 38 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine spsource_SA_comp(ndom, nitcfg, param_int, param_real,
     &                            ind_loop, 
     &                            xmut,rop,coe, ti, tj, tk, vol,dlng, 
     &                            drodm)
c***********************************************************************
c_P                          O N E R A
c     ACT
c_A    Calcul de la contribution des termes sources 
c_A    pour l'equation de Spalart Allmaras 
c
c     VAL
c_V    Gaz parfait mono-espece
c_V    Navier-Stokes
c
c     INP
c_I    ndom      : numero du domaine calcule
c_I    dvardc(6) : gradients de nutild primitif aux centres
c_I                des cellules
c_I    dvardc(5) : gradients de ronutild conservatif aux centres
c_I                des cellules
c_I    vol       : volumes
c_I    rotn      : norme du rotationnel
c_I    dlng      : distance a la paroi
c
c     OUT
c
c     I/O
c_I    drodm    : terme source de l'equation de Spalart Allmaras 
c***********************************************************************
      implicit  none

#include "FastS/param_solver.h"

      INTEGER_E ndom,ind_loop(6), nitcfg, param_int(0:*)

      REAL_E xmut( param_int(NDIMDX) )
      REAL_E rop( param_int(NDIMDX) , param_int(NEQ) )
      REAL_E coe( param_int(NDIMDX) , param_int(NEQ_COE) )
      REAL_E drodm( param_int(NDIMDX), param_int(NEQ) )
      REAL_E ti( param_int(NDIMDX_MTR) , param_int(NEQ_IJ) ),
     &       tj( param_int(NDIMDX_MTR) , param_int(NEQ_IJ) ),
     &       tk( param_int(NDIMDX_MTR) , param_int(NEQ_K ) )
      REAL_E dlng(param_int(NDIMDX)),vol(param_int(NDIMDX_MTR))

      REAL_E param_real(0:*)


c Var loc 
      INTEGER_E incmax,i,j,k,l,icoe_pos,inci,incj,inck,inci_mtr,
     & incj_mtr,inck_mtr, l1,l2,l3,l4,l5,l6,ltij,lij,lt,ndimdx,lvo

      REAL_E adtild,adtild1,ad1,adelta1,adelta2,
     1  testzg2,testzg,fw,fv1,fv2,r1_1,amulam,anulam,
     &  amutild,anutild,voldes,tci,tcj,tck,sph2,amut,xmuprov,
     &  rotx,roty,rotz,rot,chi,prod,stild,stides,r,r2,g,fwg1,fwg,
     & aseuil,destruc,dest2,anvisc,tsource,tsourceb,
     & tsourcenu,f1,f2,dist,s,t,df2dchi,dstidnu,
     & dpdnu,dfwdg,dgdr,drdnu,dfwdnu,ddesdnu,
     & temp01,cmus1,coesut,t1,t1_1,
     & auijuij,adcut,variable2,fa,testfa,ra,c1,cw1,c1d,c2d,
     & tjx, tjy,tjz,tjx1,tjy1,tjz1,si,sj,sk,
     & tix, tiy,tiz,tix1,tiy1,tiz1,
     & tkx,tky,tkz,tkx1,tky1,tkz1,u1,u2,u3,u4,u5,u6,xvol,dudx,dudy,
     & dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,gam,rgp,cc

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

c.....formulation originelle
      fw(s)      = s*((1.+SA_CW3)/(s**6+SA_CW3))**(1./6.)
      fv1(s)     = (s**3)/(s**3+SA_CV1)
      fv2(s,t)   = 1.-s/(1.+s*t)

      cmus1  = param_real(VISCO+4)
      temp01 = 1./param_real(VISCO+3)
      coesut =  param_real(VISCO+2) * (1.+cmus1*temp01)
      c1     =  SA_CB2/SA_SIGMA
      cw1    = (SA_CB1/SA_CKARM/SA_CKARM)+(1.+SA_CB2)/SA_SIGMA

      gam    = param_real( GAMMA )
      rgp    = param_real( CVINF )*(gam-1.)

      incmax=param_int(NIJK)*param_int(NIJK+1)*param_int(NIJK+4)

      inci = 1
      incj = param_int(NIJK)
      inck = param_int(NIJK)*param_int(NIJK+1)

      inci_mtr = param_int(NIJK_MTR)
      incj_mtr = param_int(NIJK_MTR+1)
      inck_mtr = param_int(NIJK_MTR+2)

      ndimdx   = param_int(NDIMDX)

      IF(param_int(ITYPZONE).eq.0) THEN  !domaine 3d:

        If(param_int(ITYPCP).le.1.and.nitcfg.le.1) then !calcul implicite, on stocke coe pour ssor SA

          icoe_pos = 6
          if(icoe_pos.gt.param_int(NEQ_COE)) then
!$OMP SINGLE 
             write(*,*)'erreur dim coe: spsource'
!$OMP END SINGLE 
             stop
          endif

#include  "FastS/Compute/loop_begin.for"

#include       "FastS/Compute/SA/sourceSA_comp_grad_3d.for"
#include       "FastS/Compute/SA/sourceSA_prod_dest.for"
#include       "FastS/Compute/SA/sourceSA_LU.for"
               drodm(l,6)= drodm(l,6) + vol(lvo)*tsource

#include  "FastS/Compute/loop_end.for"

       Else  !calcul explicit, Stockage terme source inutile

#include  "FastS/Compute/loop_begin.for"

#include       "FastS/Compute/SA/sourceSA_comp_grad_3d.for"
#include       "FastS/Compute/SA/sourceSA_prod_dest.for"
               drodm(l,6)= drodm(l,6) + vol(lvo)*tsource

#include  "FastS/Compute/loop_end.for"


       endif!explicite/implicite

      !!!!!
      !!!!!
      !!!!!
      !!!!!
      !!!!!
      !3dhomo
      ELSEIF(param_int(ITYPZONE).eq.1)  THEN 


        If(param_int(ITYPCP).le.1.and.nitcfg.le.1) then !calcul implicite, on stocke coe pour ssor SA

          icoe_pos = 6
          if(icoe_pos.gt.param_int(NEQ_COE)) then
!$OMP SINGLE 
             write(*,*)'erreur dim coe: spsource'
!$OMP END SINGLE 
             stop
          endif


#include  "FastS/Compute/loop_begin.for"

#include       "FastS/Compute/SA/sourceSA_comp_grad_3dhomo.for"
#include       "FastS/Compute/SA/sourceSA_prod_dest.for"
#include       "FastS/Compute/SA/sourceSA_LU.for"
               drodm(l,6)= drodm(l,6) + vol(lvo)*tsource
#include  "FastS/Compute/loop_end.for"

       Else  !calcul explicit, Stockage terme source inutile

#include  "FastS/Compute/loop_begin.for"

#include       "FastS/Compute/SA/sourceSA_comp_grad_3dhomo.for"
#include       "FastS/Compute/SA/sourceSA_prod_dest.for"
               drodm(l,6)= drodm(l,6) + vol(lvo)*tsource
#include  "FastS/Compute/loop_end.for"


       endif!explicite/implicite 3dhomo


      !!!!!
      !!!!!
      !!!!!
      !!!!!
      !!!!!
      !3dcart
      ELSEIF(param_int(ITYPZONE).eq.2)  THEN 


      !metric
      lt  = indmtr(1 , 1, 1)
      lvo = lt
      tix = ti(lt,1)
      tjy = tj(lt,1)
      tkz = tk(lt,1)
      si      = abs (tix)
      sj      = abs (tjy)
      sk      = abs (tkz)
      sph2    = min(si,sj,sk)
      sph2    = vol(lvo)/sph2
      voldes  = vol(lvo) 
      xvol    = 0.5/vol(lvo)

        If(param_int(ITYPCP).le.1.and.nitcfg.le.1) then !calcul implicite, on stocke coe pour ssor SA

          icoe_pos = 6
          if(icoe_pos.gt.param_int(NEQ_COE)) then
!$OMP SINGLE 
             write(*,*)'erreur dim coe: spsource'
!$OMP END SINGLE 
             stop
          endif

          do k = ind_loop(5), ind_loop(6)
           do j = ind_loop(3), ind_loop(4)
                lij  =       inddm( ind_loop(1) , j, k)
!DEC$ IVDEP
             do l = lij, lij +  ind_loop(2)- ind_loop(1)

#include       "FastS/Compute/SA/sourceSA_comp_grad_3dcart.for"
#include       "FastS/Compute/SA/sourceSA_prod_dest.for"
#include       "FastS/Compute/SA/sourceSA_LU.for"
               drodm(l,6)= drodm(l,6) + vol(lvo)*tsource
             enddo
           enddo
          enddo

       Else  !calcul explicit, Stockage terme source inutile

          do k = ind_loop(5), ind_loop(6)
           do j = ind_loop(3), ind_loop(4)
                lij  =       inddm( ind_loop(1) , j, k)
!DEC$ IVDEP
             do l = lij, lij +  ind_loop(2)- ind_loop(1)

#include       "FastS/Compute/SA/sourceSA_comp_grad_3dcart.for"
#include       "FastS/Compute/SA/sourceSA_prod_dest.for"
               drodm(l,6)= drodm(l,6) + vol(lvo)*tsource
             enddo
           enddo
          enddo

       endif!explicite/implicite 3dcart

      !!!!!
      !!!!!
      !!!!!
      !!!!!
      !!!!!
      !2D     
      ELSE !2d: SA force

        If(param_int(ITYPCP).le.1.and.nitcfg.le.1) then !calcul implicite, on stocke coe pour ssor SA

          icoe_pos = 6
          if(icoe_pos.gt.param_int(NEQ_COE)) then
!$OMP SINGLE 
          write(*,*)'erreur dim coe: spsource'
!$OMP END SINGLE 
          stop
          endif

#include  "FastS/Compute/loop_begin.for"

#include       "FastS/Compute/SA/sourceSA_comp_grad_2d.for"
#include       "FastS/Compute/SA/sourceSA_prod_dest.for"
#include       "FastS/Compute/SA/sourceSA_LU.for"
               drodm(l,6)= drodm(l,6) + vol(lvo)*tsource
#include  "FastS/Compute/loop_end.for"

       Else  !calcul explicit, Stockage terme source inutile

#include  "FastS/Compute/loop_begin.for"

#include       "FastS/Compute/SA/sourceSA_comp_grad_2d.for"
#include       "FastS/Compute/SA/sourceSA_prod_dest.for"
               drodm(l,6)= drodm(l,6) + vol(lvo)*tsource
#include  "FastS/Compute/loop_end.for"


       endif!explicite/implicite




      ENDIF! typezone (2d,3d, cart,...0
      end