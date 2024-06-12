!C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!c Partir de ce fichier pour les débits
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!c***********************************************************************
!c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
!c     $Revision: 64 $
!c     $Author: MehmetDemirci $
!c***********************************************************************
      subroutine dpJ_dpW_fluausm_euler_o3_2d(ndom, ithread, 
     &                             param_int, param_real, param_int_eff, 
     &                             ind_loop, pos, 
     &                             rop, wig, 
     &                             x, y, z,
     &                             dpCdp_dpW, dpClp_dpW,
     &                             venti, ventj, ventk, 
     &                             ti, tj, tk, vol, xmut,
     &                             cosAoA, sinAoA, surfinv)
!c*****************************************************************************
!c_U   USER : MEHMETDEMIRCI
!c
!c     ACT
!c_A    Calcul du membre de droite de l equation adjointe
!c
!c     VAL
!c_V    gaz parfait monoespece
!c_V    processeur domaine
!c_V    steady/unsteady
!c
!c     INP
!c_I    tijk     : vecteur normale aux facettes des mailles
!c_I    ventijk     : vitesses d entrainement aux facettes preced.
!c_I    qm,qp    : etats droit et gauche aux interfaces d une maille
!c
!c     LOC
!c_L    flu      : flux convectifs dans une direction de maillage
!c
!c     I/O
!c_/    dpClp_dpW,dpCdp_dpW : Membres de droite des equations adjointes
!c                                respectives pour Clp et Cdp en tenant
!c                                compte des termes dpPg_dpP
!c
!c      Hypothese : ventijk est ici suppose independant de rop
!c
!c
!c     !!! Les constantes ci_mjr, cj_mjr et ck_mjr donnees ici sont des
!c         valeurs imprimees sur l ecran dans l unique but de valider cette
!c         routine. Il reste a les extraire en utilisant shape_tab_mtr, comme
!c         dans bvbs_wall_inviscid.for 
!c*****************************************************************************
      implicit none

      real souszero
      parameter(souszero=-1e-12)

#include "FastS/param_solver.h"

!**** VARIABLES LOCALES A FAIRE PASSER DANS PARAM_REAL ****
      REAL*8 cosAoA, sinAoA, surfinv
      !real*8 AoA
      !parameter(AoA=1.25*3.1415926535897932384/180.)
      !real*8 surfinv
      !parameter(surfinv = 100.)


!# 1 "./FastS/param_solver.h" 1


!# 44 "build/x86_r8/FastS/Compute/AUSM/2d/eff_fluausm_euler_o3_2d.for" 2


      INTEGER*4 ndom, ithread, param_int(0:*), 
     & param_int_eff(0:*), iter
      INTEGER_E ind_loop(6)

      REAL*8   xmut( param_int(41) )
      REAL*8   rop( param_int(41) * param_int(36) )
      REAL*8   wig( param_int(41) * 3 )
      REAL*8   venti( param_int(44)* param_int(40))
      REAL*8   ventj( param_int(44)* param_int(40))
      REAL*8   ventk( param_int(44)* param_int(40))

      REAL*8  ti( param_int(43) * param_int(37) ), 
     &              tj( param_int(43) * param_int(37) ), 
     &              tk( param_int(43) * param_int(38 ) )
      REAL*8   vol( param_int(43) )
      REAL*8   x( param_int(42) )
      REAL*8   y( param_int(42) )
      REAL*8   z( param_int(42) )

      REAL*8 param_real(0:*), pos(3)


      REAL*8  dpClp_dpW(param_int(41) * param_int(36))
      REAL*8  dpCdp_dpW(param_int(41) * param_int(36))


!C Var loc
      INTEGER_E inc,incmax,l,lt,i,j,k,incmax2,nm,nm2,np,                  
     & l0,lt0,inci,incj,inck,ci,cj,lij,ltij,inci_mtr, incj_mtr,     
     & inck_mtr,icorr,jcorr,ls,v1,v2,v3,v4,v5,v6,wig_i,wig_j,wig_k, 
     & lfij,lxij,lf,lx,inc_x1,inc_x2,inc_x3,vslp,lvol,lvor,ir,il,   
     & lt200,lt100,lt010,lt210,lt020,lt110,lt002,lt012,lt102,lt001, 
     & lt021,lt201,lt120,v2flu,indflu,lvo,lvo200,lvo020,lvo002,     
     & l200,l100,l010,l020,l110,l101,l011,v1mtr,v2mtr,v3mtr,              
     & l001,l002,l210,l220,l201,l202,l021,l022,l120,l102,l012

      REAL*8 c1,c2,c3,c4,c5,c6,c4sa,c5sa,c6sa,si,sj,sk,qm,qp,                 
     & tcx,tcy,tcz,tc,r1,h1,rou1,rov1,row1,r2,h2,rou2,rov2,row2,              
     & gam,gam1,qn1,qn2,u,tdu,p1p2,roref,uref,tam,tam1,son,c,gam2,      
     & qm1,qm2,qm3,qm4,qm5,qm6,qp1,qp2,qp3,qp4,qp5,qp6,mut1,mut2,       
     & flu1,flu2,flu3,flu4,flu5,flu6,p1,p2,qen,sigma_1,ck_vent,               
     & div,f1,f2,f3,f4,f5,f6,fv,fv5,volinv,test,cmus1,temp01,coesut,    
     & tix,tiy,tiz,tix1,tiy1,tiz1,tjx,tjy,tjz,tjx1,tjy1,tjz1,tkx,       
     & tky,tkz,tkx1,tky1,tkz1,xmutvol,prandt,gam3,cvisq,rgp,                  
     & gradU_nx,gradU_ny,gradU_nz, gradV_nx,gradV_ny,gradV_nz,                
     & gradW_nx,gradW_ny,gradW_nz, gradT_nx,gradT_ny,gradT_nz,                
     & delp,delm,delq,slq,slp,roff,tmin_1,du,dv,dw,dp,dqn,s_1,nx,ny,nz, 
     & qn,r,v,w,h,q,r_1,psiroe,sens,sens1,flagi,flagj,flagk,norm


!    VARIABLES LOCALES POUR LE CALCUL DES DERIVEES *************************
      REAL*8 dpflu2_dpu, dpflu2_dptdu, dpflu2_dpp1p2,
     & dpflu2_dpqn1, dpflu2_dpqn2, dpflu2_dpp1, dpflu2_dpp2,
     & dpflu2_dprou1, dpflu2_dprou2,
     & dpflu2_dpqp1, dpflu2_dpqp2, dpflu2_dpqp3, dpflu2_dpqp5,
     & dpflu2_dpqm1, dpflu2_dpqm2, dpflu2_dpqm3, dpflu2_dpqm5

      REAL*8 dpflu3_dpu, dpflu3_dptdu, dpflu3_dpp1p2,
     & dpflu3_dpqn1, dpflu3_dpqn2, dpflu3_dpp1, dpflu3_dpp2,
     & dpflu3_dprov1, dpflu3_dprov2, 
     & dpflu3_dpqp1, dpflu3_dpqp2, dpflu3_dpqp3, dpflu3_dpqp5,
     & dpflu3_dpqm1, dpflu3_dpqm2, dpflu3_dpqm3, dpflu3_dpqm5

      REAL*8 dpflu2_dprop1_nm2, dpflu2_dprop1_nm, dpflu2_dprop1_l,
     & dpflu2_dprop1_np, dpflu2_dprop2_nm2, dpflu2_dprop2_nm,
     & dpflu2_dprop2_l, dpflu2_dprop2_np, dpflu2_dprop3_nm2,
     & dpflu2_dprop3_nm, dpflu2_dprop3_l, dpflu2_dprop3_np,
     & dpflu2_dprop5_nm2, dpflu2_dprop5_nm, dpflu2_dprop5_l,
     & dpflu2_dprop5_np

      REAL*8 dpflu3_dprop1_nm2, dpflu3_dprop1_nm, dpflu3_dprop1_l,
     & dpflu3_dprop1_np, dpflu3_dprop2_nm2, dpflu3_dprop2_nm,
     & dpflu3_dprop2_l, dpflu3_dprop2_np, dpflu3_dprop3_nm2,
     & dpflu3_dprop3_nm, dpflu3_dprop3_l, dpflu3_dprop3_np,
     & dpflu3_dprop5_nm2, dpflu3_dprop5_nm, dpflu3_dprop5_l,
     & dpflu3_dprop5_np

      REAL*8 dptdu_dpu, dpu_dpqn1, dpu_dptam1, dptam1_dptam,
     & dptam_dpson, dpson_dpqn1, dpu_dpqn2, dpp1p2_dpp1, dpp1p2_dpp2,
     & dpu_dpp1, dpu_dpp2

      REAL*8 dpqn1_dpqp2, dprou1_dpqp2, dpqn1_dpqp3, dprov1_dpqp3,
     & dpqn2_dpqm2, dprou2_dpqm2, dpqn2_dpqm3, dprov2_dpqm3

      REAL*8 dpp1_dpqp1, dprou1_dpqp1, dpp1_dpqp5, dpp2_dpqm1,
     & dprou2_dpqm1, dpp2_dpqm5, dprov1_dpqp1, dprov2_dpqm1

      REAL*8 dpqp_dprop_l, dpqp_dprop_nm, dpqp_dprop_nm2, dpqm_dprop_l,
     & dpqm_dprop_nm, dpqm_dprop_np

      REAL*8 dpson_dpc, dpflu2_dpc, dpflu3_dpc, dpc_dprop_l5

      REAL*8 cs,cc

      REAL*8 dpCdp_dprop1_l, dpCdp_dprop2_l,dpCdp_dprop3_l,
     & dpCdp_dprop5_l, dpCdp_dprop1_np, dpCdp_dprop2_np,
     & dpCdp_dprop3_np, dpCdp_dprop5_np, dpClp_dprop1_l,dpClp_dprop2_l,
     & dpClp_dprop3_l, dpClp_dprop5_l, dpClp_dprop1_np,dpClp_dprop2_np,
     & dpClp_dprop3_np, dpClp_dprop5_np

      REAL*8 dpCdp_dprop1_nm, dpCdp_dprop2_nm,dpCdp_dprop3_nm,
     & dpCdp_dprop5_nm, dpCdp_dprop1_nm2, dpCdp_dprop2_nm2,
     & dpCdp_dprop3_nm2, dpCdp_dprop5_nm2, dpClp_dprop1_nm,
     & dpClp_dprop2_nm, dpClp_dprop3_nm, dpClp_dprop5_nm,
     & dpClp_dprop1_nm2,dpClp_dprop2_nm2,
     & dpClp_dprop3_nm2, dpClp_dprop5_nm2

!    VARIABLES LOCALES POUR LE CALCUL DE dpPg_dpP **************************

      REAL*8 dpRop1_b_dpRop1_h, dpRop1_b_dpRop2_h,
     & dpRop1_b_dpRop3_h, dpRop1_b_dpRop4_h,
     & dpRop1_b_dpRop5_h

      REAL*8 dpRop2_b_dpRop1_h, dpRop2_b_dpRop2_h,
     & dpRop2_b_dpRop3_h, dpRop2_b_dpRop4_h,
     & dpRop2_b_dpRop5_h

      REAL*8 dpRop3_b_dpRop1_h, dpRop3_b_dpRop2_h,
     & dpRop3_b_dpRop3_h, dpRop3_b_dpRop4_h,
     & dpRop3_b_dpRop5_h

      REAL*8 dpRop4_b_dpRop1_h, dpRop4_b_dpRop2_h,
     & dpRop4_b_dpRop3_h, dpRop4_b_dpRop4_h,
     & dpRop4_b_dpRop_h

      REAL*8 dpRop5_b_dpRop1_h, dpRop5_b_dpRop2_h,
     & dpRop5_b_dpRop3_h, dpRop5_b_dpRop4_h,
     & dpRop5_b_dpRop5_h

      REAL*8 ncx,ncy,ncz,surf_nc

      INTEGER  neq_mtr, ic,jc,kc,kc_vent, idir
      REAL*8 ci_mtr,cj_mtr,ck_mtr,c_ale

!    VARIABLES LOCALES POUR LE CALCUL DE dpP_dpW **************************

       REAL*8 dpRop1_dpW1_l, dpRop1_dpW1_np,
     & dpRop2_dpW1_l, dpRop2_dpW1_np,
     & dpRop2_dpW2_l, dpRop2_dpW2_np,
     & dpRop3_dpW1_l, dpRop3_dpW1_np,
     & dpRop3_dpW3_l, dpRop3_dpW3_np,
     & dpRop4_dpW1_l, dpRop4_dpW1_np,
     & dpRop4_dpW4_l, dpRop4_dpW4_np,
     & dpRop5_dpW1_l, dpRop5_dpW2_l, dpRop5_dpW3_l,
     & dpRop5_dpW4_l, dpRop5_dpW5_l,
     & dpRop5_dpW1_np, dpRop5_dpW2_np, dpRop5_dpW3_np,
     & dpRop5_dpW4_np, dpRop5_dpW5_np

       REAL*8 dpRop1_dpW1_nm, dpRop1_dpW1_nm2,
     & dpRop2_dpW1_nm, dpRop2_dpW1_nm2,
     & dpRop2_dpW2_nm, dpRop2_dpW2_nm2,
     & dpRop3_dpW1_nm, dpRop3_dpW1_nm2,
     & dpRop3_dpW3_nm, dpRop3_dpW3_nm2,
     & dpRop4_dpW1_nm, dpRop4_dpW1_nm2,
     & dpRop4_dpW4_nm, dpRop4_dpW4_nm2,
     & dpRop5_dpW1_nm, dpRop5_dpW2_nm, dpRop5_dpW3_nm,
     & dpRop5_dpW4_nm, dpRop5_dpW5_nm,
     & dpRop5_dpW1_nm2, dpRop5_dpW2_nm2, dpRop5_dpW3_nm2,
     & dpRop5_dpW4_nm2, dpRop5_dpW5_nm2

       REAL*8 Cpm, Cvm, rho_e_cin_np, rho_e_cin_l,
     & W1_l, W2_l, W3_l, W4_l, W5_l,
     & W1_np, W2_np, W3_np, W4_np, W5_np,
     & rho_e_cin_nm2, rho_e_cin_nm,
     & W1_nm, W2_nm, W3_nm, W4_nm, W5_nm,
     & W1_nm2, W2_nm2, W3_nm2, W4_nm2, W5_nm2


c -----------------*------------------*-------------------*-----------------*

!C    adresse point courant pour tableau de la taille d'un domaine 
      INTEGER_E inddm, i_1,j_1,k_1

      inddm(i_1,j_1,k_1) = 1  
     &     + (i_1+param_int(0+3)-1) 
     &     + (j_1+param_int(0+3)-1)*param_int(0) 
     &     + (k_1+param_int(0+4)-1)*param_int(0)*param_int(0+1)


!C    adresse point courant pour tableau x ou de la taille d'un domaine 
      INTEGER_E indcg, i_2,j_2,k_2

      indcg(i_2,j_2,k_2) = 1  
     &     + (i_2+param_int(10+3)-1)  
     &     + (j_2+param_int(10+3)-1)*param_int(10)  
     &     + (k_2+param_int(10+4)-1)  
     *       *param_int(10)*param_int(10+1)


!C    adresse interface pour tableau metric
      INTEGER_E indmtr, i_3,j_3,k_3

      indmtr(i_3,j_3,k_3) =  1
     &                    + (i_3+param_int(5+3)-1)*param_int(5)
     &                    + (j_3+param_int(5+3)-1)*param_int(5+1)
     &                    + (k_3+param_int(5+4)-1)*param_int(5+2)

      indflu(i_1,j_1,k_1) =  1
     &             + (i_1-param_int_eff(29  ))
     &             + (j_1-param_int_eff(29+1))*param_int_eff(29+3)
     &             + (k_1-param_int_eff(29+2))*param_int_eff(29+4)

      !limiteur 'minmod'

!CC!DIR$ ASSUME_ALIGNED xmut: 16

      if (ithread.eq.1) then 
        write(*,*) ""
        write(*,*) ""
        write(*,*) " dpEff_dpW_fluausm_euler_o3_2d  " 
        write(*,*) "-----------------------"
      end if

      if(ind_loop(1).gt.ind_loop(2)) return 
      if(ind_loop(3).gt.ind_loop(4)) return 
      if(ind_loop(5).gt.ind_loop(6)) return

      inci = 1
      incj = param_int(0)
      inck = param_int(0)*param_int(0+1)
      inci_mtr = param_int(5)
      incj_mtr = param_int(5+1)
      inck_mtr = param_int(5+2)

      !metric
      lt  = indmtr(1 , 1, 1)
      lvo = lt
      tcx = ti(lt)
      tcy = tj(lt)
      si      = abs (tcx)
      sj      = abs (tcy)
      volinv  = 0.5/vol(lvo)

      !-----Variables physiques
      gam    = param_real( 1 )
      rgp    = param_real( 2 )*(gam-1.)  !Cv(gama-1)= R (gas parfait)
      prandt = param_real( 10 )
      gam1   = gam/(gam-1.)
      gam2   = 1./gam
      gam3    = gam1/ prandt*rgp

      cmus1  =    param_real( 13 )
      temp01 = 1./param_real( 12)
      coesut =    param_real( 11) * (1.+cmus1*temp01)
      sigma_1 =1./(2./3.)

      roref= param_real( 3)
      uref = param_real( 5 )

      psiroe= param_real( 34 )
      tmin_1= 100./param_real( 6 )  !!si T< 0.01Tinf, alors limiteur null

      c1     = 0.02*uref                 ! modif suite chant metrique et suppression tc dans flux final
      c2     = 0.02/(uref*roref)   ! modif suite chant metrique et suppression tc dans flux final
      c3     = -2.

      !    roff MUSCL
      c6     = 1./6.
      c4     = 5.*c6
      c5     = 2.*c6
      c6     =-1.*c6

!c      c7     = c4/c5

      cvisq = 1./3

      v1 = 0
      v2 =   param_int(41)
      v3 = 2*param_int(41)
      v4 = 3*param_int(41)
      v5 = 4*param_int(41)
      v6 = 5*param_int(41)

      v1mtr =   0
      v2mtr =   param_int(43)
      v3mtr = 2*param_int(43)

      wig_i = v1
      wig_j = v2
      wig_k = v3


      v2flu =   param_int_eff(41)

      norm =2./(param_real(3)*param_real(5)**2)

      sens = norm
      sens1=1.
      if(mod(param_int_eff(35),2).eq.1) then
         sens =-sens
         sens1=-sens1
      endif

      
      if(param_int(25).eq.0) then
              flagi = 1.
              flagj = 1.
              flagk = 1.
      elseif(param_int(25).eq.1) then
              flagi = 1.
              flagj = 1.
              flagk = 0.
      elseif(param_int(25).eq.2) then
              flagi = 0.
              flagj = 0.
              flagk = 0.
      else
              flagi = 1.
              flagj = 1.
              flagk = 0.
      endif
                

!CC!DIR$ ASSUME (mod(inck,   4) .eq. 0)
!CC!DIR$ ASSUME (mod(incj,   4) .eq. 0)
!CC!DIR$ ASSUME (mod(param_int(41), 4) .eq. 0)


      IF(param_int_eff(35).eq.1) THEN

!************************************************************************************
!*********************** CALCUL DES EFFORTS SUR LA FACE I ***************************
!************************************************************************************

              inc_x1 = param_int(10)                                   !(i,j+1,k  )
              inc_x2 = param_int(10)*param_int(10+1)       !(i,j  ,k+1)
              inc_x3 = inc_x1 + inc_x2
              iter = 0

              DO k = ind_loop(5), ind_loop(6)
                 DO j = ind_loop(3), ind_loop(4)

	      lij  =       inddm( ind_loop(1) , j, k) -1
	      ltij = lij - indmtr(ind_loop(1) , j, k) +1
	      lfij = lij - indflu(ind_loop(1) , j, k) +1
	      lxij = lij - indcg( ind_loop(1) , j, k) +1
              iter = iter + 1

!CC    !DIR$ ASSUME (mod(lij,4) .eq. 0)
!DIR$ IVDEP

	      do l = lij+1, lij+1 + ind_loop(2) - ind_loop(1)

	               lt = l  - ltij
	               lvo= lt
	               lf = l  - lfij
	               lx = l  - lxij
                              

!c...............Metrique
	               tcx = ti(lt +v1mtr)
	               tcy = ti(lt +v2mtr)
	               si  = sqrt (tcx*tcx + tcy*tcy)

	               nm  = l -  inci
	               nm2 = l -2*inci
	               np  = l +  inci

! Pente (qm) a l'interface droite et  (qp) a l'interface gauche

c    ! qm: right state,  qp: left state
        qm1=  c4*rop(l +v1) +  c5*rop(nm  +v1) +  c6*rop(np +v1)
        qp1=  c4*rop(nm+v1) +  c6*rop(nm2 +v1) +  c5*rop(l  +v1)

c    ! qm: right state,  qp: left state
        qm2=  c4*rop(l +v2) +  c5*rop(nm  +v2) +  c6*rop(np +v2)
        qp2=  c4*rop(nm+v2) +  c6*rop(nm2 +v2) +  c5*rop(l  +v2)

c    ! qm: right state,  qp: left state
        qm3=  c4*rop(l +v3) +  c5*rop(nm  +v3) +  c6*rop(np +v3)
        qp3=  c4*rop(nm+v3) +  c6*rop(nm2 +v3) +  c5*rop(l  +v3)

c    ! qm: right state,  qp: left state
        qm5=  c4*rop(l +v5) +  c5*rop(nm  +v5) +  c6*rop(np +v5)
        qp5=  c4*rop(nm+v5) +  c6*rop(nm2 +v5) +  c5*rop(l  +v5)

! Determination etat gauche (rou1) et droit (rou2): ro, roui, roe+p
        r1  =qp1
        rou1=r1*qp2
        rov1=r1*qp3
        p1  =r1*qp5*rgp
        h1  =gam1*p1 + .5*(rou1*qp2+rov1*qp3)

! Determination etat droite: ro, roui, roe+p
        r2  =qm1
        rou2=r2*qm2
        rov2=r2*qm3
        p2  =r2*qm5*rgp
        h2  =gam1*p2 + .5*(rou2*qm2+rov2*qm3)


! Determination vitesse normale interface
        qn1=qp2*tcx+qp3*tcy
        qn2=qm2*tcx+qm3*tcy

! Modification de vitesse normale par ajout
! de stabilisation de type Rhie-Chow
        c   = rgp*gam*rop(l +v5)  !c^2
        son = sqrt(qn1*qn1 / c)
        tam = c3*son+si
        tam1= max(0.,tam)*c2 ! fct amortissement: c3*Mach+1
        u   = 0.25*(qn1+qn2)-tam1*(p2-p1)
        tdu = max(abs(u),c1*si)

        p1p2= (p1+p2)*0.5


C****************************************************************************************
c************              CALCUL DE dpflu_dprop (dpeff_dprop)           ****************
c****************    LE CALCUL SE FAIT DE MANIERE ASCENDANTE    *************************
C****************************************************************************************

c*********************************** FLU2 ********************************************************

c    Equation de flu2 :
c       flu2 = u*(rou1+rou2) - tdu*(rou2-rou1) + tcx*p1p2

c************ ETAPE 1 : Calcul des derivees de flu2 par rapport a u, tdu et p1p2

        dpflu2_dpu = rou1 + rou2
        dpflu2_dptdu = rou1 - rou2
        dpflu2_dpp1p2 = tcx


c************ ETAPE 2 : Calcul de la derivee de flu2 par rapport a 
c                       r1,r2, rou1,rou2, rov1,rov2, p1,p2, qn1 et qn2.
c            En fait : Liste des variables pour flu2 : 
c                                                     ro1,rou2,p1,p2,qn1,qn2

c            Calcul des derivees partielles intermediaires necessaires
c            au calcul de la derivee de flu2.
        IF (ABS(u).GE.c1*si) THEN
           dptdu_dpu = sign(1.,u)
        ELSE
           dptdu_dpu = 0.
        END IF

        dpu_dpqn1 = 0.25
        dpu_dptam1 = -(p2-p1)
        dptam1_dptam = 0.5*(1.+sign(1.,tam))*c2
        dptam_dpson = c3
        dpson_dpqn1 = sign(1./sqrt(c),qn1)
        dpu_dpqn2 = 0.25
        dpp1p2_dpp1 = 0.5
        dpu_dpp1 = tam1
        dpp1p2_dpp2 = 0.5
        dpu_dpp2 = -tam1

c            Report des calculs precedents dans le calcul de la derivee de flu2.

        dpflu2_dpqn1 = dpflu2_dptdu * dptdu_dpu * ( dpu_dpqn1 +
     &                 dpu_dptam1 * dptam1_dptam * dptam_dpson *
     &                 dpson_dpqn1 )

        dpflu2_dpqn2 = dpflu2_dptdu * dptdu_dpu * dpu_dpqn2

        dpflu2_dpp1 = dpflu2_dpp1p2 * dpp1p2_dpp1 + dpflu2_dptdu *
     &                dptdu_dpu * dpu_dpp1

       dpflu2_dpp2 = dpflu2_dpp1p2 * dpp1p2_dpp2 + dpflu2_dptdu *
     &               dptdu_dpu * dpu_dpp2

       dpflu2_dprou1 = u+tdu
       dpflu2_dprou2 = u-tdu

       !  D'apres l'equation de flu2
       !     dpflu2_dpr1 = 0
       !     dpflu2_dpr2 = 0
       !     dpflu2_dprov1 = 0
       !     dpflu2_dprov2 = 0

c********* ETAPE 3 : Calcul des derivees de flu2 par rapport a
c                    qp1 qp2, qp3, qp5 et qm1, qm2,qm3, qm5
c  
c            Calcul des derivees partielles intermediaires necessaires
c            au calcul de la derivee de flu2.

c     ! qp : left state
       dpp1_dpqp1 = qp5*rgp
       dpp1_dpqp5 = r1*rgp
       dprou1_dpqp1 = qp2
       dprou1_dpqp2 = r1
       dpqn1_dpqp2 = tcx
       dpqn1_dpqp3 = tcy

c     ! qm : right state
       dpp2_dpqm1 = qm5*rgp
       dpp2_dpqm5 = r2*rgp
       dprou2_dpqm1 = qm2
       dprou2_dpqm2 = r2
       dpqn2_dpqm2 = tcx
       dpqn2_dpqm3 = tcy

c            Report des calculs precedents dans le calcul de la derivee de flu2.

c     ! qp : left state
       dpflu2_dpqp1 = dpflu2_dpp1 * dpp1_dpqp1 + dpflu2_dprou1 *
     &                dprou1_dpqp1
       dpflu2_dpqp2 = dpflu2_dpqn1 * dpqn1_dpqp2 + dpflu2_dprou1 *
     &                dprou1_dpqp2
       dpflu2_dpqp3 = dpflu2_dpqn1 * dpqn1_dpqp3
       dpflu2_dpqp5 = dpflu2_dpp1 * dpp1_dpqp5

c     ! qm : right state
       dpflu2_dpqm1 = dpflu2_dpp2 * dpp2_dpqm1 + dpflu2_dprou2 *
     &                dprou2_dpqm1
       dpflu2_dpqm2 = dpflu2_dpqn2 * dpqn2_dpqm2 + dpflu2_dprou2 *
     &                dprou2_dpqm2
       dpflu2_dpqm3 = dpflu2_dpqn2 * dpqn2_dpqm3
       dpflu2_dpqm5 = dpflu2_dpp2 * dpp2_dpqm5

c********* ETAPE 4 : Calcul de la derivee de flu2 par rapport a rop

c            Calcul des derivees partielles intermediaires necessaires
c            au calcul de la derivee de flu2.

       dpson_dpc = -0.5 * sqrt(qn1**2)/sqrt(c**3)
       dpflu2_dpc = dpflu2_dptdu * dptdu_dpu * dpu_dptam1 * dptam1_dptam 
     &            * dptam_dpson * dpson_dpc
       dpc_dprop_l5 = rgp * gam

c            Report des calculs precedents dans le calcul de la derivee de flu2.

       dpflu2_dprop1_np = dpflu2_dpqm1 * c6
       dpflu2_dprop1_l = dpflu2_dpqp1 * c5 + dpflu2_dpqm1 * c4
       dpflu2_dprop1_nm = dpflu2_dpqp1 * c4 + dpflu2_dpqm1 * c5
       dpflu2_dprop1_nm2 = dpflu2_dpqp1 * c6

       dpflu2_dprop2_np = dpflu2_dpqm2 * c6
       dpflu2_dprop2_l = dpflu2_dpqp2 * c5 + dpflu2_dpqm2 * c4
       dpflu2_dprop2_nm = dpflu2_dpqp2 * c4 + dpflu2_dpqm2 * c5
       dpflu2_dprop2_nm2 = dpflu2_dpqp2 * c6

       dpflu2_dprop3_np = dpflu2_dpqm3 * c6
       dpflu2_dprop3_l = dpflu2_dpqp3 * c5 + dpflu2_dpqm3 * c4
       dpflu2_dprop3_nm = dpflu2_dpqp3 * c4 + dpflu2_dpqm3 * c5
       dpflu2_dprop3_nm2 = dpflu2_dpqp3 * c6

       dpflu2_dprop5_np = dpflu2_dpqm5 * c6
       dpflu2_dprop5_l = dpflu2_dpqp5 * c5 + dpflu2_dpqm5 * c4 + 
     &                          dpflu2_dpc * dpc_dprop_l5 ! Prise en compte de c dans dpflu2_dprop5_l
       dpflu2_dprop5_nm = dpflu2_dpqp5 * c4 + dpflu2_dpqm5 * c5
       dpflu2_dprop5_nm2 = dpflu2_dpqp5 * c6


c******************************************************************************
!------------ Calcul des composants de dpRopl_dpRopldjr (dpPg_dpP)
c******************************************************************************

!      NORMALISATION de tcx, tcy et tcz
c       ci_mtr = 1. ! Valeur a integrer depuis BCWallInviscid.for 
c       cj_mtr = 1. ! Valeur a integrer depuis BCWallInviscid.for
c       ck_mtr = 0. ! Valeur a integrer depuis BCWallInviscid.for

      idir = param_int_eff(EFF_IDIR)
      neq_mtr = param_int(NEQ_IJ)

      call shape_tab_mtr(neq_mtr, param_int, idir,
     &                   ic,jc,kc,kc_vent,
     &                   ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale)

       ncx = tcx*ci_mtr
       ncy = tcy*cj_mtr
       ncz = tcz*ck_mtr

       surf_nc = sqrt(ncx*ncx + ncy*ncy + ncz*ncz)
       surf_nc = max(surf_nc,1e-30)

       ncx = ncx/surf_nc
       ncy = ncy/surf_nc
       ncz = ncz/surf_nc

!      Derivees calculees a partir des fonctions se trouvant dans : .FastS/BC/BCWallInviscid.for
       dpRop1_b_dpRop1_h = 1.

       dpRop2_b_dpRop2_h = 1. - 2.*(ncx**2)
       dpRop2_b_dpRop3_h = - 2.*ncy*ncx
       dpRop2_b_dpRop4_h = - 2.*ncz*ncx

       dpRop3_b_dpRop2_h = - 2.*ncx*ncy
       dpRop3_b_dpRop3_h = 1. - 2.*(ncy**2)
       dpRop3_b_dpRop4_h = - 2.*ncz*ncy

       dpRop4_b_dpRop2_h = - 2.*ncx*ncz
       dpRop4_b_dpRop3_h = - 2.*ncy*ncz
       dpRop4_b_dpRop4_h = 1. - 2.*(ncz**2)

       dpRop5_b_dpRop5_h = 1.


c******************************************************************************
!------------ Integration de dpPg_dpP dans dpflu2_dprop
c******************************************************************************

       dpflu2_dprop1_np = dpflu2_dprop1_np +
     &                    dpflu2_dprop1_nm2 * dpRop1_b_dpRop1_h
       dpflu2_dprop1_l = dpflu2_dprop1_l +
     &                    dpflu2_dprop1_nm * dpRop1_b_dpRop1_h

       dpflu2_dprop2_np = dpflu2_dprop2_np +
     &                    dpflu2_dprop2_nm2 * dpRop2_b_dpRop2_h +
     &                    dpflu2_dprop3_nm2 * dpRop3_b_dpRop2_h
       dpflu2_dprop2_l = dpflu2_dprop2_l +
     &                   dpflu2_dprop2_nm * dpRop2_b_dpRop2_h +
     &                   dpflu2_dprop3_nm * dpRop3_b_dpRop2_h

       dpflu2_dprop3_np = dpflu2_dprop3_np +
     &                    dpflu2_dprop2_nm2 * dpRop2_b_dpRop3_h +
     &                    dpflu2_dprop3_nm2 * dpRop3_b_dpRop3_h
       dpflu2_dprop3_l = dpflu2_dprop3_l +
     &                   dpflu2_dprop2_nm * dpRop2_b_dpRop3_h +
     &                   dpflu2_dprop3_nm * dpRop3_b_dpRop3_h

       dpflu2_dprop5_np = dpflu2_dprop5_np +
     &                    dpflu2_dprop5_nm2 * dpRop5_b_dpRop5_h
       dpflu2_dprop5_l = dpflu2_dprop5_l +
     &                   dpflu2_dprop5_nm * dpRop5_b_dpRop5_h


c****************************************************************************************************
c*********************************** FLU3 ***********************************************************

c    Equation de flu3 :
c       flu3 = u*(rov1+rov2) - tdu*(rov2-rov1) + tcx*p1p2

c********* ETAPE 1 : Calcul des derivees de flu3 par rapport a u, tdu et p1p2

       dpflu3_dpu = rov1 + rov2
       dpflu3_dptdu = rov1 - rov2
       dpflu3_dpp1p2 = tcy


c********* ETAPE 2 : Calcul de la derivee de flu3 par rapport a
c                    r1,r2, rou1,rou2, rov1,rov2, p1,p2, qn1 et qn2.
c            En fait : Liste des variables pour flu3 :
c                                                     rov1,rov2, p1,p2, qn1 et qn2.

c            Toutes les derivees partielles intermediaires pour le calcul
c            de la derivee de flu3 sont identiques a celles de flu2.

c            Report des derivees partielles intermediaires dans le calcul de la derivee de flu3.

       dpflu3_dpqn1 = dpflu3_dptdu * dptdu_dpu * ( dpu_dpqn1 +
     &                dpu_dptam1 * dptam1_dptam * dptam_dpson *
     &                dpson_dpqn1 )

       dpflu3_dpqn2 = dpflu3_dptdu * dptdu_dpu * dpu_dpqn2

       dpflu3_dpp1 = dpflu3_dpp1p2 * dpp1p2_dpp1 + dpflu3_dptdu *
     &               dptdu_dpu * dpu_dpp1

       dpflu3_dpp2 = dpflu3_dpp1p2 * dpp1p2_dpp2 + dpflu3_dptdu *
     &               dptdu_dpu * dpu_dpp2

       dpflu3_dprov1 = u+tdu
       dpflu3_dprov2 = u-tdu

       !  D'apres l'equation de flu3
       !     dpflu3_dpr1 = 0
       !     dpflu3_dpr2 = 0
       !     dpflu3_dprou1 = 0
       !     dpflu3_dprou2 = 0


c********* ETAPE 3 : Calcul des dérivées de flu2 par rapport à qp et qm

c            Calcul des derivees partielles intermediaires necessaires au calcul de la derivee de flu3.

c            Toutes les derivees partielles intermediaires pour le calcul de la derivee de flu3
c            sont identiques a celles de flu2, sauf pour les quatres suivantes :

c     ! qp : left state
       dprov1_dpqp1 = qp3
       dprov1_dpqp3 = r1

c     ! qm : right state
       dprov2_dpqm1 = qm3
       dprov2_dpqm3 = r2


c            Report des derivees partielles intermediaires dans le calcul de la derivee de flu3.

c     ! qp : left state
       dpflu3_dpqp1 = dpflu3_dpp1 * dpp1_dpqp1 + dpflu3_dprov1 *
     &                dprov1_dpqp1
       dpflu3_dpqp2 = dpflu3_dpqn1 * dpqn1_dpqp2
       dpflu3_dpqp3 = dpflu3_dpqn1 * dpqn1_dpqp3 + dpflu3_dprov1 *
     &                dprov1_dpqp3
       dpflu3_dpqp5 = dpflu3_dpp1 * dpp1_dpqp5

c    ! qm : right state
       dpflu3_dpqm1 = dpflu3_dpp2 * dpp2_dpqm1 + dpflu3_dprov2 *
     &                dprov2_dpqm1
       dpflu3_dpqm2 = dpflu3_dpqn2 * dpqn2_dpqm2
       dpflu3_dpqm3 = dpflu3_dpqn2 * dpqn2_dpqm3 + dpflu3_dprov2 *
     &                dprov2_dpqm3
       dpflu3_dpqm5 = dpflu3_dpp2 * dpp2_dpqm5


c********* ETAPE 4 : Calcul de la derivee de flu3 par rapport a rop

c            Toutes les derivees partielles intermediaires pour le calcul de la derivee de flu3
c            sont identiques a celles de flu2, sauf pour la suivante :

       dpflu3_dpc = dpflu3_dptdu * dptdu_dpu * dpu_dptam1 * dptam1_dptam 
     &            * dptam_dpson * dpson_dpc


c            Report des derivees partielles intermediaires dans le calcul de la derivee de flu3.

       dpflu3_dprop1_np = dpflu3_dpqm1 * c6
       dpflu3_dprop1_l = dpflu3_dpqp1 * c5 + dpflu3_dpqm1 * c4
       dpflu3_dprop1_nm = dpflu3_dpqp1 * c4 + dpflu3_dpqm1 * c5
       dpflu3_dprop1_nm2 = dpflu3_dpqp1 * c6

       dpflu3_dprop2_np = dpflu3_dpqm2 * c6
       dpflu3_dprop2_l = dpflu3_dpqp2 * c5 + dpflu3_dpqm2 * c4
       dpflu3_dprop2_nm = dpflu3_dpqp2 * c4 + dpflu3_dpqm2 * c5
       dpflu3_dprop2_nm2 = dpflu3_dpqp2 * c6

       dpflu3_dprop3_np = dpflu3_dpqm3 * c6
       dpflu3_dprop3_l = dpflu3_dpqp3 * c5 + dpflu3_dpqm3 * c4
       dpflu3_dprop3_nm = dpflu3_dpqp3 * c4 + dpflu3_dpqm3 * c5
       dpflu3_dprop3_nm2 = dpflu3_dpqp3 * c6

       dpflu3_dprop5_np = dpflu3_dpqm5 * c6
       dpflu3_dprop5_l = dpflu3_dpqp5 * c5 + dpflu3_dpqm5 * c4 + 
     &                          dpflu3_dpc * dpc_dprop_l5 ! Prise en compte de c dans dpflu3_dprop5_l
       dpflu3_dprop5_nm = dpflu3_dpqp5 * c4 + dpflu3_dpqm5 * c5
       dpflu3_dprop5_nm2 = dpflu3_dpqp5 * c6


c******************************************************************************
!------------ Calcul des composants de dpRopl_dpRopldjr (dpPg_dpP)
c******************************************************************************

!      Derivees calculees depuis BCWallInviscid.for sont ici identiques.


c******************************************************************************
!------------ Integration de dpPg_dpP dans dpflu3_dprop
c******************************************************************************

       dpflu3_dprop1_np = dpflu3_dprop1_np +
     &                    dpflu3_dprop1_nm2 * dpRop1_b_dpRop1_h
       dpflu3_dprop1_l = dpflu3_dprop1_l +
     &                   dpflu3_dprop1_nm * dpRop1_b_dpRop1_h

       dpflu3_dprop2_np = dpflu3_dprop2_np +
     &                    dpflu3_dprop2_nm2 * dpRop2_b_dpRop2_h +
     &                    dpflu3_dprop3_nm2 * dpRop3_b_dpRop2_h
       dpflu3_dprop2_l = dpflu3_dprop2_l +
     &                   dpflu3_dprop2_nm * dpRop2_b_dpRop2_h +
     &                   dpflu3_dprop3_nm * dpRop3_b_dpRop2_h

       dpflu3_dprop3_np = dpflu3_dprop3_np +
     &                    dpflu3_dprop2_nm2 * dpRop2_b_dpRop3_h +
     &                    dpflu3_dprop3_nm2 * dpRop3_b_dpRop3_h
       dpflu3_dprop3_l = dpflu3_dprop3_l +
     &                   dpflu3_dprop2_nm * dpRop2_b_dpRop3_h +
     &                   dpflu3_dprop3_nm * dpRop3_b_dpRop3_h

       dpflu3_dprop5_np = dpflu3_dprop5_np +
     &                   dpflu3_dprop5_nm2 * dpRop5_b_dpRop5_h
       dpflu3_dprop5_l = dpflu3_dprop5_l +
     &                   dpflu3_dprop5_nm * dpRop5_b_dpRop5_h


C****************************************************************************************
c***********         FIN DE CALCUL DE dpflu_dprop (dpeff_dprop)          ****************
C****************************************************************************************


C****************************************************************************************
c***********            CALCUL DE dpClp_dprop et dpCdp_dprop             ****************
C****************************************************************************************

      cc=sens*surfinv*cosAoA
      cs=sens*surfinv*sinAoA

c -------- Calcul de dpCdp_dprop

!        drag = (eff[0]*cosAoA+eff[1]*sinAoA)/(corde*span)

       dpCdp_dprop1_np = cc*dpflu2_dprop1_np +
     &                   cs*dpflu3_dprop1_np
       dpCdp_dprop1_l  = cc*dpflu2_dprop1_l +
     &                   cs*dpflu3_dprop1_l


       dpCdp_dprop2_np = cc*dpflu2_dprop2_np +
     &                   cs*dpflu3_dprop2_np
       dpCdp_dprop2_l  = cc*dpflu2_dprop2_l +
     &                   cs*dpflu3_dprop2_l

       dpCdp_dprop3_np = cc*dpflu2_dprop3_np +
     &                   cs*dpflu3_dprop3_np
       dpCdp_dprop3_l  = cc*dpflu2_dprop3_l +
     &                   cs*dpflu3_dprop3_l

       dpCdp_dprop5_np = cc*dpflu2_dprop5_np +
     &                   cs*dpflu3_dprop5_np
       dpCdp_dprop5_l  = cc*dpflu2_dprop5_l +
     &                   cs*dpflu3_dprop5_l


c -------- Calcul de dpClp_dprop

!        lift = (eff[1]*cosAoA-eff[0]*sinAoA)/(corde*span)

       dpClp_dprop1_np = cc*dpflu3_dprop1_np -
     &                   cs*dpflu2_dprop1_np
       dpClp_dprop1_l  = cc*dpflu3_dprop1_l  -
     &                   cs*dpflu2_dprop1_l

       dpClp_dprop2_np = cc*dpflu3_dprop2_np -
     &                   cs*dpflu2_dprop2_np
       dpClp_dprop2_l  = cc*dpflu3_dprop2_l  -
     &                   cs*dpflu2_dprop2_l

       dpClp_dprop3_np = cc*dpflu3_dprop3_np -
     &                   cs*dpflu2_dprop3_np
       dpClp_dprop3_l  = cc*dpflu3_dprop3_l -
     &                   cs*dpflu2_dprop3_l

       dpClp_dprop5_np = cc*dpflu3_dprop5_np -
     &                   cs*dpflu2_dprop5_np
       dpClp_dprop5_l  = cc*dpflu3_dprop5_l -
     &                   cs*dpflu2_dprop5_l

c ---------------- Pour VERIFICATION (comparaison avec les differences finies calculees).

      if ((iter.eq.29)) then


c------ Affichage des resultats pour les points (i=84,j=1) et (i=100,j=1)

       write(*,*) " Point (j= ", iter,"i= 1)"
       write(*,*) ""
       write(*,*) " dpCdp_dprop(l+v1) ",
     &            dpCdp_dprop1_l, " j ",
     & iter

       write(*,*) " dpClp_dprop(l+v1) ",
     &            dpClp_dprop1_l, " j ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dprop(l+v2) ",
     &            dpCdp_dprop2_l, " j ",
     & iter

       write(*,*) " dpClp_dprop(l+v2) ",
     &            dpClp_dprop2_l, " j ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dprop(l+v3) ",
     &            dpCdp_dprop3_l, " j ",
     & iter

       write(*,*) " dpClp_dprop(l+v3) ",
     &            dpClp_dprop3_l, " j ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dprop(l+v5) ",
     &            dpCdp_dprop5_l, " j ",
     & iter

       write(*,*) " dpClp_dprop(l+v5) ",
     &            dpClp_dprop5_l, " j ",
     & iter
       write(*,*) ""
       write(*,*) "-----------------------"

      endif


C****************************************************************************************
c***********        FIN DE CALCUL DE dpClp_dprop et dpCdp_dprop          ****************
C****************************************************************************************

#include "FastS/ADJOINT/AUSM/2d/dpP_dpW_fluausm_euler_2d_idir_min.for" 

C****************************************************************************************
c***********              CALCUL DE dpClp_dpW et dpCdp_dpW               ****************
C****************************************************************************************

c -------- Calcul de dpCdp_dpW
       dpCdp_dpW(np+v1) = dpCdp_dprop1_np * dprop1_dpW1_np + 
     &                    dpCdp_dprop2_np * dprop2_dpW1_np +
     &                    dpCdp_dprop3_np * dprop3_dpW1_np +
     &                    dpCdp_dprop5_np * dprop5_dpW1_np
       dpCdp_dpW(l+v1)  = dpCdp_dprop1_l * dprop1_dpW1_l + 
     &                    dpCdp_dprop2_l * dprop2_dpW1_l +
     &                    dpCdp_dprop3_l * dprop3_dpW1_l +
     &                    dpCdp_dprop5_l * dprop5_dpW1_l

       dpCdp_dpW(np+v2) = dpCdp_dprop2_np * dprop2_dpW2_np +
     &                    dpCdp_dprop5_np * dprop5_dpW2_np
       dpCdp_dpW(l+v2)  = dpCdp_dprop2_l * dprop2_dpW2_l +
     &                    dpCdp_dprop5_l * dprop5_dpW2_l

       dpCdp_dpW(np+v3) = dpCdp_dprop3_np * dprop3_dpW3_np +
     &                    dpCdp_dprop5_np * dprop5_dpW3_np
       dpCdp_dpW(l+v3)  = dpCdp_dprop3_l * dprop3_dpW3_l +
     &                    dpCdp_dprop5_l * dprop5_dpW3_l

       dpCdp_dpW(np+v5) = dpCdp_dprop5_np * dprop5_dpW5_np
       dpCdp_dpW(l+v5)  = dpCdp_dprop5_l * dprop5_dpW5_l

c -------- Calcul de dpClp_dpW
       dpClp_dpW(np+v1) = dpClp_dprop1_np * dprop1_dpW1_np + 
     &                    dpClp_dprop2_np * dprop2_dpW1_np +
     &                    dpClp_dprop3_np * dprop3_dpW1_np +
     &                    dpClp_dprop5_np * dprop5_dpW1_np
       dpClp_dpW(l+v1)  = dpClp_dprop1_l * dprop1_dpW1_l + 
     &                    dpClp_dprop2_l * dprop2_dpW1_l +
     &                    dpClp_dprop3_l * dprop3_dpW1_l +
     &                    dpClp_dprop5_l * dprop5_dpW1_l

       dpClp_dpW(np+v2) = dpClp_dprop2_np * dprop2_dpW2_np +
     &                    dpClp_dprop5_np * dprop5_dpW2_np
       dpClp_dpW(l+v2)  = dpClp_dprop2_l * dprop2_dpW2_l +
     &                    dpClp_dprop5_l * dprop5_dpW2_l

       dpClp_dpW(np+v3) = dpClp_dprop3_np * dprop3_dpW3_np +
     &                    dpClp_dprop5_np * dprop5_dpW3_np
       dpClp_dpW(l+v3)  = dpClp_dprop3_l * dprop3_dpW3_l +
     &                    dpClp_dprop5_l * dprop5_dpW3_l

       dpClp_dpW(np+v5) = dpClp_dprop5_np * dprop5_dpW5_np
       dpClp_dpW(l+v5)  = dpClp_dprop5_l * dprop5_dpW5_l




c ---------------- Pour VERIFICATION (comparaison avec les differences finies calculees)

      if ((iter.eq.29).or.(iter.eq.45)) then

c------ Affichage des resultats pour les points (i=84,j=1) et (i=100,j=1)

       write(*,*) " Point (j= ", iter,"i= 1)"
       write(*,*) ""
       write(*,*) " dpCdp_dpW(l+v1) ",
     &            dpCdp_dpW(l+v1), " j ",
     & iter

       write(*,*) " dpClp_dpW(l+v1) ",
     &            dpClp_dpW(l+v1), " j ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dpW(l+v2) ",
     &            dpCdp_dpW(l+v2), " j ",
     & iter

       write(*,*) " dpClp_dpW(l+v2) ",
     &            dpClp_dpW(l+v2), " j ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dpW(l+v3) ",
     &            dpCdp_dpW(l+v3), " j ",
     & iter

       write(*,*) " dpClp_dpW(l+v3) ",
     &            dpClp_dpW(l+v3), " j ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dpW(l+v5) ",
     &            dpCdp_dpW(l+v5), " j ",
     & iter

       write(*,*) " dpClp_dpW(l+v5) ",
     &            dpClp_dpW(l+v5), " j ",
     & iter
       write(*,*) ""
       write(*,*) "-----------------------"

c------ Affichage des resultats pour les points (i=84,j=2) et (i=100,j=2)
       write(*,*) " Point (j= ", iter,"i= 2)"
       write(*,*) ""
       write(*,*) " dpCdp_dpW(np+v1) ",
     &            dpCdp_dpW(np+v1), " j ",
     & iter

       write(*,*) " dpClp_dpW(np+v1) ",
     &            dpClp_dpW(np+v1), " j ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dpW(np+v2) ",
     &            dpCdp_dpW(np+v2), " j ",
     & iter

       write(*,*) " dpClp_dpW(np+v2) ",
     &            dpClp_dpW(np+v2), " j ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dpW(np+v3) ",
     &            dpCdp_dpW(np+v3), " j ",
     & iter

       write(*,*) " dpClp_dpW(np+v3) ",
     &            dpClp_dpW(np+v3), " j ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dpW(np+v5) ",
     &            dpCdp_dpW(np+v5), " j ",
     & iter

       write(*,*) " dpClp_dpW(np+v5) ",
     &            dpClp_dpW(np+v5), " j ",
     & iter
       write(*,*) ""
       write(*,*) "-----------------------"

      endif


C****************************************************************************************
c***********          FIN DE CALCUL DE dpClp_dpW et dpCdp_dpW            ****************
C****************************************************************************************

                    enddo
                 ENDDO
              ENDDO


      ELSEIF(param_int_eff(35).eq.2) THEN

              inc_x1 = param_int(10)                                   !(i,j+1,k  )
              inc_x2 = param_int(10)*param_int(10+1)       !(i,j  ,k+1)
              inc_x3 = inc_x1 + inc_x2
              iter = 0

              DO k = ind_loop(5), ind_loop(6)
                 DO j = ind_loop(3), ind_loop(4)

	      lij  =       inddm( ind_loop(1) , j, k) -1
	      ltij = lij - indmtr(ind_loop(1) , j, k) +1
	      lfij = lij - indflu(ind_loop(1) , j, k) +1
	      lxij = lij - indcg( ind_loop(1) , j, k) +1
              iter = iter + 1

!CC    !DIR$ ASSUME (mod(lij,4) .eq. 0)
!DIR$ IVDEP

	      do l = lij+1, lij+1 + ind_loop(2) - ind_loop(1)

	               lt = l  - ltij
	               lvo= lt
	               lf = l  - lfij
	               lx = l  - lxij
                              

!c...............Metrique
	               tcx = ti(lt +v1mtr)
	               tcy = ti(lt +v2mtr)
	               si  = sqrt (tcx*tcx + tcy*tcy)

	               nm  = l -  inci
	               nm2 = l -2*inci
	               np  = l +  inci

! Pente (qm) a l'interface droite et  (qp) a l'interface gauche

c    ! qm: right state,  qp: left state
        qm1=  c4*rop(l +v1) +  c5*rop(nm  +v1) +  c6*rop(np +v1)
        qp1=  c4*rop(nm+v1) +  c6*rop(nm2 +v1) +  c5*rop(l  +v1)

c    ! qm: right state,  qp: left state
        qm2=  c4*rop(l +v2) +  c5*rop(nm  +v2) +  c6*rop(np +v2)
        qp2=  c4*rop(nm+v2) +  c6*rop(nm2 +v2) +  c5*rop(l  +v2)

c    ! qm: right state,  qp: left state
        qm3=  c4*rop(l +v3) +  c5*rop(nm  +v3) +  c6*rop(np +v3)
        qp3=  c4*rop(nm+v3) +  c6*rop(nm2 +v3) +  c5*rop(l  +v3)

c    ! qm: right state,  qp: left state
        qm5=  c4*rop(l +v5) +  c5*rop(nm  +v5) +  c6*rop(np +v5)
        qp5=  c4*rop(nm+v5) +  c6*rop(nm2 +v5) +  c5*rop(l  +v5)

! Determination etat gauche (rou1) et droit (rou2): ro, roui, roe+p
        r1  =qp1
        rou1=r1*qp2
        rov1=r1*qp3
        p1  =r1*qp5*rgp
        h1  =gam1*p1 + .5*(rou1*qp2+rov1*qp3)

! Determination etat droite: ro, roui, roe+p
        r2  =qm1
        rou2=r2*qm2
        rov2=r2*qm3
        p2  =r2*qm5*rgp
        h2  =gam1*p2 + .5*(rou2*qm2+rov2*qm3)


! Determination vitesse normale interface
        qn1=qp2*tcx+qp3*tcy
        qn2=qm2*tcx+qm3*tcy

! Modification de vitesse normale par ajout
! de stabilisation de type Rhie-Chow
        c   = rgp*gam*rop(l +v5)  !c^2
        son = sqrt(qn1*qn1 / c)
        tam = c3*son+si
        tam1= max(0.,tam)*c2 ! fct amortissement: c3*Mach+1
        u   = 0.25*(qn1+qn2)-tam1*(p2-p1)
        tdu = max(abs(u),c1*si)

        p1p2= (p1+p2)*0.5


C****************************************************************************************
c************              CALCUL DE dpflu_dprop (dpeff_dprop)           ****************
c****************    LE CALCUL SE FAIT DE MANIERE ASCENDANTE    *************************
C****************************************************************************************

c*********************************** FLU2 ********************************************************

c    Equation de flu2 :
c       flu2 = u*(rou1+rou2) - tdu*(rou2-rou1) + tcx*p1p2

c************ ETAPE 1 : Calcul des derivees de flu2 par rapport a u, tdu et p1p2

        dpflu2_dpu = rou1 + rou2
        dpflu2_dptdu = rou1 - rou2
        dpflu2_dpp1p2 = tcx


c************ ETAPE 2 : Calcul de la derivee de flu2 par rapport a 
c                       r1,r2, rou1,rou2, rov1,rov2, p1,p2, qn1 et qn2.
c            En fait : Liste des variables pour flu2 : 
c                                                     ro1,rou2,p1,p2,qn1,qn2

c            Calcul des derivees partielles intermediaires necessaires
c            au calcul de la derivee de flu2.
        IF (ABS(u).GE.c1*si) THEN
           dptdu_dpu = sign(1.,u)
        ELSE
           dptdu_dpu = 0.
        END IF

        dpu_dpqn1 = 0.25
        dpu_dptam1 = -(p2-p1)
        dptam1_dptam = 0.5*(1.+sign(1.,tam))*c2
        dptam_dpson = c3
        dpson_dpqn1 = sign(1./sqrt(c),qn1)
        dpu_dpqn2 = 0.25
        dpp1p2_dpp1 = 0.5
        dpu_dpp1 = tam1
        dpp1p2_dpp2 = 0.5
        dpu_dpp2 = -tam1

c            Report des calculs precedents dans le calcul de la derivee de flu2.

        dpflu2_dpqn1 = dpflu2_dptdu * dptdu_dpu * ( dpu_dpqn1 +
     &                 dpu_dptam1 * dptam1_dptam * dptam_dpson *
     &                 dpson_dpqn1 )

        dpflu2_dpqn2 = dpflu2_dptdu * dptdu_dpu * dpu_dpqn2

        dpflu2_dpp1 = dpflu2_dpp1p2 * dpp1p2_dpp1 + dpflu2_dptdu *
     &                dptdu_dpu * dpu_dpp1

       dpflu2_dpp2 = dpflu2_dpp1p2 * dpp1p2_dpp2 + dpflu2_dptdu *
     &               dptdu_dpu * dpu_dpp2

       dpflu2_dprou1 = u+tdu
       dpflu2_dprou2 = u-tdu

       !  D'apres l'equation de flu2
       !     dpflu2_dpr1 = 0
       !     dpflu2_dpr2 = 0
       !     dpflu2_dprov1 = 0
       !     dpflu2_dprov2 = 0

c********* ETAPE 3 : Calcul des derivees de flu2 par rapport a
c                    qp1 qp2, qp3, qp5 et qm1, qm2,qm3, qm5
c  
c            Calcul des derivees partielles intermediaires necessaires
c            au calcul de la derivee de flu2.

c     ! qp : left state
       dpp1_dpqp1 = qp5*rgp
       dpp1_dpqp5 = r1*rgp
       dprou1_dpqp1 = qp2
       dprou1_dpqp2 = r1
       dpqn1_dpqp2 = tcx
       dpqn1_dpqp3 = tcy

c     ! qm : right state
       dpp2_dpqm1 = qm5*rgp
       dpp2_dpqm5 = r2*rgp
       dprou2_dpqm1 = qm2
       dprou2_dpqm2 = r2
       dpqn2_dpqm2 = tcx
       dpqn2_dpqm3 = tcy

c            Report des calculs precedents dans le calcul de la derivee de flu2.

c     ! qp : left state
       dpflu2_dpqp1 = dpflu2_dpp1 * dpp1_dpqp1 + dpflu2_dprou1 *
     &                dprou1_dpqp1
       dpflu2_dpqp2 = dpflu2_dpqn1 * dpqn1_dpqp2 + dpflu2_dprou1 *
     &                dprou1_dpqp2
       dpflu2_dpqp3 = dpflu2_dpqn1 * dpqn1_dpqp3
       dpflu2_dpqp5 = dpflu2_dpp1 * dpp1_dpqp5

c     ! qm : right state
       dpflu2_dpqm1 = dpflu2_dpp2 * dpp2_dpqm1 + dpflu2_dprou2 *
     &                dprou2_dpqm1
       dpflu2_dpqm2 = dpflu2_dpqn2 * dpqn2_dpqm2 + dpflu2_dprou2 *
     &                dprou2_dpqm2
       dpflu2_dpqm3 = dpflu2_dpqn2 * dpqn2_dpqm3
       dpflu2_dpqm5 = dpflu2_dpp2 * dpp2_dpqm5

c********* ETAPE 4 : Calcul de la derivee de flu2 par rapport a rop

c            Calcul des derivees partielles intermediaires necessaires
c            au calcul de la derivee de flu2.

       dpson_dpc = -0.5 * sqrt(qn1**2)/sqrt(c**3)
       dpflu2_dpc = dpflu2_dptdu * dptdu_dpu * dpu_dptam1 * dptam1_dptam 
     &            * dptam_dpson * dpson_dpc
       dpc_dprop_l5 = rgp * gam

c            Report des calculs precedents dans le calcul de la derivee de flu2.

       dpflu2_dprop1_np = dpflu2_dpqm1 * c6
       dpflu2_dprop1_l = dpflu2_dpqp1 * c5 + dpflu2_dpqm1 * c4
       dpflu2_dprop1_nm = dpflu2_dpqp1 * c4 + dpflu2_dpqm1 * c5
       dpflu2_dprop1_nm2 = dpflu2_dpqp1 * c6

       dpflu2_dprop2_np = dpflu2_dpqm2 * c6
       dpflu2_dprop2_l = dpflu2_dpqp2 * c5 + dpflu2_dpqm2 * c4
       dpflu2_dprop2_nm = dpflu2_dpqp2 * c4 + dpflu2_dpqm2 * c5
       dpflu2_dprop2_nm2 = dpflu2_dpqp2 * c6

       dpflu2_dprop3_np = dpflu2_dpqm3 * c6
       dpflu2_dprop3_l = dpflu2_dpqp3 * c5 + dpflu2_dpqm3 * c4
       dpflu2_dprop3_nm = dpflu2_dpqp3 * c4 + dpflu2_dpqm3 * c5
       dpflu2_dprop3_nm2 = dpflu2_dpqp3 * c6

       dpflu2_dprop5_np = dpflu2_dpqm5 * c6
       dpflu2_dprop5_l = dpflu2_dpqp5 * c5 + dpflu2_dpqm5 * c4 + 
     &                          dpflu2_dpc * dpc_dprop_l5 ! Prise en compte de c dans dpflu2_dprop5_l
       dpflu2_dprop5_nm = dpflu2_dpqp5 * c4 + dpflu2_dpqm5 * c5
       dpflu2_dprop5_nm2 = dpflu2_dpqp5 * c6


c******************************************************************************
!------------ Calcul des composants de dpRopl_dpRopldjr (dpPg_dpP)
c******************************************************************************

!      NORMALISATION de tcx, tcy et tcz
c       ci_mtr = 1. ! Valeur a integrer depuis BCWallInviscid.for 
c       cj_mtr = 1. ! Valeur a integrer depuis BCWallInviscid.for
c       ck_mtr = 0. ! Valeur a integrer depuis BCWallInviscid.for

      idir = param_int_eff(EFF_IDIR)
      neq_mtr = param_int(NEQ_IJ)

      call shape_tab_mtr(neq_mtr, param_int, idir,
     &                   ic,jc,kc,kc_vent,
     &                   ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale)

       ncx = tcx*ci_mtr
       ncy = tcy*cj_mtr
       ncz = tcz*ck_mtr

       surf_nc = sqrt(ncx*ncx + ncy*ncy + ncz*ncz)
       surf_nc = max(surf_nc,1e-30)

       ncx = ncx/surf_nc
       ncy = ncy/surf_nc
       ncz = ncz/surf_nc

!      Derivees calculees a partir des fonctions se trouvant dans : .FastS/BC/BCWallInviscid.for
       dpRop1_b_dpRop1_h = 1.

       dpRop2_b_dpRop2_h = 1. - 2.*(ncx**2)
       dpRop2_b_dpRop3_h = - 2.*ncy*ncx
       dpRop2_b_dpRop4_h = - 2.*ncz*ncx

       dpRop3_b_dpRop2_h = - 2.*ncx*ncy
       dpRop3_b_dpRop3_h = 1. - 2.*(ncy**2)
       dpRop3_b_dpRop4_h = - 2.*ncz*ncy

       dpRop4_b_dpRop2_h = - 2.*ncx*ncz
       dpRop4_b_dpRop3_h = - 2.*ncy*ncz
       dpRop4_b_dpRop4_h = 1. - 2.*(ncz**2)

       dpRop5_b_dpRop5_h = 1.


c******************************************************************************
!------------ Integration de dpPg_dpP dans dpflu2_dprop
c******************************************************************************

       dpflu2_dprop1_nm2 = dpflu2_dprop1_nm2 +
     &                    dpflu2_dprop1_np * dpRop1_b_dpRop1_h
       dpflu2_dprop1_nm = dpflu2_dprop1_nm +
     &                    dpflu2_dprop1_l * dpRop1_b_dpRop1_h

       dpflu2_dprop2_nm2 = dpflu2_dprop2_nm2 +
     &                    dpflu2_dprop2_np * dpRop2_b_dpRop2_h +
     &                    dpflu2_dprop3_np * dpRop3_b_dpRop2_h
       dpflu2_dprop2_nm = dpflu2_dprop2_nm +
     &                   dpflu2_dprop2_l * dpRop2_b_dpRop2_h +
     &                   dpflu2_dprop3_l * dpRop3_b_dpRop2_h

       dpflu2_dprop3_nm2 = dpflu2_dprop3_nm2 +
     &                    dpflu2_dprop2_np * dpRop2_b_dpRop3_h +
     &                    dpflu2_dprop3_np * dpRop3_b_dpRop3_h
       dpflu2_dprop3_nm = dpflu2_dprop3_nm +
     &                   dpflu2_dprop2_l * dpRop2_b_dpRop3_h +
     &                   dpflu2_dprop3_l * dpRop3_b_dpRop3_h

       dpflu2_dprop5_nm2 = dpflu2_dprop5_nm2 +
     &                    dpflu2_dprop5_np * dpRop5_b_dpRop5_h
       dpflu2_dprop5_nm = dpflu2_dprop5_nm +
     &                   dpflu2_dprop5_l * dpRop5_b_dpRop5_h


c****************************************************************************************************
c*********************************** FLU3 ***********************************************************

c    Equation de flu3 :
c       flu3 = u*(rov1+rov2) - tdu*(rov2-rov1) + tcx*p1p2

c********* ETAPE 1 : Calcul des derivees de flu3 par rapport a u, tdu et p1p2

       dpflu3_dpu = rov1 + rov2
       dpflu3_dptdu = rov1 - rov2
       dpflu3_dpp1p2 = tcy


c********* ETAPE 2 : Calcul de la derivee de flu3 par rapport a
c                    r1,r2, rou1,rou2, rov1,rov2, p1,p2, qn1 et qn2.
c            En fait : Liste des variables pour flu3 :
c                                                     rov1,rov2, p1,p2, qn1 et qn2.

c            Toutes les derivees partielles intermediaires pour le calcul
c            de la derivee de flu3 sont identiques a celles de flu2.

c            Report des derivees partielles intermediaires dans le calcul de la derivee de flu3.

       dpflu3_dpqn1 = dpflu3_dptdu * dptdu_dpu * ( dpu_dpqn1 +
     &                dpu_dptam1 * dptam1_dptam * dptam_dpson *
     &                dpson_dpqn1 )

       dpflu3_dpqn2 = dpflu3_dptdu * dptdu_dpu * dpu_dpqn2

       dpflu3_dpp1 = dpflu3_dpp1p2 * dpp1p2_dpp1 + dpflu3_dptdu *
     &               dptdu_dpu * dpu_dpp1

       dpflu3_dpp2 = dpflu3_dpp1p2 * dpp1p2_dpp2 + dpflu3_dptdu *
     &               dptdu_dpu * dpu_dpp2

       dpflu3_dprov1 = u+tdu
       dpflu3_dprov2 = u-tdu

       !  D'apres l'equation de flu3
       !     dpflu3_dpr1 = 0
       !     dpflu3_dpr2 = 0
       !     dpflu3_dprou1 = 0
       !     dpflu3_dprou2 = 0


c********* ETAPE 3 : Calcul des dérivées de flu2 par rapport à qp et qm

c            Calcul des derivees partielles intermediaires necessaires au calcul de la derivee de flu3.

c            Toutes les derivees partielles intermediaires pour le calcul de la derivee de flu3
c            sont identiques a celles de flu2, sauf pour les quatres suivantes :

c     ! qp : left state
       dprov1_dpqp1 = qp3
       dprov1_dpqp3 = r1

c     ! qm : right state
       dprov2_dpqm1 = qm3
       dprov2_dpqm3 = r2


c            Report des derivees partielles intermediaires dans le calcul de la derivee de flu3.

c     ! qp : left state
       dpflu3_dpqp1 = dpflu3_dpp1 * dpp1_dpqp1 + dpflu3_dprov1 *
     &                dprov1_dpqp1
       dpflu3_dpqp2 = dpflu3_dpqn1 * dpqn1_dpqp2
       dpflu3_dpqp3 = dpflu3_dpqn1 * dpqn1_dpqp3 + dpflu3_dprov1 *
     &                dprov1_dpqp3
       dpflu3_dpqp5 = dpflu3_dpp1 * dpp1_dpqp5

c    ! qm : right state
       dpflu3_dpqm1 = dpflu3_dpp2 * dpp2_dpqm1 + dpflu3_dprov2 *
     &                dprov2_dpqm1
       dpflu3_dpqm2 = dpflu3_dpqn2 * dpqn2_dpqm2
       dpflu3_dpqm3 = dpflu3_dpqn2 * dpqn2_dpqm3 + dpflu3_dprov2 *
     &                dprov2_dpqm3
       dpflu3_dpqm5 = dpflu3_dpp2 * dpp2_dpqm5


c********* ETAPE 4 : Calcul de la derivee de flu3 par rapport a rop

c            Toutes les derivees partielles intermediaires pour le calcul de la derivee de flu3
c            sont identiques a celles de flu2, sauf pour la suivante :

       dpflu3_dpc = dpflu3_dptdu * dptdu_dpu * dpu_dptam1 * dptam1_dptam 
     &            * dptam_dpson * dpson_dpc


c            Report des derivees partielles intermediaires dans le calcul de la derivee de flu3.

       dpflu3_dprop1_np = dpflu3_dpqm1 * c6
       dpflu3_dprop1_l = dpflu3_dpqp1 * c5 + dpflu3_dpqm1 * c4
       dpflu3_dprop1_nm = dpflu3_dpqp1 * c4 + dpflu3_dpqm1 * c5
       dpflu3_dprop1_nm2 = dpflu3_dpqp1 * c6

       dpflu3_dprop2_np = dpflu3_dpqm2 * c6
       dpflu3_dprop2_l = dpflu3_dpqp2 * c5 + dpflu3_dpqm2 * c4
       dpflu3_dprop2_nm = dpflu3_dpqp2 * c4 + dpflu3_dpqm2 * c5
       dpflu3_dprop2_nm2 = dpflu3_dpqp2 * c6

       dpflu3_dprop3_np = dpflu3_dpqm3 * c6
       dpflu3_dprop3_l = dpflu3_dpqp3 * c5 + dpflu3_dpqm3 * c4
       dpflu3_dprop3_nm = dpflu3_dpqp3 * c4 + dpflu3_dpqm3 * c5
       dpflu3_dprop3_nm2 = dpflu3_dpqp3 * c6

       dpflu3_dprop5_np = dpflu3_dpqm5 * c6
       dpflu3_dprop5_l = dpflu3_dpqp5 * c5 + dpflu3_dpqm5 * c4 + 
     &                          dpflu3_dpc * dpc_dprop_l5 ! Prise en compte de c dans dpflu3_dprop5_l
       dpflu3_dprop5_nm = dpflu3_dpqp5 * c4 + dpflu3_dpqm5 * c5
       dpflu3_dprop5_nm2 = dpflu3_dpqp5 * c6


c******************************************************************************
!------------ Calcul des composants de dpRopl_dpRopldjr (dpPg_dpP)
c******************************************************************************

!      Derivees calculees depuis BCWallInviscid.for sont ici identiques.


c******************************************************************************
!------------ Integration de dpPg_dpP dans dpflu3_dprop
c******************************************************************************

       dpflu3_dprop1_nm2 = dpflu3_dprop1_nm2 +
     &                    dpflu3_dprop1_np * dpRop1_b_dpRop1_h
       dpflu3_dprop1_nm = dpflu3_dprop1_nm +
     &                   dpflu3_dprop1_l * dpRop1_b_dpRop1_h

       dpflu3_dprop2_nm2 = dpflu3_dprop2_nm2 +
     &                    dpflu3_dprop2_np * dpRop2_b_dpRop2_h +
     &                    dpflu3_dprop3_np * dpRop3_b_dpRop2_h
       dpflu3_dprop2_nm = dpflu3_dprop2_nm +
     &                   dpflu3_dprop2_l * dpRop2_b_dpRop2_h +
     &                   dpflu3_dprop3_l * dpRop3_b_dpRop2_h

       dpflu3_dprop3_nm2 = dpflu3_dprop3_nm2 +
     &                    dpflu3_dprop2_np * dpRop2_b_dpRop3_h +
     &                    dpflu3_dprop3_np * dpRop3_b_dpRop3_h
       dpflu3_dprop3_nm = dpflu3_dprop3_nm +
     &                   dpflu3_dprop2_l * dpRop2_b_dpRop3_h +
     &                   dpflu3_dprop3_l * dpRop3_b_dpRop3_h

       dpflu3_dprop5_nm2 = dpflu3_dprop5_nm2 +
     &                   dpflu3_dprop5_np * dpRop5_b_dpRop5_h
       dpflu3_dprop5_nm = dpflu3_dprop5_nm +
     &                   dpflu3_dprop5_l * dpRop5_b_dpRop5_h


C****************************************************************************************
c***********         FIN DE CALCUL DE dpflu_dprop (dpeff_dprop)          ****************
C****************************************************************************************


C****************************************************************************************
c***********            CALCUL DE dpClp_dprop et dpCdp_dprop             ****************
C****************************************************************************************

      cc=sens*surfinv*cosAoA
      cs=sens*surfinv*sinAoA

c -------- Calcul de dpCdp_dprop

!        drag = (eff[0]*cosAoA+eff[1]*sinAoA)/(corde*span)

       dpCdp_dprop1_nm2 = cc*dpflu2_dprop1_nm2 +
     &                   cs*dpflu3_dprop1_nm2
       dpCdp_dprop1_nm  = cc*dpflu2_dprop1_nm +
     &                   cs*dpflu3_dprop1_nm


       dpCdp_dprop2_nm2 = cc*dpflu2_dprop2_nm2 +
     &                   cs*dpflu3_dprop2_nm2
       dpCdp_dprop2_nm  = cc*dpflu2_dprop2_nm +
     &                   cs*dpflu3_dprop2_nm

       dpCdp_dprop3_nm2 = cc*dpflu2_dprop3_nm2 +
     &                   cs*dpflu3_dprop3_nm2
       dpCdp_dprop3_nm  = cc*dpflu2_dprop3_nm +
     &                   cs*dpflu3_dprop3_nm

       dpCdp_dprop5_nm2 = cc*dpflu2_dprop5_nm2 +
     &                   cs*dpflu3_dprop5_nm2
       dpCdp_dprop5_nm  = cc*dpflu2_dprop5_nm +
     &                   cs*dpflu3_dprop5_nm


c -------- Calcul de dpClp_dprop

!        lift = (eff[1]*cosAoA-eff[0]*sinAoA)/(corde*span)

       dpClp_dprop1_nm2 = cc*dpflu3_dprop1_nm2 -
     &                   cs*dpflu2_dprop1_nm2
       dpClp_dprop1_nm  = cc*dpflu3_dprop1_nm  -
     &                   cs*dpflu2_dprop1_nm

       dpClp_dprop2_nm2 = cc*dpflu3_dprop2_nm2 -
     &                   cs*dpflu2_dprop2_nm2
       dpClp_dprop2_nm  = cc*dpflu3_dprop2_nm  -
     &                   cs*dpflu2_dprop2_nm

       dpClp_dprop3_nm2 = cc*dpflu3_dprop3_nm2 -
     &                   cs*dpflu2_dprop3_nm2
       dpClp_dprop3_nm  = cc*dpflu3_dprop3_nm -
     &                   cs*dpflu2_dprop3_nm

       dpClp_dprop5_nm2 = cc*dpflu3_dprop5_nm2 -
     &                   cs*dpflu2_dprop5_nm2
       dpClp_dprop5_nm  = cc*dpflu3_dprop5_nm -
     &                   cs*dpflu2_dprop5_nm

c ---------------- Pour VERIFICATION (comparaison avec les differences finies calculees).

      if (iter.eq.100) then

c------ Affichage des resultats pour les points (i=84,j=1) et (i=100,j=1)

       write(*,*) " Point (j= ", iter,"j= 128)"
       write(*,*) ""
       write(*,*) " dpCdp_dprop(nm+v1) ",
     &            dpCdp_dprop1_nm, " j ",
     & iter

       write(*,*) " dpClp_dprop(nm+v1) ",
     &            dpClp_dprop1_nm, " j ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dprop(nm+v2) ",
     &            dpCdp_dprop2_nm, " j ",
     & iter

       write(*,*) " dpClp_dprop(nm+v2) ",
     &            dpClp_dprop2_nm, " j ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dprop(nm+v3) ",
     &            dpCdp_dprop3_nm, " j ",
     & iter

       write(*,*) " dpClp_dprop(nm+v3) ",
     &            dpClp_dprop3_nm, " j ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dprop(nm+v5) ",
     &            dpCdp_dprop5_nm, " j ",
     & iter

       write(*,*) " dpClp_dprop(nm+v5) ",
     &            dpClp_dprop5_nm, " j ",
     & iter
       write(*,*) ""
       write(*,*) "-----------------------"

      endif


C****************************************************************************************
c***********        FIN DE CALCUL DE dpClp_dprop et dpCdp_dprop          ****************
C****************************************************************************************

#include "FastS/ADJOINT/AUSM/2d/dpP_dpW_fluausm_euler_2d_idir_max.for" 

C****************************************************************************************
c***********              CALCUL DE dpClp_dpW et dpCdp_dpW               ****************
C****************************************************************************************

c -------- Calcul de dpCdp_dpW
       dpCdp_dpW(nm2+v1) = dpCdp_dprop1_nm2 * dprop1_dpW1_nm2 + 
     &                    dpCdp_dprop2_nm2 * dprop2_dpW1_nm2 +
     &                    dpCdp_dprop3_nm2 * dprop3_dpW1_nm2 +
     &                    dpCdp_dprop5_nm2 * dprop5_dpW1_nm2
       dpCdp_dpW(nm+v1)  = dpCdp_dprop1_nm * dprop1_dpW1_nm + 
     &                    dpCdp_dprop2_nm * dprop2_dpW1_nm +
     &                    dpCdp_dprop3_nm * dprop3_dpW1_nm +
     &                    dpCdp_dprop5_nm * dprop5_dpW1_nm

       dpCdp_dpW(nm2+v2) = dpCdp_dprop2_nm2 * dprop2_dpW2_nm2 +
     &                    dpCdp_dprop5_nm2 * dprop5_dpW2_nm2
       dpCdp_dpW(nm+v2)  = dpCdp_dprop2_nm * dprop2_dpW2_nm +
     &                    dpCdp_dprop5_nm * dprop5_dpW2_nm

       dpCdp_dpW(nm2+v3) = dpCdp_dprop3_nm2 * dprop3_dpW3_nm2 +
     &                    dpCdp_dprop5_nm2 * dprop5_dpW3_nm2
       dpCdp_dpW(nm+v3)  = dpCdp_dprop3_nm * dprop3_dpW3_nm +
     &                    dpCdp_dprop5_nm * dprop5_dpW3_nm

       dpCdp_dpW(nm2+v5) = dpCdp_dprop5_nm2 * dprop5_dpW5_nm2
       dpCdp_dpW(nm+v5)  = dpCdp_dprop5_nm * dprop5_dpW5_nm


c -------- Calcul de dpClp_dpW
       dpClp_dpW(nm2+v1) = dpClp_dprop1_nm2 * dprop1_dpW1_nm2 + 
     &                    dpClp_dprop2_nm2 * dprop2_dpW1_nm2 +
     &                    dpClp_dprop3_nm2 * dprop3_dpW1_nm2 +
     &                    dpClp_dprop5_nm2 * dprop5_dpW1_nm2
       dpClp_dpW(nm+v1)  = dpClp_dprop1_nm * dprop1_dpW1_nm + 
     &                    dpClp_dprop2_nm * dprop2_dpW1_nm +
     &                    dpClp_dprop3_nm * dprop3_dpW1_nm +
     &                    dpClp_dprop5_nm * dprop5_dpW1_nm

       dpClp_dpW(nm2+v2) = dpClp_dprop2_nm2 * dprop2_dpW2_nm2 +
     &                    dpClp_dprop5_nm2 * dprop5_dpW2_nm2
       dpClp_dpW(nm+v2)  = dpClp_dprop2_nm * dprop2_dpW2_nm +
     &                    dpClp_dprop5_nm * dprop5_dpW2_nm

       dpClp_dpW(nm2+v3) = dpClp_dprop3_nm2 * dprop3_dpW3_nm2 +
     &                    dpClp_dprop5_nm2 * dprop5_dpW3_nm2
       dpClp_dpW(nm+v3)  = dpClp_dprop3_nm * dprop3_dpW3_nm +
     &                    dpClp_dprop5_nm * dprop5_dpW3_nm

       dpClp_dpW(nm2+v5) = dpClp_dprop5_nm2 * dprop5_dpW5_nm2
       dpClp_dpW(nm+v5)  = dpClp_dprop5_nm * dprop5_dpW5_nm




c ---------------- Pour VERIFICATION (comparaison avec les differences finies calculees)

      if ((iter.eq.84).or.(iter.eq.100)) then

c------ Affichage des resultats pour les points (i=84,j=1) et (i=100,j=1)

       write(*,*) " Point (j= ", iter,"i= 127)"
       write(*,*) ""
       write(*,*) " dpCdp_dpW(nm+v1) ",
     &            dpCdp_dpW(nm+v1), " j ",
     & iter

       write(*,*) " dpClp_dpW(nm+v1) ",
     &            dpClp_dpW(nm+v1), " j ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dpW(nm+v2) ",
     &            dpCdp_dpW(nm+v2), " j ",
     & iter

       write(*,*) " dpClp_dpW(nm+v2) ",
     &            dpClp_dpW(nm+v2), " j ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dpW(nm+v3) ",
     &            dpCdp_dpW(nm+v3), " j ",
     & iter

       write(*,*) " dpClp_dpW(nm+v3) ",
     &            dpClp_dpW(nm+v3), " j ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dpW(nm+v5) ",
     &            dpCdp_dpW(nm+v5), " j ",
     & iter

       write(*,*) " dpClp_dpW(nm+v5) ",
     &            dpClp_dpW(nm+v5), " j ",
     & iter
       write(*,*) ""
       write(*,*) "-----------------------"

c------ Affichage des resultats pour les points (i=84,j=2) et (i=100,j=2)
       write(*,*) " Point (j= ", iter,"i= 128)"
       write(*,*) ""
       write(*,*) " dpCdp_dpW(nm2+v1) ",
     &            dpCdp_dpW(nm2+v1), " j ",
     & iter

       write(*,*) " dpClp_dpW(nm2+v1) ",
     &            dpClp_dpW(nm2+v1), " j ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dpW(nm2+v2) ",
     &            dpCdp_dpW(nm2+v2), " j ",
     & iter

       write(*,*) " dpClp_dpW(nm2+v2) ",
     &            dpClp_dpW(nm2+v2), " j ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dpW(nm2+v3) ",
     &            dpCdp_dpW(nm2+v3), " j ",
     & iter

       write(*,*) " dpClp_dpW(nm2+v3) ",
     &            dpClp_dpW(nm2+v3), " j ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dpW(nm2+v5) ",
     &            dpCdp_dpW(nm2+v5), " j ",
     & iter

       write(*,*) " dpClp_dpW(nm2+v5) ",
     &            dpClp_dpW(nm2+v5), " j ",
     & iter
       write(*,*) ""
       write(*,*) "-----------------------"

      endif


C****************************************************************************************
c***********          FIN DE CALCUL DE dpClp_dpW et dpCdp_dpW            ****************
C****************************************************************************************

                    enddo
                 ENDDO
              ENDDO


!************** FIN DE LA BOUCLE CALCULANT LES EFFORTS SUR LA FACE I ****************
!************************************************************************************

      ELSEIF(param_int_eff(35).eq.3) THEN

!************************************************************************************
!*********************** CALCUL DES EFFORTS SUR LA FACE J ***************************
!************************************************************************************

              inc_x1 = 1                             !(i+1,j,k  )
              inc_x2 = param_int(10)*param_int(10+1) !(i  ,j,k+1)
              inc_x3 = inc_x1 + inc_x2

              DO k = ind_loop(5), ind_loop(6)
                 DO j = ind_loop(3), ind_loop(4)

	      lij  =       inddm( ind_loop(1) , j, k) -1
	      ltij = lij - indmtr(ind_loop(1) , j, k) +1
	      lfij = lij - indflu(ind_loop(1) , j, k) +1
	      lxij = lij - indcg( ind_loop(1) , j, k) +1

!CC    !DIR$ ASSUME (mod(lij,4) .eq. 0)
!DIR$ IVDEP

	      do l = lij+1, lij+1 + ind_loop(2) - ind_loop(1)

	               lt = l  - ltij
	               lvo= lt
	               lf = l  - lfij
	               lx = l  - lxij


!c...............Metrique
	               tcx = tj(lt +v1mtr)
	               tcy = tj(lt +v2mtr)
	               sj  = sqrt (tcx*tcx + tcy*tcy)

	               nm  = l -  incj
	               nm2 = l -2*incj
	               np  = l +  incj

! Pente (qm) a l'interface droite et  (qp) a l'interface gauche

c    ! qm: right state,  qp: left state
        qm1=  c4*rop(l +v1) +  c5*rop(nm  +v1) +  c6*rop(np +v1)
        qp1=  c4*rop(nm+v1) +  c6*rop(nm2 +v1) +  c5*rop(l  +v1)

c    ! qm: right state,  qp: left state
        qm2=  c4*rop(l +v2) +  c5*rop(nm  +v2) +  c6*rop(np +v2)
        qp2=  c4*rop(nm+v2) +  c6*rop(nm2 +v2) +  c5*rop(l  +v2)

c    ! qm: right state,  qp: left state
        qm3=  c4*rop(l +v3) +  c5*rop(nm  +v3) +  c6*rop(np +v3)
        qp3=  c4*rop(nm+v3) +  c6*rop(nm2 +v3) +  c5*rop(l  +v3)

c    ! qm: right state,  qp: left state
        qm5=  c4*rop(l +v5) +  c5*rop(nm  +v5) +  c6*rop(np +v5)
        qp5=  c4*rop(nm+v5) +  c6*rop(nm2 +v5) +  c5*rop(l  +v5)

c ---------------- Pour VERIFICATION (comparaison avec les differences finies calculees).

      if(j.eq.1.and.l-lij+ind_loop(1)-1.eq.100) then


c------ Affichage des resultats pour les points (i=84,j=1) et (i=100,j=1)

       write(*,*) " Point (i= ", l-lij+ind_loop(1)-1,"j= 1)"
       write(*,*) ""

       write(*,*) "tcx ", tcx
       write(*,*) "tcy ", tcy
       write(*,*) "-----------------------"

       write(*,*) " qm1 ",
     &            qm1, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " qm2 ",
     &            qm2, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " qm3 ",
     &            qm3, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " qm5 ",
     &            qm5, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) "-----------------------"

       write(*,*) " qp1 ",
     &            qp1, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " qp2 ",
     &            qp2, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " qp3 ",
     &            qp3, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " qp5 ",
     &            qp5, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) "-----------------------"

      endif

! Determination etat gauche (rou1) et droit (rou2): ro, roui, roe+p
        r1   = qp1
        rou1 = r1*qp2
        rov1 = r1*qp3
        p1   = r1*qp5*rgp
        h1   = gam1*p1 + .5*(rou1*qp2+rov1*qp3)

! Determination etat droite: ro, roui, roe+p
        r2  =qm1
        rou2=r2*qm2
        rov2=r2*qm3
        p2  =r2*qm5*rgp
        h2  =gam1*p2 + .5*(rou2*qm2+rov2*qm3)

! Determination vitesse normale interface
        qn1=qp2*tcx+qp3*tcy
        qn2=qm2*tcx+qm3*tcy

! Modification de vitesse normale par ajout
! de stabilisation de type Rhie-Chow
        c   = rgp*gam*rop(l +v5)  !c^2
        son = sqrt(qn1*qn1 / c)
        tam = c3*son+sj
        tam1= max(0.,tam)*c2 ! fct amortissement: c3*Mach+1
        u   = 0.25*(qn1+qn2)-tam1*(p2-p1)
        tdu = max(abs(u),c1*sj)

        p1p2= (p1+p2)*0.5


C****************************************************************************************
c************              CALCUL DE dpflu_dprop (dpeff_dprop)           ****************
c****************    LE CALCUL SE FAIT DE MANIERE ASCENDANTE    *************************
C****************************************************************************************

c*********************************** FLU2 ********************************************************

c    Equation de flu2 :
c       flu2 = u*(rou1+rou2) - tdu*(rou2-rou1) + tcx*p1p2

c************ ETAPE 1 : Calcul des derivees de flu2 par rapport a u, tdu et p1p2

        dpflu2_dpu = rou1 + rou2
        dpflu2_dptdu = -(rou2 - rou1)
        dpflu2_dpp1p2 = tcx


c************ ETAPE 2 : Calcul de la derivee de flu2 par rapport a 
c                       r1,r2, rou1,rou2, rov1,rov2, p1,p2, qn1 et qn2.
c            En fait : Liste des variables pour flu2 : 
c                                                     ro1,rou2,p1,p2,qn1,qn2

c            Calcul des derivees partielles intermediaires necessaires
c            au calcul de la derivee de flu2.

        IF (ABS(u).GT.c1*sj) THEN
           dptdu_dpu = sign(1.,u)
        ELSE
           dptdu_dpu = 0.
        END IF

        dpu_dpqn1 = 0.25
        dpu_dptam1 = -(p2-p1)
        dptam1_dptam = 0.5*(1.+sign(1.,tam))*c2
        dptam_dpson = c3
        dpson_dpqn1 = sign(1./sqrt(c),qn1)
        dpu_dpqn2 = 0.25
        dpp1p2_dpp1 = 0.5
        dpu_dpp1 = tam1
        dpp1p2_dpp2 = 0.5
        dpu_dpp2 = -tam1


c            Report des calculs precedents dans le calcul de la derivee de flu2.

        dpflu2_dpqn1 = ( dpflu2_dptdu * dptdu_dpu + dpflu2_dpu ) *
     &                 ( dpu_dpqn1 + dpu_dptam1 * dptam1_dptam * 
     &                 dptam_dpson * dpson_dpqn1 )

        dpflu2_dpqn2 = ( dpflu2_dptdu * dptdu_dpu + dpflu2_dpu ) *
     &                 dpu_dpqn2

        dpflu2_dpp1 = dpflu2_dpp1p2 * dpp1p2_dpp1 + ( dpflu2_dptdu *
     &                dptdu_dpu + dpflu2_dpu ) * dpu_dpp1

       dpflu2_dpp2 = dpflu2_dpp1p2 * dpp1p2_dpp2 + ( dpflu2_dptdu *
     &               dptdu_dpu + dpflu2_dpu ) * dpu_dpp2

       dpflu2_dprou1 = u+tdu
       dpflu2_dprou2 = u-tdu

       !  D'apres l'equation de flu2
       !     dpflu2_dpr1 = 0
       !     dpflu2_dpr2 = 0
       !     dpflu2_dprov1 = 0
       !     dpflu2_dprov2 = 0

c********* ETAPE 3 : Calcul des derivees de flu2 par rapport a
c                    qp1 qp2, qp3, qp5 et qm1, qm2,qm3, qm5
c  
c            Calcul des derivees partielles intermediaires necessaires
c            au calcul de la derivee de flu2.

c     ! qp : left state
       dpp1_dpqp1 = qp5*rgp
       dpp1_dpqp5 = r1*rgp
       dprou1_dpqp1 = qp2
       dprou1_dpqp2 = r1
       dpqn1_dpqp2 = tcx
       dpqn1_dpqp3 = tcy

c     ! qm : right state
       dpp2_dpqm1 = qm5*rgp
       dpp2_dpqm5 = r2*rgp
       dprou2_dpqm1 = qm2
       dprou2_dpqm2 = r2
       dpqn2_dpqm2 = tcx
       dpqn2_dpqm3 = tcy

c            Report des calculs precedents dans le calcul de la derivee de flu2.

c     ! qp : left state
       dpflu2_dpqp1 = dpflu2_dpp1 * dpp1_dpqp1 + dpflu2_dprou1 *
     &                dprou1_dpqp1
       dpflu2_dpqp2 = dpflu2_dpqn1 * dpqn1_dpqp2 + dpflu2_dprou1 *
     &                dprou1_dpqp2
       dpflu2_dpqp3 = dpflu2_dpqn1 * dpqn1_dpqp3
       dpflu2_dpqp5 = dpflu2_dpp1 * dpp1_dpqp5

c     ! qm : right state
       dpflu2_dpqm1 = dpflu2_dpp2 * dpp2_dpqm1 + dpflu2_dprou2 *
     &                dprou2_dpqm1
       dpflu2_dpqm2 = dpflu2_dpqn2 * dpqn2_dpqm2 + dpflu2_dprou2 *
     &                dprou2_dpqm2
       dpflu2_dpqm3 = dpflu2_dpqn2 * dpqn2_dpqm3
       dpflu2_dpqm5 = dpflu2_dpp2 * dpp2_dpqm5

c********* ETAPE 4 : Calcul de la derivee de flu2 par rapport a rop

c            Calcul des derivees partielles intermediaires necessaires
c            au calcul de la derivee de flu2.

       dpson_dpc = -0.5 * sqrt(qn1**2)*(c**(-3/2))
       dpflu2_dpc = ( dpflu2_dptdu * dptdu_dpu + dpflu2_dpu ) *  
     &              dpu_dptam1 * dptam1_dptam * dptam_dpson * 
     &              dpson_dpc
       dpc_dprop_l5 = rgp * gam

c            Report des calculs precedents dans le calcul de la derivee de flu2.

       dpflu2_dprop1_np = dpflu2_dpqm1 * c6
       dpflu2_dprop1_l = dpflu2_dpqp1 * c5 + dpflu2_dpqm1 * c4
       dpflu2_dprop1_nm = dpflu2_dpqp1 * c4 + dpflu2_dpqm1 * c5
       dpflu2_dprop1_nm2 = dpflu2_dpqp1 * c6

       dpflu2_dprop2_np = dpflu2_dpqm2 * c6
       dpflu2_dprop2_l = dpflu2_dpqp2 * c5 + dpflu2_dpqm2 * c4
       dpflu2_dprop2_nm = dpflu2_dpqp2 * c4 + dpflu2_dpqm2 * c5
       dpflu2_dprop2_nm2 = dpflu2_dpqp2 * c6

       dpflu2_dprop3_np = dpflu2_dpqm3 * c6
       dpflu2_dprop3_l = dpflu2_dpqp3 * c5 + dpflu2_dpqm3 * c4
       dpflu2_dprop3_nm = dpflu2_dpqp3 * c4 + dpflu2_dpqm3 * c5
       dpflu2_dprop3_nm2 = dpflu2_dpqp3 * c6

       dpflu2_dprop5_np = dpflu2_dpqm5 * c6
       dpflu2_dprop5_l = dpflu2_dpqp5 * c5 + dpflu2_dpqm5 * c4 + 
     &                          dpflu2_dpc * dpc_dprop_l5 ! Prise en compte de c dans dpflu2_dprop5_l
       dpflu2_dprop5_nm = dpflu2_dpqp5 * c4 + dpflu2_dpqm5 * c5
       dpflu2_dprop5_nm2 = dpflu2_dpqp5 * c6

c ---------------- Pour VERIFICATION (comparaison avec les differences finies calculees).

      if(j.eq.1.and.l-lij+ind_loop(1)-1.eq.100) then


c------ Affichage des resultats pour les points (i=84,j=1) et (i=100,j=1)

       write(*,*) " AV. Point (i= ", l-lij+ind_loop(1)-1,"j= 1)"
       write(*,*) ""

       write(*,*) " dpflu2_dprop(l+v1) ",
     &            dpflu2_dprop1_l, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " dpflu2_dprop(l+v2) ",
     &            dpflu2_dprop2_l, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " dpflu2_dprop(l+v3) ",
     &            dpflu2_dprop3_l, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " dpflu2_dprop(l+v5) ",
     &            dpflu2_dprop5_l, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) "-----------------------"

      endif

c******************************************************************************
!------------ Calcul des composants de dpRopl_dpRopldjr (dpPg_dpP)
c******************************************************************************

!      NORMALISATION de tcx, tcy et tcz
   !    ci_mtr = 1. ! Valeur a integrer depuis BCWallInviscid.for 
   !    cj_mtr = 1. ! Valeur a integrer depuis BCWallInviscid.for
   !    ck_mtr = 0. ! Valeur a integrer depuis BCWallInviscid.for

       idir = param_int_eff(EFF_IDIR)
       neq_mtr = param_int(NEQ_IJ)

      call shape_tab_mtr(neq_mtr, param_int, idir,
     &                   ic,jc,kc,kc_vent,
     &                   ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale)

       ncx = tcx*ci_mtr
       ncy = tcy*cj_mtr
       ncz = tcz*ck_mtr

       surf_nc = sqrt(ncx*ncx + ncy*ncy + ncz*ncz)
       surf_nc = max(surf_nc,1e-30)

       ncx = ncx/surf_nc
       ncy = ncy/surf_nc
       ncz = ncz/surf_nc


!      Derivees calculees a partir des fonctions se trouvant dans : .FastS/BC/BCWallInviscid.for
       dpRop1_b_dpRop1_h = 1.

       dpRop2_b_dpRop2_h = 1. - 2.*(ncx**2)
       dpRop2_b_dpRop3_h = - 2.*ncy*ncx
       dpRop2_b_dpRop4_h = - 2.*ncz*ncx

       dpRop3_b_dpRop2_h = - 2.*ncx*ncy
       dpRop3_b_dpRop3_h = 1. - 2.*(ncy**2)
       dpRop3_b_dpRop4_h = - 2.*ncz*ncy

       dpRop4_b_dpRop2_h = - 2.*ncx*ncz
       dpRop4_b_dpRop3_h = - 2.*ncy*ncz
       dpRop4_b_dpRop4_h = 1. - 2.*(ncz**2)

       dpRop5_b_dpRop5_h = 1.


c******************************************************************************
!------------ Integration de dpPg_dpP dans dpflu2_dprop
c******************************************************************************

       dpflu2_dprop1_np = dpflu2_dprop1_np +
     &                    dpflu2_dprop1_nm2 * dpRop1_b_dpRop1_h
       dpflu2_dprop1_l = dpflu2_dprop1_l +
     &                    dpflu2_dprop1_nm * dpRop1_b_dpRop1_h

       dpflu2_dprop2_np = dpflu2_dprop2_np +
     &                    dpflu2_dprop2_nm2 * dpRop2_b_dpRop2_h +
     &                    dpflu2_dprop3_nm2 * dpRop3_b_dpRop2_h
       dpflu2_dprop2_l = dpflu2_dprop2_l +
     &                   dpflu2_dprop2_nm * dpRop2_b_dpRop2_h +
     &                   dpflu2_dprop3_nm * dpRop3_b_dpRop2_h

       dpflu2_dprop3_np = dpflu2_dprop3_np +
     &                    dpflu2_dprop2_nm2 * dpRop2_b_dpRop3_h +
     &                    dpflu2_dprop3_nm2 * dpRop3_b_dpRop3_h
       dpflu2_dprop3_l = dpflu2_dprop3_l +
     &                   dpflu2_dprop2_nm * dpRop2_b_dpRop3_h +
     &                   dpflu2_dprop3_nm * dpRop3_b_dpRop3_h

       dpflu2_dprop5_np = dpflu2_dprop5_np +
     &                    dpflu2_dprop5_nm2 * dpRop5_b_dpRop5_h
       dpflu2_dprop5_l = dpflu2_dprop5_l +
     &                   dpflu2_dprop5_nm * dpRop5_b_dpRop5_h

c ---------------- Pour VERIFICATION (comparaison avec les differences finies calculees).

      if(j.eq.1.and.l-lij+ind_loop(1)-1.eq.100) then


c------ Affichage des resultats pour les points (i=84,j=1) et (i=100,j=1)

       write(*,*) " AP. Point (i= ", l-lij+ind_loop(1)-1,"j= 1)"
       write(*,*) ""

       write(*,*) " dpflu2_dprop(l+v1) ",
     &            dpflu2_dprop1_l, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " dpflu2_dprop(l+v2) ",
     &            dpflu2_dprop2_l, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " dpflu2_dprop(l+v3) ",
     &            dpflu2_dprop3_l, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " dpflu2_dprop(l+v5) ",
     &            dpflu2_dprop5_l, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) "-----------------------"

      endif


c****************************************************************************************************
c*********************************** FLU3 ***********************************************************

c    Equation de flu3 :
c       flu3 = u*(rov1+rov2) - tdu*(rov2-rov1) + tcx*p1p2

c********* ETAPE 1 : Calcul des derivees de flu3 par rapport a u, tdu et p1p2

       dpflu3_dpu = rov1 + rov2
       dpflu3_dptdu = -(rov2 - rov1)
       dpflu3_dpp1p2 = tcy


c********* ETAPE 2 : Calcul de la derivee de flu3 par rapport a
c                    r1,r2, rou1,rou2, rov1,rov2, p1,p2, qn1 et qn2.
c            En fait : Liste des variables pour flu3 :
c                                                     rov1,rov2, p1,p2, qn1 et qn2.

c            Toutes les derivees partielles intermediaires pour le calcul
c            de la derivee de flu3 sont identiques a celles de flu2.

c            Report des derivees partielles intermediaires dans le calcul de la derivee de flu3.

       dpflu3_dpqn1 = ( dpflu3_dptdu * dptdu_dpu + dpflu3_dpu ) * 
     &                ( dpu_dpqn1 + dpu_dptam1 * dptam1_dptam *
     &                dptam_dpson * dpson_dpqn1 )

       dpflu3_dpqn2 = ( dpflu3_dptdu * dptdu_dpu + dpflu3_dpu ) *
     &                dpu_dpqn2

       dpflu3_dpp1 = dpflu3_dpp1p2 * dpp1p2_dpp1 + ( dpflu3_dptdu *
     &               dptdu_dpu + dpflu3_dpu ) * dpu_dpp1

       dpflu3_dpp2 = dpflu3_dpp1p2 * dpp1p2_dpp2 + ( dpflu3_dptdu *
     &               dptdu_dpu + dpflu3_dpu ) * dpu_dpp2

       dpflu3_dprov1 = u+tdu
       dpflu3_dprov2 = u-tdu

       !  D'apres l'equation de flu3
       !     dpflu3_dpr1 = 0
       !     dpflu3_dpr2 = 0
       !     dpflu3_dprou1 = 0
       !     dpflu3_dprou2 = 0


c********* ETAPE 3 : Calcul des dérivées de flu2 par rapport à qp et qm

c            Calcul des derivees partielles intermediaires necessaires au calcul de la derivee de flu3.

c            Toutes les derivees partielles intermediaires pour le calcul de la derivee de flu3
c            sont identiques a celles de flu2, sauf pour les quatres suivantes :

c     ! qp : left state
       dprov1_dpqp1 = qp3
       dprov1_dpqp3 = r1

c     ! qm : right state
       dprov2_dpqm1 = qm3
       dprov2_dpqm3 = r2


c            Report des derivees partielles intermediaires dans le calcul de la derivee de flu3.

c     ! qp : left state
       dpflu3_dpqp1 = dpflu3_dpp1 * dpp1_dpqp1 + dpflu3_dprov1 *
     &                dprov1_dpqp1
       dpflu3_dpqp2 = dpflu3_dpqn1 * dpqn1_dpqp2
       dpflu3_dpqp3 = dpflu3_dpqn1 * dpqn1_dpqp3 + dpflu3_dprov1 *
     &                dprov1_dpqp3
       dpflu3_dpqp5 = dpflu3_dpp1 * dpp1_dpqp5

c    ! qm : right state
       dpflu3_dpqm1 = dpflu3_dpp2 * dpp2_dpqm1 + dpflu3_dprov2 *
     &                dprov2_dpqm1
       dpflu3_dpqm2 = dpflu3_dpqn2 * dpqn2_dpqm2
       dpflu3_dpqm3 = dpflu3_dpqn2 * dpqn2_dpqm3 + dpflu3_dprov2 *
     &                dprov2_dpqm3
       dpflu3_dpqm5 = dpflu3_dpp2 * dpp2_dpqm5


c********* ETAPE 4 : Calcul de la derivee de flu3 par rapport a rop

c            Toutes les derivees partielles intermediaires pour le calcul de la derivee de flu3
c            sont identiques a celles de flu2, sauf pour la suivante :

       dpflu3_dpc = ( dpflu3_dptdu * dptdu_dpu + dpflu3_dpu ) *  
     &              dpu_dptam1 * dptam1_dptam * dptam_dpson * 
     &              dpson_dpc


c            Report des derivees partielles intermediaires dans le calcul de la derivee de flu3.

       dpflu3_dprop1_np = dpflu3_dpqm1 * c6
       dpflu3_dprop1_l = dpflu3_dpqp1 * c5 + dpflu3_dpqm1 * c4
       dpflu3_dprop1_nm = dpflu3_dpqp1 * c4 + dpflu3_dpqm1 * c5
       dpflu3_dprop1_nm2 = dpflu3_dpqp1 * c6

       dpflu3_dprop2_np = dpflu3_dpqm2 * c6
       dpflu3_dprop2_l = dpflu3_dpqp2 * c5 + dpflu3_dpqm2 * c4
       dpflu3_dprop2_nm = dpflu3_dpqp2 * c4 + dpflu3_dpqm2 * c5
       dpflu3_dprop2_nm2 = dpflu3_dpqp2 * c6

       dpflu3_dprop3_np = dpflu3_dpqm3 * c6
       dpflu3_dprop3_l = dpflu3_dpqp3 * c5 + dpflu3_dpqm3 * c4
       dpflu3_dprop3_nm = dpflu3_dpqp3 * c4 + dpflu3_dpqm3 * c5
       dpflu3_dprop3_nm2 = dpflu3_dpqp3 * c6

       dpflu3_dprop5_np = dpflu3_dpqm5 * c6
       dpflu3_dprop5_l = dpflu3_dpqp5 * c5 + dpflu3_dpqm5 * c4 + 
     &                          dpflu3_dpc * dpc_dprop_l5 ! Prise en compte de c dans dpflu3_dprop5_l
       dpflu3_dprop5_nm = dpflu3_dpqp5 * c4 + dpflu3_dpqm5 * c5
       dpflu3_dprop5_nm2 = dpflu3_dpqp5 * c6


c ---------------- Pour VERIFICATION (comparaison avec les differences finies calculees).

      if(j.eq.1.and.l-lij+ind_loop(1)-1.eq.100) then


c------ Affichage des resultats pour les points (i=84,j=1) et (i=100,j=1)

       write(*,*) " AV. Point (i= ", l-lij+ind_loop(1)-1,"j= 1)"
       write(*,*) ""

       write(*,*) " dpflu3_dprop(l+v1) ",
     &            dpflu3_dprop1_l, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " dpflu3_dprop(l+v2) ",
     &            dpflu3_dprop2_l, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " dpflu3_dprop(l+v3) ",
     &            dpflu3_dprop3_l, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " dpflu3_dprop(l+v5) ",
     &            dpflu3_dprop5_l, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) "-----------------------"

      endif


c******************************************************************************
!------------ Calcul des composants de dpRopl_dpRopldjr (dpPg_dpP)
c******************************************************************************

!      Derivees calculees depuis BCWallInviscid.for sont ici identiques.


c******************************************************************************
!------------ Integration de dpPg_dpP dans dpflu3_dprop
c******************************************************************************

       dpflu3_dprop1_np = dpflu3_dprop1_np +
     &                    dpflu3_dprop1_nm2 * dpRop1_b_dpRop1_h
       dpflu3_dprop1_l = dpflu3_dprop1_l +
     &                   dpflu3_dprop1_nm * dpRop1_b_dpRop1_h

       dpflu3_dprop2_np = dpflu3_dprop2_np +
     &                    dpflu3_dprop2_nm2 * dpRop2_b_dpRop2_h +
     &                    dpflu3_dprop3_nm2 * dpRop3_b_dpRop2_h
       dpflu3_dprop2_l = dpflu3_dprop2_l +
     &                   dpflu3_dprop2_nm * dpRop2_b_dpRop2_h +
     &                   dpflu3_dprop3_nm * dpRop3_b_dpRop2_h

       dpflu3_dprop3_np = dpflu3_dprop3_np +
     &                    dpflu3_dprop2_nm2 * dpRop2_b_dpRop3_h +
     &                    dpflu3_dprop3_nm2 * dpRop3_b_dpRop3_h
       dpflu3_dprop3_l = dpflu3_dprop3_l +
     &                   dpflu3_dprop2_nm * dpRop2_b_dpRop3_h +
     &                   dpflu3_dprop3_nm * dpRop3_b_dpRop3_h

       dpflu3_dprop5_np = dpflu3_dprop5_np +
     &                   dpflu3_dprop5_nm2 * dpRop5_b_dpRop5_h
       dpflu3_dprop5_l = dpflu3_dprop5_l +
     &                   dpflu3_dprop5_nm * dpRop5_b_dpRop5_h


c ---------------- Pour VERIFICATION (comparaison avec les differences finies calculees).

      if(j.eq.1.and.l-lij+ind_loop(1)-1.eq.100) then


c------ Affichage des resultats pour les points (i=84,j=1) et (i=100,j=1)

       write(*,*) " AP. Point (i= ", l-lij+ind_loop(1)-1,"j= 1)"
       write(*,*) ""

       write(*,*) " dpflu3_dprop(l+v1) ",
     &            dpflu3_dprop1_l, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " dpflu3_dprop(l+v2) ",
     &            dpflu3_dprop2_l, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " dpflu3_dprop(l+v3) ",
     &            dpflu3_dprop3_l, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " dpflu3_dprop(l+v5) ",
     &            dpflu3_dprop5_l, " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) "-----------------------"

      endif


C****************************************************************************************
c***********         FIN DE CALCUL DE dpflu_dprop (dpeff_dprop)          ****************
C****************************************************************************************


C****************************************************************************************
c***********           CALCUL DE dpClp_dprop et dpCdp_dprop              ****************
C****************************************************************************************

      cc=sens*surfinv*cosAoA
      cs=sens*surfinv*sinAoA

c -------- Calcul de dpCdp_dprop

!        drag = (eff[0]*cosAoA+eff[1]*sinAoA)/(corde*span)

       dpCdp_dprop1_np = cc*dpflu2_dprop1_np +
     &                   cs*dpflu3_dprop1_np
       dpCdp_dprop1_l  = cc*dpflu2_dprop1_l +
     &                   cs*dpflu3_dprop1_l


       dpCdp_dprop2_np = cc*dpflu2_dprop2_np +
     &                   cs*dpflu3_dprop2_np
       dpCdp_dprop2_l  = cc*dpflu2_dprop2_l +
     &                   cs*dpflu3_dprop2_l

       dpCdp_dprop3_np = cc*dpflu2_dprop3_np +
     &                   cs*dpflu3_dprop3_np
       dpCdp_dprop3_l  = cc*dpflu2_dprop3_l +
     &                   cs*dpflu3_dprop3_l

       dpCdp_dprop5_np = cc*dpflu2_dprop5_np +
     &                   cs*dpflu3_dprop5_np
       dpCdp_dprop5_l  = cc*dpflu2_dprop5_l +
     &                   cs*dpflu3_dprop5_l


c -------- Calcul de dpClp_dprop

!        lift = (eff[1]*cosAoA-eff[0]*sinAoA)/(corde*span)

       dpClp_dprop1_np = cc*dpflu3_dprop1_np -
     &                   cs*dpflu2_dprop1_np
       dpClp_dprop1_l  = cc*dpflu3_dprop1_l  -
     &                   cs*dpflu2_dprop1_l

       dpClp_dprop2_np = cc*dpflu3_dprop2_np -
     &                   cs*dpflu2_dprop2_np
       dpClp_dprop2_l  = cc*dpflu3_dprop2_l  -
     &                   cs*dpflu2_dprop2_l

       dpClp_dprop3_np = cc*dpflu3_dprop3_np -
     &                   cs*dpflu2_dprop3_np
       dpClp_dprop3_l  = cc*dpflu3_dprop3_l -
     &                   cs*dpflu2_dprop3_l

       dpClp_dprop5_np = cc*dpflu3_dprop5_np -
     &                   cs*dpflu2_dprop5_np
       dpClp_dprop5_l  = cc*dpflu3_dprop5_l -
     &                   cs*dpflu2_dprop5_l

C****************************************************************************************
c***********        FIN DE CALCUL DE dpClp_dprop et dpCdp_dprop          ****************
C****************************************************************************************

#include "FastS/ADJOINT/AUSM/2d/dpP_dpW_fluausm_euler_2d_idir_min.for" 

C****************************************************************************************
c***********              CALCUL DE dpClp_dpW et dpCdp_dpW               ****************
C****************************************************************************************
c
c -------- Calcul de dpCdp_dpW
       dpCdp_dpW(np+v1) = dpCdp_dprop1_np * dprop1_dpW1_np + 
     &                    dpCdp_dprop2_np * dprop2_dpW1_np +
     &                    dpCdp_dprop3_np * dprop3_dpW1_np +
     &                    dpCdp_dprop5_np * dprop5_dpW1_np
       dpCdp_dpW(l+v1)  = dpCdp_dprop1_l * dprop1_dpW1_l + 
     &                    dpCdp_dprop2_l * dprop2_dpW1_l +
     &                    dpCdp_dprop3_l * dprop3_dpW1_l +
     &                    dpCdp_dprop5_l * dprop5_dpW1_l

       dpCdp_dpW(np+v2) = dpCdp_dprop2_np * dprop2_dpW2_np +
     &                    dpCdp_dprop5_np * dprop5_dpW2_np
       dpCdp_dpW(l+v2)  = dpCdp_dprop2_l * dprop2_dpW2_l +
     &                    dpCdp_dprop5_l * dprop5_dpW2_l

       dpCdp_dpW(np+v3) = dpCdp_dprop3_np * dprop3_dpW3_np +
     &                    dpCdp_dprop5_np * dprop5_dpW3_np
       dpCdp_dpW(l+v3)  = dpCdp_dprop3_l * dprop3_dpW3_l +
     &                    dpCdp_dprop5_l * dprop5_dpW3_l

       dpCdp_dpW(np+v5) = dpCdp_dprop5_np * dprop5_dpW5_np
       dpCdp_dpW(l+v5)  = dpCdp_dprop5_l * dprop5_dpW5_l


c -------- Calcul de dpClp_dpW
       dpClp_dpW(np+v1) = dpClp_dprop1_np * dprop1_dpW1_np + 
     &                    dpClp_dprop2_np * dprop2_dpW1_np +
     &                    dpClp_dprop3_np * dprop3_dpW1_np +
     &                    dpClp_dprop5_np * dprop5_dpW1_np
       dpClp_dpW(l+v1)  = dpClp_dprop1_l * dprop1_dpW1_l + 
     &                    dpClp_dprop2_l * dprop2_dpW1_l +
     &                    dpClp_dprop3_l * dprop3_dpW1_l +
     &                    dpClp_dprop5_l * dprop5_dpW1_l

       dpClp_dpW(np+v2) = dpClp_dprop2_np * dprop2_dpW2_np +
     &                    dpClp_dprop5_np * dprop5_dpW2_np
       dpClp_dpW(l+v2)  = dpClp_dprop2_l * dprop2_dpW2_l +
     &                    dpClp_dprop5_l * dprop5_dpW2_l

       dpClp_dpW(np+v3) = dpClp_dprop3_np * dprop3_dpW3_np +
     &                    dpClp_dprop5_np * dprop5_dpW3_np
       dpClp_dpW(l+v3)  = dpClp_dprop3_l * dprop3_dpW3_l +
     &                    dpClp_dprop5_l * dprop5_dpW3_l

       dpClp_dpW(np+v5) = dpClp_dprop5_np * dprop5_dpW5_np
       dpClp_dpW(l+v5)  = dpClp_dprop5_l * dprop5_dpW5_l


c ---------------- Pour VERIFICATION (comparaison avec les differences finies calculees)
c
      if(j.eq.1.and.l-lij+ind_loop(1)-1.eq.100) then

c------ Affichage des resultats pour les points (i=84,j=1) et (i=100,j=1)

       write(*,*) " Point (i= ", l-lij+ind_loop(1)-1,"j= 1)"
       write(*,*) ""
       write(*,*) " dpCdp_dpW(l+v1) ",
     &            dpCdp_dpW(l+v1), " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " dpClp_dpW(l+v1) ",
     &            dpClp_dpW(l+v1), " i ",
     & l-lij+ind_loop(1)-1
       write(*,*) ""
       write(*,*) " dpCdp_dpW(l+v2) ",
     &            dpCdp_dpW(l+v2), " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " dpClp_dpW(l+v2) ",
     &            dpClp_dpW(l+v2), " i ",
     & l-lij+ind_loop(1)-1
       write(*,*) ""
       write(*,*) " dpCdp_dpW(l+v3) ",
     &            dpCdp_dpW(l+v3), " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " dpClp_dpW(l+v3) ",
     &            dpClp_dpW(l+v3), " i ",
     & l-lij+ind_loop(1)-1
       write(*,*) ""
       write(*,*) " dpCdp_dpW(l+v5) ",
     &            dpCdp_dpW(l+v5), " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " dpClp_dpW(l+v5) ",
     &            dpClp_dpW(l+v5), " i ",
     & l-lij+ind_loop(1)-1
       write(*,*) ""
       write(*,*) "-----------------------"

c------ Affichage des resultats pour les points (i=84,j=2) et (i=100,j=2)
       write(*,*) " Point (i= ", l-lij+ind_loop(1)-1,"j= 2)"
       write(*,*) ""
       write(*,*) " dpCdp_dpW(np+v1) ",
     &            dpCdp_dpW(np+v1), " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " dpClp_dpW(np+v1) ",
     &            dpClp_dpW(np+v1), " i ",
     & l-lij+ind_loop(1)-1
       write(*,*) ""
       write(*,*) " dpCdp_dpW(np+v2) ",
     &            dpCdp_dpW(np+v2), " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " dpClp_dpW(np+v2) ",
     &            dpClp_dpW(np+v2), " i ",
     & l-lij+ind_loop(1)-1
       write(*,*) ""
       write(*,*) " dpCdp_dpW(np+v3) ",
     &            dpCdp_dpW(np+v3), " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " dpClp_dpW(np+v3) ",
     &            dpClp_dpW(np+v3), " i ",
     & l-lij+ind_loop(1)-1
       write(*,*) ""
       write(*,*) " dpCdp_dpW(np+v5) ",
     &            dpCdp_dpW(np+v5), " i ",
     & l-lij+ind_loop(1)-1

       write(*,*) " dpClp_dpW(np+v5) ",
     &            dpClp_dpW(np+v5), " i ",
     & l-lij+ind_loop(1)-1
       write(*,*) ""
       write(*,*) "-----------------------"

      endif


C****************************************************************************************
c***********          FIN DE CALCUL DE dpClp_dpW et dpCdp_dpW            ****************
C****************************************************************************************

                    enddo
                 ENDDO
              ENDDO


      ELSEIF(param_int_eff(35).eq.4) THEN

              inc_x1 = 1                             !(i+1,j,k  )
              inc_x2 = param_int(10)*param_int(10+1) !(i  ,j,k+1)
              inc_x3 = inc_x1 + inc_x2
              iter = 0

          write(*,*) 'ind_loop(2) ',ind_loop(2)
          write(*,*) 'ind_loop(1) ',ind_loop(1)
          write(*,*) '-------------------------'
          write(*,*) 'ind_loop(4) ',ind_loop(4)
          write(*,*) 'ind_loop(3) ',ind_loop(3)
          write(*,*) '-------------------------'

              DO k = ind_loop(5), ind_loop(6)
                 DO j = ind_loop(3), ind_loop(4)

	      lij  =       inddm( ind_loop(1) , j, k) -1
	      ltij = lij - indmtr(ind_loop(1) , j, k) +1
	      lfij = lij - indflu(ind_loop(1) , j, k) +1
	      lxij = lij - indcg( ind_loop(1) , j, k) +1

!CC    !DIR$ ASSUME (mod(lij,4) .eq. 0)
!DIR$ IVDEP

	      do l = lij+1, lij+1 + ind_loop(2) - ind_loop(1)

                       iter = iter + 1

	               lt = l  - ltij
	               lvo= lt
	               lf = l  - lfij
	               lx = l  - lxij


!c...............Metrique
	               tcx = tj(lt +v1mtr)
	               tcy = tj(lt +v2mtr)
	               sj  = sqrt (tcx*tcx + tcy*tcy)

	               nm  = l -  incj
	               nm2 = l -2*incj
	               np  = l +  incj

       if(l-lij+ind_loop(1)-1.eq.29) then
        write(*,*) 'tcx ', tcx
        write(*,*) 'tcy ', tcy
      endif

! Pente (qm) a l'interface droite et  (qp) a l'interface gauche

c    ! qm: right state,  qp: left state
        qm1=  c4*rop(l +v1) +  c5*rop(nm  +v1) +  c6*rop(np +v1)
        qp1=  c4*rop(nm+v1) +  c6*rop(nm2 +v1) +  c5*rop(l  +v1)

c    ! qm: right state,  qp: left state
        qm2=  c4*rop(l +v2) +  c5*rop(nm  +v2) +  c6*rop(np +v2)
        qp2=  c4*rop(nm+v2) +  c6*rop(nm2 +v2) +  c5*rop(l  +v2)

c    ! qm: right state,  qp: left state
        qm3=  c4*rop(l +v3) +  c5*rop(nm  +v3) +  c6*rop(np +v3)
        qp3=  c4*rop(nm+v3) +  c6*rop(nm2 +v3) +  c5*rop(l  +v3)

c    ! qm: right state,  qp: left state
        qm5=  c4*rop(l +v5) +  c5*rop(nm  +v5) +  c6*rop(np +v5)
        qp5=  c4*rop(nm+v5) +  c6*rop(nm2 +v5) +  c5*rop(l  +v5)

! Determination etat gauche (rou1) et droit (rou2): ro, roui, roe+p
        r1   = qp1
        rou1 = r1*qp2
        rov1 = r1*qp3
        p1   = r1*qp5*rgp
        h1   = gam1*p1 + .5*(rou1*qp2+rov1*qp3)

! Determination etat droite: ro, roui, roe+p
        r2  =qm1
        rou2=r2*qm2
        rov2=r2*qm3
        p2  =r2*qm5*rgp
        h2  =gam1*p2 + .5*(rou2*qm2+rov2*qm3)

! Determination vitesse normale interface
        qn1=qp2*tcx+qp3*tcy
        qn2=qm2*tcx+qm3*tcy

! Modification de vitesse normale par ajout
! de stabilisation de type Rhie-Chow
        c   = rgp*gam*rop(l +v5)  !c^2
        son = sqrt(qn1*qn1 / c)
        tam = c3*son+sj
        tam1= max(0.,tam)*c2 ! fct amortissement: c3*Mach+1
        u   = 0.25*(qn1+qn2)-tam1*(p2-p1)
        tdu = max(abs(u),c1*sj)

        p1p2= (p1+p2)*0.5


C****************************************************************************************
c************              CALCUL DE dpflu_dprop (dpeff_dprop)           ****************
c****************    LE CALCUL SE FAIT DE MANIERE ASCENDANTE    *************************
C****************************************************************************************

c*********************************** FLU2 ********************************************************

c    Equation de flu2 :
c       flu2 = u*(rou1+rou2) - tdu*(rou2-rou1) + tcx*p1p2

c************ ETAPE 1 : Calcul des derivees de flu2 par rapport a u, tdu et p1p2

        dpflu2_dpu = rou1 + rou2
        dpflu2_dptdu = -(rou2 - rou1)
        dpflu2_dpp1p2 = tcx


c************ ETAPE 2 : Calcul de la derivee de flu2 par rapport a 
c                       r1,r2, rou1,rou2, rov1,rov2, p1,p2, qn1 et qn2.
c            En fait : Liste des variables pour flu2 : 
c                                                     ro1,rou2,p1,p2,qn1,qn2

c            Calcul des derivees partielles intermediaires necessaires
c            au calcul de la derivee de flu2.

        IF (ABS(u).GT.c1*sj) THEN
           dptdu_dpu = sign(1.,u)
        ELSE
           dptdu_dpu = 0.
        END IF

        dpu_dpqn1 = 0.25
        dpu_dptam1 = -(p2-p1)
        dptam1_dptam = 0.5*(1.+sign(1.,tam))*c2
        dptam_dpson = c3
        dpson_dpqn1 = sign(1./sqrt(c),qn1)
        dpu_dpqn2 = 0.25
        dpp1p2_dpp1 = 0.5
        dpu_dpp1 = tam1
        dpp1p2_dpp2 = 0.5
        dpu_dpp2 = -tam1


c            Report des calculs precedents dans le calcul de la derivee de flu2.

        dpflu2_dpqn1 = ( dpflu2_dptdu * dptdu_dpu + dpflu2_dpu ) *
     &                 ( dpu_dpqn1 + dpu_dptam1 * dptam1_dptam * 
     &                 dptam_dpson * dpson_dpqn1 )

        dpflu2_dpqn2 = ( dpflu2_dptdu * dptdu_dpu + dpflu2_dpu ) *
     &                 dpu_dpqn2

        dpflu2_dpp1 = dpflu2_dpp1p2 * dpp1p2_dpp1 + ( dpflu2_dptdu *
     &                dptdu_dpu + dpflu2_dpu ) * dpu_dpp1

       dpflu2_dpp2 = dpflu2_dpp1p2 * dpp1p2_dpp2 + ( dpflu2_dptdu *
     &               dptdu_dpu + dpflu2_dpu ) * dpu_dpp2

       dpflu2_dprou1 = u+tdu
       dpflu2_dprou2 = u-tdu

       !  D'apres l'equation de flu2
       !     dpflu2_dpr1 = 0
       !     dpflu2_dpr2 = 0
       !     dpflu2_dprov1 = 0
       !     dpflu2_dprov2 = 0

c********* ETAPE 3 : Calcul des derivees de flu2 par rapport a
c                    qp1 qp2, qp3, qp5 et qm1, qm2,qm3, qm5
c  
c            Calcul des derivees partielles intermediaires necessaires
c            au calcul de la derivee de flu2.

c     ! qp : left state
       dpp1_dpqp1 = qp5*rgp
       dpp1_dpqp5 = r1*rgp
       dprou1_dpqp1 = qp2
       dprou1_dpqp2 = r1
       dpqn1_dpqp2 = tcx
       dpqn1_dpqp3 = tcy

c     ! qm : right state
       dpp2_dpqm1 = qm5*rgp
       dpp2_dpqm5 = r2*rgp
       dprou2_dpqm1 = qm2
       dprou2_dpqm2 = r2
       dpqn2_dpqm2 = tcx
       dpqn2_dpqm3 = tcy

c            Report des calculs precedents dans le calcul de la derivee de flu2.

c     ! qp : left state
       dpflu2_dpqp1 = dpflu2_dpp1 * dpp1_dpqp1 + dpflu2_dprou1 *
     &                dprou1_dpqp1
       dpflu2_dpqp2 = dpflu2_dpqn1 * dpqn1_dpqp2 + dpflu2_dprou1 *
     &                dprou1_dpqp2
       dpflu2_dpqp3 = dpflu2_dpqn1 * dpqn1_dpqp3
       dpflu2_dpqp5 = dpflu2_dpp1 * dpp1_dpqp5

c     ! qm : right state
       dpflu2_dpqm1 = dpflu2_dpp2 * dpp2_dpqm1 + dpflu2_dprou2 *
     &                dprou2_dpqm1
       dpflu2_dpqm2 = dpflu2_dpqn2 * dpqn2_dpqm2 + dpflu2_dprou2 *
     &                dprou2_dpqm2
       dpflu2_dpqm3 = dpflu2_dpqn2 * dpqn2_dpqm3
       dpflu2_dpqm5 = dpflu2_dpp2 * dpp2_dpqm5

c********* ETAPE 4 : Calcul de la derivee de flu2 par rapport a rop

c            Calcul des derivees partielles intermediaires necessaires
c            au calcul de la derivee de flu2.

       dpson_dpc = -0.5 * sqrt(qn1**2)*(c**(-3/2))
       dpflu2_dpc = ( dpflu2_dptdu * dptdu_dpu + dpflu2_dpu ) *  
     &              dpu_dptam1 * dptam1_dptam * dptam_dpson * 
     &              dpson_dpc
       dpc_dprop_l5 = rgp * gam

c            Report des calculs precedents dans le calcul de la derivee de flu2.

       dpflu2_dprop1_np = dpflu2_dpqm1 * c6
       dpflu2_dprop1_l = dpflu2_dpqp1 * c5 + dpflu2_dpqm1 * c4
       dpflu2_dprop1_nm = dpflu2_dpqp1 * c4 + dpflu2_dpqm1 * c5
       dpflu2_dprop1_nm2 = dpflu2_dpqp1 * c6

       dpflu2_dprop2_np = dpflu2_dpqm2 * c6
       dpflu2_dprop2_l = dpflu2_dpqp2 * c5 + dpflu2_dpqm2 * c4
       dpflu2_dprop2_nm = dpflu2_dpqp2 * c4 + dpflu2_dpqm2 * c5
       dpflu2_dprop2_nm2 = dpflu2_dpqp2 * c6

       dpflu2_dprop3_np = dpflu2_dpqm3 * c6
       dpflu2_dprop3_l = dpflu2_dpqp3 * c5 + dpflu2_dpqm3 * c4
       dpflu2_dprop3_nm = dpflu2_dpqp3 * c4 + dpflu2_dpqm3 * c5
       dpflu2_dprop3_nm2 = dpflu2_dpqp3 * c6

       dpflu2_dprop5_np = dpflu2_dpqm5 * c6
       dpflu2_dprop5_l = dpflu2_dpqp5 * c5 + dpflu2_dpqm5 * c4 + 
     &                          dpflu2_dpc * dpc_dprop_l5 ! Prise en compte de c dans dpflu2_dprop5_l
       dpflu2_dprop5_nm = dpflu2_dpqp5 * c4 + dpflu2_dpqm5 * c5
       dpflu2_dprop5_nm2 = dpflu2_dpqp5 * c6


c******************************************************************************
!------------ Calcul des composants de dpRopl_dpRopldjr (dpPg_dpP)
c******************************************************************************

!      NORMALISATION de tcx, tcy et tcz
   !    ci_mtr = 1. ! Valeur a integrer depuis BCWallInviscid.for 
   !    cj_mtr = 1. ! Valeur a integrer depuis BCWallInviscid.for
   !    ck_mtr = 0. ! Valeur a integrer depuis BCWallInviscid.for

       idir = param_int_eff(EFF_IDIR)
       neq_mtr = param_int(NEQ_IJ)

      call shape_tab_mtr(neq_mtr, param_int, idir,
     &                   ic,jc,kc,kc_vent,
     &                   ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale)

       ncx = tcx*ci_mtr
       ncy = tcy*cj_mtr
       ncz = tcz*ck_mtr

       surf_nc = sqrt(ncx*ncx + ncy*ncy + ncz*ncz)
       surf_nc = max(surf_nc,1e-30)

       ncx = ncx/surf_nc
       ncy = ncy/surf_nc
       ncz = ncz/surf_nc


!      Derivees calculees a partir des fonctions se trouvant dans : .FastS/BC/BCWallInviscid.for
       dpRop1_b_dpRop1_h = 1.

       dpRop2_b_dpRop2_h = 1. - 2.*(ncx**2)
       dpRop2_b_dpRop3_h = - 2.*ncy*ncx
       dpRop2_b_dpRop4_h = - 2.*ncz*ncx

       dpRop3_b_dpRop2_h = - 2.*ncx*ncy
       dpRop3_b_dpRop3_h = 1. - 2.*(ncy**2)
       dpRop3_b_dpRop4_h = - 2.*ncz*ncy

       dpRop4_b_dpRop2_h = - 2.*ncx*ncz
       dpRop4_b_dpRop3_h = - 2.*ncy*ncz
       dpRop4_b_dpRop4_h = 1. - 2.*(ncz**2)

       dpRop5_b_dpRop5_h = 1.


c******************************************************************************
!------------ Integration de dpPg_dpP dans dpflu2_dprop
c******************************************************************************

       dpflu2_dprop1_nm2 = dpflu2_dprop1_nm2 +
     &                    dpflu2_dprop1_np * dpRop1_b_dpRop1_h
       dpflu2_dprop1_nm = dpflu2_dprop1_nm +
     &                    dpflu2_dprop1_l * dpRop1_b_dpRop1_h

       dpflu2_dprop2_nm2 = dpflu2_dprop2_nm2 +
     &                    dpflu2_dprop2_np * dpRop2_b_dpRop2_h +
     &                    dpflu2_dprop3_np * dpRop3_b_dpRop2_h
       dpflu2_dprop2_nm = dpflu2_dprop2_nm +
     &                   dpflu2_dprop2_l * dpRop2_b_dpRop2_h +
     &                   dpflu2_dprop3_l * dpRop3_b_dpRop2_h

       dpflu2_dprop3_nm2 = dpflu2_dprop3_nm2 +
     &                    dpflu2_dprop2_np * dpRop2_b_dpRop3_h +
     &                    dpflu2_dprop3_np * dpRop3_b_dpRop3_h
       dpflu2_dprop3_nm = dpflu2_dprop3_nm +
     &                   dpflu2_dprop2_l * dpRop2_b_dpRop3_h +
     &                   dpflu2_dprop3_l * dpRop3_b_dpRop3_h

       dpflu2_dprop5_nm2 = dpflu2_dprop5_nm2 +
     &                    dpflu2_dprop5_np * dpRop5_b_dpRop5_h
       dpflu2_dprop5_nm = dpflu2_dprop5_nm +
     &                   dpflu2_dprop5_l * dpRop5_b_dpRop5_h



c****************************************************************************************************
c*********************************** FLU3 ***********************************************************

c    Equation de flu3 :
c       flu3 = u*(rov1+rov2) - tdu*(rov2-rov1) + tcx*p1p2

c********* ETAPE 1 : Calcul des derivees de flu3 par rapport a u, tdu et p1p2

       dpflu3_dpu = rov1 + rov2
       dpflu3_dptdu = -(rov2 - rov1)
       dpflu3_dpp1p2 = tcy


c********* ETAPE 2 : Calcul de la derivee de flu3 par rapport a
c                    r1,r2, rou1,rou2, rov1,rov2, p1,p2, qn1 et qn2.
c            En fait : Liste des variables pour flu3 :
c                                                     rov1,rov2, p1,p2, qn1 et qn2.

c            Toutes les derivees partielles intermediaires pour le calcul
c            de la derivee de flu3 sont identiques a celles de flu2.

c            Report des derivees partielles intermediaires dans le calcul de la derivee de flu3.

       dpflu3_dpqn1 = ( dpflu3_dptdu * dptdu_dpu + dpflu3_dpu ) * 
     &                ( dpu_dpqn1 + dpu_dptam1 * dptam1_dptam *
     &                dptam_dpson * dpson_dpqn1 )

       dpflu3_dpqn2 = ( dpflu3_dptdu * dptdu_dpu + dpflu3_dpu ) *
     &                dpu_dpqn2

       dpflu3_dpp1 = dpflu3_dpp1p2 * dpp1p2_dpp1 + ( dpflu3_dptdu *
     &               dptdu_dpu + dpflu3_dpu ) * dpu_dpp1

       dpflu3_dpp2 = dpflu3_dpp1p2 * dpp1p2_dpp2 + ( dpflu3_dptdu *
     &               dptdu_dpu + dpflu3_dpu ) * dpu_dpp2

       dpflu3_dprov1 = u+tdu
       dpflu3_dprov2 = u-tdu

       !  D'apres l'equation de flu3
       !     dpflu3_dpr1 = 0
       !     dpflu3_dpr2 = 0
       !     dpflu3_dprou1 = 0
       !     dpflu3_dprou2 = 0


c********* ETAPE 3 : Calcul des dérivées de flu2 par rapport à qp et qm

c            Calcul des derivees partielles intermediaires necessaires au calcul de la derivee de flu3.

c            Toutes les derivees partielles intermediaires pour le calcul de la derivee de flu3
c            sont identiques a celles de flu2, sauf pour les quatres suivantes :

c     ! qp : left state
       dprov1_dpqp1 = qp3
       dprov1_dpqp3 = r1

c     ! qm : right state
       dprov2_dpqm1 = qm3
       dprov2_dpqm3 = r2


c            Report des derivees partielles intermediaires dans le calcul de la derivee de flu3.

c     ! qp : left state
       dpflu3_dpqp1 = dpflu3_dpp1 * dpp1_dpqp1 + dpflu3_dprov1 *
     &                dprov1_dpqp1
       dpflu3_dpqp2 = dpflu3_dpqn1 * dpqn1_dpqp2
       dpflu3_dpqp3 = dpflu3_dpqn1 * dpqn1_dpqp3 + dpflu3_dprov1 *
     &                dprov1_dpqp3
       dpflu3_dpqp5 = dpflu3_dpp1 * dpp1_dpqp5

c    ! qm : right state
       dpflu3_dpqm1 = dpflu3_dpp2 * dpp2_dpqm1 + dpflu3_dprov2 *
     &                dprov2_dpqm1
       dpflu3_dpqm2 = dpflu3_dpqn2 * dpqn2_dpqm2
       dpflu3_dpqm3 = dpflu3_dpqn2 * dpqn2_dpqm3 + dpflu3_dprov2 *
     &                dprov2_dpqm3
       dpflu3_dpqm5 = dpflu3_dpp2 * dpp2_dpqm5


c********* ETAPE 4 : Calcul de la derivee de flu3 par rapport a rop

c            Toutes les derivees partielles intermediaires pour le calcul de la derivee de flu3
c            sont identiques a celles de flu2, sauf pour la suivante :

       dpflu3_dpc = ( dpflu3_dptdu * dptdu_dpu + dpflu3_dpu ) *  
     &              dpu_dptam1 * dptam1_dptam * dptam_dpson * 
     &              dpson_dpc


c            Report des derivees partielles intermediaires dans le calcul de la derivee de flu3.

       dpflu3_dprop1_np = dpflu3_dpqm1 * c6
       dpflu3_dprop1_l = dpflu3_dpqp1 * c5 + dpflu3_dpqm1 * c4
       dpflu3_dprop1_nm = dpflu3_dpqp1 * c4 + dpflu3_dpqm1 * c5
       dpflu3_dprop1_nm2 = dpflu3_dpqp1 * c6

       dpflu3_dprop2_np = dpflu3_dpqm2 * c6
       dpflu3_dprop2_l = dpflu3_dpqp2 * c5 + dpflu3_dpqm2 * c4
       dpflu3_dprop2_nm = dpflu3_dpqp2 * c4 + dpflu3_dpqm2 * c5
       dpflu3_dprop2_nm2 = dpflu3_dpqp2 * c6

       dpflu3_dprop3_np = dpflu3_dpqm3 * c6
       dpflu3_dprop3_l = dpflu3_dpqp3 * c5 + dpflu3_dpqm3 * c4
       dpflu3_dprop3_nm = dpflu3_dpqp3 * c4 + dpflu3_dpqm3 * c5
       dpflu3_dprop3_nm2 = dpflu3_dpqp3 * c6

       dpflu3_dprop5_np = dpflu3_dpqm5 * c6
       dpflu3_dprop5_l = dpflu3_dpqp5 * c5 + dpflu3_dpqm5 * c4 + 
     &                          dpflu3_dpc * dpc_dprop_l5 ! Prise en compte de c dans dpflu3_dprop5_l
       dpflu3_dprop5_nm = dpflu3_dpqp5 * c4 + dpflu3_dpqm5 * c5
       dpflu3_dprop5_nm2 = dpflu3_dpqp5 * c6


c******************************************************************************
!------------ Calcul des composants de dpRopl_dpRopldjr (dpPg_dpP)
c******************************************************************************

!      Derivees calculees depuis BCWallInviscid.for sont ici identiques.


c******************************************************************************
!------------ Integration de dpPg_dpP dans dpflu3_dprop
c******************************************************************************

       dpflu3_dprop1_nm2 = dpflu3_dprop1_nm2 +
     &                    dpflu3_dprop1_np * dpRop1_b_dpRop1_h
       dpflu3_dprop1_nm = dpflu3_dprop1_nm +
     &                   dpflu3_dprop1_l * dpRop1_b_dpRop1_h

       dpflu3_dprop2_nm2 = dpflu3_dprop2_nm2 +
     &                    dpflu3_dprop2_np * dpRop2_b_dpRop2_h +
     &                    dpflu3_dprop3_np * dpRop3_b_dpRop2_h
       dpflu3_dprop2_nm = dpflu3_dprop2_nm +
     &                   dpflu3_dprop2_l * dpRop2_b_dpRop2_h +
     &                   dpflu3_dprop3_l * dpRop3_b_dpRop2_h

       dpflu3_dprop3_nm2 = dpflu3_dprop3_nm2 +
     &                    dpflu3_dprop2_np * dpRop2_b_dpRop3_h +
     &                    dpflu3_dprop3_np * dpRop3_b_dpRop3_h
       dpflu3_dprop3_nm = dpflu3_dprop3_nm +
     &                   dpflu3_dprop2_l * dpRop2_b_dpRop3_h +
     &                   dpflu3_dprop3_l * dpRop3_b_dpRop3_h

       dpflu3_dprop5_nm2 = dpflu3_dprop5_nm2 +
     &                   dpflu3_dprop5_np * dpRop5_b_dpRop5_h
       dpflu3_dprop5_nm = dpflu3_dprop5_nm +
     &                   dpflu3_dprop5_l * dpRop5_b_dpRop5_h


C****************************************************************************************
c***********         FIN DE CALCUL DE dpflu_dprop (dpeff_dprop)          ****************
C****************************************************************************************


C****************************************************************************************
c***********           CALCUL DE dpClp_dprop et dpCdp_dprop              ****************
C****************************************************************************************

      cc=sens*surfinv*cosAoA
      cs=sens*surfinv*sinAoA

c -------- Calcul de dpCdp_dprop

!        drag = (eff[0]*cosAoA+eff[1]*sinAoA)/(corde*span)

       dpCdp_dprop1_nm2 = cc*dpflu2_dprop1_nm2 +
     &                   cs*dpflu3_dprop1_nm2
       dpCdp_dprop1_nm  = cc*dpflu2_dprop1_nm +
     &                   cs*dpflu3_dprop1_nm


       dpCdp_dprop2_nm2 = cc*dpflu2_dprop2_nm2 +
     &                   cs*dpflu3_dprop2_nm2
       dpCdp_dprop2_nm  = cc*dpflu2_dprop2_nm +
     &                   cs*dpflu3_dprop2_nm

       dpCdp_dprop3_nm2 = cc*dpflu2_dprop3_nm2 +
     &                   cs*dpflu3_dprop3_nm2
       dpCdp_dprop3_nm  = cc*dpflu2_dprop3_nm +
     &                   cs*dpflu3_dprop3_nm

       dpCdp_dprop5_nm2 = cc*dpflu2_dprop5_nm2 +
     &                   cs*dpflu3_dprop5_nm2
       dpCdp_dprop5_nm  = cc*dpflu2_dprop5_nm +
     &                   cs*dpflu3_dprop5_nm


c -------- Calcul de dpClp_dprop

!        lift = (eff[1]*cosAoA-eff[0]*sinAoA)/(corde*span)

       dpClp_dprop1_nm2 = cc*dpflu3_dprop1_nm2 -
     &                   cs*dpflu2_dprop1_nm2
       dpClp_dprop1_nm  = cc*dpflu3_dprop1_nm  -
     &                   cs*dpflu2_dprop1_nm

       dpClp_dprop2_nm2 = cc*dpflu3_dprop2_nm2 -
     &                   cs*dpflu2_dprop2_nm2
       dpClp_dprop2_nm  = cc*dpflu3_dprop2_nm  -
     &                   cs*dpflu2_dprop2_nm

       dpClp_dprop3_nm2 = cc*dpflu3_dprop3_nm2 -
     &                   cs*dpflu2_dprop3_nm2
       dpClp_dprop3_nm  = cc*dpflu3_dprop3_nm -
     &                   cs*dpflu2_dprop3_nm

       dpClp_dprop5_nm2 = cc*dpflu3_dprop5_nm2 -
     &                   cs*dpflu2_dprop5_nm2
       dpClp_dprop5_nm  = cc*dpflu3_dprop5_nm -
     &                   cs*dpflu2_dprop5_nm

c ---------------- Pour VERIFICATION (comparaison avec les differences finies calculees).

      if (iter.eq.29) then
c      if(l-lij+ind_loop(1)-1.eq.29) then

c------ Affichage des resultats pour les points (i=84,j=1) et (i=100,j=1)

       write(*,*) " Point (i= ", iter,"j= 128)"
       write(*,*) ""
       write(*,*) " dpCdp_dprop(nm+v1) ",
     &            dpCdp_dprop1_nm, " i ",
     & iter

       write(*,*) " dpClp_dprop(nm+v1) ",
     &            dpClp_dprop1_nm, " i ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dprop(nm+v2) ",
     &            dpCdp_dprop2_nm, " i ",
     & iter

       write(*,*) " dpClp_dprop(nm+v2) ",
     &            dpClp_dprop2_nm, " i ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dprop(nm+v3) ",
     &            dpCdp_dprop3_nm, " i ",
     & iter

       write(*,*) " dpClp_dprop(nm+v3) ",
     &            dpClp_dprop3_nm, " i ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dprop(nm+v5) ",
     &            dpCdp_dprop5_nm, " i ",
     & iter

       write(*,*) " dpClp_dprop(nm+v5) ",
     &            dpClp_dprop5_nm, " i ",
     & iter
       write(*,*) ""
       write(*,*) "-----------------------"

      endif


C****************************************************************************************
c***********        FIN DE CALCUL DE dpClp_dprop et dpCdp_dprop          ****************
C****************************************************************************************

#include "FastS/ADJOINT/AUSM/2d/dpP_dpW_fluausm_euler_2d_idir_max.for" 

C****************************************************************************************
c***********              CALCUL DE dpClp_dpW et dpCdp_dpW               ****************
C****************************************************************************************

c -------- Calcul de dpCdp_dpW
       dpCdp_dpW(nm2+v1) = dpCdp_dprop1_nm2 * dprop1_dpW1_nm2 + 
     &                    dpCdp_dprop2_nm2 * dprop2_dpW1_nm2 +
     &                    dpCdp_dprop3_nm2 * dprop3_dpW1_nm2 +
     &                    dpCdp_dprop5_nm2 * dprop5_dpW1_nm2
       dpCdp_dpW(nm+v1)  = dpCdp_dprop1_nm * dprop1_dpW1_nm + 
     &                    dpCdp_dprop2_nm * dprop2_dpW1_nm +
     &                    dpCdp_dprop3_nm * dprop3_dpW1_nm +
     &                    dpCdp_dprop5_nm * dprop5_dpW1_nm

       dpCdp_dpW(nm2+v2) = dpCdp_dprop2_nm2 * dprop2_dpW2_nm2 +
     &                    dpCdp_dprop5_nm2 * dprop5_dpW2_nm2
       dpCdp_dpW(nm+v2)  = dpCdp_dprop2_nm * dprop2_dpW2_nm +
     &                    dpCdp_dprop5_nm * dprop5_dpW2_nm

       dpCdp_dpW(nm2+v3) = dpCdp_dprop3_nm2 * dprop3_dpW3_nm2 +
     &                    dpCdp_dprop5_nm2 * dprop5_dpW3_nm2
       dpCdp_dpW(nm+v3)  = dpCdp_dprop3_nm * dprop3_dpW3_nm +
     &                    dpCdp_dprop5_nm * dprop5_dpW3_nm

       dpCdp_dpW(nm2+v5) = dpCdp_dprop5_nm2 * dprop5_dpW5_nm2
       dpCdp_dpW(nm+v5)  = dpCdp_dprop5_nm * dprop5_dpW5_nm


c -------- Calcul de dpClp_dpW
       dpClp_dpW(nm2+v1) = dpClp_dprop1_nm2 * dprop1_dpW1_nm2 + 
     &                    dpClp_dprop2_nm2 * dprop2_dpW1_nm2 +
     &                    dpClp_dprop3_nm2 * dprop3_dpW1_nm2 +
     &                    dpClp_dprop5_nm2 * dprop5_dpW1_nm2
       dpClp_dpW(nm+v1)  = dpClp_dprop1_nm * dprop1_dpW1_nm + 
     &                    dpClp_dprop2_nm * dprop2_dpW1_nm +
     &                    dpClp_dprop3_nm * dprop3_dpW1_nm +
     &                    dpClp_dprop5_nm * dprop5_dpW1_nm

       dpClp_dpW(nm2+v2) = dpClp_dprop2_nm2 * dprop2_dpW2_nm2 +
     &                    dpClp_dprop5_nm2 * dprop5_dpW2_nm2
       dpClp_dpW(nm+v2)  = dpClp_dprop2_nm * dprop2_dpW2_nm +
     &                    dpClp_dprop5_nm * dprop5_dpW2_nm

       dpClp_dpW(nm2+v3) = dpClp_dprop3_nm2 * dprop3_dpW3_nm2 +
     &                    dpClp_dprop5_nm2 * dprop5_dpW3_nm2
       dpClp_dpW(nm+v3)  = dpClp_dprop3_nm * dprop3_dpW3_nm +
     &                    dpClp_dprop5_nm * dprop5_dpW3_nm

       dpClp_dpW(nm2+v5) = dpClp_dprop5_nm2 * dprop5_dpW5_nm2
       dpClp_dpW(nm+v5)  = dpClp_dprop5_nm * dprop5_dpW5_nm


c ---------------- Pour VERIFICATION (comparaison avec les differences finies calculees)

      if ((iter.eq.29)) then

      write(*,*) 'iter ', iter
c      if(l-lij+ind_loop(1)-1.eq.127) then

c------ Affichage des resultats pour les points (i=84,j=127) et (i=100,j=127)

       write(*,*) " Point (i= ", iter,"j= 127)"
       write(*,*) ""
       write(*,*) " dpCdp_dpW(nm+v1) ",
     &            dpCdp_dpW(nm+v1), " i ",
     & iter

       write(*,*) " dpClp_dpW(nm+v1) ",
     &            dpClp_dpW(nm+v1), " i ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dpW(nm+v2) ",
     &            dpCdp_dpW(nm+v2), " i ",
     & iter

       write(*,*) " dpClp_dpW(nm+v2) ",
     &            dpClp_dpW(nm+v2), " i ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dpW(nm+v3) ",
     &            dpCdp_dpW(nm+v3), " i ",
     & iter

       write(*,*) " dpClp_dpW(nm+v3) ",
     &            dpClp_dpW(nm+v3), " i ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dpW(nm+v5) ",
     &            dpCdp_dpW(nm+v5), " i ",
     & iter

       write(*,*) " dpClp_dpW(nm+v5) ",
     &            dpClp_dpW(nm+v5), " i ",
     & iter
       write(*,*) ""
       write(*,*) "-----------------------"

c------ Affichage des resultats pour les points (i=84,j=128) et (i=100,j=128)
       write(*,*) " Point (i= ", iter,"j= 128)"
       write(*,*) ""
       write(*,*) " dpCdp_dpW(nm2+v1) ",
     &            dpCdp_dpW(nm2+v1), " i ",
     & iter

       write(*,*) " dpClp_dpW(nm2+v1) ",
     &            dpClp_dpW(nm2+v1), " i ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dpW(nm2+v2) ",
     &            dpCdp_dpW(nm2+v2), " i ",
     & iter

       write(*,*) " dpClp_dpW(nm2+v2) ",
     &            dpClp_dpW(nm2+v2), " i ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dpW(nm2+v3) ",
     &            dpCdp_dpW(nm2+v3), " i ",
     & iter

       write(*,*) " dpClp_dpW(nm2+v3) ",
     &            dpClp_dpW(nm2+v3), " i ",
     & iter
       write(*,*) ""
       write(*,*) " dpCdp_dpW(nm2+v5) ",
     &            dpCdp_dpW(nm2+v5), " i ",
     & iter

       write(*,*) " dpClp_dpW(nm2+v5) ",
     &            dpClp_dpW(nm2+v5), " i ",
     & iter
       write(*,*) ""
       write(*,*) "-----------------------"

      endif


C****************************************************************************************
c***********          FIN DE CALCUL DE dpClp_dpW et dpCdp_dpW            ****************
C****************************************************************************************

                    enddo
                 ENDDO
              ENDDO



!************** FIN DE LA BOUCLE CALCULANT LES EFFORTS SUR LA FACE J ****************
!************************************************************************************



      ELSE

!************************************************************************************
!************************ CALCUL DES EFFORTS SUR LA FACE K **************************
!************************************************************************************

              inc_x1 = 1                                               !(i+1,j  , k)
              inc_x2 = param_int(10)                                   !(i  ,j+1, k)
              inc_x3 = inc_x1 + inc_x2

              DO k = ind_loop(5), ind_loop(6)
                 DO j = ind_loop(3), ind_loop(4)

	      lij  =       inddm( ind_loop(1) , j, k) -1
	      ltij = lij - indmtr(ind_loop(1) , j, k) +1
	      lfij = lij - indflu(ind_loop(1) , j, k) +1
	      lxij = lij - indcg( ind_loop(1) , j, k) +1

!CC    !DIR$ ASSUME (mod(lij,4) .eq. 0)
!DIR$ IVDEP

	      do l = lij+1, lij+1 + ind_loop(2) - ind_loop(1)

	               lt = l  - ltij
	               lvo= lt
	               lf = l  - lfij
	               lx = l  - lxij


       idir = param_int_eff(EFF_IDIR)
       neq_mtr = param_int(NEQ_K)

                    enddo

                 ENDDO
              ENDDO

cc
!******************* FIN DU CALCUL DES EFFORTS SUR LA FACE K ************************
!************************************************************************************
                            
      Endif

      end subroutine dpJ_dpW_fluausm_euler_o3_2d


