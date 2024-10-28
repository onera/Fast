c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 35 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine bvbs_wallmodel(idir,lrhs, nstep,
     &                         neq_mtr, mobile_coef,
     &                         param_int, ind_loop,
     &                         param_real, c4, c5, c6,
     &                         x , y , z,
     &                         ventijk, tijk, rop, xmut, state !insert0
     &                        )
c***********************************************************************
c_U   USER : MARY 
c     ACT
c     Paroi glissante
c     VAL
c_V    Optimisation NEC
c
c     COM
c_C   
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E idir,lrhs,nstep,neq_mtr
      INTEGER_E ind_loop(6), param_int(0:*)

      REAL_E x( param_int(NDIMDX_XYZ) ),y( param_int(NDIMDX_XYZ) ),
     &       z( param_int(NDIMDX_XYZ) )

      REAL_E rop    (param_int(NDIMDX     ), param_int(NEQ)      )
      REAL_E xmut   (param_int(NDIMDX     ), 2                   )
      REAL_E ventijk(param_int(NDIMDX_VENT), param_int(NEQ_VENT) )
      REAL_E tijk   (param_int(NDIMDX_MTR ), neq_mtr             )
      REAL_E state(param_int(NEQ))
      REAL_E mobile_coef, c4,c5,c6 
      REAL_E param_real(0:*)

C...  Valeurs des grandeurs thermo imposees  (insert1)

C Var local 
      INTEGER_E  inc2,inc3,li1,li2,l,iref,jref,kref,lij,lr,lp,sens_int,
     & njf,nkf,ldjr,ic,jc,kc,i,j,k,li3,ldp,kc_vent,ijkplanrec,idirlund,
     & ijk_amor,kshift,jshift, ishift,ci_amor,cj_amor,ck_amor,nijk,
     & iamor,jamor,inci,incj,inck,lmtr,li,l1,l2,l0,nitnwt,newtonmax,
     & ldx,lx,ldnm,ldnp,ldl,ir,jr,kr,lsample,sample,m,iter,lghost,incp,
     & shift_loo,exchange, shift,shiftvent,incijk, sampling,
     & loo,ltg,iadrf, mvent


      REAL_E ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale,tcx,tcy,tcz,
     & ventx,venty,ventz,r,u,v,w,ut,vt,wt,ua,va,wa,s_1,qn
      REAL_E ro,t,nut,c6inv,c0,c1,c2,c3,c7,roe,rue,rve,rwe,ete,
     & qref1,qref2,qref3,qref4,qref5,p,c,snorm, nue, gam1,gam2,
     & gam1_1,gam,cvinv,rho_w,vmax,sens,
     & sni,tn,qen,eigv1,eigv2,eigv3,eigv4,
     & eigv5,qvar1,qvar2,qvar3,qvar4,qvar5,svar1,svar2,svar3,
     & svar4,svar5,rvar1,rvar2,rvar3,rvar4,rvar5,
     & roext,ruext,rvext,rwext,etext,roint,ruint,rvint,rwint,etint,
     & roi,rui,rvi,rwi,eti, tnx,tny, tnz, ri, roinv, sn, roe_inv
      REAL_E pi,ci,roci,dqn,rog,ug,vg,wg,Tg,pg,rgp,roc0 
      REAL_E ui,vi,wi,Ti,ro0,u0,v0,w0,T0,p0,vni,usd0n,usd0n2,
     & wni,wn0,wng,residug,tolnewton,b,f,df,dwng,vng,
     & zlong_transi, freq_transi, lambdaz_transi, ampli_transi,
     & rnd,xlong,zlong,x0,y0,omega,lambda,c_pertu,surf,nx,ny,nz,
     & u_int, un,vn,wn,wall_int, wall_ut,wall_vt,wall_wt, unorm,dist,
     & utau, yplus, aa, g1,g2,g3,tp,fp, umod,pwall,twall,prandlt,cv,
     & cv1,tol,cp,ro_state,u_state,v_state,w_state,t_state,nut_state,
     & expy,nutcible,nu,a,q,ap,bp,dp,delta,delta0,delta1,
     & y1,y2,z1,p1,p2,kappa,superdelta,root,b0,cmus1,temp01,coesut,
     & c0_2,coefa,coefb,coefc,mach,clund, amor, text, uext,tt_ext,
     & cpinv,t1,utau2,tauw,e1,e2,e3,  mu_w


#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"
#include "FastS/formule_vent_param.h"
#include "FastS/formule_xyz_param.h"



c      write(*,*)'idir=', idir,nijk(4),nijk(5),ndimdx
c      write(*,*)'nijk=', nijk
c      write(*,*)'loop vis=', ind_loop


c......determine la forme des tableau metrique en fonction de la nature du domaine
      !Seule la valeur de k_vent et ck_vent a un sens dans cet appel
      call shape_tab_mtr(neq_mtr, param_int, idir,
     &                   ic,jc,kc,kc_vent,
     &                   ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale)

      c_ale = c_ale*mobile_coef
      if(lrhs.eq.1) c_ale = 0.

      snorm =-1.
      sens  =1
      if(mod(idir,2).eq.0) then
       snorm = 1.
       sens  =-1
      endif
      
      gam     = param_real(GAMMA)
      gam1    = param_real(GAMMA) - 1.
      gam1_1  = 1. / gam1
      gam2    = gam / gam1
      rgp     = param_real(CVINF)*gam1
      cv      = param_real( CVINF )
      cvinv   = 1./cv
      cpinv   = 1./(gam*cv)

      prandlt = param_real( PRANDT )

      cmus1  =    param_real( CS )
      temp01 = 1./param_real( TEMP0 )
      coesut =    param_real( XMUL0 ) * (1.+cmus1*temp01)*sqrt(temp01)

      tol   = 1.e-12
      cv1   = 7.1
      kappa = 0.4

      newtonmax    = param_int(NIT_INFLOW)
      tolnewton    = param_real(EPSI_INFLOW)

      ro_state    = state(1)
      u_state     = state(2)/state(1)
      v_state     = state(3)/state(1)
      w_state     = state(4)/state(1)
      t_state     = ( state(5) -0.5*ro_state*( u_state*u_state
     &                                        +v_state*v_state
     &                                        +w_state*w_state) )
     &             /(ro_state*param_real(CVINF))
      if (param_int(NEQ).eq.6) nut_state = state(6)/state(1)

      c0   = 1./c6
      c1   =-(c4 + c5)*c0
      c2   =- c6*c0
      c3   = (2.- c5- c4)*c0
      c7   = (c4-c5)*c0

C...  Parametre specific a la CL(insert2)
#include       "FastS/BC/BCWallModel_init.for"  

      IF (idir.eq.1) THEN 

       iref = 2*ind_loop(2) + 1
       inci = 1
       incp = param_int(IJKV)*inci
       ir   =   ind_loop(2) + 1


       if(param_int(NEQ).eq.5) then

          do k = ind_loop(5), ind_loop(6)
          do j = ind_loop(3), ind_loop(4)

             i    = ind_loop(2)
             l    = inddm( i              , j,  k ) 
             ldp  = indven( ind_loop(2)+1 , j,  k )
             lmtr = indmtr( ind_loop(2)+1 , j,  k )
#include     "FastS/BC/BCWallModel_i.for"
             do i = ind_loop(1), ind_loop(2)-1
               l    = inddm( i              , j,  k ) 
               ldp  = indven( ind_loop(2)+1 , j,  k )  !next rank
               lmtr = indmtr( ind_loop(2)+1 , j,  k )  !next rank
#include       "FastS/BC/BCWallModel_i.for"              
             enddo
          enddo
          enddo

       else

          do k = ind_loop(5), ind_loop(6)
          do j = ind_loop(3), ind_loop(4)

             i    = ind_loop(2)
             l    = inddm( i              , j,  k ) 
             ldp  = indven( ind_loop(2)+1 , j,  k )
             lmtr = indmtr( ind_loop(2)+1 , j,  k )
#include     "FastS/BC/BCWallModelSA_i.for"
             do i = ind_loop(1), ind_loop(2)-1
               l    = inddm( i              , j,  k ) 
               ldp  = indven( ind_loop(2)+1 , j,  k )  !next rank
               lmtr = indmtr( ind_loop(2)+1 , j,  k )  !next rank
#include       "FastS/BC/BCWallModelSA_i.for"            
             enddo
          enddo
          enddo

        endif !param_int(NEQ)

      ELSEIF (idir.eq.2) THEN

       iref = 2*ind_loop(1) - 1
       inci = -1
       incp = param_int(IJKV)*inci
       ir   =   ind_loop(1) - 1

       if(param_int(NEQ).eq.5) then

          do k = ind_loop(5), ind_loop(6)
          do j = ind_loop(3), ind_loop(4)

             i    = ind_loop(1)
             l    = inddm( i            , j, k ) 
             ldp  = indven( ind_loop(1) , j, k )
             lmtr = indmtr( ind_loop(1) , j, k )
#include     "FastS/BC/BCWallModel_i.for"
             do i = ind_loop(1)+1, ind_loop(2)
               l    = inddm( i            , j, k ) 
               ldp  = indven( ind_loop(1) , j, k )  !next rank
               lmtr = indmtr( ind_loop(1) , j, k )  !next rank
#include       "FastS/BC/BCWallModel_i.for"           
             enddo
          enddo
          enddo
       else

          do k = ind_loop(5), ind_loop(6)
          do j = ind_loop(3), ind_loop(4)

             i    = ind_loop(1)
             l    = inddm( i            , j, k ) 
             ldp  = indven( ind_loop(1) , j, k )
             lmtr = indmtr( ind_loop(1) , j, k )
#include     "FastS/BC/BCWallModelSA_i.for"
             do i = ind_loop(1)+1, ind_loop(2)
               l    = inddm( i            , j, k ) 
               ldp  = indven( ind_loop(1) , j, k )  !next rank
               lmtr = indmtr( ind_loop(1) , j, k )  !next rank
#include       "FastS/BC/BCWallModelSA_i.for"         
             enddo
          enddo
          enddo
        endif !param_int(NEQ)


      ELSEIF (idir.eq.3) THEN

       jref = 2*ind_loop(4) + 1
       incj = param_int(NIJK)
       incp = param_int(IJKV+1)*incj
       jr   =   ind_loop(4) +  1

       if(param_int(NEQ).eq.5) then

          do k = ind_loop(5), ind_loop(6)

            j    = ind_loop(4)
!DEC$ IVDEP
              do i = ind_loop(1), ind_loop(2)
                l    = inddm(  i,          j    , k ) 
                ldp  = indven( i, ind_loop(4)+1 , k )
                lmtr = indmtr( i, ind_loop(4)+1 , k )
#include        "FastS/BC/BCWallModel_j.for"
              enddo 
            do j = ind_loop(3), ind_loop(4)-1
!DEC$ IVDEP
              do i = ind_loop(1), ind_loop(2)
                l    = inddm(  i,          j    , k ) 
                ldp  = indven( i, ind_loop(4)+1 , k )  !next rank
                lmtr = indmtr( i, ind_loop(4)+1 , k )  !next rank
#include        "FastS/BC/BCWallModel_j.for"             
              enddo 
            enddo !j
          enddo !k

       else
          do k = ind_loop(5), ind_loop(6)


            j    = ind_loop(4)
!DEC$ IVDEP
              do i = ind_loop(1), ind_loop(2)
                l    = inddm(  i,          j    , k ) 
                ldp  = indven( i, ind_loop(4)+1 , k )
                lmtr = indmtr( i, ind_loop(4)+1 , k )
#include        "FastS/BC/BCWallModelSA_j.for"
              enddo 
            do j = ind_loop(3), ind_loop(4)-1
!DEC$ IVDEP
              do i = ind_loop(1), ind_loop(2)
                l    = inddm(  i,          j    , k ) 
                ldp  = indven( i, ind_loop(4)+1 , k )  !next rank
                lmtr = indmtr( i, ind_loop(4)+1 , k )  !next rank
#include        "FastS/BC/BCWallModelSA_j.for"           
              enddo 
            enddo !j
          enddo !k

       endif !param_int(NEQ)

      ELSEIF (idir.eq.4) THEN

       jref = 2*ind_loop(3) - 1
       incj = -param_int(NIJK)
       incp =  param_int(IJKV+1)*incj
       jr   =   ind_loop(3) -1

       if(param_int(NEQ).eq.5) then

          do k = ind_loop(5), ind_loop(6)

            j    = ind_loop(3)
!DEC$ IVDEP
              do i = ind_loop(1), ind_loop(2)
                l    = inddm(  i,          j    , k ) 
                ldp  = indven( i, ind_loop(3) , k )
                lmtr = indmtr( i, ind_loop(3) , k )
#include        "FastS/BC/BCWallModel_j.for"
              enddo 
            do j = ind_loop(3)+1, ind_loop(4)
!DEC$ IVDEP
              do i = ind_loop(1), ind_loop(2)
                l    = inddm(  i,          j  , k ) 
                ldp  = indven( i, ind_loop(3) , k )  !next rank
                lmtr = indmtr( i, ind_loop(3) , k )  !next rank
#include        "FastS/BC/BCWallModel_j.for"           
              enddo 
            enddo !j
          enddo !k

       else


          do k = ind_loop(5), ind_loop(6)

            j    = ind_loop(3)
!DEC$ IVDEP
              do i = ind_loop(1), ind_loop(2)
                l    = inddm(  i,          j    , k ) 
                ldp  = indven( i, ind_loop(3) , k )
                lmtr = indmtr( i, ind_loop(3) , k )
#include        "FastS/BC/BCWallModelSA_j.for"
              enddo 
            do j = ind_loop(3)+1, ind_loop(4)
!DEC$ IVDEP
              do i = ind_loop(1), ind_loop(2)
                l    = inddm(  i,          j  , k ) 
                ldp  = indven( i, ind_loop(3) , k )  !next rank
                lmtr = indmtr( i, ind_loop(3) , k )  !next rank
#include        "FastS/BC/BCWallModelSA_j.for"         
              enddo 
            enddo !j
          enddo !k
       endif !param_int(NEQ)


      ELSEIF (idir.eq.5) THEN

       kref = 2*ind_loop(6) + 1
       inck = param_int(NIJK)*param_int(NIJK+1)
       incp = param_int(IJKV+2)*inck
       kr   =   ind_loop(6) + 1

       if(param_int(NEQ).eq.5) then

          k = ind_loop(6)
           do j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
              do i = ind_loop(1), ind_loop(2)
                l    = inddm(  i,      j, k              ) 
                ldp  = indven( i,      j, ind_loop(6)+1  )
                lmtr = indmtr( i,      j, ind_loop(6)+1  )
#include        "FastS/BC/BCWallModel_k.for"
              enddo 
           enddo 
          do k = ind_loop(5), ind_loop(6)-1
            do j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
              do i = ind_loop(1), ind_loop(2)
                l    = inddm(  i,      j, k              ) 
                ldp  = indven( i,      j, ind_loop(6)+1  )  !next rank
                lmtr = indmtr( i,      j, ind_loop(6)+1  )  !next rank
#include        "FastS/BC/BCWallModel_k.for"                  
              enddo 
            enddo !j
          enddo !k

       else
          k = ind_loop(6)
           do j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
              do i = ind_loop(1), ind_loop(2)
                l    = inddm(  i,      j, k              ) 
                ldp  = indven( i,      j, ind_loop(6)+1  )
                lmtr = indmtr( i,      j, ind_loop(6)+1  )
#include        "FastS/BC/BCWallModelSA_k.for"
              enddo 
           enddo 
          do k = ind_loop(5), ind_loop(6)-1
            do j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
              do i = ind_loop(1), ind_loop(2)
                l    = inddm(  i,      j, k              ) 
                ldp  = indven( i,      j, ind_loop(6)+1  )  !next rank
                lmtr = indmtr( i,      j, ind_loop(6)+1  )  !next rank
#include        "FastS/BC/BCWallModelSA_k.for"                
              enddo 
            enddo !j
          enddo !k
       endif !param_int(NEQ)


      ELSE 

       kref = 2*ind_loop(5) - 1
       inck =-param_int(NIJK)*param_int(NIJK+1)
       incp = param_int(IJKV+2)*inck
       kr   =   ind_loop(5) - 1

       if(param_int(NEQ).eq.5) then

          k = ind_loop(5)
           do j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
              do i = ind_loop(1), ind_loop(2)
                l    = inddm(  i,      j, k            ) 
                ldp  = indven( i,      j, ind_loop(5)  )
                lmtr = indmtr( i,      j, ind_loop(5)  )
#include        "FastS/BC/BCWallModel_k.for"
              enddo 
           enddo 
          do k = ind_loop(5)+1, ind_loop(6)
            do j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
              do i = ind_loop(1), ind_loop(2)
                l    = inddm(  i,      j, k            ) 
                ldp  = indven( i,      j, ind_loop(5)  )  !next rank
                lmtr = indmtr( i,      j, ind_loop(5)  )  !next rank
#include        "FastS/BC/BCWallModel_k.for"                
              enddo 
            enddo !j
          enddo !k

       else

          k = ind_loop(5)
           do j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
              do i = ind_loop(1), ind_loop(2)
                l    = inddm(  i,      j, k            ) 
                ldp  = indven( i,      j, ind_loop(5)  )
                lmtr = indmtr( i,      j, ind_loop(5)  )
#include        "FastS/BC/BCWallModelSA_k.for"
              enddo 
           enddo 
          do k = ind_loop(5)+1, ind_loop(6)
            do j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
              do i = ind_loop(1), ind_loop(2)
                l    = inddm(  i,      j, k            ) 
                ldp  = indven( i,      j, ind_loop(5)  )  !next rank
                lmtr = indmtr( i,      j, ind_loop(5)  )  !next rank
#include        "FastS/BC/BCWallModelSA_k.for"              
              enddo 
            enddo !j
          enddo !k

       endif !param_int(NEQ)

      ENDIF !idir

      END
