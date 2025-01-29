c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine src_term(ndom, ithread, nitcfg, nb_pulse,
     &                    param_int, param_real,
     &                    ind_sdm, ind_rhs, ind_ssa, ind_grad,
     &                    temps, nitrun,cycl,
     &                    rop, xmut, drodm, coe, x,y,z,cellN_IBC,
     &            ti, tj, tk, vol, delta, ro_src, wig, rop_p1, rop_m1)
c***********************************************************************
c_P                          O N E R A
c
c_DC  DATE_C : Novembre 2005 - M. Terracol 
c
c     HISTORIQUE
c
c     ACT
c_A    Ajout d un terme source aux membre de droite des equations
c     VAL
c     INP
c_I    ndom   : numero du domaine
c     OUT
c     I/O
c      drodm  : bilan de flux
C-----------------------------------------------------------------------
      implicit  none

      INTEGER_E ndom,nitcfg,nb_pulse,ind_sdm(6),ind_rhs(6),ind_ssa(6),
     & ind_grad(6),param_int(0:*), nitrun,cycl, ithread
c

      REAL_E rop(*),xmut(*),drodm(*),coe(*)
      REAL_E ti(*),tj(*),tk(*), vol(*)
      REAL_E x(*),y(*),z(*), cellN_IBC(*)
      REAL_E param_real(0:*), temps
      REAL_E delta(*), ro_src(*), wig(*),rop_m1(*),rop_p1(*)

C Var loc
      INTEGER_E nacp,l,idirx,idirz,i,j,k,n,n2,i1,j1,k1,i2m1,j2m1,k2m1,
     & i2,j2,k2,ind1,ind2,lt,n_tot,ndt,nsp,ix,iy,iz,incr,indf,ci,cj,ck,
     & im,jm,km,indmy,imin,jmin,kmin,imax,jmax,kmax,l1,l2,l3,lm,lp,
     & iptribc,niloc,njloc,m,igetadr,ind2d,ierr,neq_wig,iwig

      REAL_E coefa,coefb,x0,y0,z0,amp,per,phi,qnp1,dvzparoi,
     & q_tot,qf_tot,
     & Lx,Ly,Lzz,alpha,beta,Syz,Re,voltot,vmoy,dzparoi,ectot,ran1,f,
     & xiadm2,xpente,dxinew,dxiprev,q_inv,q_res,q_res2,fn,yp,yplus,
     & rho,vu,vv,vw,q,eps,ampli,xmod,pi,xrelax,sro,su,sv,sw,sroe,seuil,
     & xc,yc,rayon,rayon1,rayon2,rayon3,rayon4,nuwall,p,roe,dphi,aphi,
     & gam1,constA,uext,dpdx,dy,dym,up,utp2,nu,nut,vara,varb,rho0,p0,
     & rhou0,rhov0,rhow0,roe0

#include "FastS/param_solver.h"


      !Initilalisation systematique de drodm
      if (param_int(EXPLOC).ne.0) then !! explicite local instationnaire
 
            call init_rhs_dtloc(ndom, nitcfg, param_int, 
     &       param_int(NDIMDX),param_int( NEQ ), ind_ssa,  drodm )
               

      else  !! pas explicite local
 
         call init_rhs(ndom, nitcfg, param_int, param_int( NDIMDX ),
     &        param_int( NEQ ), ind_rhs,  drodm )

      end if


      IF (param_int(SRC).ge.1) THEN
        call update_src(ndom, nitcfg, param_int, param_real, ind_rhs,
     &                  drodm , rop, coe, vol, ro_src)
      ENDIF

      ! Si dtloc + ALE 
      ! => on resoud la conservation du moment dans le repere relatif
      ! => terme source volumique ajoute au bilan de flux (force centrifuge)
      IF (param_int(DTLOC).eq.1 .and. param_real(ROT_TETAP).ne.0.) THEN
        call spsource_MEANFLOW(ndom, nitcfg, param_int, param_real,
     &                         ind_ssa, rop ,coe, vol, drodm)
      ENDIF
c---------------------------------------------------------------------
c----- Calcul des termes sources de l''equation de Spalart Allmaras et de mut pour ZDES
c---------------------------------------------------------------------
c
      IF (param_int(IFLOW).eq.3.and.param_int(ILES).eq.0) THEN

            if(param_int(SA_INT+ SA_IDES-1).eq.0) then !SA

             call spsource_SA(ndom, nitcfg, param_int, param_real,
     &                        ind_ssa,
     &                        xmut, rop, coe, ti, tj, tk, vol,
     &                        xmut(1+param_int(NDIMDX)),       !dist
     &                        drodm)

            elseif(param_int(SA_INT+ SA_IDES-1).eq.1) then !SA_comp

             call spsource_SA_comp(ndom, nitcfg, param_int, param_real,
     &                             ind_ssa,
     &                             xmut, rop, coe, ti, tj, tk, vol,
     &                             xmut(1+param_int(NDIMDX)),       !dist
     &                             drodm)

            elseif(param_int(SA_INT+ SA_IDES-1).eq.8) then !SA_comp

             call spsource_SA_diff(ndom, nitcfg, param_int, param_real,
     &                             ind_ssa,
     &                             xmut, rop, coe, ti, tj, tk, vol,
     &                             xmut(1+param_int(NDIMDX)),       !dist
     &                             drodm)

            elseif(param_int(SA_INT+ SA_IDES-1).eq.2) then !ZDES mode 1, delta_vol

             if(param_int(SA_DEBUG).eq.1) then  !sauvegarde de delta

               call spsource_ZDES1_vol_debug(ndom,param_int, param_real,
     &                                      ind_ssa,
     &                                      xmut,rop,coe, ti,tj,tk,vol,
     &                                      xmut(1+param_int(NDIMDX)),!dist
     &                                      drodm, delta)
             else
               call spsource_ZDES1_vol(ndom,param_int, param_real,
     &                                 ind_ssa,
     &                                 xmut,rop,coe, ti,tj,tk,vol,
     &                                 xmut(1+param_int(NDIMDX)),!dist
     &                                 drodm)
             endif
 
            elseif(param_int(SA_INT+ SA_IDES-1).eq.3) then !ZDES mode 1, delta_rot

             if(param_int(SA_DEBUG).eq.1) then  !sauvegarde de delta

               call spsource_ZDES1_rot_debug(ndom,param_int, param_real,
     &                                      ind_ssa,
     &                                      xmut,rop,coe, ti,tj,tk,vol,
     &                                      xmut(1+param_int(NDIMDX)),!dist
     &                                      drodm, delta)
             else
               call spsource_ZDES1_rot(ndom,param_int, param_real,
     &                                 ind_ssa,
     &                                 xmut,rop,coe, ti,tj,tk,vol,
     &                                 xmut(1+param_int(NDIMDX)),!dist
     &                                 drodm)
             endif

            elseif(param_int(SA_INT+ SA_IDES-1).eq.4) then !ZDES mode 2, delta_vol

             if(param_int(SA_DEBUG).eq.1) then  !sauvegarde de delta

               call spsource_ZDES2_vol_debug(ndom,param_int, param_real,
     &                                       ind_ssa,
     &                                       xmut,rop,coe, ti,tj,tk,vol,
     &                                       xmut(1+param_int(NDIMDX)), !dist
     &                                       delta, 
     &                                       delta(1+param_int(NDIMDX)),!fd  
     &                                       drodm)
             else
             call spsource_ZDES2_vol(ndom, param_int, param_real,
     &                               ind_ssa,
     &                               xmut,rop,coe, ti,tj,tk,vol,
     &                               xmut(1+param_int(NDIMDX)),       !dist
     &                               drodm)
             endif
 
            elseif(param_int(SA_INT+ SA_IDES-1).eq.5) then !ZDES mode 2, delta_rot

             if(param_int(SA_DEBUG).eq.1) then  !sauvegarde de delta

               call spsource_ZDES2_rot_debug(ndom,param_int, param_real,
     &                                       ind_ssa,
     &                                       xmut,rop,coe, ti,tj,tk,vol,
     &                                       xmut(1+param_int(NDIMDX)), !dist
     &                                       delta, 
     &                                       delta(1+param_int(NDIMDX)),!fd  
     &                                       drodm)
             else
             call spsource_ZDES2_rot(ndom, param_int, param_real,
     &                               ind_ssa,
     &                               xmut,rop,coe, ti,tj,tk,vol,
     &                               xmut(1+param_int(NDIMDX)),       !dist
     &                               drodm)
             endif

            elseif(param_int(SA_INT+ SA_IDES-1).eq.6) then !ZDES mode 3, delta_vol

             if(param_int(SA_DEBUG).eq.1) then  !sauvegarde de delta

               call spsource_ZDES3_vol_debug(ndom,param_int, param_real,
     &                                 ind_ssa,
     &                                 xmut,rop,coe, ti,tj,tk,vol,
     &                                 xmut(1+param_int(NDIMDX)),       !dist
     &                                 xmut(1+param_int(NDIMDX)*2),     !zgris
     &                                 delta, drodm)
             else
               call spsource_ZDES3_vol(ndom, param_int, param_real,
     &                                 ind_ssa,
     &                                 xmut,rop,coe, ti,tj,tk,vol,
     &                                 xmut(1+param_int(NDIMDX)),       !dist
     &                                 xmut(1+param_int(NDIMDX)*2),     !zgris
     &                                 drodm)
             endif

            elseif(param_int(SA_INT+ SA_IDES-1).eq.7) then !ZDES mode 3, delta_rot

             if(param_int(SA_DEBUG).eq.1) then  !sauvegarde de delta

               call spsource_ZDES3_rot_debug(ndom,param_int, param_real,
     &                                      ind_ssa,
     &                                      xmut,rop,coe, ti,tj,tk,vol,
     &                                      xmut(1+param_int(NDIMDX)),  !dist
     &                                      xmut(1+param_int(NDIMDX)*2),!zgris
     &                                      delta, drodm)
             else
               call spsource_ZDES3_rot(ndom, param_int, param_real,
     &                                 ind_ssa,
     &                                 xmut,rop,coe, ti,tj,tk,vol,
     &                                 xmut(1+param_int(NDIMDX)),       !dist
     &                                 xmut(1+param_int(NDIMDX)*2),     !zgris
     &                                 drodm)
             endif

            else
!$OMP SINGLE
              write(*,*)'erreur source SA: pas code...'
!$OMP END SINGLE
              stop
            endif

      ENDIF

c***********************************************************************
c**   Source term IBC a la funk ordre 0

      if(param_int(IBC).eq.1) then 

        call ibcsource(ndom,param_int,ind_rhs,
     &                 rop, cellN_IBC, coe, drodm)
      endif


c***********************************************************************
c**    Proto senseur de choc Sciacovelli directionnel

      !init wig(:, 4-6) a chaque ssiter: o3sc ou o5sc
      IF(param_int(SLOPE).eq.6.or.param_int(SLOPE).eq.7) THEN 

          neq_wig=3
          if(param_int(KFLUDOM).eq.2) neq_wig = 6

          call sciacovelli(ndom,neq_wig, param_int, param_real,
     &               nitrun, ind_sdm, rop, ti, tj, tk, vol, wig )
      ENDIF

      !!senseur scheme: init sensor
      iwig=0
      if(   (nitcfg.eq.1.and.param_int(EXPLOC).eq.0)
     &  .or.(mod(nitcfg,cycl).eq.1.and.param_int(EXPLOC).ne.0)) iwig=1

      IF(param_int(KFLUDOM).eq.2.and.iwig.eq.1) THEN

       if(param_int(SLOPE).ne.5.and.param_int(SLOPE).ne.7) then
         call wiggle2001(ndom,param_int, param_real, ind_sdm, rop, wig)

       else !!wiggle pour o5 ou o5sc
         call wiggle2023(ndom, param_int, param_real,
     &              ind_grad, ind_sdm, rop, rop_p1, rop_m1, wig, temps)
       endif
        
c          call corr_wiggle2023(ndom, param_int, param_real,
c     &              ind_rhs, ti,tj,tk, rop, wig, drodm)

      ENDIF
        

c***********************************************************************
c***********************************************************************
c***********************************************************************
c
c     Source harmonique (pulse acoustique)
c
      !IF ((iflagpert(ndom).eq.1).and.(nb_pulse.gt.0)) then
      nb_pulse=0
      !if (temps< 0.00000) nb_pulse=1

      IF (nb_pulse.gt.0) then

         do nacp = 1,nb_pulse

           !coefa = acp_a(nacp)
           !coefb = acp_b(nacp) 
           !x0    = acp_x0(nacp)
           !y0    = acp_y0(nacp)
           !z0    = acp_z0(nacp)
           !amp   = acp_amp(nacp)
           !per   = acp_per(nacp)
           !phi   = acp_phi(nacp)
           !coefa = 2.71128
           !coefb = 2.4
           coefa = 2.71128
           !coefb = 0.1
           !coefb = 0.01
           coefb = 0.03
           x0    = 0.5
           y0    = 2.0
           z0    = 0.
           amp   = 10.01
           per   = 0.002
           phi   = 0.

           !temps = param_real(TEMPS)
           
           !print*, 'temps= ',temps

           call ac_pulse(ndom, param_int, ind_rhs,
     &                   temps, param_real(STAREF),
     &                   drodm,rop,x,y,z,vol, 
     &                   x0,y0,z0,coefa,coefb,amp,per,phi)

         end do

      ENDIF

      end
