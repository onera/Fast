c***********************************************************************
c     $Date: 2010-11-04 13:25:50 +0100 (Thu, 04 Nov 2010) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine navier_stokes_struct( ndo, Nbre_thread_actif,
     &        ithread, ithread_io, 
     &        omp_mode, layer_mode, Nbre_socket, socket, mx_synchro, 
     &        lssiter_verif, lexit_lu,
     &        nptpsi, nitcfg, nssiter, nitrun, first_it,
     &        nb_pulse, flagCellN,
     &        param_int, param_real,
     &        temps, tot,
     &        ijkv_sdm,
     &        ind_dm_zone, ind_dm_socket, ind_dm_omp,
     &        socket_topology, lok , topo_s, timer_omp,
     &        krylov, norm_kry,
     &        cfl,
     &        x , y , z, cellN, cellN_IBC,
     &        rop , rop_m1     , rop_tmp , rop_ssiter,
     &        xmut  ,
     &        ti, tj, tk, vol,  ti_df, tj_df, tk_df, vol_df,
     &        venti , ventj , ventk ,
     &        wig , stat_wig, rot,
     &        drodm , coe, delta, ro_res, ro_src)

c***********************************************************************
c_U   USER : TERRACOL
c
c     ACT
c_A    Appel du calcul des flux explicites
c
c     VAL
c_V    gaz parfait monoespece
c_V    processeur domaine
c_V    steady/unsteady
c
c     INP
c_I    tijk     : vecteur param_int( IO_THREAD)rmale aux facettes des mailles
c_I    vent     : vitesses d'entrainement aux facettes preced.
c
c     LOC
c_L    flu      : flux convectifs dans une direction de maillage
c
c     I/O
c_/    drodm    : increment de la solution
c
c***********************************************************************
      implicit none
#include "FastS/param_solver.h"

      INTEGER_E ndo, Nbre_thread_actif , mx_synchro, first_it,
     & ithread, ithread_io, Nbre_socket, socket, nitrun, nptpsi, nitcfg,
     & nb_pulse,lexit_lu,
     & lssiter_verif,flagCellN,omp_mode,layer_mode,nssiter
c
      INTEGER_E  ijkv_sdm(3),ind_dm_zone(6), topo_s(3),
     & ind_dm_omp(6), ind_dm_socket(6), socket_topology(3),
     & param_int(0:*), lok(*), cycl

      REAL_E rop(*),rop_m1(*),rop_tmp(*),rop_ssiter(*),xmut(*),drodm(*),
     & coe(*), ti(*),tj(*),tk(*),vol(*),x(*),y(*),z(*),
     & venti(*),ventj(*),ventk(*), wig(*),stat_wig(*), rot(*), celln(*),
     & ti_df(*),tj_df(*),tk_df(*),vol_df(*), krylov(*), cellN_IBC(*), 
     & ro_src(*),timer_omp(*)

      REAL_E delta(*),ro_res(*)

      REAL_E psi(nptpsi)

      REAL_E temps, norm_kry, cfl(3), param_real(0:*)
c      REAL_E drodmstk(20000,param_int(NEQ)), rostk(20000,param_int(NEQ))

C Var loc
#include "FastS/HPC_LAYER/LOC_VAR_DECLARATION.for"



      INTEGER_E tot(6,Nbre_thread_actif), totf(6), glob(4), 
     & ind_loop(6),neq_rot,depth,nb_bc,thmax,th,shift1,shift2,
     & flag_wait,cells, it_dtloc


      REAL_E rhs_begin,rhs_end


#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      include 'omp_lib.h'

       rhs_begin = OMP_GET_WTIME()

#include "FastS/HPC_LAYER/SIZE_MIN.for"
#include "FastS/HPC_LAYER/WORK_DISTRIBUTION_BEGIN.for"
#include "FastS/HPC_LAYER/LOOP_CACHE_BEGIN.for"
#include "FastS/HPC_LAYER/INDICE_RANGE.for"

          flag_wait = 0

c         if(ithread_io.eq.24.and.nitcfg.le.1.and.icache.eq.1
c     &    .and.jcache.eq.1) then
c         if(ithread.eq.param_int( IO_THREAD).and.nitcfg.le.1) then
c           write(*,'(a,7i6)')'ijkv_sdm  =',ijkv_sdm,icache,jcache,
c     &                                     kcache,ithread
c           write(*,'(a,9i6)')'ind_dm_zone=',ind_dm_zone,
c     &    ithread,jbloc,ndo
c           write(*,'(a,3i6)')'topo=',topo_omp
c           write(*,'(a,3i6)')'loop_patern=',shift
c           write(*,'(a,6i6)')'ind_dm_soc=',ind_dm_socket
c           write(*,'(a,6i6)')'inddm_omp =',inddm_omp
c           write(*,'(a,6i6)')'ind_dm_thr=',ind_dm_omp
c           write(*,'(a,6i6)')'ind_sdm   =',ind_sdm
c           write(*,'(a,6i6)')'ind_ssa   =',ind_ssa
c           write(*,'(a,6i6)')'ind_coe   =',ind_coe
c           write(*,'(a,6i6)')'ind_grad  =',ind_grad
c           write(*,'(a,6i6)')'ind_rhs   =',ind_rhs
c           write(*,'(a,6i6)')'ind_mjr   =',ind_mjr
c         endif


c      if(ndo.eq.-9) then
c         thmax = OMP_get_num_threads()
c         th = OMP_get_thread_num()+1

c         do i = 1, thmax
c            if (i == th) then
c        write(*,'(a,7i5)')'ind_dm  =',ind_dm_omp(1:4),th,icache,jcache
c            write(*,'(a,4i5)')'ind_sdm =',ind_sdm(1:4)
c            write(*,'(a,4i5)')'ind_coe =',ind_coe(1:4)
c            write(*,'(a,4i5)')'ind_ssa =',ind_ssa(1:4)
c            write(*,'(a,4i5)')'ind_gra =',ind_grad(1:4)
c            endif
c!$OMP BARRIER
c         enddo 
c       endif

          
         !!! variable pour l'explicit local !!!!
         !!! utilisee pour savoir si on est a la premiere ss-ite des zones !!!
         cycl=param_int(NSSITER)/param_int(LEVEL) 

 
         IF (param_int(EXPLOC) .eq. 0) then !!! explicit global ou implicit


             IF(nitcfg.eq.1) then



              if(param_int(LALE).eq.1) then ! mise a jour Vent et tijk si mvt corps solide
                call mjr_ale(ndo,nitcfg, ithread,
     &                        param_int, param_real,
     &                        ind_dm_zone, ind_sdm,ijkv_thread,ijkv_sdm,
     &                        synchro_send_sock, synchro_send_th,
     &                        synchro_receive_sock, synchro_receive_th,
     &                        ibloc , jbloc , kbloc ,
     &                        icache, jcache, kcache,
     &                        x,y,z,ti,ti_df,tj,tj_df,tk,tk_df,vol,
     &                        venti, ventj, ventk)

                endif

                !Calcul de la viscosite laminaire si nslaminar ou (nsles + dom 2D)
                if(param_int(IFLOW).eq.2) then
                 


                 if(param_int(ILES).eq.0.or.param_int(NIJK+4).eq.0) then
                
                   if(param_int(ITYPCP).le.1) then 
                     call invist(ndo, param_int, param_real, ind_coe,
     &                           rop_ssiter, xmut )
                   else 
                     call invist(ndo, param_int, param_real, ind_grad,
     &                           rop_ssiter, xmut )
                   endif


                 !LES selective mixed scale model
                 else
                   neq_rot = 3
                   depth   = 1 !pour extrapolation mut, on travaille sur 1 seule rangee
                   call lesvist(ndo, param_int,param_real,neq_rot,depth,
     &                          ithread, nitrun,
     &                          ind_grad, ind_coe, ind_dm_zone,
     &                          xmut, rop_ssiter, ti,tj,tk, vol, rot)
                 endif

                elseif(param_int(IFLOW).eq.3) then

                !! remplissage tableau xmut si SA uniquememnt. 
                !! Pour ZDES, remplissage dans terme source 
                !call vispalart(ndo, param_int, param_real, ind_grad,
                call vispalart(ndo, param_int, param_real, ind_coe,
     &                         xmut,rop_ssiter)

              endif

             !!sinon cfl foireux
             IF(param_int(IFLOW).eq.3.and.param_int(ITYPCP).le.1) then
#include      "FastS/HPC_LAYER/SYNCHRO_WAIT.for"
               flag_wait = 1
             ENDIF

             ! Calcul du pas de temps
             call cptst3(ndo, nitcfg, nitrun, first_it, lssiter_verif,
     &                   flagCellN, param_int, param_real,
     &                   ind_sdm, ind_grad, ind_coe,
     &                   cfl, xmut,rop_ssiter, cellN, coe,
     &                   ti,tj,tk, vol,venti)

           ENDIF!!1ere sous-iteration



         ELSE !!!  explicite local instationnaire
         

             IF(mod(nitcfg,cycl).eq.1) then !!! On teste si on est a la premiere ss-ite de la zone



                if(param_int(LALE).ge.1) then ! mise a jour Vent et tijk
                   call mjr_ale(ndo,nitcfg, ithread,
     &                        param_int, param_real,
     &                        ind_dm_zone, ind_sdm,ijkv_thread,ijkv_sdm,
     &                        synchro_send_sock, synchro_send_th,
     &                        synchro_receive_sock, synchro_receive_th,
     &                        ibloc , jbloc , kbloc ,
     &                        icache, jcache, kcache,
     &                        x,y,z,ti,ti_df,tj,tj_df,tk,tk_df,vol,
     &                        venti, ventj, ventk)

                endif

                !Calcul de la viscosite laminaire si nslaminar ou (nsles + dom 2D)
                if(param_int(IFLOW).eq.2) then
                 


                if(param_int(ILES).eq.0.or.param_int(NIJK+4).eq.0) then
                
                   if(param_int(ITYPCP).le.1) then 
                     call invist(ndo, param_int, param_real, ind_coe,
     &                           rop_ssiter, xmut )
                   else 
                     call invist(ndo, param_int, param_real, ind_grad,
     &                           rop_ssiter, xmut )
                   endif


                !LES selective mixed scale model
                else
#include          "FastS/HPC_LAYER/SYNCHRO_WAIT.for"
                   flag_wait = 1

                   neq_rot = 3
                   depth   = 1 !pour extrapolation mut, on travaille sur 1 seule rangee
                   call lesvist(ndo, param_int,param_real,neq_rot,depth,
     &                          ithread,
     &                          ind_grad, ind_coe, ind_dm_zone,
     &                          xmut, rop_ssiter, ti,tj,tk, vol, rot)
                endif

                elseif(param_int(IFLOW).eq.3) then

                !! remplissage tableau xmut si SA uniquememnt. 
                !! Pour ZDES, remplissage dans terme source 
                !call vispalart(ndo, param_int, param_real, ind_grad,
                call vispalart(ndo, param_int, param_real, ind_coe,
     &                         xmut,rop_ssiter)

              endif

             !!sinon cfl foireux
             IF(param_int(IFLOW).eq.3.and.param_int(ITYPCP).le.1) then
#include      "FastS/HPC_LAYER/SYNCHRO_WAIT.for"
               flag_wait = 1
             ENDIF

             ! Calcul du pas de temps
             call cptst3(ndo, nitcfg, nitrun, first_it, lssiter_verif,
     &                   flagCellN, param_int, param_real,
     &                   ind_sdm, ind_grad, ind_coe,
     &                   cfl, xmut,rop_ssiter, cellN, coe,
     &                   ti,tj,tk, vol,venti)

           ENDIF !!1ere sous-ite de la zone


         ENDIF !! fin test explicit local

           
           !SI SA implicit, verrou ici car dependence entre coe(5)
           !calculee sur ind_coe dans cptst3 et coe(6) calculee sur
           !ind_ssa dans src_term
           ! pareil pour ibc0
           IF( (    (param_int(IFLOW).eq.3.and.param_int(ITYPCP).le.1)
     &          .or.param_int(IBC).eq.1
     &         )
     &                             .and.flag_wait.eq.0) then
#include "FastS/HPC_LAYER/SYNCHRO_WAIT.for"
               flag_wait = 1
           ENDIF





           ! - Ajout d'un eventuel terme source au second membre
           ! - initialisation drodm
           call src_term(ndo, nitcfg, nb_pulse, param_int, param_real,
     &                   ind_sdm, ind_rhs, ind_ssa,
     &                   temps,
     &                   rop_ssiter, xmut, drodm, coe, x,y,z,cellN_IBC,
     &                   ti,tj,tk,vol, delta, ro_src, wig)


           IF(flag_wait.eq.0) then
#include "FastS/HPC_LAYER/SYNCHRO_WAIT.for"
           ENDIF


           !print*, param_int(KFLUDOM)

           ! -----assemblage drodm euler+visqueux
              if(param_int(KFLUDOM).eq.1) then

                call fluausm_select(ndo,nitcfg, ithread,nptpsi,
     &                        param_int, param_real,
     &                        ind_dm_zone, ind_sdm,ijkv_thread,ijkv_sdm,
     &                        synchro_send_sock, synchro_send_th,
     &                        synchro_receive_sock, synchro_receive_th,
     &                        ibloc , jbloc , kbloc ,
     &                        icache, jcache, kcache,
     &                        psi,wig,stat_wig, rop_ssiter, drodm,
     &                        ti,ti_df,tj,tj_df,tk,tk_df, vol,vol_df,
     &                        venti, ventj, ventk, xmut)

              elseif(param_int(KFLUDOM).eq.2.and.nitcfg.eq.1 .and.
     & param_int(EXPLOC) .eq. 0 ) then
   
             !elseif(param_int(KFLUDOM).eq.2.and.nitcfg.eq.1 .and.
     & !param_int(EXPLOC) .eq. 0 .or. param_int(KFLUDOM).eq.2
     & !.and.nitcfg.eq.5 .and.param_int(EXPLOC) .eq. 0 ) then


              
                call flusenseur_init_select(ndo,nitcfg,ithread, nptpsi,
     &                        param_int, param_real,
     &                        ind_dm_zone, ind_sdm,ijkv_thread,ijkv_sdm,
     &                        synchro_send_sock, synchro_send_th,
     &                        synchro_receive_sock, synchro_receive_th,
     &                        ibloc , jbloc , kbloc ,
     &                        icache, jcache, kcache,
     &                        psi,wig,stat_wig, rop_ssiter, drodm,
     &                        ti,ti_df,tj,tj_df,tk,tk_df, vol,vol_df,
     &                        venti, ventj, ventk, xmut)


     

              elseif(param_int(KFLUDOM).eq.2.and.nitcfg.gt.1 .and.
     & param_int(EXPLOC) .eq. 0 ) then

              !elseif(param_int(KFLUDOM).eq.2.and.nitcfg.ne.1 .and.
     & !param_int(EXPLOC) .eq. 0 .and.nitcfg.ne.5 ) then
 


                call flusenseur_select(ndo,nitcfg,ithread, nptpsi,
     &                        param_int, param_real,
     &                        ind_dm_zone, ind_sdm,ijkv_thread,ijkv_sdm,
     &                        synchro_send_sock, synchro_send_th,
     &                        synchro_receive_sock, synchro_receive_th,
     &                        ibloc , jbloc , kbloc ,
     &                        icache, jcache, kcache,
     &                        psi,wig,stat_wig, rop_ssiter, drodm,
     &                        ti,ti_df,tj,tj_df,tk,tk_df, vol,vol_df,
     &                        venti, ventj, ventk, xmut)


              elseif(param_int(KFLUDOM).eq.2.and.mod(nitcfg,cycl).eq.1
     & .and. param_int(EXPLOC) .ne. 0 ) then !!! senseur pour l'explicit local
                 
                 
                call flusenseur_init_select(ndo,nitcfg,ithread, nptpsi,
     &                        param_int, param_real,
     &                        ind_dm_zone, ind_sdm,ijkv_thread,ijkv_sdm,
     &                        synchro_send_sock, synchro_send_th,
     &                        synchro_receive_sock, synchro_receive_th,
     &                        ibloc , jbloc , kbloc ,
     &                        icache, jcache, kcache,
     &                        psi,wig,stat_wig, rop_ssiter, drodm,
     &                        ti,ti_df,tj,tj_df,tk,tk_df, vol,vol_df,
     &                        venti, ventj, ventk, xmut)


     

              elseif(param_int(KFLUDOM).eq.2.and.mod(nitcfg,cycl).ne.1
     & .and. param_int(EXPLOC) .ne. 0 ) then !!! senseur pour l'explicit local

 
                call flusenseur_select(ndo,nitcfg,ithread, nptpsi,
     &                        param_int, param_real,
     &                        ind_dm_zone, ind_sdm,ijkv_thread,ijkv_sdm,
     &                        synchro_send_sock, synchro_send_th,
     &                        synchro_receive_sock, synchro_receive_th,
     &                        ibloc , jbloc , kbloc ,
     &                        icache, jcache, kcache,
     &                        psi,wig,stat_wig, rop_ssiter, drodm,
     &                        ti,ti_df,tj,tj_df,tk,tk_df, vol,vol_df,
     &                        venti, ventj, ventk, xmut)





              elseif(param_int(KFLUDOM).eq.5) then

                call fluroe_select(ndo,nitcfg,ithread, nptpsi,
     &                        param_int, param_real,
     &                        ind_dm_zone, ind_sdm,ijkv_thread,ijkv_sdm,
     &                        synchro_send_sock, synchro_send_th,
     &                        synchro_receive_sock, synchro_receive_th,
     &                        ibloc , jbloc , kbloc ,
     &                        icache, jcache, kcache,
     &                        psi,wig,stat_wig, rop_ssiter, drodm,
     &                        ti,ti_df,tj,tj_df,tk,tk_df, vol,vol_df,
     &                        venti, ventj, ventk, xmut)


             elseif(param_int(KFLUDOM).eq.8) then

                call flushyper_select(ndo,nitcfg,ithread, nptpsi,
     &                        param_int, param_real,
     &                        ind_dm_zone, ind_sdm,ijkv_thread,ijkv_sdm,
     &                        synchro_send_sock, synchro_send_th,
     &                        synchro_receive_sock, synchro_receive_th,
     &                        ibloc , jbloc , kbloc ,
     &                        icache, jcache, kcache,
     &                        psi,wig,stat_wig, rop_ssiter, drodm,
     &                        ti,ti_df,tj,tj_df,tk,tk_df, vol,vol_df,
     &                        venti, ventj, ventk, xmut)



              else

                 if(ithread.eq.1) then
                   write(*,*)'Unknown flux',param_int(KFLUDOM)
                 endif
             endif

#include "FastS/HPC_LAYER/SYNCHRO_GO.for"

           !correction flux roe au CL si pas de mvt ALE,....
           nb_bc = param_int( param_int(PT_BC) )
           if(param_int(KFLUDOM).eq.5.and.nb_bc.ne.0) then

              call bfl3(ndo, ithread, param_int, param_real, 
     &                  ind_dm_zone, ind_sdm,
     &                  psi,wig,stat_wig, rop_ssiter, drodm,
     &                  ti,ti_df,tj,tj_df,tk,tk_df, vol,vol_df,
     &                  venti, ventj, ventk, xmut)
           endif

           !boundary flux correction for LES wall model
           if(param_int(NEQ).eq.5.and.nb_bc.ne.0) then
              call wall_model_flux(ndo,ithread, param_int, param_real,
     &                  ind_dm_zone, ind_sdm, nitcfg, nitrun, cycl,
     &                  psi,wig,stat_wig, rop_ssiter, drodm,x,y,z,
     &                  ti,ti_df,tj,tj_df,tk,tk_df, vol,vol_df,
     &                  venti, ventj, ventk, xmut)
           endif

           !Extraction tableau residu
           if(param_int(EXTRACT_RES).eq.1) then
              call extract_res(ndo, param_int, param_real,
     &                         ind_mjr,
     &                         drodm, vol, ro_res)
           endif
 
           !! impicit krylov             
           if(param_int(ITYPCP).le.1.and.
     &        (param_int(IMPLICITSOLVER).eq.1.and.layer_mode.eq.1)) then
              !Assemble Residu Newton; 3q(n+1)-4Q(n)+q(n-1)) + dt (flu(i+1)-(flu(i)) =0
              if(flagCellN.eq.0) then
               call core3as2_kry(ndo,nitcfg, first_it,
     &                           param_int,param_real,
     &                           ind_mjr,
     &                           krylov, norm_kry,
     &                           rop_ssiter, rop, rop_m1, drodm, coe)
               else
               call core3as2_chim_kry(ndo,nitcfg, first_it, 
     &                                param_int,param_real,
     &                                ind_mjr, cellN,
     &                                krylov, norm_kry,
     &                                rop_ssiter, rop, rop_m1,drodm,coe)
               endif


           !! implicit Lu                 
           elseif(param_int(ITYPCP).le.1) then

              !Assemble Residu Newton; 3q(n+1)-4Q(n)+q(n-1)) + dt (flu(i+1)-(flu(i)) =0
              if(flagCellN.eq.0) then
               call core3as2(ndo,nitcfg, first_it, param_int,param_real,
     &                        ind_mjr,
     &                       rop_ssiter, rop, rop_m1, drodm, coe)
              else
               call core3as2_chim(ndo,nitcfg, first_it, 
     &                            param_int,param_real,
     &                            ind_mjr, cellN,
     &                            rop_ssiter, rop, rop_m1, drodm, coe)
              endif

             !Extraction tableau residu
             if( param_int(EXTRACT_RES).eq.2.and.
     &           (nitcfg.eq.1.or.nitcfg.eq.nssiter) ) then

                shift2=1
                if(param_int(IFLOW).eq.3) shift2=2

                if(nitcfg.eq.1)  then
                   shift1 = param_int(NDIMDX)*shift2
                else
                   shift1 = param_int(NDIMDX)*(param_int(NEQ)+shift2)
                endif

                call extract_res(ndo, param_int, param_real,
     &                           ind_mjr,
     &                         drodm, vol, xmut(1+shift1) )
             endif
           !! explicit Lu                 
           else
             !c--------------------------------------------------
             !calcul param_int( IO_THREAD)uveau champ en explicite rk3
             !c-----Mise a jour de iptdrodm0 et de la solution
             !c     iptrotmp(nitcfg+1)=iptrotmp(nitcfg)+CoefW *drodm (rk3)
             !c--------------------------------------------------



              if (param_int(EXPLOC)==0) then ! Explicit global


                 if(mod(nitcfg,2)==1) then
                       
                     call core3ark3(ndo,nitcfg, param_int, param_real,
     &                      ind_mjr,
     &                   rop_tmp, rop, drodm, coe)
                 else

                     call core3ark3(ndo,nitcfg, param_int, param_real,
     &                   ind_mjr,
     &                   rop, rop_tmp, drodm, coe)

                 end if
           
    
             else if(param_int(EXPLOC)==1) then ! Explicit local instationnaire
             !else if(param_int(EXPLOC)==2) then ! Explicit local instationnaire

                 if (MOD(nitcfg,2)==0) then

                     call core3_dtloc(ndo,nitcfg, param_int, param_real,
     &                    ind_mjr,
     &                    rop, rop_tmp, drodm, coe)
                 else

                     call core3_dtloc(ndo,nitcfg, param_int, param_real,
     &                    ind_mjr,
     &                   rop_tmp, rop, drodm, coe)
                     
                 endif

                  
              endif
          end if

#include "FastS/HPC_LAYER/LOOP_CACHE_END.for"
#include "FastS/HPC_LAYER/WORK_DISTRIBUTION_END.for"



           !!! omp barriere
      if(lexit_lu.eq.0.and.nitrun*nitcfg.gt.15
     &                .and.(nitcfg.lt.3.or.nitcfg.eq.nssiter-1)) then 
        rhs_end = OMP_GET_WTIME()
        cells = (ind_dm_omp(2)-ind_dm_omp(1)+1)
     &         *(ind_dm_omp(4)-ind_dm_omp(3)+1)
     &         *(ind_dm_omp(6)-ind_dm_omp(5)+1)

        timer_omp(1)=timer_omp(1)+(rhs_end-rhs_begin)/float(cells)
        timer_omp(2)=float(cells)
      endif



      end
