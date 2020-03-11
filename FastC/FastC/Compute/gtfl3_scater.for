c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 aout 2013) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine gtfl3_scater(ndom, nitcfg, ithread, nptpsi,
     &                        param_int, param_real,
     &                        ind_dm, ind_loop,ijkv_thread, ijkv_sdm,
     &                        synchro_send_sock, synchro_send_th,
     &                        synchro_receive_sock, synchro_receive_th,
     &                        ibloc , jbloc , kbloc ,
     &                        icache, jcache, kcache,
     &                        psi, wig, stat_wig, rop, drodm,
     &                        ti,ti_df,tj,tj_df,tk,tk_df, vol,vol_df,
     &                        venti, ventj, ventk, xmut)
c***********************************************************************
c_U   USER : PECHIER
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
c_I    tijk     : vecteur normale aux facettes des mailles
c_I    ventijk     : vitesses d entrainement aux facettes preced.
c_I    qm,qp    : etats droit et gauche aux interfaces d une maille
c
c     LOC
c_L    flu      : flux convectifs dans une direction de maillage
c
c     I/O
c_/    drodm    : increment de la solution
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ijkv_thread(3), ijkv_sdm(3), ind_loop(6),ind_dm(6),
     & ndom, nitcfg,ithread, nptpsi,
     & icache, jcache, kcache,ibloc, jbloc, kbloc,
     & synchro_send_sock(3),synchro_send_th(3),
     & synchro_receive_sock(3), synchro_receive_th(3), param_int(0:*)


      REAL_E rop(*),xmut(*),drodm(*), ti(*),tj(*),tk(*),vol(*),
     & venti(*),ventj(*),ventk(*), wig(*),stat_wig(*),
     & ti_df(*),tj_df(*),tk_df(*),vol_df(*)

      REAL_E param_real(0:*)
      REAL_E psi(nptpsi)

C Var loc
      INTEGER_E option

      if(ind_loop(1).gt.ind_loop(2)) return 
      if(ind_loop(3).gt.ind_loop(4)) return 
      if(ind_loop(5).gt.ind_loop(6)) return

      option =  1000*param_int(LALE)
     &        +  100*param_int(KFLUDOM)
     &        +   10*param_int(IFLOW)
     &        +      param_int(ITYPZONE)


      IF (option.eq.110) THEN

          call fluthese_euler_3dfull(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSEIF (option.eq.111) THEN

          call fluthese_euler_3dhomo(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSEIF (option.eq.112) THEN

          call fluthese_euler_3dcart(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSEIF(option.eq.113) THEN

          call fluthese_euler_2d(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSEIF (option.eq.120) THEN

          call fluthese_lamin_3dfull(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSEIF (option.eq.121) THEN

          call fluthese_lamin_3dhomo(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSEIF (option.eq.122) THEN

          call fluthese_lamin_3dcart(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSEIF (option.eq.123) THEN

          call fluthese_lamin_2d(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSEIF (option.eq.130) THEN

          call fluthese_SA_3dfull(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSEIF (option.eq.131) THEN

          call fluthese_SA_3dhomo(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSEIF (option.eq.132) THEN

          call fluthese_SA_3dcart(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSEIF (option.eq.133) THEN

          call fluthese_SA_2d(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
!!
!!!!
!!!!  schema senseur
!!!!
!!
      ELSEIF (option.eq.210) THEN

          if(nitcfg.eq.1) then
          call flusenseur_init_euler_3dfull(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
          else
          call flusenseur_euler_3dfull(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
          endif

      ELSEIF (option.eq.211) THEN

          if(nitcfg.eq.1) then
          call flusenseur_init_euler_3dhomo(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
          else
          call flusenseur_euler_3dhomo(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
          endif

      ELSEIF (option.eq.212) THEN

     
          if(nitcfg.eq.1) then
          call flusenseur_init_euler_3dcart(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
          else
          call flusenseur_euler_3dcart(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
          endif


      ELSEIF (option.eq.213) THEN

     
           if(nitcfg.eq.1) then

          call flusenseur_init_euler_2d(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
          else
          call flusenseur_euler_2d(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
          endif

      ELSEIF (option.eq.220) THEN

     
           if(nitcfg.eq.1) then
          call flusenseur_init_lamin_3dfull(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
           else
          call flusenseur_lamin_3dfull(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
           endif

      ELSEIF (option.eq.221) THEN

     
           if(nitcfg.eq.1) then

          call flusenseur_init_lamin_3dhomo(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
           else
          call flusenseur_lamin_3dhomo(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
           endif


      ELSEIF (option.eq.222) THEN

     
          if(nitcfg.le.1) then

          call flusenseur_init_lamin_3dcart(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
          else

            call flusenseur_lamin_3dcart(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
           endif


      ELSEIF (option.eq.223) THEN

     
          if(nitcfg.eq.1) then
          call flusenseur_init_lamin_2d(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
          else
          call flusenseur_lamin_2d(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
          endif

      ELSEIF (option.eq.230) THEN

     
           if(nitcfg.eq.1) then
          call flusenseur_init_SA_3dfull(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
           else
          call flusenseur_SA_3dfull(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
           endif

      ELSEIF (option.eq.231) THEN

     
           if(nitcfg.eq.1) then

          call flusenseur_init_SA_3dhomo(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
           else
          call flusenseur_SA_3dhomo(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
           endif


      ELSEIF (option.eq.232) THEN

     
          if(nitcfg.le.1) then

          call flusenseur_init_SA_3dcart(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
          else

            call flusenseur_SA_3dcart(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
           endif


      ELSEIF (option.eq.233) THEN

     
          if(nitcfg.eq.1) then
          call flusenseur_init_SA_2d(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
          else
          call flusenseur_SA_2d(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
          endif

!!
!!!!
!!!!  ROE minmod
!!!!
!!
      ELSEIF (option.eq.510) THEN

          call fluroe_minmod_euler_3dfull(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSEIF (option.eq.511) THEN

          call fluroe_minmod_euler_3dhomo(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSEIF (option.eq.512) THEN

          call fluroe_minmod_euler_3dcart(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSEIF(option.eq.513) THEN

          call fluroe_minmod_euler_2d(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSEIF (option.eq.520) THEN

          call fluroe_minmod_lamin_3dfull(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSEIF (option.eq.521) THEN

          call fluroe_minmod_lamin_3dhomo(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSEIF (option.eq.522) THEN

          call fluroe_minmod_lamin_3dcart(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSEIF (option.eq.523) THEN

          call fluroe_minmod_lamin_2d(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSEIF (option.eq.530) THEN

          call fluroe_minmod_SA_3dfull(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSEIF (option.eq.531) THEN

          call fluroe_minmod_SA_3dhomo(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSEIF (option.eq.532) THEN

          call fluroe_minmod_SA_3dcart(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSEIF (option.eq.533) THEN

          call fluroe_minmod_SA_2d(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, drodm, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)

      ELSE
         write(*,*) ' option = ' , option 
            write(*,*)'fluroe_minmod 3d non code'
           call error('gtfl3$',70,1)

      ENDIF
 
      end

