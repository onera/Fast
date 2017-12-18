c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 aout 2013) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine mjr_ale(ndom, nitcfg, ithread,
     &                        param_int, param_real,
     &                        ind_dm, ind_loop,ijkv_thread, ijkv_sdm,
     &                        synchro_send_sock, synchro_send_th,
     &                        synchro_receive_sock, synchro_receive_th,
     &                        ibloc , jbloc , kbloc ,
     &                        icache, jcache, kcache,
     &                        x,y,z,ti,ti0,tj,tj0,tk,tk0, vol,
     &                        venti, ventj, ventk)
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

      INTEGER_E ndom, nitcfg,ithread,
     & icache, jcache, kcache,ibloc, jbloc, kbloc,
     & ijkv_thread(3), ijkv_sdm(3), ind_loop(6),ind_dm(6),
     & synchro_send_sock(3),synchro_send_th(3),
     & synchro_receive_sock(3), synchro_receive_th(3), param_int(0:*)


      REAL_E ti(*),tj(*),tk(*),vol(*),venti(*),ventj(*),ventk(*),
     & ti0(*),tj0(*),tk0(*),x(*),y(*),z(*)

      REAL_E param_real(0:*)

C Var loc
      INTEGER_E translation_pur
      REAL_E rot(4,3), vtrans(3)

      if(ind_loop(1).gt.ind_loop(2)) return 
      if(ind_loop(3).gt.ind_loop(4)) return 
      if(ind_loop(5).gt.ind_loop(6)) return

       vtrans=0.
       !om modifie la position en fonction de la loi horaire
       call move( rot, param_real, translation_pur) 

       !mvt de rotation: modification de normale*surf
       if(param_int(ITYPZONE).eq.0) then

          call mjrtijk_3dfull(ndom, ithread,
     &                        param_int, param_real,
     &                    ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                    synchro_send_sock, synchro_send_th,
     &                    synchro_receive_sock, synchro_receive_th,
     &                    ibloc , jbloc , kbloc ,
     &                    icache, jcache, kcache,
     &                    rot, ti0,tj0,tk0,ti,tj,tk )

       elseif(param_int(ITYPZONE).eq.3) then

          call mjrtijk_2d(ndom, ithread,
     &                    param_int, param_real,
     &                    ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                    synchro_send_sock, synchro_send_th,
     &                    synchro_receive_sock, synchro_receive_th,
     &                    ibloc , jbloc , kbloc ,
     &                    icache, jcache, kcache,
     &                    rot, ti0,tj0,tk0,ti,tj,tk )
       else   !mvt metric traite avant zone paralel open, car shape metric imcompatible avec cache bloking: race....
        continue
       endif

       if(param_int(ITYPVENT).eq.0) then

          call mjrvent_3dfull(ndom, ithread,
     &                        param_int, param_real,
     &                    ind_dm, ind_loop, ijkv_thread, ijkv_sdm, 
     &                    synchro_send_sock, synchro_send_th,
     &                    synchro_receive_sock, synchro_receive_th,
     &                    ibloc , jbloc , kbloc ,
     &                    icache, jcache, kcache,
     &                    vtrans, rot, x,y,z,venti,ventj,ventk )

       elseif(param_int(ITYPVENT).eq.3) then

          call mjrvent_2d(ndom, ithread,
     &                    param_int, param_real,
     &                    ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                    synchro_send_sock, synchro_send_th,
     &                    synchro_receive_sock, synchro_receive_th,
     &                    ibloc , jbloc , kbloc ,
     &                    icache, jcache, kcache,
     &                    vtrans, rot, x,y,z,venti,ventj,ventk )

       else   ! Vent traite avant zone paralel open, car shape metric imcompatible avec cache bloking: race....
        continue
       endif



      end

