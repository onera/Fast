c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 aout 2013) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine mjr_ale_3dhomocart(ndom, param_int, param_real,
     &                        socket,  Nbre_socket,
     &                        ithread_sock, thread_parsock,
     &                        ind_dm_socket, socket_topology,
     &                        x,y,z,ti,tj,tk,ti0,tj0,tk0, vol,
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

      INTEGER_E ndom, socket, Nbre_socket, ithread_sock, thread_parsock,
     & ind_dm_socket(6), socket_topology(3), param_int(0:*)


      REAL_E ti(*),tj(*),tk(*),vol(*),venti(*),ventj(*),ventk(*),
     & ti0(*),tj0(*),tk0(*),x(*),y(*),z(*)

      REAL_E param_real(0:*)

C Var loc
      INTEGER_E translation_pur,lskip,ind_mtr(6),ind_coe(6),lmin
      REAL_E rot(4,3), vtrans(3)

      vtrans=0.
      !om modifie la position en fonction de la loi horaire
      call move( rot, param_real, translation_pur) 

      !!traitement metrique
      lskip = 0

      ind_mtr(1)= 1 
      ind_mtr(3)= 1 
      ind_mtr(5)= 1 
      if(param_int( ITYPZONE ).eq.1) then
           ind_mtr(2)   = param_int( IJKV   )
           ind_mtr(4)   = param_int( IJKV+1 )
           ind_mtr(6)   = 1
      elseif(param_int( ITYPZONE ).eq.2) then
           ind_mtr(2)   = 1
           ind_mtr(4)   = 1 
           ind_mtr(6)   = 1 
      else
          lskip =1
      endif

      IF(lskip.eq.0) THEN

         lmin = 10
         if(param_int(ITYPCP).eq.2) lmin = 4

         call indice_boucle_lu(ndom, socket, Nbre_socket, lmin, ind_mtr,
     &                         socket_topology, ind_dm_socket )

         call indice_boucle_lu(ndom, ithread_sock, thread_parsock, lmin,
     &                        ind_dm_socket, 
     &                        socket_topology, ind_coe )

         if(ind_coe(1).eq.ind_mtr(1)) 
     &    ind_coe(1)= ind_coe(1) -param_int( NIJK_MTR+3)
         if(ind_coe(3).eq.ind_mtr(3)) 
     &    ind_coe(3)= ind_coe(3) -param_int( NIJK_MTR+3)
         if(ind_coe(2).eq.ind_mtr(2)) 
     &    ind_coe(2)= ind_coe(2) +param_int( NIJK_MTR+3)
         if(ind_coe(4).eq.ind_mtr(4)) 
     &    ind_coe(4)= ind_coe(4) +param_int( NIJK_MTR+3)


          if(ind_coe(1).gt.ind_coe(2)) goto 100  
          if(ind_coe(3).gt.ind_coe(4)) goto 100
          if(ind_coe(5).gt.ind_coe(6)) goto 100

          !mvt de rotation: modific(ndom)ation de normale*surf
          if(translation_pur.ne.0) call cptijk_ale(ndom, param_int, 
     &                                             ind_coe, rot,
     &                                             ti0,tj0,tk0,ti,tj,tk)
      ENDIF

100   continue

      !!traitement Vitesse entrainement
      lskip = 0

       ind_mtr(1)=1
       ind_mtr(3)=1
       ind_mtr(5)=1
       if(param_int(ITYPVENT).eq.1.and.param_int(ITYPZONE).eq.1)then
           ind_mtr(2) = param_int( IJKV   )
           ind_mtr(4) = param_int( IJKV+1 )
           ind_mtr(6) = 1
       elseif(param_int(ITYPZONE).eq.2)then
           ind_mtr(2) = 1
           ind_mtr(4) = 1
           ind_mtr(6) = 1
       else
          lskip =1
       endif   

      IF(lskip.eq.0) THEN

         call indice_boucle_lu(ndom, socket, Nbre_socket, lmin, ind_mtr,
     &                         socket_topology, ind_dm_socket )

         call indice_boucle_lu(ndom, ithread_sock, thread_parsock, lmin,
     &                        ind_dm_socket, 
     &                        socket_topology, ind_coe )

         if(ind_coe(1).eq.ind_mtr(1)) 
     &    ind_coe(1)= ind_coe(1) -param_int( NIJK_VENT+3)
         if(ind_coe(3).eq.ind_mtr(3)) 
     &    ind_coe(3)= ind_coe(3) -param_int( NIJK_VENT+3)
         if(ind_coe(2).eq.ind_mtr(2)) 
     &    ind_coe(2)= ind_coe(2) +param_int( NIJK_VENT+3)
         if(ind_coe(4).eq.ind_mtr(4)) 
     &    ind_coe(4)= ind_coe(4) +param_int( NIJK_VENT+3)

          if(ind_coe(1).gt.ind_coe(2)) return  
          if(ind_coe(3).gt.ind_coe(4)) return
          if(ind_coe(5).gt.ind_coe(6)) return

         !mvt de rotation: modific(ndom)ation de normale*surf
         call cpventijk_ale(ndom, param_int, param_real, ind_coe, 
     &                      x,y,z, rot, venti, ventj, ventk )

      ENDIF

      end

