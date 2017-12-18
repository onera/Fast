c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine mjrvent_3dfull(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop1, ijkv_bloc, ijkv_cache,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 vtrans, rot, x,y,z, venti,ventj,ventk )
c***********************************************************************
c_U   USER : PECHIER
c
c     ACT
c_A    Appel du calcul metric corps mvr rigide
c
c     VAL
c_V    gaz parfait monoespece
c_V    processeur domaine
c_V    steady/unsteady
c
c     INP
c_I    ventijk     : vecteur normale aux facettes des mailles
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom, ithread, nptpsi,
     & icache, jcache, kcache,
     & ibloc, jbloc, kbloc, 
     & ijkv_bloc(3), ijkv_cache(3),ind_loop1(6),ind_dm(6),
     & synchro_send_sock(3),synchro_send_th(3),
     & synchro_receive_sock(3), synchro_receive_th(3), param_int(0:*)

      REAL_E x(param_int(NDIMDX_XYZ)),y(param_int(NDIMDX_XYZ)),
     &       z(param_int(NDIMDX_XYZ))

      REAL_E venti( param_int(NDIMDX_VENT)* param_int(NEQ_VENT))
      REAL_E ventj( param_int(NDIMDX_VENT)* param_int(NEQ_VENT))
      REAL_E ventk( param_int(NDIMDX_VENT)* param_int(NEQ_VENT))
      REAL_E rot(4,3), vtrans(3), param_real(0:*)

C Var loc
      INTEGER_E i,j,k,l,l1,l2,l3,lij,lxij,inci,incj,inck,icorr,jcorr,
     & kcorr,v1ven,v2ven,v3ven,lx1,lx2,lx3,lx4,l111,ind_loop(6)
      REAL_E cenx, ceny,cenz, cax,cay,caz


#include "FastS/formule_xyz_param.h"
#include "FastS/formule_vent_param.h"

      inci =1
      incj = param_int(NIJK_XYZ)
      inck = param_int(NIJK_XYZ)*param_int(NIJK_XYZ+1)

      icorr = 0 !correction indice boucle i pour traiter l'interface ind_loop(2)+1 si necessaire
      jcorr = 0 
      kcorr = 0 
      If(ibloc .eq.ijkv_bloc(1) .and.synchro_receive_sock(1).eq.0.and.
     &   icache.eq.ijkv_cache(1).and.synchro_receive_th(1).eq.0) icorr=1
      If(jbloc.eq.ijkv_bloc(2).and.synchro_receive_sock(2).eq.0.and.
     &   jcache.eq.ijkv_cache(2).and.synchro_receive_th(2).eq.0) jcorr=1
      If( kbloc.eq.ijkv_bloc(3).and.synchro_receive_sock(3).eq.0.and.
     &   kcache.eq.ijkv_cache(3).and.synchro_receive_th(3).eq.0) kcorr=1

      ind_loop = ind_loop1

      if(ind_loop1(1).eq.1) ind_loop(1)= 1-param_int(NIJK_VENT+3)
      if(ind_loop1(3).eq.1) ind_loop(3)= 1-param_int(NIJK_VENT+3)
      if(ind_loop1(5).eq.1) ind_loop(5)= 1-param_int(NIJK_VENT+4)
      if(icorr.eq.1) ind_loop(2)= ind_loop1(2)+param_int(NIJK_VENT+3)
      if(jcorr.eq.1) ind_loop(4)= ind_loop1(4)+param_int(NIJK_VENT+3)
      if(kcorr.eq.1) ind_loop(6)= ind_loop1(6)+param_int(NIJK_VENT+4)


      v1ven =   0
      v2ven =   param_int(NDIMDX_VENT)
      v3ven = 2*param_int(NDIMDX_VENT)

#include "FastS/Compute/pragma_align.for"

      DO k = ind_loop(5), ind_loop(6)
       DO j = ind_loop(3), ind_loop(4)

#include    "FastS/Compute/ALE/loopIvent_begin.for"
  
#include          "FastS/Compute/ALE/ventijk_3dfull_k.for"
            enddo

#include    "FastS/Compute/ALE/loopIvent_begin.for"
                
#include          "FastS/Compute/ALE/ventijk_3dfull_j.for"
            enddo

#include    "FastS/Compute/ALE/loopIvent_begin.for"
#include          "FastS/Compute/ALE/ventijk_3dfull_i.for"
            enddo

                if(icorr.eq.1) then !flux manquant en I

                   i = ind_loop(2) + 1
                   l = indven(  i, j, k)

#include           "FastS/Compute/ALE/ventijk_3dfull_i.for"
                endif ! 
       ENDDO !do j

       !Complement fluj en Jmax
       If(jcorr.eq.1) then

               j    = ind_loop(4)+1

#include   "FastS/Compute/ALE/loopIvent_begin.for"

#include           "FastS/Compute/ALE/ventijk_3dfull_j.for"
           enddo
        Endif

c        Endif  !Kmin ou pas

      ENDDO !do k


      !Complement fluk en Kmax
      If(kcorr.eq.1) then
           k    = ind_loop(6)+1
           do j = ind_loop(3),ind_loop(4)

#include   "FastS/Compute/ALE/loopIvent_begin.for"

#include        "FastS/Compute/ALE/ventijk_3dfull_k.for"
            enddo
          enddo
      Endif
      end

