c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine mjrtijk_3dfull(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop1, ijkv_cache,
     &                 synchro_send_th, synchro_receive_th,
     &                 icache, jcache, kcache,
     &                 rot, ti0,tj0,tk0,ti,tj,tk )
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
c_I    tijk     : vecteur normale aux facettes des mailles
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom, ithread, nptpsi,
     & icache, jcache, kcache,
     & ijkv_cache(3),ind_loop1(6),ind_dm(6),
     & synchro_send_th(3), synchro_receive_th(3), param_int(0:*)


      REAL_E  ti( param_int(NDIMDX_MTR) * param_int(NEQ_IJ) ),
     &        tj( param_int(NDIMDX_MTR) * param_int(NEQ_IJ) ),
     &        tk( param_int(NDIMDX_MTR) * param_int(NEQ_K ) )
      REAL_E  ti0( param_int(NDIMDX_MTR) * param_int(NEQ_IJ) ),
     &        tj0( param_int(NDIMDX_MTR) * param_int(NEQ_IJ) ),
     &        tk0( param_int(NDIMDX_MTR) * param_int(NEQ_K ) )

      REAL_E rot(4,3),param_real(0:*)


C Var loc
      INTEGER_E i,j,k,lt,lt2,lt3,ltij,icorr,jcorr,kcorr,
     & v1mtr,v2mtr,v3mtr, ind_loop(6) 

#include "FastS/formule_mtr_param.h"


      icorr = 0 !correction indice boucle i pour traiter l'interface ind_loop(2)+1 si necessaire
      jcorr = 0 
      kcorr = 0 
      If(icache.eq.ijkv_cache(1).and.synchro_receive_th(1).eq.0) icorr=1
      If(jcache.eq.ijkv_cache(2).and.synchro_receive_th(2).eq.0) jcorr=1
      If(kcache.eq.ijkv_cache(3).and.synchro_receive_th(3).eq.0) kcorr=1

      ind_loop = ind_loop1

      if(ind_loop1(1).eq.1) ind_loop(1)= 1-param_int(NIJK_MTR+3)
      if(ind_loop1(3).eq.1) ind_loop(3)= 1-param_int(NIJK_MTR+3)
      if(ind_loop1(5).eq.1) ind_loop(5)= 1-param_int(NIJK_MTR+4)
      if(icorr.eq.1) ind_loop(2)= ind_loop1(2)+param_int(NIJK_MTR+3)
      if(jcorr.eq.1) ind_loop(4)= ind_loop1(4)+param_int(NIJK_MTR+3)
      if(kcorr.eq.1) ind_loop(6)= ind_loop1(6)+param_int(NIJK_MTR+4)

      v1mtr =   0
      v2mtr =   param_int(NDIMDX_MTR)
      v3mtr = 2*param_int(NDIMDX_MTR)

#include "FastS/Compute/pragma_align.for"

      DO k = ind_loop(5), ind_loop(6)
       DO j = ind_loop(3), ind_loop(4)

#include    "FastS/Compute/ALE/loopI_begin.for"

#include          "FastS/Compute/ALE/tijk_3dfull_k.for"
            enddo

#include    "FastS/Compute/ALE/loopI_begin.for"
                
#include          "FastS/Compute/ALE/tijk_3dfull_j.for"
            enddo

#include    "FastS/Compute/ALE/loopI_begin.for"
#include          "FastS/Compute/ALE/tijk_3dfull_i.for"
            enddo

                if(icorr.eq.1) then !flux manquant en I

                   i  = ind_loop(2) + 1
                   lt = indmtr(  i, j, k)

#include           "FastS/Compute/ALE/tijk_3dfull_i.for"
                endif ! 
       ENDDO !do j

       !Complement fluj en Jmax
       If(jcorr.eq.1) then

               j    = ind_loop(4)+1

#include   "FastS/Compute/ALE/loopI_begin.for"

#include           "FastS/Compute/ALE/tijk_3dfull_j.for"
           enddo
        Endif

c        Endif  !Kmin ou pas

      ENDDO !do k


      !Complement fluk en Kmax
      If(kcorr.eq.1) then
           k    = ind_loop(6)+1
           do j = ind_loop(3),ind_loop(4)

#include   "FastS/Compute/ALE/loopI_begin.for"

#include        "FastS/Compute/ALE/tijk_3dfull_k.for"
            enddo
          enddo
      Endif
      end

