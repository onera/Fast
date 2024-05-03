c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine cp_gradu_3dfull(ndom, ithread, neq_grad,
     &                         param_int, c1,c2,
     &                         ind_loop,
     &                         ind_dm, ijkv_bloc, ijkv_cache,
     &                         synchro_send_sock, synchro_send_th,
     &                         synchro_receive_sock, synchro_receive_th,
     &                         ibloc , jbloc , kbloc ,
     &                         icache, jcache, kcache,
     &                         rop, dvardc,
     &                         ti, tj, tk)
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

      INTEGER_E ndom,neq_grad,ndimdx,ind_loop(6),
     & ithread, icache, jcache, kcache, ibloc, jbloc, kbloc, 
     & ijkv_bloc(3), ijkv_cache(3),ind_dm(6),
     & synchro_send_sock(3),synchro_send_th(3),
     & synchro_receive_sock(3), synchro_receive_th(3), param_int(0:*)

      REAL_E    rop( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E dvardc( param_int(NDIMDX) * neq_grad*3 )

      REAL_E  ti( param_int(NDIMDX_MTR) , param_int(NEQ_IJ) ),
     &        tj( param_int(NDIMDX_MTR) , param_int(NEQ_IJ) ),
     &        tk( param_int(NDIMDX_MTR) , param_int(NEQ_K ) )

      REAL_E c1,c2

C Var loc
      INTEGER_E inc,incmax,l,lt,i,j,k,incmax2,l5,l6,
     & l0,lt0,inci,incj,inck,ci,cj,lij,ltij,inci_mtr, incj_mtr,
     & inck_mtr,icorr,jcorr,ls,v1,v2,v3,v4,v5,
     & inc2i,inc2j,inc2k,lvo,
     & Vdudx,Vdudy,Vdudz, Vdvdx,Vdvdy,Vdvdz,Vdwdx,Vdwdy,Vdwdz

       REAL_E tix,tiy,tiz,tjx,tjy,tjz,tkx,tky,tkz,u

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

CC!DIR$ ASSUME_ALIGNED xmut: CACHELINE

      if(ind_loop(1).gt.ind_loop(2)) return 
      if(ind_loop(3).gt.ind_loop(4)) return 
      if(ind_loop(5).gt.ind_loop(6)) return


      inci = 1
      incj = param_int(NIJK)
      inck = param_int(NIJK)*param_int(NIJK+1)

      !fomule classique sur 2 point (bord domaine)
      if(c2.le.0.00001) then
       inc2i = 0
       inc2j = 0
       inc2k = 0
      else ! formule 4 points
       inc2i = 2*inci
       inc2j = 2*incj
       inc2k = 2*inck
      endif

      !metric
      lt  = indmtr(1 , 1, 1)
      tix = ti(lt,1)
      tjy = tj(lt,1)
      tkz = tk(lt,1)

      icorr = 0 !correction indice boucle i pour traiter l'interface ind_loop(2)+1 si necessaire
      jcorr = 0 
      If(ibloc .eq.ijkv_bloc(1) .and.synchro_receive_sock(1).eq.0.and.
     &   icache.eq.ijkv_cache(1).and.synchro_receive_th(1).eq.0) icorr=1
      If(jbloc.eq.ijkv_bloc(2).and.synchro_receive_sock(2).eq.0.and.
     &   jcache.eq.ijkv_cache(2).and.synchro_receive_th(2).eq.0) jcorr=1

      v1 = 0
      v2 =   param_int(NDIMDX)
      v3 = 2*param_int(NDIMDX)
      v4 = 3*param_int(NDIMDX)
      v5 = 4*param_int(NDIMDX)

      Vdudx = 0
      Vdudy =   param_int(NDIMDX)
      Vdudz = 2*param_int(NDIMDX)
      Vdvdx = 3*param_int(NDIMDX)
      Vdvdy = 4*param_int(NDIMDX)
      Vdvdz = 5*param_int(NDIMDX)
      Vdwdx = 6*param_int(NDIMDX)
      Vdwdy = 7*param_int(NDIMDX)
      Vdwdz = 8*param_int(NDIMDX)

#include "FastS/Compute/pragma_align.for"

      DO k = ind_loop(5), ind_loop(6)
       DO j = ind_loop(3), ind_loop(4)

#include        "FastS/POST/loopI_begin.for"

                  l0= l  - inck                    
#include          "FastS/POST/grad_k_3dfull.for"
                enddo

#include        "FastS/POST/loopI_begin.for"
                
                  l0= l  - incj  !on modifie rhs(j et j-1)
#include          "FastS/POST/grad_j_3dfull.for"
                enddo

#include       "FastS/POST/loopI_begin.for"
                
                  l0= l  - inci
#include          "FastS/POST/grad_i_3dfull.for"
                enddo

                if(icorr.eq.1) then !flux manquant en I

                   i = ind_loop(2) + 1
                   l = inddm( i, j, k)
                   ls= l -inci

#include          "FastS/POST/grad_i_3dfull_send.for"
                endif ! 
       ENDDO !do j

       !Complement fluj en Jmax
       If(jcorr.eq.1) then

               j    = ind_loop(4)+1

#include       "FastS/POST/loopI_begin.for"

                  ls = l -incj
#include          "FastS/POST/grad_j_3dfull_send.for"
               enddo
        Endif

c        Endif  !Kmin ou pas

      ENDDO !do k



      !Complement fluk en Kmax
      If( kbloc.eq.ijkv_bloc(3).and.synchro_receive_sock(3).eq.0.and.
     &   kcache.eq.ijkv_cache(3).and.synchro_receive_th(3).eq.0) then

           k    = ind_loop(6)+1
           do j = ind_loop(3),ind_loop(4)

#include    "FastS/POST/loopI_begin.for"

                 ls = l -inck
#include        "FastS/POST/grad_k_3dfull_send.for"
            enddo
          enddo
      Endif
      end

