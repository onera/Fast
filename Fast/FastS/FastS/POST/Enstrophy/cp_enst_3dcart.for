c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine cp_enst_3dcart(ndom, ithread, dim_grad,
     &                         param_int, c1,c2,
     &                         ind_loop,
     &                         ind_dm, ijkv_bloc, ijkv_cache,
     &                         synchro_send_sock, synchro_send_th,
     &                         synchro_receive_sock, synchro_receive_th,
     &                         ibloc , jbloc , kbloc ,
     &                         icache, jcache, kcache,
     &                         dvardc,
     &                         rop, enst, 
     &                         ti, tj, tk, vol)
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

      INTEGER_E ndom,neq_grad,ndimdx,ind_loop(6),dim_grad,
     & ithread, icache, jcache, kcache, ibloc, jbloc, kbloc, 
     & ijkv_bloc(3), ijkv_cache(3),ind_dm(6),
     & synchro_send_sock(3),synchro_send_th(3),
     & synchro_receive_sock(3), synchro_receive_th(3), param_int(0:*)

      REAL_E    rop( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E   enst( param_int(NDIMDX) )

      REAL_E  ti( param_int(NDIMDX_MTR) * param_int(NEQ_IJ) ),
     &        tj( param_int(NDIMDX_MTR) * param_int(NEQ_IJ) ),
     &        tk( param_int(NDIMDX_MTR) * param_int(NEQ_K ) ),
     &       vol( param_int(NDIMDX_MTR) )

      REAL_E dvardc(dim_grad*9)

      REAL_E c1,c2

C Var loc
      INTEGER_E inc,incmax,l,lt,i,j,k,incmax2,l5,l6,
     & l0,lt0,inci,incj,inck,ci,cj,lij,ltij,inci_mtr, incj_mtr,
     & inck_mtr,icorr,jcorr,ls,v1,v2,v3,v4,v5,
     & inc2i,inc2j,inc2k,incdir,v1mtr,v2mtr,v3mtr,
     & Vdudx,Vdvdx,Vdwdx, Vvar, Vgrad, lg, lvo


       REAL_E tcx,tcy,tcz,tcx1,tcy1,tcz1,u1,u3,volinv

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

CC!DIR$ ASSUME_ALIGNED xmut: CACHELINE

      if(ind_loop(1).gt.ind_loop(2)) return 
      if(ind_loop(3).gt.ind_loop(4)) return 
      if(ind_loop(5).gt.ind_loop(6)) return


      inci = 1
      incj = param_int(NIJK)
      inck = param_int(NIJK)*param_int(NIJK+1)

      inci_mtr = param_int(NIJK_MTR)
      incj_mtr = param_int(NIJK_MTR+1)
      inck_mtr = param_int(NIJK_MTR+2)

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
      tcx = ti(lt)
      tcy = tj(lt)
      tcz = tk(lt)

      volinv= 0.25/(vol(1)*vol(1))

      v1 = 0
      v2 =   param_int(NDIMDX)
      v3 = 2*param_int(NDIMDX)
      v4 = 3*param_int(NDIMDX)
      v5 = 4*param_int(NDIMDX)

      v1mtr =   0
      v2mtr =   param_int(NDIMDX_MTR)
      v3mtr = 2*param_int(NDIMDX_MTR)

      Vdudx = 0
      Vdvdx = 3*dim_grad
      Vdwdx = 6*dim_grad


#include "FastS/Compute/pragma_align.for"

      DO k = ind_loop(5), ind_loop(6)
       DO j = ind_loop(3), ind_loop(4)

          incdir=inc2k
#include "FastS/POST/loopI_begin.for"
           lg = l -lij
           l5    = l+inck
           l6    = l-inck
             
           !
           !contribution k 
           !
           !pour gradU
           Vgrad = Vdudx
           Vvar  = v2
#include   "FastS/POST/GradU/gradvar_3dcart_init.for"
           !pour gradV
           Vgrad = Vdvdx
           Vvar  = v3
#include   "FastS/POST/GradU/gradvar_3dcart_init.for"
           !pour gradW
           Vgrad = Vdwdx
           Vvar  = v4
#include   "FastS/POST/GradU/gradvar_3dcart_init.for"
         enddo

          incdir=inc2j
#include "FastS/POST/loopI_begin.for"
           lg = l -lij
           l5    = l+incj
           l6    = l-incj

           !
           !contribution j 
           !
           !pour gradU
           Vgrad = Vdudx
           Vvar  = v2
#include   "FastS/POST/GradU/gradvar_3dcart_j.for"
           !pour gradV
           Vgrad = Vdvdx
           Vvar  = v3
#include   "FastS/POST/GradU/gradvar_3dcart_j.for"
           !pour gradW
           Vgrad = Vdwdx
           Vvar  = v4
#include   "FastS/POST/GradU/gradvar_3dcart_j.for"
         enddo

          incdir=inc2i
#include"FastS/POST/loopI_begin.for"
           lg = l -lij
           l5    = l+inci
           l6    = l-inci

           !
           !contribution i 
           !
           !pour gradU
           Vgrad = Vdudx
           Vvar  = v2
#include   "FastS/POST/GradU/gradvar_3dcart_i.for"
           !pour gradV
           Vgrad = Vdvdx
           Vvar  = v3
#include   "FastS/POST/GradU/gradvar_3dcart_i.for"
           !pour gradW
           Vgrad = Vdwdx
           Vvar  = v4
#include   "FastS/POST/GradU/gradvar_3dcart_i.for"

#include   "FastS/POST/Enstrophy/enstrophy.for"

          !valeur interface green 2 fois trop grande(gradient auu carree
          !4 fois trop garnd
        enddo

       ENDDO !do j


      ENDDO !do k



      end

