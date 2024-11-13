c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine cp_gradu(ndom, ithread, neq_grad,
     &                    param_int, c1,c2,
     &                    ind_sdm, ind_loop,
     &                    ind_dm, ijkv_cache,
     &                    synchro_send_th, synchro_receive_th,
     &                    icache, jcache, kcache,
     &                    rop, dvardc,
     &                    ti, tj, tk, vol)
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

      INTEGER_E ndom,neq_grad,ndimdx,ind_loop(6),ind_sdm(6),
     & ithread, icache, jcache, kcache,
     & ijkv_cache(3),ind_dm(6),
     & synchro_send_th(3), synchro_receive_th(3), param_int(0:*)

      REAL_E    rop( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E dvardc( param_int(NDIMDX) * neq_grad*3 )

      REAL_E vol( param_int(NDIMDX_MTR) )
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

       REAL_E tix,tjy,tkz,u,volinv

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

            !!calcul Grad*vol
            if(param_int(ITYPZONE).eq.0) then

               call cp_gradu_3dfull(ndom, ithread, neq_grad,
     &                        param_int, c1,c2,
     &                        ind_sdm,
     &                        ind_dm, ijkv_cache,
     &                        synchro_send_th, synchro_receive_th,
     &                        icache, jcache, kcache,
     &                        rop, dvardc,ti,tj,tk)
            
            elseif(param_int(ITYPZONE).eq.1) then

               call cp_gradu_3dhomo(ndom, ithread, neq_grad,
     &                        param_int, c1,c2,
     &                        ind_sdm,
     &                        ind_dm, ijkv_cache,
     &                        synchro_send_th, synchro_receive_th,
     &                        icache, jcache, kcache,
     &                        rop, dvardc,ti,tj,tk)
              
            elseif(param_int(ITYPZONE).eq.2) then

              
               call cp_gradu_3dcart(ndom, ithread, neq_grad,
     &                        param_int, c1,c2,
     &                        ind_sdm,
     &                        ind_dm, ijkv_cache,
     &                        synchro_send_th, synchro_receive_th,
     &                        icache, jcache, kcache,
     &                        rop, dvardc,ti,tj,tk)

            else
               call cp_gradu_2d(ndom, ithread, neq_grad,
     &                        param_int, c1,c2,
     &                        ind_sdm,
     &                        ind_dm, ijkv_cache,
     &                        synchro_send_th, synchro_receive_th,
     &                        icache, jcache, kcache,
     &                        rop, dvardc,ti,tj,tk)
            endif

CCC
CCC  On acheve le calcul gradient en divisant par le volume sur la list
CCC  de mise a jour
CCC
CCC

      Vdudx = 0
      Vdudy =   param_int(NDIMDX)
      Vdudz = 2*param_int(NDIMDX)
      Vdvdx = 3*param_int(NDIMDX)
      Vdvdy = 4*param_int(NDIMDX)
      Vdvdz = 5*param_int(NDIMDX)
      Vdwdx = 6*param_int(NDIMDX)
      Vdwdy = 7*param_int(NDIMDX)
      Vdwdz = 8*param_int(NDIMDX)

      if(param_int(ITYPZONE).eq.2) then

       volinv=1./vol(1)

#include "FastS/Compute/loop_begin.for"

          dvardc(l + Vdudx) = dvardc(l + Vdudx)*volinv
          dvardc(l + Vdudy) = dvardc(l + Vdudy)*volinv
          dvardc(l + Vdudz) = dvardc(l + Vdudz)*volinv
          dvardc(l + Vdvdx) = dvardc(l + Vdvdx)*volinv
          dvardc(l + Vdvdy) = dvardc(l + Vdvdy)*volinv
          dvardc(l + Vdvdz) = dvardc(l + Vdvdz)*volinv
          dvardc(l + Vdwdx) = dvardc(l + Vdwdx)*volinv
          dvardc(l + Vdwdy) = dvardc(l + Vdwdy)*volinv
          dvardc(l + Vdwdz) = dvardc(l + Vdwdz)*volinv
#include "FastS/Compute/loop_end.for"

      else

#include "FastS/Compute/loop_begin.for"

          volinv=1./vol(lvo)

          dvardc(l + Vdudx) = dvardc(l + Vdudx)*volinv
          dvardc(l + Vdudy) = dvardc(l + Vdudy)*volinv
          dvardc(l + Vdudz) = dvardc(l + Vdudz)*volinv
          dvardc(l + Vdvdx) = dvardc(l + Vdvdx)*volinv
          dvardc(l + Vdvdy) = dvardc(l + Vdvdy)*volinv
          dvardc(l + Vdvdz) = dvardc(l + Vdvdz)*volinv
          dvardc(l + Vdwdx) = dvardc(l + Vdwdx)*volinv
          dvardc(l + Vdwdy) = dvardc(l + Vdwdy)*volinv
          dvardc(l + Vdwdz) = dvardc(l + Vdwdz)*volinv
#include "FastS/Compute/loop_end.for"
 
      endif

      end

