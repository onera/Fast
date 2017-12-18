c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 aout 2013) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine eff_flusenseur_init_select(ndom, ithread,
     &                        param_int, param_real,param_int_eff,
     &                        ind_loop, effort, pos,
     &                        rop, flu  , wig,
     &                        x, y, z,
     &                        venti, ventj, ventk,
     &                        ti,tj,tk,vol, xmut)
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

      INTEGER_E ndom, ithread,ind_loop(6), param_int(0:*),
     & param_int_eff(0:*)


      REAL_E rop(*),xmut(*),flu(*), ti(*),tj(*),tk(*),vol(*),
     & venti(*),ventj(*),ventk(*), wig(*),effort(*), pos(*),
     & x(*),y(*),z(*)

      REAL_E param_real(0:*)

C Var loc
      INTEGER_E option

      if(ind_loop(1).gt.ind_loop(2)) return 
      if(ind_loop(3).gt.ind_loop(4)) return 
      if(ind_loop(5).gt.ind_loop(6)) return


       !Si ALE et/ou SA, on force le passage par non ALE et NSLam
       !Si  Roe, on passe par AUSM pour avoir bon flux paroi

       !kflu_loc = param_int(KFLUDOM)
       !if(param_int(KFLUDOM).eq.5) kflu_loc =1

       option = 100*param_int(SLOPE)
     &        +  10*min(2,param_int(IFLOW))
     &        +      param_int(ITYPZONE)

       IF  (option.eq.220) THEN
                                               
           call eff_flusenseur_init_lamin_o3_3dfull(ndom, ithread, 
     &                 param_int, param_real,param_int_eff,
     &                 ind_loop, effort, pos,
     &                 rop, flu , wig,
     &                 x, y, z,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
                                               
       ELSEIF (option.eq.221) THEN
                                               
           call eff_flusenseur_init_lamin_o3_3dhomo(ndom, ithread, 
     &                 param_int, param_real,param_int_eff,
     &                 ind_loop, effort, pos,
     &                 rop, flu , wig,
     &                 x, y, z,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
                                               
       ELSEIF (option.eq.222) THEN
                                               
           call eff_flusenseur_init_lamin_o3_3dcart(ndom, ithread, 
     &                 param_int, param_real,param_int_eff,
     &                 ind_loop, effort, pos,
     &                 rop, flu , wig,
     &                 x, y, z,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
                                               
       ELSEIF (option.eq.223) THEN
                                               
           call eff_flusenseur_init_lamin_o3_2d(ndom, ithread, 
     &                 param_int, param_real,param_int_eff,
     &                 ind_loop, effort, pos,
     &                 rop, flu , wig,
     &                 x, y, z,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
                                               
       ELSEIF (option.eq.210) THEN
                                               
           call eff_flusenseur_init_euler_o3_3dfull(ndom, ithread, 
     &                 param_int, param_real,param_int_eff,
     &                 ind_loop, effort, pos,
     &                 rop, flu , wig,
     &                 x, y, z,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
                                               
       ELSEIF (option.eq.211) THEN
                                               
           call eff_flusenseur_init_euler_o3_3dhomo(ndom, ithread, 
     &                 param_int, param_real,param_int_eff,
     &                 ind_loop, effort, pos,
     &                 rop, flu , wig,
     &                 x, y, z,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
                                               
       ELSEIF (option.eq.212) THEN
                                               
           call eff_flusenseur_init_euler_o3_3dcart(ndom, ithread, 
     &                 param_int, param_real,param_int_eff,
     &                 ind_loop, effort, pos,
     &                 rop, flu , wig,
     &                 x, y, z,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
                                               
       ELSEIF (option.eq.213) THEN
                                               
           call eff_flusenseur_init_euler_o3_2d(ndom, ithread, 
     &                 param_int, param_real,param_int_eff,
     &                 ind_loop, effort, pos,
     &                 rop, flu , wig,
     &                 x, y, z,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
                                               
      ELSE
         write(*,*) ' option = ' , option 
            write(*,*)'Unknown flux options'
           call error('gtfl3$',70,1)

      ENDIF
 
      end
