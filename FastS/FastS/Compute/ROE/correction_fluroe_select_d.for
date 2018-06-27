c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 aout 2013) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine corr_fluroe_select_d(ndom, ithread, idir,
     &                        param_int, param_real,
     &                        ind_loop, 
     &                        rop, ropd, drodm, drodmd , wig,
     &                        venti, ventj, ventk,
     &                        ti,tj,tk,vol)
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

      INTEGER_E ndom, ithread, idir, ind_loop(6), param_int(0:*)


      REAL_E rop(*),drodm(*), ti(*),tj(*),tk(*),vol(*),
     & venti(*),ventj(*),ventk(*), wig(*)

      REAL_E param_real(0:*)

      REAL_E ropd(*),drodmd(*)
C Var loc
      INTEGER_E option, iflow_loc

      if(ind_loop(1).gt.ind_loop(2)) return 
      if(ind_loop(3).gt.ind_loop(4)) return 
      if(ind_loop(5).gt.ind_loop(6)) return

      iflow_loc = param_int(IFLOW)
      if(iflow_loc.eq.2) iflow_loc = 1

       option =1000*param_int(LALE)
     &        + 100*param_int(SLOPE)
     &        +  10*iflow_loc
     &        +      param_int(ITYPZONE)

       IF  (option.eq.223) THEN
                                               
         call corr_fluroe_lamin_o3_2d_d(ndom, ithread, idir,
     &                 param_int, param_real,
     &                 ind_loop,
     &                 rop, ropd, drodm, drodmd, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol)
                                               
       ELSEIF (option.eq.233) THEN
                                               
         call corr_fluroe_SA_o3_2d_d(ndom, ithread, idir,
     &                 param_int, param_real,
     &                 ind_loop,
     &                 rop, ropd, drodm, drodmd, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol)
                                               
       ELSEIF (option.eq.213) THEN
                                               
         call corr_fluroe_euler_o3_2d_d(ndom, ithread, idir,
     &                 param_int, param_real,
     &                 ind_loop,
     &                 rop, ropd, drodm, drodmd, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol)
                                               
      ELSE
         write(*,*) ' option = ' , option 
            write(*,*)'Unknown flux options'
           call error('correction_flu$',70,1)

      ENDIF
 
      end
