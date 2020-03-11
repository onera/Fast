c***********************************************************************
c     $Date: 2010-12-22 11:30:40 +0100 (Wed, 22 Dec 2010) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine shape_tab_mtr(neq_mtr, param_int, idir,
     &                         ic,jc,kc,kc_vent,
     &                         ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale)
c***********************************************************************
c_P                          O N E R A
c
c     ACT
c_A    determine la forme des tableuz metrique en fonction de la nature du domaine
c          2D, 3D homogene, #D general
c
c     VAL
c
c     INP
c_I    neq_mtr   : nombre de composant du tableau des normales
c
c     OUT
c=======================================================================
      implicit none

#include "FastP/param_solver.h"

      INTEGER_E neq_mtr, idir, ic,jc,kc,kc_vent, param_int(0:*)

      REAL_E ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale

c...ne nombre de normale varie selon le type de domaine

      if(neq_mtr.eq.3) then !Domaine 3D quelconque
       ic     = 1
       jc     = 2
       kc     = 3
       ci_mtr = 1.
       cj_mtr = 1.
       ck_mtr = 1.
      elseif(neq_mtr.eq.2) then !Domaine 3D avec une direction homogene k: traitement facette i et j
       ic     = 1
       jc     = 2
       kc     = 2
       kc_vent= 1
       ci_mtr = 1.
       cj_mtr = 1.
       ck_mtr = 0.
      elseif(neq_mtr.eq.1) then !Domaine 3D avec une direction homogene k ou dom cartesien:

        ic     = 1
        jc     = 1
        kc     = 1

        if(idir.le.2) then
          ci_mtr = 1.
          cj_mtr = 0.
          ck_mtr = 0.
        elseif(idir.le.4) then
          ci_mtr = 0.
          cj_mtr = 1.
          ck_mtr = 0.
        else
          ci_mtr = 0.
          cj_mtr = 0.
          ck_mtr = 1.
        endif
      
      else   !Domaine 2D 
       ic     = 1
       jc     = 1
       kc     = 1
       ci_mtr = 0.
       cj_mtr = 0.
       ck_mtr = 0.
      endif

      if(param_int(NEQ_VENT).ne.3) then
       ck_vent= 0.
       kc_vent= 1
      else
       ck_vent= 1
       kc_vent= 3
      endif

      c_ale = 1.*min(param_int(LALE),1)

      end
