c***********************************************************************
c     $Date: 2010-12-22 11:30:40 +0100 (Wed, 22 Dec 2010) $
c     $Revision: 59 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine mj_lund_planrecyl(ndom,idir, param_int,
     &                             ind_fen,indth, inc_bc,
     &                             size_data, lund_param,
     &                             rop,rof2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom,idir,size_data
      INTEGER_E ind_fen(6), inc_bc(3), indth(6), param_int(0:*)

      REAL_E rop(param_int(NDIMDX), param_int(NEQ) )
      REAL_E rof2( size_data      , param_int(NEQ) )
      REAL_E lund_param(5)

c  Var loc
      INTEGER_E i,j,k,l,nif,ninjf,ci,cj,ck,li,lt,lvo,ltij,lij,
     &        njf,nkf,ic,jc,kc,l1,ind_loop(6),indbci,
     &        ldjrec,iplanrec,lrec,jrec,krec,irec
      REAL_E cnm,cn,nechant,c1,c2

#include "FastS/formule_mtr_param.h"
#include "FastS/formule_param.h"

      indbci(j_1,k_1) = 1 + (j_1-inc_bc(2)) + (k_1-inc_bc(3))*inc_bc(1)


      !!mise a jour echantillon lund
      lund_param(1) = lund_param(1) + 1.
      nechant  = lund_param(1)
      iplanrec = int(lund_param(2))


      ! coeff pour moyenne glissante temporelle
      ! c1 = nbr echantillon instant N
      ! c2 = nbr echantillon instant N-1
      !!! nbr echantillon
      c1 = nechant
      c2 = c1-1.

      cnm = c2/c1
      cn  = 1./c1

      ind_loop(:) = ind_fen(:)
      ci      = 0
      cj      = 0
      ck      = 0
      if (idir.eq.1.or.idir.eq.2) then
          ci  = 1
          ind_loop(1) = iplanrec
          ind_loop(2) = iplanrec
      elseif (idir.eq.3.or.idir.eq.4) then
          cj  = 1
          ind_loop(3) = iplanrec
          ind_loop(4) = iplanrec
      else 
          ck  = 1
          ind_loop(5) = iplanrec
          ind_loop(6) = iplanrec
      endif

      !write(*,*)'fen',ind_fen
      !write(*,*)'lop',ind_loop
      !write(*,*)'inbc',inc_bc


      IF(   (indth(1).le.ind_loop(2)).and.(indth(2).ge.ind_loop(1))
     & .and.(indth(3).le.ind_loop(4)).and.(indth(4).ge.ind_loop(3)) 
     & .and.(indth(5).le.ind_loop(6)).and.(indth(6).ge.ind_loop(5)))
     &    THEN

       ind_loop(1)=max(ind_loop(1),indth(1))
       ind_loop(2)=min(ind_loop(2),indth(2))
       ind_loop(3)=max(ind_loop(3),indth(3))
       ind_loop(4)=min(ind_loop(4),indth(4))
       ind_loop(5)=max(ind_loop(5),indth(5))
       ind_loop(6)=min(ind_loop(6),indth(6))

      ELSE
        return
      ENDIF

      IF (idir.le.2) THEN

         i = ind_loop(2)
#include "FastS/Compute/loopPlanI_begin.for"
            li   = indbci(j,  k )

            rof2(li,1) = rop(l,1)*cn + cnm*rof2(li,1)
            rof2(li,2) = rop(l,2)*cn + cnm*rof2(li,2)
            rof2(li,3) = rop(l,3)*cn + cnm*rof2(li,3)
            rof2(li,4) = rop(l,4)*cn + cnm*rof2(li,4)
            rof2(li,5) = rop(l,5)*cn + cnm*rof2(li,5)
#include "FastS/Compute/loopPlan_end.for"

      ELSEIF (idir.le.4) THEN
         j = ind_loop(4)
#include "FastS/Compute/loopPlanJ_begin.for"
            li   = indbci(i,  k )

            rof2(li,1) = rop(l,1)*cn + cnm*rof2(li,1)
            rof2(li,2) = rop(l,2)*cn + cnm*rof2(li,2)
            rof2(li,3) = rop(l,3)*cn + cnm*rof2(li,3)
            rof2(li,4) = rop(l,4)*cn + cnm*rof2(li,4)
            rof2(li,5) = rop(l,5)*cn + cnm*rof2(li,5)
#include "FastS/Compute/loopPlan_end.for"
      ELSE
         k = ind_loop(6)
#include "FastS/Compute/loopPlanK_begin.for"
            li   = indbci(i,  j )

            rof2(li,1) = rop(l,1)*cn + cnm*rof2(li,1)
            rof2(li,2) = rop(l,2)*cn + cnm*rof2(li,2)
            rof2(li,3) = rop(l,3)*cn + cnm*rof2(li,3)
            rof2(li,4) = rop(l,4)*cn + cnm*rof2(li,4)
            rof2(li,5) = rop(l,5)*cn + cnm*rof2(li,5)
#include "FastS/Compute/loopPlan_end.for"
      ENDIF

      end
