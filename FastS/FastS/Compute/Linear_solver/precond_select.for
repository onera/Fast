c***********************************************************************
c     $Date: 2018-11-04 13:25:50 +0100 (Thu, 04 Nov 2018) $
c     $Revision: 64 $
c     $Author: Thibaut $
c***********************************************************************
      subroutine precond_select(param_int, param_real, ind_loop,
     &     step, nb_prec,
     &     rop, coe, ti, tj,
     &     vectin, vectout, ssor)
c***********************************************************************
c_U   USER : GUEGAN
c
c     ACT
c_A    selection preconditionneur a droite pour Krylov
c
c     VAL
c
c     I/O
c_/    vectin    : vecteur kryloc IN
c_/    vectout   : vecteur kryloc OUT
c
c***********************************************************************
      implicit none

      INTEGER_E ind_loop(6), param_int(0:*), nb_prec, step, nobl

      REAL_E param_real(0:*), vectin(*), vectout(*), coe(*),
     &     ti(*), tj(*), rop(*), ssor(*),norm


      if (nb_prec == 1) then
c
c         call precond(param_int, param_real, ind_loop, coe,
c     &        ti, tj, rop, vectin, vectout)
c
      elseif (nb_prec == 2) then
c
c         nobl = min(ind_loop(2), ind_loop(4))
c         call new_precond(param_int, param_real, ind_loop, step,
c     &        nobl, vectin, vectout, rop, coe, ti, tj)

      elseif (nb_prec == 3) then

         nobl = min(ind_loop(2), ind_loop(4))
         call LDURelaxScal(param_int, param_real, ind_loop,
     &                     step, nobl, 
     &                     vectin, vectout, rop, coe, ti, tj, ssor)

      endif

      end
