c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 40 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine invlu(ndo, nitcfg, nitrun, param_int, param_real,
     &     ind_loop_lu, mjrnewton,
     &     rotmp,rop_ssiter,
     &     drodm_in, drodm_out,
     &     ti,tj,tk,
     &     venti,ventj,ventk,
     &     coe, ssor, ssortmp)
c***********************************************************************
c     _U   USER : PECHIER 
c     _U   USER : DECK
c     ACT
c     _A      Factorisation LU
c     LU+SSOR pour Spalart Allmaras
c     VAL
c     _V      LCI + Jameson-Turkel
c     I/O
c     _/    rop_ssiter
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndo, nitcfg, nitrun, ind_loop_lu(6), param_int(0:*),
     &     mjrnewton

      REAL_E rotmp(*),rop_ssiter(*),coe(*),drodm_out(*), drodm_in(*)
      REAL_E ti(*),tj(*),tk(*),venti(*), ventj(*),ventk(*)
      REAL_E ssor(*), ssortmp(*)
      REAL_E param_real(0:*)

c     Var loc
      INTEGER_E lSA, lussor_end, i, k, j, lij, l,
     &     v2, v3, v4, v5, v6
      logical llower

#include "FastS/formule_param.h"
#include "FastS/formule_ssor_param.h"

! on blinde si pas assez de travail pour tous les threads

      if(ind_loop_lu(2).lt.ind_loop_lu(1)) return
      lSA           = 0
      if(param_int(IFLOW).eq.3.and.param_int(ILES).eq.0) lSA= 1

      if (param_int(NB_RELAX) == 1) then
         
! domaine fixe
         IF(param_int(LALE).eq.0) THEN

            if(lSA.eq.1) Then

!lower
               call invlu_l_SA(ndo, param_int, param_real, ind_loop_lu,
     &              drodm_in,drodm_out,rop_ssiter,
     &              ti,tj,tk,
     &              coe)
!Diag + upper + mjr_newton
               call invlu_u_SA(ndo, param_int, param_real,
     &              param_real(VISCO),param_real(SA_REAL),
     &              ind_loop_lu,
     &              drodm_out,rop_ssiter, rotmp,
     &              ti,tj,tk,
     &              coe, mjrnewton)

            else

!     lower
               call invlu_l(ndo, param_int, param_real, ind_loop_lu,
     &              drodm_in,drodm_out,rop_ssiter,
     &              ti,tj,tk,
     &              coe)
!     Diag + upper + mjr_newton
               call invlu_u(ndo, param_int, param_real, ind_loop_lu,
     &              drodm_out,rop_ssiter, rotmp,
     &              ti,tj,tk,
     &              coe, mjrnewton)
c     
            endif               !5 ou 6 Eq
         ELSE

            if(lSA.eq.1) Then

               call invlu_ale_l_SA(ndo, param_int, param_real, 
     &              ind_loop_lu,
     &              drodm_in,drodm_out,rop_ssiter,
     &              ti,tj,tk,venti,ventj,ventk,
     &              coe)
!     Diag + upper + mjr_newton
               call invlu_ale_u_SA(ndo, param_int, param_real,
     &              param_real(VISCO),param_real(SA_REAL),
     &              ind_loop_lu,
     &              drodm_out,rop_ssiter, rotmp,
     &              ti,tj,tk,venti,ventj,ventk,
     &              coe, mjrnewton)

            else

!     lower 
               call invlu_ale_l(ndo, param_int, param_real, ind_loop_lu,
     &              drodm_in,drodm_out,rop_ssiter,
     &              ti,tj,tk,venti,ventj,ventk,
     &              coe)

!     Diag + upper + mjr_newton
               call invlu_ale_u(ndo, param_int, param_real, ind_loop_lu,
     &              drodm_out,rop_ssiter, rotmp,
     &              ti,tj,tk,venti,ventj,ventk,
     &              coe, mjrnewton)


            endif               !5 ou 6 Eq
         ENDIF                  !maillge fixe/mobile

      elseif (param_int(NB_RELAX) .GE. 2) then

         v2 = param_int(NDIMDX)
         v3 = 2 * param_int(NDIMDX)
         v4 = 3 * param_int(NDIMDX)
         v5 = 4 * param_int(NDIMDX)
         v6 = 5 * param_int(NDIMDX)

         do k = ind_loop_lu(5), ind_loop_lu(6)
            do j = ind_loop_lu(3) - 1, ind_loop_lu(4) + 1

               lij = indssor(ind_loop_lu(1), j, k)

               do l = lij-1, lij + ind_loop_lu(2) - ind_loop_lu(1)+1

                  ssor(l) = 0.
                  ssor(l + v2) = 0.
                  ssor(l + v3) = 0.
                  ssor(l + v4) = 0.
                  ssor(l + v5) = 0.

               enddo
            enddo
         enddo

         do i = 1, param_int(NB_RELAX)

            lussor_end = 0
            if ((i == param_int(NB_RELAX)) .AND. (mjrnewton == 1)) then
               lussor_end = 1
            endif

#include "FastS/Compute/LU/lussor_mjrtmp.for"

!     domaine fixe
            IF(param_int(LALE).eq.0) THEN

               if(lSA.eq.1) Then

!     lower
                  call invlussor_l_SA(ndo, param_int, param_real, 
     &                 ind_loop_lu,
     &                 ssortmp,rop_ssiter,
     &                 ti,tj,tk,
     &                 coe, ssor)

#include "FastS/Compute/LU/lussor_mjrtmp.for"

!     Diag + upper + mjr_newton
                  call invlussor_u_SA(ndo, param_int, param_real,
     &                 param_real(VISCO),param_real(SA_REAL),
     &                 ind_loop_lu,
     &                 ssortmp,rop_ssiter, rotmp,
     &                 ti,tj,tk,
     &                 coe, ssor, lussor_end)

               else

                  call invlussor_l(ndo, param_int, param_real, 
     &                 ind_loop_lu,
     &                 ssortmp,rop_ssiter,
     &                 ti,tj,tk,
     &                 coe, ssor)

#include "FastS/Compute/LU/lussor_mjrtmp.for"

!     Diag + upper + mjr_newton
                  call invlussor_u(ndo, param_int, param_real,
     &                 ind_loop_lu,
     &                 ssortmp,rop_ssiter, rotmp,
     &                 ti,tj,tk,
     &                 coe, ssor, lussor_end)
c     
               endif            !5 ou 6 Eq
            ELSE

               if(lSA.eq.1) Then

                  call invlussor_ale_l_SA(ndo, param_int, param_real, 
     &                 ind_loop_lu,
     &                 ssortmp,rop_ssiter,
     &                 ti,tj,tk,venti,ventj,ventk,
     &                 coe, ssor)

#include "FastS/Compute/LU/lussor_mjrtmp.for"

!     Diag + upper + mjr_newton
                  call invlussor_ale_u_SA(ndo, param_int, param_real,
     &                 param_real(VISCO),param_real(SA_REAL),
     &                 ind_loop_lu,
     &                 ssortmp,rop_ssiter, rotmp,
     &                 ti,tj,tk,venti,ventj,ventk,
     &                 coe, ssor, lussor_end)

               else

!     lower 
                  call invlussor_ale_l(ndo, param_int, param_real,
     &                 ind_loop_lu,
     &                 ssortmp,rop_ssiter,
     &                 ti,tj,tk,venti,ventj,ventk,
     &                 coe, ssor)

#include "FastS/Compute/LU/lussor_mjrtmp.for"

!     Diag + upper + mjr_newton
                  call invlussor_ale_u(ndo, param_int, param_real, 
     &                 ind_loop_lu,
     &                 ssortmp,rop_ssiter, rotmp,
     &                 ti,tj,tk,venti,ventj,ventk,
     &                 coe, ssor, lussor_end)


               endif            !5 ou 6 Eq
            ENDIF               !maillge fixe/mobile

         enddo

      endif
      
      end
