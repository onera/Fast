c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 40 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine invlu(ndo, nitcfg, nitrun, param_int, param_real,
     &     ind_loop_lu, ind_loop_sdm, mjrnewton,
     &     rotmp,rop_ssiter,
     &     drodm_in, drodm_out,
     &     ti,tj,tk,
     &     venti,ventj,ventk,
     &     coe, ssor, ssortmp, ssor_size)
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
     &     mjrnewton, ssor_size, ind_loop_sdm(6)

      REAL_E rotmp(*),rop_ssiter(*),coe(*),drodm_out(*)
      REAL_E ti(*),tj(*),tk(*),venti(*), ventj(*),ventk(*)
      REAL_E ssor(ssor_size, param_int(NEQ)), 
     &     ssortmp(ssor_size, param_int(NEQ)),
     &     drodm_in(param_int(NDIMDX), param_int(NEQ))
      REAL_E param_real(0:*)

c     Var loc
      INTEGER_E lSA, lussor_end, m, k, j, lij, ls, i, l, depth
      logical llower

#include "FastS/formule_param.h"
#include "FastS/formule_ssor_param.h"
! on blinde si pas assez de travail pour tous les threads

      if(ind_loop_lu(2).lt.ind_loop_lu(1)) return
      lSA           = 0
      depth         = 1
      if(param_int(IFLOW).eq.3.and.param_int(ILES).eq.0) lSA= 1

      if (param_int(NB_RELAX) == 1) then
         
! domaine fixe
         IF(param_int(LALE).eq.0) THEN

            if(lSA.eq.1) Then

               !extrap coe(6) sur la 2eme rangee fictive
               if(param_int(LU_MATCH).eq.1) then
                 call extrap_coe(ndo, param_int, depth, 
     &                           ind_loop_lu,coe(1+5*param_int(NDIMDX)))
               endif
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

               if(param_int(LU_MATCH).eq.1) then
                 call extrap_coe(ndo, param_int, depth, 
     &                           ind_loop_lu,coe(1+5*param_int(NDIMDX)))
               endif

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

         i_size = ind_loop_sdm(2) - ind_loop_sdm(1) + 1 +
     &        2 * param_int(NIJK + 3) !taille de la fenetre + ghostcells
         j_size = ind_loop_sdm(4) - ind_loop_sdm(3) + 1 +
     &        2 * param_int(NIJK + 3)

#include "FastS/Compute/LU/lussor_initssor.for"

         do m = 1, param_int(NB_RELAX)
         !do m =1,1
            lussor_end = 0
            if ((m == param_int(NB_RELAX)) .AND. (mjrnewton == 1)) then
               lussor_end = 1
            endif

#include "FastS/Compute/LU/lussor_mjrtmp.for"

!     domaine fixe
            IF(param_int(LALE).eq.0) THEN

               if(lSA.eq.1) Then

               if(param_int(LU_MATCH).eq.1) then
                 call extrap_coe(ndo, param_int, depth, 
     &                           ind_loop_lu,coe(1+5*param_int(NDIMDX)))
               endif
!     lower
                  call invlussor_l_SA(ndo, param_int, param_real, 
     &                 ind_loop_lu, ind_loop_sdm,
     &                 ssortmp,rop_ssiter,
     &                 ti,tj,tk,
     &                 coe, ssor, ssor_size)

#include "FastS/Compute/LU/lussor_mjrtmp.for"

!     Diag + upper + mjr_newton
                  call invlussor_u_SA(ndo, param_int, param_real,
     &                 param_real(VISCO),param_real(SA_REAL),
     &                 ind_loop_lu, ind_loop_sdm,
     &                 ssortmp,rop_ssiter, rotmp,
     &                 ti,tj,tk,
     &                 coe, ssor, lussor_end, ssor_size)

               else

                  call invlussor_l(ndo, param_int, param_real, 
     &                 ind_loop_lu, ind_loop_sdm,
     &                 ssortmp,rop_ssiter,
     &                 ti,tj,tk,
     &                 coe, ssor, ssor_size)

#include "FastS/Compute/LU/lussor_mjrtmp.for"

!     Diag + upper + mjr_newton
                  call invlussor_u(ndo, param_int, param_real,
     &                 ind_loop_lu, ind_loop_sdm,
     &                 ssortmp,rop_ssiter, rotmp,
     &                 ti,tj,tk,
     &                 coe, ssor, lussor_end, ssor_size)
c     
               endif            !5 ou 6 Eq
            ELSE

               if(lSA.eq.1) Then

               if(param_int(LU_MATCH).eq.1) then
                 call extrap_coe(ndo, param_int, depth, 
     &                           ind_loop_lu,coe(1+5*param_int(NDIMDX)))
               endif

                  call invlussor_ale_l_SA(ndo, param_int, param_real, 
     &                 ind_loop_lu, ind_loop_sdm,
     &                 ssortmp,rop_ssiter,
     &                 ti,tj,tk,venti,ventj,ventk,
     &                 coe, ssor, ssor_size)

#include "FastS/Compute/LU/lussor_mjrtmp.for"

!     Diag + upper + mjr_newton
                  call invlussor_ale_u_SA(ndo, param_int, param_real,
     &                 param_real(VISCO),param_real(SA_REAL),
     &                 ind_loop_lu, ind_loop_sdm,
     &                 ssortmp,rop_ssiter, rotmp,
     &                 ti,tj,tk,venti,ventj,ventk,
     &                 coe, ssor, lussor_end, ssor_size)

               else

!     lower 
                  call invlussor_ale_l(ndo, param_int, param_real,
     &                 ind_loop_lu, ind_loop_sdm,
     &                 ssortmp,rop_ssiter,
     &                 ti,tj,tk,venti,ventj,ventk,
     &                 coe, ssor, ssor_size)

#include "FastS/Compute/LU/lussor_mjrtmp.for"

!     Diag + upper + mjr_newton
                  call invlussor_ale_u(ndo, param_int, param_real, 
     &                 ind_loop_lu, ind_loop_sdm,
     &                 ssortmp,rop_ssiter, rotmp,
     &                 ti,tj,tk,venti,ventj,ventk,
     &                 coe, ssor, lussor_end, ssor_size)


               endif            !5 ou 6 Eq
            ENDIF               !maillge fixe/mobile

         enddo

      endif
      
      end
