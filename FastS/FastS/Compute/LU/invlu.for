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
     &     coe, ssor, ssor_size)
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
     &     drodm_in(param_int(NDIMDX), param_int(NEQ))
      REAL_E param_real(0:*)

c     Var loc
      INTEGER_E lSA, lussor_end, m, k, j, lij, ls, i, l, depth, mjr_ssor
      logical llower

#include "FastS/formule_param.h"
#include "FastS/formule_ssor_param.h"
! on blinde si pas assez de travail pour tous les threads

      if(ind_loop_lu(2).lt.ind_loop_lu(1)) return
      lSA           = 0
      depth         = 1
      if(param_int(IFLOW).eq.3.and.param_int(ILES).eq.0) lSA= 1

      IF (param_int(NB_RELAX) == 1) THEN
         
       !write(*,*)'lu', ind_loop_lu
       !write(*,*)'sdm',ind_loop_sdm
! domaine fixe
         If(param_int(LALE).eq.0) Then

            if(lSA.eq.1) Then

               !extrap coe(6) sur la 2eme rangee fictive
               if(param_int(LU_MATCH).eq.1) then
                 call extrap_coe(ndo, param_int, depth, 
     &                           ind_loop_lu,coe(1+5*param_int(NDIMDX)))

                 !lower
                 call invlu_overlap_l_SA(ndo, param_int, param_real, 
     &                   ind_loop_lu, ind_loop_sdm,
     &                   drodm_out,drodm_in, rop_ssiter,
     &                   ti,tj,tk,
     &                   coe, ssor_size)

                  !Diag + upper + mjr_newton
                  call invlu_overlap_u_SA(ndo, param_int, param_real,
     &                 param_real(VISCO),param_real(SA_REAL),
     &                 ind_loop_lu, ind_loop_sdm, mjrnewton,
     &                 drodm_out, rop_ssiter, rotmp,
     &                 ti,tj,tk,
     &                 coe, ssor_size)

               else !sans overlap

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
               endif

            else !5 Eqs


               if(param_int(LU_MATCH).eq.1) then
                 !lower
                 call invlu_overlap_l(ndo, param_int, param_real, 
     &                   ind_loop_lu, ind_loop_sdm,
     &                   drodm_out,drodm_in, rop_ssiter,
     &                   ti,tj,tk, coe, ssor_size)

                  !Diag + upper + mjr_newton
                  call invlu_overlap_u(ndo, param_int, param_real,
     &                             ind_loop_lu, ind_loop_sdm, mjrnewton,
     &                             drodm_out, rop_ssiter, rotmp,
     &                             ti,tj,tk, coe, ssor_size)
               else !sans overlap
                !lower
                 call invlu_l(ndo, param_int, param_real, ind_loop_lu,
     &                        drodm_in,drodm_out,rop_ssiter,
     &                        ti,tj,tk, coe)
                !Diag + upper + mjr_newton
                 call invlu_u(ndo, param_int, param_real, ind_loop_lu,
     &                        drodm_out,rop_ssiter, rotmp,
     &                        ti,tj,tk, coe, mjrnewton)
     
               endif ! overlap          
            endif    !5 ou 6 Eq

         ELSE   !ALE

            if(lSA.eq.1) Then

               if(param_int(LU_MATCH).eq.1) then
                 call extrap_coe(ndo, param_int, depth, 
     &                           ind_loop_lu,coe(1+5*param_int(NDIMDX)))

                 !lower
                 call invlu_overlap_ale_l_SA(ndo, param_int,param_real, 
     &                   ind_loop_lu, ind_loop_sdm,
     &                   drodm_out,drodm_in, rop_ssiter,
     &                   ti,tj,tk,venti,ventj,ventk,
     &                   coe, ssor_size)

                  !Diag + upper + mjr_newton
                  call invlu_overlap_ale_u_SA(ndo, param_int,param_real,
     &                 param_real(VISCO),param_real(SA_REAL),
     &                 ind_loop_lu, ind_loop_sdm, mjrnewton,
     &                 drodm_out, rop_ssiter, rotmp,
     &                 ti,tj,tk,venti,ventj,ventk,
     &                 coe, ssor_size)
               else
                 !lower
                 call invlu_ale_l_SA(ndo, param_int, param_real, 
     &                               ind_loop_lu,
     &                               drodm_in,drodm_out,rop_ssiter,
     &                               ti,tj,tk,venti,ventj,ventk,
     &                               coe)
                 !Diag + upper + mjr_newton
                 call invlu_ale_u_SA(ndo, param_int, param_real,
     &                            param_real(VISCO),param_real(SA_REAL),
     &                            ind_loop_lu,
     &                            drodm_out, rop_ssiter, rotmp,
     &                            ti,tj,tk,venti,ventj,ventk,
     &                            coe, mjrnewton)

               endif !!recouvrememnt
            else

               if(param_int(LU_MATCH).eq.1) then

                 !lower
                 call invlu_overlap_ale_l(ndo, param_int,param_real, 
     &                   ind_loop_lu, ind_loop_sdm,
     &                   drodm_out,drodm_in, rop_ssiter,
     &                   ti,tj,tk,venti,ventj,ventk,
     &                   coe, ssor_size)

                  !Diag + upper + mjr_newton
                  call invlu_overlap_ale_u(ndo, param_int,param_real,
     &                 ind_loop_lu, ind_loop_sdm, mjrnewton,
     &                 drodm_out, rop_ssiter, rotmp,
     &                 ti,tj,tk,venti,ventj,ventk,
     &                 coe, ssor_size)

               else
                 !lower 
                  call invlu_ale_l(ndo,param_int,param_real,ind_loop_lu,
     &                             drodm_in, drodm_out, rop_ssiter,
     &                             ti,tj,tk,venti,ventj,ventk, coe)

                 !Diag + upper + mjr_newton
                 call invlu_ale_u(ndo, param_int,param_real,ind_loop_lu,
     &                            drodm_out, rop_ssiter, rotmp,
     &                            ti,tj,tk,venti,ventj,ventk,
     &                            coe, mjrnewton)
               endif !!recouvrememnt

            endif               !5 ou 6 Eq
         ENDIF                  !maillge fixe/mobile

      elseif (param_int(NB_RELAX) .GE. 2) then


         mjr_ssor =1 

         !loop relaxation
         do m = 1, param_int(NB_RELAX)

            lussor_end = 0
            if (m.eq.param_int(NB_RELAX)) then

              mjr_ssor =0 
              if(mjrnewton.eq.1)lussor_end=1

            endif

            !domaine fixe
            IF(param_int(LALE).eq.0) THEN

               if(lSA.eq.1) Then

                  if(m.eq.1.and.param_int(LU_MATCH).eq.1) Then

                     call extrap_coe(ndo, param_int, depth, 
     &                           ind_loop_lu,coe(1+5*param_int(NDIMDX)))
                  endif

                  !lower
                  if(m.eq.1) then

                    call invlussor1_l_SA(ndo, param_int, param_real, 
     &                   ind_loop_lu, ind_loop_sdm,
     &                   drodm_out,drodm_in, rop_ssiter,
     &                   ti,tj,tk,
     &                   coe, ssor, ssor_size)
                  else

                    call invlussor_l_SA(ndo, param_int, param_real, 
     &                   ind_loop_lu, ind_loop_sdm,
     &                   drodm_out, drodm_in, rop_ssiter,
     &                   ti,tj,tk,
     &                   coe, ssor, ssor_size)
                  endif

                  
                  !Diag + upper + mjr_newton
                  call invlussor_u_SA(ndo, param_int, param_real,
     &                 param_real(VISCO),param_real(SA_REAL),
     &                 ind_loop_lu, ind_loop_sdm, mjr_ssor,
     &                 drodm_out, drodm_in, rop_ssiter, rotmp,
     &                 ti,tj,tk,
     &                 coe, ssor, lussor_end, ssor_size)

               else

                  !lower
                  if(m.eq.1) then

                  call invlussor1_l(ndo, param_int, param_real, 
     &                 ind_loop_lu, ind_loop_sdm,
     &                 drodm_out,drodm_in,rop_ssiter,
     &                 ti,tj,tk,
     &                 coe, ssor, ssor_size)
                  else
                  call invlussor_l(ndo, param_int, param_real, 
     &                 ind_loop_lu, ind_loop_sdm,
     &                 drodm_out, drodm_in, rop_ssiter,
     &                 ti,tj,tk,
     &                 coe, ssor, ssor_size)
                  endif

!#include "FastS/Compute/LU/lussor_mjrtmp.for"

!     Diag + upper + mjr_newton
                  call invlussor_u(ndo, param_int, param_real,
     &                 ind_loop_lu, ind_loop_sdm, mjr_ssor,
     &                 drodm_out, drodm_in, rop_ssiter, rotmp,
     &                 ti,tj,tk,
     &                 coe, ssor, lussor_end, ssor_size)
c     
               endif            !5 ou 6 Eq
            ELSE

               if(lSA.eq.1) Then

                 if(m.eq.1.and.param_int(LU_MATCH).eq.1) Then

                 call extrap_coe(ndo, param_int, depth, 
     &                           ind_loop_lu,coe(1+5*param_int(NDIMDX)))
                 endif

                  !lower
                  if(m.eq.1) then

                  call invlussor1_ale_l_SA(ndo, param_int, param_real, 
     &                 ind_loop_lu, ind_loop_sdm,
     &                 drodm_out,drodm_in,rop_ssiter,
     &                 ti,tj,tk,venti,ventj,ventk,
     &                 coe, ssor, ssor_size)
                  else
                  call invlussor_ale_l_SA(ndo, param_int, param_real, 
     &                 ind_loop_lu, ind_loop_sdm,
     &                 drodm_out, drodm_in, rop_ssiter,
     &                 ti,tj,tk,venti,ventj,ventk,
     &                 coe, ssor, ssor_size)
                  endif

!     Diag + upper + mjr_newton
                  call invlussor_ale_u_SA(ndo, param_int, param_real,
     &                 param_real(VISCO),param_real(SA_REAL),
     &                 ind_loop_lu, ind_loop_sdm, mjr_ssor,
     &                 drodm_out, drodm_in,rop_ssiter, rotmp,
     &                 ti,tj,tk,venti,ventj,ventk,
     &                 coe, ssor, lussor_end, ssor_size)

               else

                  !lower
                  if(m.eq.1) then
                  call invlussor1_ale_l(ndo, param_int, param_real,
     &                 ind_loop_lu, ind_loop_sdm,
     &                 drodm_out, drodm_in, rop_ssiter,
     &                 ti,tj,tk,venti,ventj,ventk,
     &                 coe, ssor, ssor_size)

                  else
                  call invlussor_ale_l(ndo, param_int, param_real,
     &                 ind_loop_lu, ind_loop_sdm,
     &                 drodm_out, drodm_in, rop_ssiter,
     &                 ti,tj,tk,venti,ventj,ventk,
     &                 coe, ssor, ssor_size)
                  endif


!     Diag + upper + mjr_newton
                  call invlussor_ale_u(ndo, param_int, param_real, 
     &                 ind_loop_lu, ind_loop_sdm, mjr_ssor,
     &                 drodm_out, drodm_in, rop_ssiter, rotmp,
     &                 ti,tj,tk,venti,ventj,ventk,
     &                 coe, ssor, lussor_end, ssor_size)


               endif            !5 ou 6 Eq
            ENDIF               !maillge fixe/mobile

         enddo !loop relaxation

      endif
      
      end
