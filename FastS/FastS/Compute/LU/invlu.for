c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 40 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine invlu(ndo, nitcfg, nitrun, param_int, param_real,
     &                 ind_loop_lu,
     &                 rotmp,rop_ssiter,
     &                 drodm,
     &                 ti,tj,tk,
     &                 venti,ventj,ventk,
     &                 coe)
c***********************************************************************
c_U   USER : PECHIER 
c_U   USER : DECK
c     ACT
c_A      Factorisation LU
c        LU+SSOR pour Spalart Allmaras
c     VAL
c_V      LCI + Jameson-Turkel
c     I/O
c_/    rop_ssiter
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndo, nitcfg, nitrun, ind_loop_lu(6), param_int(0:*)

      REAL_E rotmp(*),rop_ssiter(*),drodm(*),coe(*)
      REAL_E ti(*),tj(*),tk(*),venti(*), ventj(*),ventk(*)
      REAL_E param_real(0:*)

c Var loc
      INTEGER_E lSA
      logical llower

      ! on blinde si pas assez de travail pour tous les threads
      if(ind_loop_lu(2).lt.ind_loop_lu(1)) return

      lSA           = 0
      if(param_int(IFLOW).eq.3.and.param_int(ILES).eq.0) lSA= 1

      ! domaine fixe
      IF(param_int(LALE).eq.0) THEN

        if(lSA.eq.1) Then

        !lower
        call invlu_l_SA(ndo, param_int, param_real, ind_loop_lu,
     &                  drodm,rop_ssiter,
     &                  ti,tj,tk,
     &                  coe)
        !Diag + upper + mjr_newton
        call invlu_u_SA(ndo, param_int, param_real,
     &                  param_real(VISCO),param_real(SA_REAL),
     &                  ind_loop_lu,
     &                  drodm,rop_ssiter, rotmp,
     &                  ti,tj,tk,
     &                  coe)

        else

c        llower =.true. !flag pour parcourir la matrice L (lower) ou U (upper)
c        call invlu_lu(ndo,param_int(NEQ),param_int(NEQ_IJ),param_int(NEQ_K),param_int(NEQ_COE),
c     &                param_int(NDIMDX), param_int(NDIMDX_MTR), param_int(NIJK), param_int(NIJK_MTR),
c     &                param_int(ITYPZONE), param_real(GAMMA), param_real(CVINF),
c     &                ind_loop_lu,
c     &                llower,
c     &                drodm,rop_ssiter,
c     &                ti,tj,tk,
c     &                coe)
c        call invlu_d(ndo,param_int(NEQ),param_int(NEQ_COE),param_int(NDIMDX),param_int(NIJK),lSA,param_int(ITYPZONE),
c     &               ind_loop_lu,
c     &               drodm,
c     &               coe)
c       llower =.false. 
c       call invlu_lu(ndo,param_int(NEQ),param_int(NEQ_IJ),param_int(NEQ_K),param_int(NEQ_COE),
c     &                param_int(NDIMDX), param_int(NDIMDX_MTR), param_int(NIJK), param_int(NIJK_MTR),
c     &                param_int(ITYPZONE), param_real(GAMMA), param_real(CVINF),
c     &                ind_loop_lu,
c     &                llower,
c     &                drodm,rop_ssiter,
c     &                ti,tj,tk,
c     &                coe)
c       ! Mise a jour de la solution roptmp(P+1) = rop_ssiter(P) + w(idrodm)
c        call mjro_newton(ndo, nitcfg,  param_int(NEQ), param_int(NDIMDX), param_int(NIJK), ijkv,
c     &                   lSA, param_int(ITYPZONE), param_real(CVINF),  param_real(VISCO),param_real(RATIOM),
c     &                   ind_loop_lu,
c     &                   rotmp,rop_ssiter,
c     &                   drodm)
c
        !lower 
        call invlu_l(ndo, param_int, param_real, ind_loop_lu,
     &                drodm,rop_ssiter,
     &                ti,tj,tk,
     &                coe)

        !Diag + upper + mjr_newton
        call invlu_u(ndo, param_int, param_real, ind_loop_lu,
     &                drodm,rop_ssiter, rotmp,
     &                ti,tj,tk,
     &                coe)
c
        endif!5 ou 6 Eq
      ELSE

        if(lSA.eq.1) Then

        call invlu_ale_l_SA(ndo, param_int, param_real, ind_loop_lu,
     &                      drodm,rop_ssiter,
     &                      ti,tj,tk,venti,ventj,ventk,
     &                      coe)
        !Diag + upper + mjr_newton
        call invlu_ale_u_SA(ndo, param_int, param_real,
     &                      param_real(VISCO),param_real(SA_REAL),
     &                      ind_loop_lu,
     &                      drodm,rop_ssiter, rotmp,
     &                      ti,tj,tk,venti,ventj,ventk,
     &                      coe)



        else

        !lower 
        call invlu_ale_l(ndo, param_int, param_real, ind_loop_lu,
     &                drodm,rop_ssiter,
     &                ti,tj,tk,venti,ventj,ventk,
     &                coe)

        !Diag + upper + mjr_newton
        call invlu_ale_u(ndo, param_int, param_real, ind_loop_lu,
     &                drodm,rop_ssiter, rotmp,
     &                ti,tj,tk,venti,ventj,ventk,
     &                coe)


        endif !5 ou 6 Eq
      ENDIF !maillge fixe/mobile
c

c      call check(ndom,param_int(NEQ),param_int(NDIMDX),
c     &           param_int(NIJK)(1),param_int(NIJK)(1)*param_int(NIJK)(2),2,2,
c     &           w(itro_ssiter),'champ 1  p',1,-1,3,88,92,-1,3)
 
      end
