c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 59 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine cptst3(ndom,nitcfg,nitrun,first_it, lssiter_verif,
     &                  flagCellN, param_int, param_real,
     &                  ind_sdm, ind_grad, ind_coe,
     &                  cfl, xmut,rop, cellN, coe, ti,tj,tk, vol,venti)
c***********************************************************************
c_U   USER : PECHIER
c
c     ACT
c_A    Calcul du pas de temps 
c
c     VAL
c_V    Valide pour un G.P euler/vis/turb en STEADY/UNSTEADY
c
c     OUT
c_O    coe(:,1)   : tableau des pas de temps locaux  = dt/Vol
c_O    coe(:,2)   : tableau resolution LU
c_O    .....
c_O    coe(:,6)   : 
c

c     COM
c_C    STEADY   => dt nodal /dt fixe (common param_real(DTC))
c_C    UNSTEADY => dt fixe (common param_real(DTC))
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ind_sdm(6),ind_grad(6),ind_coe(6), param_int(0:*),
     & ndom,nitcfg,nitrun,lssiter_verif,first_it, flagCellN
 
      REAL_E xmut( param_int(NDIMDX) )
      REAL_E  rop( param_int(NDIMDX) , param_int(NEQ) )
      REAL_E  coe( param_int(NDIMDX) , param_int(NEQ_COE) )
      REAL_E ti( param_int(NDIMDX_MTR) , param_int(NEQ_IJ) ),
     &       tj( param_int(NDIMDX_MTR) , param_int(NEQ_IJ) ),
     &       tk( param_int(NDIMDX_MTR) , param_int(NEQ_K ) ),
     &      vol( param_int(NDIMDX_MTR) ) 

      REAL_E cellN( param_int(NDIMDX) )
      REAL_E venti( param_int(NDIMDX_VENT) * param_int(NEQ_VENT) )

      REAL_E cfl(3), param_real(0:*)

C Var loc
      INTEGER_E  l,i,j,k,lij,lt, ltij,lvo
      REAL_E volinv

#include "FastS/formule_mtr_param.h"
#include "FastS/formule_param.h"

      ! Calcul pas de temps local
      IF (param_int(DTLOC).eq.1) THEN

           if(param_int(ITYPCP).le.1.and.flagCellN.ne.1) then
 
           call tstb3_impli(ndom, first_it, param_int, param_real,
     &                      ind_coe ,
     &                      rop,coe,xmut,venti,cellN,
     &                      ti,tj,tk,vol)

           elseif(param_int(ITYPCP).le.1.and.flagCellN.eq.1) then

           call tstb3_impli_chimer(ndom,first_it, param_int, param_real,
     &                             ind_coe ,
     &                             rop,coe,xmut,venti,cellN,
     &                             ti,tj,tk,vol)

           elseif(param_int(ITYPCP).eq.2.and.flagCellN.ne.1) then

           call tstb3_expli(ndom, first_it, param_int, param_real,
     &                      ind_sdm ,
     &                      rop,coe,xmut,venti,cellN,
     &                      ti,tj,tk,vol)

           else
           call tstb3_expli_chimer(ndom,first_it, param_int, param_real,
     &                      ind_sdm ,
     &                      rop,coe,xmut,venti,cellN,
     &                      ti,tj,tk,vol)
           endif

      ! Calcul pas de temps global
      ELSE

       if(param_int(ITYPCP).le.1) then  !implicit

         !! calcul dt/vol
         call tstb3_global(ndom,flagCellN, param_int, param_real,
     &                     ind_coe,
     &                     coe, vol, cellN)

         !! calcul coef matrice LU
         call tstb3c( ndom, first_it, param_int, param_real,
     &                    ind_coe,
     &                    rop, coe, xmut, venti, ti, tj, tk, vol)
       else
         !! calcul dt/vol
         call tstb3_global(ndom,flagCellN, param_int, param_real,
     &                     ind_sdm,
     &                     coe, vol, cellN)
       endif

      ENDIF  !param_int(DTLOC) ou pas


      !calcul CFL en expli et implicit
      IF (nitcfg.eq.1.and.lssiter_verif.ne.0) then

            call cpcfl3(ndom , param_int, param_real,
     &                  ind_sdm,
     &                  cfl, rop,xmut,venti, ti,tj,tk,vol)
      ENDIF
      end
