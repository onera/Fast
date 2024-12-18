c***********************************************************************
c     $Date: 2013-04-30 16:01:06 +0200 (mar. 30 avril 2013) $
c     $Revision: 56 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine init_rhs(ndom, nitcfg, param_int, ndimdx, neq,
     &                   ind_loop,
     &                   drodm )
c***********************************************************************
c_P                          O N E R A
c
c_PR  PRODUIT :  mjdro3.sf/pechier1/e1/i1
c
c     ACT
c_A    Mise a jour de drodm apres calcul des flux.
c
c     INP
c_I    ndom   : numero-utilisateur du domaine
c_I    ipara  : direction du calcul des flux
c_I    neq    : nombre d equations du systeme
c_I    ndimdx : nombre de points maximal dans un domaine
c_I    flu    : flux aux interfaces
c
c     OUT
c
c     I/O
c_/    drodm  : increment des variables conservatives
c***********************************************************************
      implicit none

      REAL_E cxsi1, cxsi2
      parameter( cxsi1   = -0.6830127018922193 )
      parameter( cxsi2   = - 4./3.             )

#include "FastS/param_solver.h"

      INTEGER_E ndom, nitcfg, ndimdx, neq, param_int(0:*), ind_loop(6)

      REAL_E drodm( ndimdx * neq )
 

C Var loc
      INTEGER_E incmax, l, i,j,k,ne,lij,ltij,lt,b,ind, vg, lvo
      REAL_E ratio,coefH,xmut(1),rop(1) !!ajout pour feinter option de vecto 

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"


      b   = param_int(RK)
      ind = param_int(NSSITER)/param_int(LEVEL)
     
c***********************************************************************
Calcul des coeffcients : Explicit global
c***********************************************************************


      if (b.eq.3) then ! Runge-Kutta ordre 3
       
         if (nitcfg.eq.2) then
             coefH = cxsi1
         else if (nitcfg.eq.3) then
             coefH = cxsi2
         else if (nitcfg.eq.1) then
             coefH = 0.
         end if

      else if (b.eq.2) then ! Runge-Kutta ordre 2

         if (nitcfg.eq.2) then
            coefH =-1. 
         else if (nitcfg.eq.1) then 
            coefH = 0.
         end if

      else if (b.eq.1) then ! Runge-Kutta ordre 1

         if (nitcfg.eq.1) then
             coefH = 0.
         end if
      end if



      IF(param_int(ITYPCP).le.1.or.nitcfg.eq.1) THEN

          if(neq.eq.5) then
#include      "FastC/HPC_LAYER/loop_begin.for"
                 drodm(l          )=  0.
                 drodm(l +ndimdx  )=  0.
                 drodm(l +ndimdx*2)=  0.
                 drodm(l +ndimdx*3)=  0.
                 drodm(l +ndimdx*4)=  0.
#include      "FastC/HPC_LAYER/loop_end.for"
          else
#include      "FastC/HPC_LAYER/loop_begin.for"
                 drodm(l          )=  0.
                 drodm(l +ndimdx  )=  0.
                 drodm(l +ndimdx*2)=  0.
                 drodm(l +ndimdx*3)=  0.
                 drodm(l +ndimdx*4)=  0.
                 drodm(l +ndimdx*5)=  0.
#include      "FastC/HPC_LAYER/loop_end.for"
          endif

      ELSE


          if(neq.eq.5) then
#include      "FastC/HPC_LAYER/loop_begin.for"
                 drodm(l          )=  drodm(l          )*coefH
                 drodm(l +ndimdx  )=  drodm(l +ndimdx  )*coefH
                 drodm(l +ndimdx*2)=  drodm(l +ndimdx*2)*coefH
                 drodm(l +ndimdx*3)=  drodm(l +ndimdx*3)*coefH
                 drodm(l +ndimdx*4)=  drodm(l +ndimdx*4)*coefH
#include      "FastC/HPC_LAYER/loop_end.for"
          else
#include      "FastC/HPC_LAYER/loop_begin.for"
                 drodm(l          )=  drodm(l          )*coefH
                 drodm(l +ndimdx  )=  drodm(l +ndimdx  )*coefH
                 drodm(l +ndimdx*2)=  drodm(l +ndimdx*2)*coefH
                 drodm(l +ndimdx*3)=  drodm(l +ndimdx*3)*coefH
                 drodm(l +ndimdx*4)=  drodm(l +ndimdx*4)*coefH
                 drodm(l +ndimdx*5)=  drodm(l +ndimdx*5)*coefH
#include      "FastC/HPC_LAYER/loop_end.for"
          endif


      ENDIF

      end
