c***********************************************************************
c     $Date: 2013-04-30 16:01:06 +0200 (mar. 30 avril 2013) $
c     $Revision: 56 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine init_rhs_dtloc(ndom, nitcfg, param_int, ndimdx, neq,
     &                   ind_loop,
     &                   drodm)
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

      REAL_E cxsi1, cxsi2,cxsi3,cxsi1_, cxsi2_,cxsi3_,cxsi4_,cxsi5_
      REAL_E cxsi6_,cxsi7_,cxsi8_,cxsi9_,cxsi10_,cxsi11_,cxsi12_,cxsi13_

      parameter( cxsi1   = -0.6830127018922193 )
      parameter( cxsi2   = - 4./3.             )
      parameter( cxsi3   = -0.96037658773652734 )


      parameter( cxsi1_   = -0.4178904745)
      parameter( cxsi2_   = -1.192151694643)
      parameter( cxsi3_   = -1.697784692471)
      parameter( cxsi4_   = -1.514183444257)



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

      if (MOD(nitcfg,ind)==1) then
            coefH=0.
            !print*,'coeff= ',coefH
      elseif(MOD(nitcfg,ind)==ind/2) then
            coefH = cxsi1
            !print*,'coeff= ',coefH
      elseif(MOD(nitcfg,ind)==ind-1) then
            coefH = cxsi2
            !print*,'coeff= ',coefH
      end if

c*****************************************************************************



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
