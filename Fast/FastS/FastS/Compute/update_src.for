c***********************************************************************
c     $Date: 2013-04-30 16:01:06 +0200 (mar. 30 avril 2013) $
c     $Revision: 56 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine update_src(ndom, nitcfg, param_int, param_real,
     &                   ind_loop, drodm, rop, coe, vol, ro_src)
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


#include "FastS/param_solver.h"

      INTEGER_E ndom, nitcfg,  param_int(0:*), ind_loop(6)

      REAL_E drodm( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E ro_src( param_int(NDIMDX) * (param_int(NEQ)+1) )
      REAL_E rop( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E coe( param_int(NDIMDX) * param_int(NEQ_COE) )
      REAL_E vol(param_int(NDIMDX_MTR))
      REAL_E param_real(0:*)
 

C Var loc
      INTEGER_E incmax, l, i,j,k,ne,lij,ltij,lt,b,ind, lvo
      INTEGER_E v1,v2,v3,v4,v5,v6,v7
      !REAL_E ratio,coefH,xmut(1),rop(1) !!ajout pour feinter option de vecto 
      REAL_E ratio,coefH,xmut(1), vol_dt,ec,ec_tg

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      v1 = 0
      v2 = param_int(NDIMDX)
      v3 = param_int(NDIMDX)*2
      v4 = param_int(NDIMDX)*3
      v5 = param_int(NDIMDX)*4
      v6 = param_int(NDIMDX)*5
      v7 = param_int(NDIMDX)*6

      IF(param_int(SRC).eq.1) THEN
        If(param_int(NEQ).eq.6) Then

#include "FastC/HPC_LAYER/loop_begin.for"

            !vol_dt     = 1./max(coe(l+v1),1.e-15)
            !vol_dt     = vol(lvo)
            vol_dt     = 1.

            drodm(l+v1)= drodm( l +v1) + ro_src( l +v1) *vol_dt
            drodm(l+v2)= drodm( l +v2) + ro_src( l +v2) *vol_dt
            drodm(l+v3)= drodm( l +v3) + ro_src( l +v3) *vol_dt
            drodm(l+v4)= drodm( l +v4) + ro_src( l +v4) *vol_dt
            drodm(l+v5)= drodm( l +v5) + ro_src( l +v5) *vol_dt
            drodm(l+v6)= drodm( l +v6) + ro_src( l +v6) *vol_dt
#include "FastC/HPC_LAYER/loop_end.for"
        Else
#include "FastC/HPC_LAYER/loop_begin.for"
            drodm(l +v1)=  drodm(l +v1) + ro_src(l +v1)
            drodm(l +v2)=  drodm(l +v2) + ro_src(l +v2)
            drodm(l +v3)=  drodm(l +v3) + ro_src(l +v3)
            drodm(l +v4)=  drodm(l +v4) + ro_src(l +v4)
            drodm(l +v5)=  drodm(l +v5) + ro_src(l +v5)
#include "FastC/HPC_LAYER/loop_end.for"
        Endif

      ELSE
        If(param_int(NEQ).eq.6) Then

          !ro_src( l +v7 ou v6): cellN_src
          !ro_src( l +v1): density_target
          !ro_src( l +v2): VelocityX_target
          !ro_src( l +v5): Temperature_target
          !ro_src( l +v6): nutildeSA_target
#include "FastC/HPC_LAYER/loop_begin.for"

            vol_dt     = 1./max(coe(l+v1),1.e-15)
            !vol_dt     = vol(lvo)
            !vol_dt     = 1.

            ec    = rop(l+v2)*rop(l+v2)+rop(l+v3)*rop(l+v3)
     &             +rop(l+v4)*rop(l+v4)
            ec_tg = ro_src(l+v2)*ro_src(l+v2)+ro_src(l+v3)*ro_src(l+v3)
     &             +ro_src(l+v4)*ro_src(l+v4)

            drodm(l+v1)= drodm( l +v1)+  ro_src( l +v7)*vol_dt
     &                                 *(ro_src( l +v1)-rop(l+v1) )
     &                                 *param_real(NUDGING_EQ1)
            drodm(l+v2)= drodm( l +v2)+  ro_src( l +v7)*vol_dt
     &                                 *(ro_src( l +v2)-rop(l+v2) )
     &                                 *param_real(NUDGING_EQ2)
     &                                 *rop(l+v1)
            drodm(l+v3)= drodm( l +v3)+  ro_src( l +v7)*vol_dt
     &                                 *(ro_src( l +v3)-rop(l+v3) )
     &                                 *param_real(NUDGING_EQ3)
     &                                 *rop(l+v1)
            drodm(l+v4)= drodm( l +v4)+  ro_src( l +v7)*vol_dt
     &                                 *(ro_src( l +v4)-rop(l+v4) )
     &                                 *param_real(NUDGING_EQ4)
     &                                 *rop(l+v1)
            drodm(l+v5)= drodm( l +v5)+  ro_src( l +v7)*vol_dt
     &                                 *(ec_tg-ec)
     &                                 *param_real(NUDGING_EQ5)
     &                                 *rop(l+v1)
            drodm(l+v6)= drodm( l +v6)+  ro_src( l +v7)*vol_dt
     &                                 *(ro_src( l +v6)-rop(l+v6) )
     &                                 *param_real(NUDGING_EQ6)
     &                                 *rop(l+v1)
#include "FastC/HPC_LAYER/loop_end.for"
        Else
#include "FastC/HPC_LAYER/loop_begin.for"

            vol_dt     = 1./max(coe(l+v1),1.e-15)
            !vol_dt     = param_real(NUDGING_AMPLI)/max(coe(l+v1),1.e-15)
            !vol_dt     = vol(lvo)
            !vol_dt     = 1.

            ec    = rop(l+v2)*rop(l+v2)+rop(l+v3)*rop(l+v3)
     &             +rop(l+v4)*rop(l+v4)
            ec_tg = ro_src(l+v2)*ro_src(l+v2)+ro_src(l+v3)*ro_src(l+v3)
     &             +ro_src(l+v4)*ro_src(l+v4)

            drodm(l+v1)= drodm( l +v1)+  ro_src( l +v6)*vol_dt
     &                                 *(ro_src( l +v1)-rop(l+v1) )
     &                                 *param_real(NUDGING_EQ1)
            drodm(l+v2)= drodm( l +v2)+  ro_src( l +v6)*vol_dt
     &                                 *(ro_src( l +v2)-rop(l+v2) )
     &                                 *param_real(NUDGING_EQ2)
     &                                 *rop(l+v1)
            drodm(l+v3)= drodm( l +v3)+  ro_src( l +v6)*vol_dt
     &                                 *(ro_src( l +v3)-rop(l+v3) )
     &                                 *param_real(NUDGING_EQ3)
     &                                 *rop(l+v1)
            drodm(l+v4)= drodm( l +v4)+  ro_src( l +v6)*vol_dt
     &                                 *(ro_src( l +v4)-rop(l+v4) )
     &                                 *param_real(NUDGING_EQ4)
     &                                 *rop(l+v1)
            drodm(l+v5)= drodm( l +v5)+  ro_src( l +v6)*vol_dt
     &                                 *(ec_tg-ec)
     &                                 *param_real(NUDGING_EQ5)
     &                                 *rop(l+v1)
#include "FastC/HPC_LAYER/loop_end.for"

        Endif!Neq
      ENDIF! source =1 ou 2

      end
