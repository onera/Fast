c***********************************************************************
c     $Date: 2013-04-30 16:01:06 +0200 (mar. 30 avril 2013) $
c     $Revision: 56 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine update_src(ndom, nitcfg, param_int,
     &                   ind_loop,
     &                   drodm, ro_src)
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
      REAL_E ro_src( param_int(NDIMDX) * param_int(NEQ) )
 

C Var loc
      INTEGER_E incmax, l, i,j,k,ne,lij,ltij,lt,b,ind, lvo
      INTEGER_E v1,v2,v3,v4,v5,v6
      REAL_E ratio,coefH,xmut(1),rop(1) !!ajout pour feinter option de vecto 

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      v1 = 0
      v2 = param_int(NDIMDX)
      v3 = param_int(NDIMDX)*2
      v4 = param_int(NDIMDX)*3
      v5 = param_int(NDIMDX)*4
      v6 = param_int(NDIMDX)*5

      IF(param_int(NEQ).eq.6) THEN

          do  k = ind_loop(5), ind_loop(6)
          do  j = ind_loop(3), ind_loop(4)
#include    "FastS/Compute/loopI_begin.for"

           drodm(l +v1)=  drodm(l +v1) + ro_src(l +v1)
           drodm(l +v2)=  drodm(l +v2) + ro_src(l +v2)
           drodm(l +v3)=  drodm(l +v3) + ro_src(l +v3)
           drodm(l +v4)=  drodm(l +v4) + ro_src(l +v4)
           drodm(l +v5)=  drodm(l +v5) + ro_src(l +v5)
           drodm(l +v6)=  drodm(l +v6) + ro_src(l +v6)
        enddo
        enddo
        enddo

      ELSE

          do  k = ind_loop(5), ind_loop(6)
          do  j = ind_loop(3), ind_loop(4)
#include    "FastS/Compute/loopI_begin.for"

           drodm(l +v1)=  drodm(l +v1) + ro_src(l +v1)
           drodm(l +v2)=  drodm(l +v2) + ro_src(l +v2)
           drodm(l +v3)=  drodm(l +v3) + ro_src(l +v3)
           drodm(l +v4)=  drodm(l +v4) + ro_src(l +v4)
           drodm(l +v5)=  drodm(l +v5) + ro_src(l +v5)
        enddo
        enddo
        enddo

      ENDIF

      end
