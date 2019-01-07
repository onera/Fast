c***********************************************************************
c     $Date: 2013-04-30 16:01:06 +0200 (mar. 30 avril 2013) $
c     $Revision: 56 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine init_tab(ndom, val, param_int, ndimdx, neq,
     &                   ind_loop,  drodm )
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

      INTEGER_E ndom, ndimdx, neq, param_int(0:*), ind_loop(6)

      REAL_E drodm( ndimdx * neq ), val
 

C Var loc
      INTEGER_E incmax, l, i,j,k,ne,lij,ltij,lt,b,ind, vg, lvo
      REAL_E ratio,coefH,xmut(1),rop(1) !!ajout pour feinter option de vecto 

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"


      do  ne = 1, neq
         vg = (ne-1)*ndimdx
         do  k = ind_loop(5), ind_loop(6)
         do  j = ind_loop(3), ind_loop(4)
#include      "FastS/Compute/loopI_begin.for"
                 drodm(l +vg      )= val
               enddo
         enddo
         enddo
      enddo


      end
