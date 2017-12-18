c***********************************************************************
c     $Date: 2013-04-30 16:01:06 +0200 (mar. 30 avril 2013) $
c     $Revision: 56 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine copynuma(ind_loop, ni,nj,shiftvar,ific,jfic,kfic,
     &                    cible, source )
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
c_I    param_int(NDIMDX) : nombre de points maximal dans un domaine
c_I    flu    : flux aux interfaces
c
c     OUT
c
c     I/O
c_/    drodm  : increment des variables conservatives
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ind_loop(6), ni,nj,shiftvar,ific,jfic,kfic

      REAL_E cible(*),source(*)
 

C Var loc
      INTEGER_E incmax, l, i,j,k,ne,lij,ltij,lt,lsrc,it
      REAL_E xmut(1),rop(1)  !!ajout pour feinter option de vecto 

         do  k = ind_loop(5), ind_loop(6)
         do  j = ind_loop(3), ind_loop(4)
         do  i = ind_loop(1), ind_loop(2)

               lsrc = (i+ific-1) + (j+jfic-1)*ni +(k+kfic-1)*ni*nj +1
               l    = lsrc +shiftvar

               cible(l) = source(lsrc)
        enddo
        enddo
        enddo

C        do it=1,3000
C
C         do  k = ind_loop(5), ind_loop(6)
C         do  j = ind_loop(3), ind_loop(4)
C         do  i = ind_loop(1), ind_loop(2)
C
C               lsrc = (i+ific-1) + (j+jfic-1)*ni +(k+kfic-1)*ni*nj +1
C               l    = lsrc +shiftvar
C          cible(l)=  cible(l)*1.0000001*it
CC
C        enddo
C        enddo
C        enddo
C
C       enddo

      end
