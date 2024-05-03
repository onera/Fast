c***********************************************************************
c     $Date: 2013-02-04 19:27:20 +0100 (lun. 04 fevr. 2013) $
c     $Revision: 40 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine invlu_d(ndom,neq,neq_coe,ndimdx,nijk,lSA,imtr,
     &                   ind_loop,
     @                   drodm,
     @                   coe)
c***********************************************************************
c                              O N E R A
c
c_D   DATE_C/M : 1996 
c
c_U   USER : DARRACQ 
c
c     ACT
c_A    Inversion de la matrice diagonale  D-1
c
c     VAL
c_V
c
c     INP
c_I    ndom,neq,ndimdx
c
c     OUT
c
c     I/O
c_/    drodm
c
c***********************************************************************
      implicit none

      INTEGER_E ndom,neq,neq_coe,ndimdx,lSA,imtr,nijk(5),ind_loop(6)

      REAL_E drodm(ndimdx,neq),coe(ndimdx,neq_coe)

c Var loc
      INTEGER_E inck,l,lij,i,j,k

#include "FastS/formule.h"


      inck=nijk(1)*nijk(2)*nijk(5)

      IF(lSA.eq.1) THEN

        if(imtr.ne.3)then
#ifndef E_SCALAR_COMPUTER
CDIR$ IVDEP
!CDIR NODEP
CVD$  NODEPCHK
          do 10 l=1+inck,ndimdx-inck
#else
          do 10 k = ind_loop(5), ind_loop(6)
          do 10 j = ind_loop(3), ind_loop(4)
             lij  = inddm(ind_loop(1), j, k)
!DEC$ IVDEP
           do 10 l = lij, lij + ind_loop(2) - ind_loop(1)
#endif
c 
            drodm(l,1)=drodm(l,1)*coe(l,5)
            drodm(l,2)=drodm(l,2)*coe(l,5)
            drodm(l,3)=drodm(l,3)*coe(l,5)
            drodm(l,4)=drodm(l,4)*coe(l,5)
            drodm(l,5)=drodm(l,5)*coe(l,5)
            drodm(l,6)=drodm(l,6)*coe(l,6)
 
10    continue

        else
#ifndef E_SCALAR_COMPUTER
CDIR$ IVDEP
!CDIR NODEP
CVD$  NODEPCHK
          do 11 l=1+inck,ndimdx-inck
#else
          do 11 j = ind_loop(3), ind_loop(4)
             lij  = inddm(ind_loop(1), j, 1)
!DEC$ IVDEP
           do 11 l = lij, lij + ind_loop(2) - ind_loop(1)
#endif
            drodm(l,1)=drodm(l,1)*coe(l,5)
            drodm(l,2)=drodm(l,2)*coe(l,5)
            drodm(l,3)=drodm(l,3)*coe(l,5)
            drodm(l,5)=drodm(l,5)*coe(l,5)
            drodm(l,6)=drodm(l,6)*coe(l,6)

11    continue
        endif

      ELSE  !on traite 5 equation

        if(imtr.ne.3)then
#ifndef E_SCALAR_COMPUTER
CDIR$ IVDEP
!CDIR NODEP
CVD$  NODEPCHK
          do 100  l=1+inck,ndimdx-inck
#else
          do 100  k = ind_loop(5), ind_loop(6)
          do 100  j = ind_loop(3), ind_loop(4)
             lij  = inddm(ind_loop(1), j, k)
CCCCC!DEC$ IVDEP
           do 100 l = lij, lij + ind_loop(2) - ind_loop(1)
#endif
c     
            drodm(l,1)=drodm(l,1)*coe(l,5)
            drodm(l,2)=drodm(l,2)*coe(l,5)
            drodm(l,3)=drodm(l,3)*coe(l,5)
            drodm(l,4)=drodm(l,4)*coe(l,5)
            drodm(l,5)=drodm(l,5)*coe(l,5)
100   continue

        else
#ifndef E_SCALAR_COMPUTER
CDIR$ IVDEP
!CDIR NODEP
CVD$  NODEPCHK
          do 101 l=1+inck,ndimdx-inck
#else
          do 101 j = ind_loop(3), ind_loop(4)
             lij  = inddm(ind_loop(1), j, 1)
!DEC$ IVDEP
           do 101 l = lij, lij + ind_loop(2) - ind_loop(1)
#endif
            drodm(l,1)=drodm(l,1)*coe(l,5)
            drodm(l,2)=drodm(l,2)*coe(l,5)
            drodm(l,3)=drodm(l,3)*coe(l,5)
            drodm(l,5)=drodm(l,5)*coe(l,5)

101   continue
        endif
     
      ENDIF  !5 ou 6 equation

      end
