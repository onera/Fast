c***********************************************************************
c     $Date: 2010-04-06 11:08:16 +0200 (Tue, 06 Apr 2010) $
c     $Revision: 56 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine derdf_xyz(nijk_xyz, ndimt_xyz,
     &                     f0, f1,
     &                     iordre, ipara)
c***********************************************************************
c_P                          O N E R A
c
c_PR  PRODUIT :  der1_x.sf
c
c_DC  DATE_C : Novembre 2001 - M. Terracol 
c
c     HISTORIQUE
c
c     ACT
c_A    Derivation dans la direction x 
c
c     VAL
c
c     INP
c_I    ndom   : numero du domaine
c_I    f0     : "fonction" a deriver
c
c     OUT
c_O    f1     : derivee
c***********************************************************************
      implicit none

      INTEGER_E nijk_xyz(5), ndimt_xyz,iordre,ipara
      REAL_E    f0(ndimt_xyz), f1(ndimt_xyz)

C var loc
      INTEGER_E i, j, k, l,li,ci,cj,ck,ishift,
     &        deb, ifin, jfin, kfin, kdeb,
     &        indcr, indav, indap,iv,jv,kv,
     &        indav2, indap2,
     &        indav3, indap3,
     &        indav4, indap4,
     &        indav5, indap5,
     &        indav6, indap6,
     &        ind1,ind2,inc,ind_loop(6)

      REAL_E    c6,c1,c2,c3,c4,c5

#include "FastS/formule_xyz.h"

      !si dom2d et derive en k: on annule tableau sortie et on sort
      if(nijk_xyz(3).eq.2.and.ipara.eq.3) then

        f1(:) = 0.

        return
      endif

      ishift = iordre/2
      ci    = 0
      cj    = 0
      ck    = 0
      if(ipara.eq.1) then
        inc = 1
        ci    = ishift
      elseif(ipara.eq.2) then
        inc = nijk_xyz(1)
        cj  = ishift
      else
        inc = nijk_xyz(2)*nijk_xyz(1) 
        ck  = ishift
      endif

      iv = nijk_xyz(1)-2*nijk_xyz(4)-1
      jv = nijk_xyz(2)-2*nijk_xyz(4)-1
      kv = nijk_xyz(3)-2*nijk_xyz(5)-1

#ifndef E_SCALAR_COMPUTER
      deb   =  1 - nijk_xyz(4)
      kdeb  =  1 - nijk_xyz(5)
      ifin  = iv + nijk_xyz(4)
      jfin  = jv + nijk_xyz(4)
      kfin  = kv + nijk_xyz(5)

      ind1 = indcg( deb + ci, deb + cj, kdeb+ ck)
      ind2 = indcg( ifin- ci, jfin- cj, kfin- ck)
#else

      ind_loop(1) = 1  - nijk_xyz(4) + ci
      ind_loop(2) = iv + nijk_xyz(4) - ci
      ind_loop(3) = 1  - nijk_xyz(4) + cj
      ind_loop(4) = jv + nijk_xyz(4) - cj
      ind_loop(5) = 1  - nijk_xyz(5) + ck
      ind_loop(6) = kv + nijk_xyz(5) - ck
#endif


      c1 = 45./60.
      c2 = -9./60.
      c3 =  1./60.
           
        if (iordre.eq.2) then

#ifndef E_SCALAR_COMPUTER
CDIR$ IVDEP
!CDIR NODEP
CVD$  NODEPCHK
         do 10 l=ind1,ind2
#else
         do 10 k = ind_loop(5), ind_loop(6)
         do 10 j = ind_loop(3), ind_loop(4)
         do 10 i = ind_loop(1), ind_loop(2)
      
           l     = indcg(i,j,k)
#endif
           indav = l - inc
           indap = l + inc

           f1(l) = (f0(indap) - f0(indav))*c1
  10     continue

        elseif(iordre.eq.4) then

#ifndef E_SCALAR_COMPUTER
CDIR$ IVDEP
!CDIR NODEP
CVD$  NODEPCHK
         do 20 l=ind1,ind2
#else
         do 20 k = ind_loop(5), ind_loop(6)
         do 20 j = ind_loop(3), ind_loop(4)
         do 20 i = ind_loop(1), ind_loop(2)
      
           l     = indcg(i,j,k)
#endif
           indav = l     - inc
           indap = l     + inc
           indav2= indav - inc 
           indap2= indap + inc 

           f1(l) = (f0(indap ) - f0(indav ))*c1
     &            +(f0(indap2) - f0(indav2))*c2
 20      continue

        elseif(iordre.eq.6) then

#ifndef E_SCALAR_COMPUTER
CDIR$ IVDEP
!CDIR NODEP
CVD$  NODEPCHK
         do 30 l=ind1,ind2
#else
         do 30 k = ind_loop(5), ind_loop(6)
         do 30 j = ind_loop(3), ind_loop(4)
         do 30 i = ind_loop(1), ind_loop(2)
      
           l     = indcg(i,j,k)
#endif
           indav = l     - inc
           indap = l     + inc
           indav2= indav - inc 
           indap2= indap + inc
           indav3= indav2- inc 
           indap3= indap2+ inc

           f1(l) = (f0(indap ) - f0(indav ))*c1
     &            +(f0(indap2) - f0(indav2))*c2
     &            +(f0(indap3) - f0(indav3))*c3
 30      continue

        else 
         call error('derdf_xyz$',70,1)
        endif

      end
