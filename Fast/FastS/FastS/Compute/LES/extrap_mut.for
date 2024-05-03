c***********************************************************************
c     $Date: 2011-10-10 16:10:53 +0200 (lun 10 oct 2011) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine extrap_mut(ndom, param_int,depth, ind_loop,ind_dm,tab)
c***********************************************************************
c_U   USER : C. laurent
c
c     ACT
c_A    calcul grandeur moyenne pour bilan rij
c_A    commande de calcul des statistiques 
c
c     VAL
c
c     INP
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom,depth, ind_loop(6), ind_dm(6), param_int(0:*)
c
      REAL_E    tab( param_int(NDIMDX) )

c Var loc
      INTEGER_E i,j,k,l,ldjr,ific,shift
      INTEGER_E i1,i2,j1,j2,k1,k2


#include "FastS/formule_param.h"

      i1   =0
      i2   =0
      j1   =0
      j2   =0
      k1   =0
      k2   =0

      !extrap imin
      if( ind_loop(1).eq.1) then

          i = 1
          i1= depth

          do  ific = 1, depth 
            do  k = ind_loop(5), ind_loop(6)
            do  j = ind_loop(3), ind_loop(4)

               ldjr   = inddm(  i      ,j, k)
               l      = inddm(  i-ific ,j, k)
               tab(l) = tab(ldjr)
             enddo
             enddo
          enddo
      endif

      !extrap imax
      if(ind_loop(2).eq.ind_dm(2)) then

          i = ind_dm(2)
          i2= depth

          do  ific = 1, depth
            do  k = ind_loop(5), ind_loop(6)
            do  j = ind_loop(3), ind_loop(4)

               ldjr  = inddm(  i      ,j, k)
               l     = inddm(  i+ific ,j, k)
               tab(l)= tab(ldjr)
            enddo
            enddo
          enddo
      endif

      !extrap jmin
      if( ind_loop(3).eq.1) then

          j = 1
          j1= depth

          do  ific = 1, depth 
            do  k = ind_loop(5)   , ind_loop(6)
            do  i = ind_loop(1)-i1, ind_loop(2)+i2

               ldjr   = inddm(i,  j     , k)
               l      = inddm(i,  j-ific, k)
               tab(l) = tab(ldjr)
            enddo
            enddo
          enddo
      endif
      !extrap jmax
      if( ind_loop(4).eq.ind_dm(4)) then


          j = ind_dm(4)
          j2= depth

          do  ific = 1, depth 
            do  k = ind_loop(5)   , ind_loop(6)
            do  i = ind_loop(1)-i1, ind_loop(2)+i2

               ldjr   = inddm(i,  j     , k)
               l      = inddm(i,  j+ific, k)
               tab(l) = tab(ldjr)
            enddo
            enddo
          enddo
      endif


      !extrap kmin
      if( ind_loop(5).eq.1) then
          k = 1
          do  ific = 1, depth 
            do  j = ind_loop(3)-j1, ind_loop(4)+j2
            do  i = ind_loop(1)-i1, ind_loop(2)+i2

               ldjr   = inddm(i,  j, k      )
               l      = inddm(i,  j, k -ific)
               tab(l) = tab(ldjr)
            enddo
            enddo
          enddo
      endif
      !extrap kmax
      if(ind_loop(6).eq.ind_dm(6) ) then
          k = ind_dm(6)
          do  ific = 1, depth 
            do  j = ind_loop(3)-j1, ind_loop(4)+j2
            do  i = ind_loop(1)-i1, ind_loop(2)+i2

               ldjr   = inddm(i,  j, k     )
               l      = inddm(i,  j, k+ific)
               tab(l) = tab(ldjr)
            enddo
            enddo
          enddo
      endif

      end
