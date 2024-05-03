c***********************************************************************
c     $Date: 2011-10-10 16:10:53 +0200 (lun 10 oct 2011) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine extrap_coe(ndom, param_int,depth, ind_loop,tab)
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

      INTEGER_E ndom,depth, ind_loop(6), param_int(0:*)
c
      REAL_E    tab( param_int(NDIMDX) )

c Var loc
      INTEGER_E i,j,k,l,ldjr,ific,shift
      INTEGER_E i1,i2,j1,j2,k1,k2,jdeb,jfin,kdeb,kfin


#include "FastS/formule_param.h"

      i1   =0
      i2   =0
      j1   =0
      j2   =0
      k1   =0
      k2   =0

      jdeb = max(ind_loop(3), 0                   )
      jfin = min(ind_loop(4), param_int(IJKV+1) +1)
      kdeb = max(ind_loop(5), 0                   )
      kfin = min(ind_loop(6), param_int(IJKV+2) +1)
      !extrap imin
      if( ind_loop(1).eq.1-param_int(NIJK+3) ) then

          i = 0
          i1= depth

          do  ific = 1, depth 
            do  k = kdeb, kfin
            do  j = jdeb, jfin

               ldjr   = inddm(  i      ,j, k)
               l      = inddm(  i-ific ,j, k)
               tab(l) = tab(ldjr)
             enddo
             enddo
          enddo
      endif

      !extrap imax
      i = param_int(IJKV) +1
      if(ind_loop(2).eq. i-1 + param_int(NIJK+3) ) then

          i2= depth

          do  ific = 1, depth
            do  k = kdeb, kfin
            do  j = jdeb, jfin

               ldjr  = inddm(  i      ,j, k)
               l     = inddm(  i+ific ,j, k)
               tab(l)= tab(ldjr)
            enddo
            enddo
          enddo
      endif

      !extrap jmin
      if( ind_loop(3).eq.1-param_int(NIJK+3) ) then

          j = 0
          j1= depth

          do  ific = 1, depth 
            do  k = kdeb, kfin
            do  i = ind_loop(1), ind_loop(2)

               ldjr   = inddm(i,  j     , k)
               l      = inddm(i,  j-ific, k)
               tab(l) = tab(ldjr)
            enddo
            enddo
          enddo
      endif
      !extrap jmax
      j = param_int(IJKV+1) +1
      if(ind_loop(4).eq. j-1 + param_int(NIJK+3) ) then

          j2= depth

          do  ific = 1, depth 
            do  k = kdeb, kfin
            do  i = ind_loop(1), ind_loop(2)

               ldjr   = inddm(i,  j     , k)
               l      = inddm(i,  j+ific, k)
               tab(l) = tab(ldjr)
            enddo
            enddo
          enddo
      endif


      if(param_int(ITYPZONE).ne.3) then

        !extrap kmin
        if( ind_loop(5).eq.1-param_int(NIJK+4) ) then
            k = 1
            do  ific = 1, depth 
              do  j = ind_loop(3), ind_loop(4)
              do  i = ind_loop(1), ind_loop(2)

                 ldjr   = inddm(i,  j, k      )
                 l      = inddm(i,  j, k -ific)
                 tab(l) = tab(ldjr)
              enddo
              enddo
            enddo
        endif
        !extrap kmax
        k = param_int(IJKV+2) +1
        if(ind_loop(6).eq. k-1 + param_int(NIJK+4) ) then
            do  ific = 1, depth 
              do  j = ind_loop(3), ind_loop(4)
              do  i = ind_loop(1), ind_loop(2)

                 ldjr   = inddm(i,  j, k     )
                 l      = inddm(i,  j, k+ific)
                 tab(l) = tab(ldjr)
              enddo
              enddo
            enddo
        endif
      endif !2d

c      !do  k = ind_loop(5), ind_loop(6)
c      do  k = 4,4
c       !do  j = ind_loop(3), ind_loop(4)
c       do  j = 0,0
c         do  i = ind_loop(1), ind_loop(2)
c                 l      = inddm(i,  j, k)
c                 write(*,*)tab(l),i,j,k
c         enddo
c       enddo
c      enddo

             
      end
