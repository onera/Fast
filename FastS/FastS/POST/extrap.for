c***********************************************************************
c     $Date: 2011-10-10 16:10:53 +0200 (lun 10 oct 2011) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine extrap(ndom, param_int, c1, ind_loop, ind_dm, tab)
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

      INTEGER_E ndom, ind_loop(6), ind_dm(6), param_int(0:*)
c
      REAL_E    tab( param_int(NDIMDX) ), c1

c Var loc
      INTEGER_E i,j,k,l,ldjr,ific,shift, inck


#include "FastS/formule_param.h"

      shift = 0
      if(c1.eq.1.) shift=1

      inck  = 0

      if(param_int(ITYPZONE).ne.3) then

        inck  = 1

              if(    (ind_loop(5).eq.0.and.c1.eq.1.)
     &           .or.(ind_loop(5).eq.1.and.c1.ne.1.) ) then
              !extrap kmin
              k = 1-shift
                do  ific = 1,param_int(NIJK+4)-shift
                 do  j = ind_loop(3), ind_loop(4)
                   do  i = ind_loop(1), ind_loop(2)

                     ldjr      = inddm(i,  j, k      )
                     l         = inddm(i,  j, k -ific)
                      tab(l)   = tab(ldjr)
                  enddo
                 enddo
                enddo
              endif
              if(    (ind_loop(6).eq.ind_dm(6)+1.and.c1.eq.1.)
     &           .or.(ind_loop(6).eq.ind_dm(6)  .and.c1.ne.1.) ) then
              !extrap kmax
              k = ind_dm(6)+shift
                do  ific = 1,param_int(NIJK+4)-shift
                 do  j = ind_loop(3), ind_loop(4)
                   do  i = ind_loop(1), ind_loop(2)

                     ldjr      = inddm(i,  j, k     )
                     l         = inddm(i,  j, k+ific)
                      tab(l)   = tab(ldjr)
                  enddo
                 enddo
                enddo
              endif
      endif

      if(    (ind_loop(3).eq.0.and.c1.eq.1.)
     &   .or.(ind_loop(3).eq.1.and.c1.ne.1.) ) then
      !extrap jmin
          j = 1-shift
                do  ific = 1,param_int(NIJK+3)-shift
                 do  k = ind_loop(5)-inck, ind_loop(6)+inck
                   do  i = ind_loop(1), ind_loop(2)

                     ldjr      = inddm(i,  j     , k)
                     l         = inddm(i,  j-ific, k)
                      tab(l)   = tab(ldjr)
                  enddo
                 enddo
                enddo
      endif

      if(    (ind_loop(4).eq.ind_dm(4)+1.and.c1.eq.1.)
     &   .or.(ind_loop(4).eq.ind_dm(4)  .and.c1.ne.1.) ) then
      !extrap jmax
            j = ind_dm(4)+shift
                do  ific = 1,param_int(NIJK+3)-1
                 do  k = ind_loop(5)-inck, ind_loop(6)+inck
                   do  i = ind_loop(1), ind_loop(2)

                     ldjr      = inddm(i,  j     , k)
                     l         = inddm(i,  j+ific, k)
                      tab(l)   = tab(ldjr)
                  enddo
                 enddo
                enddo
      endif

      if(    (ind_loop(1).eq.0.and.c1.eq.1.)
     &   .or.(ind_loop(1).eq.1.and.c1.ne.1.) ) then
      !extrap imin
          i = 1-shift
                do  ific = 1,param_int(NIJK+3)-shift
                 do  k = ind_loop(5)-inck, ind_loop(6)+inck
                   do  j = ind_loop(3)-1, ind_loop(4)+1

                     ldjr      = inddm(  i      ,j, k)
                     l         = inddm(  i-ific ,j, k)
                      tab(l)   = tab(ldjr)
                  enddo
                 enddo
                enddo
      endif

      if(    (ind_loop(2).eq.ind_dm(2)+1.and.c1.eq.1.) 
     &   .or.(ind_loop(2).eq.ind_dm(2)  .and.c1.ne.1.) ) then
      !extrap imax
          i = ind_dm(2)+shift
                do  ific = 1,param_int(NIJK+3)-shift
                 do  k = ind_loop(5)-inck, ind_loop(6)+inck
                   do  j = ind_loop(3)-1, ind_loop(4)+1

                     ldjr      = inddm(  i      ,j, k)
                     l         = inddm(  i+ific ,j, k)
                      tab(l)   = tab(ldjr)
                  enddo
                 enddo
                enddo
      endif

      end
