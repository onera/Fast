c***********************************************************************
c     $Date: 2010-12-22 11:30:40 +0100 (Wed, 22 Dec 2010) $
c     $Revision: 59 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine mj_wallmodel_plan(ndom,idir, param_int,
     &                            ind_fen,indth, inc_bc,
     &                            size_data, lprint, wmles_param,
     &                            rop,moy, snap)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom,idir,size_data,lprint
      INTEGER_E ind_fen(6), inc_bc(3), indth(6), param_int(0:*)

      REAL_E rop(param_int(NDIMDX), param_int(NEQ) )
      REAL_E moy( size_data       , param_int(NEQ) )
      REAL_E snap(size_data       , param_int(NEQ) )
      REAL_E wmles_param(4)

c  Var loc
      INTEGER_E i,j,k,l,nif,ninjf,ci,cj,ck,li,li1,
     &        njf,nkf,ic,jc,kc,l1,ind_loop(6),indbci,
     &        ldjrec,iplanrec,lrec,jrec,krec,irec
      REAL_E cnm,cn,nechant,c1,c2,c3,c4, vmoy1,vmoy2,vmoy3,vmoy4

#include "FastS/formule_param.h"

      indbci(j_1,k_1) = 1 + (j_1-inc_bc(2)) + (k_1-inc_bc(3))*inc_bc(1)


      !!mise a jour echantillon lund
      !wmles_param(1): N+1 plan stockes: [ moyenne, snap1, snap2,... snapN]
      !wmles_param(2): nbre de snap stock�s
      !wmles_param(3): nO du snap courant

      wmles_param(2) = wmles_param(2)+1
      nechant  = int( wmles_param(2))
      c4=1
      if(nechant.le.wmles_param(1) -1) c4=0
      wmles_param(2) = min(wmles_param(2), wmles_param(1) -1)

      !wmles_param(2) = min(wmles_param(2)+1, 1.)
      !wmles_param(2) = wmles_param(2)+1

      ! coeff pour moyenne glissante temporelle
      ! c1 = nbr echantillon N
      !!! nbr echantillon
      c1 = nechant
      cn = 1./c1
      c2 = (nechant-1)*cn
      c3 = 1./param_int(IJKV+2)


      !write(*,*)'c1,c2',c1,cn,c2

      !if(c1.eq.wmles_param(2)-1) c2=1.

      ind_loop(:) = ind_fen(:)
      if (idir.eq.1.or.idir.eq.2) then
          ind_loop(1) = param_int(WM_SAMPLING)
          ind_loop(2) = param_int(WM_SAMPLING)
      elseif (idir.eq.3.or.idir.eq.4) then
          ind_loop(3) = param_int(WM_SAMPLING)
          ind_loop(4) = param_int(WM_SAMPLING)
      else 
          ind_loop(5) = param_int(WM_SAMPLING)
          ind_loop(6) = param_int(WM_SAMPLING)
      endif

      !write(*,*)'fen',ind_fen
      !write(*,*)'lop',ind_loop
      !write(*,*)'inbc',inc_bc


      IF(   (indth(1).le.ind_loop(2)).and.(indth(2).ge.ind_loop(1))
     & .and.(indth(3).le.ind_loop(4)).and.(indth(4).ge.ind_loop(3)) 
     & .and.(indth(5).le.ind_loop(6)).and.(indth(6).ge.ind_loop(5)))
     &    THEN

       ind_loop(1)=max(ind_loop(1),indth(1))
       ind_loop(2)=min(ind_loop(2),indth(2))
       ind_loop(3)=max(ind_loop(3),indth(3))
       ind_loop(4)=min(ind_loop(4),indth(4))
       ind_loop(5)=max(ind_loop(5),indth(5))
       ind_loop(6)=min(ind_loop(6),indth(6))

      ELSE
        return
      ENDIF

      IF (idir.le.2) THEN


         do  k = ind_loop(5), ind_loop(6)
          do  j = ind_loop(3), ind_loop(4)

            l    = inddm( ind_loop(2), j,  k )
            li   = indbci(j,  k )

            moy(li,1) = rop(l,1)*cn + c2*moy(li,1)
            moy(li,2) = rop(l,2)*cn + c2*moy(li,2)
            moy(li,3) = rop(l,3)*cn + c2*moy(li,3)
            moy(li,4) = rop(l,4)*cn + c2*moy(li,4)
            moy(li,5) = rop(l,5)*cn + c2*moy(li,5)
           enddo
         enddo

       return


      ELSEIF (idir.le.4) THEN

         !calcul nouvelle moyenne et sauvegarde du snap
         k = 1
c          do  i = ind_loop(1), ind_loop(2)
c
c            li   = indbci(i,  k )
c
c            moy(li,1) = c2*moy(li,1) 
c            moy(li,2) = c2*moy(li,2)
c            moy(li,3) = c2*moy(li,3)
c            moy(li,4) = c2*moy(li,4)
c            moy(li,5) = c2*moy(li,5)
c           enddo


       if(lprint.eq.1) then
         do  k = ind_loop(5), ind_loop(6)
          do  i = ind_loop(1), ind_loop(2)

            l    = inddm( i, ind_loop(4),  k )
            li1  = indbci(i,  1 )
            li   = indbci(i,  k )
            !li   = indbci(i,  k )

            moy(li,1) = c2*moy(li,1) + (rop(l,1)-c4*snap(li,1))*cn
            moy(li,2) = c2*moy(li,2) + (rop(l,2)-c4*snap(li,2))*cn
            moy(li,3) = c2*moy(li,3) + (rop(l,3)-c4*snap(li,3))*cn
            moy(li,4) = c2*moy(li,4) + (rop(l,4)-c4*snap(li,4))*cn
            moy(li,5) = c2*moy(li,5) + (rop(l,5)-c4*snap(li,5))*cn
            !moy(li,1) = c2*moy(li,1) + (rop(l,1))*cn
            !moy(li,2) = c2*moy(li,2) + (rop(l,2))*cn
            !moy(li,3) = c2*moy(li,3) + (rop(l,3))*cn
            !moy(li,4) = c2*moy(li,4) + (rop(l,4))*cn
            !moy(li,5) = c2*moy(li,5) + (rop(l,5))*cn
            !moy(li1,1) = moy(li1,1) + (rop(l,1))*cn*c3
            !moy(li1,2) = moy(li1,2) + (rop(l,2))*cn*c3
            !moy(li1,3) = moy(li1,3) + (rop(l,3))*cn*c3
            !moy(li1,4) = moy(li1,4) + (rop(l,4))*cn*c3
            !moy(li1,5) = moy(li1,5) + (rop(l,5))*cn*c3

            snap(li,1) = rop(l,1)
            snap(li,2) = rop(l,2)
            snap(li,3) = rop(l,3)
            snap(li,4) = rop(l,4)
            snap(li,5) = rop(l,5)
           enddo
         enddo

         if(ndom.eq.1.and.lprint.eq.1) then
           i = 185
           k = 1
           l    = inddm( i, ind_loop(4),  k )
           li1  = indbci(i,  1 )
           write(*,'(a,8f15.10)')'mjr old new',moy(li1,2),
     &     moy(li1,3),moy(li1,4),rop(l,2),rop(l,3),rop(l,4),c2,cn

           vmoy1 =0
           vmoy2 =0
           vmoy3 =0
           vmoy4 =0
           do  k = ind_loop(5), ind_loop(6)
!DEC$ IVDEP
             do i = ind_loop(1), ind_loop(2) 

               l    = inddm( i ,  ind_loop(4)-2    , k )
               vmoy1 = vmoy1+rop(l,3)
               l    = inddm( i ,  ind_loop(4)-1    , k )
               vmoy2 = vmoy2+rop(l,3)
               l    = inddm( i ,  ind_loop(4)      , k )
               vmoy3 = vmoy3 +  rop(l,3)
               l    = inddm( i ,  ind_loop(4)+1    , k )
               vmoy4 = vmoy4+rop(l,3)
             enddo !i
            enddo !k
            l =(ind_loop(6)-ind_loop(5)+1)*(ind_loop(2)-ind_loop(1)+1)
            write(*,'(a,4f15.10)')'Vmoy',vmoy1/float(l),vmoy2/float(l),
     &                           vmoy3/float(l),vmoy4/float(l)
         endif
       endif!lprint

       vmoy3 =0
       do  k = ind_loop(5), ind_loop(6)
!DEC$ IVDEP
         do i = ind_loop(1), ind_loop(2) 
            l    = inddm( i ,  ind_loop(4)      , k )
            vmoy3 = vmoy3 +  rop(l,3)
         enddo !i
        enddo !k

        wmles_param(4)= vmoy3


c         vmoy =0.
c         k=1
c          do  i = ind_loop(1), ind_loop(2)
c            li   = indbci(i,  k )
c            vmoy = vmoy +  moy(li,3)
c          enddo
c         vmoy = vmoy/float(ind_loop(2))

c         wmles_param(3)= vmoy
 
c         do  k = ind_loop(5), ind_loop(6)
c          do  i = ind_loop(1), ind_loop(2)
c
c            li   = indbci(i,  k )
c            li1  = indbci(i,  1 )
c            l    = inddm( i, ind_loop(4),  k )
c            !l    = inddm( i, ind_loop(4)-1,  k )
c            !l1   = inddm( i, ind_loop(4)-2,  k )

            !rop(l,3) = rop(l,3) -vmoy
            !rop(l,3) = rop(l, 3) -vmoy
            !rop(l1,3)= rop(l1,3) -vmoy

c            moy(li,1) = moy(li1,1) 
c            moy(li,2) = moy(li1,2)
c            moy(li,3) = moy(li1,3) -vmoy
c            moy(li,4) = moy(li1,4)
c            moy(li,5) = moy(li1,5)

c           enddo
c         enddo

c         vmoy =0.
c         k=1
c          do  i = ind_loop(1), ind_loop(2)
c            li   = indbci(i,  k )
c            vmoy = vmoy +  moy(li,3)
c          enddo
c      write(*,'(a,8f20.12)')'verif moy',vmoy/float(ind_loop(2))

      ELSE

         do  j = ind_loop(3), ind_loop(4)
          do  i = ind_loop(1), ind_loop(2)

            l    = inddm( i, j, ind_loop(6) )
            li   = indbci(i,  j )

            moy(li,1) = rop(l,1)*cn + c2*moy(li,1)
            moy(li,2) = rop(l,2)*cn + c2*moy(li,2)
            moy(li,3) = rop(l,3)*cn + c2*moy(li,3)
            moy(li,4) = rop(l,4)*cn + c2*moy(li,4)
            moy(li,5) = rop(l,5)*cn + c2*moy(li,5)
           enddo
         enddo

      ENDIF

      end
