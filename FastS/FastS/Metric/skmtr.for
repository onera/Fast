C  
C    Copyright 2013-2019 Onera.
C
C    This file is part of Cassiopee.
C
C    Cassiopee is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    Cassiopee is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
C ============================================================================
      SUBROUTINE skmtr(ndom, param_int, param_real,
     &                 x, y, z, degener, dist,
     &                 ti, tj, tk,ti0, tj0, tk0, vol,venti,ventj,ventk,
     &                 ijkv_sdm,
     &                 ind_sdm, ind_mtr, ind_grad, 
     &                 ind_dm_zone, ind_dm_socket, ind_dm_omp,
     &                 socket_topology,
     &                 ithread_sock, thread_parsock,
     &                 Nbre_thread_actif, Nbre_socket, socket,
     &                 ithread)
C==============================================================================
C_IN
      IMPLICIT NONE

      include "omp_lib.h"
#include "FastS/param_solver.h"

      INTEGER_E ndom                                 
      INTEGER_E ithread, ijkv_sdm(3), ind_sdm(6),
     & ind_mtr(6),ind_grad(6), degener(*),
     & ind_dm_zone(6), ind_dm_omp(6), ind_dm_socket(6), 
     & socket_topology(3),Nbre_thread_actif,Nbre_socket,
     & thread_parsock, ithread_sock,
     & socket,param_int(0:*)

      REAL_E x(param_int( NDIMDX_XYZ )),y(param_int( NDIMDX_XYZ )),
     &       z(param_int( NDIMDX_XYZ )), param_real(0:*)
C_OUT
      REAL_E  ti(param_int( NDIMDX_MTR ),param_int( NEQ_IJ )) 
      REAL_E ti0(param_int( NDIMDX_MTR ),param_int( NEQ_IJ )) 
      REAL_E  tj(param_int( NDIMDX_MTR ),param_int( NEQ_IJ )) 
      REAL_E tj0(param_int( NDIMDX_MTR ),param_int( NEQ_IJ )) 
      REAL_E  tk(param_int( NDIMDX_MTR ),param_int( NEQ_K  )) 
      REAL_E tk0(param_int( NDIMDX_MTR ),param_int( NEQ_K  )) 
      REAL_E vol(param_int( NDIMDX_MTR ))

      REAL_E  dist( param_int( NDIMDX ) )

      REAL_E venti( param_int(NDIMDX_VENT) * param_int(NEQ_VENT) )
      REAL_E ventj( param_int(NDIMDX_VENT) * param_int(NEQ_VENT) )
      REAL_E ventk( param_int(NDIMDX_VENT) * param_int(NEQ_VENT) )

C_LOCAL
      logical l_initmtr
      INTEGER_E ind,nd,
     & i,j,k, ndo,thread_topology(3),ind_dm_thread(6),ind_loop(6),
     & imin,jmin,kmin,ind_rhs(6),ind_ven(6),
     & icache,jcache,kcache, ind0, lmtr,lmtr0, ne, l,l0,lmtrj,lmtr0j,
     & lmtri,lmtr0i,lmtrk,lmtr0k,li,lj,lk,ind1,ind2,ind3,no,
     & jmax,kmax,imax,translation_pur,lmin
      REAL_E ix,iy,iz,eps

#include "FastS/param_solver.h"
#include "FastS/formule_mtr_param.h"
#include "FastS/formule_xyz_param.h"
 
      l_initmtr=.false.
!!!
!on initialise les metriques en debut de calcul a partir du maillage lu dans init grid
!!!
      IF(.not.l_initmtr) THEN

        lmin = 10;
        if (param_int(ITYPCP).eq.2) lmin = 4

        call indice_boucle_lu(ndom, ithread_sock, thread_parsock, lmin,
     &                        ind_dm_socket, 
     &                        thread_topology, ind_dm_omp )

        ind_mtr = ind_dm_omp

        if(ind_dm_omp(1).eq.ind_dm_zone(1)) 
     &    ind_mtr(1)= ind_mtr(1) -param_int( NIJK+3)
        if(ind_dm_omp(3).eq.ind_dm_zone(3)) 
     &    ind_mtr(3)= ind_mtr(3) -param_int( NIJK+3)
        if(ind_dm_omp(5).eq.ind_dm_zone(5)) 
     &    ind_mtr(5)= ind_mtr(5) -param_int( NIJK+4)
        if(ind_dm_omp(2).eq.ind_dm_zone(2)) 
     &    ind_mtr(2)= ind_mtr(2) +param_int( NIJK+3)
        if(ind_dm_omp(4).eq.ind_dm_zone(4)) 
     &    ind_mtr(4)= ind_mtr(4) +param_int( NIJK+3)
        if(ind_dm_omp(6).eq.ind_dm_zone(6)) 
     &    ind_mtr(6)= ind_mtr(6) +param_int( NIJK+4)

           !tableau metric d'epaisseur k =1 si 3Dhomogene, cartesien ou 2D
           if(param_int( ITYPZONE ).ge.1) then
            ind_mtr(5)= 1
            ind_mtr(6)= 1
           endif
           if(param_int( ITYPZONE ).eq.2) then
            if(ithread.eq.1) then 
                ind_mtr(1:6) = 1
            else
                ind_mtr(2:6:2) = -2
            endif
           endif

        !if(ithread.eq.param_int( IO_THREAD)) then
        !   write(*,*)'thread_parsoc', thread_parsock ,Nbre_socket
        !   write(*,*)'zone',ind_dm_zone
        !   write(*,*)'sock',ind_dm_socket
        !   write(*,*)'mtr ',ind_mtr
        !endif

        !calcul normales facette I,J,K , surface et volume(DF et VF pour volume) from x,y,z
        call cp_tijk( param_int, x,y,z,ti,tj,tk,ti0,tj0,tk0, ind_mtr)

!$OMP BARRIER

!!!
!!! Calcul volume
!!! 
!!!
        !calcul volume(DF et VF pour volume) from x,y,z
        call cp_vol(param_int, x,y,z,ti,tj,tk,ti0,tj0,tk0, vol, ind_mtr)


!$OMP BARRIER
!! Barrier car besoin normale et volume OK pour extrapol

      !extrapol metric si zone pas cartesienne
      !! Bande OK, coin a verifier pour implict recouvert....
      if(param_int( ITYPZONE ).ne.2) then

!$OMP SINGLE

        call tijk_extrap(param_int( NDIMDX_MTR ),param_int(NDIMDX_XYZ ),
     &                    param_int( NIJK_XYZ ), param_int( NIJK_MTR ),
     &                    param_int( NEQ_IJ ),param_int( NEQ_K ),
     &                    ind_dm_zone,
     &                    degener,
     &                    ti,tj,tk, ti0,tj0,tk0, vol)

        eps =1e-11 
        if(param_int( ITYPZONE ).eq.0) then

          do k= ind_dm_zone(5), ind_dm_zone(6)
          do j= ind_dm_zone(3), ind_dm_zone(4)
          do i= ind_dm_zone(1), ind_dm_zone(2)

            l = indmtr(i  ,j  ,k)
            li= indmtr(i+1,j  ,k)
            lj= indmtr(i  ,j+1,k)
            lk= indmtr(i  ,j  ,k+1)

            ix =(ti(l,1)-ti(li,1))+(tj(l,1)-tj(lj,1))+(tk(l,1)-tk(lk,1))
            iy =(ti(l,2)-ti(li,2))+(tj(l,2)-tj(lj,2))+(tk(l,2)-tk(lk,2))
            iz =(ti(l,3)-ti(li,3))+(tj(l,3)-tj(lj,3))+(tk(l,3)-tk(lk,3))
          
            if(abs(ix).ge.eps.or.abs(iy).ge.eps.or.abs(iz).ge.eps) 
     &        write(*,'(a,3f25.20,4i5)')'cons',ix,iy,iz,i,j,k,ndom
          enddo 
          enddo 
          enddo

        elseif(param_int( ITYPZONE ).eq.1) then

          do k= ind_dm_zone(5), ind_dm_zone(6)
          do j= ind_dm_zone(3), ind_dm_zone(4)
          do i= ind_dm_zone(1), ind_dm_zone(2)

            l = indmtr(i  ,j  ,k)
            li= indmtr(i+1,j  ,k)
            lj= indmtr(i  ,j+1,k)
            lk= indmtr(i  ,j  ,k+1)

            ix =(ti(l,1)-ti(li,1))+(tj(l,1)-tj(lj,1))
            iy =(ti(l,2)-ti(li,2))+(tj(l,2)-tj(lj,2))
            iz =(tk(l,1)-tk(lk,1))
          
            if(abs(ix).ge.eps.or.abs(iy).ge.eps.or.abs(iz).ge.eps) 
     &        write(*,'(a,3f25.20,4i5)')'cons',ix,iy,iz,i,j,k,ndom

          enddo 
          enddo 
          enddo

       elseif(param_int( ITYPZONE ).eq.3) then

          do k= ind_dm_zone(5), ind_dm_zone(6)
          do j= ind_dm_zone(3), ind_dm_zone(4)
          do i= ind_dm_zone(1), ind_dm_zone(2)

            l = indmtr(i  ,j  ,k)
            li= indmtr(i+1,j  ,k)
            lj= indmtr(i  ,j+1,k)
            lk= indmtr(i  ,j  ,k+1)

            ix =(ti(l,1)-ti(li,1))+(tj(l,1)-tj(lj,1))
            iy =(ti(l,2)-ti(li,2))+(tj(l,2)-tj(lj,2))
          
            if(abs(ix).ge.eps.or.abs(iy).ge.eps) 
     &        write(*,'(a,3f25.20,4i5)')'cons',ix,iy,iz,i,j,k,ndom
          enddo 
          enddo 
          enddo
        endif

        !Face k
        do k= 0,  ind_dm_zone(6)+1, ind_dm_zone(6)+1
        do j= ind_dm_zone(3), ind_dm_zone(4)
        do i= ind_dm_zone(1), ind_dm_zone(2)

            l = indmtr(i  ,j  ,k)
            li= indmtr(i+1,j  ,k)
            lj= indmtr(i  ,j+1,k)
            lk= indmtr(i  ,j  ,k+1)

            if(param_int( ITYPZONE ).eq.0) then
            ix =(ti(l,1)-ti(li,1))+(tj(l,1)-tj(lj,1))+(tk(l,1)-tk(lk,1))
            iy =(ti(l,2)-ti(li,2))+(tj(l,2)-tj(lj,2))+(tk(l,2)-tk(lk,2))
            iz =(ti(l,3)-ti(li,3))+(tj(l,3)-tj(lj,3))+(tk(l,3)-tk(lk,3))
            elseif(param_int( ITYPZONE ).eq.1) then
            ix =(ti(l,1)-ti(li,1))+(tj(l,1)-tj(lj,1))
            iy =(ti(l,2)-ti(li,2))+(tj(l,2)-tj(lj,2))
            iz =(tk(l,1)-tk(lk,1))

            elseif(param_int( ITYPZONE ).eq.3) then
            ix =(ti(l,1)-ti(li,1))+(tj(l,1)-tj(lj,1))
            iy =(ti(l,2)-ti(li,2))+(tj(l,2)-tj(lj,2))
            iz =0.
            endif

            if(abs(ix).ge.eps.or.abs(iy).ge.eps.or.abs(iz).ge.eps) 
     &        write(*,'(a,3f25.20,4i3)')'consk',ix,iy,iz,i,j,k,ndom
        enddo 
        enddo 
        enddo

        !Face j
        do k= ind_dm_zone(5), ind_dm_zone(6)
        do j= 0,  ind_dm_zone(4)+1, ind_dm_zone(4)+1
        do i= ind_dm_zone(1), ind_dm_zone(2)
            l = indmtr(i  ,j  ,k)
            li= indmtr(i+1,j  ,k)
            lj= indmtr(i  ,j+1,k)
            lk= indmtr(i  ,j  ,k+1)

            if(param_int( ITYPZONE ).eq.0) then
            ix =(ti(l,1)-ti(li,1))+(tj(l,1)-tj(lj,1))+(tk(l,1)-tk(lk,1))
            iy =(ti(l,2)-ti(li,2))+(tj(l,2)-tj(lj,2))+(tk(l,2)-tk(lk,2))
            iz =(ti(l,3)-ti(li,3))+(tj(l,3)-tj(lj,3))+(tk(l,3)-tk(lk,3))
            elseif(param_int( ITYPZONE ).eq.1) then
            ix =(ti(l,1)-ti(li,1))+(tj(l,1)-tj(lj,1))
            iy =(ti(l,2)-ti(li,2))+(tj(l,2)-tj(lj,2))
            iz =(tk(l,1)-tk(lk,1))
            elseif(param_int( ITYPZONE ).eq.3) then
            ix =(ti(l,1)-ti(li,1))+(tj(l,1)-tj(lj,1))
            iy =(ti(l,2)-ti(li,2))+(tj(l,2)-tj(lj,2))
            iz =0.
            endif

            if(abs(ix).ge.eps.or.abs(iy).ge.eps.or.abs(iz).ge.eps) 
     &        write(*,'(a,3f25.20,4i5)')'consj',ix,iy,iz,i,j,k,ndom

        enddo 
        enddo 
        enddo

        !Face i
        do k= ind_dm_zone(5), ind_dm_zone(6)
        do j= ind_dm_zone(3), ind_dm_zone(4)
        do i= 0,  ind_dm_zone(2)+1, ind_dm_zone(2)+1
            l = indmtr(i  ,j  ,k)
            li= indmtr(i+1,j  ,k)
            lj= indmtr(i  ,j+1,k)
            lk= indmtr(i  ,j  ,k+1)

            if(param_int( ITYPZONE ).eq.0) then
            ix =(ti(l,1)-ti(li,1))+(tj(l,1)-tj(lj,1))+(tk(l,1)-tk(lk,1))
            iy =(ti(l,2)-ti(li,2))+(tj(l,2)-tj(lj,2))+(tk(l,2)-tk(lk,2))
            iz =(ti(l,3)-ti(li,3))+(tj(l,3)-tj(lj,3))+(tk(l,3)-tk(lk,3))
            elseif(param_int( ITYPZONE ).eq.1) then
            ix =(ti(l,1)-ti(li,1))+(tj(l,1)-tj(lj,1))
            iy =(ti(l,2)-ti(li,2))+(tj(l,2)-tj(lj,2))
            iz =(tk(l,1)-tk(lk,1))
            elseif(param_int( ITYPZONE ).eq.3) then
            ix =(ti(l,1)-ti(li,1))+(tj(l,1)-tj(lj,1))
            iy =(ti(l,2)-ti(li,2))+(tj(l,2)-tj(lj,2))
            iz =0
            endif

            if(abs(ix).ge.eps.or.abs(iy).ge.eps.or.abs(iz).ge.eps) 
     &        write(*,'(a,3f25.20,4i5)')'consi',ix,iy,iz,i,j,k,ndom
        enddo 
        enddo 
        enddo

        !extap distance paroi
        if(param_int(IFLOW).eq.3) then

            ind_dm_zone(2) = param_int(IJKV  )
            ind_dm_zone(4) = param_int(IJKV+1)
            ind_dm_zone(6) = param_int(IJKV+2)
            call dist_extrap( param_int(NDIMDX), param_int(NDIMDX_XYZ),
     &                        param_int(NIJK), param_int( NIJK_XYZ ),
     &                        ind_dm_zone, degener , dist)
        endif



!$OMP END SINGLE

      endif !extrap typezone

        l_initmtr=.true.
      ENDIF


      END
