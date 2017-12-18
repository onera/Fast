c***********************************************************************
c     $Date: 2014-03-19 20:08:08 +0100 (mer. 19 mars 2014) $
c     $Revision: 38 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine cprdu3s1(ndom,nitcfg,nssiter,
     &                    neq, ndimdx,ndim_rdm, nitrun,
     &                    iflw, les, ithread, omp_mode,
     &                    mx_ssdom_lu, it_bloc, epsi,
     &                    size_ssdom, nisdom_residu,
     &                    nijk, ijkv, ijkv_lu, 
     &                    ind_loop,
     &                    it_lu_ssdom, it_target_ssdom,
     &                    it_target_old, it_temp_ssdom, no_lu, 
     &                    drodm,rdm_sdm)
c***********************************************************************
c_U   USER : PECHIER
c
c     ACT
c_A    Calcul des residus moyens par sous iteration,
c_A    variable et enfin domaine !
c
c     VAL
c_V    Processeur domaine; NS 5 equations
c
c     INP
c_I    nssiter : nbre total de sous-iterations
c
c     OUT
c_O    rdm_sdm : residus newton
c***********************************************************************
      implicit none

      INTEGER_E ndom,nssiter,neq, ndimdx ,ndim_rdm,nitcfg,nitrun,
     & iflw, les, mx_ssdom_lu, ithread, omp_mode,
     & size_ssdom(3), nijk(5), ijkv(3), ijkv_lu(3), ind_loop(6)

      REAL_E drodm(ndimdx,neq)
      REAL_E rdm_sdm(nssiter,neq,2,ndim_rdm)
      REAL_E epsi


      INTEGER_E nisdom_residu(nssiter), it_temp_ssdom(mx_ssdom_lu),
     & it_lu_ssdom(mx_ssdom_lu),it_target_ssdom(mx_ssdom_lu),
     & it_target_old(mx_ssdom_lu), no_lu(mx_ssdom_lu), it_bloc

c Var loc
      logical lcomput
      INTEGER_E l,i,j,k,ne,no_rdm,nd_rdm,i_lu,j_lu,k_lu,it,iunit,
     & ind_loop_lu(6),imin_lu,jmin_lu,kmin_lu,imax_lu,jmax_lu,kmax_lu
c     & ii,jj,kk
      REAL_E xinterm(6,2),xp,rmax,rmoy,xro1,xrou1,xrov1,xrow1,xroe1,
     & conv_loo, cut0x

      character *200 cfich,fmt

#include "FastS/formule.h"

      !write(*,'(a,7i7)')'ssdom=',ind_loop,ndom

      cut0x = 1e-20

      if(nisdom_residu(nitcfg).ne.0) then
        if(omp_mode .eq.0) then
!$OMP SINGLE
          it_bloc= it_bloc+1
!$OMP END SINGLE
        else
          it_bloc= it_bloc+1
        endif
      else
        goto 1000
      endif

      conv_loo = ALOG10(epsi)


      imin_lu = 1 + ind_loop(1)/size_ssdom(1)
      imax_lu =  max( 1, ind_loop(2)/size_ssdom(1))
      jmin_lu = 1 + ind_loop(3)/size_ssdom(2)
      jmax_lu =  max( 1, ind_loop(4)/size_ssdom(2))
      kmin_lu = 1 + ind_loop(5)/size_ssdom(3)
      kmax_lu =  max( 1, ind_loop(6)/size_ssdom(3))

c       write(*,'(a,9i7)')'ssdom=',imin_lu,jmin_lu,kmin_lu,
c     & imax_lu,jmax_lu,kmax_lu,ndo,nitcfg,nisdom_residu(nitcfg)

      if(omp_mode .eq.0) then
!$OMP DO SCHEDULE(DYNAMIC,1)
        DO no_rdm=1,ndim_rdm

#include "FastS/Compute/cprdu3s1_incl.for"

        ENDDO
!$OMP END  DO
      
      else
        DO no_rdm=1,ndim_rdm
#include "FastS/Compute/cprdu3s1_incl.for"
        ENDDO
      endif



 1000 continue


      !if(carte_residu.eq.1) then
      if(nitrun.eq.-10) then
!$OMP SINGLE

        do i=1,200
          cfich(i:i)=' '
        enddo

         cfich(1:11)='svt_ijk_dom'

         if(ndom.lt.10) then
           write(cfich,'(a,i1.1)')trim(cfich)//'0000',ndom
         elseif(ndom.ge.10.and.ndom.lt.100) then
           write(cfich,'(a,i2.2)')trim(cfich)//'000',ndom
         elseif(ndom.ge.100.and.ndom.lt.1000) then
           write(cfich,'(a,i3.3)')trim(cfich)//'00',ndom
         elseif(ndom.ge.1000.and.ndom.lt.10000) then
           write(cfich,'(a,i4.4)')trim(cfich)//'0',ndom
         elseif(ndom.ge.10000.and.ndom.lt.100000) then
           write(cfich,'(a,i5.5)')trim(cfich),ndom
         else
           write(*,*)'erreur definition cfich'
         endif

         cfich= cfich(1:16)//'_iter='
         if(nitrun.lt.10) then
           write(cfich,'(a,i1.1)')trim(cfich)//'00',nitrun
         elseif(nitrun.ge.10.and.nitrun.lt.100) then
           write(cfich,'(a,i2.2)')trim(cfich)//'0',nitrun
         elseif(nitrun.ge.100.and.nitrun.lt.1000) then
           write(cfich,'(a,i3.3)')trim(cfich),nitrun
         else
           write(*,*)'erreur definition cfich'
         endif

         cfich= cfich(1:25)//'_ssiter='

         if(nitcfg.lt.10) then
           write(cfich,'(a,i1.1)')trim(cfich)//'00',nitcfg
         elseif(nitcfg.ge.10.and.nitcfg.lt.100) then
           write(cfich,'(a,i2.2)')trim(cfich)//'0',nitcfg
         elseif(nitcfg.ge.100.and.nitcfg.lt.1000) then
           write(cfich,'(a,i3.3)')trim(cfich),nitcfg
         else
           write(*,*)'erreur definition cfich'
         endif

         cfich= cfich(1:36)//'.dat'


      iunit = 1
c      fmt   ='(3i5)'
      fmt   ='(3i5, 6(1x,e23.16))       '
      open(unit=iunit,file=cfich, form='formatted',status='new', 
     &            iostat=i)
      if(i.ne.0) open(unit=iunit,file=cfich, form='formatted',
     &                status='old', iostat=i)


      write(iunit,'(a)')'ZONE T="join.11"'
      write(iunit,'(a,i4,a,i4,a,i4,a)')'I=',ind_loop(2),
     &       'J=', ind_loop(4),  'K=',ind_loop(6), 'ZONETYPE=Ordered'
      write(iunit,'(a)')'DATAPACKING=POINT'

      do k = ind_loop(5), ind_loop(6)
      do j = ind_loop(3), ind_loop(4)
      do i = ind_loop(1), ind_loop(2)

        l = inddm(i,j,k)
        write(iunit,fmt) i,j,k,abs(drodm(l,1)),abs(drodm(l,2)),
     &                         abs(drodm(l,3)),abs(drodm(l,4)),
     &                         abs(drodm(l,5)),abs(drodm(l,6))
c
      enddo
      enddo
      enddo
      close(iunit)

!$OMP END SINGLE
      endif


      END
