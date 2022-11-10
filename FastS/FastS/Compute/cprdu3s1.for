c***********************************************************************
c     $Date: 2014-03-19 20:08:08 +0100 (mer. 19 mars 2014) $
c     $Revision: 38 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine cprdu3s1(ndom,nitcfg,nssiter,
     &                    neq, ndimdx,ndim_rdm, nitrun,
     &                    iflw, les, ithread,Nbre_thread_actif,omp_mode,
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
     & iflw, les, mx_ssdom_lu, ithread, Nbre_thread_actif, omp_mode,
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
     & ind_loop_lu(6), size_rdm, no_start,no_end 
c     & ii,jj,kk
      REAL_E xinterm(6,2),xp,rmax,rmoy,xro1,xrou1,xrov1,xrow1,xroe1,
     & conv_loo, cut0x

      character *200 cfich,fmt

#include "FastS/formule.h"

c      if (ithread.eq.1.and.ndom.eq.0)
c     &  write(*,'(a,6i4,a,4i4)')'ssdom=',ind_loop,' nstep ',nitcfg,
c     & ndim_rdm,Nbre_thread_actif,nisdom_residu(nitcfg)

      cut0x = 1e-20

      !if(nisdom_residu(nitcfg).eq.0) goto 1000

      conv_loo = ALOG10(epsi)

        if(ndim_rdm.le.Nbre_thread_actif) then
           no_start = ithread
           no_end   = ithread
           if(ithread.gt.ndim_rdm) goto 1000       !on skippe
        else
           size_rdm = ndim_rdm/Nbre_thread_actif
           no_start = 1+ (ithread-1)*size_rdm
           no_end   =    ithread*size_rdm
           if(ithread.eq.Nbre_thread_actif) no_end = ndim_rdm
        endif

        DO no_rdm=no_start,no_end
#include   "FastS/Compute/cprdu3s1_incl.for"
        ENDDO

 1000 continue


      !if(carte_residu.eq.1) then
c      if(nitrun.eq.-7.and.ndom.eq.41.and.ithread.eq.1) then
CCC!$OMP SINGLE

c        do i=1,200
c          cfich(i:i)=' '
c        enddo

c         cfich(1:11)='svt_ijk_dom'

c         if(ndom.lt.10) then
c           write(cfich,'(a,i1.1)')trim(cfich)//'0000',ndom
c         elseif(ndom.ge.10.and.ndom.lt.100) then
c           write(cfich,'(a,i2.2)')trim(cfich)//'000',ndom
c         elseif(ndom.ge.100.and.ndom.lt.1000) then
c           write(cfich,'(a,i3.3)')trim(cfich)//'00',ndom
c         elseif(ndom.ge.1000.and.ndom.lt.10000) then
c           write(cfich,'(a,i4.4)')trim(cfich)//'0',ndom
c         elseif(ndom.ge.10000.and.ndom.lt.100000) then
c           write(cfich,'(a,i5.5)')trim(cfich),ndom
c         else
c           write(*,*)'erreur definition cfich'
c         endif

c         cfich= cfich(1:16)//'_iter='
c         if(nitrun.lt.10) then
c           write(cfich,'(a,i1.1)')trim(cfich)//'00',nitrun
c         elseif(nitrun.ge.10.and.nitrun.lt.100) then
c           write(cfich,'(a,i2.2)')trim(cfich)//'0',nitrun
c         elseif(nitrun.ge.100.and.nitrun.lt.1000) then
c           write(cfich,'(a,i3.3)')trim(cfich),nitrun
c         else
c           write(*,*)'erreur definition cfich'
c         endif

c         cfich= cfich(1:25)//'_ssiter='

c         if(nitcfg.lt.10) then
c           write(cfich,'(a,i1.1)')trim(cfich)//'00',nitcfg
c         elseif(nitcfg.ge.10.and.nitcfg.lt.100) then
c           write(cfich,'(a,i2.2)')trim(cfich)//'0',nitcfg
c         elseif(nitcfg.ge.100.and.nitcfg.lt.1000) then
c           write(cfich,'(a,i3.3)')trim(cfich),nitcfg
c         else
c           write(*,*)'erreur definition cfich'
c         endif

c         cfich= cfich(1:36)//'.dat'


c      iunit = 1
cc      fmt   ='(3i5)'
cc      fmt   ='(3i5, 6(1x,e23.16))       '
c      fmt   ='(3i5, 5(1x,e23.16))       '
c      open(unit=iunit,file=cfich, form='formatted',status='new', 
c     &            iostat=i)
c      if(i.ne.0) open(unit=iunit,file=cfich, form='formatted',
c     &                status='old', iostat=i)
c

c      write(iunit,'(a)')'ZONE T="join.11"'
c      write(iunit,'(a,i4,a,i4,a,i4,a)')'I=',ind_loop(2),
c     &       'J=', ind_loop(4),  'K=',ind_loop(6), 'ZONETYPE=Ordered'
c      write(iunit,'(a)')'DATAPACKING=POINT'

c      do k = ind_loop(5), ind_loop(6)
c      do j = ind_loop(3), ind_loop(4)
c      do i = ind_loop(1), ind_loop(2)

c        l = inddm(i,j,k)
c        write(iunit,fmt) i,j,k,abs(drodm(l,1)),abs(drodm(l,2)),
c     &                         abs(drodm(l,3)),abs(drodm(l,4)),
c     &                         abs(drodm(l,5))
c
c      enddo
c      enddo
c      enddo
c      close(iunit)
c
cCCC!$OMP END SINGLE
c      endif


      END
