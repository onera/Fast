c***********************************************************************
c     $Date: 2010-07-16 17:20:59 +0200 (Fri, 16 Jul 2010) $
c     $Revision: 59 $
c     $Author: MarcTerracol $
c***********************************************************************
      subroutine init_ssiter_bloc(nd , nitcfg, nssiter, nitrun,
     &     lssiter_loc, itypcp, flag_res,
     &     ijkv, ijkv_lu, ijk_lu, size_ssdom,
     &     mx_ssdom_lu, iskip_lu, dtloc,
     &     ind_dm, nidom_loc, it_bloc, nisdom_residu,
     &     it_lu_ssdom,  it_target_ssdom, it_target_old, no_lu,
     &     param_int)
c***********************************************************************
c     _P                          O N E R A
c     
c     ACT
c     _A    nitcfg     = No sous-iteration courante
c     _A    nssiter    = Nombre de sous-iteration Maximale
c     _A    ind_dm     = Borne du sous-domaine calculee a l iteration courante si Newton adaptatif
c     _A    ssiter_loc = list des domaines reellement calculee a la ssiteration nitcfg
c***********************************************************************
      implicit none

#include "FastC/param_solver.h"

      INTEGER_E nd , nitcfg, nssiter, nitrun, itypcp,
     &     lssiter_loc, mx_ssdom_lu,flag_res,
     &     ijkv(3), ijkv_lu(3), ijk_lu(3,mx_ssdom_lu),
     &     size_ssdom(3), iskip_lu(nssiter,2), dtloc(0:*) 

      INTEGER_E ind_dm(6,mx_ssdom_lu,nssiter), 
     &     nidom_loc(nssiter), nisdom_residu(nssiter),
     &     it_lu_ssdom(mx_ssdom_lu),it_target_ssdom(mx_ssdom_lu),
     &     it_target_old(mx_ssdom_lu), no_lu(mx_ssdom_lu), it_bloc,
     &     param_int(0:*), it_cycl_lbm,level_tg


c     Var loc
      INTEGER_E i,j,k,ndo, icp,nitmin,iseuil,ndfin,b,cycl,count,no_sdm

      !print*,'mx_ssdom_lu= ' , mx_ssdom_lu
      !write(*,*)'lssiter_loc',lssiter_loc,param_int(IFLOW),nd
      
!     initialisation compteur ssiter par domain
      if (nitcfg.eq.1.and.flag_res.eq.1) then

         it_bloc                     = 0
         it_lu_ssdom(1:mx_ssdom_lu)  = 0  

         if(lssiter_loc.eq.1.and.nd.eq.0) then
            call skip_lu(nssiter  , nssiter, iskip_lu )
            call skip_lu(nssiter-1, nssiter, iskip_lu )
         endif
      endif
      !!  si on declenche le calcul residu, on ne modifie pas les souszones 
      if (flag_res.eq.1) return

      !on evalue les domaine a zapper a cette iteration du newton
      ndfin   = 0

      IF(lssiter_loc.eq.1 ) then

         ndo      = nd

         !on estime le nbr de sous-iteration necessaire en fonction du CFL
         nitmin    = 3          !nombre minimal de ss-iteration
         if(nitcfg.ne.nssiter) then
            iseuil      = iskip_lu( nitcfg, 2) + nitmin
         else
            iseuil      = iskip_lu( nitcfg, 1) + nitmin
         endif
         
         !write(*,*)'iseuil',iseuil, nitcfg

         !on cherche les sous-domain ou it_lu_ssdom(,ndo) <= iseuil
         call ssdom_lu_ijk(nd, ndfin, nssiter, mx_ssdom_lu, nitcfg,
     &        iseuil,
     &        size_ssdom, ijkv, ijkv_lu, it_target_ssdom,
     &        no_lu, nisdom_residu,
     &        ind_dm, ijk_lu )

      !ELSE IF (param_int(EXPLOC).eq.1) THEN 
      ELSE IF (lssiter_loc.eq.2.and.param_int(EXPLOC).eq.1) THEN ! Explicite local instationnaire d'ordre 3 

         
         cycl    = param_int(NSSITER)/param_int(LEVEL)      
      
         IF (MOD(nitcfg,cycl)==1.OR.MOD(nitcfg,cycl)==cycl/2.OR.
     & MOD(nitcfg,cycl)==cycl-1) THEN

         ndfin = ndfin + 1

         !write(*,*)'ndfin, nitcfg', ndfin, nitcfg

         ind_dm(1, ndfin, nitcfg)   = 1
         ind_dm(3, ndfin, nitcfg)   = 1
         ind_dm(5, ndfin, nitcfg)   = 1
         ind_dm(2, ndfin, nitcfg)   = ijkv(1)
         ind_dm(4, ndfin, nitcfg)   = ijkv(2)
         ind_dm(6, ndfin, nitcfg)   = ijkv(3)
         
         nisdom_residu( nitcfg)     = 1

         END IF               

      ELSE                      !explict global, LBM  ou implicit ss-iter constant 
     
         ndfin = ndfin + 1

         ind_dm(1, ndfin, nitcfg)   = 1
         ind_dm(3, ndfin, nitcfg)   = 1
         ind_dm(5, ndfin, nitcfg)   = 1
         ind_dm(2, ndfin, nitcfg)   = ijkv(1)
         ind_dm(4, ndfin, nitcfg)   = ijkv(2)
         ind_dm(6, ndfin, nitcfg)   = ijkv(3)
 
         !!!  necessaire pour calcul residu OK
         count = 0
         do k=1,ijkv_lu(3)
         do j=1,ijkv_lu(2)
         do i=1,ijkv_lu(1)
           no_sdm= i + (j-1)*ijkv_lu(1) +(k-1)*ijkv_lu(1)*ijkv_lu(2)

           count            = count + 1
           no_lu(  count  ) = no_sdm
         enddo
         enddo
         enddo

         nisdom_residu( nitcfg ) = count
         !nisdom_residu( nitcfg)     = 1

C     !$       write(*,'(a16,9i4)')'loop sousdom RK3 ',
C     !$   &    nd,ind_dm(1,ndfin,nitcfg),ind_dm(2,ndfin,nitcfg),
C     !$   &       ind_dm(3,ndfin,nitcfg),ind_dm(4,ndfin,nitcfg),
C     !$   &       ind_dm(5,ndfin,nitcfg),ind_dm(6,ndfin,nitcfg),
C     !$   &       ndfin,nitcfg

      ENDIF

      nidom_loc(nitcfg) = ndfin

!     print*, 'it : ',nitcfg,ind_dm(1, ndfin, nitcfg)

!     write(*,*)'nidom_loc(nitcfg)',ndfin,nitcfg,nd

      end


