c***********************************************************************
c     $Date: 2010-11-04 13:25:50 +0100 (Thu, 04 Nov 2010) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine ssdom_lu_ijk(ndo, ndfin, nssiter, mx_ssdom_lu,  nitcfg,
     &                       iseuil,
     &                       size_ssdom, ijkv, ijkv_lu, it_target_ssdom,
     &                       no_lu, nisdom_residu,
     &                       ind_dm, ijk_lu)
c***********************************************************************
c_U   USER : TERRACOL
c
c     ACT
c_A    determine les indices a parcourir sur un sous domaine pour le calcul de navier-stokes
c      par sous-domaine
c      nsdm: No interne sous-domaine
c
c     I/O
c_/    IN:  ndom, ndsdm
c_/    OUT: ind_dm,ind_coe,ind_grad,....,ind_fluk
c***********************************************************************
      implicit none

      INTEGER_E ndo, nidom, ndfin, nssiter,
     &  mx_ssdom_lu, nitcfg, iseuil,
     & size_ssdom(3), ijkv(3), ijkv_lu(3),ijk_lu(3,mx_ssdom_lu)

      INTEGER_E it_target_ssdom( mx_ssdom_lu), no_lu(mx_ssdom_lu),
     & nisdom_residu( nssiter), ind_dm(6, mx_ssdom_lu , nssiter)


c Var loc
      logical ldone(mx_ssdom_lu),lsearch_global,lsearch
      INTEGER_E i,j,k,l,i0,j0,k0,ishift,count,count1,nd,no_sdm,
     & ipara,k2,inc, imin_lu,jmin_lu,kmin_lu,imax_lu,jmax_lu,kmax_lu,
     & ind_lu(6,mx_ssdom_lu),c1,c2,nd2,nd3,iverbs


      iverbs = 0

        count = 0
        do k=1,ijkv_lu(3)
        do j=1,ijkv_lu(2)
        do i=1,ijkv_lu(1)

           no_sdm= i + (j-1)*ijkv_lu(1) +(k-1)*ijkv_lu(1)*ijkv_lu(2)

          if (it_target_ssdom( no_sdm ).ge.iseuil) then
           ldone(no_sdm)    =.true.
           count            = count + 1
           ijk_lu(1,count)  = i
           ijk_lu(2,count)  = j
           ijk_lu(3,count)  = k
           no_lu(  count  ) = no_sdm
          else
           ldone(no_sdm)=.false.
          endif

        enddo
        enddo
        enddo


        nisdom_residu( nitcfg ) = count

        !CCCCCCCCCCCCCCCCCCC
        !determine taille sous-domaine LU dans chaque direction I,J,K
        !CCCCCCCCCCCCCCCCCCC

        count = 0
        c1    = 0

        DO 10 nd=1,nisdom_residu( nitcfg )

         no_sdm =  no_lu(nd) 

         if(.not.ldone(no_sdm)) goto 10 !si bloc deja traite par merging, on skip

         c1  = c1+1
 
         do ipara = 1,6,1

          if(ipara.eq.1) then
            k0      = ijk_lu(1, nd)
            k2      = ijk_lu(1, nd)
            inc     = 1
            ishift  =-1
          elseif(ipara.eq.2) then
            k0      = ijk_lu(1,nd)
            k2      = ijkv_lu(1)- ijk_lu(1,nd) + 1
            inc     = 1
            ishift  = 1
          elseif(ipara.eq.3) then
            k0      = ijk_lu(2, nd)
            k2      = ijk_lu(2, nd)
            inc     = ijkv_lu(1)
            ishift  =-1
          elseif(ipara.eq.4) then
            k0      = ijk_lu(2, nd)
            k2      = ijkv_lu(2)-ijk_lu(2,nd)  + 1
            inc     = ijkv_lu(1)
            ishift  = 1
          elseif(ipara.eq.5) then
            k0      = ijk_lu(3, nd)
            k2      = ijk_lu(3, nd)
            inc     = ijkv_lu(1)*ijkv_lu(2)
            ishift  =-1
          else
            k0      = ijk_lu(3, nd)
            k2      = ijkv_lu(3)- ijk_lu(3, nd) + 1
            inc     = ijkv_lu(1)*ijkv_lu(2)
            ishift  = 1
          endif

          ind_lu(ipara,c1) = k0
          j             = ind_lu(ipara,c1) - ishift
          count1        = 0
          l             = no_sdm 

          do while( (ind_lu(ipara,c1).eq.k0).and.(count1.lt.k2)
     &                                   .and.ldone(l)) 
            j        = j + ishift
            l        = no_sdm +  count1*inc*ishift
            count1   = count1 + 1

            if(.not.ldone(l))             ind_lu(ipara,c1) = j-ishift
            if(count1.ge.k2.and.ldone(l)) ind_lu(ipara,c1) = j
          enddo

         enddo !ipara

         !test existance de sous-domaine qui se touche
         ! Si oui,  alors elargit le sous-domaine
         ! car problem potentiel calcul coe et Cond limite en multithreade.
         lsearch_global   = .true.
         count1        = 0
         do while(lsearch_global)

         count1         = count1 + 1
         lsearch_global = .false.

         do ipara=1,2

            if(ipara.eq.1) then
              k0      = 1
              ishift  =-1
            else
              k0      = ijkv_lu(1)
              ishift  = 1
            endif

            lsearch   = .true.
            do while( ind_lu(ipara,c1).ne.k0.and.lsearch)

             lsearch   = .false.
             i         = ind_lu(ipara,c1) +ishift

             do k=max(1, ind_lu(5,c1)-1), min(ind_lu(6,c1)+1,ijkv_lu(3))
             do j=max(1, ind_lu(3,c1)-1), min(ind_lu(4,c1)+1,ijkv_lu(2))

               l    = i + (j-1)*ijkv_lu(1) +(k-1)*ijkv_lu(1)*ijkv_lu(2)
     
               if(ldone(l)) then
                 !ldone(l)       = .false.
                 lsearch_global = .true.
                 lsearch        = .true.
               endif
             enddo
             enddo
             if(lsearch) ind_lu(ipara,c1) = ind_lu(ipara,c1) +ishift
           enddo !while
         enddo !ipara

         do ipara=3,4

            if(ipara.eq.3) then
              k0      = 1
              ishift  =-1
            else
              k0      = ijkv_lu(2)
              ishift  = 1
            endif

            lsearch   = .true.
            do while( ind_lu(ipara,c1).ne.k0.and.lsearch)

             lsearch   = .false.
                j      = ind_lu(ipara,c1) +ishift
             do k=max(1, ind_lu(5,c1)-1), min(ind_lu(6,c1)+1,ijkv_lu(3))
             do i=max(1, ind_lu(1,c1)-1), min(ind_lu(2,c1)+1,ijkv_lu(1))

               l    = i + (j-1)*ijkv_lu(1) +(k-1)*ijkv_lu(1)*ijkv_lu(2)
     
               if(ldone(l)) then
                 !ldone(l)       = .false.
                 lsearch_global = .true.
                 lsearch        = .true.
               endif
             enddo
             enddo
             if(lsearch) ind_lu(ipara,c1) = ind_lu(ipara,c1) +ishift
           enddo !while
         enddo !ipara

         do ipara=5,6

            if(ipara.eq.5) then
              k0      = 1
              ishift  =-1
            else
              k0      = ijkv_lu(3)
              ishift  = 1
            endif

            lsearch   = .true.
            do while( ind_lu(ipara,c1).ne.k0.and.lsearch)

             lsearch   = .false.
               k       = ind_lu(ipara,c1) +ishift
             do j=max(1, ind_lu(3,c1)-1), min(ind_lu(4,c1)+1,ijkv_lu(2))
             do i=max(1, ind_lu(1,c1)-1), min(ind_lu(2,c1)+1,ijkv_lu(1))

               l    = i + (j-1)*ijkv_lu(1) +(k-1)*ijkv_lu(1)*ijkv_lu(2)
     
               if(ldone(l)) then
                 !ldone(l)       = .false.
                 lsearch_global = .true.
                 lsearch        = .true.
               endif
             enddo
             enddo
             if(lsearch) ind_lu(ipara,c1) = ind_lu(ipara,c1) +ishift
            enddo !while
         enddo !ipara
        !if(count1.ne.1)write(*,'(a,4i6)')'passijk',ndom,nitcfg,nd,count1
        enddo !while search global

        do k=ind_lu(5,c1),ind_lu(6,c1)
        do j=ind_lu(3,c1),ind_lu(4,c1)
        do i=ind_lu(1,c1),ind_lu(2,c1)

            count    = count + 1
            l        = i + (j-1)*ijkv_lu(1) +(k-1)*ijkv_lu(1)*ijkv_lu(2)
            !no_lu(count,ndo) = l
            ldone(l)         = .false.
        enddo
        enddo
        enddo

 10   continue

      !Protection recouvrement
      !On recherche eventuels recouvrement entre sous-domiane. Si oui: MERGE
      !
      lsearch =.true.
      do while(lsearch.and.(c1.ne.1) )
      lsearch=.false.

      nd   = 1
      do while( (nd.le.c1).and.(c1.ne.1) )

        nd2     = 1
        do while( nd2.le.c1)

         IF( (nd2.ne.nd.and.(c1.ne.1))
     &      .and.(ind_lu(1,nd)-1.le.ind_lu(2,nd2) )
     &      .and.(ind_lu(2,nd)+1.ge.ind_lu(1,nd2) )
     &      .and.(ind_lu(3,nd)-1.le.ind_lu(4,nd2) )
     &      .and.(ind_lu(4,nd)+1.ge.ind_lu(3,nd2) )
     &      .and.(ind_lu(5,nd)-1.le.ind_lu(6,nd2) )
     &      .and.(ind_lu(6,nd)+1.ge.ind_lu(5,nd2) ))THEN


            lsearch=.true.

            ind_lu(1,nd)=min(ind_lu(1,nd2),ind_lu(1,nd))
            ind_lu(2,nd)=max(ind_lu(2,nd2),ind_lu(2,nd))
            ind_lu(3,nd)=min(ind_lu(3,nd2),ind_lu(3,nd))
            ind_lu(4,nd)=max(ind_lu(4,nd2),ind_lu(4,nd))
            ind_lu(5,nd)=min(ind_lu(5,nd2),ind_lu(5,nd))
            ind_lu(6,nd)=max(ind_lu(6,nd2),ind_lu(6,nd))
            
            do nd3 = nd2,c1-1
             ind_lu(:,nd3)= ind_lu(:,nd3+1)
            enddo

            c1      = c1 - 1

         ELSE
           nd2 = nd2+1
         ENDIF
        enddo !while
        nd = nd+1
      enddo !while
      enddo !while

      !Mise a jour finale tableau sous-domain 
      count = 0
      do nd=1,c1

        do k=ind_lu(5,nd),ind_lu(6,nd)
        do j=ind_lu(3,nd),ind_lu(4,nd)
        do i=ind_lu(1,nd),ind_lu(2,nd)

            count    = count + 1
            l        = i + (j-1)*ijkv_lu(1) +(k-1)*ijkv_lu(1)*ijkv_lu(2)

            ldone(l)     =.false.
            no_lu(count) = l
        enddo
        enddo
        enddo

        ndfin = ndfin + 1

        ind_dm(1,ndfin,nitcfg) = (ind_lu(1,nd)-1)*size_ssdom(1) + 1
        ind_dm(2,ndfin,nitcfg) =     ind_lu(2,nd)*size_ssdom(1)
        if(ind_lu(2,nd).eq.ijkv_lu(1)) ind_dm(2,ndfin,nitcfg) =ijkv(1)
        ind_dm(3,ndfin,nitcfg) = (ind_lu(3,nd)-1)*size_ssdom(2) + 1
        ind_dm(4,ndfin,nitcfg) =     ind_lu(4,nd)*size_ssdom(2)
        if(ind_lu(4,nd).eq.ijkv_lu(2)) ind_dm(4,ndfin,nitcfg) =ijkv(2)
        ind_dm(5,ndfin,nitcfg) = (ind_lu(5,nd)-1)*size_ssdom(3) + 1
        ind_dm(6,ndfin,nitcfg) =     ind_lu(6,nd)*size_ssdom(3)
        if(ind_lu(6,nd).eq.ijkv_lu(3)) ind_dm(6,ndfin,nitcfg) =ijkv(3)

      enddo

      if(count.lt.  nisdom_residu(nitcfg) ) then
!$OMP SINGLE
       write(*,'(a,4i6)')'WARNING MERGE ssdom_Lu',nd,nitcfg,
     &            nisdom_residu(nitcfg),count
!$OMP END SINGLE
      endif

       do k=1,ijkv_lu(3)
       do j=1,ijkv_lu(2)
       do i=1,ijkv_lu(1)

           l= i + (j-1)*ijkv_lu(1) +(k-1)*ijkv_lu(1)*ijkv_lu(2)

          if (it_target_ssdom(no_sdm).ge.iseuil.and.ldone(l)) then
!$OMP SINGLE
            write(*,'(a,5i6)')'WARNING MERGE ssdom_Lu',nd,nitcfg,i,j,k
!$OMP END SINGLE
          endif

       enddo
       enddo
       enddo


      nisdom_residu(nitcfg) = count

      end
