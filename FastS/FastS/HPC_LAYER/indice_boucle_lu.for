c***********************************************************************
c     $Date: 2010-11-04 13:25:50 +0100 (Thu, 04 Nov 2010) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine indice_boucle_lu(ndom,ndsdm,nsdom_lu,icp,
     &                            ind_dm,
     &                            thread_topology, ind_sdm)
c***********************************************************************
c_U   USER : TERRACOL
c
c     ACT
c_A    determine les indices a parcourir sur un sous domaine pour le calcul de navir-stokes
c      par sous-domaine
c      ndom: No domaine
c      nsdm: No interne sous-domaine
c
c     I/O
c_/    IN:  ndom, ndsdm
c_/    OUT: ind_sdm,ind_coe,ind_grad,....,ind_fluk
c***********************************************************************
      implicit none

      include "omp_lib.h"

      INTEGER_E ndom,ndsdm,nsdom_lu,icp,ind_dm(6),ind_sdm(6),
     & thread_topology(3)

c Var loc
      logical l1,l2,l21,l12
      INTEGER_E p1,p2,r1,i,j,k,ivsdm,jvsdm,kvsdm,
!     & lmin(3),Imin,Kmin,Jmin,ncombi,first(4),loptim,l,
     & lmin,Imin,Kmin,Jmin,ncombi,first(7),loptim,l,
     & Imax,Kmax,Jmax,ivloc,jvloc,kvloc,iseuil, iverbs
      REAL_E cost, cout, surcout

      !on assure une decoupe minimale
      Kmax = 1
      Jmax = 1
      Imax = 1

      If(nsdom_lu.eq.1) Then  !pas de sous domaine; un seul bloc!

      ind_sdm(1)  = ind_dm(1)
      ind_sdm(2)  = ind_dm(2)
      ind_sdm(3)  = ind_dm(3)
      ind_sdm(4)  = ind_dm(4)
      ind_sdm(5)  = ind_dm(5)
      ind_sdm(6)  = ind_dm(6)

      Else   ! calcul indice sous-domaine
      
 
        first(1)=1
        ncombi  =1
        if(mod(nsdom_lu,2).eq.0) then
          ncombi        = ncombi + 1
          first(ncombi) = 2
        endif
        if(mod(nsdom_lu,3).eq.0) then
          ncombi        = ncombi + 1
          first(ncombi) = 3
        endif
        if(mod(nsdom_lu,4).eq.0) then
         ncombi        = ncombi + 1
          first(ncombi) = 4
        endif
        if(mod(nsdom_lu,5).eq.0) then
          ncombi        = ncombi + 1
          first(ncombi) = 5
        endif
        if(mod(nsdom_lu,6).eq.0) then
          ncombi        = ncombi + 1
          first(ncombi) = 6
        endif
        if(mod(nsdom_lu,7).eq.0) then
          ncombi        = ncombi + 1
          first(ncombi) = 7
        endif
        if(mod(nsdom_lu,8).eq.0) then
          ncombi        = ncombi + 1
          first(ncombi) = 8
        endif

           

        ! on determine l'epaisseur minimal du sous domaine
        if(icp.eq.2) then   !integration explicite: on minimise l'epaisseur pour optimiser la repartition
          !lmin(1) = 1000    !(ific+2)
          lmin =  4         !(ific+2)
          !lmin(1) =   4     !(ific+2)
          !lmin(2) =   4     !(ific+2)
          !lmin(3) =   4     !(ific+2)
        else                !integration implicite: on evite les domaines trop fin pour preserver la convergence du LU
          !lmin(:) = 10
          lmin = 10
        endif
 

        ivloc = ind_dm(2)- ind_dm(1) + 1
        jvloc = ind_dm(4)- ind_dm(3) + 1
        kvloc = ind_dm(6)- ind_dm(5) + 1

        ! determination des tailles a decouper suivant type domaine 2d/3d
        Imin = max(1,ivloc/lmin)
        Jmin = max(1,jvloc/lmin)
        Kmin = max(1,kvloc/lmin)

        !Imin = ivloc/lmin(1)
        !if(lmin(1).eq.1000) Imin = 1
        !Jmin = jvloc/lmin(2)
        !Kmin = kvloc/lmin(3)
        if(kvloc.eq.1) Kmin = 1
        iseuil= Kmin*Jmin

        !if(ndsdm.eq.12) write(*,*)'size',ivloc,jvloc,kvloc
        !if(ndsdm.eq.1) write(*,*)'min',Imin,Jmin,Kmin

        IF( Kmin*Jmin.ge.nsdom_lu) THEN !dom 3d: decoupe direction k (et j eventuellement)

           Imax   = 1 
           cout   = 100000.
           surcout= 2.           !!on privilegie la decoupe en K versus J en 3D
           loptim = 0
           do l =1,ncombi
             p1 = first(l)
             p2 = nsdom_lu/p1
        !if(ndsdm.eq.1) write(*,*)'p0',p1,p2,l
             if(p2.le.Kmin.and.p1.le.Jmin) then
                cost = mod(kvloc,p2)/float(kvloc/p2)
     &               + mod(jvloc,p1)/float(jvloc/p1)*surcout

                if(cost.lt.cout) then
                       loptim =-l
                       cout   =  cost
        !if(ndsdm.eq.1) write(*,*)'p2',cost,loptim,l
                endif

             endif
             if(p1.le.Kmin.and.p2.le.Jmin) then
                cost = mod(kvloc,p1)/float(kvloc/p1)
     &               + mod(jvloc,p2)/float(jvloc/p2)*surcout

                if(cost.lt.cout) then
                       loptim = l
                       cout   =  cost
        !if(ndsdm.eq.1) write(*,*)'p1',cost,loptim,l
                endif
             endif
         
           enddo

           if(loptim.eq.0) then
                !pas assez de travail a partager, le thread 1 prend tout
                 Jmax =1 
                 Kmax =1 
                 if(ndsdm.eq.1) then
                     ind_sdm(1) = ind_dm(1)
                     ind_sdm(2) = ind_dm(2)
                     ind_sdm(3) = ind_dm(3)
                     ind_sdm(4) = ind_dm(4)
                     ind_sdm(5) = ind_dm(5)
                     ind_sdm(6) = ind_dm(6)
                 else
                     ind_sdm(1) = ind_dm(1)
                     ind_sdm(2) = ind_dm(1)-1
                     ind_sdm(3) = ind_dm(3)
                     ind_sdm(4) = ind_dm(3)-1
                     ind_sdm(5) = ind_dm(5)
                     ind_sdm(6) = ind_dm(5)-1
                 endif
                 goto 1000
            else
                 if(loptim.gt.0) then
                     Kmax = first(loptim)
                     Jmax = nsdom_lu/first(loptim)
                 else
                     Jmax = first(-loptim)
                     Kmax = nsdom_lu/first(-loptim)
                 endif
            endif


        ELSEIF(Imin*Jmin.ge.nsdom_lu) then !dom 2d: decoupe direction k (et j eventuellement)


           Kmax   = 1 
           cout   = 100000.
           surcout= 5.           !!on privilegie la decoupe en J versus I en 2D
           loptim = 0
           do l =1,ncombi
             p1 = first(l)
             p2 = nsdom_lu/p1
             !if(ndsdm.eq.12) write(*,*)'p1,p2',p1,p2,l
             if(p1.le.Jmin.and.p2.le.Imin) then
                cost = mod(jvloc,p1)/float(jvloc/p1)
     &               + mod(ivloc,p2)/float(ivloc/p2)*surcout

             !   if(ndsdm.eq.12) write(*,*)'cout+',cost,l
                if(cost.lt.cout) then
                       loptim = l
                       cout   =  cost
                     !if(ndsdm.eq.12) write(*,*)'opt+',cout,l,loptim
                endif
             endif
             if(p2.le.Jmin.and.p1.le.Imin) then
                cost = mod(jvloc,p2)/float(jvloc/p2)
     &               + mod(ivloc,p1)/float(ivloc/p1)*surcout

                !if(ndsdm.eq.12) write(*,*)'cout-',cost,l
                if(cost.lt.cout) then
                       loptim =-l
                       cout   =  cost
                       !if(ndsdm.eq.12) write(*,*)'opt-',cout,l,loptim
                endif
             endif
           enddo
           if(loptim.eq.0) then
                !pas assez de travail a partager, le thread 1 prend tout
                 Imax =1 
                 Jmax =1 
                 if(ndsdm.eq.1) then
                     ind_sdm(1) = ind_dm(1)
                     ind_sdm(2) = ind_dm(2)
                     ind_sdm(3) = ind_dm(3)
                     ind_sdm(4) = ind_dm(4)
                     ind_sdm(5) = ind_dm(5)
                     ind_sdm(6) = ind_dm(6)
                 else
                     ind_sdm(1) = ind_dm(1)
                     ind_sdm(2) = ind_dm(1)-1
                     ind_sdm(3) = ind_dm(3)
                     ind_sdm(4) = ind_dm(3)-1
                     ind_sdm(5) = ind_dm(5)
                     ind_sdm(6) = ind_dm(5)-1
                 endif
                 goto 1000
            else
                 if(loptim.gt.0) then
                     Jmax = first(loptim)
                     Imax = nsdom_lu/first(loptim)
                 else
                     Imax = first(-loptim)
                     Jmax = nsdom_lu/first(-loptim)
                 endif
            endif
            !if(ndsdm.eq.12) write(*,*)'IJKMax',Imax,Jmax,Kmax

        ELSE  !on blinde les 1er threads prennent tout, les autres dorment....

            p1 = Imin*Jmin
            p2 = Jmin*Kmin
            
            if(p1.gt.p2) Then
              Imax =Imin
              Jmax =Jmin
              Kmax = 1
            else
              Imax = 1  
              Jmax =Jmin
              Kmax =Kmin
            endif

            !if(Imax.gt.Kmax) then
            ! Kmax = 1
            !else
            ! Imax = 1
            !endif

            if(ndsdm.gt.Imax*Jmax*Kmax) then  !thread dormant
             ind_sdm(1) = ind_dm(1)
             ind_sdm(2) = ind_dm(1)-1
             ind_sdm(3) = ind_dm(3)
             ind_sdm(4) = ind_dm(3)-1
             ind_sdm(5) = ind_dm(5)
             ind_sdm(6) = ind_dm(5)-1
            goto 1000
            endif
            !goto 1000
         ENDIF


         !if(iseuil.ge.nsdom_lu) then !dom 3d: decoupe direction k (et j eventuellement)
         if(Imax.eq.1) then !dom 3d: decoupe direction k (et j eventuellement)

           ind_sdm(1) = ind_dm(1)
           ind_sdm(2) = ind_dm(2)
 
           k = 1 +(ndsdm-1)/Jmax
           j = ndsdm -(k-1)*Jmax

           if(j.le.mod(jvloc,Jmax)) then !test pour repartir le residu sur les 1er sous-doamine
             jvsdm       = jvloc/Jmax  + 1
             ind_sdm(3)  = ind_dm(3)     + jvsdm*(j-1) 
             ind_sdm(4)  = ind_dm(3) -1  + jvsdm*j    
           else
             jvsdm       = jvloc/Jmax
             ind_sdm(3)  = ind_dm(3)     + jvsdm*(j-1) + mod(jvloc,Jmax)
             ind_sdm(4)  = ind_dm(3) -1  + jvsdm*j     + mod(jvloc,Jmax)
           endif
           if(k.le.mod(kvloc,Kmax)) then
             kvsdm       = kvloc/Kmax  + 1
             ind_sdm(5)  = ind_dm(5)     + kvsdm*(k-1) 
             ind_sdm(6)  = ind_dm(5) -1  + kvsdm*k    
           else
             kvsdm       = kvloc/Kmax
             ind_sdm(5)  = ind_dm(5)     + kvsdm*(k-1) + mod(kvloc,Kmax)
             ind_sdm(6)  = ind_dm(5) -1  + kvsdm*k     + mod(kvloc,Kmax)
           endif

         else

           !Imax = Jmax
           !Jmax = Kmax
           !Kmax = 1

           ind_sdm(5) = ind_dm(5)
           ind_sdm(6) = ind_dm(6)
 
           j = 1 +(ndsdm-1)/Imax
           i = ndsdm -(j-1)*Imax

           if(i.le.mod(ivloc,Imax)) then
             ivsdm       = ivloc/Imax  + 1
             ind_sdm(1)  = ind_dm(1)     + ivsdm*(i-1) 
             ind_sdm(2)  = ind_dm(1) -1  + ivsdm*i    
           else
             ivsdm       = ivloc/Imax
             ind_sdm(1)  = ind_dm(1)     + ivsdm*(i-1) + mod(ivloc,Imax)
             ind_sdm(2)  = ind_dm(1) -1  + ivsdm*i     + mod(ivloc,Imax)
           endif
           if(j.le.mod(jvloc,Jmax)) then
             jvsdm       = jvloc/Jmax  + 1
             ind_sdm(3)  = ind_dm(3)     + jvsdm*(j-1) 
             ind_sdm(4)  = ind_dm(3) -1  + jvsdm*j    
           else
             jvsdm       = jvloc/Jmax
             ind_sdm(3)  = ind_dm(3)     + jvsdm*(j-1) + mod(jvloc,Jmax)
             ind_sdm(4)  = ind_dm(3) -1  + jvsdm*j     + mod(jvloc,Jmax)
           endif

         endif
      ENDIF  ! calcul indice sous-domaine

1000  continue 

      thread_topology(1) = Imax
      thread_topology(2) = Jmax
      thread_topology(3) = Kmax
     
      iverbs = 0

!$    if(iverbs.ge.1.and.ndsdm.eq.6) write(*,'(a10,3i4)')'topologie ',
!$   &  thread_topology

!$    if(iverbs.ge.1.and.ndsdm.eq.1) write(*,'(a9,6i4)')'loop lu  ',
!$   & ind_sdm(1),ind_sdm(2),
!$   & ind_sdm(3),ind_sdm(4),ind_sdm(5),ind_sdm(6)
      end
