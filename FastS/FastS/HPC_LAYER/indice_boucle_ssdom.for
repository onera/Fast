c***********************************************************************
c     $Date: 2010-11-04 13:25:50 +0100 (Thu, 04 Nov 2010) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine indice_boucle_ssdom(ndo, extended_range,
     &                              ith,jth,kth, ic,jc,kc,
     &                               topo_th, ithread, thread_pos,
     &                               size_cache, 
     &                               synchro_receive_sock, 
     &                               synchro_receive_th,
     &                               synchro_send_sock, 
     &                               synchro_send_th,
     &                               nijk, ijkv,
     &                               ind_dm_zone, ind_dm_sock,
     &                               ind_dm_thread,ijkv_sdm,
     &                               ind_sdm, ind_coe,
     &                               ind_grad, ind_rhs, ind_mjr)
c***********************************************************************
c_U   USER : TERRACOL
c

c_A    determine les indices a parcourir sur un sous domaine pour le calcul de navir-stokes
c      par sous-domaine
c      ndom: No domaine
c      ndo : No interne domaine
c      nsdm: No interne sous-domaine
c
c     I/O
c_/    IN:  ndom, ndsdm, ndo
c_/    OUT: ind_sdm,ind_coe,ind_grad,....,ind_fluk
c***********************************************************************
      implicit none

      INTEGER_E ndo, extended_range,ith,jth,kth, ic,jc,kc, ithread,
     & topo_th(3), size_cache(3), thread_pos(3),
     & synchro_receive_sock(3), synchro_receive_th(3), ijkv(3),
     & synchro_send_sock(3), synchro_send_th(3),
     & nijk(5), ind_dm_zone(6), ind_dm_sock(6),ind_dm_thread(6),
     & ijkv_sdm(3),ind_sdm(6),
     & ind_coe(6),ind_grad(6), ind_rhs(6), ind_mjr(6)

c Var loc
      INTEGER_E inck, shift_coe,shift_grad, shift_rhs,c1,c2,i,j,k,
     & ideb,ifin,l,inc,ific,l20,
     & shift_coe2, shift_gra2,shift_rhs2,shift_mjr2, 
     & shift_coe1, shift_gra1,shift_rhs1,shift_mjr1

      inck = 1
      if(nijk(5).eq.0) inck = 0

      ind_sdm(1)= ind_dm_thread(1)+ (ic-1)*size_cache(1)
      ind_sdm(2)= ind_dm_thread(1)+     ic*size_cache(1)-1
      ind_sdm(3)= ind_dm_thread(3)+ (jc-1)*size_cache(2)
      ind_sdm(4)= ind_dm_thread(3)+     jc*size_cache(2)-1
      ind_sdm(5)= ind_dm_thread(5)+ (kc-1)*size_cache(3)
      ind_sdm(6)= ind_dm_thread(5)+     kc*size_cache(3)-1

      if(ic.eq.ijkv_sdm(1)) ind_sdm(2) = ind_dm_thread(2)
      if(jc.eq.ijkv_sdm(2)) ind_sdm(4) = ind_dm_thread(4)
      if(kc.eq.ijkv_sdm(3)) ind_sdm(6) = ind_dm_thread(6)

  
      !correction eventuelle au bord pour calcul dans 1 rangee maille fictive
      if(extended_range.eq.1) then

       if( ind_sdm(1).eq.1) ind_sdm(1)= 1-nijk(4)+1
       if( ind_sdm(3).eq.1) ind_sdm(3)= 1-nijk(4)+1
       if( ind_sdm(5).eq.1) ind_sdm(5)= 1-nijk(5)+inck

      if(ind_sdm(2).eq.ind_dm_zone(2))ind_sdm(2)=ind_sdm(2)+nijk(4)-1
      if(ind_sdm(4).eq.ind_dm_zone(4))ind_sdm(4)=ind_sdm(4)+nijk(4)-1
      if(ind_sdm(6).eq.ind_dm_zone(6))ind_sdm(6)=ind_sdm(6)+nijk(5)-inck
      endif

         l20=ic*jc*kc
         do i=1,3  !! boucle dir I,J,K

           ideb = 2*i-1
           ifin = 2*i

           inc = 1
           if(i.eq.1) then
            ific = nijk(4)
            l    = ic
           elseif(i.eq.2) then
            ific = nijk(4)
            l    = jc
           else
            ific = nijk(5)
            l    = kc
            inc  = inck
           endif

           If(topo_th(i).eq.1) Then


             if(ind_dm_thread(ideb).eq.ind_dm_sock(ideb).and.l.eq.1)then

                 c1 =  -synchro_send_sock(i)*inc

                 shift_coe1  =-ific
                 shift_gra1  =-inc
                 shift_rhs1  = c1
                 shift_mjr1  = 0

             else
                 shift_coe1  = 0
                 shift_gra1  = 0
                 shift_rhs1  = 0
                 shift_mjr1  =-inc
             endif

             if(ind_dm_thread(ifin).eq.ind_dm_sock(ifin).and.
     &                        l.eq.ijkv_sdm(i)         ) then
          

                 c1 = -synchro_receive_sock(i)*inc
                 c2 =  inc
                 if(synchro_receive_sock(i).eq.1) c2 = -inc

                 shift_coe2  = ific*c2
                 shift_gra2  = c2
                 shift_rhs2  = c1
                 shift_mjr2  = 0
             else
                 shift_coe2  = 0
                 shift_gra2  = 0
                 shift_rhs2  = 0
                 shift_mjr2  =-inc
             endif


           Else ! depend du thread, de jcache,(en plus de debut et fin de bloc)
          
            if(thread_pos(i).eq.1) then

                 !if(ithread.eq.1.and.i.eq.2) write(*,*)'coucou1'

             if(ind_dm_thread(ideb).eq.ind_dm_sock(ideb).and.l.eq.1)then

                 c1 =  -synchro_send_sock(i)*inc

      !if(ithread.eq.1.and.i.eq.2) write(*,*)'coucou11',synchro_send_sock

                 shift_coe1  =-ific
                 shift_coe2  = 0
                 shift_gra1  =-inc
                 shift_gra2  = 0
                 shift_rhs1  = c1
                 shift_rhs2  = 0
                 shift_mjr1  = 0
                 shift_mjr2  =-inc
             elseif( l .eq.ijkv_sdm(i) ) then
                 !if(ithread.eq.1.and.i.eq.2) write(*,*)'coucou12'
                 shift_coe1  = 0
                 shift_coe2  = 0
                 shift_gra1  = 0
                 shift_gra2  =-inc
                 shift_rhs1  = 0
                 shift_rhs2  =-inc
                 shift_mjr1  =-inc
                 shift_mjr2  = 0
             else
                 !if(ithread.eq.1.and.i.eq.2) write(*,*)'coucou13'
                 shift_coe1  = 0
                 shift_coe2  = 0
                 shift_gra1  = 0
                 shift_gra2  = 0
                 shift_rhs1  = 0
                 shift_rhs2  = 0
                 shift_mjr1  =-inc
                 shift_mjr2  =-inc
             endif

            elseif(thread_pos(i).eq.topo_th(i)) then 


              if(l.eq.1 ) then

                 shift_coe1  = 0
                 shift_coe2  = 0
                 shift_gra1  =-inc
                 shift_gra2  = 0
                 shift_rhs1  =-inc
                 shift_rhs2  = 0
                 shift_mjr1  = 0
                 shift_mjr2  =-inc

              elseif(ind_dm_thread(ifin).eq.ind_dm_sock(ifin).and.
     &                             l.eq.ijkv_sdm(i)         ) then
          
                 c1 =  -synchro_receive_sock(i)*inc
                 c2 = -inc
                 if(synchro_receive_sock(i).eq.0) c2 =  inc

                 shift_coe1  = 0
                 shift_coe2  = ific*c2
                 shift_gra1  = 0
                 shift_gra2  = c2
                 shift_rhs1  = 0
                 shift_rhs2  = c1
                 shift_mjr1  =-inc
                 shift_mjr2  = 0
              else
                 shift_coe1  = 0
                 shift_coe2  = 0
                 shift_gra1  = 0
                 shift_gra2  = 0
                 shift_rhs1  = 0
                 shift_rhs2  = 0
                 shift_mjr1  =-inc
                 shift_mjr2  =-inc
              endif

            else

               if(l.eq.1 ) then
                 shift_coe1  = 0
                 shift_coe2  = 0
                 shift_gra1  =-inc
                 shift_gra2  = 0
                 shift_rhs1  =-inc
                 shift_rhs2  = 0
                 shift_mjr1  = 0
                 shift_mjr2  =-inc
              elseif(l.eq.ijkv_sdm(i) ) then
                shift_coe1  = 0
                shift_coe2  = 0
                shift_gra1  = 0
                shift_gra2  =-inc
                shift_rhs1  = 0
                shift_rhs2  =-inc
                shift_mjr1  =-inc
                shift_mjr2  = 0
              else
                shift_coe1  = 0
                shift_coe2  = 0
                shift_gra1  = 0
                shift_gra2  = 0
                shift_rhs1  = 0
                shift_rhs2  = 0
                shift_mjr1  =-inc
                shift_mjr2  =-inc
              endif
            endif

         Endif

         ind_coe(ideb)  = ind_sdm(ideb) + shift_coe1
         ind_grad(ideb) = ind_sdm(ideb) + shift_gra1
         ind_rhs(ideb)  = ind_sdm(ideb) + shift_rhs1
         ind_mjr(ideb)  = ind_sdm(ideb) + shift_mjr1
         ind_coe(ifin)  = ind_sdm(ifin) + shift_coe2
         ind_grad(ifin) = ind_sdm(ifin) + shift_gra2
         ind_rhs(ifin)  = ind_sdm(ifin) + shift_rhs2
         ind_mjr(ifin)  = ind_sdm(ifin) + shift_mjr2

        enddo !! boucle dir I,J,K

      end
