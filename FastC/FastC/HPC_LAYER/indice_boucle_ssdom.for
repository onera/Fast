c***********************************************************************
c     $Date: 2010-11-04 13:25:50 +0100 (Thu, 04 Nov 2010) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine indice_boucle_ssdom(ndo, extended_range,
     &                              ith,jth,kth, ic,jc,kc,kfludom,
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
     &                               ind_grad, ind_rhs, 
     &                               ind_mjr, ind_ssa)
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
c_/    OUT: ind_sdm,ind_coe,ind_grad,....,ind_flt 
c***********************************************************************
      implicit none

      INTEGER_E ndo, extended_range,ith,jth,kth, ic,jc,kc, ithread,
     & topo_th(3), size_cache(3), thread_pos(3),
     & synchro_receive_sock(3), synchro_receive_th(3), ijkv(3),
     & synchro_send_sock(3), synchro_send_th(3),
     & nijk(5), ind_dm_zone(6), ind_dm_sock(6),ind_dm_thread(6),
     & ijkv_sdm(3),ind_sdm(6),ind_ssa(6),
     & ind_coe(6),ind_grad(6), ind_rhs(6), ind_mjr(6),
     & kfludom 

c Var loc
      INTEGER_E inck, shift_coe,shift_grad, shift_rhs,c1,c2,i,j,k,
     & ideb,ifin,l,inc,ific,
     & shift_coe2, shift_gra2,shift_rhs2,shift_mjr2, shift_flt2,
     & shift_coe1, shift_gra1,shift_rhs1,shift_mjr1, shift_flt1,
     & shift_ssa1, shift_ssa2,
     & delta

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

  
c      if(ithread.eq.2) then
c       write(*,*)'sockS', synchro_send_sock
c       write(*,*)'sockR', synchro_receive_sock
c       write(*,*)'threS', synchro_send_th
c       write(*,*)'threR', synchro_receive_th
c       write(*,'(a,6i4)')'dm_th',ind_dm_thread
c       write(*,'(a,6i4)')'dm_sk',ind_dm_sock
c      endif

      !correction eventuelle au bord pour calcul dans 1 rangee maille fictive
      if(extended_range.eq.1) then

       if( ind_sdm(1).eq.1) ind_sdm(1)= 1-nijk(4)+1
       if( ind_sdm(3).eq.1) ind_sdm(3)= 1-nijk(4)+1
       if( ind_sdm(5).eq.1) ind_sdm(5)= 1-nijk(5)+inck

      if(ind_sdm(2).eq.ind_dm_zone(2))ind_sdm(2)=ind_sdm(2)+nijk(4)-1
      if(ind_sdm(4).eq.ind_dm_zone(4))ind_sdm(4)=ind_sdm(4)+nijk(4)-1
      if(ind_sdm(6).eq.ind_dm_zone(6))ind_sdm(6)=ind_sdm(6)+nijk(5)-inck
      endif

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

           !1er bloc cache
           If(l.eq.1) Then

             shift_coe1  =-ific
             shift_coe2  = ific
             shift_gra1  =-inc
             shift_gra2  = inc
             shift_mjr1  = 0
             shift_mjr2  =-inc
             shift_ssa1  =-inc
             shift_ssa2  = 0

             !test pour eventuel optim  du init RHS
             !(suppression debordement a gauche au bord)
             c1 =-inc
             if(ind_dm_thread(ideb).eq.ind_dm_sock(ideb)) 
     &           c1 =  -synchro_send_sock(i)*inc

             shift_rhs1  = c1
             shift_rhs2  = 0

             !Si 1er bloc cache = dernier
             if( l.eq.ijkv_sdm(i)) then

               !Si  bloc cache recoit des donnee du bloc a droite
               if(    (synchro_receive_sock(i).eq.1)
     &            .or.(topo_th(i).ne.thread_pos(i))   )  then 
                 shift_coe2  =-ific
                 shift_gra2  =-inc
                 shift_rhs2  =-inc
                 shift_mjr2  = 0
                 shift_ssa2  =-inc
               !Sinon: 
               !  init rhs sur derniere cellule reelle en laminaire et calcul coe 
               !  init rhs sur 1er ghostcell en SA +calcul coe(6) 
               else
                 shift_mjr2=0
                 shift_ssa2=inc
               endif

               !if(topo_th(i).eq.1) then
               !  shift_mjr2=0
               !  shift_ssa2=inc
               !endif

             endif

           !Si dernier bloc cache
           Elseif(l.eq.ijkv_sdm(i)) then
 
             shift_coe1  = ific
             shift_coe2  = ific
             shift_gra1  = inc
             shift_gra2  = inc
             shift_rhs1  =   0
             shift_rhs2  =   0
             shift_mjr1  =-inc
             shift_mjr2  =   0
             shift_ssa1  =   0
             shift_ssa2  =   0
            
             if(   (synchro_receive_sock(i).eq.1)
     &         .or.(topo_th(i).ne.thread_pos(i) )     )  then 
                 shift_coe2  =-ific
                 shift_gra2  =-inc
                 shift_rhs2  =-inc
                 shift_ssa2  =-inc
             else
                 shift_mjr2 = 0
                 shift_ssa2 = inc
             endif
             !if(topo_th(i).eq.1) then
             !     shift_mjr2 = 0
             !     shift_ssa2 = inc
             !endif

           Else
             shift_coe1  = ific
             shift_coe2  = ific
             shift_gra1  = inc
             shift_gra2  = inc
             shift_rhs1  =   0
             shift_rhs2  =   0
             shift_mjr1  =-inc
             shift_mjr2  =-inc
             shift_ssa1  =   0
             shift_ssa2  =   0

             if(   (synchro_receive_sock(i).eq.1)
     &         .or.(topo_th(i).ne.thread_pos(i) )     )  then 

                 delta       = ind_dm_thread(ifin) - ind_sdm(ifin)
                 if(delta.lt.4) shift_coe2  = delta - ific
                 if(delta.lt.2) shift_gra2  = delta - inc

             endif
           Endif


           ind_coe(ideb)  = ind_sdm(ideb) + shift_coe1
           ind_grad(ideb) = ind_sdm(ideb) + shift_gra1
           ind_rhs(ideb)  = ind_sdm(ideb) + shift_rhs1
           ind_mjr(ideb)  = ind_sdm(ideb) + shift_mjr1
           ind_ssa(ideb)  = ind_sdm(ideb) + shift_ssa1
           ind_coe(ifin)  = ind_sdm(ifin) + shift_coe2
           ind_grad(ifin) = ind_sdm(ifin) + shift_gra2
           ind_rhs(ifin)  = ind_sdm(ifin) + shift_rhs2
           ind_mjr(ifin)  = ind_sdm(ifin) + shift_mjr2
           ind_ssa(ifin)  = ind_sdm(ifin) + shift_ssa2

        enddo !! boucle dir I,J,K

         if(kfludom.eq.7)then
         if(jc.eq.1)          ind_mjr(3)=ind_sdm(3)   
         if(jc.eq.ijkv_sdm(2))ind_mjr(4)=ind_sdm(4)
         if(kc.eq.1)          ind_mjr(5)=ind_sdm(5)   
         if(kc.eq.ijkv_sdm(3))ind_mjr(6)=ind_sdm(6)
         endif



      end
