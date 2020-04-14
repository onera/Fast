c***********************************************************************
c     $Date: 2010-11-04 13:25:50 +0100 (Thu, 04 Nov 2010) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine topo_scater(ndo, ithread, socket, omp_mode, lmin,
     &                       kfludom,
     &                       thread_parsock, thread_parsock_actif,
     &                       ithread_sock, socket_topology,
     &                       size_target, ind_dm_zone, topo_omp,
     &                       topo_s, ijkvloc, thread_pos, socket_pos, 
     &                       thread_topology, size_thread, ijkv_thread)

      implicit none

      INTEGER_E ndo,ithread, socket,thread_parsock, lmin,
     & ithread_sock, socket_topology(3),thread_parsock_actif,
     & size_target(3), ind_dm_zone(6), topo_s(3),topo_omp(3),
     & ijkvloc(3), thread_pos(3), socket_pos(3), thread_topology(3),
     & size_thread(3), ijkv_thread(3),kfludom, omp_mode

C Var loc 
      logical new_try,lskip,lpair
      INTEGER_E compteur, i,l,IJKmin(3),topo_target(3),res,iv,
     & ind_dm_loc(6)

      !Par defaut tous les threads bossent
      thread_parsock_actif = thread_parsock
      


      !!Bypass de la routine pour forcer le mode mono-bloc si Nombre de
      !Socket =1
      i = socket_topology(1)*socket_topology(2)*socket_topology(3)
      if(kfludom.eq.7) i=2
      If(i.eq.1) then

         if(omp_mode.eq.0) then
           call indice_boucle_lu(ndo, ithread, thread_parsock, lmin,
     &                           ind_dm_zone, 
     &                           topo_s, ind_dm_loc )
         else
             topo_s(1:3)     = topo_omp(1:3)
             !ind_dm_loc[1:6] = inddm_omp[1:6] tableau local: affectation inutile. voir indice_boucle_scater
         endif

         thread_parsock_actif = topo_s(1)*topo_s(2)*topo_s(3)

         do i =1,3
         iv             =  ind_dm_zone(2*i)-ind_dm_zone(2*i-1) +1
         res            = mod( iv , topo_s(i) )
         ijkvloc(i)     = ( iv - res)/ topo_s(i)
         size_target(i) =  ijkvloc(i)
         IJKmin(i)      = topo_s(i)

         enddo


      Else



      compteur   = 0
      lpair=.false.
      if(mod(thread_parsock,2).eq.0) lpair=.true.


      !!on determine la socket la plus petite afin de garantir que
      !toutes les socket auront la meme topology de thread
      do i =1,3
         iv         =  ind_dm_zone(2*i)-ind_dm_zone(2*i-1) +1
         res        = mod( iv , socket_topology(i) )
         ijkvloc(i) = ( iv - res)/ socket_topology(i)
      enddo

 1000 continue
      new_try=.false.
      lskip  =.false.

      IJKmin(1) =  ijkvloc(1)/size_target(1)
      IJKmin(2) =  ijkvloc(2)/size_target(2)
      IJKmin(3) =  ijkvloc(3)/size_target(3)

      !if(ithread.eq.1) write(*,*)'IJKmin',IJKmin, compteur
      !!on cherche si les parametre target sont compatibles avec la
      !taille de la zone

      !priorite decoupe K
      if( IJKmin(3).ge.thread_parsock) then

           topo_s(1)= 1
           topo_s(2)= 1
           topo_s(3)= thread_parsock

      !priorite decoupe K et J si possible
      elseif( IJKmin(3).ge.thread_parsock/2.and.
     &        IJKmin(2).ge.2.and.lpair)     then

           topo_s(1)= 1
           topo_s(2)= 2
           topo_s(3)= thread_parsock/2

      !decoupe J
      elseif( IJKmin(2).ge.thread_parsock.and.compteur.ge.1) then

           topo_s(1)      = 1
           topo_s(2)      = thread_parsock
           topo_s(3)      = 1

      !priorite decoupe J et I si possible
      elseif( IJKmin(2).ge.thread_parsock/2.and.
     &        IJKmin(1).ge.2.and.lpair.and.compteur.ge.1)     then

           topo_s(1)= 2
           topo_s(2)= thread_parsock/2
           topo_s(3)= 1

      elseif( IJKmin(1).ge.thread_parsock.and.compteur.gt.1) then

           topo_s(1)      = thread_parsock
           topo_s(2)      = 1
           topo_s(3)      = 1
        
      else
           if(compteur .lt. 6) new_try =.true.
           lskip = .true.
      endif

      !if(ithread.eq.1) write(*,*)'size_target',size_target,compteur

      compteur = compteur + 1
      if(new_try) then
         size_target(1) = max(32, size_target(1)/2)   !taille minimale vecto
         size_target(2) = max(8 , size_target(2)/2)
         size_target(3) = max(8 , size_target(3)/2)     !taille minimal pour synchro omp
         !size_target(2) = max(2 , size_target(2)/2)
         !size_target(3) = max(2 , size_target(3)/2)     !taille minimal pour synchro omp
        goto 1000
      endif


      if(lskip) then

         size_target(1) = min(32, ijkvloc(1) )   !taille minimale vecto
         size_target(2) = 2
         size_target(3) = 2

         IJKmin(1) =  ijkvloc(1)/size_target(1)
         IJKmin(2) =  ijkvloc(2)/size_target(2)
         IJKmin(3) =  ijkvloc(3)/size_target(3)

         topo_s(1) = ijkvloc(1)/size_target(1)
         topo_s(2) = ijkvloc(2)/size_target(2)
         topo_s(3) = ijkvloc(3)/size_target(3)

         if(topo_s(3).gt.topo_s(1)) then

           topo_s(1) = 1
           thread_parsock_actif = topo_s(3)*topo_s(2)

           res= thread_parsock - thread_parsock_actif
           if(res.lt.0) then
            
              if(topo_s(3).gt.topo_s(2)) then
                
                topo_s(3)      = max(1,topo_s(3) + (res/topo_s(3)-1) )
                size_target(3) = ijkvloc(3) / topo_s(3)
                IJKmin(3)      = ijkvloc(3)/size_target(3)
              else
                topo_s(2)      = max(1,topo_s(2) + (res/topo_s(2)-1) )
                size_target(2) = ijkvloc(2) / topo_s(2)
                IJKmin(2)      = ijkvloc(2)/size_target(2)
              endif
              thread_parsock_actif = topo_s(1)*topo_s(2)
           endif

         else
           topo_s(3) = 1
           thread_parsock_actif = topo_s(1)*topo_s(2)

           res= thread_parsock - thread_parsock_actif
           if(res.lt.0) then
            
              if(topo_s(2).gt.topo_s(2)) then
                
                topo_s(2)      = max(1,topo_s(2) + (res/topo_s(2)-1) )
                size_target(2) = ijkvloc(2) / topo_s(2)
                IJKmin(2)      = ijkvloc(2)/size_target(2)
              else
                topo_s(1)      = max(1,topo_s(1) + (res/topo_s(1)-1) )
                size_target(1) = ijkvloc(1) / topo_s(1)
                IJKmin(1)      = ijkvloc(1)/size_target(1)
              endif
              thread_parsock_actif = topo_s(1)*topo_s(2)
           endif
         endif

      endif

      Endif !!socket =1 ou > 1


c      if(ithread.eq.1)then
c        write(*,*)'topooo',topo_s,thread_parsock_actif,ndo
c       write(*,*)'topo',topo_s
c       write(*,*)'size', size_target
c       write(*,*)'IJK', IJKmin
c      endif

      do i =1,3
            thread_topology(i) =topo_s(i)
            size_thread(i)     = size_target(i)
            ijkv_thread(i)     = max(1,IJKmin(i)/topo_s(i))
      enddo
      
      i              = socket_topology(1)*socket_topology(2)
      socket_pos(3)  = 1+ (socket-1)/i
      l              = socket -(socket_pos(3)-1)*i
      socket_pos(2)  = 1 + (l-1)/socket_topology(1)
      socket_pos(1)  = l - (socket_pos(2)-1)*socket_topology(1)

      ithread_sock   = ithread-(socket-1)*thread_parsock

      thread_pos(3) = 1 + (ithread_sock-1)/(topo_s(1)*topo_s(2))
      l    = ithread_sock -(thread_pos(3)-1)*topo_s(1)*topo_s(2)
      thread_pos(2) = 1 + (l-1)/topo_s(1)
      thread_pos(1) = l - (thread_pos(2)-1)*topo_s(1)

      end
